

struct BLSPeriodogram{TT,VT<:AbstractVector{TT},YT<:AbstractVector,ET<:AbstractVector,PT,DIT,PWT,DT,TOT}
    t::VT
    y::YT
    err::ET
    period::PT
    duration_in::DIT
    power::PWT
    duration::DT
    t0::TOT
    method::Symbol
end

period(bls::BLSPeriodogram) = bls.period
duration(bls::BLSPeriodogram) = bls.duration
power(bls::BLSPeriodogram) = bls.power
t0(bls::BLSPeriodogram) = bls.t0

"""
    params(::BLSPeriodogram)

Return the transit parameters for the best fitting period. Returns period, duration, t0, and power.
"""
function params(bls::BLSPeriodogram)
    best_ind = argmax(power(bls))
    best_pow = power(bls)[best_ind]
    best_per = period(bls)[best_ind]
    best_dur = duration(bls)[best_ind]
    best_t0 = t0(bls)[best_ind]
    return (power=best_pow, period=best_per, duration=best_dur, t0=best_t0)
end

function Base.show(io::IO, ::MIME"text/plain", bls::BLSPeriodogram)
    p = params(bls)
    method = string(bls.method)
    n = length(bls.t)
    m = length(bls.power)
    println(io, "BLSPeriodogram\n==============")
    println(io, "input dim: ", n)
    println(io, "output dim: ", m)
    println(io, "period range: ", range_str(extrema(bls.period)))
    println(io, "duration range: ", range_str(extrema(bls.duration_in)))
    println(io, "\nparameters\n----------")
    println(io, "period: ", p.period)
    println(io, "duration: ", p.duration)
    println(io, "t0: ", p.t0)
    print(io, "$method: ", p.power)
end

range_str((min, max)) = "$min - $max"


function autoperiod(t::AbstractVector{T},
                    duration;
                    minimum_n_transit=3,
                    minimum_period=default_minimum_period(duration),
                    maximum_period=default_maximum_period(t, minimum_n_transit),
                    frequency_factor=1.0) where T
    baseline = maximum(t) - minimum(t)
    df = frequency_factor * minimum(duration) / baseline^2

    if maximum_period < minimum_period
        minimum_period, maximum_period = maximum_period, minimum_period
    end
    minimum_period ≥ zero(T) || error("minimum period must be positive")

    minimum_frequency = inv(maximum_period)
    maximum_frequency = inv(minimum_period)

    nf = 1 + (maximum_frequency - minimum_frequency) ÷ df
    return @. convert(T, inv(maximum_frequency - df * (0:nf)))
end

default_minimum_period(duration) = 2 * maximum(duration)
function default_maximum_period(t, minimum_n_transit)
    minimum_n_transit > 1 || error("minimum number of transits must be greater than 1")
    baseline = maximum(t) - minimum(t)
    return baseline / (minimum_n_transit - 1)
end

@inline wrap(x, period) = x - period * floor(x / period)

function BLS(t, y, yerr;
    duration,
    objective=:likelihood,
    oversample=10,
    period=autoperiod(t, duration))

    # set up arrays
    powers = similar(period, Float64)
    durations = similar(period, eltype(duration))
    t0s = similar(period)

    # find parameter ranges
    max_P = maximum(period)
    min_τ = minimum(duration)
    bin_duration = min_τ / oversample
    max_bins = ceil(Int, max_P / bin_duration) + oversample
    # work arrays
    mean_y = similar(y, typeof(inv(first(y))), max_bins + 1)
    mean_ivar = similar(yerr, typeof(inv(first(yerr)^2)), max_bins + 1)

    # pre-accumulate some factors
    ymed = median(y)
    min_t = minimum(t)

    sum_y = zero(eltype(mean_y))
    sum_ivar = zero(eltype(mean_ivar))
    for i in eachindex(y)
        iv = inv(yerr[i]^2)
        sum_y += (y[i] - ymed) * iv
        sum_ivar += iv
    end

    # loop over periods in grid search
    @inbounds for idx in eachindex(period)
        P = period[idx]
        n_bins = ceil(Int, P / bin_duration) + oversample
        for n in 0:n_bins
            mean_y[begin + n] = zero(eltype(mean_y))
            mean_ivar[begin + n] = zero(eltype(mean_ivar))
        end

        for n in 0:length(t) - 1
            ind = round(Int, wrap(t[begin + n] - min_t, P) / bin_duration) + 1
            iv = inv(yerr[begin + n]^2)
            mean_y[begin + ind] += (y[begin + n] - ymed) * iv
            mean_ivar[begin + ind] += iv
        end

        # simplify calcualations by wrapping binned values
        ind = n_bins - oversample
        for n in 1:oversample
            mean_y[begin + ind] = mean_y[begin + n]
            mean_ivar[begin + ind] = mean_ivar[begin + n]
            ind += 1
        end

        for n in 1:n_bins
            mean_y[begin + n] += mean_y[begin + n - 1]
            mean_ivar[begin + n] += mean_ivar[begin + n - 1]
        end

        # now, loop over phases
        best = (power=-Inf, duration=zero(P), t0=zero(P))
        for k in 0:length(duration) - 1
            τ = round(Int, duration[begin + k] / bin_duration)
            n_max = n_bins - τ
            for n in 0:n_max
                # estimate in- and out-of-transit flux
                y_in = mean_y[begin + n + τ] - mean_y[begin + n]
                ivar_in = mean_ivar[begin + n + τ] - mean_ivar[begin + n]
                y_out = sum_y - y_in
                ivar_out = sum_ivar - ivar_in
                # check if no points in transit
                if iszero(ivar_in) || iszero(ivar_out)
                    continue
                end
                # normalize
                y_in /= ivar_in
                y_out /= ivar_out

                # compute best-fit depth
                if objective === :snr
                    depth = y_out - y_in
                    depth_err = sqrt(inv(ivar_in) + inv(ivar_out))
                    power = depth / depth_err
                elseif objective === :likelihood
                    power = 0.5 * ivar_in * (y_out - y_in)^2
                else
                    error("invalid objective $objective. Should be one of `:snr` or `:likelihood`")
                end

                if power > best.power
                    dur = τ * bin_duration
                    t0 = mod(n * bin_duration + 0.5 * dur, P)
                    best = (power=power, duration=dur, t0=t0)
                end
            end
        end
        powers[idx] = best.power
        durations[idx] = best.duration
        t0s[idx] = best.t0 + min_t
    end
    return BLSPeriodogram(t, y, yerr, period, duration, powers, durations, t0s, objective)
end


function BLS_slow(t, y, yerr;
    duration,
    objective=:likelihood,
    oversample=10,
    period=autoperiod(t, duration))
    invvar = @. inv(yerr^2)
    y_res = y .- median(y)
    min_t = minimum(t)

    # set up arrays
    powers = similar(period, Float64)
    durations = similar(period, eltype(duration))
    t0s = similar(period)

    @inbounds for idx in eachindex(period)
        P = period[idx]
        best = (power=-Inf, duration=zero(P), t0=zero(P))
        half_P = 0.5 * P
        for τ in duration
            # compute phase grid
            dp = τ / oversample
            phase = zero(P):dp:P

            for t0 in phase
                # find in-transit data points
                invvar_in = invvar_out = zero(eltype(invvar))
                y_in = y_out = zero(eltype(y_res))
                for i in eachindex(t)
                    transiting = abs((t[i] - min_t - t0 + half_P) % P - half_P) < 0.5 * τ
                    if transiting
                        invvar_in += invvar[i]
                        y_in += y_res[i] * invvar[i]
                    else
                        invvar_out += invvar[i]
                        y_out += y_res[i] * invvar[i]
                    end
                end
                y_in /= invvar_in
                y_out /= invvar_out

                # compute best-fit depth
                if objective === :snr
                    depth = y_out - y_in
                    depth_err = sqrt(inv(invvar_in) + inv(invvar_out))
                    snr = depth / depth_err
                    if snr > best.power
                        best = (power=snr, duration=τ, t0=t0)
                    end
                elseif objective === :likelihood
                    loglike = 0.5 * invvar_in * (y_out - y_in)^2
                    if loglike > best.power
                        best = (power=loglike, duration=τ, t0=t0)
                    end
                else
                    error("invalid objective $objective. Should be one of `:snr` or `:likelihood`")
                end
            end
        end
        powers[idx] = best.power
        durations[idx] = best.duration
        t0s[idx] = best.t0 + min_t
    end
    return BLSPeriodogram(t, y, yerr, period, duration, powers, durations, t0s, objective)
end


function model(bls::BLSPeriodogram)
    # compute depth
    pars = params(bls)
    return model(bls.t, bls.y, bls.err; pars...)
end

function model(t, y, yerr; period, duration, t0, kwargs...)
    # compute depth
    hp = 0.5 * period
    y_in = y_out = zero(typeof(inv(first(y))))
    ivar_in = ivar_out = zero(typeof(inv(first(yerr)^2)))
    for idx in eachindex(t, y, yerr)
        transiting = abs((t[idx] - t0 + hp) % period - hp) < 0.5 * duration
        iv = inv(yerr[idx]^2)
        χ = y[idx] * iv
        if transiting
            y_in += χ
            ivar_in += iv
        else
            y_out += χ
            ivar_out += iv
        end
    end
    y_in /= ivar_in
    y_out /= ivar_out
    T = eltype(y)
    y_model = map(t) do t
        transiting = abs((t - t0 + hp) % period - hp) < 0.5 * duration
        convert(T, ifelse(transiting, y_in, y_out))
    end
    return y_model
end