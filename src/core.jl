

struct BLSPeriodogram{TT,VT<:AbstractVector{TT},YT<:AbstractVector,ET<:AbstractVector,PT,PWT,DT,TOT}
    t::VT
    y::YT
    err::ET
    period::PT
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
    print(io, "BLSPeriodogram\n==============\ninput dim: ", n, "\noutput dim: ", m, "\n\nparameters\n----------\n", "period: ", p.period, "\nduration: ", p.duration, "\nt0: ", p.t0, "\n$method: ", p.power)
end


function autoperiod(t,
                    duration;
                    minimum_n_transit=3,
                    minimum_period=default_minimum_period(duration),
                    maximum_period=default_maximum_period(t, minimum_n_transit),
                    frequency_factor=1.0)

    lo, hi = extrema(t)
    baseline = hi - lo
    df = frequency_factor * minimum(duration) / baseline^2

    if maximum_period < minimum_period
        minimum_period, maximum_period = maximum_period, minimum_period
    end
    minimum_period ≥ zero(minimum_period) || error("minimum period must be positive")

    minimum_frequency = inv(maximum_period)
    maximum_frequency = inv(minimum_period)

    nf = 1 + (maximum_frequency - minimum_frequency) ÷ df
    return @. inv(maximum_frequency - df * (0:nf))
end

default_minimum_period(duration) = 2 * maximum(duration)
function default_maximum_period(t, minimum_n_transit)
    minimum_n_transit > 1 || error("minimum number of transits must be greater than 1")
    lo, hi = extrema(t)
    return (hi - lo) / (minimum_n_transit - 1)
end

function BLS(t, y, yerr;
               duration,
               objective=:likelihood,
               oversample=10,
               period=autoperiod(t, duration))
    invvar = @. inv(yerr^2)
    y_res = y .- median(y)
    min_t = minimum(t)

    # set up arrays
    powers = similar(period, Float64)
    durations = similar(period, typeof(duration))
    t0s = similar(period)

    @inbounds Threads.@threads for idx in eachindex(period)
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
    return BLSPeriodogram(t, y, yerr, period, powers, durations, t0s, objective)
end
