using LoopVectorization

"""
    BLSPeriodogram

A convenient wrapper for outputs from [`BLS`](@ref).

# Methods
* [`BoxLeastSquares.params`](@ref)
* [`BoxLeastSquares.power`](@ref)
* [`BoxLeastSquares.periods`](@ref)

# Attributes
* `t` - input time grid
* `y` - input data
* `yerr` - input data uncertainty
* `periods` - the input periods
* `duration_in` - the input durations
* `objective` - the objective that was maximized
* `power` - the power calculated at each period
* `duration` - the best duration at each period
* `t0` - the best transit time at each period
* `depth` - the best transit depth at each period
* `snr` - the signal-to-noise ratio at each period
* `loglike` - the log-likeilhood at each period

# Plotting

Plotting recipes are provided for `BLSPeriodogram` which automatically plots the period and the power
"""
struct BLSPeriodogram{TT,FT,FET,PT,DT,VTT<:AbstractVector{TT},VFT<:AbstractVector{FT},VFET<:AbstractVector{FET},T1,T2,T3,T4,T5,T6}
    t::VTT
    y::VFT
    yerr::VFET
    periods::PT
    duration_in::DT
    objective::Symbol

    power::T1
    duration::T2
    t0::T3
    depth::T4
    snr::T5
    loglike::T6
end

"""
    BoxLeastSquares.periods(::BLSPeriodogram)

Return the period grid for the periodogram
"""
periods(bls::BLSPeriodogram) = bls.periods
"""
    BoxLeastSquares.power(::BLSPeriodogram)

Return the power calculated for each period for the periodogram
"""
power(bls::BLSPeriodogram) = bls.power

"""
    BoxLeastSquares.params(::BLSPeriodogram)

Return the transit parameters for the best fitting period. Returns period, duration, t0, and power as well as the index of the max-power period.
"""
function params(bls::BLSPeriodogram)
    ind = argmax(power(bls))
    pow = power(bls)[ind]
    per = periods(bls)[ind]
    dur = bls.duration[ind]
    t0 = bls.t0[ind]
    depth = bls.depth[ind]
    snr = bls.snr[ind]
    depth_err = depth / snr
    loglike = bls.loglike[ind]
    return (index=ind, power=pow, period=per, duration=dur, t0=t0,
            depth=depth, depth_err=depth_err, snr=snr, loglike=loglike)
end

"""
    autoperiod(t, duration;
        minimum_n_transit=3, frequency_factor=1.0,
        [minimum_period, maximum_period])

Automatically determine a period grid from the given times and duration(s). Periods are selected such that at least `minimum_n_trasnit` transits occur. The default minimum period is twice the maximum duration. The default maximum period is `(maximum(t) - minimum(t)) / (minimum_n_transit - 1)`. The frequency factor changes the granularity in frequency space- a smaller frequency factor will create a finer period grid.
"""
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

@inline wrap(x, period) = x - period * trunc(x / period)

@inline function wrap_index(t::Real, min_t, P, bin_duration)
    unsafe_trunc(Int, wrap(t - min_t, P) / bin_duration) + 1
end
# this method will catch e.g. AbstractQuantity <: Number
# which is important because unsafe_trunc not defined
@inline function wrap_index(t, min_t, P, bin_duration)
    trunc(Int, wrap(t - min_t, P) / bin_duration) + 1
end

"""
    BLS(t, y, [yerr];
        duration, periods=autoperiod(t, duration, kwargs...), 
        objective=:likelihood, oversample=10, kwargs...)

Compute the box-least-squares periodogram.

# Parameters
* `t` - the time for each observation. Units are irrelevant, except that they must be consistent for all temporal parameters (e.g., `duration`). `Unitful.jl` units work seamlessly without needing to convert.
* `y` - the flux value for each observation
* `yerr`, optional - the uncertainty for each observation, if not provided, will default to ones
* `duration` - The duration or durations to consider. Same units as `t`
* `periods`, optional - The period grid to computer the BLS power over. If not provided, [`autoperiod`](@ref) will be called along with any extra keyword arguments (like `minimum_period`)
* `objective`, optional - Choose between maximizing the likeilhood (`:likeilhood`, default) or the signal-to-noise ratio (`:snr`).
* `oversample`, optional - The number of bins per duration that should be used. Larger values of `oversample` will lead to a finer grid.

The returned values are wrapped into a convenience type [`BoxLeastSquares.BLSPeriodogram`](@ref)
"""
function BLS(t, y, yerr=fill!(similar(y), one(eltype(y)));
    duration,
    periods=nothing,
    objective=:likelihood,
    oversample=10,
    kwargs...)
    if isnothing(periods)
        periods = autoperiod(t, duration; kwargs...)
    end

    # set up arrays
    powers = similar(periods, Float64)
    durations = similar(powers, eltype(duration))
    t0s = similar(powers, eltype(t))
    depths = similar(powers, eltype(y))
    snrs = similar(powers)
    loglikes = similar(powers)

    # find parameter ranges
    max_P = maximum(periods)
    min_τ = minimum(duration)
    bin_duration = min_τ / oversample
    max_bins = ceil(Int, max_P / bin_duration) + oversample
    # work arrays
    mean_y = similar(y, typeof(inv(first(y))), max_bins + 1)
    mean_ivar = similar(yerr, typeof(inv(first(yerr)^2)), max_bins + 1)
    _begin = firstindex(mean_y) # need this because @turbo can't use implicit begin

    # pre-accumulate some factors
    ymed = median(y)
    min_t = minimum(t)

    sum_y = zero(eltype(mean_y))
    sum_ivar = zero(eltype(mean_ivar))
    @inbounds for i in eachindex(y, yerr)
        iv = inv(yerr[i]^2)
        sum_y += (y[i] - ymed) * iv
        sum_ivar += iv
    end

    # loop over periods in grid search
    @inbounds for idx in eachindex(periods, powers, durations, t0s, depths, snrs, loglikes)
        P = periods[idx]
        n_bins = ceil(Int, P / bin_duration) + oversample
        @turbo for n in 0:n_bins
            mean_y[_begin + n] = zero(eltype(mean_y))
            mean_ivar[_begin + n] = zero(eltype(mean_ivar))
        end

        @turbo for idx in eachindex(t, y, yerr)
            ind = wrap_index(t[idx], min_t, P, bin_duration)
            iv = inv(yerr[idx]^2)
            mean_y[_begin + ind] += (y[idx] - ymed) * iv
            mean_ivar[_begin + ind] += iv
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
        powers[idx] = -Inf
        for dur in duration
            τ = round(Int, dur / bin_duration)
            for n in 0:n_bins - τ
                # estimate in- and out-of-transit flux
                y_in = mean_y[begin + n + τ] - mean_y[begin + n]
                ivar_in = mean_ivar[begin + n + τ] - mean_ivar[begin + n]
                y_out = sum_y - y_in
                ivar_out = sum_ivar - ivar_in
                # check if no points in transit
                (iszero(ivar_in) || iszero(ivar_out)) &&  continue
                # normalize
                y_in /= ivar_in
                y_out /= ivar_out
                # check if negative SNR
                y_out < y_in && continue

                # compute best-fit depth
                depth = y_out - y_in
                depth_err = sqrt(inv(ivar_in) + inv(ivar_out))
                # calculate objectives
                snr = depth / depth_err
                loglike = 0.5 * ivar_in * (y_out - y_in)^2

                if objective === :snr
                    power = snr
                elseif objective === :likelihood
                    power = loglike
                else
                    error("invalid objective $objective. Should be one of `:snr` or `:likelihood`")
                end

                if power > powers[idx]
                    powers[idx] = power
                    durations[idx] = τ * bin_duration
                    t0s[idx] = mod(n * bin_duration + 0.5 * durations[idx] + min_t, P)
                    depths[idx] = depth
                    snrs[idx] = snr
                    loglikes[idx] = loglike
                end
            end
        end
    end

    return BLSPeriodogram(t, y, yerr, periods, duration, objective, powers, durations, t0s, depths, snrs, loglikes)
end

