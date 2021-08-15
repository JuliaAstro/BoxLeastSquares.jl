"""
    BoxLeastSquares.model(::BLSPeriodogram; kwargs...)

Create a transit model using the data and best-fitting parameters from the given BLS periodogram. Any keyword parameters can be overriden.
"""
function model(bls::BLSPeriodogram; kwargs...)
    # compute depth
    pars = params(bls)
    return model(bls.t, bls.y, bls.yerr; pars..., kwargs...)
end

"""
    BoxLeastSquares.model(t, y, [yerr]; period, duration, t0)

Evaluate the transit model on the given time grid. If `yerr` is not provided, it will default to 1. The following transit parameters must be set:
* `period` orbital period in the same units as `t`
* `duration` the transit duration in the same units as `t`
* `t0` the transit time (middle of transit) in the same units as `t`

If you are using `Unitful.jl`, the unit conversions will be made automatically.
"""
function model(t, y, yerr=fill!(similar(y), one(eltype(y))); period, duration, t0, kwargs...)
    # compute depth
    hp = 0.5 * period
    y_in = y_out = zero(typeof(first(y)))
    ivar_in = ivar_out = zero(typeof(inv(first(yerr)^2)))
    @inbounds for idx in eachindex(t, y, yerr)
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
    y_model = similar(y)
    @inbounds for idx in eachindex(t)
        transiting = abs((t[idx] - t0 + hp) % period - hp) < 0.5 * duration
        y_model[idx] = ifelse(transiting, y_in, y_out)
    end
    return y_model
end