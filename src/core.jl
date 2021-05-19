

struct BLSPeriodogram{TT,VT<:AbstractVector{TT},YT<:AbstractVector,ET<:AbstractVector}
    t::VT
    y::YT
    err::ET
end


function autoperiod(bls::BLSPeriodogram,
                    duration;
                    minimum_period=default_minimum_period(duration),
                    maximum_period=default_maximum_period(bls, minimum_n_transit),
                    minimum_n_transit=3,
                    frequency_factor=1.0)

    lo, hi = extrema(bls.t)
    baseline = hi - lo
    df = frequency_factor * minimum(duration) / baseline^2

    if maximum_period < minimum_period
        minimum_period, maximum_period = maximum_period, minimum_period
    end
    minimum_period ≥ 0 || error("minimum period must be positive")

    minimum_frequency = inv(maximum_period)
    maximum_frequency = inv(minimum_period)

    nf = (maximum_frequency - minimum_frequency) ÷ df
    return @. maximum_frequency - df * (0:nf)
end

default_minimum_period(duration) = 2 * maximum(duration)
function default_maximum_period(bls, minimum_n_transit)
    lo, hi = extrema(bls.t)
    return (hi - lo) / (minimum_n_transit - 1)
end