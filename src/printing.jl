
function Base.show(io::IO, bls::BLSPeriodogram{T,S}) where {T,S}
    N = length(bls.power)
    print(io, "$N-element BLSPeriodogram{$T,$S}")
end

function Base.show(io::IO, ::MIME"text/plain", bls::BLSPeriodogram)
    p = params(bls)
    println(io, "BLSPeriodogram\n==============")
    println(io, "input dim: ", length(bls.t))
    println(io, "output dim: ", length(bls.power))
    println(io, "period range: ", range_str(extrema(bls.period)))
    println(io, "duration range: ", range_str(extrema(bls.duration_in)))
    println(io, "objective: ", bls.method)
    println(io, "\nparameters\n----------")
    println(io, "index: ", p.index)
    println(io, "period: ", p.period)
    println(io, "duration: ", p.duration)
    println(io, "t0: ", p.t0)
    println(io, "depth: ", p.depth, " ± ", p.depth_err)
    println(io, "snr: ", p.snr)
    print(io, "log-likelihood: ", p.loglike)
end

range_str((min, max)) = "$min - $max"