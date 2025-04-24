function containsnan(K::T) where {T<:AbstractArray}
    @simd for i in eachindex(K)
        isnan(K[i]) && return true
    end
    return false
end

function nancheck(chunk::AbstractBodyChunk, t)
    if containsnan(chunk.storage.b_int)
        msg = "NaN's found in field `b_int` at simulation time $(t)!\n"
        error(msg)
    end
    return nothing
end
