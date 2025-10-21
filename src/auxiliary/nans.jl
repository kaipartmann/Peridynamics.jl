function containsnan(K::T) where {T<:AbstractArray}
    @simd for i in eachindex(K)
        isnan(K[i]) && return true
    end
    return false
end

function nancheck(storage::AbstractStorage, t, Δt)
    if containsnan(storage.b_int)
        n = (t > 0 && Δt > 0) ? Int(t ÷ Δt) : 0
        msg = "NaN's found in field `b_int` at simulation time $(t), step $(n)!\n"
        error(msg)
    end
    return nothing
end
