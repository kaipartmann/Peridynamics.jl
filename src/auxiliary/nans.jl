function containsnan(K::T) where {T<:AbstractArray}
    @simd for i in eachindex(K)
        isnan(K[i]) && return true
    end
    return false
end

containsnan(storage::AbstractStorage) = containsnan(storage.b_int)

function check_for_nans(chunk::AbstractBodyChunk, t, Δt)
    check_for_nans(chunk.storage, t, Δt)
    return nothing
end

function check_for_nans(storage::AbstractStorage, t, Δt)
    containsnan(storage) && throw(NaNError(t, calculate_current_step(t, Δt)))
    return nothing
end

function check_for_nans_mpi(chunk::AbstractBodyChunk, t, Δt)
    check_for_nans_mpi(chunk.storage, t, Δt)
    return nothing
end

function check_for_nans_mpi(storage::AbstractStorage, t, Δt)
    rank_containsnan = containsnan(storage)
    someone_containsnan = MPI.Allreduce(rank_containsnan, |, mpi_comm())
    someone_containsnan && throw(NaNError(t, calculate_current_step(t, Δt)))
    return nothing
end

calculate_current_step(t, Δt) = (t > 0 && Δt > 0) ? Int(t ÷ Δt) : 0
