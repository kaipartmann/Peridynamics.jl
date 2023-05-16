struct MultiMaterial{T<:AbstractPDMaterial,N,I<:Integer}
    materials::NTuple{N,T}
    matofpoint::Vector{I}

    function MultiMaterial(
        materials::NTuple{N,T}, matofpoint::AbstractVector{I}
    ) where {T<:AbstractPDMaterial,N,I<:Integer}
        if N == 1
            error("N = 1\nUsage of MultiMaterial only makes sense for N > 1")
        end

        if all(1 .<= matofpoint .<= N)
            error("Index in `matofpoint` out of bounds!")
        end

        # save memory for small N
        if N <= typemax(UInt8)
            _matofpoint = convert.(UInt8, matofpoint)
        else
            _matofpoint = matofpoint
        end

        return new{T,N,eltype(_matofpoint)}(materials, _matofpoint)
    end
end

function Base.show(io::IO, ::MIME"text/plain", mm::MultiMaterial)
    print(io, typeof(mm), " for ", length(mm.matofpoint), " points")
    return nothing
end

const PDMaterial{T} = Union{T,MultiMaterial{T,N}} where {T<:AbstractPDMaterial,N}

Base.getindex(mm::MultiMaterial, i::Int) = mm.materials[mm.matofpoint[i]]
function Base.getindex(mm::MultiMaterial, I::AbstractVector{T}) where {T<:Integer}
    return mm.materials[mm.matofpoint[I]]
end
Base.getindex(mm::MultiMaterial, ::Colon) = mm.materials[mm.matofpoint]
Base.firstindex(mm::MultiMaterial) = 1
Base.lastindex(mm::MultiMaterial) = length(mm.matofpoint)

Base.eltype(::MultiMaterial{T}) where {T} = T

Base.getindex(mat::AbstractPDMaterial, ::Int) = mat
Base.getindex(mat::AbstractPDMaterial, ::AbstractVector{T}) where {T<:Integer} = mat
Base.getindex(mat::AbstractPDMaterial, ::Colon) = mat
Base.firstindex(::AbstractPDMaterial) = 1
Base.lastindex(::AbstractPDMaterial) = 1

Base.eltype(::T) where {T<:AbstractPDMaterial} = T
