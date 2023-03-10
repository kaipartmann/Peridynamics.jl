struct MultiMaterial{T<:AbstractPDMaterial, N} # TODO: rename to MaterialMap?
    materials::NTuple{N,T}
    matofpoint::Vector{UInt8}

    function MultiMaterial(
        materials::NTuple{N,T},
        matofpoint::AbstractVector{I},
    ) where {N, T<:AbstractPDMaterial, I<:Integer}

        if N == 1
            return new{T,N}(materials, Vector{UInt8}())
        end

        if !allequal(1 .<= matofpoint .<= N)
            error("Index in `matofpoint` out of bounds!")
        end
        if N <= typemax(UInt8)
            _matofpoint = convert.(UInt8, matofpoint)
        elseif N <= typemax(UInt16)
            _matofpoint = convert.(UInt16, matofpoint)
        elseif N <= typemax(UInt32)
            _matofpoint = convert.(UInt32, matofpoint)
        elseif N <= typemax(UInt64)
            _matofpoint = convert.(UInt64, matofpoint)
        end
        new{T,N}(materials, _matofpoint)
    end
end

function Base.show(io::IO, ::MIME"text/plain", mm::MultiMaterial)
    print(io, typeof(mm), " for ", length(mm.matofpoint), " points")
    return nothing
end

const PDMaterial{T} = Union{T, MultiMaterial{T,N}} where {T<:AbstractPDMaterial,N}

Base.getindex(mm::MultiMaterial, i::Int) = mm.materials[mm.matofpoint[i]]
function Base.getindex(mm::MultiMaterial, I::AbstractVector{T}) where {T<:Integer}
    return mm.materials[mm.matofpoint[I]]
end
Base.getindex(mm::MultiMaterial, ::Colon) = mm.materials[mm.matofpoint]
Base.getindex(mm::MultiMaterial{T,1}, ::Int) where {T} = mm.materials[1]
Base.firstindex(mm::MultiMaterial) = 1
Base.lastindex(mm::MultiMaterial) = length(mm.matofpoint)

Base.eltype(::MultiMaterial{T}) where {T} = T

Base.getindex(mat::AbstractPDMaterial, ::Int) = mat
Base.getindex(mat::AbstractPDMaterial, ::AbstractVector{T}) where {T<:Integer} = mat
Base.getindex(mat::AbstractPDMaterial, ::Colon) = mat
Base.firstindex(::AbstractPDMaterial) = 1
Base.lastindex(::AbstractPDMaterial) = 1

Base.eltype(::T) where{T<:AbstractPDMaterial} = T
