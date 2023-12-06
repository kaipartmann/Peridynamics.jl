"""
    MultiMaterial{T <: AbstractMaterial, N, I <: Integer}

Multiple material properties for one body.

# Fields
- `materials::NTuple{N, T}`: different materials, all of the same type!
- `matofpoint::Vector{I}`: material index in `materials` of each point

# Example
```julia-repl
julia> mat1 = BBMaterial(horizon=1, rho=8e-6, E=2.1e5, Gc=2)
BBMaterial:
  δ:   1.0
  rho: 8.0e-6
  E:   210000.0
  nu:  0.25
  G:   84000.0
  K:   140000.0
  bc:  802140.9131831526
  Gc:  2.0
  εc:  0.0028171808490950554

julia> mat2 = BBMaterial(horizon=1.1, rho=3e-6, E=2.7e4, Gc=1)
BBMaterial:
  δ:   1.1
  rho: 3.0e-6
  E:   27000.0
  nu:  0.25
  G:   10800.0
  K:   18000.0
  bc:  70440.81901751804
  Gc:  1.0
  εc:  0.0052970143846977355

julia> matofpoint = [1, 1, 2, 2] # points 1:2 -> mat1; points 3:4 -> mat2
4-element Vector{Int64}:
 1
 1
 2
 2

julia> mm = MultiMaterial((mat1, mat2), matofpoint)
MultiMaterial{BBMaterial, 2, UInt8} for 4 points
```
"""
struct MultiMaterial{T <: AbstractMaterial, N, I <: Integer}
    materials::NTuple{N, T}
    matofpoint::Vector{I}

    function MultiMaterial(materials::NTuple{N, T},
                           matofpoint::AbstractVector{I}) where {T <: AbstractMaterial, N,
                                                                 I <: Integer}
        if N == 1
            error("N = 1\nUsage of MultiMaterial only makes sense for N > 1")
        end

        if !all(1 .<= matofpoint .<= N)
            error("Index in `matofpoint` out of bounds!")
        end

        # save memory for small N
        if N <= typemax(UInt8)
            _matofpoint = convert.(UInt8, matofpoint)
        else
            _matofpoint = matofpoint
        end

        return new{T, N, eltype(_matofpoint)}(materials, _matofpoint)
    end
end

function Base.show(io::IO, ::MIME"text/plain", mm::MultiMaterial)
    print(io, typeof(mm), " for ", length(mm.matofpoint), " points")
    return nothing
end

const PDMaterial{T} = Union{T, MultiMaterial{T, N}} where {T <: AbstractMaterial, N}

Base.getindex(mm::MultiMaterial, i::Int) = mm.materials[mm.matofpoint[i]]
function Base.getindex(mm::MultiMaterial, I::AbstractVector{T}) where {T <: Integer}
    return mm.materials[mm.matofpoint[I]]
end
Base.getindex(mm::MultiMaterial, ::Colon) = mm.materials[mm.matofpoint]
Base.firstindex(mm::MultiMaterial) = 1
Base.lastindex(mm::MultiMaterial) = length(mm.matofpoint)

Base.eltype(::MultiMaterial{T}) where {T} = T

Base.getindex(mat::AbstractMaterial, ::Int) = mat
Base.getindex(mat::AbstractMaterial, ::AbstractVector{T}) where {T <: Integer} = mat
Base.getindex(mat::AbstractMaterial, ::Colon) = mat
Base.firstindex(::AbstractMaterial) = 1
Base.lastindex(::AbstractMaterial) = 1

Base.eltype(::T) where {T <: AbstractMaterial} = T
