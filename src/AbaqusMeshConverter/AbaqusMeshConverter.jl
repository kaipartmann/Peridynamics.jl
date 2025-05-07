module AbaqusMeshConverter

using StaticArrays
using LinearAlgebra
using AbaqusReader: abaqus_read_mesh

const SUPPORTED_ELEMENT_TYPES = [:Tet4, :Hex8]

export read_inp

function get_points(nodes::Dict{Int, Vector{Float64}}, elements::Dict{Int, Vector{Int}},
                    element_types::Dict{Int, Symbol})
    nel = length(elements)
    midpoints = zeros(3, nel)
    volumes = zeros(nel)
    for i in 1:nel
        T = element_types[i]
        if T in SUPPORTED_ELEMENT_TYPES
            if T == :Tet4
                nodeids = elements[i]
                if length(nodeids) == 4
                    a = nodes[nodeids[1]]
                    b = nodes[nodeids[2]]
                    c = nodes[nodeids[3]]
                    d = nodes[nodeids[4]]
                    @views midpoints[:, i] = midpoint(a, b, c, d)
                    volumes[i] = tetvol(a, b, c, d)
                else
                    throw(DimensionMismatch("4 nodes needed for element of type Tet4!"))
                end
            end
            if T == :Hex8
                nodeids = elements[i]
                if length(nodeids) == 8
                    n1 = nodes[nodeids[1]]
                    n2 = nodes[nodeids[2]]
                    n3 = nodes[nodeids[3]]
                    n4 = nodes[nodeids[4]]
                    n5 = nodes[nodeids[5]]
                    n6 = nodes[nodeids[6]]
                    n7 = nodes[nodeids[7]]
                    n8 = nodes[nodeids[8]]
                    @views midpoints[:, i] = midpoint(n1, n2, n3, n4, n5, n6, n7, n8)
                    volume1 = tetvol(n1, n2, n4, n5)
                    volume2 = tetvol(n2, n3, n4, n7)
                    volume3 = tetvol(n2, n5, n6, n7)
                    volume4 = tetvol(n4, n5, n7, n8)
                    volume5 = tetvol(n2, n4, n5, n7)
                    volumes[i] = volume1 + volume2 + volume3 + volume4 + volume5
                else
                    throw(DimensionMismatch("8 nodes needed for element of type Hex8"))
                end
            end
        else
            msg = "Element of type $T not supported!\n"
            msg *= "Supported types: $(SUPPORTED_ELEMENT_TYPES)"
            throw(DomainError(T, msg))
        end
    end
    return midpoints, volumes
end

"""
    read_inp(file::String)

Read Abaqus .inp-file and convert meshes to a point cloud with the help of the
[`AbaqusReader.jl`](https://github.com/JuliaFEM/AbaqusReader.jl) package.
Every element is converted to a point. The center of the element becomes the position of the
point and the element volume becomes the point volume.
Element sets defined in Abaqus are converted to corresponding point sets.

Currently supported mesh elements: $SUPPORTED_ELEMENT_TYPES

# Arguments
- `file::String`: Path to Abaqus .inp-file.

# Returns
- `position::Matrix{Float64}`: Point position (midpoint of every element).
- `volume::Vector{Float64}`: Point volume (volume of every element).
- `point_sets`: Element sets defined in the .inp-file.
"""
function read_inp(file::String)
    if !endswith(file, ".inp")
        _, extension = splitext(basename(file))
        throw(ArgumentError("cannot read file with extension $(extension)!\n"))
    end
    am = abaqus_read_mesh(file)
    position, volume = get_points(am["nodes"], am["elements"], am["element_types"])
    point_sets = am["element_sets"]
    return position, volume, point_sets
end

function midpoint(a::T, b::T, c::T, d::T) where {T <: AbstractVector}
    mp = 1 / 4 .* (a + b + c + d)
    return mp
end

function midpoint(a::T, b::T, c::T, d::T, e::T, f::T, g::T,
                  h::T) where {T <: AbstractVector}
    mp = 1 / 8 .* (a + b + c + d + e + f + g + h)
    return mp
end

function tetvol(a::T, b::T, c::T, d::T) where {T <: AbstractVector}
    v1 = SVector{3}(b - a)
    v2 = SVector{3}(c - a)
    v3 = SVector{3}(d - a)
    vol = 1 / 6 * abs(dot(v1, cross(v2, v3)))
    return vol
end

end # module AbaqusMeshConverter
