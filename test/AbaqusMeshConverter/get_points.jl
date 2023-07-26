using Test
using Peridynamics.AbaqusMeshConverter: get_points

nodes1 = Dict(
    1 => [0.0, 0.0, 0.0],
    2 => [1.0, 0.0, 0.0],
    3 => [0.0, 1.0, 0.0],
    4 => [0.0, 0.0, 1.0],
)
elements1 = Dict(
    1 => [1, 2, 3, 4]
)
element_types1 = Dict(
    1 => :Tet4
)
midpoints1, volumes1 = get_points(nodes1, elements1, element_types1)
@test midpoints1 == [0.25; 0.25; 0.25;;]
@test volumes1 == [1/6]

wrong_elem1 = Dict(
    1 => [1, 2, 3, 4, 5]
)
@test_throws DimensionMismatch get_points(nodes1, wrong_elem1, element_types1)

nodes2 = Dict(
    1 => [0.0, 1.0, 1.0],
    2 => [0.0, 0.0, 1.0],
    3 => [0.0, 0.0, 0.0],
    4 => [0.0, 1.0, 0.0],
    5 => [1.0, 1.0, 1.0],
    6 => [1.0, 0.0, 1.0],
    7 => [1.0, 0.0, 0.0],
    8 => [1.0, 1.0, 0.0],
)
elements2 = Dict(
    1 => [1, 2, 3, 4, 5, 6, 7, 8]
)
element_types2 = Dict(
    1 => :Hex8
)
midpoints2, volumes2 = get_points(nodes2, elements2, element_types2)
@test midpoints2 == [0.5; 0.5; 0.5;;]
@test volumes2 == [1]

wrong_elem2 = Dict(
    1 => [1, 2, 3, 4, 5, 6, 7]
)
@test_throws DimensionMismatch get_points(nodes1, wrong_elem1, element_types1)

wrong_elemtype = Dict(
    1 => :Tet8
)
@test_throws DomainError get_points(nodes1, elements1, wrong_elemtype)
