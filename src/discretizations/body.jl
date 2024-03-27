"""
    Body{M<:AbstractMaterial,P<:AbstractPointParameters}

Body for use in peridynamic calculation

# Fields

- `mat<:AbstractMaterial`: specified material model
- `n_points::Int`: number of material points that represent the body
- `position::Matrix{Float64}`: 3×n_points matrix with position for each point
- `volume::Vector{Float64}`: vector with volume for each point
- `fail_permit::Vector{Bool}`: vector that describes if failure is allowed for each point
- `point_sets::Dict{Symbol,Vector{Int}}`: dictionary containing the defined point sets
- `point_params::Vector{P}`: vector with material parameter sets
- `params_map::Vector{Int}`: vector that assigns a material parameter set to each point
- `single_dim_bcs::Vector{SingleDimBC}`: vector with defined boundary conditions
- `single_dim_ics::Vector{SingleDimIC}`: vector with defined initial conditions
- `point_sets_precracks::Vector{PointSetsPreCrack}`: vector with defined cracks

---

Constructors:

    Body(mat<:AbstractMaterial, position::AbstractMatrix, volume::AbstractVector)

Creates a body for use in peridynamic calculation

# Arguments

- `mat<:AbstractMaterial`: specifies which material model is used
- `position::AbstractMatrix`: 3×n matrix with position of each point
- `volume::AbstractVector`: vector with volume of each point

# Throws

- error if n_points = 0
- `DimensionMismatch`: error if dimension of position != 3

# Example

```julia-repl
julia> l, Δx, = 1.0, 1/50;
julia> pos, vol = uniform_box(l, l, 0.1l, Δx);
julia> b = Body(BBMaterial(), pos, vol);
julia> b
Body{BBMaterial, Peridynamics.BBPointParameters}(BBMaterial(), 12500,
[-0.49 -0.47 … 0.47 0.49; -0.49 -0.49 … 0.49 0.49; -0.04 -0.04 … 0.04 0.04],
[8.000000000000001e-6, 8.000000000000001e-6,  …  8.000000000000001e-6],
Bool[1, 1, 1, 1, 1, 1, 1, 1, 1, 1  …  1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
Dict{Symbol, Vector{Int64}}(), Peridynamics.BBPointParameters[],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0  …  0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
Peridynamics.SingleDimBC[], Peridynamics.SingleDimIC[], Peridynamics.PointSetsPreCrack[])
```
"""
struct Body{M<:AbstractMaterial,P<:AbstractPointParameters}
    mat::M
    n_points::Int
    position::Matrix{Float64}
    volume::Vector{Float64}
    fail_permit::Vector{Bool}
    point_sets::Dict{Symbol,Vector{Int}}
    point_params::Vector{P}
    params_map::Vector{Int}
    single_dim_bcs::Vector{SingleDimBC}
    posdep_single_dim_bcs::Vector{PosDepSingleDimBC}
    single_dim_ics::Vector{SingleDimIC}
    point_sets_precracks::Vector{PointSetsPreCrack}

    function Body(mat::M, position::AbstractMatrix, volume::AbstractVector) where {M}
        n_points = length(volume)
        check_pos_and_vol(n_points, position, volume)
        fail_permit = fill(true, length(volume))
        point_sets = Dict{Symbol,Vector{Int}}()

        P = point_param_type(mat)
        point_params = Vector{P}()
        params_map = zeros(Int, n_points)

        single_dim_bcs = Vector{SingleDimBC}()
        posdep_single_dim_bcs = Vector{PosDepSingleDimBC}()
        single_dim_ics = Vector{SingleDimIC}()
        point_sets_precracks = Vector{PointSetsPreCrack}()

        new{M,P}(mat, n_points, position, volume, fail_permit, point_sets, point_params,
                 params_map, single_dim_bcs, posdep_single_dim_bcs, single_dim_ics,
                 point_sets_precracks)
    end
end

@inline material_type(::Body{M,P}) where {M,P} = M

function check_pos_and_vol(n_points::Int, position::AbstractMatrix, volume::AbstractVector)
    # check if n_points is greater than zero
    n_points > 0 || error("number of points `n_points` must be greater than zero!\n")

    # check dimension of position
    dim_position, n_points_position = size(position)
    if dim_position != 3 || n_points_position != n_points
        err_msg = "incorrect dimensions of `position`!\n"
        err_msg *= @sprintf("  should be: (%d, %d)\n", 3, n_points)
        err_msg *= @sprintf("  evaluated: (%d, %d)\n", dim_position, n_points_position)
        throw(DimensionMismatch(err_msg))
    end

    # check if they contain NaN's
    sum(isnan.(position)) > 0 && error("matrix `position` contains NaN values!\n")
    sum(isnan.(volume)) > 0 && error("vector `volume` contains NaN values!\n")

    return nothing
end

function check_if_set_is_defined(point_sets::Dict{Symbol,V}, name::Symbol) where {V}
    if !haskey(point_sets, name)
        error("there is no point set with name $(name)!")
    end
    return nothing
end

"""
    point_set!(b::Body, name::Symbol, points::V) where {V<:AbstractVector}
    point_set!(f::F, b::Body, name::Symbol) where {F<:Function}

Creates the point set `name` in Body `b` containing all points either defined in `V` or
described by `F`

# Arguments

- `b::Body`: peridynamic body
- `name::Symbol`: name of the point set
- `points<:AbstractVector`: vector of point indices
- `f<:Function`: function that describes points contained in point set

# Throws

- error if a point set called `name` is already defined
- `BoundsError`: if points in `V` do not exist

# Example

```julia-repl
julia> point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, b, :set_a)
julia> point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, b, :set_b)
julia> b.point_sets
Dict{Symbol, Vector{Int64}} with 4 entries:
  :set_a      => [1251, 1252, 1253, 1254, 1255, 1256, 1257, 1258, 1259, …
  :set_b      => [951, 952, 953, 954, 955, 956, 957, 958, 959, 960  …  1…
```
"""
function point_set! end

function point_set!(b::Body, name::Symbol, points::V) where {V<:AbstractVector}
    checkbounds(b.volume, points)
    _point_set!(b.point_sets, name, points)
    return nothing
end

function point_set!(f::F, b::Body, name::Symbol) where {F<:Function}
    points = find_points(f, b.position)
    point_set!(b, name, points)
    return nothing
end

function _point_set!(point_sets::Dict{Symbol,Vector{Int}}, name::Symbol,
                     points::V) where {V<:AbstractVector{<:Integer}}
    if haskey(point_sets, name)
        error("there is already a point set with name $(name)!")
    end
    point_sets[name] = points
    return nothing
end

"""
    failure_permit!(b::Body, fail_permit::Bool)
    failure_permit!(b::Body, name::Symbol, fail_permit::Bool)

determines whether failure is permitted `fail_permit = true` or prohibited
`fail_permit = false` in the body `b` or the point set `name` of body `b`

# Arguments

- `b::Body`: peridynamic body
- `name::Symbol`: name of the point set on body `b`
- `fail_permit::Bool`: decides if failure is allowed on considered body or point set

# Throws

- error if no point set called `name` exists

# Example

```julia-repl
julia> failure_permit!(b, :set_bottom, false)
julia> b.fail_permit
12500-element Vector{Bool}:
 0
 0
 0
 ⋮
 1
 1
```
"""
function failure_permit! end

function failure_permit!(b::Body, fail_permit::Bool)
    b.fail_permit .= fail_permit
    return nothing
end

function failure_permit!(b::Body, name::Symbol, fail_permit::Bool)
    check_if_set_is_defined(b.point_sets, name)
    b.fail_permit[b.point_sets[name]] .= fail_permit
    return nothing
end

function check_material_kwargs(mat::M, p::Dict{Symbol,Any}) where {M<:AbstractMaterial}
    allowed_kwargs = allowed_material_kwargs(mat)
    check_kwargs(p, allowed_kwargs)
    return nothing
end

"""
    material!(b::Body{M,P}, name::Symbol; kwargs...)
    material!(b::Body{M,P}; kwargs...)
        where {M<:AbstractMaterial,P<:AbstractPointParameters}

specifies material parameters used for body `b` or point set `name`

# Arguments

- `b::Body{M<:AbstractMaterial,P<:AbstractPointParameters}`: peridynamic body
- `name::Symbol`: point set on body `b`

# Keywords

Allowed keywords depend on selected material model. See material type documentation.
Default material parameter keywords:

- `:horizon::Float64`: radius of point interactions
- `:rho::Float64`: density
- `:E::Float64`: Young's modulus
- `:nu::Float64`: Poisson's ratio
- `:Gc::Float64`: critical energy release rate
- `:epsilon_c::Float64`: critical strain

# Throws

- error if parameter is not eligible for specification in selected material model

# Example

```julia-repl
julia> material!(b; horizon=δ, E=2.1e5, rho=8e-6, Gc=2.7)
julia> b.point_params
1-element Vector{Peridynamics.BBPointParameters}:
 Peridynamics.BBPointParameters(0.0603, 8.0e-6, 210000.0, 0.25, 84000.0, 140000.0, 84000.0,
84000.0, 2.7, 0.013329779199368195, 6.0671037207022026e10)
```
"""
function material! end

function material!(b::Body{M,P}, name::Symbol; kwargs...) where {M,P}
    check_if_set_is_defined(b.point_sets, name)

    p = Dict{Symbol,Any}(kwargs)
    check_material_kwargs(b.mat, p)

    points = b.point_sets[name]
    params = get_point_params(b.mat, p)

    _material!(b, points, params)

    return nothing
end

function material!(b::Body{M,P}; kwargs...) where {M,P}
    p = Dict{Symbol,Any}(kwargs)
    check_material_kwargs(b.mat, p)

    isempty(b.point_params) || empty!(b.point_params)

    points = 1:b.n_points
    params = get_point_params(b.mat, p)

    _material!(b, points, params)
    return nothing
end

function _material!(b::Body{M,P}, points::V, params::P) where {M,P,V}
    push!(b.point_params, params)
    id = length(b.point_params)
    b.params_map[points] .= id
    return nothing
end

function get_dim(dim::I) where {I<:Integer}
    in(dim, 1:3) || error("specified dimension should be 1=x, 2=y, or 3=z!\n")
    isa(UInt8, I) && return dim
    return convert(UInt8, dim)
end

function get_dim(dim::Symbol)
    if !haskey(SYMBOL_TO_DIM, dim)
        error("unknown dimension symbol $(dim)! Should be :x, :y, or :z!\n")
    end
    return SYMBOL_TO_DIM[dim]
end

function _condition!(conditions::Vector{B}, condition::B) where {B<:AbstractCondition}
    # check if conditions override each other!
    if is_duplicate(condition, conditions)
        error("duplicate conditions for point set $(condition.point_set)!\n")
    end
    push!(conditions, condition)
    return nothing
end

"""
    velocity_bc!(f::F, b::Body, name::Symbol, d::DimensionSpec) where {F<:Function}

specifies velocity boundary conditions for point set `name` on body `b`

# Arguments

- `f<:Function`: velocity condition function
- `b::Body`: peridynamic body
- `name::Symbol`: point set on body `b`
- `d::DimensionSpec`: direction of velocity

# Throws

- error if no point set called `name` exists
- error if dimension is not correctly specified
- error if function is not suitable as condition function

# Example

```julia-repl
julia> velocity_bc!(t -> -9.81 * t, b, :set_bottom, :y)
julia> b.single_dim_bcs
1-element Vector{Peridynamics.SingleDimBC}:
 Peridynamics.SingleDimBC{var"#15#16"}(var"#15#16"(), :velocity_half, :set_bottom, 0x02)
julia> velocity_bc!(t -> 40, b, :set_a, 1)
julia> b.single_dim_bcs
2-element Vector{Peridynamics.SingleDimBC}:
 Peridynamics.SingleDimBC{var"#15#16"}(var"#15#16"(), :velocity_half, :set_bottom, 0x02)
 Peridynamics.SingleDimBC{var"#17#18"}(var"#17#18"(), :velocity_half, :set_a, 0x01)
```
"""
function velocity_bc!(f::F, b::Body, name::Symbol, d::DimensionSpec) where {F<:Function}
    check_if_set_is_defined(b.point_sets, name)
    type = check_condition_function(f)
    dim = get_dim(d)
    if type === :sdbc
        sdbc = SingleDimBC(f, :velocity_half, name, dim)
        _condition!(b.single_dim_bcs, sdbc)
    elseif type === :pdsdbc
        pdsdbc = PosDepSingleDimBC(f, :velocity_half, name, dim)
        _condition!(b.posdep_single_dim_bcs, pdsdbc)
    end
    return nothing
end

"""
    velocity_ic!(b::Body, name::Symbol, d::DimensionSpec, value::Real)

specifies initital conditions for the velocity of points in point set `name` on body `b`

# Arguments

- `b::Body`: peridynamic body
- `name::Symbol`: point set on body `b`
- `d::DimensionSpec`: direction of velocity
- `value::Real`: initial velocity value

# Throws

- error if no point set called `name` exists
- error if dimension is not correctly specified

# Example

```julia-repl
julia> velocity_ic!(b, :set_b, :y, 20)
julia> b.single_dim_ics
1-element Vector{Peridynamics.SingleDimIC}:
 Peridynamics.SingleDimIC(20.0, :velocity, :set_b, 0x02)
```
"""
function velocity_ic!(b::Body, name::Symbol, d::DimensionSpec, value::Real)
    check_if_set_is_defined(b.point_sets, name)
    dim = get_dim(d)
    sdic = SingleDimIC(convert(Float64, value), :velocity, name, dim)
    _condition!(b.single_dim_ics, sdic)
    return nothing
end

"""
    forcedensity_bc!(f::F, b::Body, name::Symbol, d::DimensionSpec) where {F<:Function}

specifies boundary conditions for force density on points of point set `name` on body `b`

# Arguments

- `f<:Function`: condition function
- `b::Body`: peridynamic body
- `name::Symbol`: point set on body `b`
- `d::DimensionSpec`: direction of force density

# Throws

- error if no point set called `name` exists
- error if dimension is not correctly specified
- error if function is not suitable as condition function

# Example

```julia-repl
julia> forcedensity_bc!(t -> 40, b, :set_a, 1)
julia> b.single_dim_bcs
1-element Vector{Peridynamics.SingleDimBC}:
 Peridynamics.SingleDimBC{var"#25#26"}(var"#25#26"(), :b_ext, :set_a, 0x01)
```
"""
function forcedensity_bc!(f::F, b::Body, name::Symbol, d::DimensionSpec) where {F<:Function}
    check_if_set_is_defined(b.point_sets, name)
    type = check_condition_function(f)
    dim = get_dim(d)
    if type === :sdbc
        sdbc = SingleDimBC(f, :b_ext, name, dim)
        _condition!(b.single_dim_bcs, sdbc)
    elseif type === :pdsdbc
        pdsdbc = PosDepSingleDimBC(f, :b_ext, name, dim)
        _condition!(b.posdep_single_dim_bcs, pdsdbc)
    end
    return nothing
end

function check_if_sets_intersect(point_sets::Dict{Symbol,Vector{Int}}, key_a::Symbol,
                                 key_b::Symbol)
    set_a, set_b = point_sets[key_a], point_sets[key_b]
    if !isempty(set_a ∩ set_b)
        msg = "set :$key_a and :$key_b intersect!\n"
        msg *= "No point of the first set is allowed in the second!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

"""
    precrack!(b::Body, set_a::Symbol, set_b::Symbol)

creates a crack between two point sets by prohibiting interaction between points of
different point sets

# Arguments

- `b::Body`: peridynamic body
- `set_a::Symbol`: first point set
- `set_b::Symbol`: second point set

# Throws

- error if point set `set_a` or `set_b` does not exist
- error if point sets contain common points

# Example

```julia-repl
julia> precrack!(b, :set_a, :set_b)
julia> b.point_sets_precracks
1-element Vector{Peridynamics.PointSetsPreCrack}:
 Peridynamics.PointSetsPreCrack(:set_a, :set_b)
```
"""
function precrack!(b::Body, set_a::Symbol, set_b::Symbol)
    check_if_set_is_defined(b.point_sets, set_a)
    check_if_set_is_defined(b.point_sets, set_b)
    check_if_sets_intersect(b.point_sets, set_a, set_b)
    push!(b.point_sets_precracks, PointSetsPreCrack(set_a, set_b))
    return nothing
end

function pre_submission_check(b::Body)
    #TODO: check if everything is defined for job submission!
    return nothing
end

@inline function get_point_param(b::Body, key::Symbol, i::Int)
    return getfield(b.point_params[b.params_map[i]], key)
end

@inline storage_type(b::Body, ts::AbstractTimeSolver) = storage_type(b.mat, ts)
