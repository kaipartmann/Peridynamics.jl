
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
        single_dim_ics = Vector{SingleDimIC}()
        point_sets_precracks = Vector{PointSetsPreCrack}()

        new{M,P}(mat, n_points, position, volume, fail_permit, point_sets, point_params,
                 params_map, single_dim_bcs, single_dim_ics, point_sets_precracks)
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

function _point_set!(point_sets::Dict{Symbol,Vector{Int}}, name::Symbol,
                     points::V) where {V<:AbstractVector{<:Integer}}
    if haskey(point_sets, name)
        error("there is already a point set with name $(name)!")
    end
    point_sets[name] = points
    return nothing
end

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

function velocity_bc!(f::F, b::Body, name::Symbol, d::DimensionSpec) where {F<:Function}
    check_if_set_is_defined(b.point_sets, name)
    check_condition_function(f)
    dim = get_dim(d)
    _condition!(b.single_dim_bcs, SingleDimBC(f, :velocity_half, name, dim))
    return nothing
end

function velocity_ic!(b::Body, name::Symbol, d::DimensionSpec, value::Real)
    check_if_set_is_defined(b.point_sets, name)
    dim = get_dim(d)
    _condition!(b.single_dim_ics, SingleDimIC(value, :velocity, name, dim))
end

function forcedensity_bc!(f::F, b::Body, name::Symbol, d::DimensionSpec) where {F<:Function}
    check_if_set_is_defined(b.point_sets, name)
    check_condition_function(f)
    dim = get_dim(d)
    _condition!(b.single_dim_bcs, SingleDimBC(f, :b_ext, name, dim))
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
