
struct Body
    n_points::Int
    position::Matrix{Float64}
    volume::Vector{Float64}
    failure_allowed::Vector{Bool}
    point_sets::Dict{Symbol,AbstractVector{<:Integer}}
    materials::Dict{Symbol,AbstractMaterial}
    single_dim_bcs::Vector{SingleDimBC}
    single_dim_ics::Vector{SingleDimIC}
    point_sets_precracks::Vector{PointSetsPreCrack}

    function Body(position::AbstractMatrix, volume::AbstractVector)
        # check if n_points is greater than zero
        n_points = length(volume)
        n_points > 0 || error("number of points `n_points` must be greater than zero!\n")

        # check dimension of position
        dim_position, n_points_position = size(position)
        if dim_position != 3 || n_points_position != n_points
            err_msg = "incorrect dimensions of `position`!\n"
            err_msg *= @sprintf("  should be: (%d, %d)\n", 3, n_points)
            err_msg *= @sprintf("  evaluated: (%d, %d)\n", dim_position, n_points_position)
            throw(DimensionMismatch(err_msg))
        end

        sum(isnan.(position)) > 0 && error("matrix `position` contains NaN values!\n")
        sum(isnan.(volume)) > 0 && error("vector `volume` contains NaN values!\n")

        failure_allowed = BitVector(fill(true, length(volume)))
        point_sets = Dict{Symbol,AbstractVector{<:Integer}}()
        materials = Dict{Symbol,AbstractMaterial}()
        single_dim_bcs = Vector{SingleDimBC}()
        single_dim_ics = Vector{SingleDimIC}()
        point_sets_precracks = Vector{PointSetsPreCrack}()

        new(n_points, position, volume, failure_allowed, point_sets, materials,
            single_dim_bcs, single_dim_ics, point_sets_precracks)
    end
end

function check_if_set_is_defined(point_sets::Dict{Symbol,V}, name::Symbol) where {V}
    if !haskey(point_sets, name)
        error("there is no point set with name $(name)!")
    end
    return nothing
end

function _point_set!(point_sets::Dict{Symbol,V}, name::Symbol,
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

function failure_allowed!(b::Body, name::Symbol, failure_allowed::Bool)
    check_if_set_is_defined(b.point_sets, name)
    b.failure_allowed[b.point_sets[name]] .= failure_allowed
    return nothing
end

function get_material_type(materials::Dict{Symbol,AbstractMaterial})
    if isempty(materials)
        return AbstractMaterial
    end
    mat_types = typeof.(values(materials))
    if !allequal(mat_types)
        error("Materials contains multiple material types!")
    end
    mat_type = first(mat_types)
    return mat_type
end

function _material!(materials::Dict{Symbol,AbstractMaterial}, name::Symbol,
                    mat::M) where {M<:AbstractMaterial}
    mat_type = get_material_type(materials)
    if mat_type != M && isconcretetype(mat_type)
        msg = "body is already assigned with material type $(mat_type))), "
        msg *= "cannot add material with type $(M)!\n"
        throw(ArgumentError(msg))
    end
    if haskey(materials, name)
        msg = "material for set $(name) already defined!\n"
        throw(ArgumentError(msg))
    end
    materials[name] = mat
    return nothing
end

function material!(b::Body, mat::M) where {M<:AbstractMaterial}
    if !haskey(b.point_sets, :__all__)
        _point_set!(b.point_sets, :__all__, 1:b.n_points)
    end
    material!(b, :__all__, mat)
    return nothing
end

function material!(b::Body, name::Symbol, mat::M) where {M<:AbstractMaterial}
    check_if_set_is_defined(b.point_sets, name)
    _material!(b.materials, name, mat)
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

function sets_intersect(set_a::U, set_b::V) where {U,V<:AbstractVector{<:Integer}}
    isempty(set_a ∩ set_b) && return false
    return true
end

function precrack!(b::Body, set_a::Symbol, set_b::Symbol)
    check_if_set_is_defined(b.point_sets, set_a)
    check_if_set_is_defined(b.point_sets, set_b)
    if sets_intersect(b.point_sets[set_a], b.point_sets[set_b])
        msg = "set :$set_a and :$set_b intersect!\n"
        msg *= "No point of the first set is allowed in the second!\n"
        throw(ArgumentError(msg))
    end
    push!(b.point_sets_precracks, PointSetsPreCrack(set_a, set_b))
    return nothing
end
