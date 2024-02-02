
mutable struct Body
    const n_points::Int
    const position::Matrix{Float64}
    const volume::Vector{Float64}
    const failure_allowed::Vector{Bool}
    const single_dim_bcs::Vector{SingleDimBC}
    const single_dim_ics::Vector{SingleDimIC}
    const point_sets_precracks::Vector{PointSetsPreCrack}
    psh::PointSetHandler{<:AbstractMaterial}

    function Body(position::AbstractMatrix, volume::AbstractVector,
                  failure_allowed::BitVector)
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

        # check dimension of failure_allowed
        n_points_failure_allowed = length(failure_allowed)
        if n_points_failure_allowed != n_points
            err_msg = "Incorrect length of `failure_allowed`!\n"
            err_msg *= @sprintf("  should be: %d\n", n_points)
            err_msg *= @sprintf("  evaluated: %d\n", n_points_failure_allowed)
            throw(DimensionMismatch(err_msg))
        end

        sum(isnan.(position)) > 0 && error("matrix `position` contains NaN values!\n")
        sum(isnan.(volume)) > 0 && error("vector `volume` contains NaN values!\n")

        single_dim_bcs = Vector{SingleDimBC}()
        single_dim_ics = Vector{SingleDimIC}()
        point_sets_precracks = Vector{PointSetsPreCrack}()
        psh = PointSetHandler(AbstractMaterial)

        new(n_points, position, volume, failure_allowed, single_dim_bcs, single_dim_ics,
            point_sets_precracks, psh)
    end
end

function Body(position::AbstractMatrix, volume::AbstractVector)
    failure_allowed = BitVector(fill(true, length(volume)))
    return Body(position, volume, failure_allowed)
end

function point_set!(b::Body, name::Symbol, points::V) where {V<:AbstractVector}
    checkbounds(b.volume, points)
    _point_set!(b.psh, name, points)
    return nothing
end

function point_set!(f::F, b::Body, name::Symbol) where {F<:Function}
    points = find_points(f, b.position)
    point_set!(b, name, points)
    return nothing
end

function material!(b::Body, mat::M) where {M<:AbstractMaterial}
    if !haskey(b.psh.point_sets, :__all__)
        _point_set!(b.psh, :__all__, 1:b.n_points)
    end
    material!(b, :__all__, mat)
    return nothing
end

function material!(b::Body, name::Symbol, mat::M) where {M<:AbstractMaterial}
    check_if_set_is_defined(b.psh, name)
    if material_is_undefined(b.psh)
        b.psh = PointSetHandler(M, b.psh)
    end
    _material!(b.psh, name, mat)
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
    check_if_set_is_defined(b.psh, name)
    check_condition_function(f)
    dim = get_dim(d)
    _condition!(b.single_dim_bcs, SingleDimBC(f, :velocity_half, name, dim))
    return nothing
end

function velocity_ic!(b::Body, name::Symbol, d::DimensionSpec, value::Real)
    check_if_set_is_defined(b.psh, name)
    dim = get_dim(d)
    _condition!(b.single_dim_ics, SingleDimIC(value, :velocity, name, dim))
end

function forcedensity_bc!(f::F, b::Body, name::Symbol, d::DimensionSpec) where {F<:Function}
    check_if_set_is_defined(b.psh, name)
    check_condition_function(f)
    dim = get_dim(d)
    _condition!(b.single_dim_bcs, SingleDimBC(f, :b_ext, name, dim))
end

function sets_intersect(set_a::U, set_b::V) where {U,V<:AbstractVector{<:Integer}}
    isempty(set_a ∩ set_b) && return false
    return true
end

function precrack!(b::Body, set_a::Symbol, set_b::Symbol)
    check_if_set_is_defined(b.psh, set_a)
    check_if_set_is_defined(b.psh, set_b)
    if sets_intersect(b.psh.point_sets[set_a], b.psh.point_sets[set_b])
        msg = "set :$set_a and :$set_b intersect!\n"
        msg *= "No point of the first set is allowed in the second!\n"
        throw(ArgumentError(msg))
    end
    push!(b.point_sets_precracks, PointSetsPreCrack(set_a, set_b))
    return nothing
end
