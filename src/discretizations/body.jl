
mutable struct Body
    const n_points::Int
    const position::Matrix{Float64}
    const volume::Vector{Float64}
    const failure_allowed::Vector{Bool}
    psh::AbstractPointSetHandler

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

        psh = UndefPointSetHandler()

        new(n_points, position, volume, failure_allowed, psh)
    end
end

function Body(position::AbstractMatrix, volume::AbstractVector)
    failure_allowed = BitVector(fill(true, length(volume)))
    return Body(position, volume, failure_allowed)
end

function init_no_mat_point_set_handler!(b::Body)
    b.psh = NoMatPointSetHandler()
    return nothing
end

function init_point_set_handler!(::Type{M}, b::Body) where {M<:AbstractMaterial}
    b.psh = PointSetHandler(M)
    return nothing
end

function convert_to_point_set_handler!(::Type{M}, b::Body) where {M<:AbstractMaterial}
    b.psh = PointSetHandler(M, b.psh)
end

function check_if_undef_psh(psh::H, name::Symbol) where {H<:AbstractPointSetHandler}
    if isa(psh, UndefPointSetHandler)
        error("there is no point set with name $(name)! Define the point set first!")
    end
    return nothing
end

function point_set!(b::Body, name::Symbol, points::V) where {V<:AbstractVector}
    isa(b.psh, UndefPointSetHandler) && init_no_mat_point_set_handler!(b)
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
    isa(b.psh, UndefPointSetHandler) && init_point_set_handler!(M, b)
    isa(b.psh, NoMatPointSetHandler) && convert_to_point_set_handler!(M, b)
    if haskey(b.psh.point_sets, :__all__)
        error("body already has a material for all points defined! Choose a point set!")
    end
    _point_set!(b.psh, :__all__, 1:b.n_points)
    _material!(b.psh, :__all__, mat)
    return nothing
end

function material!(b::Body, name::Symbol, mat::M) where {M<:AbstractMaterial}
    check_if_undef_psh(b.psh, name)
    isa(b.psh, NoMatPointSetHandler) && convert_to_point_set_handler!(M, b)
    _material!(b.psh, name, mat)
    return nothing
end

function velocity_bc!(f::F, b::Body, name::Symbol) where {F<:Function}
    # add velocity boundary condition
    check_if_undef_psh(b.psh, name)
end

function velocity_ic!(f::F, b::Body, name::Symbol) where {F<:Function}
    # add velocity initial condition
    check_if_undef_psh(b.psh, name)
end

function forcedensity_bc!(f::F, b::Body, name::Symbol) where {F<:Function}
    # add force density boundary condition
    check_if_undef_psh(b.psh, name)
end

function forcedensity_ic!(f::F, b::Body, name::Symbol) where {F<:Function}
    # add force density initial condition
    check_if_undef_psh(b.psh, name)
end
