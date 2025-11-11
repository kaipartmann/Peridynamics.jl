"""
    Body(material, position, volume)
    Body(material, inp_file)

Construct a `Body` for a peridynamics simulation.

# Arguments
- `material::AbstractMaterial`: The material which is defined for the whole body.
    Available material models:
    - [`BBMaterial`](@ref): Bond-based peridynamics.
    - [`DHBBMaterial`](@ref): Dual-horizon bond-based peridynamics.
    - [`GBBMaterial`](@ref): Generalized bond-based peridynamics.
    - [`OSBMaterial`](@ref): Ordinary state-based peridynamics, also called linear
        peridynamic solid (LPS).
    - [`CMaterial`](@ref): Correspondence formulation.
    - [`CRMaterial`](@ref): Correspondence formulation with stress rotation for objectivity
        enforcement.
    - [`RKCMaterial`](@ref): Reproducing kernel peridynamics with bond-associated
        higher-order integration.
    - [`RKCRMaterial`](@ref): Reproducing kernel peridynamics with bond-associated
        higher-order integration with stress rotation for objectivity enforcement.
    - [`BACMaterial`](@ref): Bond-associated correspondence formulation of Chen and Spencer.
    - [`CKIMaterial`](@ref): Continuum-kinematics-inspired peridynamics.
- `position::AbstractMatrix`: A `3×n` matrix with the point position of the `n` points.
- `volume::AbstractVector`: A vector with the volume of each point.
- `inp_file::AbstractString`: An Abaqus input file containing meshes, imported with
    [`read_inp`](@ref).

# Throws
- Error if the number of points is not larger than zero.
- Error if `position` is not a `3×n` matrix and has the same length as `volume`.
- Error if `position` or `volume` contain `NaN` values.

# Example

```julia-repl
julia> Body(BBMaterial(), rand(3, 10), rand(10))
10-point Body{BBMaterial{NoCorrection}}:
  1 point set(s):
    10-point set `all_points`
```

---

!!! warning "Internal use only"
    Please note that the fields are intended for internal use only. They are *not* part of
    the public API of Peridynamics.jl, and thus can be altered (or removed) at any time
    without it being considered a breaking change.

```julia
Body{Material,PointParameters}
```

# Type Parameters

- `Material <: AbstractMaterial`: Type of the specified material model.
- `PointParameters <: AbstractPointParameters`: Type of the point parameters.

# Fields

- `mat::Material`: The material formulation.
- `n_points::Int`: The number of points that in the body.
- `position::Matrix{Float64}`: A `3×n_points` matrix with the position of the points.
- `volume::Vector{Float64}`: A vector with the volume of each point.
- `fail_permit::Vector{Bool}`: A vector that describes if failure is allowed for each point.
- `point_sets::Dict{Symbol,Vector{Int}}`: A dictionary containing point sets.
- `point_params::Vector{PointParameters}`: A vector containing all different point parameter
    instances of the body. Each point can have its own `PointParameters` instance.
- `params_map::Vector{Int}`: A vector that maps each point index to a parameter instance in
    `point_params`.
- `single_dim_bcs::Vector{SingleDimBC}`: A vector with boundary conditions on a single
    dimension.
- `posdep_single_dim_bcs::Vector{PosDepSingleDimBC}`: A vector with position dependent
    boundary conditions on a single dimension.
- `single_dim_ics::Vector{SingleDimIC}`: A vector with initial conditions on a single
    dimension.
- `posdep_single_dim_ics::Vector{PosDepSingleDimIC}`: A vector with position dependent
    initial conditions on a single dimension.
- `data_bcs::Vector{DataBC}`: A vector with data boundary conditions.
- `point_sets_precracks::Vector{PointSetsPreCrack}`: A vector with predefined point set
    cracks.
"""
struct Body{M<:AbstractMaterial,P<:AbstractPointParameters} <: AbstractBody{M}
    name::Ref{Symbol}
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
    pos_single_dim_bcs::Vector{PosSingleDimBC}
    data_bcs::Vector{DataBC}
    single_dim_ics::Vector{SingleDimIC}
    posdep_single_dim_ics::Vector{PosDepSingleDimIC}
    point_sets_precracks::Vector{PointSetsPreCrack}

    function Body(mat::M, position::AbstractMatrix, volume::AbstractVector) where {M}
        name = Symbol("")
        n_points = length(volume)
        check_pos_and_vol(n_points, position, volume)
        fail_permit = fill(true, length(volume))
        point_sets = Dict{Symbol,Vector{Int}}(:all_points => 1:length(volume))

        P = point_param_type(mat)
        point_params = Vector{P}()
        params_map = zeros(Int, n_points)

        single_dim_bcs = Vector{SingleDimBC}()
        posdep_single_dim_bcs = Vector{PosDepSingleDimBC}()
        pos_single_dim_bcs = Vector{PosSingleDimBC}()
        single_dim_ics = Vector{SingleDimIC}()
        posdep_single_dim_ics = Vector{PosDepSingleDimIC}()
        data_bcs = Vector{DataBC}()
        point_sets_precracks = Vector{PointSetsPreCrack}()

        new{M,P}(name, mat, n_points, position, volume, fail_permit, point_sets,
                 point_params, params_map, single_dim_bcs, posdep_single_dim_bcs,
                 pos_single_dim_bcs, data_bcs, single_dim_ics, posdep_single_dim_ics,
                 point_sets_precracks)
    end
end

function Base.show(io::IO, @nospecialize(body::AbstractBody))
    print(io, body.n_points, "-point Body{", material_type(body), "}")
    if has_name(body)
        print(io, " with name `", get_name(body), "`")
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(body::AbstractBody))
    if get(io, :compact, false)
        show(io, body)
        return nothing
    end
    print(io, body.n_points, "-point Body{", material_type(body), "}")
    if has_name(body)
        print(io, " with name `", get_name(body), "`")
    end
    if has_point_sets(body) || has_params(body) || has_conditions(body)
        print(io, ":")
    end
    if has_point_sets(body)
        print(io, "\n  ", length(keys(body.point_sets)), " point set(s):")
        for (name, points) in body.point_sets
            print(io, "\n    ", length(points), "-point set `", name, "`")
        end
    end
    if has_params(body)
        print(io, "\n  ", n_params(body), " point parameter(s):")
        for params in body.point_params
            print(io, "\n    ")
            show(io, params)
        end
    end
    if has_bcs(body)
        print(io, "\n  ", n_bcs(body), " boundary condition(s):")
        for bc in body.single_dim_bcs
            print(io, "\n    ")
            show(io, bc)
        end
        for bc in body.posdep_single_dim_bcs
            print(io, "\n    ")
            show(io, bc)
        end
        for bc in body.pos_single_dim_bcs
            print(io, "\n    ")
            show(io, bc)
        end
        for bc in body.data_bcs
            print(io, "\n    ")
            show(io, bc)
        end
    end
    if has_ics(body)
        print(io, "\n  ", n_ics(body), " initial condition(s):")
        for bc in body.single_dim_ics
            print(io, "\n    ")
            show(io, bc)
        end
        for bc in body.posdep_single_dim_ics
            print(io, "\n    ")
            show(io, bc)
        end
    end
    if has_precracks(body)
        print(io, "\n  ", n_precracks(body), " predefined crack(s)")
    end
    n_failure_permit_false = body.n_points - sum(body.fail_permit)
    if n_failure_permit_false > 0
        print(io, "\n  ", n_failure_permit_false,
              " points with failure prohibited")
    end
    return nothing
end

function Body(mat::AbstractMaterial, inp_file::AbstractString)
    position, volume, point_sets = read_inp(inp_file)
    body = Body(mat, position, volume)
    for (name, point_ids) in point_sets
        point_set!(body, Symbol(clean_point_set_name(name)), point_ids)
    end
    return body
end

@inline material_type(::AbstractBody{M}) where {M} = M.name.wrapper

"""
    check_pos_and_vol(n_points, position, volume)

$(internal_api_warning())

Check if the positions and volumes for the points are correctly specified in the fields of
a [`Body`](@ref).
"""
function check_pos_and_vol(n_points::Int, position::AbstractMatrix, volume::AbstractVector)
    # check if n_points is greater than zero
    n_points > 0 || error("the number of points must be greater than zero!\n")

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

"""
    pre_submission_check(body::Body; body_in_multibody_setup::Bool=false)
    pre_submission_check(ms::AbstractMultibodySetup)

$(internal_api_warning())

Check if necessary material parameters and conditions are defined when defining a
[`Job`](@ref).
"""
function pre_submission_check(body::Body; body_in_multibody_setup::Bool=false)
    # the body should have material properties
    if isempty(body.point_params)
        msg = "no material parameters found!\n"
        msg *= "Bodies without material parameters are not ready for job submission!\n"
        error(msg)
    end
    # all points should have a material property defined
    points_without_material = findfirst(x -> x == 0, body.params_map)
    if points_without_material !== nothing
        msg = "not all points have material parameters!\n"
        msg *= "You probably just used the `material!(body, set; kwargs...)` function\n"
        msg *= "on a set that does not include all material points.\n"
        error(msg)
    end
    # if the body is not part of a multibody setup, then there has to be any condition
    if !body_in_multibody_setup
        if n_bcs(body) + n_ics(body) == 0
            msg = "no initial or boundary condition specified!\n"
            msg *= "Bodies that are not part of a `MultibodySetup` need at least one\n"
            msg *= "initial or boundary condition! Otherwise nothing will happen during\n"
            msg *= "the simulation!\n"
            error(msg)
        end
    end
    return nothing
end

@inline function get_point_param(b::AbstractBody, i::Int)
    return b.point_params[b.params_map[i]]
end

@inline function get_point_param(b::AbstractBody, key::Symbol, i::Int)
    return getfield(b.point_params[b.params_map[i]], key)
end

@inline storage_type(b::AbstractBody) = storage_type(b.mat)

function log_msg_body(body::AbstractBody)
    msg = "BODY"
    body_name = string(get_name(body))
    isempty(body_name) || (msg *= " `" * body_name * "`")
    msg *= "\n"
    msg *= "  POINT CLOUD\n"
    msg *= msg_qty("number of points", body.n_points; indentation=4)
    @views min_x, max_x = minimum(body.position[1, :]), maximum(body.position[1, :])
    @views min_y, max_y = minimum(body.position[2, :]), maximum(body.position[2, :])
    @views min_z, max_z = minimum(body.position[3, :]), maximum(body.position[3, :])
    minmax_x = @sprintf("%.7g, %.7g", min_x, max_x)
    minmax_y = @sprintf("%.7g, %.7g", min_y, max_y)
    minmax_z = @sprintf("%.7g, %.7g", min_z, max_z)
    msg *= msg_qty("min, max values x-direction", minmax_x; indentation=4)
    msg *= msg_qty("min, max values y-direction", minmax_y; indentation=4)
    msg *= msg_qty("min, max values z-direction", minmax_z; indentation=4)
    msg *= "  POINT SETS\n"
    for (key, points) in body.point_sets
        descr = @sprintf("number of points in set `%s`", string(key))
        msg *= msg_qty(descr, length(points); indentation=4)
    end
    has_ics(body) && (msg *= "  INITIAL CONDITIONS\n")
    for ic in body.single_dim_ics
        descr = @sprintf("%s condition", field_to_name(ic.field))
        settings = @sprintf("set `%s`, dimension %d", ic.point_set, ic.dim)
        msg *= msg_qty(descr, settings; indentation=4)
    end
    for ic in body.posdep_single_dim_ics
        descr = @sprintf("%s condition", field_to_name(ic.field))
        settings = @sprintf("set `%s`, dimension %d", ic.point_set, ic.dim)
        msg *= msg_qty(descr, settings; indentation=4)
    end
    has_bcs(body) && (msg *= "  BOUNDARY CONDITIONS\n")
    for bc in body.single_dim_bcs
        descr = @sprintf("%s condition", field_to_name(bc.field))
        settings = @sprintf("set `%s`, dimension %d", bc.point_set, bc.dim)
        msg *= msg_qty(descr, settings; indentation=4)
    end
    for bc in body.posdep_single_dim_bcs
        descr = @sprintf("%s condition", field_to_name(bc.field))
        settings = @sprintf("set `%s`, dimension %d", bc.point_set, bc.dim)
        msg *= msg_qty(descr, settings; indentation=4)
    end
    for bc in body.pos_single_dim_bcs
        descr = @sprintf("%s condition", field_to_name(bc.field))
        settings = @sprintf("set `%s`, dimension %d", bc.point_set, bc.dim)
        msg *= msg_qty(descr, settings; indentation=4)
    end
    msg *= "  MATERIAL\n"
    msg *= log_material(body.mat; indentation=4)
    n_point_params = length(body.point_params)
    if n_point_params == 1
        msg *= log_material_parameters(first(body.point_params); indentation=4)
    elseif n_point_params > 1
        for (i, params) in enumerate(body.point_params)
            msg *= @sprintf("    MATERIAL PROPERTIES #%d\n", i)
            msg *= log_material_parameters(params; indentation=6)
        end
    end
    return msg
end

function log_spatial_setup(options::AbstractJobOptions, body::AbstractBody)
    msg = log_msg_body(body)
    log_it(options, msg)
    return nothing
end

function maximum_horizon(b::AbstractBody)
    n_params = length(b.point_params)
    n_params == 1 && return get_horizon(first(b.point_params))
    n_params == 0 && error("body has no material parameters!\n")
    δmax = 0.0
    for point_id in eachindex(b.volume)
        δ::Float64 = get_horizon(b.point_params[b.params_map[point_id]])
        if δ > δmax
            δmax = δ
        end
    end
    return δmax
end

@inline get_horizon(param::AbstractPointParameters) = getfield(param, :δ)

@inline get_name(body::AbstractBody) = body.name[]

@inline function change_name!(body::AbstractBody, name::Symbol)
    body.name[] = name
    return nothing
end

@inline has_name(body::AbstractBody) = body.name[] !== Symbol("")

@inline has_point_sets(body::AbstractBody) = !isempty(body.point_sets)

@inline n_params(body::AbstractBody) = length(body.point_params)
@inline has_params(body::AbstractBody) = !isempty(body.point_params)

@inline function n_bcs(body::AbstractBody)
    n_sdbcs = length(body.single_dim_bcs)
    n_pdsdbcs = length(body.posdep_single_dim_bcs)
    n_psdbcs = length(body.pos_single_dim_bcs)
    n_databcs = length(body.data_bcs)
    return n_sdbcs + n_pdsdbcs + n_psdbcs + n_databcs
end
@inline function n_ics(body::AbstractBody)
    n_sdics = length(body.single_dim_ics)
    n_pdsdics = length(body.posdep_single_dim_ics)
    return n_sdics + n_pdsdics
end

@inline has_bcs(body::AbstractBody) = n_bcs(body) > 0 ? true : false
@inline has_ics(body::AbstractBody) = n_ics(body) > 0 ? true : false
@inline has_conditions(body::AbstractBody) = has_bcs(body) || has_ics(body)

@inline n_precracks(body::AbstractBody) = length(body.point_sets_precracks)
@inline has_precracks(body::AbstractBody) = n_precracks(body) > 0 ? true : false

"""
    n_points(body)

Return the total number of points in a body.

# Arguments
- `body::Body`: [`Body`](@ref).

# Returns
- `n_points::Int`: The number of points in the body.

# Examples
```julia-repl
julia> body = Body(BBMaterial(), pos, vol)
1000-point Body{BBMaterial{NoCorrection}}:
  1 point set(s):
    1000-point set `all_points`

julia> n_points(body)
1000
```

---

    n_points(multibody_setup)

Return the total number of points in a multibody setup.

# Arguments
- `multibody_setup::MultibodySetup`: [`MultibodySetup`](@ref).

# Returns
- `n_points::Int`: The sum of all points from all bodies in the multibody setup.

# Examples
```julia-repl
julia> ms = MultibodySetup(:a => body_a, :b => body_b)
2000-point MultibodySetup:
  1000-point Body{BBMaterial{NoCorrection}} with name `a`
  1000-point Body{BBMaterial{NoCorrection}} with name `b`

julia> n_points(ms)
2000
```

"""
function n_points end

@inline n_points(body::AbstractBody) = body.n_points
