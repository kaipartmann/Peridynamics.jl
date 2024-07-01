"""
    Body(material, position, volume)

Creates a body for use in peridynamic calculation

# Arguments

- `material::AbstractMaterial`: Specifies which material model is used
- `position::AbstractMatrix`: 3×n matrix with position of each point
- `volume::AbstractVector`: Vector with volume of each point

# Throws

- Error if n_points = 0
- `DimensionMismatch`: Error if dimension of position != 3

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

---

!!! warning "Internal use only"
    Please note that the fields are intended for internal use only. They are *not* part of
    the public API of Peridynamics.jl, and thus can be altered (or removed) at any time
    without it being considered a breaking change.

```julia
Body{Material,PointParameters}
```

# Type Parameters

- `Material <: AbstractMaterial`: Type of the specified material model
- `PointParameters <: AbstractPointParameters`: Type of the point parameters

# Fields

- `mat::Material`: Specified material model
- `n_points::Int`: Number of material points that represent the body
- `position::Matrix{Float64}`: 3×n_points matrix with position for each point
- `volume::Vector{Float64}`: Vector with volume for each point
- `fail_permit::Vector{Bool}`: Vector that describes if failure is allowed for each point
- `point_sets::Dict{Symbol,Vector{Int}}`: Dictionary containing the defined point sets
- `point_params::Vector{PointParameters}`: Vector with material parameter sets
- `params_map::Vector{Int}`: Vector that assigns a material parameter set to each point
- `single_dim_bcs::Vector{SingleDimBC}`: Vector with defined boundary conditions
- `single_dim_ics::Vector{SingleDimIC}`: Vector with defined initial conditions
- `point_sets_precracks::Vector{PointSetsPreCrack}`: Vector with defined cracks
"""
struct Body{M<:AbstractMaterial,P<:AbstractPointParameters,PC} <: AbstractBody{M}
    name::Ref{Symbol}
    mat::M
    point_cloud::PC
    fail_permit::Vector{Bool}
    point_params::Vector{P}
    params_map::Vector{Int}
    single_dim_bcs::Vector{SingleDimBC}
    posdep_single_dim_bcs::Vector{PosDepSingleDimBC}
    single_dim_ics::Vector{SingleDimIC}
    point_sets_precracks::Vector{PointSetsPreCrack}

    function Body(mat::M, point_cloud::PC) where {M,PC}
        name = Symbol("")

        fail_permit = fill(true, length(point_cloud.n_points))

        P = point_param_type(mat)
        point_params = Vector{P}()
        params_map = zeros(Int, point_cloud.n_points)

        single_dim_bcs = Vector{SingleDimBC}()
        posdep_single_dim_bcs = Vector{PosDepSingleDimBC}()
        single_dim_ics = Vector{SingleDimIC}()
        point_sets_precracks = Vector{PointSetsPreCrack}()

        new{M,P,PC}(name, mat, point_cloud, fail_permit, point_params, params_map,
                    single_dim_bcs, posdep_single_dim_bcs, single_dim_ics,
                    point_sets_precracks)
    end
end

function Base.show(io::IO, body::AbstractBody)
    print(io, n_points(body), "-point Body{", material_type(body), "}")
    if has_name(body)
        print(io, " with name `", get_name(body), "`")
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", body::AbstractBody)
    if get(io, :compact, false)
        show(io, body)
        return nothing
    end
    print(io, n_points(body), "-point Body{", material_type(body), "}")
    if has_name(body)
        print(io, " with name `", get_name(body), "`")
    end
    if has_point_sets(body) || has_params(body) || has_conditions(body)
        print(io, ":")
    end
    for (name, points) in body.point_cloud.point_sets
        print(io, "\n  ", length(points), "-point set `", name, "`")
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
    end
    if has_ics(body)
        print(io, "\n  ", n_ics(body), " initial condition(s):")
        for bc in body.single_dim_ics
            print(io, "\n    ")
            show(io, bc)
        end
    end
    if has_precracks(body)
        print(io, "\n  ", n_precracks(body), " predefined crack(s)")
    end
    n_failure_permit_false = body.n_points - sum(body.fail_permit)
    if n_failure_permit_false > 0
        print(io, "\n  ", n_failure_permit_false, " points with no failure permission")
    end
    return nothing
end

@inline material_type(::AbstractBody{M}) where {M} = M

function pre_submission_check(b::Body)
    #TODO: check if everything is defined for job submission!
    return nothing
end

@inline function get_point_param(b::AbstractBody, key::Symbol, i::Int)
    return getfield(b.point_params[b.params_map[i]], key)
end

@inline storage_type(b::AbstractBody, ts::AbstractTimeSolver) = storage_type(b.mat, ts)

function log_spatial_setup(options::AbstractJobOptions, body::AbstractBody;
                           bodyname::AbstractString="")
    (; position, point_sets) = body.point_cloud
    msg = "BODY"
    isempty(bodyname) || (msg *= " `" * bodyname * "`")
    msg *= "\n"
    msg *= "  POINT CLOUD\n"
    msg *= msg_qty("number of points", body.n_points; indentation=4)
    @views min_x, max_x = minimum(position[1, :]), maximum(position[1, :])
    @views min_y, max_y = minimum(position[2, :]), maximum(position[2, :])
    @views min_z, max_z = minimum(position[3, :]), maximum(position[3, :])
    minmax_x = @sprintf("%.7g, %.7g", min_x, max_x)
    minmax_y = @sprintf("%.7g, %.7g", min_y, max_y)
    minmax_z = @sprintf("%.7g, %.7g", min_z, max_z)
    msg *= msg_qty("min, max values x-direction", minmax_x; indentation=4)
    msg *= msg_qty("min, max values y-direction", minmax_y; indentation=4)
    msg *= msg_qty("min, max values z-direction", minmax_z; indentation=4)
    msg *= "  POINT SETS\n"
    for (key, points) in point_sets
        descr = @sprintf("number of points in set `%s`", string(key))
        msg *= msg_qty(descr, length(points); indentation=4)
    end
    has_conditions(body) && (msg *= "  CONDITIONS\n")
    has_bcs(body) && (msg *= msg_qty("number of BC's", n_bcs(body); indentation=4))
    has_ics(body) && (msg *= msg_qty("number of IC's", n_ics(body); indentation=4))
    n_point_params = length(body.point_params)
    if length(body.point_params) == 1
        msg *= "  MATERIAL\n"
        msg *= msg_qty("material type", material_type(body); indentation=4)
        msg *= log_material_parameters(first(body.point_params); indentation=4)
    elseif n_point_params > 1
        for (i, params) in enumerate(body.point_params)
            msg *= @sprintf("  MATERIAL %d\n", i)
            msg *= log_material_parameters(params; indentation=4)
        end
    end
    log_it(options, msg)
    return nothing
end

function maximum_horizon(b::AbstractBody)
    n_params = length(b.point_params)
    n_params == 1 && return get_horizon(first(b.point_params))
    n_params == 0 && error("body has no material parameters!\n")
    δmax = 0.0
    for point_id in eachindex(b.point_cloud.volume)
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

@inline n_points(body::AbstractBody) = n_points(body.point_cloud)

@inline has_point_sets(body::AbstractBody) = has_point_sets(body.point_cloud)
@inline each_user_point_set(body::AbstractBody) = each_user_point_set(body.point_cloud)

@inline n_params(body::AbstractBody) = length(body.point_params)
@inline has_params(body::AbstractBody) = !isempty(body.point_params)

@inline function n_bcs(body::AbstractBody)
    return length(body.single_dim_bcs) + length(body.posdep_single_dim_bcs)
end
@inline n_ics(body::AbstractBody) = length(body.single_dim_ics)

@inline has_bcs(body::AbstractBody) = n_bcs(body) > 0 ? true : false
@inline has_ics(body::AbstractBody) = n_ics(body) > 0 ? true : false
@inline has_conditions(body::AbstractBody) = has_bcs(body) || has_ics(body)

@inline n_precracks(body::AbstractBody) = length(body.point_sets_precracks)
@inline has_precracks(body::AbstractBody) = n_precracks(body) > 0 ? true : false
