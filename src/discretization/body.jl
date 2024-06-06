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
struct Body{M<:AbstractMaterial,P<:AbstractPointParameters} <: AbstractBody{M}
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

@inline material_type(::AbstractBody{M}) where {M} = M

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

function pre_submission_check(b::Body)
    #TODO: check if everything is defined for job submission!
    return nothing
end

@inline function get_point_param(b::AbstractBody, key::Symbol, i::Int)
    return getfield(b.point_params[b.params_map[i]], key)
end

@inline storage_type(b::AbstractBody, ts::AbstractTimeSolver) = storage_type(b.mat, ts)

function log_spatial_setup(options::AbstractOptions, body::AbstractBody)
    msg = "BODY\n"
    msg *= log_msg("number of points", body.n_points)
    msg *= log_msg("number of point sets", length(keys(body.point_sets)))
    msg *= log_msg("number of parameters", length(body.point_params))
    n_bcs = length(body.single_dim_bcs) + length(body.posdep_single_dim_bcs)
    n_bcs > 0 && (msg *= log_msg("number of BC's", n_bcs))
    n_ics = length(body.single_dim_ics)
    n_ics > 0 && (msg *= log_msg("number of IC's", length(body.single_dim_ics)))
    log_it(options, msg)
    return nothing
end
