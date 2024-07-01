struct PointCloud <: AbstractPointCloud
    n_points::Int
    position::Matrix{Float64}
    volume::Vector{Float64}
    point_sets::Dict{Symbol,Vector{Int}}

    function PointCloud(position::AbstractMatrix, volume::AbstractVector)
        n_points = length(volume)
        check_pos_and_vol(n_points, position, volume)
        point_sets = Dict{Symbol,Vector{Int}}(:all_points => 1:length(volume))
        return new(n_points, position, volume, point_sets)
    end
end

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

@inline has_point_sets(point_cloud::AbstractPointCloud) = !isempty(point_cloud.point_sets)

@inline function each_user_point_set(point_cloud::AbstractPointCloud)
    return filter(x -> x !== :all_points, point_cloud.point_sets)
end

@inline n_points(point_cloud::AbstractPointCloud) = point_cloud.n_points
