using Peridynamics

RESPATH_ROOT::String = length(ARGS) ≥ 1 ? ARGS[1] : joinpath("results")
RESPATH::String = joinpath(RESPATH_ROOT, "mode_I_perf_test")

function get_job(N::Int)
    l, Δx, δ, a = 1.0, 1/N, 3.015/N, 0.5
    pos, vol = uniform_box(l, l, 0.1l, Δx)
    ids = sortperm(pos[2,:])
    mat = BBMaterial()
    b = Body(mat, pos[:, ids], vol[ids])
    material!(b; horizon=3.015Δx, E=2.1e5, rho=8e-6, Gc=2.7)
    point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, b, :set_a)
    point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, b, :set_b)
    precrack!(b, :set_a, :set_b)
    point_set!(p -> p[2] > l/2-Δx, b, :set_top)
    point_set!(p -> p[2] < -l/2+Δx, b, :set_bottom)
    velocity_bc!(t -> -30, b, :set_bottom, :y)
    velocity_bc!(t -> 30, b, :set_top, :y)
    vv = VelocityVerlet(steps=2000)
    job = Job(b, vv; path=RESPATH)
    return job
end

function main()
    # @rootdo println("-- compilation --")
    # job = get_job(40)
    # @mpitime submit(job)

    @rootdo println("-- work --")
    job = get_job(40)
    @mpitime submit(job)
    return nothing
end

using BenchmarkTools, PointNeighbors, Plots

function benchmark_peridynamics(neighborhood_search, job)
    point_decomp = Peridynamics.PointDecomposition(job.spatial_setup, 1)
    tdh = Peridynamics.ThreadsDataHandler(job.spatial_setup, job.time_solver, point_decomp)

    return @belapsed Peridynamics.calc_force_density!($(tdh.chunks[1]), $neighborhood_search);
end

function plot_benchmarks(N_start, iterations; title = "")
    neighborhood_searches_names = [
        "Original Code";;
        "GridNeighborhoodSearch";;
        "NeighborListsNeighborhoodSearch";;
        "NeighborListsNHS contiguous"
    ]

    # Multiply number of points in each iteration (roughly) by this factor
    scaling_factor = 4
    per_dimension_factor = scaling_factor^(1 / 3)
    N = [round(Int, N_start * per_dimension_factor^(iter - 1))
             for iter in 1:iterations]
    n_points = zero(N)

    times = zeros(iterations, length(neighborhood_searches_names))

    for iter in 1:iterations
        job = get_job(N[iter])
        n_points[iter] = job.spatial_setup.n_points
        search_radius = Peridynamics.maximum_horizon(job.spatial_setup)

        neighborhood_searches = [
            nothing,
            GridNeighborhoodSearch{3}(search_radius, job.spatial_setup.n_points),
            NeighborListsNeighborhoodSearch{3}(search_radius, job.spatial_setup.n_points,
                                               backend = Vector{Vector{Int}}),
            NeighborListsNeighborhoodSearch{3}(search_radius, job.spatial_setup.n_points),
        ]

        for i in eachindex(neighborhood_searches)
            neighborhood_search = neighborhood_searches[i]
            if !isnothing(neighborhood_search)
                initialize!(neighborhood_search,
                            job.spatial_setup.position, job.spatial_setup.position)
            end

            time = benchmark_peridynamics(neighborhood_search, job)
            times[iter, i] = time
            time_string = BenchmarkTools.prettytime(time * 1e9)
            println("$(neighborhood_searches_names[i])")
            println("with $(n_points[iter]) points finished in $time_string\n")
        end
    end

    plot(n_points, times,
         xaxis = :log, yaxis = :log,
         xticks = (n_points, n_points),
         xlabel = "#particles", ylabel = "Runtime [s]",
         legend = :outerright, size = (750, 400), dpi = 200,
         label = neighborhood_searches_names, title = title)
end
