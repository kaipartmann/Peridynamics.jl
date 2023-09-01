using CairoMakie
using DataFrames

TITLE_TL = Dict("bbvv" => "bond-based, velocity verlet",
                "bbdr" => "bond-based, dynamic relaxation",
                "bbmmvv" => "bond-based multi-material, velocity verlet",
                "cpdvv" => "continuum-kinematics-based, velocity verlet",
                "contact" => "contact analysis")

function get_data(logfiles)
    simsets = parse_filename.(logfiles)
    display(simsets)
    timings = get_timings.(logfiles, simsets)
    display(timings)
    return DataFrame(merge.(simsets, timings))
end

function parse_filename(file)
    parts = split(first(splitext(basename(file))), "_")
    @assert length(parts) == 3
    @assert startswith(parts[2], "t=")
    @assert startswith(parts[3], "e=")
    benchmark::String = parts[1]
    t = parse(Int, last(split(parts[2], "=")))
    e = parse(Bool, last(split(parts[3], "=")))
    return (; file = file, benchmark = benchmark, t = t, e = e)
end

function get_timings(file, simset)
    @info "analyzing file..." file=basename(file)
    lines = readlines(file)
    time = extract_total_time(lines)
    timings = extract_timings(lines, simset.benchmark)
    return (; walltime = time, timings...)
end

function extract_total_time(lines)
    timing_line = findfirst(x -> startswith(x, "✓ Simulation completed after"), lines)
    time = parse(Float64, split(lines[timing_line])[end - 1])
    return time
end

function extract_timer(lines, section; offset::Int = 3)
    timing_line = findfirst(x -> contains(x, section), lines)
    isnothing(timing_line) && error("section $section not found in logfile!\n")
    parts = split(lines[timing_line])
    # time
    time_str, unit_str = parse_time_with_unit(parts[offset])
    if unit_str == "ms"
        @assert time_str * "ms" == parts[offset]
        time = parse(Float64, time_str) * 1e-3 # convert to seconds
    elseif unit_str == "μs"
        @assert time_str * "μs" == parts[offset]
        time = parse(Float64, time_str) * 1e-6 # convert to seconds
    elseif endswith(parts[offset], "s")
        @assert time_str * "s" == parts[offset]
        time = parse(Float64, time_str)
    else
        error("time in section $section ends with unknown unit $unit_str\n")
    end
    # %tot
    @assert endswith(parts[offset + 1], "%")
    percenttot_string = parts[offset + 1][begin:(end - 1)]
    @assert percenttot_string * "%" == parts[offset + 1]
    percenttot = parse(Float64, percenttot_string)
    return time, percenttot
end

function parse_time_with_unit(input_string)
    # Define a regular expression pattern to match the numeric value and unit
    pattern = r"([\d\.]+)(\p{L}+)"

    # Use the match function to extract the parts of the string
    match_result = match(pattern, input_string)

    # Check if a match was found
    if match_result !== nothing
        # Extract the numeric value and unit from the match
        numeric_value = match_result.captures[1]
        unit = match_result.captures[2]

        # Create a tuple with the extracted values
        result_tuple = (numeric_value, unit)

        return result_tuple
    else
        error("No match found in the input string!\n")
    end
end

function extract_timings(lines, benchmark)
    timer1 = extract_timer(lines, "init_body")
    timer2 = extract_timer(lines, "compute_forcedensity!")
    timer3 = extract_timer(lines, "export_vtk")
    if benchmark in ("bbvv", "bbdr", "cpdvv", "bbmmvv")
        timer4 = extract_timer(lines, "update_thread_cache!")
        timer5 = extract_timer(lines, "define cracks"; offset = 4)
        timer6 = (0.0, 0.0)
    elseif benchmark in ("contact",)
        timer4 = extract_timer(lines, "update_thread_cache_contact!")
        timer5 = (0.0, 0.0)
        timer6 = extract_timer(lines, "compute_contactforcedensity!")
    end
    return (t_init_body = timer1, t_comp_fd = timer2, t_export_vtk = timer3,
            t_update_cache = timer4, t_def_cracks = timer5, t_comp_cfd = timer6)
end

function walltime_plots(df_all, benchmark)
    df_etrue = subset(df_all, :benchmark => x -> x .== benchmark, :e => x -> x .== true)
    df_efalse = subset(df_all, :benchmark => x -> x .== benchmark, :e => x -> x .== false)
    @assert allequal(df_etrue.benchmark)
    @assert allequal(df_efalse.benchmark)
    @assert allequal(df_etrue.e)
    @assert allequal(df_efalse.e)
    @assert df_etrue.benchmark == df_efalse.benchmark
    sort!(df_etrue, :t)
    sort!(df_efalse, :t)
    fig = walltime_plot(df_etrue, df_efalse)
    display(fig)
    return nothing
end

function walltime_plot(df_etrue, df_efalse)
    benchmark = df_etrue.benchmark[1]
    colors = Makie.wong_colors()
    ymax = maximum([df_etrue.walltime; df_efalse.walltime])
    fig = Figure()
    ax = Axis(fig[1, 1];
              title = TITLE_TL[benchmark],
              xlabel = "number of threads",
              ylabel = "walltime [s]",
              topspinevisible = false,
              rightspinevisible = false,
              xautolimitmargin = (0.02f0, 0.02f0),
              limits = (nothing, nothing, -0.03ymax, 1.03ymax),
              yticks = LinearTicks(6),
              xticks = [1, 2, 4, 8, 16, 32, 64])
    scatterlines!(ax, df_etrue.t, df_etrue.walltime, color = colors[1],
                  label = "with vtk export")
    scatterlines!(ax, df_efalse.t, df_efalse.walltime, color = colors[2], linestyle = :dash,
                  marker = :star4, label = "no vtk export")
    axislegend(ax; position = :rt)
    save(joinpath(IMG_PATH, "pdf", "$(benchmark)_walltime.pdf"), fig)
    save(joinpath(IMG_PATH, "png", "$(benchmark)_walltime.png"), fig; px_per_unit = 3)
    return fig
end

function speedup_plots(df_all, benchmark)
    df_etrue = subset(df_all, :benchmark => x -> x .== benchmark, :e => x -> x .== true)
    df_efalse = subset(df_all, :benchmark => x -> x .== benchmark, :e => x -> x .== false)
    @assert allequal(df_etrue.benchmark)
    @assert allequal(df_efalse.benchmark)
    @assert allequal(df_etrue.e)
    @assert allequal(df_efalse.e)
    @assert df_etrue.benchmark == df_efalse.benchmark
    sort!(df_etrue, :t)
    sort!(df_efalse, :t)
    @assert df_etrue[1, :t] == 1
    @assert df_efalse[1, :t] == 1
    df_etrue[!, "speedup"] = df_etrue[1, :walltime] ./ df_etrue[!, :walltime]
    df_efalse[!, "speedup"] = df_efalse[1, :walltime] ./ df_efalse[!, :walltime]
    fig = speedup_plot(df_etrue, df_efalse)
    display(fig)
    return nothing
end

function speedup_plot(df_etrue, df_efalse)
    benchmark = df_etrue.benchmark[1]
    colors = Makie.wong_colors()
    fig = Figure()
    ax = Axis(fig[1, 1];
              title = TITLE_TL[benchmark],
              xlabel = "number of threads",
              ylabel = "speedup",
              topspinevisible = false,
              rightspinevisible = false,
              xautolimitmargin = (0.02f0, 0.02f0),
              yautolimitmargin = (0.03f0, 0.03f0),
              yticks = 1:15,
              xticks = [1, 2, 4, 8, 16, 32, 64])
    # ideal
    max_speedup = maximum([df_etrue.speedup; df_efalse.speedup])
    lines!(ax, 1:0.01:max_speedup, 1:0.01:max_speedup, label = "ideal", linestyle = :dash,
           color = :black)
    scatterlines!(ax, df_etrue.t, df_etrue.speedup, color = colors[1],
                  label = "with vtk export")
    scatterlines!(ax, df_efalse.t, df_efalse.speedup, color = colors[2], linestyle = :dash,
                  marker = :star4, label = "no vtk export")
    axislegend(ax; position = :rb)
    save(joinpath(IMG_PATH, "pdf", "$(benchmark)_speedup.pdf"), fig)
    save(joinpath(IMG_PATH, "png", "$(benchmark)_speedup.png"), fig; px_per_unit = 3)
    return fig
end

function profiling_plots(df_all, benchmark)
    df_etrue = subset(df_all, :benchmark => x -> x .== benchmark, :e => x -> x .== true)
    df_efalse = subset(df_all, :benchmark => x -> x .== benchmark, :e => x -> x .== false)
    @assert allequal(df_etrue.benchmark)
    @assert allequal(df_efalse.benchmark)
    @assert allequal(df_etrue.e)
    @assert allequal(df_efalse.e)
    @assert df_etrue.benchmark == df_efalse.benchmark
    sort!(df_etrue, :t)
    sort!(df_efalse, :t)
    @assert df_etrue[1, :t] == 1
    @assert df_efalse[1, :t] == 1
    fig = profiling_plot(df_etrue, df_efalse)
    display(fig)
    return nothing
end

function profiling_plot(df_etrue, df_efalse)
    benchmark = df_etrue.benchmark[1]
    colors = Makie.wong_colors()
    fig = Figure(resolution = (1000, 600))
    ax = Axis(fig[1, 1];
              title = TITLE_TL[benchmark],
              xlabel = "number of threads",
              ylabel = "percentage of walltime [%]",
              topspinevisible = false,
              rightspinevisible = false,
              xautolimitmargin = (0.02f0, 0.02f0),
              yautolimitmargin = (0.03f0, 0.03f0),
              yticks = LinearTicks(6),
              xticks = [1, 2, 4, 8, 16, 32, 64])

    sce1 = scatterlines!(ax, df_etrue.t, last.(df_etrue.t_init_body), color = colors[1])
    sce2 = scatterlines!(ax, df_etrue.t, last.(df_etrue.t_comp_fd), color = colors[2])
    sce3 = scatterlines!(ax, df_etrue.t, last.(df_etrue.t_update_cache), color = colors[3])
    lg1 = [sce1, sce2, sce3]
    lg1_lbl = ["init. body", "comp. forcedensity", "update cache"]
    if is_timed(df_etrue.t_def_cracks)
        sce41 = scatterlines!(ax, df_etrue.t, last.(df_etrue.t_def_cracks),
                              color = colors[4])
        push!(lg1, sce41)
        push!(lg1_lbl, "def. cracks")
    end
    if is_timed(df_etrue.t_comp_cfd)
        sce42 = scatterlines!(ax, df_etrue.t, last.(df_etrue.t_comp_cfd), color = colors[4])
        push!(lg1, sce42)
        push!(lg1_lbl, "comp. contact")
    end
    sce5 = scatterlines!(ax, df_etrue.t, last.(df_etrue.t_export_vtk), color = colors[5])
    push!(lg1, sce5)
    push!(lg1_lbl, "vtk export")

    sc1 = scatterlines!(ax, df_efalse.t, last.(df_efalse.t_init_body), color = colors[1],
                        linestyle = :dash, marker = :star4)
    sc2 = scatterlines!(ax, df_efalse.t, last.(df_efalse.t_comp_fd), color = colors[2],
                        linestyle = :dash, marker = :star4)
    sc3 = scatterlines!(ax, df_efalse.t, last.(df_efalse.t_update_cache), color = colors[3],
                        linestyle = :dash, marker = :star4)
    lg2 = [sc1, sc2, sc3]
    lg2_lbl = ["init. body", "comp. forcedensity", "update cache"]
    if is_timed(df_etrue.t_def_cracks)
        sc41 = scatterlines!(ax, df_etrue.t, last.(df_etrue.t_def_cracks),
                             color = colors[4], linestyle = :dash, marker = :star4)
        push!(lg2, sc41)
        push!(lg2_lbl, "def. cracks")
    end
    if is_timed(df_etrue.t_comp_cfd)
        sc42 = scatterlines!(ax, df_etrue.t, last.(df_etrue.t_comp_cfd), color = colors[4],
                             linestyle = :dash, marker = :star4)
        push!(lg2, sc42)
        push!(lg2_lbl, "comp. contact")
    end

    Legend(fig[1, 2], [lg1, lg2], [lg1_lbl, lg2_lbl], ["with vtk export", "no vtk export"];
           titlehalign = :left)

    save(joinpath(IMG_PATH, "pdf", "$(benchmark)_profiling.pdf"), fig)
    save(joinpath(IMG_PATH, "png", "$(benchmark)_profiling.png"), fig; px_per_unit = 3)

    return fig
end

function is_timed(t::Vector{Tuple{Float64, Float64}})
    if isapprox(first(first(t)), 0) && isapprox(first(last(t)), 0)
        return false
    else
        return true
    end
end

#==========================================================================================#
#--- IMPORT DATA
LOGS_PATH = joinpath(@__DIR__, "..", "benchmark_results_logs")
@assert isdir(LOGS_PATH)
LOGFILES = filter(x -> endswith(x, ".log"), readdir(LOGS_PATH; join = true))

#--- CREATE IMG DIRECTORY
const IMG_PATH = joinpath(@__DIR__, "..", "img")
rm(IMG_PATH, recursive = true, force = true)
mkpath(joinpath(IMG_PATH, "pdf"))
mkpath(joinpath(IMG_PATH, "png"))

#--- CREATE DATAFRAME
df = get_data(LOGFILES)

#--- MAKE THE PLOTS
benchmarks = unique(df.benchmark)
for b in benchmarks
    walltime_plots(df, b)
    speedup_plots(df, b)
    profiling_plots(df, b)
end
