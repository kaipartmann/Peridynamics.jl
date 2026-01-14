# Perform a parameter study

The following section explains how to perform a parameter study. 

### 1. Define the simulation conditions

Here is an example function that creates a peridynamic body and defines the boundary conditions of each simulation.
This simulation is intended to simulate wave propagation in the defined body. 

```julia
function job_xwave(setup, root)
    (; E, nu, rho, T, vmax, npyz, m) = setup
    lx, lyz = 0.2, 0.002
    Δx = lyz / npyz
    pos, vol = uniform_box(lx, lyz, lyz, Δx)
    ids = sortperm(pos[1,:])
    body = Body(BBMaterial(), pos[:, ids], vol[ids])
    horizon = m * Δx
    material!(body; horizon, rho, E, nu)
    point_set!(x -> x < -lx / 2 + 1.2Δx, body, :left)
    velocity_bc!(t -> t < T ? vmax * sin(2π / T * t) : 0.0, body, :left, :x)
    vv = VelocityVerlet(time=1e-4, safety_factor=0.7)
    simname = @sprintf("xwave_npyz-%d_m-%d", npyz, round(Int, m))
    path = joinpath(root, simname)
    job = Job(body, vv; path, freq=5, fields=(:displacement,))
    return job
end
```

The line 
```julia
(; E, nu, rho, T, vmax, npyz, m) = setup
```
creates a setup in which the material parameters are defined.


Now various setups are created in which the respective parameters are defined as desired. A separate simulation is performed for each setup. Here the parameter `m`, which controls the peridynamic horizon `h = m * Δx`, varies. 

```julia
setups = [
    (; E=210e9, nu=0.25, rho=7850.0, T=1.0e-5, vmax=2.0, npyz=4, m=3.015),
    (; E=210e9, nu=0.25, rho=7850.0, T=1.0e-5, vmax=2.0, npyz=4, m=4.015),
    (; E=210e9, nu=0.25, rho=7850.0, T=1.0e-5, vmax=2.0, npyz=4, m=5.015),
]
```

The location where the simulation results are to be stored can be specified with `joinpath`
```julia
root = joinpath("results", "xwave_study")
```

### 2. Submit study

To execute the previously defined setups one after the other, first create a study object and then call `submit!`:
```julia
study = Study(job_xwave, setups; root)
submit!(study)
```


### 3. Postprocessing

Now that the simulations are complete, post-processing can begin. This is necessary in order to calculate and display the actual values needed from the raw data generated. In this example the needed values are the theoretical wave velocity `c_0`, the measured wave velocity from the simulation `c_w` and the difference between them two `dc` as well as the percentage deviation `dcp`.

The complete postprocessing code is shown below and will be discussed in the following sections. 



```julia
default_result = (; c_0=NaN, c_w=NaN, dc=NaN, dcp=NaN)
results = process_each_job(study, default_result) do job, setup
    (; E, rho, T, npyz) = setup
    Δx = 0.002 / npyz

    post_path = joinpath(job.options.root, "post") 
    mkpath(post_path)
    wave_position_data_file = joinpath(post_path, "wavepos.csv")

    process_each_export(job; serial=true) do r0, r, id 
        t = first(r[:time])
        t < 1.5T && return nothing 
        ymax = @views maximum(r0[:position][2, :]) 
        pids = findall(eachcol(r0[:position])) do point 
            isapprox(point[2], ymax; atol=1.01Δx/2) 
        end
        û, pids_id = @views findmax(r[:displacement][1, pids])
        x̂ = r[:position][1, pids[pids_id]] 
        open(wave_position_data_file, "a+") do io 
            @printf(io, "%.12f,%.12f,%.12f\n", t, x̂, û)
        end
        return nothing
    end

    results = readdlm(wave_position_data_file, ',', Float64)
    t, x̂, û = results[:,1], results[:,2], results[:,3]
    c_0 = sqrt(E / rho)
    c_w = calc_velocity(t, x̂, û)
    dc = c_w - c_0
    dcp = 100 * dc / c_0
    printstyled("\nRESULTS FOR SIMULATION: `$(basename(job.options.root))`:\n", bold=true)
    printstyled(@sprintf("c_0 = %8.2f m/s\n", c_0); color=:green, bold=true)
    printstyled(@sprintf("c_w = %8.2f m/s\n", c_w); color=:blue, bold=true)
    printstyled(@sprintf("Δc  = %8.2f m/s (%.3f %%)\n", dc, dcp); color=:red, bold=true)

    return (; c_0, c_w, dc, dcp)
end
```

#### 3.1 Postprocessing function

First, a standard result structure for post-processing results is created. Then the `process_each_job` function iterates over all simulations (study) and executes the postprocessing code shown below for each individual setup. 

```julia
default_result = (; c_0=NaN, c_w=NaN, dc=NaN, dcp=NaN)
results = process_each_job(study, default_result) do job, setup
    (; E, rho, T, npyz) = setup
    Δx = 0.002 / npyz
    
    # ... postprocessing code for each job (see sections 3.2-3.3) ...
    
    return (; c_0, c_w, dc, dcp)
end
```

The `do` block receives `job` and `setup` as parameters and should return the computed results for each simulation.


#### 3.2 Create an output file and write the results to it

To create a CSV file (here `wavepos.csv` resp. `wave_position_data_file`) in a desired path where the data can be saved, the following can be executed.
```julia
post_path = joinpath(job.options.root, "post")
mkpath(post_path) 
wave_position_data_file = joinpath(post_path, "wavepos.csv")
```

```julia
open(wave_position_data_file, "a+") do io
    @printf(io, "%.12f,%.12f,%.12f\n", t, x̂, û)
end
```


#### 3.3 Calculate and display the values

With the `readdlm` command the created CSV file is read. Then the resulting columns are assigned to `t`, `x̂`, and `û`.

```julia
results = readdlm(wave_position_data_file, ',', Float64)
t, x̂, û = results[:,1], results[:,2], results[:,3]
```


Now the actual values can be calculated:
```julia
c_0 = sqrt(E / rho)
c_w = calc_velocity(t, x̂, û)
dc = c_w - c_0
dcp = 100 * dc / c_0 
```


Finally, the four determined values are displayed in the terminal:
```julia
printstyled("\nRESULTS FOR SIMULATION: `$(basename(job.options.root))`:\n", bold=true)
printstyled(@sprintf("c_0 = %8.2f m/s\n", c_0); color=:green, bold=true)
printstyled(@sprintf("c_w = %8.2f m/s\n", c_w); color=:blue, bold=true)
printstyled(@sprintf("Δc  = %8.2f m/s (%.3f %%)\n", dc, dcp); color=:red, bold=true)
```
