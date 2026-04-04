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

Now that the simulations are complete, the post-processing can begin. This is necessary in order to calculate and display the actual values needed from the raw data generated. In this example the needed values are the theoretical wave velocity `c_0`, the measured wave velocity from the simulation `c_w` and the difference between them two `dc` as well as the percentage deviation `dcp`.

The complete postprocessing code is shown below and will be discussed in the following sections. 

```julia
default_result = (; c_0=NaN, c_w=NaN, dc=NaN, dcp=NaN)
results = process_each_job(study, default_result) do job, setup
    (; E, rho, T, npyz) = setup
    Δx = 0.002 / npyz

    default_exp_ret = (; time=NaN, wave_position=NaN, amplitude_ux=NaN)
    res = process_each_export(job, default_exp_ret) do r0, r, id
        t = first(r[:time])
        t < 1.5T && return default_exp_ret  
        ymax = @views maximum(r0[:position][2, :]) 
        pids = findall(eachcol(r0[:position])) do point 
            isapprox(point[2], ymax; atol=1.01Δx/2) 
        end
        û, pids_id = @views findmax(r[:displacement][1, pids])
        x̂ = r[:position][1, pids[pids_id]] 
        return (; time=t, wave_position=x̂, amplitude_ux=û)  
    end

    allres = (; (k => getfield.(res, k) for k in keys(default_exp_ret))...)  
    (; time, wave_position, amplitude_ux) = allres  
    idxs = findall(x -> !isnan(x), time)

    c_0 = sqrt(E / rho)  
    c_w = calc_velocity(time[idxs], wave_position[idxs], amplitude_ux[idxs])  
    dc = c_w - c_0  
    dcp = 100 * dc / c_0 

    @mpiroot begin  
            bold = true  
            simname = basename(job.options.root)  
            printstyled("\nRESULTS FOR SIMULATION: `$(simname)`:\n"; bold)  
            printstyled(@sprintf("c_0 = %8.2f m/s\n", c_0); color=:green, bold)  
            printstyled(@sprintf("c_w = %8.2f m/s\n", c_w); color=:blue, bold)  
            printstyled(@sprintf("Δc  = %8.2f m/s (%.3f %%)\n", dc, dcp); color=:red, bold)  
    end  

    return (; c_0, c_w, dc, dcp)  
end  
```


#### 3.1 Postprocessing function

First, `default_result` is used to define a default result structure (a NamedTuple). 
The function `process_each_job` then runs through all simulations (`study`) and executes the postprocessing code listed in sections 3.2-3 for each individual setup. The results are stored in `res`.
The `do` block receives `job` and `setup` as parameters and returns the calculated results for each simulation.


```julia
default_result = (; c_0=NaN, c_w=NaN, dc=NaN, dcp=NaN)
results = process_each_job(study, default_result) do job, setup
    (; E, rho, T, npyz) = setup
    Δx = 0.002 / npyz
    
    # ... postprocessing code for each job (see sections 3.2-3.3) ...
    
    return (; c_0, c_w, dc, dcp)
end
```


#### 3.2 Identify irrelevant data 

Now, invalid or irrelevant data that is not to be processed further is collected.
For this purpose, another named tuple `default_exp_ret` is defined as the result structure. The function `process_each_export` iterates over the stored export data of the job and then stores only the relevant results in `res`. 

```julia
default_exp_ret = (; time=NaN, wave_position=NaN, amplitude_ux=NaN)
res = process_each_export(job, default_exp_ret) do r0, r, id
        t = first(r[:time])
        t < 1.5T && return default_exp_ret  
        ymax = @views maximum(r0[:position][2, :]) 
        pids = findall(eachcol(r0[:position])) do point 
            isapprox(point[2], ymax; atol=1.01Δx/2) 
        end
        û, pids_id = @views findmax(r[:displacement][1, pids])
        x̂ = r[:position][1, pids[pids_id]] 
        return (; time=t, wave_position=x̂, amplitude_ux=û)  
end
```


#### 3.3 Transformation of the results

Now, from the results stored in `res`, the vectors of the fields (`time`, `wave_position`, `amplitude_ux`) are collected in a NamedTuple, unpacked into local vectors, and the indices of the non-NaN entries are determined with `findall` so that only valid measurement points are used later.

```julia
allres = (; (k => getfield.(res, k) for k in keys(default_exp_ret))...)  
(; time, wave_position, amplitude_ux) = allres  
idxs = findall(x -> !isnan(x), time)
```



#### 3.4 Calculate and display the values

Now the actual values can be calculated and displayed:

```julia
    c_0 = sqrt(E / rho)  
    c_w = calc_velocity(time[idxs], wave_position[idxs], amplitude_ux[idxs])  
    dc = c_w - c_0  
    dcp = 100 * dc / c_0 

    @mpiroot begin  
            bold = true  
            simname = basename(job.options.root)  
            printstyled("\nRESULTS FOR SIMULATION: `$(simname)`:\n"; bold)  
            printstyled(@sprintf("c_0 = %8.2f m/s\n", c_0); color=:green, bold)  
            printstyled(@sprintf("c_w = %8.2f m/s\n", c_w); color=:blue, bold)  
            printstyled(@sprintf("Δc  = %8.2f m/s (%.3f %%)\n", dc, dcp); color=:red, bold)  
    end  
```
