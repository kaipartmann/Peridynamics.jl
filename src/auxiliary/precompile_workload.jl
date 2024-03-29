@setup_workload begin
    if mpi_run()
        msg = "precompilation triggered during mpirun!\n"
        msg *= "The precompilation is not safe to use with MPI!\n"
        msg *= "Trigger package precompilation manually and then restart the mpirun!\n"
        error(msg)
    end
    l, Δx, δ, a = 1.0, 1 / 4, 3.015 / 4, 0.5
    root = joinpath(@__DIR__, "temp_precompilation")
    pos, vol = uniform_box(l, l, l, Δx)
    path_bb = joinpath(root, "BB")
    path_osb = joinpath(root, "OSB")
    path_nosb = joinpath(root, "NOSB")

    @compile_workload begin
        b1 = Body(BBMaterial(), pos, vol)
        material!(b1; horizon=3.015Δx, E=2.1e5, rho=8e-6, Gc=2.7)
        point_set!(p -> p[1] ≤ -l / 2 + a && 0 ≤ p[2] ≤ 2δ, b1, :set_a)
        point_set!(p -> p[1] ≤ -l / 2 + a && -2δ ≤ p[2] < 0, b1, :set_b)
        precrack!(b1, :set_a, :set_b)
        point_set!(y -> y > l / 2 - Δx, b1, :set_top)
        point_set!(y -> y < -l / 2 + Δx, b1, :set_bottom)
        velocity_ic!(b1, :set_bottom, :x, 0)
        velocity_ic!(b1, :set_bottom, :y, 0)
        velocity_ic!(b1, :set_bottom, :z, 0)
        velocity_ic!(b1, :set_top, 1, 0)
        velocity_ic!(b1, :set_top, 2, 0)
        velocity_ic!(b1, :set_top, 3, 0)
        velocity_bc!(t -> 0, b1, :set_bottom, :x)
        velocity_bc!(t -> -30, b1, :set_bottom, :y)
        velocity_bc!(t -> 0, b1, :set_bottom, :z)
        velocity_bc!(t -> 0, b1, :set_top, 1)
        velocity_bc!(t -> 30, b1, :set_top, 2)
        velocity_bc!(t -> 0, b1, :set_top, 3)
        forcedensity_bc!(t -> 0, b1, :set_top, :x)
        forcedensity_bc!(t -> 0, b1, :set_top, :y)
        forcedensity_bc!(t -> 0, b1, :set_top, :z)
        forcedensity_bc!(t -> 0, b1, :set_bottom, 1)
        forcedensity_bc!(t -> 0, b1, :set_bottom, 2)
        forcedensity_bc!(t -> 0, b1, :set_bottom, 3)

        b2 = Body(OSBMaterial(), pos, vol)
        material!(b2; horizon=3.015Δx, E=2.1e5, nu=0.25, rho=8e-6, Gc=2.7)
        point_set!(p -> p[1] ≤ -l / 2 + a && 0 ≤ p[2] ≤ 2δ, b2, :set_a)
        point_set!(p -> p[1] ≤ -l / 2 + a && -2δ ≤ p[2] < 0, b2, :set_b)
        failure_permit!(b2, :set_a, false)
        precrack!(b2, :set_a, :set_b)
        point_set!(p -> p[2] > l / 2 - Δx, b2, :set_top)
        point_set!(p -> p[2] < -l / 2 + Δx, b2, :set_bottom)
        velocity_bc!(t -> -30, b2, :set_bottom, :y)
        velocity_bc!(t -> 30, b2, :set_top, 2)

        b3 = Body(NOSBMaterial(), pos, vol)
        material!(b3; horizon=3.015Δx, E=2.1e5, nu=0.25, rho=8e-6, Gc=2.7)
        point_set!(p -> p[1] ≤ -l / 2 + a && 0 ≤ p[2] ≤ 2δ, b3, :set_a)
        point_set!(p -> p[1] ≤ -l / 2 + a && -2δ ≤ p[2] < 0, b3, :set_b)
        failure_permit!(b3, :set_a, false)
        precrack!(b3, :set_a, :set_b)
        point_set!(p -> p[2] > l / 2 - Δx, b3, :set_top)
        point_set!(p -> p[2] < -l / 2 + Δx, b3, :set_bottom)
        velocity_bc!(t -> -30, b3, :set_bottom, :y)
        velocity_bc!(t -> 30, b3, :set_top, 2)

        vv = VelocityVerlet(steps=2)

        submit(Job(b1, vv; path=path_bb, freq=1); quiet=true)
        submit(Job(b2, vv; path=path_osb, freq=1); quiet=true)
        submit(Job(b3, vv; path=path_nosb, freq=1); quiet=true)
    end

    rm(root; recursive=true, force=true)
end
