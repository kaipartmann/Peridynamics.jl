@setup_workload begin
    if mpi_run()
        msg = "precompilation triggered during mpirun!\n"
        msg *= "The precompilation is not safe to use with MPI!\n"
        msg *= "Trigger package precompilation manually and then restart the mpirun!\n"
        error(msg)
    end
    root = joinpath(@__DIR__, "temp_precompilation")
    pos1, vol1 = uniform_box(1, 1, 1, 0.5; center=(0.5, 0.5, 0.5))
    pos2, vol2 = uniform_box(1, 1, 1, 0.5; center=(-0.5, 0.5, 0.5))
    path_bb = joinpath(root, "BB")
    path_osb = joinpath(root, "OSB")
    path_cc = joinpath(root, "CC")
    path_ms = joinpath(root, "MS")

    @compile_workload begin
        b1 = Body(BBMaterial(), pos1, vol1)
        material!(b1; horizon=2, E=2.1e5, rho=8e-6, Gc=2.7)
        point_set!(p -> p[1] ≤ 0.5, b1, :set_a)
        point_set!(x -> x > 0.5, b1, :set_b)
        precrack!(b1, :set_a, :set_b)
        velocity_ic!(b1, :set_a, :x, 0)
        velocity_ic!(b1, :set_a, :y, 0)
        velocity_ic!(b1, :set_a, :z, 0)
        velocity_ic!(b1, :set_b, 1, 0)
        velocity_ic!(b1, :set_b, 2, 0)
        velocity_ic!(b1, :set_b, 3, 0)
        velocity_bc!(t -> -1, b1, :set_a, :x)
        velocity_bc!(t -> 1, b1, :set_b, 1)
        forcedensity_bc!(t -> 0, b1, :set_a, :x)
        forcedensity_bc!(t -> 0, b1, :set_a, :y)
        forcedensity_bc!(t -> 0, b1, :set_a, :z)
        forcedensity_bc!(t -> 0, b1, :set_b, 1)
        forcedensity_bc!(t -> 0, b1, :set_b, 2)
        forcedensity_bc!(t -> 0, b1, :set_b, 3)

        b2 = Body(OSBMaterial{EnergySurfaceCorrection}(), pos1, vol1)
        material!(b2; horizon=2, E=2.1e5, nu=0.25, rho=8e-6, Gc=2.7)
        point_set!(p -> p[1] ≤ 0.5, b2, :set_a)
        point_set!(x -> x > 0.5, b2, :set_b)
        failure_permit!(b2, :set_a, false)
        precrack!(b2, :set_a, :set_b)
        velocity_bc!(t -> -1, b2, :set_a, :x)
        velocity_bc!(t -> 1, b2, :set_b, 1)

        b3 = Body(CMaterial(), pos1, vol1)
        material!(b3; horizon=2, E=2.1e5, nu=0.25, rho=8e-6, Gc=2.7)
        point_set!(p -> p[1] ≤ 0.5, b3, :set_a)
        point_set!(x -> x > 0.5, b3, :set_b)
        failure_permit!(b3, :set_a, false)
        precrack!(b3, :set_a, :set_b)
        velocity_bc!(t -> -1, b3, :set_a, :x)
        velocity_bc!(t -> 1, b3, :set_b, 1)

        b4 = Body(BBMaterial{EnergySurfaceCorrection}(), pos1, vol1)
        material!(b4; horizon=2, E=2.1e5, rho=8e-6, Gc=2.7)
        b5 = Body(OSBMaterial{EnergySurfaceCorrection}(), pos2, vol2)
        material!(b5; horizon=2, E=2.1e5, nu=0.25, rho=8e-6, Gc=2.7)
        velocity_ic!(b5, :all_points, :x, 100)

        ms = MultibodySetup(:b4 => b4, :b5 => b5)
        contact!(ms, :b4, :b5; radius=0.6)

        vv = VelocityVerlet(steps=2)

        submit(Job(b1, vv; path=path_bb, freq=1); quiet=true)
        submit(Job(b2, vv; path=path_osb, freq=1); quiet=true)
        submit(Job(b3, vv; path=path_cc, freq=1); quiet=true)
        submit(Job(ms, vv; path=path_ms, freq=1); quiet=true)
    end

    rm(root; recursive=true, force=true)
end
