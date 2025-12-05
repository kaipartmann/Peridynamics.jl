@testitem "Uniform tension DynamicRelaxation" begin
    using Peridynamics.Printf
    l, w, h, Δx = 1.0, 0.1, 0.1, 1/30
    F = 2e6 # Force in N applied at the right boundary
    E = 200e9 # Young's modulus
    steps = 2000 # Number of time steps for the simulation

    # Create a uniform box with the specified dimensions and discretization
    pos, vol = uniform_box(l+3Δx, w, h, Δx; center=(-1.5Δx, 0, 0))
    body = Body(CMaterial(), pos, vol)
    material!(body; horizon=3.015Δx, rho=7850.0, E, nu=0.25, epsilon_c=1.0)
    point_set!(x -> x <-l/2, body, :left)
    point_set!(x -> x > l/2-Δx, body, :right)

    velocity_bc!(t -> 0.0, body, :left, :x)
    velocity_bc!(t -> 0.0, body, :left, :y)
    velocity_bc!(t -> 0.0, body, :left, :z)

    volume_right = sum(vol[body.point_sets[:right]])
    forcedensity_bc!(t -> F / volume_right, body, :right, :x)

    # run the simulation
    vv = DynamicRelaxation(; steps)
    path = mktempdir()
    rm(path; recursive=true, force=true)
    job = Job(body, vv; freq=steps, path, fields=(:displacement, :b_ext))
    submit(job; quiet=true)

    # Analytical solution
    A = w * h # Cross-sectional area
    Δl = F / (E * A) * l # Displacement in x-direction

    # Post-processing of the simulation results
    vtk_path = joinpath(path, "vtk")
    @test isdir(vtk_path)
    last_vtk_file = last(Peridynamics.find_vtk_files(vtk_path))
    results = read_vtk(last_vtk_file)
    u = results[:displacement]
    ux = sum(u[1, body.point_sets[:right]]) / length(body.point_sets[:right])
    err_ux = abs(ux - Δl) / Δl
    print_msg = @sprintf("x-displacement error: %.2f %%\n", 100err_ux)
    printstyled(print_msg; color=:red, bold=true)

    # Verify that the error in the x-direction displacement is within a 10% tolerance
    # This should be okay, since a very coarse discretization is used.
    # It can be shown that the error gets smaller for a finer discretization and better
    # material models.
    @test_broken err_ux < 0.02
end

@testitem "Uniform tension NewtonRaphson" begin
    using Peridynamics.Printf
    using Peridynamics: NewtonRaphson, displacement_bc!

    l, w, h, Δx = 1.0, 0.1, 0.1, 1/30
    F = 2e6 # Force in N applied at the right boundary
    E = 200e9 # Young's modulus
    steps = 20 # Number of time steps for the simulation

    # Create a uniform box with the specified dimensions and discretization
    pos, vol = uniform_box(l+3Δx, w, h, Δx; center=(-1.5Δx, 0, 0))
    body = Body(CMaterial(), pos, vol)
    material!(body; horizon=3.015Δx, rho=7850.0, E, nu=0.25)
    point_set!(x -> x <-l/2, body, :left)
    point_set!(x -> x > l/2-Δx, body, :right)

    displacement_bc!(p -> 0.0, body, :left, :x)
    displacement_bc!(p -> 0.0, body, :left, :y)
    displacement_bc!(p -> 0.0, body, :left, :z)

    volume_right = sum(vol[body.point_sets[:right]])
    forcedensity_bc!(p -> F / volume_right, body, :right, :x)

    # run the simulation
    nr = NewtonRaphson(; steps, tol=1e-3, maxiter=50)
    path = mktempdir()
    rm(path; recursive=true, force=true)
    job = Job(body, nr; freq=steps, path, fields=(:displacement,))
    submit(job; quiet=false)

    # Analytical solution
    A = w * h # Cross-sectional area
    Δl = F / (E * A) * l # Displacement in x-direction

    # Post-processing of the simulation results
    vtk_path = joinpath(path, "vtk")
    @test isdir(vtk_path)
    last_vtk_file = last(Peridynamics.find_vtk_files(vtk_path))
    results = read_vtk(last_vtk_file)
    u = results[:displacement]
    ux = sum(u[1, body.point_sets[:right]]) / length(body.point_sets[:right])
    err_ux = abs(ux - Δl) / Δl
    print_msg = @sprintf("x-displacement error: %.2f %%\n", 100err_ux)
    printstyled(print_msg; color=:red, bold=true)

    # Verify that the error in the x-direction displacement is within a 10% tolerance
    # This should be okay, since a very coarse discretization is used.
    # It can be shown that the error gets smaller for a finer discretization and better
    # material models.
    @test err_ux < 0.02
end
