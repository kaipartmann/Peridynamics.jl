@testitem "data boundary conditions" begin
    using Peridynamics: velocity_databc!, forcedensity_databc!
    using Peridynamics.Printf

    l, w, h, Δx = 1.0, 0.1, 0.1, 1/40
    F = 2e6 # Force in N applied at the right boundary
    E = 200e9 # Young's modulus
    steps = 2000 # Number of time steps for the simulation

    # Create a uniform box with the specified dimensions and discretization
    pos, vol = uniform_box(l+3Δx, w, h, Δx; center=(-1.5Δx, 0, 0))
    body = Body(BBMaterial{EnergySurfaceCorrection}(), pos, vol)
    material!(body; horizon=3.015Δx, rho=7850.0, E, nu=0.25, epsilon_c=1.0)
    point_set!(x -> x <-l/2, body, :left)
    point_set!(x -> x > l/2-Δx, body, :right)

    # Declare the dimensions involved in boundary conditions.
    # 1 represents the x-direction, and 2, 3 represent the yz directions.
    # If all three dimensions are involved, dims = UInt8[1, 2, 3].
    dims_f = [:x]
    dims_v = [1, 2, 3]

    # Define a matrix to describe the boundary conditions.
    # The number of rows is equal to the number of dimensions,
    # and the columns are initialized for all points.
    f_matrix = zeros(length(dims_f), body.n_points)
    v_matrix = zeros(length(dims_v), body.n_points)
    # Assign values to the boundary points in the matrix according to the boundary
    # conditions.
    volume_right = sum(vol[body.point_sets[:right]])
    for i in body.point_sets[:all_points]
        if i ∈ body.point_sets[:right]
            f_matrix[1, i] = F / volume_right
        elseif i ∈ body.point_sets[:left]
            v_matrix[:, i] .= 0.0 # Set boundary condition for fix layer
        end
    end

    velocity_databc!(body, v_matrix, :left, dims_v)
    forcedensity_databc!(body, f_matrix, :right, dims_f)

    # test show method for DataBCs
    io = IOBuffer()
    show(IOContext(io, :compact=>false), MIME("text/plain"), body)
    msg = String(take!(io))
    @test contains(msg, "2 boundary condition(s):")
    @test contains(msg, "Data BC on velocity: point_set=left, dims=UInt8[0x01, 0x02, 0x03]")
    @test contains(msg, "Data BC on force density: point_set=right, dims=UInt8[0x01]")

    # run the simulation
    vv = DynamicRelaxation(; steps)
    path = mktempdir()
    rm(path; recursive=true, force=true)
    job = Job(body, vv; freq=steps, path, fields=(:displacement, :b_ext))
    submit(job)

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
    @test err_ux < 0.1
end
