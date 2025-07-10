@testitem "DynamicRelaxation wrong input" begin
    @test_throws ArgumentError DynamicRelaxation(steps=0)
    @test_throws ArgumentError DynamicRelaxation(steps=10, stepsize=0)
    @test_throws ArgumentError DynamicRelaxation(steps=10, damping_factor=0)
end

@testitem "dynamic_relaxation_check" begin
    dr = DynamicRelaxation(steps=1)
    dr.n_steps = -1
    msg = "`n_steps` of DynamicRelaxation smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.dynamic_relaxation_check(dr)

    dr = DynamicRelaxation(steps=1)
    dr.Δt = -1
    msg = "`Δt` of DynamicRelaxation smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.dynamic_relaxation_check(dr)

    dr = DynamicRelaxation(steps=1)
    dr.Λ = -1
    msg = "`Λ` of DynamicRelaxation smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.dynamic_relaxation_check(dr)
end

@testitem "_init_density_matrix" begin
    using Peridynamics
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 0.0 1.0
                0.0 0.0 1.0 0.0]
    volume = fill(1.0, 4)
    body = Body(BBMaterial(), position, volume)
    material!(body; horizon=1.5, rho=8e-6, E=210e3)
    point_set!(body, :a, 1:2)
    material!(body, :a; horizon=1.5, rho=7.8e-6, E=200e3)

    dr = DynamicRelaxation(steps=1)
    dh = Peridynamics.threads_data_handler(body, dr, 1)
    chunk = dh.chunks[1]
    Peridynamics._init_density_matrix!(chunk,dr,chunk.paramsetup)

    @test chunk.storage.density_matrix ≈  [4.0e6  4.0e6  4.2e6  4.2e6
                                           4.0e6  4.0e6  4.2e6  4.2e6
                                           4.0e6  4.0e6  4.2e6  4.2e6]
end

@testitem "show DynamicRelaxation" begin
    io = IOBuffer()

    dr = DynamicRelaxation(steps=1)

    show(IOContext(io, :compact=>true), MIME("text/plain"), dr)
    msg = String(take!(io))
    @test contains(msg, "DynamicRelaxation(n_steps=1, Δt=1.0, Λ=1.0)")

    show(IOContext(io, :compact=>false), MIME("text/plain"), dr)
    msg = String(take!(io))
    @test contains(msg, "DynamicRelaxation:\n  n_steps  1\n  Δt       1\n  Λ        1")
end
