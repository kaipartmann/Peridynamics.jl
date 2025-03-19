@testitem "LinearElastic" begin
    using Peridynamics.StaticArrays

    model = LinearElastic()
    mat = CMaterial(; model)
    pos, vol = uniform_box(1.0, 1.0, 1.0, 0.5)
    body = Body(mat, pos, vol)
    material!(body, horizon=1, rho=1, E=1, nu=0.25)
    params = body.point_params[1]
    decomp = Peridynamics.PointDecomposition(body, 1)
    system = Peridynamics.get_system(body, decomp, 1)
    solver = VelocityVerlet(steps=2)
    storage = Peridynamics.CStorage(mat, solver, system)

    F = @SMatrix [1 0 0; 0 1 0; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ zero(SMatrix{3,3,Float64,9}) atol=eps()

    F = @SMatrix [2 0 0; 0 2 0; 0 0 2]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [3 0 0; 0 3 0; 0 0 3]

    F = @SMatrix [1 0.5 0; 0 1 0; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [0.05 0.2 0; 0.2 0.15 0; 0 0 0.05]

    F = @SMatrix [1 0 0; 0 1 0.5; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [0.05 0 0; 0 0.05 0.2; 0 0.2 0.15]

    F = @SMatrix [1 0 0.5; 0 1 0; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [0.05 0 0.2; 0 0.05 0.0; 0.2 0 0.15]

    F = @SMatrix [1 0 0; 0.5 1 0; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [0.15 0.2 0; 0.2 0.05 0; 0 0 0.05]

end

@testitem "NeoHooke" begin
    using Peridynamics.StaticArrays

    model = NeoHooke()
    mat = CMaterial(; model)
    pos, vol = uniform_box(1.0, 1.0, 1.0, 0.5)
    body = Body(mat, pos, vol)
    material!(body, horizon=1, rho=1, E=1, nu=0.25)
    params = body.point_params[1]
    decomp = Peridynamics.PointDecomposition(body, 1)
    system = Peridynamics.get_system(body, decomp, 1)
    solver = VelocityVerlet(steps=2)
    storage = Peridynamics.CStorage(mat, solver, system)

    F = @SMatrix [1 0 0; 0 1 0; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ zero(SMatrix{3,3,Float64,9}) atol=eps()

    F = @SMatrix [2 0 0; 0 2 0; 0 0 2]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [1.0158883083359673 0.0 0.0;
                        0.0 1.0158883083359673 0.0;
                        0.0 0.0 1.0158883083359673]

    F = @SMatrix [1 0.5 0; 0 1 0; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [0.0 0.2 0.0; 0.2 0.0 0.0; 0.0 0.0 0.0]

    F = @SMatrix [1 0 0; 0 1 0.5; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [0.0 0.0 0.0; 0.0 0.0 0.2; 0.0 0.2 0.0]

    F = @SMatrix [1 0 0.5; 0 1 0; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [0.0 0.0 0.2; 0.0 0.0 0.0; 0.2 0.0 0.0]

    F = @SMatrix [1 0 0; 0.5 1 0; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [0.0 0.2 0.0; 0.2 0.0 0.0; 0.0 0.0 0.0]
end

@testitem "MooneyRivlin" begin
    using Peridynamics.StaticArrays

    model = MooneyRivlin()
    mat = CMaterial(; model)
    pos, vol = uniform_box(1.0, 1.0, 1.0, 0.5)
    body = Body(mat, pos, vol)
    material!(body, horizon=1, rho=1, E=1, nu=0.25)
    params = body.point_params[1]
    decomp = Peridynamics.PointDecomposition(body, 1)
    system = Peridynamics.get_system(body, decomp, 1)
    solver = VelocityVerlet(steps=2)
    storage = Peridynamics.CStorage(mat, solver, system)

    F = @SMatrix [1 0 0; 0 1 0; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ zero(SMatrix{3,3,Float64,9}) atol=eps()

    F = @SMatrix [2 0 0; 0 2 0; 0 0 2]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [5.33203125 0.0 0.0; 0.0 5.33203125 0.0; 0.0 0.0 5.33203125]

    F = @SMatrix [1 0.5 0; 0 1 0; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [-0.03333333333333327 0.2 0.0;
                        0.21666666666666667 -0.033333333333333305 0.0;
                        0.0 0.0 -0.033333333333333305]

    F = @SMatrix [1 0 0; 0 1 0.5; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [-0.033333333333333305 0.0 0.0;
                        0.0 -0.03333333333333327 0.2;
                        0.0 0.21666666666666667 -0.033333333333333305]

    F = @SMatrix [1 0 0.5; 0 1 0; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [-0.03333333333333327 0.0 0.2;
                        0.0 -0.033333333333333305 0.0;
                        0.21666666666666667 0.0 -0.033333333333333305]

    F = @SMatrix [1 0 0; 0.5 1 0; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [-0.033333333333333305 0.21666666666666667 0.0;
                        0.2 -0.03333333333333327 0.0;
                        0.0 0.0 -0.033333333333333305]
end

@testitem "SaintVenantKirchhoff" begin
    using Peridynamics.StaticArrays

    model = SaintVenantKirchhoff()
    mat = CMaterial(; model)
    pos, vol = uniform_box(1.0, 1.0, 1.0, 0.5)
    body = Body(mat, pos, vol)
    material!(body, horizon=1, rho=1, E=1, nu=0.25)
    params = body.point_params[1]
    decomp = Peridynamics.PointDecomposition(body, 1)
    system = Peridynamics.get_system(body, decomp, 1)
    solver = VelocityVerlet(steps=2)
    storage = Peridynamics.CStorage(mat, solver, system)

    F = @SMatrix [1 0 0; 0 1 0; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ zero(SMatrix{3,3,Float64,9}) atol=eps()

    F = @SMatrix [2 0 0; 0 2 0; 0 0 2]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [6.0 0.0 0.0; 0.0 6.0 0.0; 0.0 0.0 6.0]

    F = @SMatrix [1 0.5 0; 0 1 0; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [0.15 0.275 0.0; 0.2 0.15 0.0; 0.0 0.0 0.05]

    F = @SMatrix [1 0 0; 0 1 0.5; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [0.05 0.0 0.0; 0.0 0.15 0.275; 0.0 0.2 0.15]

    F = @SMatrix [1 0 0.5; 0 1 0; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [0.15 0.0 0.275; 0.0 0.05 0.0; 0.2 0.0 0.15]

    F = @SMatrix [1 0 0; 0.5 1 0; 0 0 1]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test P ≈ @SMatrix [0.15 0.2 0.0; 0.275 0.15 0.0; 0.0 0.0 0.05]
end
