# The custom constitutive model tutorial is included rather than copied, so that it cannot
# rot. It defines a Mooney-Rivlin model with parameters of its own and a J2 plasticity model
# with linear isotropic hardening that carries its plastic strain as a state.

@testmodule CustomConstitutiveModelTutorial begin
    using Peridynamics
    using Peridynamics.LinearAlgebra
    using Peridynamics.StaticArrays

    const TUTORIAL = normpath(@__DIR__, "..", "..", "docs", "src", "literate",
                              "tutorial_custom_constitutive_model.jl")
    include(TUTORIAL)

    "A steel bar of `mat` pulled apart at both ends with `v`."
    function steel_bar(mat; l=0.1, Δx=0.004, v=5.0, kwargs...)
        pos, vol = uniform_box(l, 0.2l, 0.2l, Δx)
        b = Body(mat, pos, vol)
        material!(b; horizon=3.015Δx, rho=7850, E=210e9, nu=0.3, kwargs...)
        point_set!(x -> x < -0.4l, b, :left)
        point_set!(x -> x > 0.4l, b, :right)
        velocity_bc!(t -> -v, b, :left, :x)
        velocity_bc!(t -> v, b, :right, :x)
        return b
    end

    "The storage and the parameters of a one-chunk body of `mat`."
    function chunk_of(body)
        chunk = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1).chunks[1]
        return chunk.storage, only(body.point_params)
    end
end

@testitem "tutorial constitutive model: Mooney-Rivlin parameters and stress" setup=[CustomConstitutiveModelTutorial] begin
    T = CustomConstitutiveModelTutorial
    # the model registered its constants as keywords of `material!`
    @test Peridynamics.constitutive_param_kwargs(T.MooneyRivlin()) == (:C10, :C01)
    params = only(T.rubber.point_params)
    @test params.C10 == 0.3e6
    @test params.C01 == 0.1e6

    # `C01` is optional, and the model is not tied to one material family
    pos, vol = uniform_box(0.1, 0.02, 0.02, 0.004)
    for mat in (RKCMaterial(; model=T.MooneyRivlin()), CMaterial(; model=T.MooneyRivlin()),
                BACMaterial(; model=T.MooneyRivlin()))
        @test Peridynamics.get_constitutive_model(mat) isa T.MooneyRivlin
        @test !Peridynamics.is_history_dependent(T.MooneyRivlin())
        body = Body(mat, pos, vol)
        material!(body; horizon=0.012, rho=1100, E=2.4e6, nu=0.49, C10=0.3e6)
        @test only(body.point_params).C01 == 0
    end

    # the stress of an undeformed configuration vanishes, and the Cauchy stress of a
    # deformed one is symmetric
    F0 = T.SMatrix{3,3,Float64,9}(1, 0, 0, 0, 1, 0, 0, 0, 1)
    P0 = Peridynamics.first_piola_kirchhoff(T.MooneyRivlin(), nothing, params, F0)
    @test isapprox(maximum(abs, P0), 0; atol=1e-8)
    F = T.SMatrix{3,3,Float64,9}(1.1, 0.02, 0, 0, 0.97, 0, 0.01, 0, 0.98)
    P = Peridynamics.first_piola_kirchhoff(T.MooneyRivlin(), nothing, params, F)
    σ = P * F' / T.det(F)
    @test σ ≈ σ'
    @test maximum(abs, σ) > 0
end

@testitem "tutorial constitutive model: J2 plasticity below and above the yield stress" setup=[CustomConstitutiveModelTutorial] begin
    T = CustomConstitutiveModelTutorial
    body = T.steel_bar(RKCMaterial(; model=T.J2Plasticity()); sigma_y=250e6, H=1e9)
    storage, params = T.chunk_of(body)
    @test Peridynamics.constitutive_param_kwargs(T.J2Plasticity()) == (:sigma_y, :H)
    @test params.sigma_y == 250e6
    @test Peridynamics.is_history_dependent(T.J2Plasticity())
    state = Peridynamics.constitutive_state(storage)
    idx = 1

    # below the yield stress the model is the elastic model in logarithmic strain space and
    # the state stays untouched
    e = 0.5 * params.sigma_y / params.E
    F = T.SMatrix{3,3,Float64,9}(1 + e, 0, 0, 0, 1 - 0.3e, 0, 0, 0, 1 - 0.3e)
    P = Peridynamics.first_piola_kirchhoff(T.J2Plasticity(), storage, params, F, idx, 1e-6)
    ε, Uinv = Peridynamics.hencky_and_invstretch(F' * F)
    τ = params.λ * T.tr(ε) * T.I + 2 * params.μ * ε
    @test P ≈ F * (Uinv * τ * Uinv)
    @test state.bond_eqps[idx] == 0
    @test all(iszero, Peridynamics.get_sym_tensor(state.bond_plastic_strain, idx))

    # above the yield stress the stress is returned onto the yield surface and the plastic
    # strain accumulates, isochoric
    e = 4 * params.sigma_y / params.E
    F = T.SMatrix{3,3,Float64,9}(1 + e, 0, 0, 0, 1 - 0.3e, 0, 0, 0, 1 - 0.3e)
    P = Peridynamics.first_piola_kirchhoff(T.J2Plasticity(), storage, params, F, idx, 1e-6)
    eqps = state.bond_eqps[idx]
    @test eqps > 0
    εp = Peridynamics.get_sym_tensor(state.bond_plastic_strain, idx)
    @test isapprox(T.tr(εp), 0; atol=1e-12)
    # the Kirchhoff stress that belongs to the returned P lies on the yield surface
    ε, Uinv = Peridynamics.hencky_and_invstretch(F' * F)
    U = inv(Uinv)
    τ = U * (F \ P) * U
    s = τ - T.tr(τ) / 3 * T.I
    q = sqrt(1.5) * T.norm(s)
    @test q ≈ params.sigma_y + params.H * eqps
    # calling again with the same deformation is elastic up to round-off: the state stays
    Peridynamics.first_piola_kirchhoff(T.J2Plasticity(), storage, params, F, idx, 1e-6)
    @test state.bond_eqps[idx] ≈ eqps

    # the strain energy density is the elastic part and leaves the state alone
    eqps = state.bond_eqps[idx]
    Ψ = Peridynamics.strain_energy_density(T.J2Plasticity(), storage, params, F, idx)
    @test Ψ > 0
    @test state.bond_eqps[idx] == eqps
end

@testitem "tutorial constitutive model: the plastic bar" tags=[:simulation] setup=[Fixtures, CustomConstitutiveModelTutorial] begin
    T = CustomConstitutiveModelTutorial
    body = T.steel_bar(RKCMaterial(; model=T.J2Plasticity()); v=50.0, sigma_y=250e6, H=1e9)
    dh = submit(Job(body, VelocityVerlet(steps=200); path=mktempdir()); quiet=true)
    u = Fixtures.whole_body(dh, :displacement)
    @test all(isfinite, u)
    @test maximum(abs, u) > 0

    # the bar yielded, and the plastic strain reaches the export as one value per point
    eqps_max = map(dh.chunks) do c
        eqps = Peridynamics.export_field(Val(:equivalent_plastic_strain), c.mat, c.system,
                                         c.storage, c.paramsetup, 0.0)
        @test length(eqps) == Peridynamics.get_n_loc_points(c.system)
        @test all(isfinite, eqps)
        return maximum(eqps)
    end
    @test maximum(eqps_max) > 0

    # a history-dependent model is rejected by a solver that evaluates the force density
    # several times per step
    @test_throws Peridynamics.HistoryDependenceError Job(body, NewtonKrylov(steps=1);
                                                          path=mktempdir())
end
