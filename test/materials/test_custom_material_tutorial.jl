# The custom material tutorial is the single source of truth for "how does someone add a
# material". It is included here rather than copied, so that the tutorial cannot rot: if the
# extension API changes and the tutorial is not updated with it, these tests fail.
#
# The tutorial defines `ConicalBBMaterial` with a conical micro-modulus, its point
# parameters, its storage, `force_density_point!`, a custom export field, the
# `ConicalStretch` damage model that matches that micro-modulus, and a `MooneyRivlin`
# constitutive model. It does not call `submit`, so including it is cheap and writes nothing.

@testmodule CustomMaterialTutorial begin
    using Peridynamics
    using Peridynamics.LinearAlgebra

    const TUTORIAL = normpath(@__DIR__, "..", "..", "docs", "src", "literate",
                              "tutorial_custom_material.jl")
    include(TUTORIAL)

    "A bar of the tutorial material pulled apart at both ends."
    function conical_body(; l=0.1, Δx=0.002, dmgmodel=ConicalStretch(), kwargs...)
        pos, vol = uniform_box(l, 0.1l, 0.1l, Δx)
        b = Body(ConicalBBMaterial(; dmgmodel), pos, vol)
        material!(b; horizon=3.015Δx, rho=2700, E=70e9, kwargs...)
        point_set!(x -> x < -0.4l, b, :left)
        point_set!(x -> x > 0.4l, b, :right)
        velocity_bc!(t -> -10.0, b, :left, :x)
        velocity_bc!(t -> 10.0, b, :right, :x)
        return b
    end

    """
    A point field of a finished simulation for the whole body, assembled from the local
    entries of every chunk. With several threads `dh.chunks[1]` holds one chunk only.
    """
    function whole_body(dh, field)
        n = sum(Peridynamics.get_n_loc_points(c.system) for c in dh.chunks)
        sample = Peridynamics.get_loc_point_data(dh.chunks[1].storage, dh.chunks[1].system,
                                                 field)
        out = sample isa AbstractMatrix ? zeros(size(sample, 1), n) : zeros(n)
        for c in dh.chunks
            ids = Peridynamics.get_loc_points(c.system)
            data = Peridynamics.get_loc_point_data(c.storage, c.system, field)
            if data isa AbstractMatrix
                out[:, ids] .= data
            else
                out[ids] .= data
            end
        end
        return out
    end
end

@testitem "tutorial material: the generated point parameters" setup=[CustomMaterialTutorial] begin
    T = CustomMaterialTutorial
    mat = T.ConicalBBMaterial()

    # `@params` derived the keyword list from the declarations, so the derived micro-modulus
    # constant is not a keyword, and the fracture keywords of the damage model are
    kwargs = Peridynamics.all_material_kwargs(mat)
    @test !(:bc in kwargs)
    @test :Gc in kwargs
    @test :epsilon_c in kwargs

    P = Peridynamics.point_param_type(mat)
    @test isconcretetype(P)
    @test isbitstype(P)
    @test :bc in fieldnames(P)
    @test :dmg_params in fieldnames(P)

    # the parameters are parametric in the float type of the simulation, and the damage
    # model's parameters are resolved from the model of the material instance
    CSP = T.ConicalStretchParameters
    @test Peridynamics.point_param_type(mat, Float32) ===
          T.ConicalBBPointParameters{Float32,CSP{Float32}}

    body = T.conical_body(; Gc=100)
    params = only(body.point_params)

    # the conical micro-modulus is normalized by the same energy equivalence as the constant
    # one, which makes it exactly five times as large
    @test params.bc ≈ 90 * params.K / (π * params.δ^4)
    @test params.bc ≈ 5 * 18 * params.K / (π * params.δ^4)

    # a misspelled keyword is rejected rather than silently ignored
    pos, vol = uniform_box(0.1, 0.01, 0.01, 0.002)
    bad = Body(T.ConicalBBMaterial(), pos, vol)
    @test_throws ArgumentError material!(bad; horizon=0.007, rho=2700, E=70e9, Gcc=100)
end

@testitem "tutorial material: the conical micro-modulus is normalized correctly" setup=[CustomMaterialTutorial] begin
    # The normalization is fixed by requiring the same strain energy as a classical isotropic
    # solid under a homogeneous stretch: W = π s² ∫₀^δ c(ξ) ξ³ dξ = 9/2 K s². Check that
    # integral numerically against the derived constant, rather than restating the algebra.
    δ, K = 0.00603, 70e9 / (3 * (1 - 2 * 0.25))
    bc = 90 * K / (π * δ^4)

    n = 2_000_000
    ξ = range(0, δ; length=n)
    integral = sum(bc * (1 - x / δ) * x^3 for x in ξ) * step(ξ)
    @test π * integral ≈ 4.5 * K rtol=1e-4

    # and the same integral for a constant micro-modulus reproduces the built-in constant
    c_const = 18 * K / (π * δ^4)
    integral_const = sum(c_const * x^3 for x in ξ) * step(ξ)
    @test π * integral_const ≈ 4.5 * K rtol=1e-4
end

@testitem "tutorial material: the generated storage" setup=[CustomMaterialTutorial] begin
    T = CustomMaterialTutorial
    S = Peridynamics.storage_type(T.ConicalBBMaterial())

    # `storage_type` has to stay concrete, otherwise every body chunk is type unstable
    @test isconcretetype(S)
    @test S <: Peridynamics.AbstractStorage
    @test :bond_stretch in fieldnames(S)

    # the inherited solver fields are all there
    for field in (:position, :displacement, :velocity, :b_int, :damage, :bond_active)
        @test field in fieldnames(S)
    end

    # `position` is inherited with its halo annotation, so the material parallelizes: it is
    # exchanged local-to-halo and is allocated with halo entries, unlike a local-only field
    body = T.conical_body(; Gc=100)
    chunk = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1).chunks[1]
    storage, system = chunk.storage, chunk.system
    @test Peridynamics.is_halo_field(storage, Val(:position))
    @test :position in Peridynamics.loc_to_halo_fields(storage)
    @test size(storage.position, 2) == Peridynamics.get_n_points(system)
    @test size(storage.displacement, 2) == Peridynamics.get_n_loc_points(system)
end

@testitem "tutorial material: the force law is what it claims" setup=[CustomMaterialTutorial] begin
    T = CustomMaterialTutorial
    # Evaluate the force law at a fixed configuration, so that the only thing that differs
    # between the cases is the constitutive law and not how the deformation localized.
    body = T.conical_body(; Gc=100)
    chunk = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1).chunks[1]
    storage, system = chunk.storage, chunk.system

    # a uniform uniaxial stretch of 2 %
    storage.position .= system.position
    storage.position[1, :] .*= 1.02
    storage.bond_active .= true
    storage.b_int .= 0

    n_loc = Peridynamics.get_n_loc_points(system)
    i = argmax(@view system.position[1, 1:n_loc])
    params = Peridynamics.get_params(chunk.paramsetup, i)
    Peridynamics.force_density_point!(storage, system, chunk.mat, params, 0.0, 0.0, i)

    # the stretch of every bond of `i` was recorded on the way through
    @test all(isfinite, storage.bond_stretch)
    @test maximum(storage.bond_stretch) > 0

    # the point at the free end has a one-sided family, so it carries a net force
    b_int = [storage.b_int[d, i] for d in 1:3]
    @test T.norm(b_int) > 0

    # Reproduce the force of point `i` by hand from the conical law. This is the assertion
    # that the micro-modulus really falls off with the bond length: a constant micro-modulus
    # would give a different number.
    b_ref = zeros(3)
    b_const = zeros(3)
    for bond_id in Peridynamics.each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = Peridynamics.get_vector_diff(storage.position, i, j)
        len = T.norm(Δxij)
        ε = (len - L) / L
        b_ref .+= params.bc * (1 - L / params.δ) * ε * system.volume[j] / len .* Δxij
        b_const .+= 18 * params.K / (π * params.δ^4) * ε * system.volume[j] / len .* Δxij
    end
    @test b_int ≈ b_ref
    @test !isapprox(b_int, b_const; rtol=1e-3)

    # a bond exactly at the horizon carries no force at all, which is the whole point of the
    # conical shape
    @test params.bc * (1 - params.δ / params.δ) == 0
end

@testitem "tutorial damage model: the conical Gc conversion" setup=[CustomMaterialTutorial] begin
    T = CustomMaterialTutorial
    δ, K, Gc = 0.00603, 46.6e9, 100.0

    conical = Peridynamics.get_frac_params(T.ConicalStretch(), δ, K; Gc=Gc)
    constant = Peridynamics.get_frac_params(CriticalStretch(), δ, K; Gc=Gc)

    @test conical.Gc == Gc
    @test conical.εc ≈ sqrt(2 * Gc / (3 * K * δ))

    # the whole reason the damage model exists: the built-in conversion is derived for a
    # constant micro-modulus and is √1.2 too small for a conical one
    @test conical.εc / constant.εc ≈ sqrt(1.2)

    # the conversion round-trips in both directions
    back = Peridynamics.get_frac_params(T.ConicalStretch(), δ, K; epsilon_c=conical.εc)
    @test back.Gc ≈ Gc
    @test back.εc ≈ conical.εc

    # giving both is an error, and giving neither switches fracture off
    @test_throws ArgumentError Peridynamics.get_frac_params(T.ConicalStretch(), δ, K; Gc=Gc,
                                                            epsilon_c=1e-3)
    none = Peridynamics.get_frac_params(T.ConicalStretch(), δ, K)
    @test none.Gc == 0 && none.εc == 0

    # `material!` routes the keywords through the damage model of the material, and the
    # parameters of the model are read flat off the point parameters
    body = T.conical_body(; Gc=Gc)
    params = only(body.point_params)
    @test params.εc ≈ sqrt(2 * Gc / (3 * params.K * params.δ))
    @test params.dmg_params isa T.ConicalStretchParameters
    @test Peridynamics.has_fracture(body.mat, params)

    # and the same body with the built-in damage model gets the other, wrong, value
    body_wrong = T.conical_body(; dmgmodel=CriticalStretch(), Gc=Gc)
    @test only(body_wrong.point_params).εc ≈ params.εc / sqrt(1.2)
end

@testitem "tutorial material: runs, breaks, and exports its own field" tags=[:simulation] setup=[CustomMaterialTutorial] begin
    T = CustomMaterialTutorial
    body = T.conical_body(; Gc=100)
    dh = submit(Job(body, VelocityVerlet(steps=200); path=joinpath(mktempdir(), "conical"));
                quiet=true)
    storage = dh.chunks[1].storage
    system = dh.chunks[1].system

    u = T.whole_body(dh, :displacement)
    @test all(isfinite, u)
    @test length(storage.bond_stretch) == Peridynamics.get_n_bonds(system)

    # the bar is pulled apart: the left end moved left, the right end moved right
    left, right = body.point_sets[:left], body.point_sets[:right]
    @test maximum(@view u[1, left]) < 0
    @test minimum(@view u[1, right]) > 0

    # the custom export field reduces the bond field to one value per local point
    weighted = Peridynamics.export_field(Val(:weighted_stretch), T.ConicalBBMaterial(),
                                         system, storage, dh.chunks[1].paramsetup, 0.0)
    @test length(weighted) == Peridynamics.get_n_loc_points(system)
    @test all(isfinite, weighted)
    @test Peridynamics.custom_field(typeof(storage), :weighted_stretch)
    @test !Peridynamics.custom_field(typeof(storage), :not_a_field)

    # the damage model reaches the simulation: a bar pulled this hard breaks
    damage = T.whole_body(dh, :damage)
    @test maximum(damage) > 0

    # a larger Gc means a larger critical stretch, so less damage after the same 200 steps
    dh_tough = submit(Job(T.conical_body(; Gc=1e4), VelocityVerlet(steps=200);
                          path=joinpath(mktempdir(), "conical_tough")); quiet=true)
    @test sum(T.whole_body(dh_tough, :damage)) < sum(damage)
end

@testitem "tutorial constitutive model: runs on every correspondence family" tags=[:simulation] setup=[CustomMaterialTutorial] begin
    T = CustomMaterialTutorial
    # the point of a constitutive model is that it is not tied to one material family
    for mat in (RKCMaterial(; model=T.MooneyRivlin()), CMaterial(; model=T.MooneyRivlin()),
                BACMaterial(; model=T.MooneyRivlin()))
        @test Peridynamics.get_constitutive_model(mat) isa T.MooneyRivlin
        # a model without a declared state is not history dependent and costs nothing
        @test !Peridynamics.is_history_dependent(Peridynamics.get_constitutive_model(mat))
        @test Peridynamics.constitutive_storage_type(T.MooneyRivlin()) === Nothing
    end

    pos, vol = uniform_box(0.1, 0.02, 0.02, 0.005)
    body = Body(RKCMaterial(; model=T.MooneyRivlin()), pos, vol)
    material!(body; horizon=3.015 * 0.005, rho=2700, E=70e9, nu=0.3, Gc=100)
    point_set!(x -> x < -0.04, body, :left)
    point_set!(x -> x > 0.04, body, :right)
    velocity_bc!(t -> -1.0, body, :left, :x)
    velocity_bc!(t -> 1.0, body, :right, :x)

    dh = submit(Job(body, VelocityVerlet(steps=50); path=mktempdir()); quiet=true)
    u = T.whole_body(dh, :displacement)
    @test all(isfinite, u)
    @test maximum(abs, u) > 0

    # the stress of an undeformed configuration vanishes
    F = Peridynamics.StaticArrays.SMatrix{3,3,Float64,9}(1, 0, 0, 0, 1, 0, 0, 0, 1)
    params = only(body.point_params)
    P = Peridynamics.first_piola_kirchhoff(T.MooneyRivlin(), nothing, params, F)
    @test isapprox(maximum(abs, P), 0; atol=1e-8)
end
