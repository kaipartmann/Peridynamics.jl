# The custom damage model tutorial is included rather than copied, so that it cannot rot: if
# the damage model interface changes and the tutorial is not updated with it, these tests
# fail. It defines `DelayedFailure`, a model with a parameter and a state of its own, whose
# bonds break only after they were overstretched for the delay `tau`.

@testmodule CustomDamageModelTutorial begin
    using Peridynamics
    using Peridynamics.LinearAlgebra

    const TUTORIAL = normpath(@__DIR__, "..", "..", "docs", "src", "literate",
                              "tutorial_custom_damage_model.jl")
    include(TUTORIAL)

    "A bar of `mat` pulled apart at both ends, with the fracture keywords of `kwargs`."
    function bar(mat; l=0.1, Δx=0.002, kwargs...)
        pos, vol = uniform_box(l, 0.1l, 0.1l, Δx)
        b = Body(mat, pos, vol)
        material!(b; horizon=3.015Δx, rho=2700, E=70e9, kwargs...)
        point_set!(x -> x < -0.4l, b, :left)
        point_set!(x -> x > 0.4l, b, :right)
        velocity_bc!(t -> -10.0, b, :left, :x)
        velocity_bc!(t -> 10.0, b, :right, :x)
        return b
    end
end

@testitem "tutorial damage model: the parameters" setup=[CustomDamageModelTutorial] begin
    T = CustomDamageModelTutorial
    dmgmodel = T.DelayedFailure()

    # the model registered its own keyword next to the standard fracture keywords
    @test Peridynamics.damage_param_kwargs(dmgmodel) == (:Gc, :epsilon_c, :tau)
    mat = BBMaterial(; dmgmodel)
    @test :tau in Peridynamics.all_material_kwargs(mat)

    # the parameters are read flat off the point parameters, and `Gc` is converted with the
    # default relation of the material, so nothing about the conversion had to be written
    body = T.bar(mat; Gc=100, tau=2e-6)
    params = only(body.point_params)
    @test params.τ == 2e-6
    @test params.εc ≈ sqrt(5 * 100 / (9 * params.K * params.δ))
    @test params.dmg_params isa T.DelayedFailureParameters

    # fracture is enabled by the default `has_fracture`, without the model defining it
    @test Peridynamics.has_fracture(mat, params)
    @test all(body.fail_permit)

    # the delay is required
    pos, vol = uniform_box(0.1, 0.01, 0.01, 0.002)
    b = Body(mat, pos, vol)
    @test_throws UndefKeywordError material!(b; horizon=0.007, rho=2700, E=70e9, Gc=100)
end

@testitem "tutorial damage model: the state" setup=[CustomDamageModelTutorial] begin
    T = CustomDamageModelTutorial
    mat = BBMaterial(; dmgmodel=T.DelayedFailure())

    # the model brought a state with the inherited bookkeeping and its own field, and the
    # storage of the material carries it concretely
    @test Peridynamics.damage_storage_type(T.DelayedFailure(),
                                          Peridynamics.system_type(mat)) <:
          Peridynamics.AbstractDamageState
    S = Peridynamics.storage_type(mat)
    @test isconcretetype(S)
    @test Peridynamics.has_damage_state(S)
    @test fieldtype(S, :dmg_state) ===
          T.DelayedFailureState{3,Float64,Vector{Float64},Vector{Int},Vector{Bool}}

    # one entry per bond of the chunk, starting undamaged
    body = T.bar(mat; Gc=100, tau=2e-6)
    chunk = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1).chunks[1]
    state = Peridynamics.damage_state(chunk.storage)
    @test length(state.bond_damage) == Peridynamics.get_n_bonds(chunk.system)
    @test all(iszero, state.bond_damage)
end

@testitem "tutorial damage model: a bond breaks only after the delay" setup=[CustomDamageModelTutorial] begin
    T = CustomDamageModelTutorial
    mat = BBMaterial(; dmgmodel=T.DelayedFailure())
    body = T.bar(mat; Gc=100, tau=2e-6)
    chunk = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1).chunks[1]
    (; storage, system, paramsetup) = chunk
    params = Peridynamics.get_params(paramsetup, 1)
    (; εc, τ) = params

    # a uniaxial stretch of two and a half times the critical stretch: the bonds along the
    # axis are at ε = 2.5 εc, the diagonal ones below that, the transverse ones at zero
    storage.position .= system.position
    storage.position[1, :] .*= 1 + 2.5 * εc
    Δt = τ / 4
    i = argmin(abs.(@view system.position[1, 1:Peridynamics.get_n_loc_points(system)]))
    bond_ids = Peridynamics.each_bond_idx(system, i)
    state = Peridynamics.damage_state(storage)

    # one step: damage has grown, but no bond has reached one. The package fills the current
    # bond lengths before it calls the criterion, so a direct call has to do the same, which
    # is what `update_bond_lengths!` is public for.
    storage.n_active_bonds[i] = 0
    Peridynamics.update_bond_lengths!(storage, system, i)
    Peridynamics.calc_failure!(storage, system, mat, T.DelayedFailure(), paramsetup, 0.0, Δt, i)
    @test all(storage.bond_active[bond_ids])
    @test storage.n_active_bonds[i] == system.n_neighbors[i]
    @test maximum(state.bond_damage[bond_ids]) ≈ 1.5 * Δt / τ
    # a bond that is not overstretched accumulates nothing
    for bond_id in bond_ids
        j = Peridynamics.get_neighbor(system, bond_id)
        L = Peridynamics.reference_bond_length(system, bond_id)
        ε = (Peridynamics.LinearAlgebra.norm(Peridynamics.get_vector_diff(storage.position, i, j, Peridynamics.dims(system))) - L) / L
        ε > εc || @test state.bond_damage[bond_id] == 0
    end

    # after the delay the bonds along the axis are gone, the others are still there
    for _ in 1:2
        storage.n_active_bonds[i] = 0
        Peridynamics.update_bond_lengths!(storage, system, i)
        Peridynamics.calc_failure!(storage, system, mat, T.DelayedFailure(), paramsetup, 0.0,
                                   Δt, i)
    end
    active = storage.bond_active[bond_ids]
    @test !all(active)
    @test any(active)
    @test all(state.bond_damage[bond_ids][.!active] .>= 1)
    @test all(state.bond_damage[bond_ids][active] .< 1)
    @test storage.n_active_bonds[i] == count(active)

    # a point that may not fail keeps every bond, however long it is overstretched
    body_nf = T.bar(mat; Gc=100, tau=2e-6)
    no_failure!(body_nf)
    chunk_nf = Peridynamics.threads_data_handler(body_nf, VelocityVerlet(steps=1), 1).chunks[1]
    chunk_nf.storage.position .= chunk_nf.system.position
    chunk_nf.storage.position[1, :] .*= 1 + 2.5 * εc
    for _ in 1:8
        chunk_nf.storage.n_active_bonds[i] = 0
        Peridynamics.update_bond_lengths!(chunk_nf.storage, chunk_nf.system, i)
        Peridynamics.calc_failure!(chunk_nf.storage, chunk_nf.system, mat, T.DelayedFailure(),
                                   chunk_nf.paramsetup, 0.0, Δt, i)
    end
    @test all(chunk_nf.storage.bond_active[bond_ids])
    @test all(iszero, Peridynamics.damage_state(chunk_nf.storage).bond_damage[bond_ids])
end

@testitem "tutorial damage model: runs, breaks, and exports its state" tags=[:simulation] setup=[Fixtures, CustomDamageModelTutorial] begin
    T = CustomDamageModelTutorial
    delayed = submit(Job(T.bar(BBMaterial(; dmgmodel=T.DelayedFailure()); Gc=100, tau=2e-6),
                         VelocityVerlet(steps=400); path=joinpath(mktempdir(), "delayed"));
                     quiet=true)

    u = Fixtures.whole_body(delayed, :displacement)
    @test all(isfinite, u)

    # the model reaches the simulation: a bar pulled this hard breaks
    @test maximum(Fixtures.whole_body(delayed, :damage)) > 0

    # the accumulated damage reaches the export as one value per local point, and every
    # broken bond reached one on the way
    storage, system = delayed.chunks[1].storage, delayed.chunks[1].system
    out = Peridynamics.export_field(Val(:bond_damage), delayed.chunks[1].mat, system, storage,
                                    delayed.chunks[1].paramsetup, 0.0)
    @test length(out) == Peridynamics.get_n_loc_points(system)
    @test all(isfinite, out)
    @test Peridynamics.custom_field(typeof(storage), :bond_damage)
    state = Peridynamics.damage_state(storage)
    @test all(state.bond_damage[.!storage.bond_active] .>= 1)
end
