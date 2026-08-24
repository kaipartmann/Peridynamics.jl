@testitem "required_fields_fracture" begin
    @test Peridynamics.required_fields_fracture(Peridynamics.AbstractMaterial) == ()

    rff_bond_system = (:damage, :n_active_bonds, :bond_active)
    @test Peridynamics.required_fields_fracture(BBMaterial) === rff_bond_system

    rff_interaction_system = (:damage, :n_active_one_nis, :one_ni_active)
    @test Peridynamics.required_fields_fracture(CKIMaterial) === rff_interaction_system
end

@testitem "get_frac_params: the fracture keywords of CriticalStretch and the Dict bridge" begin
    import Peridynamics: get_frac_params

    # `Gc` and `epsilon_c` are converted into each other, one of them is enough
    frac = get_frac_params(CriticalStretch(), 2.0, 3.0; Gc=1.0)
    @test frac.Gc == 1.0
    @test frac.εc ≈ sqrt(5.0 * 1.0 / (9.0 * 3.0 * 2.0))
    frac = get_frac_params(CriticalStretch(), 2.0, 3.0; epsilon_c=0.1)
    @test frac.εc == 0.1
    @test frac.Gc ≈ 9.0 / 5.0 * 3.0 * 2.0 * 0.1^2
    # both is an error, none means no fracture, unknown keywords are ignored
    @test_throws ArgumentError get_frac_params(CriticalStretch(), 2.0, 3.0; Gc=1.0,
                                               epsilon_c=0.1)
    @test get_frac_params(CriticalStretch(), 2.0, 3.0) == (; Gc=0.0, εc=0.0)
    frac_gc = get_frac_params(CriticalStretch(), 2.0, 3.0; Gc=1.0)
    @test get_frac_params(CriticalStretch(), 2.0, 3.0; Gc=1.0, some_other_keyword=1.0) == frac_gc
    # a keyword that was not given arrives as `nothing`
    @test get_frac_params(CriticalStretch(), 2.0, 3.0; Gc=1.0, epsilon_c=nothing) ==
          get_frac_params(CriticalStretch(), 2.0, 3.0; Gc=1.0)

    # a damage model without parameters answers with an empty named tuple
    struct ParameterlessDamage <: Peridynamics.AbstractDamageModel end
    @test get_frac_params(ParameterlessDamage(), 1.0, 1.0; Gc=1.0) == (;)

    # the point parameters of `@params` still read the keywords from the `Dict` of
    # `material!`; the bridge passes every fracture keyword, missing ones as `nothing`
    p = Dict{Symbol,Any}(:Gc => 1.0, :horizon => 1.0)
    @test get_frac_params(ParameterlessDamage(), p, 1.0, 1.0) == (;)
    @test get_frac_params(CriticalStretch(), p, 2.0, 3.0) ==
          get_frac_params(CriticalStretch(), 2.0, 3.0; Gc=1.0)
    @test get_frac_params(CriticalStretch(), Dict{Symbol,Any}(), 2.0, 3.0) == (; Gc=0.0, εc=0.0)
end

@testitem "damage model hooks: the defaults of a model that deletes bonds" begin
    import Peridynamics: kinematic_weight, safe_degradation, degrade_bond_stress, log_dmgmodel,
                         damage_storage_type, get_dmg_storage, init_damage_state, damage_state,
                         has_damage_state, storage_type, req_storage_fields
    using Peridynamics.StaticArrays

    pos, vol = uniform_box(1, 1, 1, 0.5)
    body = Body(BBMaterial(), pos, vol)
    material!(body; horizon=1.5, rho=1, E=1, Gc=1)
    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    (; storage, system) = dh.chunks[1]
    dmg = CriticalStretch()

    # no degradation: the weight is one and the stress passes through unchanged
    @test kinematic_weight(dmg, storage, 1) === 1.0
    @test safe_degradation(dmg, storage, 1) === 1.0
    P = SMatrix{3,3,Float64,9}(1:9)
    @test degrade_bond_stress(dmg, storage, 1, P) === P

    # no state: the storage carries `nothing` and the contract asks for nothing
    @test damage_storage_type(dmg) === Nothing
    @test damage_storage_type(dmg, Float32) === Nothing
    @test isnothing(get_dmg_storage(dmg, VelocityVerlet(steps=1), system))
    @test isnothing(init_damage_state(BBMaterial(), VelocityVerlet(steps=1), system))
    @test isnothing(damage_state(storage))
    @test !has_damage_state(storage_type(BBMaterial()))
    @test req_storage_fields(BBMaterial(), dmg) == ()
    @test req_storage_fields(BBMaterial(), nothing) == ()

    # the log names the model
    @test contains(log_dmgmodel(dmg; indentation=2), "damage model type")
    @test contains(log_dmgmodel(dmg; indentation=2), "CriticalStretch")
    @test contains(Peridynamics.log_material_property(Val(:dmgmodel), BBMaterial(); indentation=2),
                   "CriticalStretch")
    @test contains(Peridynamics.log_material_property(Val(:dmgmodel), CKIMaterial(); indentation=2),
                   "CriticalStretch")
    @test contains(Peridynamics.log_material_property(Val(:dmgmodel), BACMaterial(); indentation=2),
                   "CriticalStretch")
end

@testsnippet FatigueModel begin
    using Peridynamics: damage_state, get_params, each_bond_idx, get_vector_diff
    using Peridynamics.LinearAlgebra

    # A damage model with a state of its own: a bond fails only after its stretch exceeded
    # the critical stretch `n_cycles` times. The count lives in the state of the model, so
    # the model works with every material whose storage declares `dmg_state::DamageState`.
    struct FatigueDamage <: Peridynamics.AbstractDamageModel
        n_cycles::Int
    end
    FatigueDamage() = FatigueDamage(3)

    Peridynamics.@dmg_storage FatigueDamage struct FatigueState
        bond_exceedances::BondScalar{Int}
        bond_weight::BondScalar = 1.0
    end

    function Peridynamics.get_frac_params(::FatigueDamage, δ, K; kwargs...)
        return Peridynamics.get_frac_params(CriticalStretch(), δ, K; kwargs...)
    end
    function Peridynamics.has_fracture(::FatigueDamage, params)
        return Peridynamics.has_fracture(CriticalStretch(), params)
    end

    function Peridynamics.calc_failure!(storage, system, mat, dmg::FatigueDamage, paramsetup,
                                        i)
        (; εc) = get_params(paramsetup, i)
        (; bond_exceedances) = damage_state(storage)
        for bond_id in each_bond_idx(system, i)
            bond = system.bonds[bond_id]
            j, L = bond.neighbor, bond.length
            ε = (norm(get_vector_diff(storage.position, i, j)) - L) / L
            if ε > εc && bond.fail_permit
                bond_exceedances[bond_id] += 1
                if bond_exceedances[bond_id] >= dmg.n_cycles
                    storage.bond_active[bond_id] = false
                end
            end
            storage.n_active_bonds[i] += storage.bond_active[bond_id]
        end
        return nothing
    end

    # a bar pulled apart at both ends, relaxed quasi-statically; the total stretch is
    # `v * steps`, so the bonds near the ends exceed `epsilon_c` after a few steps
    function stretched_body(mat; Δx=0.1, v=2e-4)
        pos, vol = uniform_box(1.0, 0.3, 0.3, Δx)
        body = Body(mat, pos, vol)
        material!(body; horizon=3.015Δx, rho=8e-6, E=2.1e5, nu=0.25, epsilon_c=0.005)
        point_set!(x -> x < -0.4, body, :left)
        point_set!(x -> x > 0.4, body, :right)
        velocity_bc!(t -> -v, body, :left, :x)
        velocity_bc!(t -> v, body, :right, :x)
        return body
    end
end

@testitem "@dmg_storage: the generated state of a damage model" setup=[FatigueModel] begin
    using Peridynamics: damage_storage_type, get_dmg_storage, storage_fields_expr, get_n_bonds

    dmg = FatigueDamage()

    # the state is a parametric struct, exactly like a storage
    @test FatigueState isa UnionAll
    @test FatigueState <: Peridynamics.AbstractDamageState
    @test fieldnames(FatigueState) == (:bond_exceedances, :bond_weight)
    @test damage_storage_type(dmg) === FatigueState{Float64,Vector{Int},Vector{Float64}}
    @test damage_storage_type(dmg, Float32) === FatigueState{Float32,Vector{Int},Vector{Float32}}
    @test isconcretetype(damage_storage_type(dmg))
    @test [d.name for d in storage_fields_expr(FatigueState)] == [:bond_exceedances, :bond_weight]

    # the state is allocated per chunk from the system, with the initial values of the
    # declarations
    pos, vol = uniform_box(1, 1, 1, 0.5)
    body = Body(BBMaterial(), pos, vol)
    material!(body; horizon=1.5, rho=1, E=1, Gc=1)
    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    system = dh.chunks[1].system
    state = get_dmg_storage(dmg, VelocityVerlet(steps=1), system)
    @test state isa FatigueState{Float64,Vector{Int},Vector{Float64}}
    @test length(state.bond_exceedances) == get_n_bonds(system)
    @test all(iszero, state.bond_exceedances)
    @test all(isone, state.bond_weight)
    @test Peridynamics.Adapt.adapt(Array, state) isa FatigueState
end

@testitem "@dmg_storage: the macro input checks" setup=[FatigueModel] begin
    # halo annotations are not supported: the state is chunk-local
    @test_throws LoadError @eval Peridynamics.@dmg_storage FatigueDamage struct BadHaloState
        @lth some_field::PointScalar
    end

    # every field needs a field shape, there is no `init_field` hook for a model
    @test_throws LoadError @eval Peridynamics.@dmg_storage FatigueDamage struct BadTypeState
        some_field::Vector{Float64}
    end

    # a model cannot carry the state of another model
    @test_throws LoadError @eval Peridynamics.@dmg_storage FatigueDamage struct BadNestedState
        nested::DamageState
    end

    # only a damage model can declare a state
    struct NotADamageModel end
    @test_throws ArgumentError @eval Peridynamics.@dmg_storage NotADamageModel struct NotAState
        a::BondScalar
    end
    @test_throws ArgumentError Peridynamics.typecheck_damage_model(NotADamageModel)
    # a value instead of a type is rejected as well
    @test_throws ArgumentError Peridynamics.typecheck_damage_model(FatigueDamage())
end

@testitem "DamageState: a storage carries the state of its damage model" setup=[FatigueModel] begin
    using Peridynamics: storage_type, has_damage_state, damage_state, check_storage_contract,
                        req_storage_fields, StorageContractError

    # the storage stays concrete whichever model is used, and the model fills the parameter
    for (mat, DMS) in ((RKCMaterial(), Nothing),
                       (RKCMaterial(; dmgmodel=FatigueDamage()),
                        FatigueState{Float64,Vector{Int},Vector{Float64}}))
        S = storage_type(mat)
        @test isconcretetype(S)
        @test has_damage_state(S)
        @test fieldtype(S, :dmg_state) === DMS
    end

    # a stateful model requires a storage that carries its state, a stateless one does not
    @test req_storage_fields(RKCMaterial(), FatigueDamage()) == (:dmg_state,)
    @test req_storage_fields(RKCMaterial(), CriticalStretch()) == ()
    @test isnothing(check_storage_contract(RKCMaterial(; dmgmodel=FatigueDamage()),
                                           VelocityVerlet(steps=1)))
    @test !has_damage_state(storage_type(BBMaterial()))
    err = try
        check_storage_contract(BBMaterial(; dmgmodel=FatigueDamage()), VelocityVerlet(steps=1))
    catch e
        e
    end
    @test err isa StorageContractError
    @test contains(sprint(showerror, err), "dmg_state")
    @test contains(sprint(showerror, err), "FatigueDamage")

    # the chunk of a material with the state has it, filled from the declarations
    body = stretched_body(RKCMaterial(; dmgmodel=FatigueDamage()))
    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    state = damage_state(dh.chunks[1].storage)
    @test state isa FatigueState
    @test length(state.bond_exceedances) == Peridynamics.get_n_bonds(dh.chunks[1].system)

    # a storage may declare the state only once, and never annotated for halo exchange
    @test_throws LoadError @eval Peridynamics.@storage RKCMaterial struct TwiceDmgState
        @inherit Peridynamics.VelocityVerletFields
        a::DamageState
        b::DamageState
    end
    @test_throws LoadError @eval Peridynamics.@storage RKCMaterial struct HaloDmgState
        @inherit Peridynamics.VelocityVerletFields
        @htl a::DamageState
    end
end

@testitem "stateful damage model: a simulation with a model that brings its own state" tags=[:simulation] setup=[FatigueModel] begin
    using Peridynamics: damage_state

    # every bond has to exceed the critical stretch three times before it fails, so the
    # failure happens later than with `CriticalStretch` and the count is in the state
    body = stretched_body(RKCMaterial(; dmgmodel=FatigueDamage(3), monomial=:RK1))
    dh = submit(Job(body, DynamicRelaxation(steps=100)); quiet=true)
    (; storage, system) = dh.chunks[1]
    state = damage_state(storage)
    @test maximum(state.bond_exceedances) >= 3
    @test any(!, storage.bond_active)
    @test maximum(storage.damage) > 0
    # a bond that failed exceeded the stretch three times, an active one fewer than that
    @test all(state.bond_exceedances[.!storage.bond_active] .>= 3)
    @test all(state.bond_exceedances[storage.bond_active] .< 3)
    # the same bar with a model that never lets a bond fail stays intact
    body = stretched_body(RKCMaterial(; dmgmodel=FatigueDamage(typemax(Int)), monomial=:RK1))
    dh = submit(Job(body, DynamicRelaxation(steps=100)); quiet=true)
    @test all(dh.chunks[1].storage.bond_active)
    @test maximum(damage_state(dh.chunks[1].storage).bond_exceedances) >= 3
end
