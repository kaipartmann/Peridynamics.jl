@testitem "required_fields_fracture" begin
    @test Peridynamics.required_fields_fracture(Peridynamics.AbstractMaterial) == ()

    rff_bond_system = (:damage, :n_active_bonds, :bond_active)
    @test Peridynamics.required_fields_fracture(BBMaterial) === rff_bond_system

    rff_interaction_system = (:damage, :n_active_one_nis, :one_ni_active)
    @test Peridynamics.required_fields_fracture(CKIMaterial) === rff_interaction_system
end

@testitem "get_frac_params: the fracture keywords of CriticalStretch" begin
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
end

@testitem "FractureParameters: the block resolves Gc and εc through the damage model" begin
    import Peridynamics: FractureParameters, param_fields_expr, is_provided

    # the block is one `@derived` group, resolved by `get_frac_params` of the damage model
    spec = param_fields_expr(FractureParameters)
    @test [decl.name for decl in spec.decls] == [:Gc, :εc]
    @test all(is_provided, spec.decls)
    @test spec.kwargs == [:Gc, :epsilon_c]

    # both parameters carry a simulation-log label
    labels = Dict(decl.name => decl.label for decl in spec.decls)
    @test labels[:Gc] == "critical energy release rate"
    @test labels[:εc] == "critical stretch"
end

@testitem "damage model hooks: the defaults of a model that deletes bonds" begin
    import Peridynamics: kinematic_weight, bond_integrity, log_dmgmodel,
                         damage_storage_type, get_dmg_storage, init_damage_state, damage_state,
                         has_damage_state, storage_type, req_storage_fields

    pos, vol = uniform_box(1, 1, 1, 0.5)
    body = Body(BBMaterial(), pos, vol)
    material!(body; horizon=1.5, rho=1, E=1, Gc=1)
    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    (; storage, system) = dh.chunks[1]
    dmg = CriticalStretch()

    # no softening: a model that deletes bonds has fully intact bonds and full trust
    @test kinematic_weight(dmg, storage, 1) === 1.0
    @test bond_integrity(dmg, storage, 1) === 1.0

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

    # the model owns the standard fracture keywords: inheriting `FractureParameters`
    # registers `Gc`/`epsilon_c` and resolves them through `get_frac_params` below
    Peridynamics.@dmg_params FatigueDamage struct FatigueDamageParameters
        @inherit FractureParameters
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

@testitem "bond integrity and kinematic weight: the wiring into the RKC force path" setup=[FatigueModel] begin
    import Peridynamics: bond_integrity, kinematic_weight, get_params, each_point_idx

    # a stateless damage model with constant softening factors, so every wired-in factor
    # shows up as an exact scaling relative to the unsoftened reference
    struct ConstSoftening <: Peridynamics.AbstractDamageModel
        wkin::Float64
        g::Float64
    end
    Peridynamics.@dmg_params ConstSoftening struct ConstSofteningParameters
        @inherit FractureParameters
    end
    function Peridynamics.get_frac_params(::ConstSoftening, δ, K; kwargs...)
        return Peridynamics.get_frac_params(CriticalStretch(), δ, K; kwargs...)
    end
    function Peridynamics.has_fracture(::ConstSoftening, params)
        return Peridynamics.has_fracture(CriticalStretch(), params)
    end
    function Peridynamics.calc_failure!(storage, system, mat, ::ConstSoftening, paramsetup,
                                        i)
        for bond_id in Peridynamics.each_bond_idx(system, i)
            storage.n_active_bonds[i] += storage.bond_active[bond_id]
        end
        return nothing
    end
    @inline function Peridynamics.kinematic_weight(dmg::ConstSoftening,
                                                   ::Peridynamics.AbstractStorage, bond_id)
        return dmg.wkin
    end
    @inline function Peridynamics.bond_integrity(dmg::ConstSoftening,
                                                 ::Peridynamics.AbstractStorage, bond_id)
        return dmg.g
    end

    function force_calc(dmgmodel)
        body = stretched_body(RKCMaterial(; dmgmodel, monomial=:RK1))
        dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
        chunk = dh.chunks[1]
        chunk.storage.position[1, :] .*= 1.01 # uniaxial stretch, so P and Ψ are nonzero
        Peridynamics.calc_weights_and_defgrad!(chunk, 0.0, 1e-7)
        Peridynamics.calc_force_density!(chunk, 0.0, 1e-7)
        for i in each_point_idx(chunk.system)
            params = get_params(chunk.paramsetup, i)
            Peridynamics.strain_energy_density_point!(chunk.storage, chunk.system,
                                                      chunk.mat, params, i)
        end
        return chunk.storage
    end

    ref = force_calc(ConstSoftening(1.0, 1.0))

    # a uniform kinematic weight scales the weighted volume, but cancels in the
    # least-squares fit: the gradient weights and the deformation gradient are invariant
    weighted = force_calc(ConstSoftening(0.5, 1.0))
    @test weighted.weighted_volume ≈ 0.5 .* ref.weighted_volume
    @test weighted.gradient_weight ≈ ref.gradient_weight
    @test weighted.defgrad ≈ ref.defgrad

    # the integrity scales the stress of every bond and with it everything linear in it:
    # the internal force density and the strain energy density
    softened = force_calc(ConstSoftening(1.0, 0.25))
    @test softened.weighted_volume ≈ ref.weighted_volume
    @test softened.bond_first_piola_kirchhoff ≈ 0.25 .* ref.bond_first_piola_kirchhoff
    @test softened.b_int ≈ 0.25 .* ref.b_int
    @test softened.strain_energy_density ≈ 0.25 .* ref.strain_energy_density

    # a vanishing kinematic weight leaves every point without a family for the fit: the
    # weighted volume is zero, the bonds stay active, and instead of `1 / 0` in the stress
    # integral such an isolated point transmits no stress at all, see `isolated_point`
    isolated = force_calc(ConstSoftening(0.0, 1.0))
    @test all(iszero, isolated.weighted_volume)
    @test all(isolated.bond_active)
    @test all(iszero, isolated.gradient_weight)
    @test all(iszero, isolated.b_int)
    @test all(iszero, isolated.bond_first_piola_kirchhoff)
    @test all(iszero, isolated.strain_energy_density)
    @test !any(isnan, isolated.b_int)

    # the rotated material takes the same path
    body = stretched_body(RKCRMaterial(; dmgmodel=ConstSoftening(0.0, 1.0)))
    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    chunk.storage.position[1, :] .*= 1.01
    Peridynamics.calc_weights_and_defgrad!(chunk, 0.0, 1e-7)
    Peridynamics.calc_force_density!(chunk, 0.0, 1e-7)
    @test all(iszero, chunk.storage.weighted_volume)
    @test all(iszero, chunk.storage.b_int)
    @test all(iszero, chunk.storage.bond_first_piola_kirchhoff)
end

@testitem "softening support: ignored hooks fail at Job creation" begin
    import Peridynamics: supports_bond_integrity, supports_kinematic_weight,
                         check_damage_model, SofteningSupportError

    # two models that each define one softening hook, but nothing else special
    struct SofteningNotSupported <: Peridynamics.AbstractDamageModel end
    struct WeightNotSupported <: Peridynamics.AbstractDamageModel end
    for D in (SofteningNotSupported, WeightNotSupported)
        @eval begin
            function Peridynamics.get_frac_params(::$D, δ, K; kwargs...)
                return Peridynamics.get_frac_params(CriticalStretch(), δ, K; kwargs...)
            end
            function Peridynamics.has_fracture(::$D, params)
                return Peridynamics.has_fracture(CriticalStretch(), params)
            end
        end
    end
    Peridynamics.@dmg_params SofteningNotSupported struct SofteningNotSupportedParameters
        @inherit FractureParameters
    end
    Peridynamics.@dmg_params WeightNotSupported struct WeightNotSupportedParameters
        @inherit FractureParameters
    end

    @inline function Peridynamics.bond_integrity(::SofteningNotSupported,
                                                 ::Peridynamics.AbstractStorage, bond_id)
        return 0.5
    end
    @inline function Peridynamics.kinematic_weight(::WeightNotSupported,
                                                   ::Peridynamics.AbstractStorage, bond_id)
        return 0.5
    end

    # the RKC family declares support for both hooks, every other material answers false
    @test !supports_bond_integrity(BBMaterial())
    @test !supports_kinematic_weight(BBMaterial())
    @test !supports_bond_integrity(CMaterial())
    @test supports_bond_integrity(RKCMaterial())
    @test supports_kinematic_weight(RKCMaterial())
    @test supports_bond_integrity(RKCRMaterial())
    @test supports_kinematic_weight(RKCRMaterial())

    function body_with(mat; kwargs...)
        pos, vol = uniform_box(1, 1, 1, 0.5)
        body = Body(mat, pos, vol)
        material!(body; horizon=1.5, rho=1, E=1, kwargs...)
        velocity_bc!(t -> 0.1, body, :all_points, :x)
        return body
    end

    # a material that ignores a defined hook fails at Job creation and names everything
    err = try
        Job(body_with(BBMaterial(; dmgmodel=SofteningNotSupported()); Gc=1.0),
            VelocityVerlet(steps=1))
    catch e
        e
    end
    @test err isa SofteningSupportError
    msg = sprint(showerror, err)
    @test contains(msg, "bond_integrity")
    @test contains(msg, "BBMaterial")
    @test contains(msg, "SofteningNotSupported")
    @test contains(msg, "supports_bond_integrity")

    # the kinematic weight alone triggers the check as well
    err = try
        Job(body_with(BBMaterial(; dmgmodel=WeightNotSupported()); Gc=1.0),
            VelocityVerlet(steps=1))
    catch e
        e
    end
    @test err isa SofteningSupportError
    @test contains(sprint(showerror, err), "kinematic_weight")

    # the RKC family honors the hooks, so the same models pass
    job = Job(body_with(RKCMaterial(; dmgmodel=SofteningNotSupported()); nu=0.25,
                        epsilon_c=0.01), VelocityVerlet(steps=1))
    @test job isa Job
    @test isnothing(check_damage_model(RKCMaterial(; dmgmodel=WeightNotSupported())))

    # a model without hook methods passes with every material
    @test Job(body_with(BBMaterial(); Gc=1.0), VelocityVerlet(steps=1)) isa Job
    @test isnothing(check_damage_model(BBMaterial()))
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

@testitem "CriticalStretch: Gc and εc live in the damage model parameters" begin
    using Peridynamics: CriticalStretchParameters, damage_param_type, damage_param_kwargs,
                        has_fracture

    @test damage_param_type(CriticalStretch(), Float64) === CriticalStretchParameters{Float64}
    @test damage_param_kwargs(CriticalStretch()) == (:Gc, :epsilon_c)
    @test isbitstype(CriticalStretchParameters{Float64})

    pos, vol = uniform_box(1.0, 1.0, 1.0, 0.5)
    body = Body(BBMaterial(), pos, vol)
    material!(body; horizon=1.5, rho=8e-6, E=2.1e5, Gc=2.7)
    par = only(body.point_params)
    @test par.dmg_params isa CriticalStretchParameters{Float64}
    @test has_fracture(CriticalStretch(), par)
    # failure permissions were granted through the flat reads of `Gc` and `εc`
    @test all(body.fail_permit)

    # without fracture keywords the parameters resolve to zero and failure stays prohibited
    body0 = Body(BBMaterial(), pos, vol)
    material!(body0; horizon=1.5, rho=8e-6, E=2.1e5)
    @test only(body0.point_params).Gc == 0.0
    @test !any(body0.fail_permit)

    # a damage model without parameters prohibits failure by default
    struct FPNoParamsDamage <: Peridynamics.AbstractDamageModel end
    @test !has_fracture(FPNoParamsDamage(), par)
end
