@testsnippet J2Model begin
    using Peridynamics: constitutive_state, get_tensor, update_tensor!
    using Peridynamics.StaticArrays
    using LinearAlgebra

    # A history-dependent constitutive model: J2 plasticity with linear isotropic hardening
    # and a radial return, evaluated per bond.
    struct J2Plasticity <: Peridynamics.AbstractConstitutiveModel
        sigma_y::Float64
        hardening::Float64
    end
    J2Plasticity(; sigma_y=250.0, hardening=0.0) = J2Plasticity(sigma_y, hardening)

    Peridynamics.@cm_storage J2Plasticity struct J2PlasticityState
        bond_plastic_strain::BondTensor
        bond_eqps::BondScalar
    end

    function Peridynamics.first_piola_kirchhoff(cm::J2Plasticity, storage, params, F, idx,
                                                Δt)
        state = constitutive_state(storage)
        εᵖ = get_tensor(state.bond_plastic_strain, idx)
        ε = 0.5 .* (F' * F - I)
        εᵉ = ε - εᵖ
        σ = params.λ * tr(εᵉ) * I + 2 * params.μ * εᵉ
        s = σ - tr(σ) / 3 * I
        s_norm = sqrt(sum(abs2, s))
        α = state.bond_eqps[idx]
        f = s_norm - sqrt(2 / 3) * (cm.sigma_y + cm.hardening * α)
        if f > 0
            Δγ = f / (2 * params.μ + 2 / 3 * cm.hardening)
            εᵖ = εᵖ + Δγ .* (s ./ s_norm)
            update_tensor!(state.bond_plastic_strain, idx, εᵖ)
            state.bond_eqps[idx] = α + sqrt(2 / 3) * Δγ
            εᵉ = ε - εᵖ
            σ = params.λ * tr(εᵉ) * I + 2 * params.μ * εᵉ
        end
        return F * σ
    end

    # a bar pulled apart at both ends, so that it really deforms — a uniform initial
    # velocity would only translate it and never load the material. Relaxed quasi-statically
    # with `DynamicRelaxation`, whose pseudo time step is one, so the total stretch is
    # `v * steps`, here 2 % of the length.
    function stretched_body(mat; Δx=0.1, v=1e-4)
        pos, vol = uniform_box(1.0, 0.3, 0.3, Δx)
        body = Body(mat, pos, vol)
        material!(body; horizon=3.015Δx, rho=8e-6, E=2.1e5, nu=0.25, Gc=1e12)
        point_set!(x -> x < -0.4, body, :left)
        point_set!(x -> x > 0.4, body, :right)
        velocity_bc!(t -> -v, body, :left, :x)
        velocity_bc!(t -> v, body, :right, :x)
        return body
    end

    # `:RK1` reproduces a constant field exactly, also at a free surface, so a point that
    # is not deformed really sees `F = I` and the model does not yield spuriously
    rkc_plastic(; kwargs...) = RKCMaterial(; model=J2Plasticity(; kwargs...), monomial=:RK1)
end

@testitem "@cm_storage: the generated state of a constitutive model" setup=[J2Model] begin
    using Peridynamics: constitutive_storage_type, is_history_dependent, get_cm_storage,
                        storage_fields_expr

    cm = J2Plasticity()

    # the state is a parametric struct, exactly like a storage
    @test J2PlasticityState isa UnionAll
    @test J2PlasticityState <: Peridynamics.AbstractConstitutiveState
    @test fieldnames(J2PlasticityState) == (:bond_plastic_strain, :bond_eqps)
    @test constitutive_storage_type(cm) === J2PlasticityState{Float64,Matrix{Float64},
                                                              Vector{Float64}}
    @test constitutive_storage_type(cm, Float32) === J2PlasticityState{Float32,
                                                                       Matrix{Float32},
                                                                       Vector{Float32}}
    @test isconcretetype(constitutive_storage_type(cm))

    # declaring a state makes the model history-dependent
    @test is_history_dependent(cm)

    # a model without a state costs nothing and is not history-dependent
    @test constitutive_storage_type(SaintVenantKirchhoff()) === Nothing
    @test !is_history_dependent(SaintVenantKirchhoff())
    @test isnothing(get_cm_storage(SaintVenantKirchhoff(), VelocityVerlet(steps=1), nothing))

    # the declarations are inheritable, like those of a storage
    @test [d.name for d in storage_fields_expr(J2PlasticityState)] ==
          [:bond_plastic_strain, :bond_eqps]
end

@testitem "@cm_storage: the macro input checks" setup=[J2Model] begin
    # halo annotations are not supported: the state is chunk-local
    @test_throws LoadError @eval Peridynamics.@cm_storage J2Plasticity struct BadHaloState
        @lth some_field::PointScalar
    end

    # every field needs a field shape, there is no `init_field` hook for a model
    @test_throws LoadError @eval Peridynamics.@cm_storage J2Plasticity struct BadTypeState
        some_field::Vector{Float64}
    end

    # a model cannot carry the state of another model
    @test_throws LoadError @eval Peridynamics.@cm_storage J2Plasticity struct BadNestedState
        nested::ConstitutiveState
    end

    # only a constitutive model can declare a state
    struct NotAModel end
    @test_throws ArgumentError @eval Peridynamics.@cm_storage NotAModel struct NotAState
        a::BondScalar
    end
end

@testitem "ConstitutiveState: a storage carries the state of its model" setup=[J2Model] begin
    using Peridynamics: storage_type, has_constitutive_state, constitutive_state

    # the storage stays concrete whichever model is used, and the model fills the parameter
    for (mat, CMS) in ((RKCMaterial(), Nothing),
                       (RKCMaterial(; model=J2Plasticity()),
                        J2PlasticityState{Float64,Matrix{Float64},Vector{Float64}}))
        S = storage_type(mat)
        @test isconcretetype(S)
        @test has_constitutive_state(S)
        # read the parameter through the field, so the assertion survives a storage that
        # also carries a `DamageState` (its parameter would come after `CMS`)
        @test fieldtype(S, :cm_state) === CMS
    end

    # every material family with a constitutive model carries its state
    for mat in (CMaterial(), BACMaterial())
        @test has_constitutive_state(storage_type(mat))
        @test fieldtype(storage_type(mat), :cm_state) === Nothing
    end

    # a storage that does not declare the field answers `nothing`
    @test !has_constitutive_state(storage_type(BBMaterial()))

    # a storage may declare the state only once
    @test_throws LoadError @eval Peridynamics.@storage RKCMaterial struct TwiceState
        @inherit Peridynamics.VelocityVerletFields
        a::ConstitutiveState
        b::ConstitutiveState
    end

    # ... and never annotated for halo exchange
    @test_throws LoadError @eval Peridynamics.@storage RKCMaterial struct HaloState
        @inherit Peridynamics.VelocityVerletFields
        @lth a::ConstitutiveState
    end
end

@testitem "history-dependent model: the state is integrated by a simulation" tags=[:simulation] setup=[J2Model] begin
    using Peridynamics: constitutive_state

    body = stretched_body(rkc_plastic(; sigma_y=200.0, hardening=1e4))
    dh = submit(Job(body, DynamicRelaxation(steps=200)); quiet=true)
    state = constitutive_state(dh.chunks[1].storage)

    # the model integrated a history, and it is the state of the storage that holds it
    @test any(!iszero, state.bond_plastic_strain)
    @test 0 < maximum(state.bond_eqps) < 1
    @test all(isfinite, state.bond_eqps)
    @test length(state.bond_eqps) == Peridynamics.get_n_bonds(dh.chunks[1].system)

    # a yield stress that is never reached leaves the state untouched
    body = stretched_body(rkc_plastic(; sigma_y=1e12))
    dh = submit(Job(body, DynamicRelaxation(steps=200)); quiet=true)
    state = constitutive_state(dh.chunks[1].storage)
    @test all(iszero, state.bond_plastic_strain)
    @test all(iszero, state.bond_eqps)
end

@testitem "check_constitutive_model: history dependence is checked by Job" setup=[J2Model] begin
    using Peridynamics: check_constitutive_model, supports_history_dependence

    # a solver that evaluates the force density several times per step cannot integrate a
    # history correctly
    @test !supports_history_dependence(NewtonKrylov(steps=1, stepsize=1e-6))
    @test supports_history_dependence(VelocityVerlet(steps=1))
    @test supports_history_dependence(DynamicRelaxation(steps=1))
    body = stretched_body(rkc_plastic())
    err = try
        Job(body, NewtonKrylov(steps=1, stepsize=1e-6))
    catch e
        e
    end
    @test err isa Peridynamics.HistoryDependenceError
    msg = sprint(showerror, err)
    @test occursin("J2Plasticity", msg)
    @test occursin("NewtonKrylov", msg)
    @test occursin("more than once per time step", msg)

    # a material that does not support history dependence cannot index a state
    struct NoHistoryMat <: Peridynamics.AbstractMaterial end
    Peridynamics.supports_history_dependence(::NoHistoryMat) = false
    Peridynamics.get_constitutive_model(::NoHistoryMat) = J2Plasticity()
    err = try
        check_constitutive_model(NoHistoryMat(), VelocityVerlet(steps=1))
    catch e
        e
    end
    @test err isa Peridynamics.HistoryDependenceError
    @test occursin("NoHistoryMat", sprint(showerror, err))
    @test occursin("well defined index", sprint(showerror, err))

    # a storage that does not carry the state cannot hold the history
    struct NoStateMat <: Peridynamics.AbstractMaterial end
    Peridynamics.get_constitutive_model(::NoStateMat) = J2Plasticity()
    Peridynamics.@storage NoStateMat struct NoStateStorage
        @inherit Peridynamics.VelocityVerletFields
    end
    err = try
        check_constitutive_model(NoStateMat(), VelocityVerlet(steps=1))
    catch e
        e
    end
    @test err isa Peridynamics.HistoryDependenceError
    @test occursin("NoStateStorage", sprint(showerror, err))
    @test occursin("cm_state::ConstitutiveState", sprint(showerror, err))

    # a stateless model passes every check, also in a multibody setup
    @test Job(stretched_body(RKCMaterial()), NewtonKrylov(steps=1, stepsize=1e-6)) isa Job
    @test isnothing(check_constitutive_model(stretched_body(rkc_plastic()),
                                             VelocityVerlet(steps=1)))
    ms = MultibodySetup(:a => stretched_body(rkc_plastic()),
                        :b => stretched_body(RKCMaterial()))
    @test isnothing(check_constitutive_model(ms, VelocityVerlet(steps=1)))
    @test_throws Peridynamics.HistoryDependenceError check_constitutive_model(ms,
                                                                              NewtonKrylov(steps=1,
                                                                                           stepsize=1e-6))
end

@testitem "constitutive model interface: bridge, accessors and empty state" begin
    using Peridynamics: first_piola_kirchhoff, strain_energy_density, get_constitutive_model,
                        is_history_dependent, constitutive_state
    using Peridynamics.StaticArrays

    F = SMatrix{3,3,Float64,9}(1.01, 0, 0, 0, 1, 0, 0, 0, 1)
    pos, vol = uniform_box(1, 1, 1, 0.5)
    body = Body(CMaterial(), pos, vol)
    material!(body; horizon=1.5, rho=8e-6, E=2.1e5, nu=0.25, Gc=2.7)
    params = only(body.point_params)
    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    storage = dh.chunks[1].storage

    # the six-argument form is bridged to the four-argument form of a stateless model, so a
    # model that was written before the index and the time step existed keeps working
    for model in (SaintVenantKirchhoff(), LinearElastic(), NeoHooke(), NeoHookePenalty())
        @test first_piola_kirchhoff(model, storage, params, F, 1, 1e-8) ≈
              first_piola_kirchhoff(model, storage, params, F)
        @test strain_energy_density(model, storage, params, F, 1) ≈
              strain_energy_density(model, storage, params, F)
        @test !is_history_dependent(model)
    end

    # a model without any method gets an error that says which one to define
    struct EmptyModel <: Peridynamics.AbstractConstitutiveModel end
    err = try
        first_piola_kirchhoff(EmptyModel(), storage, params, F, 1, 1e-8)
    catch e
        e
    end
    @test err isa Peridynamics.InterfaceError
    @test occursin("first_piola_kirchhoff", sprint(showerror, err))
    err = try
        strain_energy_density(EmptyModel(), storage, params, F, 1)
    catch e
        e
    end
    @test err isa Peridynamics.InterfaceError
    @test occursin("strain_energy_density", sprint(showerror, err))

    # every material that has a constitutive model answers with it
    @test get_constitutive_model(CMaterial()) === SaintVenantKirchhoff()
    @test get_constitutive_model(CRMaterial()) === SaintVenantKirchhoff()
    @test get_constitutive_model(RKCMaterial()) === SaintVenantKirchhoff()
    @test get_constitutive_model(RKCRMaterial()) === SaintVenantKirchhoff()
    @test get_constitutive_model(BACMaterial()) === SaintVenantKirchhoff()
    @test isnothing(get_constitutive_model(BBMaterial()))
    @test isnothing(get_constitutive_model(OSBMaterial()))
    @test !is_history_dependent(nothing)

    # a stateless model leaves the state of the storage empty, and a storage that does not
    # declare the field at all answers `nothing` as well
    @test isnothing(constitutive_state(storage))
    bb = Body(BBMaterial(), pos, vol)
    material!(bb; horizon=1.5, rho=8e-6, E=2.1e5, Gc=2.7)
    bb_dh = Peridynamics.threads_data_handler(bb, VelocityVerlet(steps=1), 1)
    @test isnothing(constitutive_state(bb_dh.chunks[1].storage))
end
