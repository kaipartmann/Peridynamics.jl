@testsnippet RecorderModels begin
    using Peridynamics: constitutive_state, first_piola_kirchhoff, strain_energy_density,
                        get_tensor, update_tensor!
    using Peridynamics.StaticArrays

    # History-dependent constitutive models built for these tests. Mechanically they are
    # Saint-Venant-Kirchhoff and therefore transparent to a simulation, but they record in
    # their state what a history-dependent model relies on: that its state is indexed by the
    # quantity the material evaluates the model for (a bond or a point), that every entry is
    # advanced exactly once per time step, and that the model receives the time step of the
    # solver and the deformation gradient of the evaluated quantity. A plasticity model would
    # hide all of that behind its physics; the recorders make it observable and exact.
    struct BondRecorder <: Peridynamics.AbstractConstitutiveModel end   # per bond: RKC, BAC
    struct PointRecorder <: Peridynamics.AbstractConstitutiveModel end  # per point: C

    Peridynamics.@cm_storage BondRecorder struct BondRecorderState
        bond_evals::BondScalar{Int}
        bond_time::BondScalar
        bond_defgrad::BondTensor
    end

    Peridynamics.@cm_storage PointRecorder struct PointRecorderState
        point_evals::PointScalar{Int}
        point_time::PointScalar
        point_defgrad::PointTensor
    end

    function record!(evals, time, defgrad, idx, F, Δt)
        evals[idx] += 1
        time[idx] += Δt
        update_tensor!(defgrad, idx, F)
        return nothing
    end

    function Peridynamics.first_piola_kirchhoff(::BondRecorder, storage, params, F, idx, Δt)
        state = constitutive_state(storage)
        record!(state.bond_evals, state.bond_time, state.bond_defgrad, idx, F, Δt)
        return first_piola_kirchhoff(SaintVenantKirchhoff(), storage, params, F)
    end

    function Peridynamics.first_piola_kirchhoff(::PointRecorder, storage, params, F, idx, Δt)
        state = constitutive_state(storage)
        record!(state.point_evals, state.point_time, state.point_defgrad, idx, F, Δt)
        return first_piola_kirchhoff(SaintVenantKirchhoff(), storage, params, F)
    end

    # the strain energy density must not change the state, so it is the stateless form
    function Peridynamics.strain_energy_density(::Union{BondRecorder,PointRecorder}, storage,
                                                params, F)
        return strain_energy_density(SaintVenantKirchhoff(), storage, params, F)
    end

    # a small bar pulled apart at both ends, so that every point really deforms — a uniform
    # velocity would only translate it and never load the material
    function stretched_bar(mat; Δx=0.125, v=0.01)
        pos, vol = uniform_box(1.0, 0.5, 0.5, Δx)
        body = Body(mat, pos, vol)
        material!(body; horizon=3.015Δx, rho=8e-6, E=2.1e5, nu=0.25, Gc=1e12)
        point_set!(x -> x < -0.4, body, :left)
        point_set!(x -> x > 0.4, body, :right)
        velocity_bc!(t -> -v, body, :left, :x)
        velocity_bc!(t -> v, body, :right, :x)
        return body
    end
end

@testitem "@cm_storage: the generated state of a constitutive model" setup=[RecorderModels] begin
    using Peridynamics: constitutive_storage_type, is_history_dependent, get_cm_storage,
                        storage_fields_expr

    cm = BondRecorder()

    # the state is a parametric struct, exactly like a storage: the float type of the
    # simulation first, then one parameter per array, in declaration order
    @test BondRecorderState isa UnionAll
    @test BondRecorderState <: Peridynamics.AbstractConstitutiveState
    @test fieldnames(BondRecorderState) == (:bond_evals, :bond_time, :bond_defgrad)
    @test constitutive_storage_type(cm) ===
          BondRecorderState{Float64,Vector{Int},Vector{Float64},Matrix{Float64}}
    @test constitutive_storage_type(cm, Float32) ===
          BondRecorderState{Float32,Vector{Int},Vector{Float32},Matrix{Float32}}
    @test isconcretetype(constitutive_storage_type(cm))
    @test constitutive_storage_type(PointRecorder()) ===
          PointRecorderState{Float64,Vector{Int},Vector{Float64},Matrix{Float64}}

    # declaring a state makes the model history-dependent
    @test is_history_dependent(cm)
    @test is_history_dependent(PointRecorder())

    # a model without a state costs nothing and is not history-dependent
    @test constitutive_storage_type(SaintVenantKirchhoff()) === Nothing
    @test constitutive_storage_type(SaintVenantKirchhoff(), Float32) === Nothing
    @test !is_history_dependent(SaintVenantKirchhoff())
    @test isnothing(get_cm_storage(SaintVenantKirchhoff(), VelocityVerlet(steps=1), nothing))

    # the declarations are inheritable, like those of a storage
    @test [d.name for d in storage_fields_expr(BondRecorderState)] ==
          [:bond_evals, :bond_time, :bond_defgrad]

    # a state whose fields all fix their element type does not depend on the float type of
    # the simulation, so the generated type has no float type parameter
    struct CounterModel <: Peridynamics.AbstractConstitutiveModel end
    Peridynamics.@cm_storage CounterModel struct CounterState
        bond_evals::BondScalar{Int}
    end
    @test constitutive_storage_type(CounterModel()) === CounterState{Vector{Int}}
    @test constitutive_storage_type(CounterModel(), Float32) === CounterState{Vector{Int}}
    @test is_history_dependent(CounterModel())

    # a state without fields has no type parameters at all, but still marks the model as
    # history-dependent by declaration
    struct EmptyStateModel <: Peridynamics.AbstractConstitutiveModel end
    Peridynamics.@cm_storage EmptyStateModel struct EmptyState end
    @test constitutive_storage_type(EmptyStateModel()) === EmptyState
    @test !(EmptyState isa UnionAll)
    @test is_history_dependent(EmptyStateModel())
end

@testitem "@cm_storage: the macro input checks" setup=[RecorderModels] begin
    # halo annotations are not supported: the state is chunk-local
    @test_throws LoadError @eval Peridynamics.@cm_storage BondRecorder struct BadHaloState
        @lth some_field::PointScalar
    end
    @test_throws LoadError @eval Peridynamics.@cm_storage BondRecorder struct BadHaloState2
        @htl some_field::PointScalar
    end

    # every field needs a field shape, there is no `init_field` hook for a model
    @test_throws LoadError @eval Peridynamics.@cm_storage BondRecorder struct BadTypeState
        some_field::Vector{Float64}
    end

    # a model cannot carry the state of another model
    @test_throws LoadError @eval Peridynamics.@cm_storage BondRecorder struct BadNestedState
        nested::ConstitutiveState
    end
    @test_throws LoadError @eval Peridynamics.@cm_storage BondRecorder struct BadNestedState2
        nested::DamageState
    end

    # the messages say which field and which declaration is the problem
    err = try
        decl = Peridynamics.StorageFieldDecl(:a, :lth, Peridynamics.PointScalar,
                                             Peridynamics.PointScalar(), nothing)
        Peridynamics.check_cm_storage_decls([decl])
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("`a`", err.msg) && occursin("@lth", err.msg)
    @test occursin("@cm_storage", err.msg)

    # only a constitutive model can declare a state
    struct NotAModel end
    @test_throws ArgumentError @eval Peridynamics.@cm_storage NotAModel struct NotAState
        a::BondScalar
    end
    err = try
        Peridynamics.typecheck_constitutive_model(NotAModel)
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("NotAModel is not a valid constitutive model type", err.msg)
    @test occursin("AbstractConstitutiveModel", err.msg)
    # ... and the check also rejects a value that is not a type at all
    @test_throws ArgumentError Peridynamics.typecheck_constitutive_model(1)

    # the first argument has to name the model, the second has to be a struct definition
    @test_throws LoadError @eval Peridynamics.@cm_storage BondRecorder BondRecorderState
    @test_throws LoadError @eval Peridynamics.@cm_storage 1 struct BadModelState
        a::BondScalar
    end
end

@testitem "ConstitutiveState: a storage carries the state of its model" setup=[RecorderModels] begin
    using Peridynamics: storage_type, has_constitutive_state, constitutive_state

    # the storage stays concrete whichever model is used, and the model fills the parameter
    for (mat, CMS) in ((RKCMaterial(), Nothing),
                       (RKCMaterial(; model=BondRecorder()),
                        BondRecorderState{Float64,Vector{Int},Vector{Float64},
                                          Matrix{Float64}}),
                       (CMaterial(; model=PointRecorder()),
                        PointRecorderState{Float64,Vector{Int},Vector{Float64},
                                           Matrix{Float64}}),
                       (BACMaterial(; model=BondRecorder()),
                        BondRecorderState{Float64,Vector{Int},Vector{Float64},
                                          Matrix{Float64}}))
        S = storage_type(mat)
        @test isconcretetype(S)
        @test has_constitutive_state(S)
        # read the parameter through the field, so the assertion survives a storage that
        # also carries a `DamageState` (its parameter would come after `CMS`)
        @test fieldtype(S, :cm_state) === CMS
    end

    # the state follows the float type of the simulation
    S32 = storage_type(RKCMaterial(; model=BondRecorder()), Float32)
    @test fieldtype(S32, :cm_state) ===
          BondRecorderState{Float32,Vector{Int},Vector{Float32},Matrix{Float32}}

    # every material family with a constitutive model carries its state
    for mat in (CMaterial(), BACMaterial())
        @test has_constitutive_state(storage_type(mat))
        @test fieldtype(storage_type(mat), :cm_state) === Nothing
    end

    # a storage that does not declare the field answers `nothing`
    @test !has_constitutive_state(storage_type(BBMaterial()))
    @test !has_constitutive_state(Nothing)

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

@testitem "get_cm_storage: the state is allocated per chunk and adaptable" setup=[Fixtures, RecorderModels] begin
    using Peridynamics: get_cm_storage, init_constitutive_state, constitutive_state,
                        get_n_bonds, get_n_loc_points, get_n_points

    # a minimal stand-in for the array type of another backend
    struct WrappedArray{T,N} <: AbstractArray{T,N}
        a::Array{T,N}
    end
    Base.size(x::WrappedArray) = size(x.a)
    Base.getindex(x::WrappedArray, i...) = getindex(x.a, i...)
    struct WrappedBackend end
    function Peridynamics.Adapt.adapt_storage(::WrappedBackend, a::Array{T,N}) where {T,N}
        return WrappedArray{T,N}(a)
    end

    body = stretched_bar(RKCMaterial(; model=BondRecorder()))
    solver = VelocityVerlet(steps=1)
    dh = Fixtures.handler(body, solver; n_chunks=2, init=false)
    for chunk in dh.chunks
        state = constitutive_state(chunk.storage)
        @test state isa BondRecorderState
        @test state === chunk.storage.cm_state
        # a bond field is sized by the bonds of the chunk, and the state starts at zero
        n_bonds = get_n_bonds(chunk.system)
        @test length(state.bond_evals) == n_bonds
        @test length(state.bond_time) == n_bonds
        @test size(state.bond_defgrad) == (9, n_bonds)
        @test all(iszero, state.bond_evals)
        @test all(iszero, state.bond_defgrad)
        # the state is allocated by the model and reached through the material
        direct = get_cm_storage(BondRecorder(), solver, chunk.system)
        @test typeof(direct) === typeof(state)
        @test length(direct.bond_evals) == n_bonds
        @test typeof(init_constitutive_state(body.mat, solver, chunk.system)) ===
              typeof(state)
        # the generated `Adapt.adapt_structure` rebuilds the state from its adapted arrays,
        # so the state follows the storage to another array backend
        adapted = Peridynamics.Adapt.adapt(WrappedBackend(), state)
        @test adapted isa BondRecorderState
        @test adapted.bond_evals isa WrappedArray{Int,1}
        @test adapted.bond_time isa WrappedArray{Float64,1}
        @test adapted.bond_defgrad isa WrappedArray{Float64,2}
        @test adapted.bond_defgrad == state.bond_defgrad
    end

    # a point field of the state covers the local points only: the state is chunk-local
    body = stretched_bar(CMaterial(; model=PointRecorder()))
    dh = Fixtures.handler(body, solver; n_chunks=2, init=false)
    for chunk in dh.chunks
        state = constitutive_state(chunk.storage)
        @test state isa PointRecorderState
        @test length(state.point_evals) == get_n_loc_points(chunk) < get_n_points(chunk)
        @test size(state.point_defgrad) == (9, get_n_loc_points(chunk))
    end

    # a material with a stateless model leaves the field empty
    body = stretched_bar(RKCMaterial())
    dh = Fixtures.handler(body, solver; n_chunks=2, init=false)
    @test all(isnothing(constitutive_state(c.storage)) for c in dh.chunks)
    @test isnothing(init_constitutive_state(body.mat, solver, dh.chunks[1].system))
end

@testitem "history-dependent model: the state is advanced once per step by a simulation" tags=[:simulation] setup=[Fixtures, RecorderModels] begin
    using Peridynamics: constitutive_state, get_n_bonds, get_n_loc_points

    # The contract a history-dependent model relies on: for every bond, its state is
    # advanced exactly once per time step, with the time step of the solver and the
    # deformation gradient of that bond. The data handler evaluates the force density once
    # more when it is initialized, in the reference configuration at `t = 0`, so every
    # entry sees `steps + 1` evaluations — the one thing a model must not assume is that
    # the first evaluation already belongs to the first step.
    steps = 4
    body = stretched_bar(RKCMaterial(; model=BondRecorder()); v=1e3)
    vv = VelocityVerlet(; steps)
    dh = Fixtures.run_threads(Job(body, vv); n_chunks=2)
    Δt = vv.Δt
    params = only(body.point_params)
    svk = SaintVenantKirchhoff()
    for chunk in dh.chunks
        state = constitutive_state(chunk.storage)
        @test length(state.bond_evals) == get_n_bonds(chunk.system)
        @test all(==(steps + 1), state.bond_evals)
        @test all(≈((steps + 1) * Δt), state.bond_time)
        # the bar is stretched, so the bonds near the loaded ends are really deformed ...
        F11 = @view state.bond_defgrad[1, :]
        @test maximum(F11) > 1 + 1e-4
        @test all(isfinite, state.bond_defgrad)
        # ... and the recorded deformation gradient is the one of the last step: the stress
        # the material stored in that step is the stress of the recorded gradient
        @test all(1:get_n_bonds(chunk.system)) do bond_id
            Fij = get_tensor(state.bond_defgrad, bond_id)
            Pij = get_tensor(chunk.storage.bond_first_piola_kirchhoff, bond_id)
            return Pij ≈ first_piola_kirchhoff(svk, chunk.storage, params, Fij)
        end
    end

    # the recorder is mechanically Saint-Venant-Kirchhoff, so it does not change the result
    dh_svk = Fixtures.run_threads(Job(stretched_bar(RKCMaterial(); v=1e3),
                                      VelocityVerlet(; steps)); n_chunks=2)
    for (chunk, chunk_svk) in zip(dh.chunks, dh_svk.chunks)
        @test chunk.storage.displacement ≈ chunk_svk.storage.displacement
        @test chunk.storage.bond_first_piola_kirchhoff ≈
              chunk_svk.storage.bond_first_piola_kirchhoff
    end

    # per point for the correspondence formulation, and `DynamicRelaxation` integrates
    # over its pseudo time step of one
    steps = 3
    body = stretched_bar(CMaterial(; model=PointRecorder()))
    dh = Fixtures.run_threads(Job(body, DynamicRelaxation(; steps)); n_chunks=2)
    for chunk in dh.chunks
        state = constitutive_state(chunk.storage)
        @test length(state.point_evals) == get_n_loc_points(chunk)
        @test all(==(steps + 1), state.point_evals)
        @test all(==((steps + 1) * 1.0), state.point_time)
        @test all(isfinite, state.point_defgrad)
    end

    # per bond for the bond-associated formulation as well
    steps = 2
    body = stretched_bar(BACMaterial(; model=BondRecorder()))
    dh = Fixtures.run_threads(Job(body, VelocityVerlet(; steps)); n_chunks=2)
    for chunk in dh.chunks
        state = constitutive_state(chunk.storage)
        @test length(state.bond_evals) == get_n_bonds(chunk.system)
        @test all(==(steps + 1), state.bond_evals)
    end
end

@testitem "check_constitutive_model: history dependence is checked by Job" setup=[RecorderModels] begin
    using Peridynamics: check_constitutive_model, supports_history_dependence

    # a solver that evaluates the force density several times per step cannot integrate a
    # history correctly
    @test !supports_history_dependence(NewtonKrylov(steps=1, stepsize=1e-6))
    @test supports_history_dependence(VelocityVerlet(steps=1))
    @test supports_history_dependence(DynamicRelaxation(steps=1))
    body = stretched_bar(RKCMaterial(; model=BondRecorder()))
    err = try
        Job(body, NewtonKrylov(steps=1, stepsize=1e-6))
    catch e
        e
    end
    @test err isa Peridynamics.HistoryDependenceError
    msg = sprint(showerror, err)
    @test occursin("BondRecorder", msg)
    @test occursin("RKCMaterial", msg)
    @test occursin("NewtonKrylov", msg)
    @test occursin("more than once per time step", msg)
    @test occursin("VelocityVerlet", msg)

    # a material that does not support history dependence cannot index a state
    struct NoHistoryMat <: Peridynamics.AbstractMaterial end
    Peridynamics.supports_history_dependence(::NoHistoryMat) = false
    Peridynamics.get_constitutive_model(::NoHistoryMat) = BondRecorder()
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
    Peridynamics.get_constitutive_model(::NoStateMat) = BondRecorder()
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

    # a model that declares itself history-dependent without a state only needs a solver
    # that evaluates it once per step, not a storage that carries a state
    struct HistoryNoStateModel <: Peridynamics.AbstractConstitutiveModel end
    Peridynamics.is_history_dependent(::HistoryNoStateModel) = true
    struct HistoryNoStateMat <: Peridynamics.AbstractMaterial end
    Peridynamics.get_constitutive_model(::HistoryNoStateMat) = HistoryNoStateModel()
    Peridynamics.@storage HistoryNoStateMat struct HistoryNoStateStorage
        @inherit Peridynamics.VelocityVerletFields
    end
    @test isnothing(check_constitutive_model(HistoryNoStateMat(), VelocityVerlet(steps=1)))
    @test_throws Peridynamics.HistoryDependenceError check_constitutive_model(HistoryNoStateMat(),
                                                                              NewtonKrylov(steps=1,
                                                                                           stepsize=1e-6))

    # a stateless model passes every check, also in a multibody setup
    @test Job(stretched_bar(RKCMaterial()), NewtonKrylov(steps=1, stepsize=1e-6)) isa Job
    @test isnothing(check_constitutive_model(stretched_bar(RKCMaterial(; model=BondRecorder())),
                                             VelocityVerlet(steps=1)))
    @test Job(stretched_bar(RKCMaterial(; model=BondRecorder())), VelocityVerlet(steps=1)) isa Job
    ms = MultibodySetup(:a => stretched_bar(RKCMaterial(; model=BondRecorder())),
                        :b => stretched_bar(RKCMaterial()))
    @test isnothing(check_constitutive_model(ms, VelocityVerlet(steps=1)))
    @test_throws Peridynamics.HistoryDependenceError check_constitutive_model(ms,
                                                                              NewtonKrylov(steps=1,
                                                                                           stepsize=1e-6))
    @test_throws Peridynamics.HistoryDependenceError Job(ms, NewtonKrylov(steps=1,
                                                                          stepsize=1e-6))
end

@testitem "constitutive model interface: bridge, accessors and empty state" begin
    using Peridynamics: first_piola_kirchhoff, strain_energy_density, get_constitutive_model,
                        is_history_dependent, constitutive_state, constitutive_model_hint
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
    msg = sprint(showerror, err)
    @test occursin("first_piola_kirchhoff", msg)
    @test occursin("EmptyModel", msg)
    @test occursin("function Peridynamics.first_piola_kirchhoff(::EmptyModel, storage, params, F)",
                   msg)
    @test occursin("history-dependent", msg)
    err = try
        strain_energy_density(EmptyModel(), storage, params, F, 1)
    catch e
        e
    end
    @test err isa Peridynamics.InterfaceError
    @test occursin("strain_energy_density", sprint(showerror, err))
    # the hint names the model whether it is given as an instance or as a type
    @test occursin("`EmptyModel`", constitutive_model_hint(EmptyModel(), "f"))
    @test occursin("`EmptyModel`", constitutive_model_hint(EmptyModel, "f"))

    # every material that has a constitutive model answers with it
    @test get_constitutive_model(CMaterial()) === SaintVenantKirchhoff()
    @test get_constitutive_model(CRMaterial()) === SaintVenantKirchhoff()
    @test get_constitutive_model(RKCMaterial()) === SaintVenantKirchhoff()
    @test get_constitutive_model(RKCRMaterial()) === SaintVenantKirchhoff()
    @test get_constitutive_model(BACMaterial()) === SaintVenantKirchhoff()
    @test get_constitutive_model(CMaterial(; model=NeoHooke())) === NeoHooke()
    @test isnothing(get_constitutive_model(BBMaterial()))
    @test isnothing(get_constitutive_model(OSBMaterial()))
    @test isnothing(get_constitutive_model(CKIMaterial()))
    @test !is_history_dependent(nothing)

    # a stateless model leaves the state of the storage empty, and a storage that does not
    # declare the field at all answers `nothing` as well
    @test isnothing(constitutive_state(storage))
    bb = Body(BBMaterial(), pos, vol)
    material!(bb; horizon=1.5, rho=8e-6, E=2.1e5, Gc=2.7)
    bb_dh = Peridynamics.threads_data_handler(bb, VelocityVerlet(steps=1), 1)
    @test isnothing(constitutive_state(bb_dh.chunks[1].storage))
end
