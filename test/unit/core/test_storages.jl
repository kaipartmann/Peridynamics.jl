# The storage contract and the `@storage` macro of `src/core/storages.jl`. The field
# declaration framework (shapes, blocks, allocation) is tested in `test_storage_fields.jl`.

@testitem "required_fields: the type-level part of the storage contract" begin
    # the fracture bookkeeping belongs to the damage model and is checked at Job creation,
    # so nothing is left that could be known from the material type alone
    @test Peridynamics.required_fields(Peridynamics.AbstractMaterial) === ()
    @test Peridynamics.required_fields(BBMaterial) === ()
    @test Peridynamics.required_fields(CKIMaterial) === ()
end

@testitem "req_storage_fields: material, damage model and time solver" begin
    # materials without a declared contract
    @test Peridynamics.req_storage_fields(BBMaterial()) === ()
    @test Peridynamics.req_storage_fields(CKIMaterial()) === ()

    # the RKC family declares the fields its inherited code reads
    rf_rkc = (:defgrad, :weighted_volume, :gradient_weight, :bond_first_piola_kirchhoff,
              :update_gradients)
    @test Peridynamics.req_storage_fields(RKCMaterial()) === rf_rkc
    @test Peridynamics.req_storage_fields(RKCRMaterial()) === rf_rkc

    # damage models are dispatched together with the material, and a material without a
    # damage model must not error; a model with a state needs the marker that carries it
    @test Peridynamics.req_storage_fields(BBMaterial(), CriticalStretch()) === (:dmg_state,)
    @test Peridynamics.req_storage_fields(BBMaterial(), nothing) === ()

    # time solvers are dispatched on the instance that is actually used
    @test Peridynamics.req_storage_fields(VelocityVerlet(steps=1)) ===
          Peridynamics.required_fields_timesolver(VelocityVerlet)
    @test in(:residual, Peridynamics.req_storage_fields(NewtonKrylov(steps=1)))

    # a solver that does not declare its fields cannot be checked and contributes nothing
    struct SilentSolver <: Peridynamics.AbstractTimeSolver end
    @test Peridynamics.req_storage_fields(SilentSolver()) === ()
end

@testitem "check_storage_contract: the contract is checked when a Job is created" begin
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection,
                         StorageContractError, check_storage_contract

    function testbody(mat)
        pos, vol = uniform_box(1, 1, 1, 0.5)
        body = Body(mat, pos, vol)
        material!(body; horizon=1.5, rho=1, E=1, nu=0.25, Gc=1.0)
        velocity_bc!(t -> 0.0, body, :all_points, 1)
        return body
    end
    vv = VelocityVerlet(steps=1)
    nk = NewtonKrylov(steps=1)

    # all default combinations of the package fulfill the contract
    for mat in (BBMaterial(), OSBMaterial(), CKIMaterial(), CMaterial(), CRMaterial(),
                RKCMaterial(), RKCRMaterial(), BACMaterial())
        @test isnothing(check_storage_contract(mat, vv))
        @test Job(testbody(mat), vv) isa Job
    end
    @test isnothing(check_storage_contract(RKCMaterial(), nk))

    # `RKCRStorage` does not carry the fields of the `NewtonKrylov` solver
    @test_throws StorageContractError check_storage_contract(RKCRMaterial(), nk)
    @test_throws StorageContractError Job(testbody(RKCRMaterial()), nk)

    # a material may ask for a field its storage does not have
    struct ContractMat <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::CriticalStretch
    end
    ContractMat() = ContractMat(CriticalStretch())
    Peridynamics.@params ContractMat struct ContractMatParams
        @inherit StandardParameters
    end
    Peridynamics.@storage ContractMat struct ContractStorage
        @inherit VelocityVerletFields
        dmg_state::DamageState
    end
    function Peridynamics.force_density_point!(::ContractStorage, system, ::ContractMat,
                                               params, t, Δt, i)
        return nothing
    end
    Peridynamics.req_storage_fields(::ContractMat) = (:my_field,)
    @test_throws StorageContractError check_storage_contract(ContractMat(), vv)
    @test_throws StorageContractError Job(testbody(ContractMat()), vv)

    # every body of a multibody setup is checked
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false
    b_ok = testbody(RKCMaterial())
    b_bad = testbody(ContractMat())
    @test Job(MultibodySetup(:a => testbody(BBMaterial()), :b => b_ok), vv) isa Job
    @test_throws StorageContractError begin
        Job(MultibodySetup(:a => testbody(BBMaterial()), :b => b_bad), vv)
    end
    Peridynamics.MPI_RUN[] = mpi_run_current_value

    # the error names every missing field and the reason why it is required
    err = try
        check_storage_contract(RKCRMaterial(), nk)
    catch e
        e
    end
    @test err isa StorageContractError
    @test err.storage <: Peridynamics.RKCRStorage
    @test :residual in first.(err.missing_fields)
    msg = sprint(showerror, err)
    @test contains(msg, "RKCRStorage")
    @test contains(msg, "residual")
    @test contains(msg, "required by the time solver `NewtonKrylov`")

    err = try
        check_storage_contract(ContractMat(), vv)
    catch e
        e
    end
    @test err isa StorageContractError
    @test first.(err.missing_fields) == [:my_field]
    @test contains(sprint(showerror, err), "required by the material `ContractMat`")
end

@testitem "get_storage_header: name and supertype of a storage definition" begin
    import Peridynamics: get_storage_header

    # Test case 1: Simple struct declaration
    expr1 = :(struct MyStorage end)
    type1, supertype1 = get_storage_header(expr1)
    @test type1 == :MyStorage
    @test supertype1 == :(Peridynamics.AbstractStorage)

    # Test case 2: Struct declaration with a subtype
    expr2 = :(struct MyStorage <: Peridynamics.AbstractStorage end)
    type2, supertype2 = get_storage_header(expr2)
    @test type2 == :MyStorage
    @test supertype2 == :(Peridynamics.AbstractStorage)

    # Test case 3: the type parameters are derived from the field declarations, so a storage
    # must not declare its own
    expr3 = :(struct MyStorage{A,B,C} <: Peridynamics.AbstractStorage end)
    err3 = try
        get_storage_header(expr3)
    catch e
        e
    end
    @test err3 isa ArgumentError
    @test contains(err3.msg, "cannot declare its own type parameters")

    # Test case 4: Parametric struct declaration without a subtype
    expr4 = :(struct MyStorage{A,B,C} end)
    err4 = try
        get_storage_header(expr4)
    catch e
        e
    end
    @test err4 isa ArgumentError
    @test contains(err4.msg, "cannot declare its own type parameters")

    # Test case 5: anything else is still rejected as an unsupported header
    expr5 = :(struct (a + b) end)
    err5 = try
        get_storage_header(expr5)
    catch e
        e
    end
    @test err5 isa ArgumentError
    @test contains(err5.msg, "not supported")
end

@testitem "@storage: the macro input checks" begin

    input = :(Peridynamics.BBStorage)
    @test isnothing(Peridynamics.macrocheck_input_storage_type(input))

    input = :(MyNonexistingStorage)
    @test isnothing(Peridynamics.macrocheck_input_storage_type(input))

    input = :(Base.MyNonexistingStorage)
    @test isnothing(Peridynamics.macrocheck_input_storage_type(input))

    input = :(:MyNonexistingStorageAsSymbol)
    @test_throws ArgumentError Peridynamics.macrocheck_input_storage_type(input)

    input = :(
        struct MyTestStorage <: Peridynamics.AbstractStorage
            a::Int
        end
    )
    @test isnothing(Peridynamics.macrocheck_input_storage_struct(input))

    input = :(
        struct MyTestStorage{A,B,C} <: Peridynamics.AbstractStorage
            c::Vector{Float64}
        end
    )
    @test isnothing(Peridynamics.macrocheck_input_storage_struct(input))

    input = :(
        @kwdef struct MyTestStorage <: Peridynamics.AbstractStorage
            c::Vector{Float64}
        end
    )
    @test_throws ArgumentError Peridynamics.macrocheck_input_storage_struct(input)

    input = :(Peridynamics.VelocityVerlet)
    @test isnothing(Peridynamics.macrocheck_input_timesolver(input))

    input = :(MyNonexistingSolver)
    @test isnothing(Peridynamics.macrocheck_input_timesolver(input))

    input = :(Base.MyNonexistingSolver)
    @test isnothing(Peridynamics.macrocheck_input_timesolver(input))

    input = :(:MyNonexistingSolverAsSymbol)
    @test_throws ArgumentError Peridynamics.macrocheck_input_timesolver(input)

    input = :(:testfield)
    @test isnothing(Peridynamics.macrocheck_input_field(input))

    input = :testfield
    @test_throws ArgumentError Peridynamics.macrocheck_input_field(input)
end

@testitem "@storage: generated methods and the checks at macro expansion" begin
    import Peridynamics: @storage, AbstractBondSystemMaterial, NoCorrection,
                         AbstractInteractionSystemMaterial, InterfaceError,
                         AbstractTimeSolver, AbstractSystem

    struct Mat1 <: AbstractBondSystemMaterial{NoCorrection} end
    struct Mat2 <: AbstractInteractionSystemMaterial end

    @test_throws InterfaceError Peridynamics.storage_type(Mat1())
    @test_throws InterfaceError Peridynamics.get_storage(Mat1(), VelocityVerlet,
                                                         Peridynamics.BondSystem)

    struct StorageWrong1 <: Peridynamics.AbstractStorage end
    @test_throws InterfaceError Peridynamics.point_data_fields(StorageWrong1)

    # the fracture bookkeeping belongs to the damage model and is checked when a Job is
    # created, see `check_damage_model`, so nothing about it is checked at expansion time

    # an untyped field is rejected when the macro is expanded
    try
        eval(quote
            @storage Mat1 struct StorageUntyped1
                @inherit VelocityVerletFields BondFracFields
                mycustomfield
            end
        end)
        @test false
    catch e
        @test isa(e, LoadError)
    end

    # a storage for one solver ...
    @storage Mat1 VelocityVerlet struct Storage1 <: Peridynamics.AbstractStorage
        @inherit VelocityVerletFields BondFracFields
        mycustomfield::PointScalar
    end

    @test hasmethod(Peridynamics.storage_type, Tuple{Mat1})
    @test hasmethod(Storage1, Tuple{Mat1,VelocityVerlet,Peridynamics.AbstractSystem})
    @test Peridynamics.storage_type(Mat1()) <: Storage1

    @test hasmethod(Peridynamics.loc_to_halo_fields, Tuple{Storage1})
    @test hasmethod(Peridynamics.is_halo_field, Tuple{Storage1,Val{:position}})
    @test hasmethod(Peridynamics.halo_to_loc_fields, Tuple{Storage1})
    @test Peridynamics.point_data_fields(Peridynamics.storage_type(Mat1())) ==
          (:position, :displacement, :velocity, :velocity_half, :acceleration, :b_int,
           :b_ext, :damage, :n_active_bonds, :mycustomfield)

    # ... and a storage for every solver
    @storage Mat1 struct Storage2 <: Peridynamics.AbstractStorage
        @inherit VelocityVerletFields BondFracFields
        mycustomfield::PointScalar
    end

    @test hasmethod(Storage2, Tuple{Mat1,AbstractTimeSolver,AbstractSystem})
    @test Peridynamics.storage_type(Mat1()) <: Storage2
end

@testitem "@storage: a custom material with a container-typed field and @halo_fields" begin
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection,
                         AbstractInteractionSystemMaterial, InterfaceError,
                         AbstractPointParameters, AbstractDamageModel

    struct Mat3{DM} <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::DM
        function Mat3(dmgmodel::DM) where DM
            new{DM}(dmgmodel)
        end
    end
    Mat3(; dmgmodel::AbstractDamageModel=CriticalStretch()) = Mat3(dmgmodel)
    Peridynamics.@params Mat3 struct Params3
        @inherit StandardParameters
    end

    pos, vol = uniform_box(1,1,1,0.4)
    mat = Mat3()
    body = Body(mat, pos, vol)
    material!(body, horizon=1, rho=1, E=1, nu=0.25, Gc=1)
    pd = Peridynamics.PointDecomposition(body, 1)
    system = Peridynamics.get_system(body, pd, 1)
    solver = VelocityVerlet(steps=1)

    @test_throws InterfaceError Peridynamics.storage_type(mat)
    @test_throws InterfaceError Peridynamics.get_storage(mat, solver, system)

    # `myfld` is declared with a container type, so it needs an `init_field` method and is
    # not point data
    Peridynamics.@storage Mat3 struct Storage3
        @inherit VelocityVerletFields DynamicRelaxationFields BondFracFields
        myfld::Matrix{Float64}
    end

    @test Peridynamics.storage_type(mat) <: Storage3

    pointfields = (:position, :displacement, :velocity, :velocity_half, :acceleration, :b_int,
                   :b_ext, :velocity_half_old, :b_int_old, :density_matrix, :damage,
                   :n_active_bonds)
    @test Peridynamics.point_data_fields(Storage3) === pointfields

    @test_throws InterfaceError Peridynamics.init_field(mat, solver, system, Val(:myfld))

    function Peridynamics.init_field(::Mat3, ::Peridynamics.AbstractTimeSolver,
                                     system::Peridynamics.AbstractSystem, ::Val{:myfld})
        return zeros(3, Peridynamics.get_n_points(system))
    end
    @test Peridynamics.init_field(mat, solver, system, Val(:myfld)) ≈ zeros(3, 8)

    storage = Peridynamics.get_storage(mat, solver, system)
    @test Peridynamics.loc_to_halo_fields(storage) === (:position,)
    @test Peridynamics.halo_to_loc_fields(storage) == ()

    @test Peridynamics.get_loc_to_halo_fields(storage) == (storage.position,)
    @test Peridynamics.get_halo_to_loc_fields(storage) == ()

    @test Peridynamics.is_halo_field(storage, Val(:position)) == true
    for field in (:displacement, :velocity, :velocity_half, :velocity_half_old,
                  :acceleration, :b_int, :b_int_old, :b_ext, :density_matrix, :damage,
                  :n_active_bonds, :myfld)
        @test Peridynamics.is_halo_field(storage, Val(field)) == false
    end

    # the fields of the other solver are empty arrays of the right type
    @test size(storage.velocity_half_old) == (0, 0)
    @test storage.velocity_half_old isa Matrix{Float64}

    Peridynamics.@halo_fields Storage3 :myfld

    @test Peridynamics.is_halo_field(storage, Val(:myfld)) == true

    @test_throws InterfaceError Peridynamics.point_data_field(storage, Val(:bond_active))
end

@testitem "typecheck_storage: fallbacks and missing fields" begin
    import Peridynamics: typecheck_storage, typecheck_is_storage, typecheck_storage_fields,
                         typecheck_req_fields_missing, required_fields, AbstractMaterial
    struct NotAStorage end

    @test_throws ArgumentError typecheck_storage(BBMaterial, NotAStorage)
    @test_throws ArgumentError typecheck_storage(BBMaterial, NotAStorage())
    @test_throws ArgumentError typecheck_is_storage(NotAStorage)
    @test_throws ArgumentError typecheck_is_storage(NotAStorage())
    @test typecheck_is_storage(Peridynamics.BBStorage) === nothing

    @test typecheck_storage_fields(Peridynamics.BBStorage, (:position, :b_int)) === nothing
    @test_throws ArgumentError typecheck_storage_fields(Peridynamics.BBStorage,
                                                        (:position, :not_a_field))

    # a material family can require fields from its type alone, and a storage that does not
    # declare them is rejected while `@storage` is expanded
    struct ReqFieldsMat <: AbstractMaterial end
    Peridynamics.required_fields(::Type{ReqFieldsMat}) = (:position, :my_extra_field)
    @test required_fields(ReqFieldsMat) === (:position, :my_extra_field)
    @test typecheck_req_fields_missing(Peridynamics.BBStorage, (:position, :b_int)) === false
    err = try
        typecheck_req_fields_missing(Peridynamics.BBStorage, (:my_extra_field,))
    catch e
        e
    end
    @test err isa ErrorException
    @test contains(err.msg, "required field my_extra_field not found")
    @test_throws ErrorException typecheck_storage(ReqFieldsMat, Peridynamics.BBStorage)
end

@testitem "storage interface: halo field fallbacks and local point data" setup=[Fixtures] begin
    struct BareStorage <: Peridynamics.AbstractStorage end
    @test_throws Peridynamics.InterfaceError Peridynamics.loc_to_halo_fields(BareStorage())
    @test_throws Peridynamics.InterfaceError Peridynamics.halo_to_loc_fields(BareStorage())
    @test Peridynamics.is_halo_field(BareStorage(), Val(:position)) == false

    # a halo field holds the halo points as well and its local part is a view on the local
    # points; a field without halo exchange holds only the local points and is returned as is
    body = Fixtures.line10()
    c = Fixtures.chunk(body; n_chunks=2, chunk_id=1)
    n_loc = Peridynamics.get_n_loc_points(c.system)
    n_all = Peridynamics.get_n_points(c.system)
    @test n_all > n_loc
    @test Peridynamics.is_halo_field(c.storage, Val(:position))
    @test size(c.storage.position) == (3, n_all)
    loc_position = Peridynamics.get_loc_point_data(c.storage, c.system, :position)
    @test loc_position isa SubArray && size(loc_position) == (3, n_loc)
    @test !Peridynamics.is_halo_field(c.storage, Val(:velocity))
    @test Peridynamics.get_loc_point_data(c.storage, c.system, :velocity) === c.storage.velocity
    @test size(c.storage.velocity) == (3, n_loc)
end

@testitem "@storage: the expansion of a storage that carries nested states" begin
    # the states of the constitutive model and of the damage model are not allocated by
    # `init_field` but by the models, and they are the trailing type parameters
    ex = @macroexpand Peridynamics.@storage BBMaterial struct NestedStateStorage
        @inherit Peridynamics.VelocityVerletFields Peridynamics.BondFracFields
        cm_state::ConstitutiveState
        dmg_state::DamageState
    end
    # (the names of the arguments are mangled by the macro hygiene, the functions are not)
    s = string(ex)
    @test contains(s, ").init_constitutive_state(")
    @test contains(s, ").init_damage_state(")
    @test occursin(r"\)\.constitutive_storage_type\(\(.*\)\.get_constitutive_model\(", s)
    @test occursin(r"\)\.damage_storage_type\(\(.*\)\.get_dmgmodel\(", s)
    @test contains(s, ").constitutive_state(")
    @test contains(s, ").damage_state(")
    @test contains(s, ").has_constitutive_state(")
    @test contains(s, ").has_damage_state(")
    @test occursin(r"#CMS\", var\"#\d+#DMS\"}", s)

    # each state at most once and never annotated for halo exchange, see
    # `derive_storage_type_params`
    @test_throws LoadError @eval Peridynamics.@storage BBMaterial struct TwiceState
        @inherit Peridynamics.VelocityVerletFields
        a::ConstitutiveState
        b::ConstitutiveState
    end
    @test_throws LoadError @eval Peridynamics.@storage BBMaterial struct HaloState
        @inherit Peridynamics.VelocityVerletFields
        @htl a::DamageState
    end
end

@testitem "@storage: a storage without derived type parameters" begin
    import Peridynamics: AbstractMaterial, AbstractTimeSolver, AbstractSystem, storage_type,
                         get_storage
    using Peridynamics.StaticArrays

    # a storage whose fields are all declared with types that are not arrays has nothing
    # that could follow the float type or be moved to another array backend
    struct StaticMat <: AbstractMaterial end
    Peridynamics.@storage StaticMat struct StaticStorage
        tensor::MArray{Tuple{3,3},Float64,2,9}
        counter::Int
    end
    function Peridynamics.init_field(::StaticMat, ::AbstractTimeSolver, ::AbstractSystem,
                                     ::Val{:tensor})
        return zero(MArray{Tuple{3,3},Float64,2,9})
    end
    Peridynamics.init_field(::StaticMat, ::AbstractTimeSolver, ::AbstractSystem, ::Val{:counter}) = 0

    # the number of spatial dimensions is the one parameter every storage carries
    @test StaticStorage isa UnionAll
    @test storage_type(StaticMat()) === StaticStorage{3}
    @test storage_type(StaticMat(), Float32) === StaticStorage{3}
    @test storage_type(StaticMat(), Float64, Val(2)) === StaticStorage{2}

    position = zeros(3, 4)
    position[1, :] = 0.0:3.0
    body = Body(BBMaterial(), position, ones(4))
    material!(body; horizon=1.5, rho=1.0, E=1.0, Gc=1.0)
    pd = Peridynamics.PointDecomposition(body, 1)
    ps = Peridynamics.get_param_spec(body)
    system = Peridynamics.BodyChunk(body, VelocityVerlet(steps=1), pd, 1, ps).system
    s = get_storage(StaticMat(), VelocityVerlet(steps=1), system)
    @test s isa StaticStorage{3}
    @test Peridynamics.get_n_dim(s) == Peridynamics.get_n_dim(system)
    @test iszero(s.tensor) && s.counter == 0
    # no `Adapt.adapt_structure` method is generated, so `adapt` passes the storage through
    @test Peridynamics.Adapt.adapt(Array, s) === s
end

@testitem "storage_type / init_field_hint / storage_contract: the fallbacks" begin
    import Peridynamics: AbstractMaterial, AbstractStorage, storage_type, init_field_hint,
                         storage_contract

    # a storage that is not generated by `@storage` ignores the requested float type
    struct ManualMat <: AbstractMaterial end
    struct ManualStorage <: AbstractStorage end
    Peridynamics.storage_type(::ManualMat) = ManualStorage
    @test storage_type(ManualMat(), Float32) === ManualStorage
    # ... and the number of spatial dimensions of the simulation as well
    @test storage_type(ManualMat(), Float32, Val(2)) === ManualStorage
    @test storage_type(ManualMat(), Float64, Val(3)) === ManualStorage

    # the hint of `init_field` only knows what to say for a `Val` field
    @test contains(init_field_hint(Val(:my_field)), "my_field::PointScalar")
    @test init_field_hint(:my_field) == ""

    # a material without a damage model has a contract whose damage part is `Nothing`
    @test isnothing(Peridynamics.get_dmgmodel(ManualMat()))
    contract = storage_contract(ManualMat(), VelocityVerlet(steps=1))
    @test length(contract) == 3
    @test contract[2] == ((), "the damage model `Nothing`")
    @test contract[1] == ((), "the material `ManualMat`")
end

@testitem "storage property forwarding: flat reads reach into the nested states" setup=[Fixtures] begin
    # `storage.bond_active` reads a field that lives in the state of the damage model, so
    # a kernel never sees the nesting; a flat field is read exactly as before
    body = Fixtures.cube(BBMaterial())
    storage = Fixtures.chunk(body).storage
    state = Peridynamics.damage_state(storage)
    @test storage.bond_active === state.bond_active
    @test storage.damage === state.damage
    @test storage.b_int === Base.getfield(storage, :b_int)

    # destructuring goes through `getproperty` and works for both kinds
    (; n_active_bonds, b_int) = storage
    @test n_active_bonds === state.n_active_bonds
    @test b_int === Base.getfield(storage, :b_int)

    # `propertynames` lists the fields plus what the states hold
    names = propertynames(storage)
    @test :bond_active in names
    @test :damage in names
    @test :b_int in names
    @test :dmg_state in names

    # an unknown name falls through to `getfield` and its native error
    @test_throws Exception storage.not_a_field

    # an ambiguous name reports both homes and how to read it directly
    err = try
        Peridynamics.ambiguous_storage_property(:damage, storage, (:cm_state, :dmg_state))
    catch e
        e
    end
    @test err isa ArgumentError
    @test contains(err.msg, "`damage` exists in `cm_state` and `dmg_state`")
    @test contains(err.msg, "storage.cm_state.damage")
end

@testitem "check_state_field_collisions: a flat field must not shadow a state field" begin
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection, check_storage_contract

    # a storage that declares the bookkeeping flat although the damage model carries it
    struct CollMat{D} <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::D
    end
    CollMat() = CollMat(CriticalStretch())
    Peridynamics.@params CollMat struct CollParams
        @inherit StandardParameters
    end
    Peridynamics.@storage CollMat struct CollStorage
        @inherit VelocityVerletFields BondFracFields
        dmg_state::DamageState
    end
    err = try
        check_storage_contract(CollMat(), VelocityVerlet(steps=1))
    catch e
        e
    end
    @test err isa ArgumentError
    @test contains(err.msg, "CollStorage")
    @test contains(err.msg, "Remove the flat declaration")
end

@testitem "storage_type: the nested damage state follows the dimension of the storage" begin
    import Peridynamics: storage_type, system_type, damage_storage_type, get_n_dim,
                         get_dmgmodel, BondSystem, InteractionSystem

    # `storage_type` asks `system_type` for the system the damage state is declared for, and
    # it asks with the same `FT` and `N` it was asked with itself, so a two-dimensional
    # storage can never carry a three-dimensional damage state
    for mat in (BBMaterial(), OSBMaterial(), CMaterial(), CKIMaterial())
        S3 = storage_type(mat, Float64, Val(3))
        S2 = storage_type(mat, Float64, Val(2))
        @test S3.parameters[1] === 3
        @test S2.parameters[1] === 2
        @test fieldtype(S3, :dmg_state) ===
              damage_storage_type(get_dmgmodel(mat), system_type(mat, Float64, Val(3)),
                                  Float64)
        @test fieldtype(S2, :dmg_state).parameters[1] === 2
        @test fieldtype(S3, :dmg_state).parameters[1] === 3
    end

    # the float type reaches the nested state as well
    S32 = storage_type(BBMaterial(), Float32, Val(2))
    @test fieldtype(S32, :dmg_state).parameters[2] === Float32
end

@testitem "macrocheck_input_system: the system a damage state is declared for" begin
    import Peridynamics: macrocheck_input_system, damage_storage_type, system_type,
                         AbstractDamageModel, BondSystem

    # the system is named, either bare or qualified with the module it lives in
    @test isnothing(macrocheck_input_system(:BondSystem))
    @test isnothing(macrocheck_input_system(:(Peridynamics.BondSystem)))

    # everything that cannot name a type is rejected
    @test_throws ArgumentError macrocheck_input_system(1)
    err = try
        macrocheck_input_system(:(system_type(mat)))
    catch e
        e
    end
    @test err isa ArgumentError
    @test contains(err.msg, "is not a valid system input")

    # both accepted forms expand, and the state is declared for that system family only
    struct QualifiedDmg <: AbstractDamageModel end
    Peridynamics.@dmg_storage QualifiedDmg Peridynamics.BondSystem struct QualifiedDmgState
        marker::BondScalar{Int}
    end
    struct BareDmg <: AbstractDamageModel end
    Peridynamics.@dmg_storage BareDmg BondSystem struct BareDmgState
        marker::BondScalar{Int}
    end
    for (model, State) in ((QualifiedDmg(), QualifiedDmgState), (BareDmg(), BareDmgState))
        @test damage_storage_type(model, system_type(BBMaterial())) === State{3,Vector{Int}}
        @test damage_storage_type(model, system_type(CKIMaterial())) === Nothing
    end
end

@testitem "@dmg_storage: a damage state without derived type parameters" begin
    import Peridynamics: AbstractDamageModel, AbstractDamageState, damage_storage_type,
                         system_type, get_n_dim, Adapt

    # a state without fields has nothing that could follow the float type or be moved to
    # another array backend, so the dimension is the only parameter it carries and no
    # `Adapt.adapt_structure` is generated for it
    struct MarkerDmg <: AbstractDamageModel end
    Peridynamics.@dmg_storage MarkerDmg struct MarkerDmgState end

    @test MarkerDmgState isa UnionAll
    @test MarkerDmgState <: AbstractDamageState
    @test damage_storage_type(MarkerDmg(), system_type(BBMaterial())) === MarkerDmgState{3}
    state = MarkerDmgState{3}()
    @test get_n_dim(state) == 3
    @test Adapt.adapt(Array, state) === state
end

@testitem "check_state_field_collisions: two nested states must not carry the same field" begin
    import Peridynamics: check_state_field_collisions, storage_type,
                         AbstractConstitutiveState, AbstractDamageState, AbstractStorage

    # each state is free to name its own fields, but a name that both of them carry has no
    # single home on the storage that reads them flat
    struct TwinCMState <: AbstractConstitutiveState
        damage::Vector{Float64}
    end
    struct TwinDmgState <: AbstractDamageState
        damage::Vector{Float64}
    end
    struct TwinStorage <: AbstractStorage
        cm_state::TwinCMState
        dmg_state::TwinDmgState
    end
    err = try
        check_state_field_collisions(TwinStorage)
    catch e
        e
    end
    @test err isa ArgumentError
    @test contains(err.msg, "both nested states of the storage `TwinStorage`")
    @test contains(err.msg, "`damage`")
    @test contains(err.msg, "rename the field in one of the two models")

    # a storage whose states have nothing in common passes
    @test isnothing(check_state_field_collisions(storage_type(BBMaterial())))
end

@testitem "get_storage_property: a name that both nested states carry is ambiguous" begin
    import Peridynamics: AbstractConstitutiveState, AbstractDamageState, AbstractStorage

    # the backstop of `check_state_field_collisions` for a storage that is built without a
    # `Job`: the flat read cannot decide which state the name belongs to
    struct AmbiCMState <: AbstractConstitutiveState
        damage::Vector{Float64}
        plastic_strain::Vector{Float64}
    end
    struct AmbiDmgState <: AbstractDamageState
        damage::Vector{Float64}
    end
    struct AmbiStorage <: AbstractStorage
        cm_state::AmbiCMState
        dmg_state::AmbiDmgState
    end
    s = AmbiStorage(AmbiCMState([0.0], [1.0]), AmbiDmgState([2.0]))

    # a flat field and a name that only one state carries still read
    @test s.cm_state === Base.getfield(s, :cm_state)
    @test s.plastic_strain == [1.0]

    err = try
        s.damage
    catch e
        e
    end
    @test err isa ArgumentError
    @test contains(err.msg, "`damage` exists in `cm_state` and `dmg_state`")
end

@testitem "point_data_field: a damage model without a state serves no point data" begin
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection, AbstractDamageModel,
                         point_data_field, point_data_fields, nested_point_data_fields,
                         damage_state, storage_type, get_storage

    # `Nothing` is the state of a damage model that declares none: no point data at all,
    # and a read that reaches it says why the field cannot be served
    @test point_data_fields(Nothing) === ()
    @test nested_point_data_fields(Nothing) === ()
    err = try
        point_data_field(nothing, Val(:damage))
    catch e
        e
    end
    @test err isa ArgumentError
    @test contains(err.msg, "no point data field `damage`")

    # a storage that declares `dmg_state` although its model has no state reaches exactly
    # those methods, through the forwarding the macro generates for the marker
    struct StatelessDmg <: AbstractDamageModel end
    struct StatelessDmgMat{D} <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::D
    end
    StatelessDmgMat() = StatelessDmgMat(StatelessDmg())
    Peridynamics.@params StatelessDmgMat struct StatelessDmgParams
        @inherit StandardParameters
    end
    Peridynamics.@storage StatelessDmgMat struct StatelessDmgStorage
        @inherit VelocityVerletFields
        dmg_state::DamageState
    end

    mat = StatelessDmgMat()
    @test fieldtype(storage_type(mat), :dmg_state) === Nothing
    @test point_data_fields(storage_type(mat)) ===
          (:position, :displacement, :velocity, :velocity_half, :acceleration, :b_int,
           :b_ext)

    pos, vol = uniform_box(1, 1, 1, 0.5)
    body = Body(mat, pos, vol)
    material!(body; horizon=0.8, rho=1, E=1, nu=0.25)
    pd = Peridynamics.PointDecomposition(body, 1)
    system = Peridynamics.get_system(body, pd, 1)
    storage = get_storage(mat, VelocityVerlet(steps=1), system)
    @test isnothing(damage_state(storage))
    @test point_data_field(storage, Val(:position)) === Base.getfield(storage, :position)
    @test_throws ArgumentError point_data_field(storage, Val(:damage))
end
