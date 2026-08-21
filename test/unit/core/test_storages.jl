# The storage contract and the `@storage` macro of `src/core/storages.jl`. The field
# declaration framework (shapes, blocks, allocation) is tested in `test_storage_fields.jl`.

@testitem "required_fields: the type-level part of the storage contract" begin
    # only the solver-independent, type-level part of the contract is checked by `@storage`
    @test Peridynamics.required_fields(Peridynamics.AbstractMaterial) === ()

    rf_bb = (:damage, :n_active_bonds, :bond_active)
    @test Peridynamics.required_fields(BBMaterial) === rf_bb

    rf_cki = (:damage, :n_active_one_nis, :one_ni_active)
    @test Peridynamics.required_fields(CKIMaterial) === rf_cki
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
    # damage model must not error
    @test Peridynamics.req_storage_fields(BBMaterial(), CriticalStretch()) === ()
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
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection, StandardPointParameters,
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
    Peridynamics.@params ContractMat StandardPointParameters
    Peridynamics.@storage ContractMat struct ContractStorage
        @inherit VelocityVerletFields BondFracFields
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

    # the fields of the system are checked when the macro is expanded: `n_active_bonds` and
    # `n_active_one_nis` are missing
    @test_throws ErrorException @storage Mat1 struct StorageMissing1
        @inherit VelocityVerletFields
        damage::PointScalar
        bond_active::BondScalar{Bool}
    end

    @test_throws ErrorException @storage Mat2 struct StorageMissing2
        @inherit VelocityVerletFields
        damage::PointScalar
        one_ni_active::Vector{Bool}
    end

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
    struct Params3 <: AbstractPointParameters
        δ::Float64
        rho::Float64
        E::Float64
        nu::Float64
        G::Float64
        K::Float64
        λ::Float64
        μ::Float64
        Gc::Float64
        εc::Float64
        bc::Float64
    end
    @test_throws InterfaceError Peridynamics.@params Mat3 Params3
    function Params3(mat::Mat3, p::Dict{Symbol,Any})
        (; δ, rho, E, nu, G, K, λ, μ) = Peridynamics.get_required_point_parameters(mat, p)
        (; Gc, εc) = Peridynamics.get_frac_params(mat.dmgmodel, p, δ, K)
        bc = 18 * K / (π * δ^4) # bond constant
        return Params3(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
    end
    Peridynamics.@params Mat3 Params3

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
    import Peridynamics: typecheck_storage, typecheck_is_storage, typecheck_storage_fields
    struct NotAStorage end

    @test_throws ArgumentError typecheck_storage(BBMaterial, NotAStorage)
    @test_throws ArgumentError typecheck_storage(BBMaterial, NotAStorage())
    @test_throws ArgumentError typecheck_is_storage(NotAStorage)
    @test_throws ArgumentError typecheck_is_storage(NotAStorage())
    @test typecheck_is_storage(Peridynamics.BBStorage) === nothing

    @test typecheck_storage_fields(Peridynamics.BBStorage, (:position, :b_int)) === nothing
    @test_throws ArgumentError typecheck_storage_fields(Peridynamics.BBStorage,
                                                        (:position, :not_a_field))
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
