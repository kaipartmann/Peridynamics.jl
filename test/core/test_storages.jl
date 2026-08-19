@testitem "required_fields" begin
    rf_abstractmat = (:position, :displacement, :velocity, :velocity_half, :acceleration,
                      :b_int, :b_ext, :velocity_half_old, :b_int_old, :density_matrix)
    @test Peridynamics.required_fields(Peridynamics.AbstractMaterial) === rf_abstractmat

    rf_bb = (:position, :displacement, :velocity, :velocity_half, :acceleration, :b_int,
             :b_ext, :velocity_half_old, :b_int_old, :density_matrix, :damage,
             :n_active_bonds, :bond_active)
    @test Peridynamics.required_fields(BBMaterial) === rf_bb

    rf_cki = (:position, :displacement, :velocity, :velocity_half, :acceleration, :b_int,
              :b_ext, :velocity_half_old, :b_int_old, :density_matrix, :damage,
              :n_active_one_nis, :one_ni_active)
    @test Peridynamics.required_fields(CKIMaterial) === rf_cki
end

@testitem "get_storage_header" begin
    import Peridynamics: get_storage_header

    # Test case 1: Simple struct declaration
    expr1 = :(struct MyStorage end)
    header1, type1 = get_storage_header(expr1)
    @test header1 == Expr(:(<:), :MyStorage, :(Peridynamics.AbstractStorage))
    @test type1 == :MyStorage

    # Test case 2: Struct declaration with a subtype
    expr2 = :(struct MyStorage <: Peridynamics.AbstractStorage end)
    header2, type2 = get_storage_header(expr2)
    @test header2 == Expr(:(<:), :MyStorage, :(Peridynamics.AbstractStorage))
    @test type2 == :MyStorage

    # Test case 3: Parametric struct declaration with a subtype
    expr3 = :(struct MyStorage{A,B,C} <: Peridynamics.AbstractStorage end)
    header3, type3 = get_storage_header(expr3)
    @test header3 == Expr(:(<:), :(MyStorage{A,B,C}), :(Peridynamics.AbstractStorage))
    @test type3 == :MyStorage

    # Test case 4: Parametric struct declaration without a subtype
    expr4 = :(struct MyStorage{A,B,C} end)
    header4, type4 = get_storage_header(expr4)
    @test header4 == Expr(:(<:), :(MyStorage{A,B,C}), :(Peridynamics.AbstractStorage))
    @test type4 == :MyStorage
end

@testitem "macrochecks @storage macro" begin

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

@testitem "@storage: storage declaration" begin
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

    @test_throws ErrorException @storage Mat1 struct StorageMissing1
        @lthfield position::Matrix{Float64}
        # @pointfield displacement::Matrix{Float64}
        @pointfield velocity::Matrix{Float64}
        @pointfield velocity_half::Matrix{Float64}
        @pointfield velocity_half_old::Matrix{Float64}
        @pointfield acceleration::Matrix{Float64}
        @pointfield b_int::Matrix{Float64}
        @pointfield b_int_old::Matrix{Float64}
        @pointfield b_ext::Matrix{Float64}
        @pointfield density_matrix::Matrix{Float64}
        @pointfield damage::Vector{Float64}
        bond_active::Vector{Bool}
        @pointfield n_active_bonds::Vector{Int}
    end

    @test_throws ErrorException @storage Mat1 struct StorageMissing2
        @lthfield position::Matrix{Float64}
        @pointfield displacement::Matrix{Float64}
        @pointfield velocity::Matrix{Float64}
        @pointfield velocity_half::Matrix{Float64}
        @pointfield velocity_half_old::Matrix{Float64}
        @pointfield acceleration::Matrix{Float64}
        @pointfield b_int::Matrix{Float64}
        @pointfield b_int_old::Matrix{Float64}
        @pointfield b_ext::Matrix{Float64}
        @pointfield density_matrix::Matrix{Float64}
        @pointfield damage::Vector{Float64}
        bond_active::Vector{Bool}
        # @pointfield n_active_bonds::Vector{Int}
    end

    @test_throws ErrorException @storage Mat2 struct StorageMissing3
        @lthfield position::Matrix{Float64}
        @pointfield displacement::Matrix{Float64}
        @pointfield velocity::Matrix{Float64}
        @pointfield velocity_half::Matrix{Float64}
        @pointfield velocity_half_old::Matrix{Float64}
        @pointfield acceleration::Matrix{Float64}
        @pointfield b_int::Matrix{Float64}
        @pointfield b_int_old::Matrix{Float64}
        @pointfield b_ext::Matrix{Float64}
        @pointfield density_matrix::Matrix{Float64}
        @pointfield damage::Vector{Float64}
        one_ni_active::Vector{Bool}
        # @pointfield n_active_one_nis::Vector{Int}
    end

    try
        eval(quote
            @storage Mat1 struct StorageUntyped1
                @lthfield position
                @pointfield displacement::Matrix{Float64}
                @pointfield velocity::Matrix{Float64}
                @pointfield velocity_half::Matrix{Float64}
                @pointfield velocity_half_old::Matrix{Float64}
                @pointfield acceleration::Matrix{Float64}
                @pointfield b_int::Matrix{Float64}
                @pointfield b_int_old::Matrix{Float64}
                @pointfield b_ext::Matrix{Float64}
                @pointfield density_matrix::Matrix{Float64}
                @pointfield damage::Vector{Float64}
                bond_active::Vector{Bool}
                @pointfield n_active_bonds::Vector{Int}
                @pointfield mycustomfield
            end
        end)
        @test false
    catch e
        @test isa(e, LoadError)
    end

    try
        eval(quote
            @storage Mat1 struct StorageUntyped1
                @lthfield position
                @pointfield displacement::Matrix{Float64}
                @pointfield velocity::Matrix{Float64}
                @pointfield velocity_half::Matrix{Float64}
                @pointfield velocity_half_old::Matrix{Float64}
                @pointfield acceleration::Matrix{Float64}
                @pointfield b_int::Matrix{Float64}
                @pointfield b_int_old::Matrix{Float64}
                @pointfield b_ext::Matrix{Float64}
                @pointfield density_matrix::Matrix{Float64}
                @pointfield damage::Vector{Float64}
                bond_active::Vector{Bool}
                @pointfield n_active_bonds::Vector{Int}
                @lthfield mycustomfield
            end
        end)
        @test false
    catch e
        @test isa(e, LoadError)
    end

    try
        eval(quote
            @storage Mat1 struct StorageUntyped1
                @lthfield position
                @pointfield displacement::Matrix{Float64}
                @pointfield velocity::Matrix{Float64}
                @pointfield velocity_half::Matrix{Float64}
                @pointfield velocity_half_old::Matrix{Float64}
                @pointfield acceleration::Matrix{Float64}
                @pointfield b_int::Matrix{Float64}
                @pointfield b_int_old::Matrix{Float64}
                @pointfield b_ext::Matrix{Float64}
                @pointfield density_matrix::Matrix{Float64}
                @pointfield damage::Vector{Float64}
                bond_active::Vector{Bool}
                @pointfield n_active_bonds::Vector{Int}
                @htlfield mycustomfield
            end
        end)
        @test false
    catch e
        @test isa(e, LoadError)
    end

    try
        eval(quote
            @storage Mat1 struct StorageUntyped1
                @lthfield position
                @pointfield displacement::Matrix{Float64}
                @pointfield velocity::Matrix{Float64}
                @pointfield velocity_half::Matrix{Float64}
                @pointfield velocity_half_old::Matrix{Float64}
                @pointfield acceleration::Matrix{Float64}
                @pointfield b_int::Matrix{Float64}
                @pointfield b_int_old::Matrix{Float64}
                @pointfield b_ext::Matrix{Float64}
                @pointfield density_matrix::Matrix{Float64}
                @pointfield damage::Vector{Float64}
                bond_active::Vector{Bool}
                @pointfield n_active_bonds::Vector{Int}
                mycustomfield
            end
        end)
        @test false
    catch e
        @test isa(e, LoadError)
    end

    @storage Mat1 VelocityVerlet struct Storage1 <: Peridynamics.AbstractStorage
        @lthfield position::Matrix{Float64}
        @pointfield displacement::Matrix{Float64}
        @pointfield velocity::Matrix{Float64}
        @pointfield velocity_half::Matrix{Float64}
        @pointfield velocity_half_old::Matrix{Float64}
        @pointfield acceleration::Matrix{Float64}
        @pointfield b_int::Matrix{Float64}
        @pointfield b_int_old::Matrix{Float64}
        @pointfield b_ext::Matrix{Float64}
        @pointfield density_matrix::Matrix{Float64}
        @pointfield damage::Vector{Float64}
        bond_active::Vector{Bool}
        @pointfield n_active_bonds::Vector{Int}
        @pointfield mycustomfield::Vector{Float64}
    end

    @test hasmethod(Peridynamics.storage_type, Tuple{Mat1})
    @test hasmethod(Storage1, Tuple{Mat1,VelocityVerlet,Peridynamics.AbstractSystem})
    mat, vv = Mat1(), VelocityVerlet(steps=1)
    @test Peridynamics.storage_type(mat) == Storage1

    @test hasmethod(Peridynamics.loc_to_halo_fields, Tuple{Storage1})
    @test hasmethod(Peridynamics.is_halo_field, Tuple{Storage1,Val{:position}})
    @test hasmethod(Peridynamics.is_halo_field, Tuple{Storage1,Val{:displacement}})

    @test hasmethod(Peridynamics.halo_to_loc_fields, Tuple{Storage1})
    @test hasmethod(Peridynamics.is_halo_field, Tuple{Storage1,Val{:b_int}})
    @test hasmethod(Peridynamics.is_halo_field, Tuple{Storage1,Val{:b_ext}})

    @storage Mat1 struct Storage2 <: Peridynamics.AbstractStorage
        @lthfield position::Matrix{Float64}
        @pointfield displacement::Matrix{Float64}
        @pointfield velocity::Matrix{Float64}
        @pointfield velocity_half::Matrix{Float64}
        @pointfield velocity_half_old::Matrix{Float64}
        @pointfield acceleration::Matrix{Float64}
        @pointfield b_int::Matrix{Float64}
        @pointfield b_int_old::Matrix{Float64}
        @pointfield b_ext::Matrix{Float64}
        @pointfield density_matrix::Matrix{Float64}
        @pointfield damage::Vector{Float64}
        bond_active::Vector{Bool}
        @pointfield n_active_bonds::Vector{Int}
        @pointfield mycustomfield::Vector{Float64}
    end

    @test hasmethod(Peridynamics.storage_type, Tuple{Mat1})
    @test hasmethod(Storage2, Tuple{Mat1,AbstractTimeSolver,AbstractSystem})
    mat, vv = Mat1(), VelocityVerlet(steps=1)
    @test Peridynamics.storage_type(mat) == Storage2

    @test hasmethod(Peridynamics.loc_to_halo_fields, Tuple{Storage2})
    @test hasmethod(Peridynamics.is_halo_field, Tuple{Storage2,Val{:position}})
    @test hasmethod(Peridynamics.is_halo_field, Tuple{Storage2,Val{:displacement}})

    @test hasmethod(Peridynamics.halo_to_loc_fields, Tuple{Storage2})
    @test hasmethod(Peridynamics.is_halo_field, Tuple{Storage2,Val{:b_int}})
    @test hasmethod(Peridynamics.is_halo_field, Tuple{Storage2,Val{:b_ext}})
end

@testitem "@storage: custom materials and halo fields" begin
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

    Peridynamics.@storage Mat3 struct Storage3
        @lthfield position::Matrix{Float64}
        @pointfield displacement::Matrix{Float64}
        @pointfield velocity::Matrix{Float64}
        @pointfield velocity_half::Matrix{Float64}
        @pointfield velocity_half_old::Matrix{Float64}
        @pointfield acceleration::Matrix{Float64}
        @pointfield b_int::Matrix{Float64}
        @pointfield b_int_old::Matrix{Float64}
        @pointfield b_ext::Matrix{Float64}
        @pointfield density_matrix::Matrix{Float64}
        @pointfield damage::Vector{Float64}
        bond_active::Vector{Bool}
        @pointfield n_active_bonds::Vector{Int}
        @pointfield myfld::Matrix{Float64}
    end

    @test Peridynamics.storage_type(mat) == Storage3

    pointfields = (:position, :displacement, :velocity, :velocity_half, :velocity_half_old,
                   :acceleration, :b_int, :b_int_old, :b_ext, :density_matrix, :damage,
                   :n_active_bonds, :myfld)
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
    @test Peridynamics.is_halo_field(storage, Val(:displacement)) == false
    @test Peridynamics.is_halo_field(storage, Val(:velocity)) == false
    @test Peridynamics.is_halo_field(storage, Val(:velocity_half)) == false
    @test Peridynamics.is_halo_field(storage, Val(:velocity_half_old)) == false
    @test Peridynamics.is_halo_field(storage, Val(:acceleration)) == false
    @test Peridynamics.is_halo_field(storage, Val(:b_int)) == false
    @test Peridynamics.is_halo_field(storage, Val(:b_int_old)) == false
    @test Peridynamics.is_halo_field(storage, Val(:b_ext)) == false
    @test Peridynamics.is_halo_field(storage, Val(:density_matrix)) == false
    @test Peridynamics.is_halo_field(storage, Val(:damage)) == false
    @test Peridynamics.is_halo_field(storage, Val(:n_active_bonds)) == false
    @test Peridynamics.is_halo_field(storage, Val(:myfld)) == false

    Peridynamics.@halo_fields Storage3 :myfld

    @test Peridynamics.is_halo_field(storage, Val(:myfld)) == true

    @test_throws InterfaceError Peridynamics.point_data_field(storage, Val(:bond_active))
end
