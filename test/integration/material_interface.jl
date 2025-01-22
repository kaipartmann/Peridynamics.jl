@testitem "Material declaration" begin
    import Peridynamics: NoCorrection, InterfaceError

    struct TestMaterial1 <: Peridynamics.AbstractBondSystemMaterial{NoCorrection} end
    @test isnothing(Peridynamics.typecheck_material(TestMaterial1))
    @test Peridynamics.required_point_parameters(TestMaterial1) === (:δ, :rho, :E, :nu, :G,
           :K, :λ, :μ)
    @test Peridynamics.allowed_material_kwargs(TestMaterial1()) === (:horizon, :rho, :E,
           :nu, :G, :K, :lambda, :mu, :Gc, :epsilon_c)

    struct WrongTestMaterial end
    @test_throws ArgumentError Peridynamics.typecheck_material(WrongTestMaterial)

    struct WrongTestMaterial2 <: Peridynamics.AbstractMaterial end
    @test isnothing(Peridynamics.typecheck_material(WrongTestMaterial2))
    @test_throws InterfaceError Peridynamics.required_point_parameters(WrongTestMaterial2)
    @test_throws InterfaceError Peridynamics.allowed_material_kwargs(WrongTestMaterial2())
end

@testitem "Point parameters declaration" begin
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection, AbstractPointParameters,
                         InterfaceError, typecheck_params, constructor_check,
                         point_param_type, material_type, get_point_params,
                         macrocheck_input_material, macrocheck_input_params
    struct TestMaterial2 <: AbstractBondSystemMaterial{NoCorrection} end
    struct TestPointParameters2 <: AbstractPointParameters
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
    end
    tpp2 = TestPointParameters2(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    TestPointParameters2(::TestMaterial2, ::Dict{Symbol,Any}) = nothing

    @test isnothing(typecheck_params(TestMaterial2, TestPointParameters2))
    ie1 = InterfaceError(TestMaterial2, "point_param_type")
    @test_throws ie1 point_param_type(TestMaterial2())
    ie2 = InterfaceError(TestPointParameters2, "material_type")
    @test_throws ie2 material_type(tpp2)
    ie3 = InterfaceError(TestMaterial2, "get_point_params")
    @test_throws ie3 get_point_params(TestMaterial2(), Dict{Symbol,Any}())

    struct PointParametersNoSubtype
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
    end
    @test_throws ArgumentError typecheck_params(TestMaterial2, PointParametersNoSubtype)
    @test_throws InterfaceError constructor_check(TestMaterial2, PointParametersNoSubtype)

    struct PointParametersMissingHorizon <: AbstractPointParameters
        rho::Float64
        E::Float64
        nu::Float64
        G::Float64
        K::Float64
        λ::Float64
        μ::Float64
        Gc::Float64
        εc::Float64
    end
    @test_throws ErrorException typecheck_params(TestMaterial2,
                                                 PointParametersMissingHorizon)

    @test_throws InterfaceError point_param_type(TestMaterial2())

    Peridynamics.@params TestMaterial2 TestPointParameters2
    @test hasmethod(point_param_type, Tuple{TestMaterial2})
    @test Peridynamics.point_param_type(TestMaterial2()) == TestPointParameters2
    @test hasmethod(get_point_params, Tuple{TestMaterial2,Dict{Symbol,Any}})

    @test isnothing(macrocheck_input_material(:MyMaterial))
    @test isnothing(macrocheck_input_material(:(MyModule.MyMaterial)))
    @test_throws ArgumentError macrocheck_input_material(:(1 + 1))
    @test isnothing(macrocheck_input_params(:MyParams))
    @test isnothing(macrocheck_input_params(:(MyModule.MyParams)))
    @test_throws ArgumentError macrocheck_input_params(:(1 + 1))
end

@testitem "Storage declaration" begin
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

@testitem "Custom Materials" begin
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection,
                         AbstractInteractionSystemMaterial, InterfaceError,
                         AbstractPointParameters

    struct Mat3 <: AbstractBondSystemMaterial{NoCorrection} end
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
        (; Gc, εc) = Peridynamics.get_frac_params(p, δ, K)
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

@testitem "TestMaterial halo exchange" begin
    include("test_material.jl")

    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = TestMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, nu=0.25, Gc=1)
    point_set!(body, :a, 1:2)
    point_set!(body, :b, 3:4)
    velocity_ic!(body, :a, :x, 1.0)
    velocity_bc!(t -> t, body, :a, :x)
    forcedensity_bc!(t -> t, body, :a, :x)
    precrack!(body, :a, :b)
    ts = VelocityVerlet(steps=10)
    dh = Peridynamics.threads_data_handler(body, ts, 2)

    b1 = dh.chunks[1]
    @test b1 isa Peridynamics.BodyChunk
    @test b1.storage.position == position
    @test b1.storage.displacement == zeros(3, 2)
    @test b1.storage.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.velocity_half == zeros(3, 2)
    @test b1.storage.acceleration == zeros(3, 2)
    @test b1.storage.b_int == zeros(3, 4)
    @test b1.storage.b_ext == zeros(3, 2)
    @test b1.storage.damage ≈ [2 / 3, 2 / 3]
    @test b1.storage.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.storage.n_active_bonds == [1, 1]
    b2 = dh.chunks[2]
    @test b2 isa Peridynamics.BodyChunk
    @test b2.storage.position == position[:, [3, 4, 1, 2]]
    @test b2.storage.displacement == zeros(3, 2)
    @test b2.storage.velocity == [0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test b2.storage.velocity_half == zeros(3, 2)
    @test b2.storage.acceleration == zeros(3, 2)
    @test b2.storage.b_int == zeros(3, 4)
    @test b2.storage.b_ext == zeros(3, 2)
    @test b2.storage.damage ≈ [2 / 3, 2 / 3]
    @test b2.storage.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.storage.n_active_bonds == [1, 1]

    randpos = rand(3, 4)
    dh.chunks[2].storage.position .= randpos
    Peridynamics.exchange_loc_to_halo!(dh, 1)

    @test b1.storage.position[:, 1:2] ≈ [0.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.position[:, 3:4] ≈ randpos[:, 1:2]
    @test b1.storage.displacement == zeros(3, 2)
    @test b1.storage.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.velocity_half == zeros(3, 2)
    @test b1.storage.acceleration == zeros(3, 2)
    @test b1.storage.b_int == zeros(3, 4)
    @test b1.storage.b_ext == zeros(3, 2)
    @test b1.storage.damage ≈ [2 / 3, 2 / 3]
    @test b1.storage.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.storage.n_active_bonds == [1, 1]

    @test b2.storage.position ≈ randpos
    @test b2.storage.displacement == zeros(3, 2)
    @test b2.storage.velocity == [0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test b2.storage.velocity_half == zeros(3, 2)
    @test b2.storage.acceleration == zeros(3, 2)
    @test b2.storage.b_int == zeros(3, 4)
    @test b2.storage.b_ext == zeros(3, 2)
    @test b2.storage.damage ≈ [2 / 3, 2 / 3]
    @test b2.storage.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.storage.n_active_bonds == [1, 1]

    randbint = rand(3, 4)
    b2.storage.b_int .= randbint
    b1.storage.b_int .+= 1
    Peridynamics.exchange_halo_to_loc!(dh, 1)

    @test b1.storage.position[:, 1:2] ≈ [0.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.position[:, 3:4] ≈ randpos[:, 1:2]
    @test b1.storage.displacement == zeros(3, 2)
    @test b1.storage.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.velocity_half == zeros(3, 2)
    @test b1.storage.acceleration == zeros(3, 2)
    @test b1.storage.b_int[:, 1:2] ≈ 1 .+ randbint[:, 3:4]
    @test b1.storage.b_int[:, 3:4] ≈ ones(3, 2)
    @test b1.storage.b_ext == zeros(3, 2)
    @test b1.storage.damage ≈ [2 / 3, 2 / 3]
    @test b1.storage.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.storage.n_active_bonds == [1, 1]

    @test b2.storage.position ≈ randpos
    @test b2.storage.displacement == zeros(3, 2)
    @test b2.storage.velocity == [0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test b2.storage.velocity_half == zeros(3, 2)
    @test b2.storage.acceleration == zeros(3, 2)
    @test b2.storage.b_int ≈ randbint
    @test b2.storage.b_ext == zeros(3, 2)
    @test b2.storage.damage ≈ [2 / 3, 2 / 3]
    @test b2.storage.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.storage.n_active_bonds == [1, 1]

    Peridynamics.exchange_halo_to_loc!(dh, 2)

    @test b1.storage.position[:, 1:2] ≈ [0.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.position[:, 3:4] ≈ randpos[:, 1:2]
    @test b1.storage.displacement == zeros(3, 2)
    @test b1.storage.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.velocity_half == zeros(3, 2)
    @test b1.storage.acceleration == zeros(3, 2)
    @test b1.storage.b_int[:, 1:2] ≈ 1 .+ randbint[:, 3:4]
    @test b1.storage.b_int[:, 3:4] ≈ ones(3, 2)
    @test b1.storage.b_ext == zeros(3, 2)
    @test b1.storage.damage ≈ [2 / 3, 2 / 3]
    @test b1.storage.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.storage.n_active_bonds == [1, 1]

    @test b2.storage.position ≈ randpos
    @test b2.storage.displacement == zeros(3, 2)
    @test b2.storage.velocity == [0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test b2.storage.velocity_half == zeros(3, 2)
    @test b2.storage.acceleration == zeros(3, 2)
    @test b2.storage.b_int[:, 1:2] ≈ 1 .+ randbint[:, 1:2]
    @test b2.storage.b_int[:, 3:4] ≈ randbint[:, 3:4]
    @test b2.storage.b_ext == zeros(3, 2)
    @test b2.storage.damage ≈ [2 / 3, 2 / 3]
    @test b2.storage.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.storage.n_active_bonds == [1, 1]
end

@testitem "TestMaterial full simulation" begin
    include("test_material.jl")
    root = joinpath(@__DIR__, "temp_testmaterial_full_simulation")
    path_vtk = joinpath(root, "vtk")

    N = 10
    l, Δx, δ, a = 1.0, 1 / N, 3.015 / N, 0.5
    pos, vol = uniform_box(l, l, 0.1l, Δx)
    ids = sortperm(pos[2, :])
    body = Body(TestMaterial(), pos[:, ids], vol[ids])
    material!(body; horizon=3.015Δx, E=2.1e5, nu=0.25, rho=8e-6, Gc=2.7)
    point_set!(p -> p[1] ≤ -l / 2 + a && 0 ≤ p[2] ≤ 2δ, body, :set_a)
    point_set!(p -> p[1] ≤ -l / 2 + a && -2δ ≤ p[2] < 0, body, :set_b)
    precrack!(body, :set_a, :set_b)
    point_set!(p -> p[2] > l / 2 - Δx, body, :set_top)
    point_set!(p -> p[2] < -l / 2 + Δx, body, :set_bottom)
    velocity_bc!(t -> -30, body, :set_bottom, :y)
    velocity_bc!(t -> 30, body, :set_top, :y)
    vv = VelocityVerlet(steps=2)
    job = Job(body, vv; path=root, freq=1)
    submit(job)

    @test isdir(path_vtk)
    vtk_files = Peridynamics.find_vtk_files(path_vtk)
    @test length(vtk_files) == 3

    rm(root; recursive=true, force=true)
end
