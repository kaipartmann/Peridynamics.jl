@testitem "Material declaration" begin
    import Peridynamics: NoCorrection, InterfaceError

    struct TestMaterial1 <: Peridynamics.AbstractBondSystemMaterial{NoCorrection} end
    @test isnothing(Peridynamics.typecheck_material(TestMaterial1))
    @test Peridynamics.required_point_parameters(TestMaterial1) === (:δ, :rho, :E, :nu, :G,
                                                                     :K, :λ, :μ)
    @test Peridynamics.allowed_material_kwargs(TestMaterial1()) === (:horizon, :rho, :E,
                                                                     :nu, :Gc, :epsilon_c)

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
    tpp2 = TestPointParameters2(0,0,0,0,0,0,0,0,0,0)
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
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection,
                         AbstractInteractionSystemMaterial

    struct Mat1 <: AbstractBondSystemMaterial{NoCorrection} end
    struct Mat2 <: AbstractInteractionSystemMaterial end

    struct StorageNoSubtype
        position::Matrix{Float64}
        displacement::Matrix{Float64}
        velocity::Matrix{Float64}
        velocity_half::Matrix{Float64}
        velocity_half_old::Matrix{Float64}
        acceleration::Matrix{Float64}
        b_int::Matrix{Float64}
        b_int_old::Matrix{Float64}
        b_ext::Matrix{Float64}
        density_matrix::Matrix{Float64}
        damage::Vector{Float64}
        bond_active::Vector{Bool}
        n_active_bonds::Vector{Int}
    end

    e1 = ArgumentError("StorageNoSubtype is not a valid storage type!\n")
    @test_throws e1 Peridynamics.@storage Mat1 StorageNoSubtype

    struct StorageMissing1 <: Peridynamics.AbstractStorage
        position::Matrix{Float64}
        # displacement::Matrix{Float64}
        velocity::Matrix{Float64}
        velocity_half::Matrix{Float64}
        velocity_half_old::Matrix{Float64}
        acceleration::Matrix{Float64}
        b_int::Matrix{Float64}
        b_int_old::Matrix{Float64}
        b_ext::Matrix{Float64}
        density_matrix::Matrix{Float64}
        damage::Vector{Float64}
        bond_active::Vector{Bool}
        n_active_bonds::Vector{Int}
    end
    # TODO: somehow this throws ArgumentError defined as fallback
    # errmsg2 = "required field displacement not found in TestVerletStorageMissingField1!"
    # e2 = ErrorException(errmsg2)
    @test_throws ArgumentError Peridynamics.@storage Mat1 StorageMissing1

    struct StorageMissing2 <: Peridynamics.AbstractStorage
        position::Matrix{Float64}
        displacement::Matrix{Float64}
        velocity::Matrix{Float64}
        velocity_half::Matrix{Float64}
        velocity_half_old::Matrix{Float64}
        acceleration::Matrix{Float64}
        b_int::Matrix{Float64}
        b_int_old::Matrix{Float64}
        b_ext::Matrix{Float64}
        density_matrix::Matrix{Float64}
        damage::Vector{Float64}
        bond_active::Vector{Bool}
        # n_active_bonds::Vector{Int}
    end
    # TODO: somehow this throws ArgumentError defined as fallback
    # errmsg2 = "required field displacement not found in TestVerletStorageMissingField1!"
    # e2 = ErrorException(errmsg2)
    @test_throws ArgumentError Peridynamics.@storage Mat1 StorageMissing2

    struct StorageMissing3 <: Peridynamics.AbstractStorage
        position::Matrix{Float64}
        displacement::Matrix{Float64}
        velocity::Matrix{Float64}
        velocity_half::Matrix{Float64}
        velocity_half_old::Matrix{Float64}
        acceleration::Matrix{Float64}
        b_int::Matrix{Float64}
        b_int_old::Matrix{Float64}
        b_ext::Matrix{Float64}
        density_matrix::Matrix{Float64}
        damage::Vector{Float64}
        one_ni_active::Vector{Bool}
    end
    # TODO: somehow this throws ArgumentError defined as fallback
    # errmsg2 = "required field displacement not found in TestVerletStorageMissingField1!"
    # e2 = ErrorException(errmsg2)
    @test_throws ArgumentError Peridynamics.@storage Mat2 StorageMissing3

    struct Storage1 <: Peridynamics.AbstractStorage
        position::Matrix{Float64}
        displacement::Matrix{Float64}
        velocity::Matrix{Float64}
        velocity_half::Matrix{Float64}
        velocity_half_old::Matrix{Float64}
        acceleration::Matrix{Float64}
        b_int::Matrix{Float64}
        b_int_old::Matrix{Float64}
        b_ext::Matrix{Float64}
        density_matrix::Matrix{Float64}
        damage::Vector{Float64}
        bond_active::Vector{Bool}
        n_active_bonds::Vector{Int}
    end

    Peridynamics.@storage Mat1 Storage1
    @test hasmethod(Peridynamics.storage_type, Tuple{Mat1})
    mat, vv = Mat1(), VelocityVerlet(steps=1)
    @test Peridynamics.storage_type(mat) == Storage1

    @test_throws ArgumentError Peridynamics.@loc_to_halo_fields(StorageNoSubtype, :position)
    @test_throws ArgumentError Peridynamics.@loc_to_halo_fields(Storage1, :randomfield)

    Peridynamics.@loc_to_halo_fields Storage1 :position :displacement
    @test hasmethod(Peridynamics.loc_to_halo_fields, Tuple{Storage1})
    @test hasmethod(Peridynamics.is_halo_field, Tuple{Storage1,Val{:position}})
    @test hasmethod(Peridynamics.is_halo_field,
                    Tuple{Storage1,Val{:displacement}})

    @test_throws ArgumentError Peridynamics.@halo_to_loc_fields(StorageNoSubtype,
                                                               :b_int)

    @test_throws ArgumentError Peridynamics.@halo_to_loc_fields(Storage1,
                                                                :randomfield)

    Peridynamics.@halo_to_loc_fields Storage1 :b_int :b_ext
    @test hasmethod(Peridynamics.halo_to_loc_fields, Tuple{Storage1})
    @test hasmethod(Peridynamics.is_halo_field, Tuple{Storage1,Val{:b_int}})
    @test hasmethod(Peridynamics.is_halo_field, Tuple{Storage1,Val{:b_ext}})
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
    @test b1.storage.damage ≈ [2/3, 2/3]
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
    @test b2.storage.damage ≈ [2/3, 2/3]
    @test b2.storage.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.storage.n_active_bonds == [1, 1]

    randpos = rand(3, 4)
    dh.chunks[2].storage.position .= randpos
    Peridynamics.exchange_loc_to_halo!(dh, 1)

    @test b1.storage.position[:,1:2] ≈ [0.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.position[:,3:4] ≈ randpos[:,1:2]
    @test b1.storage.displacement == zeros(3, 2)
    @test b1.storage.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.velocity_half == zeros(3, 2)
    @test b1.storage.acceleration == zeros(3, 2)
    @test b1.storage.b_int == zeros(3, 4)
    @test b1.storage.b_ext == zeros(3, 2)
    @test b1.storage.damage ≈ [2/3, 2/3]
    @test b1.storage.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.storage.n_active_bonds == [1, 1]

    @test b2.storage.position ≈ randpos
    @test b2.storage.displacement == zeros(3, 2)
    @test b2.storage.velocity == [0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test b2.storage.velocity_half == zeros(3, 2)
    @test b2.storage.acceleration == zeros(3, 2)
    @test b2.storage.b_int == zeros(3, 4)
    @test b2.storage.b_ext == zeros(3, 2)
    @test b2.storage.damage ≈ [2/3, 2/3]
    @test b2.storage.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.storage.n_active_bonds == [1, 1]

    randbint = rand(3, 4)
    b2.storage.b_int .= randbint
    b1.storage.b_int .+= 1
    Peridynamics.exchange_halo_to_loc!(dh, 1)

    @test b1.storage.position[:,1:2] ≈ [0.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.position[:,3:4] ≈ randpos[:,1:2]
    @test b1.storage.displacement == zeros(3, 2)
    @test b1.storage.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.velocity_half == zeros(3, 2)
    @test b1.storage.acceleration == zeros(3, 2)
    @test b1.storage.b_int[:,1:2] ≈ 1 .+ randbint[:,3:4]
    @test b1.storage.b_int[:,3:4] ≈ ones(3, 2)
    @test b1.storage.b_ext == zeros(3, 2)
    @test b1.storage.damage ≈ [2/3, 2/3]
    @test b1.storage.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.storage.n_active_bonds == [1, 1]

    @test b2.storage.position ≈ randpos
    @test b2.storage.displacement == zeros(3, 2)
    @test b2.storage.velocity == [0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test b2.storage.velocity_half == zeros(3, 2)
    @test b2.storage.acceleration == zeros(3, 2)
    @test b2.storage.b_int ≈ randbint
    @test b2.storage.b_ext == zeros(3, 2)
    @test b2.storage.damage ≈ [2/3, 2/3]
    @test b2.storage.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.storage.n_active_bonds == [1, 1]

    Peridynamics.exchange_halo_to_loc!(dh, 2)

    @test b1.storage.position[:,1:2] ≈ [0.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.position[:,3:4] ≈ randpos[:,1:2]
    @test b1.storage.displacement == zeros(3, 2)
    @test b1.storage.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.velocity_half == zeros(3, 2)
    @test b1.storage.acceleration == zeros(3, 2)
    @test b1.storage.b_int[:,1:2] ≈ 1 .+ randbint[:,3:4]
    @test b1.storage.b_int[:,3:4] ≈ ones(3, 2)
    @test b1.storage.b_ext == zeros(3, 2)
    @test b1.storage.damage ≈ [2/3, 2/3]
    @test b1.storage.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.storage.n_active_bonds == [1, 1]

    @test b2.storage.position ≈ randpos
    @test b2.storage.displacement == zeros(3, 2)
    @test b2.storage.velocity == [0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test b2.storage.velocity_half == zeros(3, 2)
    @test b2.storage.acceleration == zeros(3, 2)
    @test b2.storage.b_int[:,1:2] ≈ 1 .+ randbint[:,1:2]
    @test b2.storage.b_int[:,3:4] ≈ randbint[:,3:4]
    @test b2.storage.b_ext == zeros(3, 2)
    @test b2.storage.damage ≈ [2/3, 2/3]
    @test b2.storage.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.storage.n_active_bonds == [1, 1]
end

@testitem "TestMaterial full simulation" begin
    include("test_material.jl")
    root = joinpath(@__DIR__, "temp_testmaterial_full_simulation")
    path_vtk = joinpath(root, "vtk")

    N = 10
    l, Δx, δ, a = 1.0, 1/N, 3.015/N, 0.5
    pos, vol = uniform_box(l, l, 0.1l, Δx)
    ids = sortperm(pos[2,:])
    body = Body(TestMaterial(), pos[:, ids], vol[ids])
    material!(body; horizon=3.015Δx, E=2.1e5, nu=0.25, rho=8e-6, Gc=2.7)
    point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, body, :set_a)
    point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, body, :set_b)
    precrack!(body, :set_a, :set_b)
    point_set!(p -> p[2] > l/2-Δx, body, :set_top)
    point_set!(p -> p[2] < -l/2+Δx, body, :set_bottom)
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
