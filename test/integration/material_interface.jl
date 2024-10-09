@testitem "Material declaration" begin
    struct TestMaterial1 <: Peridynamics.AbstractMaterial end
    struct WrongTestMaterial end

    @test isnothing(Peridynamics.typecheck_material(TestMaterial1))
    @test_throws ArgumentError Peridynamics.typecheck_material(WrongTestMaterial)
end

@testitem "Point parameters declaration" begin
    struct TestMaterial2 <: Peridynamics.AbstractMaterial end
    Peridynamics.@params TestMaterial2 struct TestPointParameters2 <: Peridynamics.AbstractPointParameters
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
    @test isnothing(Peridynamics.typecheck_params(TestPointParameters2))

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
    @test_throws ArgumentError Peridynamics.typecheck_params(PointParametersNoSubtype)

    struct PointParametersMissingHorizon <: Peridynamics.AbstractPointParameters
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
    @test_throws ErrorException Peridynamics.typecheck_params(PointParametersMissingHorizon)

    @test hasmethod(Peridynamics.point_param_type, Tuple{TestMaterial2})
    @test Peridynamics.point_param_type(TestMaterial2()) == TestPointParameters2
    @test hasmethod(Peridynamics.get_point_params, Tuple{TestMaterial2,Dict{Symbol,Any}})
end

# TODO
# @testitem "System declaration" begin
#     struct TestMaterial3 <: Peridynamics.AbstractMaterial end
#     struct WrongMaterial2 end

#     @test_throws ArgumentError Peridynamics.@system WrongMaterial2 Peridynamics.BondSystem
#     @test_throws ArgumentError Peridynamics.@system TestMaterial3 Peridynamics.BBMaterial

#     Peridynamics.@system TestMaterial3 Peridynamics.BondSystem
#     @test hasmethod(Peridynamics.system_type, Tuple{TestMaterial3})
#     @test Peridynamics.system_type(TestMaterial3()) == Peridynamics.BondSystem
#     @test hasmethod(Peridynamics.get_system,
#                     Tuple{Peridynamics.AbstractBody{TestMaterial3}})
# end

@testitem "Storage declaration" begin
    struct TestMaterial4 <: Peridynamics.AbstractMaterial end

    struct TestVerletStorageNoSubtype
        position::Matrix{Float64}
        displacement::Matrix{Float64}
        velocity::Matrix{Float64}
        velocity_half::Matrix{Float64}
        acceleration::Matrix{Float64}
        b_int::Matrix{Float64}
        b_ext::Matrix{Float64}
        damage::Vector{Float64}
        bond_active::Vector{Bool}
        n_active_bonds::Vector{Int}
    end
    @test_throws ArgumentError Peridynamics.@storage(TestMaterial4, VelocityVerlet,
                                                     TestVerletStorageNoSubtype)

    struct TestVerletStorageMissingField1 <: Peridynamics.AbstractStorage
        position::Matrix{Float64}
        # displacement::Matrix{Float64}
        velocity::Matrix{Float64}
        velocity_half::Matrix{Float64}
        acceleration::Matrix{Float64}
        b_int::Matrix{Float64}
        b_ext::Matrix{Float64}
        damage::Vector{Float64}
        bond_active::Vector{Bool}
        n_active_bonds::Vector{Int}
    end
    @test_throws ErrorException Peridynamics.@storage(TestMaterial4, VelocityVerlet,
                                                      TestVerletStorageMissingField1)

    # TODO
    # struct TestVerletStorageMissingField2 <: Peridynamics.AbstractStorage
    #     position::Matrix{Float64}
    #     displacement::Matrix{Float64}
    #     velocity::Matrix{Float64}
    #     velocity_half::Matrix{Float64}
    #     acceleration::Matrix{Float64}
    #     b_int::Matrix{Float64}
    #     b_ext::Matrix{Float64}
    #     # damage::Vector{Float64}
    #     bond_active::Vector{Bool}
    #     n_active_bonds::Vector{Int}
    # end
    # @test_throws ErrorException Peridynamics.@storage(TestMaterial4, VelocityVerlet,
    #                                                  TestVerletStorageMissingField2)

    struct TestVerletStorage1 <: Peridynamics.AbstractStorage
        position::Matrix{Float64}
        displacement::Matrix{Float64}
        velocity::Matrix{Float64}
        velocity_half::Matrix{Float64}
        acceleration::Matrix{Float64}
        b_int::Matrix{Float64}
        b_ext::Matrix{Float64}
        damage::Vector{Float64}
        bond_active::Vector{Bool}
        n_active_bonds::Vector{Int}
    end

    Peridynamics.@storage TestMaterial4 VelocityVerlet TestVerletStorage1
    @test hasmethod(Peridynamics.storage_type, Tuple{TestMaterial4,VelocityVerlet})
    mat, vv = TestMaterial4(), VelocityVerlet(steps=1)
    @test Peridynamics.storage_type(mat, vv) == TestVerletStorage1

    @test_throws ArgumentError Peridynamics.@loc_to_halo_fields(TestVerletStorageNoSubtype,
                                                              :position)

    @test_throws ArgumentError Peridynamics.@loc_to_halo_fields(TestVerletStorage1,
                                                              :randomfield)

    Peridynamics.@loc_to_halo_fields TestVerletStorage1 :position :displacement
    @test hasmethod(Peridynamics.loc_to_halo_fields, Tuple{TestVerletStorage1})
    @test hasmethod(Peridynamics.is_halo_field, Tuple{TestVerletStorage1,Val{:position}})
    @test hasmethod(Peridynamics.is_halo_field,
                    Tuple{TestVerletStorage1,Val{:displacement}})

    @test_throws ArgumentError Peridynamics.@halo_to_loc_fields(TestVerletStorageNoSubtype,
                                                               :b_int)

    @test_throws ArgumentError Peridynamics.@halo_to_loc_fields(TestVerletStorage1,
                                                               :randomfield)

    Peridynamics.@halo_to_loc_fields TestVerletStorage1 :b_int :b_ext
    @test hasmethod(Peridynamics.halo_to_loc_fields, Tuple{TestVerletStorage1})
    @test hasmethod(Peridynamics.is_halo_field, Tuple{TestVerletStorage1,Val{:b_int}})
    @test hasmethod(Peridynamics.is_halo_field, Tuple{TestVerletStorage1,Val{:b_ext}})
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
