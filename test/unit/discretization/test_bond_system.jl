@testitem "BondSystem" begin
    # setup
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    pd = Peridynamics.PointDecomposition(body, 2)

    # 1
    system = Peridynamics.BondSystem(body, pd, 1)

    @test system.position == position
    @test system.volume == volume
    @test system.neighbor == [2, 3, 4, 1, 3, 4]
    @test system.bond_length == [1.0, 1.0, 1.0, 1.0, √2, √2]
    @test system.fail_permit == [true, true, true, true, true, true]
    @test system.n_neighbors == [3, 3]
    @test system.bond_ids == [1:3, 4:6]

    ch = system.chunk_handler
    @test ch.point_ids == [1, 2, 3, 4]
    @test ch.loc_points == [1, 2]
    @test ch.halo_points == [3, 4]
    @test ch.hidxs_by_src[2] == 3:4

    for i in 1:4
        @test ch.localizer[i] == i
    end

    # 2
    system = Peridynamics.BondSystem(body, pd, 2)

    @test system.position == position[:, [3, 4, 1, 2]]
    @test system.volume == volume[[3, 4, 1, 2]]
    @test system.neighbor == [3, 4, 2, 3, 4, 1]
    @test system.bond_length == [1.0, √2, √2, 1.0, √2, √2]
    @test system.fail_permit == [true, true, true, true, true, true]
    @test system.n_neighbors == [3, 3]
    @test system.bond_ids == [1:3, 4:6]

    ch = system.chunk_handler
    @test ch.point_ids == [3, 4, 1, 2]
    @test ch.loc_points == [3, 4]
    @test ch.halo_points == [1, 2]
    @test ch.hidxs_by_src[1] == 3:4
    @test ch.localizer[3] == 1
    @test ch.localizer[4] == 2
    @test ch.localizer[1] == 3
    @test ch.localizer[2] == 4

    # 3
    mat_incompatible = CKIMaterial()
    body_incompatible = Body(mat_incompatible, position, volume)

    @test_throws ArgumentError Peridynamics.BondSystem(body_incompatible, pd, 1)
end

@testitem "find_bonds!" begin
    using Peridynamics: PointNeighbors
    # setup
    position = [0.0 1.0
                0.0 0.0
                0.0 0.0]
    fail_permit = [true, true]
    δmax = 1.5
    nhs = PointNeighbors.GridNeighborhoodSearch{3}(search_radius=δmax, n_points=2)
    PointNeighbors.initialize!(nhs, position, position)

    # find point 2
    δ = 1.5
    neighbor = Vector{Int}()
    bond_length = Vector{Float64}()
    bp = Vector{Bool}()
    n_neighbors = Peridynamics.find_bonds!(neighbor, bond_length, bp, nhs, position,
                                           fail_permit, δ, 1)
    @test n_neighbors == 1
    @test neighbor == [2]
    @test bond_length == [1.0]
    @test bp == [true]

    # horizon too small - find nothing
    δ = 0.9
    neighbor = Vector{Int}()
    bond_length = Vector{Float64}()
    bp = Vector{Bool}()
    n_neighbors = Peridynamics.find_bonds!(neighbor, bond_length, bp, nhs, position,
                                           fail_permit, δ, 1)
    @test n_neighbors == 0
    @test neighbor == Int[]
    @test bond_length == Float64[]
    @test bp == Bool[]

    # no failure allowed for point 2
    fail_permit[2] = false
    δ = 1.5
    neighbor = Vector{Int}()
    bond_length = Vector{Float64}()
    bp = Vector{Bool}()
    n_neighbors = Peridynamics.find_bonds!(neighbor, bond_length, bp, nhs, position,
                                           fail_permit, δ, 1)
    @test n_neighbors == 1
    @test neighbor == [2]
    @test bond_length == [1.0]
    @test bp == [false]
end

@testitem "find_bonds" begin
    # setup
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1, 1, 1, 1]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)

    # all points are local points
    loc_points = 1:4
    neighbor, bond_length, fail_permit, n_neighbors = Peridynamics.find_bonds(body, loc_points)
    @test neighbor == [2, 3, 4, 1, 3, 4, 1, 2, 4, 1, 2, 3]
    @test bond_length == [1.0, 1.0, 1.0, 1.0, √2, √2, 1.0, √2, √2, 1.0, √2, √2]
    @test fail_permit == fill(true, 12)
    @test n_neighbors == [3, 3, 3, 3]

    loc_points = 1:2
    neighbor, bond_length, fail_permit, n_neighbors = Peridynamics.find_bonds(body, loc_points)
    @test neighbor == [2, 3, 4, 1, 3, 4]
    @test bond_length == [1.0, 1.0, 1.0, 1.0, √2, √2]
    @test fail_permit == fill(true, 6)
    @test n_neighbors == [3, 3]

    loc_points = 2:3
    neighbor, bond_length, fail_permit, n_neighbors = Peridynamics.find_bonds(body, loc_points)
    @test neighbor == [1, 3, 4, 1, 2, 4]
    @test bond_length == [1.0, √2, √2, 1.0, √2, √2]
    @test fail_permit == fill(true, 6)
    @test n_neighbors == [3, 3]
end

@testitem "find_halo_points" begin
    neighbor = [2, 3, 4, 1, 2, 3, 4, 1]

    # no halo point
    loc_points = 1:4
    halo_points = Peridynamics.find_halo_points(neighbor, loc_points)
    @test halo_points == Int[]

    # only 1 halo point
    loc_points = 1:3
    halo_points = Peridynamics.find_halo_points(neighbor, loc_points)
    @test halo_points == [4]

    # 2 halo points
    loc_points = 1:2
    halo_points = Peridynamics.find_halo_points(neighbor, loc_points)
    @test halo_points == [3, 4] || halo_points == [4, 3]

    # 2 halo points
    loc_points = 3:4
    halo_points = Peridynamics.find_halo_points(neighbor, loc_points)
    @test halo_points == [1, 2] || halo_points == [2, 1]
end

@testitem "find_bond_ids" begin
    n_neighbors = [1, 2]
    bond_ids = Peridynamics.find_bond_ids(n_neighbors)
    @test bond_ids == [1:1, 2:3]

    n_neighbors = [3, 4, 5]
    bond_ids = Peridynamics.find_bond_ids(n_neighbors)
    @test bond_ids == [1:3, 4:7, 8:12]
end

@testitem "log material properties" begin
    indentation = 0

    mat = BBMaterial()
    msg = Peridynamics.log_material_property(Val(:randomthing), mat; indentation)
    @test msg == ""
    msg = Peridynamics.log_material_property(Val(:dmgmodel), mat; indentation)
    @test contains(msg, "CriticalStretch")

    mat = OSBMaterial()
    msg = Peridynamics.log_material_property(Val(:randomthing), mat; indentation)
    @test msg == ""
    msg = Peridynamics.log_material_property(Val(:dmgmodel), mat; indentation)
    @test contains(msg, "CriticalStretch")
    msg = Peridynamics.log_material_property(Val(:kernel), mat; indentation)
    @test contains(msg, "linear_kernel")

    mat = CMaterial()
    msg = Peridynamics.log_material_property(Val(:randomthing), mat; indentation)
    @test msg == ""
    msg = Peridynamics.log_material_property(Val(:dmgmodel), mat; indentation)
    @test contains(msg, "CriticalStretch")
    msg = Peridynamics.log_material_property(Val(:kernel), mat; indentation)
    @test contains(msg, "linear_kernel")
    msg = Peridynamics.log_material_property(Val(:zem), mat; indentation)
    @test contains(msg, "ZEMSilling")
end

@testitem "find_bonds: duplicate points are an error" begin
    pos, vol = uniform_box(1, 1, 1, 1 / 5)
    pos[:, end] .= pos[:, end - 1] # two points at the same position
    body = Body(BBMaterial(), pos, vol)
    material!(body, horizon=4 / 100, E=1, rho=1, Gc=1)
    @test_throws ErrorException Peridynamics.find_bonds(body, 1:n_points(body))
end

# `current_bond_length` and `bond_stretch` are the one idiom for the kinematics of a bond. A
# material that inherits `BondLengthCache` reads the cache that `update_bond_lengths!` fills,
# one without the field gets the distance computed, and the two have to agree with the hand
# computed value on both.
@testitem "current_bond_length and bond_stretch: with and without the cache" setup=[Fixtures] begin
    using Peridynamics: current_bond_length, bond_stretch, update_bond_lengths!,
                        each_bond_idx, each_point_idx

    # `BBStorage` carries `bond_length`, `CStorage` does not, so the two chunks take the two
    # branches of the `hasfield` test
    cached = Fixtures.chunk(Fixtures.cube(BBMaterial(); n=4))
    uncached = Fixtures.chunk(Fixtures.cube(CMaterial(); n=4))
    @test hasfield(typeof(cached.storage), :bond_length)
    @test !hasfield(typeof(uncached.storage), :bond_length)

    for chunk in (cached, uncached)
        (; system, storage) = chunk
        # a deformation that stretches every bond by a different amount
        storage.position .= system.position .* [1.02, 0.99, 1.0]
        for i in each_point_idx(system)
            update_bond_lengths!(storage, system, i)
            for bond_id in each_bond_idx(system, i)
                j = Peridynamics.get_neighbor(system, bond_id)
                L = Peridynamics.reference_bond_length(system, bond_id)
                Δx = storage.position[:, j] .- storage.position[:, i]
                l = sqrt(Δx[1]^2 + Δx[2]^2 + Δx[3]^2)
                @test current_bond_length(storage, system, i, bond_id) ≈ l
                @test bond_stretch(storage, system, i, bond_id) ≈ (l - L) / L
            end
        end
    end

    # the cache really is what the accessor reads on a material that keeps one, and the
    # undeformed body has every bond at its reference length
    (; system, storage) = cached
    storage.position .= system.position
    for i in each_point_idx(system)
        update_bond_lengths!(storage, system, i)
    end
    for bond_id in eachindex(system.neighbor)
        @test storage.bond_length[bond_id] ≈ Peridynamics.reference_bond_length(system, bond_id)
    end
    i = 1
    for bond_id in each_bond_idx(system, i)
        @test current_bond_length(storage, system, i, bond_id) === storage.bond_length[bond_id]
        @test bond_stretch(storage, system, i, bond_id) ≈ 0 atol=1e-14
    end
end
