# The material interface: what a material has to implement to work with the package, checked
# on a minimal custom material (`TestMaterialImpl`, which implements only the documented
# interface) and on every built-in material.

@testitem "material interface: every built-in material implements it" setup=[Fixtures] begin
    # The functions the data handlers, time solvers and the export call on a material. If a
    # material misses one, the `InterfaceError` fallback is hit.
    for (name, mat) in Fixtures.MATERIALS
        @testset "$name" begin
            M = typeof(mat)
            @test Peridynamics.system_type(mat) <: Peridynamics.AbstractSystem
            @test Peridynamics.point_param_type(mat) <: Peridynamics.AbstractPointParameters
            @test Peridynamics.storage_type(mat) <: Peridynamics.AbstractStorage
            kwargs = Peridynamics.allowed_material_kwargs(mat)
            @test :horizon in kwargs && :rho in kwargs && :E in kwargs
            # the type-level part of the contract is the fracture bookkeeping of the system,
            # the rest is checked against the solver that is actually used
            required = Peridynamics.required_fields(M)
            @test :damage in required
            S = Peridynamics.storage_type(mat)
            @test Peridynamics.point_data_fields(S) ⊆ fieldnames(S)
            @test all(hasfield(S, f) for f in required)
            @test all(hasfield(S, f) for f in (:position, :displacement, :b_int))
            # a storage for every solver the material supports, with consistent halo queries:
            # the fields exchanged between chunks are declared halo fields and are point data
            body = Fixtures.line10(mat)
            for solver in (VelocityVerlet(steps=1), DynamicRelaxation(steps=1), NewtonKrylov(steps=1))
                Fixtures.supports(solver, mat) || continue
                @test isnothing(Peridynamics.check_storage_contract(mat, solver))
                (; storage) = Fixtures.chunk(body, solver)
                @test storage isa S
                lth = Peridynamics.loc_to_halo_fields(storage)
                htl = Peridynamics.halo_to_loc_fields(storage)
                @test :position in lth
                @test all(Peridynamics.is_halo_field(storage, Val(f)) for f in lth)
                @test all(Peridynamics.is_halo_field(storage, Val(f)) for f in htl)
                @test all(f in Peridynamics.point_data_fields(S) for f in lth)
                @test all(f in Peridynamics.point_data_fields(S) for f in htl)
            end
        end
    end
end

@testitem "material interface: a custom material in chunks with halo exchange" setup=[Fixtures, TestMaterialImpl] begin
    # Two chunks of the 4-point body, each with two local points and the two others as halo.
    # The custom material declares `position` as local-to-halo and `b_int` as halo-to-local
    # field, so exactly those two are exchanged.
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = TestMaterialImpl.TestMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, nu=0.25, Gc=1)
    point_set!(body, :a, 1:2)
    point_set!(body, :b, 3:4)
    velocity_ic!(body, :a, :x, 1.0)
    velocity_bc!(Fixtures.f_t, body, :a, :x)
    forcedensity_bc!(Fixtures.f_t, body, :a, :x)
    precrack!(body, :a, :b)
    ts = VelocityVerlet(steps=10)
    dh = Peridynamics.threads_data_handler(body, ts, 2)
    b1, b2 = dh.chunks

    # the initial state of both chunks
    @test b1.storage isa TestMaterialImpl.TestStorage
    @test b1.storage.position == position
    @test b1.storage.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.b_int == zeros(3, 4) # halo-to-local field: sized for all points
    @test b1.storage.b_ext == zeros(3, 2)
    @test b1.storage.damage ≈ [2 / 3, 2 / 3] # the precrack broke the bonds to the other set
    @test b1.storage.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.storage.n_active_bonds == [1, 1]
    @test b2.storage.position == position[:, [3, 4, 1, 2]]
    @test iszero(b2.storage.velocity)
    @test b2.storage.b_int == zeros(3, 4)
    @test b2.storage.damage ≈ [2 / 3, 2 / 3]
    @test b2.storage.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.storage.n_active_bonds == [1, 1]

    # local-to-halo: the positions of chunk 2 appear in the halo of chunk 1
    newpos = [0.1 0.2 0.3 0.4; 0.5 0.6 0.7 0.8; 0.9 1.0 1.1 1.2]
    b2.storage.position .= newpos
    Peridynamics.exchange_loc_to_halo!(dh, 1)
    @test b1.storage.position[:, 1:2] ≈ position[:, 1:2] # local points untouched
    @test b1.storage.position[:, 3:4] ≈ newpos[:, 1:2]
    @test b2.storage.position ≈ newpos
    @test iszero(b1.storage.b_int) && iszero(b2.storage.b_int) # no other field is touched

    # halo-to-local: the halo forces of chunk 2 are added to the local forces of chunk 1
    newbint = [0.11 0.21 0.31 0.41; 0.12 0.22 0.32 0.42; 0.13 0.23 0.33 0.43]
    b2.storage.b_int .= newbint
    b1.storage.b_int .+= 1
    Peridynamics.exchange_halo_to_loc!(dh, 1)
    @test b1.storage.b_int[:, 1:2] ≈ 1 .+ newbint[:, 3:4]
    @test b1.storage.b_int[:, 3:4] ≈ ones(3, 2)
    @test b2.storage.b_int ≈ newbint # the source is not modified
    @test b1.storage.position[:, 3:4] ≈ newpos[:, 1:2] # no other field is touched
    # and the other way round
    Peridynamics.exchange_halo_to_loc!(dh, 2)
    @test b2.storage.b_int[:, 1:2] ≈ 1 .+ newbint[:, 1:2]
    @test b2.storage.b_int[:, 3:4] ≈ newbint[:, 3:4]
    @test b1.storage.b_int[:, 1:2] ≈ 1 .+ newbint[:, 3:4]
end

@testitem "material interface: a custom material runs a simulation" tags=[:simulation] setup=[Fixtures, TestMaterialImpl] begin
    # The custom material through the public API: a small pre-cracked plate, two explicit steps,
    # with the results exported and read back.
    root = mktempdir()
    N = 10
    l, Δx, δ, a = 1.0, 1 / N, 3.015 / N, 0.5
    pos, vol = uniform_box(l, l, 0.1l, Δx)
    body = Body(TestMaterialImpl.TestMaterial(), pos, vol)
    material!(body; horizon=3.015Δx, E=2.1e5, nu=0.25, rho=8e-6, Gc=2.7)
    point_set!(p -> p[1] ≤ -l / 2 + a && 0 ≤ p[2] ≤ 2δ, body, :set_a)
    point_set!(p -> p[1] ≤ -l / 2 + a && -2δ ≤ p[2] < 0, body, :set_b)
    precrack!(body, :set_a, :set_b)
    point_set!(p -> p[2] > l / 2 - Δx, body, :set_top)
    point_set!(p -> p[2] < -l / 2 + Δx, body, :set_bottom)
    velocity_bc!(t -> -30, body, :set_bottom, :y)
    velocity_bc!(t -> 30, body, :set_top, :y)
    job = Job(body, VelocityVerlet(steps=2); path=root, freq=1)
    submit(job; quiet=true)

    vtk_files = Peridynamics.find_vtk_files(joinpath(root, "vtk"))
    @test length(vtk_files) == 3
    r = read_vtk(last(vtk_files))
    @test size(r[:displacement]) == (3, n_points(body))
    @test !iszero(r[:displacement])
    @test any(r[:damage] .> 0) # the pre-crack
    @test isfile(joinpath(root, "logfile.log"))
end
