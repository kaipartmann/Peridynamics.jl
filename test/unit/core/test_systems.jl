@testitem "DOF handling, raw functions" begin
    A = zeros(Int, 3, 10)
    for (dof, dim, i) in Peridynamics.each_dof_idx(size(A, 1), axes(A, 2))
        A[dof] = dof
    end
    @test A[1, 1] == 1
    @test A[2, 1] == 2
    @test A[3, 1] == 3
    @test A[1, 2] == 4
    @test A[2, 2] == 5
    @test A[3, 2] == 6
    @test A[1, 10] == 28
    @test A[2, 10] == 29
    @test A[3, 10] == 30

    @test Peridynamics.get_dof(3, 1, 1) == 1
    @test Peridynamics.get_dof(3, 2, 1) == 2
    @test Peridynamics.get_dof(3, 3, 1) == 3
    @test Peridynamics.get_dof(3, 1, 2) == 4
    @test Peridynamics.get_dof(3, 2, 2) == 5
    @test Peridynamics.get_dof(3, 3, 2) == 6
    @test Peridynamics.get_dof(3, 1, 10) == 28
    @test Peridynamics.get_dof(3, 2, 10) == 29
    @test Peridynamics.get_dof(3, 3, 10) == 30

    @test Peridynamics.get_point(3, 1) == 1
    @test Peridynamics.get_point(3, 2) == 1
    @test Peridynamics.get_point(3, 3) == 1
    @test Peridynamics.get_point(3, 4) == 2
    @test Peridynamics.get_point(3, 5) == 2
    @test Peridynamics.get_point(3, 6) == 2
    @test Peridynamics.get_point(3, 28) == 10
    @test Peridynamics.get_point(3, 29) == 10
    @test Peridynamics.get_point(3, 30) == 10

    @test Peridynamics.get_dim(3, 1) == 1
    @test Peridynamics.get_dim(3, 2) == 2
    @test Peridynamics.get_dim(3, 3) == 3
    @test Peridynamics.get_dim(3, 4) == 1
    @test Peridynamics.get_dim(3, 5) == 2
    @test Peridynamics.get_dim(3, 6) == 3
    @test Peridynamics.get_dim(3, 28) == 1
    @test Peridynamics.get_dim(3, 29) == 2
    @test Peridynamics.get_dim(3, 30) == 3
end

@testitem "DOF handling, BondSystem interface, -t 1" begin
    position, volume = uniform_box(1,1,1,0.25)
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    ts = VelocityVerlet(steps=10)
    pd = Peridynamics.PointDecomposition(body, 1)
    system = Peridynamics.get_system(body, pd, 1)

    @test Peridynamics.get_n_dim(system) == 3
    @test Peridynamics.get_n_points(system) == 64
    @test Peridynamics.get_n_loc_points(system) == 64
    @test Peridynamics.get_n_dof(system) == 192
    @test Peridynamics.get_n_loc_dof(system) == 192
    @test Peridynamics.get_dof(system, 1, 1) == 1
    @test Peridynamics.get_dof(system, 2, 1) == 2
    @test Peridynamics.get_dof(system, 3, 1) == 3
    @test Peridynamics.get_dof(system, 1, 2) == 4
    @test Peridynamics.get_dof(system, 2, 2) == 5
    @test Peridynamics.get_dof(system, 3, 2) == 6
    @test Peridynamics.get_dof(system, 1, 10) == 28
    @test Peridynamics.get_dof(system, 2, 10) == 29
    @test Peridynamics.get_dof(system, 3, 10) == 30
    @test Peridynamics.get_dof(system, 1, 64) == 190

    @test Peridynamics.each_dim(system) == 1:3
    all_dof_idxs = collect(Peridynamics.each_dof_idx(system))
    @test size(all_dof_idxs) == (64, 3)
    @test all_dof_idxs[1, 1] == (1, 1, 1)
    @test all_dof_idxs[1, 2] == (2, 2, 1)
    @test all_dof_idxs[1, 3] == (3, 3, 1)
    @test all_dof_idxs[2, 1] == (4, 1, 2)
    @test all_dof_idxs[2, 2] == (5, 2, 2)
    @test all_dof_idxs[2, 3] == (6, 3, 2)
    @test all_dof_idxs[10, 1] == (28, 1, 10)
    @test all_dof_idxs[10, 2] == (29, 2, 10)
    @test all_dof_idxs[10, 3] == (30, 3, 10)
    @test all_dof_idxs[64, 1] == (190, 1, 64)
    @test all_dof_idxs == collect(Peridynamics.each_loc_dof_idx(system))
    @test collect(Peridynamics.each_dof_idx(system, [1,2,10])) == [
        (1, 1, 1) (2, 2, 1) (3, 3, 1)
        (4, 1, 2) (5, 2, 2) (6, 3, 2)
        (28, 1, 10) (29, 2, 10) (30, 3, 10)
    ]
    all_dofs = collect(Peridynamics.each_dof(system))
    @test size(all_dofs) == (64, 3)
    @test all_dofs[1, 1] == 1
    @test all_dofs[1, 2] == 2
    @test all_dofs[1, 3] == 3
    @test all_dofs[2, 1] == 4
    @test all_dofs[2, 2] == 5
    @test all_dofs[2, 3] == 6
    @test all_dofs[10, 1] == 28
    @test all_dofs[10, 2] == 29
    @test all_dofs[10, 3] == 30
    @test all_dofs[64, 1] == 190
    @test all_dofs == collect(Peridynamics.each_loc_dof(system))

    @test Peridynamics.get_point(system, 1) == 1
    @test Peridynamics.get_dim(system, 1) == 1
    @test Peridynamics.get_point(system, 2) == 1
    @test Peridynamics.get_dim(system, 2) == 2
    @test Peridynamics.get_point(system, 3) == 1
    @test Peridynamics.get_dim(system, 3) == 3
    @test Peridynamics.get_point(system, 4) == 2
    @test Peridynamics.get_dim(system, 4) == 1
    @test Peridynamics.get_point(system, 28) == 10
    @test Peridynamics.get_dim(system, 28) == 1
    @test Peridynamics.get_point(system, 29) == 10
    @test Peridynamics.get_dim(system, 29) == 2
    @test Peridynamics.get_point(system, 30) == 10
    @test Peridynamics.get_dim(system, 30) == 3
end

@testitem "DOF handling, BondSystem interface, -t 2" begin
    position, volume = uniform_box(1,1,1,0.5)
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    ts = VelocityVerlet(steps=10)
    pd = Peridynamics.PointDecomposition(body, 2)

    X1 = [-0.25, -0.25, -0.25]
    X4 = [0.25, 0.25, -0.25]
    X5 = [-0.25, -0.25, 0.25]
    X8 = [0.25, 0.25, 0.25]

    # first chunk
    system = Peridynamics.get_system(body, pd, 1)

    # point ids on the first chunk similar to the body point ids
    @test system.position[:, 1] ≈ X1
    @test system.position[:, 4] ≈ X4
    @test system.position[:, 5] ≈ X5
    @test system.position[:, 8] ≈ X8

    @test Peridynamics.get_n_dim(system) == 3
    @test Peridynamics.get_n_points(system) == 8
    @test Peridynamics.get_n_loc_points(system) == 4
    @test Peridynamics.get_n_dof(system) == 24
    @test Peridynamics.get_n_loc_dof(system) == 12
    @test Peridynamics.get_dof(system, 1, 1) == 1
    @test Peridynamics.get_dof(system, 2, 1) == 2
    @test Peridynamics.get_dof(system, 3, 1) == 3
    @test Peridynamics.get_dof(system, 1, 2) == 4
    @test Peridynamics.get_dof(system, 2, 2) == 5
    @test Peridynamics.get_dof(system, 3, 2) == 6
    @test Peridynamics.get_dof(system, 1, 8) == 22
    @test Peridynamics.get_dof(system, 2, 8) == 23
    @test Peridynamics.get_dof(system, 3, 8) == 24

    @test Peridynamics.each_dim(system) == 1:3
    all_dof_idxs = collect(Peridynamics.each_dof_idx(system))
    @test size(all_dof_idxs) == (8, 3)
    @test all_dof_idxs[1, 1] == (1, 1, 1)
    @test all_dof_idxs[1, 2] == (2, 2, 1)
    @test all_dof_idxs[1, 3] == (3, 3, 1)
    @test all_dof_idxs[2, 1] == (4, 1, 2)
    @test all_dof_idxs[2, 2] == (5, 2, 2)
    @test all_dof_idxs[2, 3] == (6, 3, 2)
    @test all_dof_idxs[8, 1] == (22, 1, 8)
    @test all_dof_idxs[8, 2] == (23, 2, 8)
    @test all_dof_idxs[8, 3] == (24, 3, 8)
    @test all_dof_idxs[1:4, :] == collect(Peridynamics.each_loc_dof_idx(system))
    @test collect(Peridynamics.each_dof_idx(system, [1,2,8])) == [
        (1, 1, 1) (2, 2, 1) (3, 3, 1)
        (4, 1, 2) (5, 2, 2) (6, 3, 2)
        (22, 1, 8) (23, 2, 8) (24, 3, 8)
    ]
    all_dofs = collect(Peridynamics.each_dof(system))
    @test size(all_dofs) == (8, 3)
    @test all_dofs[1, 1] == 1
    @test all_dofs[1, 2] == 2
    @test all_dofs[1, 3] == 3
    @test all_dofs[2, 1] == 4
    @test all_dofs[2, 2] == 5
    @test all_dofs[2, 3] == 6
    @test all_dofs[8, 1] == 22
    @test all_dofs[8, 2] == 23
    @test all_dofs[8, 3] == 24
    @test all_dofs[1:4, :] == collect(Peridynamics.each_loc_dof(system))

    @test Peridynamics.get_point(system, 1) == 1
    @test Peridynamics.get_dim(system, 1) == 1
    @test Peridynamics.get_point(system, 2) == 1
    @test Peridynamics.get_dim(system, 2) == 2
    @test Peridynamics.get_point(system, 3) == 1
    @test Peridynamics.get_dim(system, 3) == 3
    @test Peridynamics.get_point(system, 4) == 2
    @test Peridynamics.get_dim(system, 4) == 1
    # no bounds checking, so the following still works although we have only 8 points
    @test Peridynamics.get_point(system, 28) == 10
    @test Peridynamics.get_dim(system, 28) == 1
    @test Peridynamics.get_point(system, 29) == 10
    @test Peridynamics.get_dim(system, 29) == 2
    @test Peridynamics.get_point(system, 30) == 10
    @test Peridynamics.get_dim(system, 30) == 3

    # second chunk
    system = Peridynamics.get_system(body, pd, 2)

    # point ids on the second chunk NOT similar to the body point ids
    @test system.position[:, 1] ≈ X5
    @test system.position[:, 4] ≈ X8
    @test system.position[:, 5] ≈ X1
    @test system.position[:, 8] ≈ X4

    @test Peridynamics.get_n_dim(system) == 3
    @test Peridynamics.get_n_points(system) == 8
    @test Peridynamics.get_n_loc_points(system) == 4
    @test Peridynamics.get_n_dof(system) == 24
    @test Peridynamics.get_n_loc_dof(system) == 12
    @test Peridynamics.get_dof(system, 1, 1) == 1
    @test Peridynamics.get_dof(system, 2, 1) == 2
    @test Peridynamics.get_dof(system, 3, 1) == 3
    @test Peridynamics.get_dof(system, 1, 2) == 4
    @test Peridynamics.get_dof(system, 2, 2) == 5
    @test Peridynamics.get_dof(system, 3, 2) == 6
    @test Peridynamics.get_dof(system, 1, 8) == 22
    @test Peridynamics.get_dof(system, 2, 8) == 23
    @test Peridynamics.get_dof(system, 3, 8) == 24

    @test Peridynamics.each_dim(system) == 1:3
    all_dof_idxs = collect(Peridynamics.each_dof_idx(system))
    @test size(all_dof_idxs) == (8, 3)
    @test all_dof_idxs[1, 1] == (1, 1, 1)
    @test all_dof_idxs[1, 2] == (2, 2, 1)
    @test all_dof_idxs[1, 3] == (3, 3, 1)
    @test all_dof_idxs[2, 1] == (4, 1, 2)
    @test all_dof_idxs[2, 2] == (5, 2, 2)
    @test all_dof_idxs[2, 3] == (6, 3, 2)
    @test all_dof_idxs[8, 1] == (22, 1, 8)
    @test all_dof_idxs[8, 2] == (23, 2, 8)
    @test all_dof_idxs[8, 3] == (24, 3, 8)
    @test all_dof_idxs[1:4, :] == collect(Peridynamics.each_loc_dof_idx(system))
    @test collect(Peridynamics.each_dof_idx(system, [1,2,8])) == [
        (1, 1, 1) (2, 2, 1) (3, 3, 1)
        (4, 1, 2) (5, 2, 2) (6, 3, 2)
        (22, 1, 8) (23, 2, 8) (24, 3, 8)
    ]
    all_dofs = collect(Peridynamics.each_dof(system))
    @test size(all_dofs) == (8, 3)
    @test all_dofs[1, 1] == 1
    @test all_dofs[1, 2] == 2
    @test all_dofs[1, 3] == 3
    @test all_dofs[2, 1] == 4
    @test all_dofs[2, 2] == 5
    @test all_dofs[2, 3] == 6
    @test all_dofs[8, 1] == 22
    @test all_dofs[8, 2] == 23
    @test all_dofs[8, 3] == 24
    @test all_dofs[1:4, :] == collect(Peridynamics.each_loc_dof(system))

    @test Peridynamics.get_point(system, 1) == 1
    @test Peridynamics.get_dim(system, 1) == 1
    @test Peridynamics.get_point(system, 2) == 1
    @test Peridynamics.get_dim(system, 2) == 2
    @test Peridynamics.get_point(system, 3) == 1
    @test Peridynamics.get_dim(system, 3) == 3
    @test Peridynamics.get_point(system, 4) == 2
    @test Peridynamics.get_dim(system, 4) == 1
    # no bounds checking, so the following still works although we have only 8 points
    @test Peridynamics.get_point(system, 28) == 10
    @test Peridynamics.get_dim(system, 28) == 1
    @test Peridynamics.get_point(system, 29) == 10
    @test Peridynamics.get_dim(system, 29) == 2
    @test Peridynamics.get_point(system, 30) == 10
    @test Peridynamics.get_dim(system, 30) == 3
end

@testitem "system interface: fallbacks and forwarding to the chunk handler" setup=[Fixtures] begin
    struct NoSystemMaterial <: Peridynamics.AbstractMaterial end
    @test_throws Peridynamics.InterfaceError Peridynamics.system_type(NoSystemMaterial())

    c = Fixtures.chunk(Fixtures.line10(); n_chunks=2, chunk_id=1)
    system, ch = c.system, c.system.chunk_handler
    @test Peridynamics.get_halo_points(system) == Peridynamics.get_halo_points(ch)
    a = reshape(collect(1.0:3 * Peridynamics.get_n_points(system)), 3, :)
    @test Peridynamics.get_loc_view(a, system) == Peridynamics.get_loc_view(a, ch)
    @test vec(collect(Peridynamics.each_dof(system, [1, 3]))) ==
          vec(collect(Peridynamics.each_dof(3, [1, 3])))
    @test vec(collect(Peridynamics.each_dof(system, [2]))) == [4, 5, 6]
end

@testitem "system_type(mat) matches the constructed system, N = 3" setup=[Fixtures] begin
    # a bond-based, a bond-associated and an interaction-system material: `system_type`
    # is what `body_chunk_type` uses to preallocate before any system is built, so it has
    # to agree with what `get_system` actually returns, `N` included
    for (mat, kwargs) in ((BBMaterial(), (;)),
                          (BACMaterial(), (;)),
                          (CKIMaterial(), Fixtures.cki_kwargs()))
        c = Fixtures.chunk(Fixtures.cube(mat; n=3, kwargs...))
        @test typeof(c.system) === Peridynamics.system_type(mat)
        @test Peridynamics.get_n_dim(c.system) == 3
    end
end

# --- the @system macro, the twin of test_storages.jl ---

@testitem "@system: the generated header and the order of its type parameters" begin
    import Peridynamics: @system, AbstractCorrection, AbstractSystem, get_n_dim, float_type,
                         host_system_type, ChunkHandler, BondSystem

    @system struct MacroSystem{Correction<:AbstractCorrection} <: AbstractSystem
        position::PointVector{Float64}
        volume::PointScalar
        neighbor::BondScalar{Int}
        bond_ids::PointScalar{UnitRange{Int}}
        kernels::BondScalar
        correction::Correction
    end

    # the declared parameters come first, then `N`, then `FT`, then `CH` for the injected
    # chunk handler, then one parameter per distinct array type in the order in which the
    # fields ask for them; `chunk_handler` itself is the last field, injected by the macro
    @test MacroSystem.body.body.body.body.body.body.body.body isa DataType
    @test fieldnames(MacroSystem) ===
          (:position, :volume, :neighbor, :bond_ids, :kernels, :correction, :chunk_handler)

    S = host_system_type(MacroSystem, NoCorrection, Val(3), Float64)
    @test S === MacroSystem{NoCorrection,3,Float64,ChunkHandler,Matrix{Float64},
                            Vector{Float64},Vector{Int},Vector{UnitRange{Int}}}
    @test S.parameters[1] === NoCorrection
    @test S.parameters[2] === 3
    @test S.parameters[3] === Float64
    @test S.parameters[4] === ChunkHandler

    # the float type of the simulation reaches every field declared without an element type,
    # and the pinned `position` keeps `Float64`
    S32 = host_system_type(MacroSystem, NoCorrection, Val(2), Float32)
    @test fieldtype(S32, :position) === Matrix{Float64}
    @test fieldtype(S32, :volume) === Vector{Float32}
    @test fieldtype(S32, :kernels) === Vector{Float32}
    @test fieldtype(S32, :neighbor) === Vector{Int}
    @test fieldtype(S32, :chunk_handler) === ChunkHandler

    # `FT` defaults to the float type of the simulation
    @test host_system_type(MacroSystem, NoCorrection, Val(3)) === S

    # a `<:` pattern that names only the declared parameter still dispatches, which is what
    # every dispatch on a system depends on
    @test S <: MacroSystem{<:AbstractCorrection}
    @test S <: MacroSystem{NoCorrection}

    # the same header order holds for a system of this package, and `chunk_handler` is
    # concrete after `host_system_type`
    names = [p.name for p in Base.unwrap_unionall(BondSystem).parameters]
    @test names == [:Correction, :N, :FT, :CH, names[5:end]...]
    SBond = host_system_type(BondSystem, NoCorrection, Val(3), Float64)
    @test fieldtype(SBond, :chunk_handler) === ChunkHandler
    @test isconcretetype(fieldtype(SBond, :chunk_handler))
end

@testitem "@system: the generated constructor, accessors and Adapt" begin
    import Peridynamics: @system, AbstractSystem, ChunkHandler, SystemSizes, alloc_field,
                         BondScalar, PointScalar, LocalPoints, get_n_dim, float_type,
                         get_n_bonds, get_n_points, get_n_loc_points, host_system_type,
                         storage_fields_expr

    @system struct CtorSystem <: AbstractSystem
        position::PointVector{Float64}
        volume::PointScalar
        neighbor::BondScalar{Int}
    end

    pos, vol = uniform_box(1, 1, 1, 0.5)
    body = Body(BBMaterial(), pos, vol)
    material!(body; horizon=0.8, rho=1, E=1, nu=0.25, Gc=1)
    pd = Peridynamics.PointDecomposition(body, 1)
    ch = Peridynamics.get_system(body, pd, 1).chunk_handler

    system = CtorSystem{3,Float64}(pos, vol, [1, 2, 3, 4], ch)

    # the array parameters and `CH` are inferred from the values, `N` and `FT` are the ones
    # named
    @test system isa CtorSystem{3,Float64,ChunkHandler,Matrix{Float64},Vector{Float64},
                                Vector{Int}}
    @test get_n_dim(system) == 3
    @test get_n_dim(typeof(system)) == 3
    @test float_type(system) === Float64
    @test get_n_bonds(system) == 4
    # the point counts come from the chunk handler through the `AbstractSystem` forwarding
    @test get_n_loc_points(system) == get_n_loc_points(ch)
    @test get_n_points(system) == get_n_points(ch)

    # the declarations are registered, so `block_table` can read them, and the injected
    # `chunk_handler` shows up last
    @test [d.name for d in storage_fields_expr(CtorSystem)] ==
          [:position, :volume, :neighbor, :chunk_handler]

    # `position` and `N` cannot disagree
    @test_throws DimensionMismatch CtorSystem{2,Float64}(pos, vol, [1, 2, 3, 4], ch)

    # moving a system to another array backend moves every field, including the chunk
    # handler, and keeps `N` and `FT`
    struct SysWrappedArray{T,N} <: AbstractArray{T,N}
        a::Array{T,N}
    end
    Base.size(x::SysWrappedArray) = size(x.a)
    Base.getindex(x::SysWrappedArray, i...) = getindex(x.a, i...)
    struct SysWrappedBackend end
    function Peridynamics.Adapt.adapt_storage(::SysWrappedBackend,
                                              a::Array{T,N}) where {T,N}
        return SysWrappedArray{T,N}(a)
    end

    moved = Peridynamics.Adapt.adapt(SysWrappedBackend(), system)
    @test moved isa CtorSystem{3,Float64}
    @test moved.position isa SysWrappedArray{Float64,2}
    @test moved.neighbor isa SysWrappedArray{Int,1}
    @test get_n_dim(moved) == 3
    @test float_type(moved) === Float64
    @test get_n_bonds(moved) == 4

    # nothing moves when the backend is the one the system already lives on
    @test Peridynamics.Adapt.adapt(Array, system) === system
end

@testitem "@system: a system without bonds and a field kept verbatim" begin
    import Peridynamics: @system, AbstractSystem, get_n_bonds, host_system_type,
                         ChunkHandler, SVector

    struct GridOffset
        offset::NTuple{3,Int}
    end

    @system struct GridSystem <: AbstractSystem
        position::PointVector{Float64}
        grid_idx::PointScalar{Int}
        offsets::Vector{GridOffset}
        origin::SVector{3,Float64}
        n_neighbors::Int
    end

    S = host_system_type(GridSystem, Val(3), Float64)
    # a concrete `Array` field takes part in the parameters, an isbits field does not
    @test fieldtype(S, :offsets) === Vector{GridOffset}
    @test fieldtype(S, :origin) === SVector{3,Float64}
    @test fieldtype(S, :n_neighbors) === Int
    @test fieldtype(S, :grid_idx) === Vector{Int}
    @test fieldtype(S, :chunk_handler) === ChunkHandler

    # no field has a bond shape, so the macro generates no `get_n_bonds` for this system
    # and only the generic fallback is left, which fails with the native no-field error
    @test which(get_n_bonds, Tuple{GridSystem}).sig ===
          Tuple{typeof(get_n_bonds),Peridynamics.AbstractSystem}
    @test_throws Exception get_n_bonds(S(zeros(3, 2), [1, 2], GridOffset[],
                                        SVector{3,Float64}(0, 0, 0), 0,
                                        Peridynamics.ChunkHandler(2, [1, 2], 1:2, Int[],
                                                                  Dict{Int,UnitRange{Int}}(),
                                                                  Dict(1 => 1, 2 => 2))))
end

@testitem "@system: what a system may not declare" begin
    import Peridynamics: get_system_header, user_param_name, macrocheck_input_system_struct

    # a system may declare type parameters of its own, unlike a storage
    name, params, super = get_system_header(:(struct MySys{P<:Real} <: MySuper end))
    @test name === :MySys
    @test params == Any[:(P <: Real)]
    @test super === :MySuper
    @test user_param_name(params[1]) === :P
    name, params, super = get_system_header(:(struct MySys end))
    @test name === :MySys
    @test isempty(params)
    @test super == :(Peridynamics.AbstractSystem)
    @test_throws ArgumentError get_system_header(:(struct (a + b) end))
    @test_throws ArgumentError user_param_name(:(P <: Real <: Q))

    @test isnothing(macrocheck_input_system_struct(:(struct S
                                                         a::Int
                                                     end)))
    @test_throws ArgumentError macrocheck_input_system_struct(:(MySystem))

    # the macro provides `chunk_handler` itself, declaring it is an error
    err = try
        @eval Peridynamics.@system struct HandlerSystem <: Peridynamics.AbstractSystem
            position::PointVector{Float64}
            chunk_handler::Peridynamics.AbstractChunkHandler
        end
    catch e
        e
    end
    @test err isa LoadError
    @test contains(err.error.msg, "provides this field itself")

    # reusable field blocks exist only for a storage
    err = try
        @eval Peridynamics.@system struct InheritSystem <: Peridynamics.AbstractSystem
            @inherit Peridynamics.VelocityVerletFields
            position::PointVector{Float64}
        end
    catch e
        e
    end
    @test err isa LoadError
    @test contains(err.error.msg, "blocks exist only for storages")

    # a field declared with an abstract type is an error, unlike before
    err = try
        @eval Peridynamics.@system struct AbstractFieldSystem <: Peridynamics.AbstractSystem
            position::PointVector{Float64}
            correction::Peridynamics.AbstractCorrection
        end
    catch e
        e
    end
    @test err isa LoadError
    @test contains(err.error.msg, "abstract type")

    # `N`, `FT` and `CH` are reserved for the macro
    err = try
        @eval Peridynamics.@system struct ReservedParamSystem{CH} <: Peridynamics.AbstractSystem
            position::PointVector{Float64}
        end
    catch e
        e
    end
    @test err isa LoadError
    @test contains(err.error.msg, "collides with the parameter")

    # the constructor of a system fills every field, so an initial value would never be read
    err = try
        @eval Peridynamics.@system struct InitSystem <: Peridynamics.AbstractSystem
            position::PointVector{Float64}
            volume::PointScalar = 1.0
        end
    catch e
        e
    end
    @test err isa LoadError
    @test contains(err.error.msg, "specifies the initial value")

    # a system is never exchanged between chunks
    err = try
        @eval Peridynamics.@system struct HaloSystem <: Peridynamics.AbstractSystem
            @lth position::PointVector{Float64}
        end
    catch e
        e
    end
    @test err isa LoadError
    @test contains(err.error.msg, "is annotated with")

    # model state belongs into the storage of the material
    err = try
        @eval Peridynamics.@system struct StateSystem <: Peridynamics.AbstractSystem
            position::PointVector{Float64}
            dmg_state::DamageState
        end
    catch e
        e
    end
    @test err isa LoadError
    @test contains(err.error.msg, "carries no model state")
end

@testitem "SystemSizes: what a constructor allocates against" begin
    import Peridynamics: SystemSizes, alloc_field, PointScalar, PointVector, BondScalar,
                         LocalPoints, HaloPoints, get_n_dim, float_type, get_n_bonds,
                         get_n_loc_points, get_n_points

    pos, vol = uniform_box(1, 1, 1, 0.5)
    body = Body(BBMaterial(), pos, vol)
    material!(body; horizon=0.8, rho=1, E=1, nu=0.25, Gc=1)
    pd = Peridynamics.PointDecomposition(body, 2)
    ch = Peridynamics.get_system(body, pd, 1).chunk_handler

    sizes = SystemSizes{2,Float32}(ch, 17)
    @test get_n_dim(sizes) == 2
    @test float_type(sizes) === Float32
    @test get_n_bonds(sizes) == 17
    @test get_n_loc_points(sizes) == get_n_loc_points(ch)
    @test get_n_points(sizes) == get_n_points(ch)

    # every shape of a storage field works on it, which is what replaces `zeros(3, n)`
    @test size(alloc_field(PointVector(), sizes, LocalPoints())) ==
          (2, get_n_loc_points(ch))
    @test size(alloc_field(PointVector(), sizes, HaloPoints())) == (2, get_n_points(ch))
    @test size(alloc_field(BondScalar(), sizes, LocalPoints())) == (17,)
    @test eltype(alloc_field(PointScalar(), sizes, LocalPoints())) === Float32
    @test all(isone, alloc_field(BondScalar(), sizes, LocalPoints(), 1))
end

@testitem "system_type and check_system_compat of the three systems" begin
    import Peridynamics: system_type, check_system_compat, host_system_type, get_n_dim,
                         float_type, BondSystem, BondAssociatedSystem, InteractionSystem,
                         NoCorrection, EnergySurfaceCorrection, ChunkHandler, InterfaceError,
                         AbstractMaterial

    # one line per system family, with the float type and the dimension of the simulation
    @test system_type(BBMaterial()) ===
          host_system_type(BondSystem, NoCorrection, Val(3), Float64)
    @test system_type(BACMaterial()) === host_system_type(BondAssociatedSystem, Val(3),
                                                          Float64)
    @test system_type(CKIMaterial()) === host_system_type(InteractionSystem, Val(3),
                                                          Float64)

    S2 = system_type(BBMaterial(), Float32, Val(2))
    @test get_n_dim(S2) == 2
    @test S2 <: BondSystem{NoCorrection,2,Float32}

    # the correction is a host type as well, so a simulation in `Float32` gets `Float32`
    # correction arrays
    SE = system_type(BBMaterial{EnergySurfaceCorrection}(), Float32, Val(3))
    @test SE <: BondSystem{EnergySurfaceCorrection{Matrix{Float32},Vector{Float32}}}

    # a material without a system says so
    struct NoSystemMat <: AbstractMaterial end
    @test_throws InterfaceError system_type(NoSystemMat())

    # one function instead of the three `check_*_compat` of before
    @test isnothing(check_system_compat(BondSystem, BBMaterial()))
    @test isnothing(check_system_compat(BondAssociatedSystem, BACMaterial()))
    @test isnothing(check_system_compat(InteractionSystem, CKIMaterial()))
    @test_throws ArgumentError check_system_compat(BondSystem, CKIMaterial())
    @test_throws ArgumentError check_system_compat(BondAssociatedSystem, BBMaterial())
    @test_throws ArgumentError check_system_compat(InteractionSystem, BBMaterial())

    # a system without a restriction accepts every material
    @test isnothing(check_system_compat(Peridynamics.AbstractSystem, BBMaterial()))
end
