# The storage field declaration framework of `src/core/storage_fields.jl`: field shapes,
# allocation, solver markers, `@storage_fields` / `@inherit`, and the derived type parameters
# of a generated storage. The storage contract and the `@storage` checks are in
# `test_storages.jl`.

@testitem "field shapes: shape, element type and container of a declaration" begin
    import Peridynamics: PointScalar, PointVector, PointTensor, PointSymTensor, PointField,
                         BondScalar, BondVector, BondTensor, BondSymTensor, BondField,
                         DofVector, SimFloat, field_shape, container_type, concrete_eltype,
                         field_eltype, storage_field_type

    # a plain container type is not a field shape and is used as declared
    @test field_shape(Matrix{Float64}) === nothing
    @test field_shape(Vector{Bool}) === nothing
    @test storage_field_type(Matrix{Float64}) === Matrix{Float64}

    # a shape without an element type follows the float type of the simulation
    @test field_shape(PointScalar) === PointScalar{SimFloat}()
    @test field_shape(PointVector) === PointVector{SimFloat}()
    @test field_shape(BondTensor) === BondTensor{SimFloat}()
    @test field_shape(DofVector) === DofVector{SimFloat}()
    @test field_shape(PointField{7}) === PointField{7,SimFloat}()
    @test field_shape(BondField{7}) === BondField{7,SimFloat}()
    @test field_eltype(PointVector{SimFloat}()) === SimFloat
    @test concrete_eltype(PointVector{SimFloat}()) === Float64

    # a shape with an element type pins it
    @test field_shape(PointScalar{Bool}) === PointScalar{Bool}()
    @test field_shape(PointScalar{Int}) === PointScalar{Int}()
    @test field_shape(PointVector{Float32}) === PointVector{Float32}()
    @test field_shape(PointField{7,Float32}) === PointField{7,Float32}()
    @test concrete_eltype(PointVector{Float32}()) === Float32

    # scalar and dof shapes are vectors, everything else is a matrix
    @test container_type(PointScalar{Float64}()) === Vector{Float64}
    @test container_type(PointScalar{Bool}()) === Vector{Bool}
    @test container_type(BondScalar{Int}()) === Vector{Int}
    @test container_type(DofVector{Float64}()) === Vector{Float64}
    @test container_type(PointVector{Float64}()) === Matrix{Float64}
    @test container_type(PointTensor{Float64}()) === Matrix{Float64}
    @test container_type(PointSymTensor{Float64}()) === Matrix{Float64}
    @test container_type(PointField{7,Float64}()) === Matrix{Float64}
    @test container_type(BondVector{Float64}()) === Matrix{Float64}
    @test container_type(BondSymTensor{Float64}()) === Matrix{Float64}

    # an unpinned shape uses the default float type as its concrete container
    @test container_type(PointVector{SimFloat}()) === Matrix{Float64}
    @test container_type(PointScalar{SimFloat}()) === Vector{Float64}

    # with an element type expression it returns the type expression of the struct field
    @test container_type(PointVector{SimFloat}(), :FT) == :(Base.Matrix{FT})
    @test container_type(PointScalar{SimFloat}(), :FT) == :(Base.Vector{FT})

    @test storage_field_type(PointVector) === Matrix{Float64}
    @test storage_field_type(PointScalar{Bool}) === Vector{Bool}
end

@testitem "alloc_field: size, element type and initial value of a shaped field" begin
    using Peridynamics.LinearAlgebra
    import Peridynamics: PointScalar, PointVector, PointTensor, PointSymTensor, BondScalar,
                         BondVector, BondTensor, BondSymTensor, LocalPoints, HaloPoints,
                         alloc_field, get_n_loc_points, get_n_points, get_n_bonds

    position = zeros(3, 10)
    position[1, :] = 0.0:9.0
    body = Body(BBMaterial(), position, ones(10))
    material!(body; horizon=1.5, rho=1.0, E=1.0, nu=0.25, Gc=1.0)
    pd = Peridynamics.PointDecomposition(body, 2)
    system = Peridynamics.get_system(body, pd, 1)
    n_loc, n_pts, n_bonds = get_n_loc_points(system), get_n_points(system),
                            get_n_bonds(system)
    @test n_loc < n_pts # the chunk really has halo points

    # the annotation decides the number of entries of a point field
    @test size(alloc_field(PointVector{Float64}(), system, LocalPoints())) == (3, n_loc)
    @test size(alloc_field(PointVector{Float64}(), system, HaloPoints())) == (3, n_pts)
    @test size(alloc_field(PointScalar{Float64}(), system, LocalPoints())) == (n_loc,)
    @test size(alloc_field(PointScalar{Float64}(), system, HaloPoints())) == (n_pts,)

    # the shape decides the number of rows
    @test size(alloc_field(PointTensor{Float64}(), system, LocalPoints())) == (9, n_loc)
    @test size(alloc_field(PointSymTensor{Float64}(), system, LocalPoints())) == (6, n_loc)

    # bond fields ignore the extent
    @test size(alloc_field(BondScalar{Float64}(), system, LocalPoints())) == (n_bonds,)
    @test size(alloc_field(BondVector{Float64}(), system, HaloPoints())) == (3, n_bonds)
    @test size(alloc_field(BondTensor{Float64}(), system, LocalPoints())) == (9, n_bonds)
    @test size(alloc_field(BondSymTensor{Float64}(), system, LocalPoints())) == (6, n_bonds)

    # element types and initial values
    @test eltype(alloc_field(PointScalar{Bool}(), system, LocalPoints())) === Bool
    @test eltype(alloc_field(PointScalar{Int}(), system, LocalPoints())) === Int
    @test all(iszero, alloc_field(PointVector{Float64}(), system, LocalPoints()))
    @test all(alloc_field(PointScalar{Bool}(), system, LocalPoints(), true))
    @test all(!, alloc_field(PointScalar{Bool}(), system, LocalPoints()))

    # a `UniformScaling` writes that tensor into every column
    R = alloc_field(BondTensor{Float64}(), system, LocalPoints(), I)
    @test size(R) == (9, n_bonds)
    @test all(isone, R[[1, 5, 9], :])
    @test sum(R) == 3 * n_bonds
    R2 = alloc_field(BondTensor{Float64}(), system, LocalPoints(), 2I)
    @test all(==(2.0), R2[[1, 5, 9], :])
    @test sum(R2) == 6 * n_bonds

    # a `UniformScaling` is only meaningful for a square tensor
    @test_throws ArgumentError alloc_field(BondScalar{Float64}(), system, LocalPoints(), I)
    @test_throws ArgumentError alloc_field(BondSymTensor{Float64}(), system, LocalPoints(),
                                           I)

    # everything that is neither a number nor a `UniformScaling` needs an `init_field`
    err = try
        alloc_field(PointVector{Float64}(), system, LocalPoints(), "one")
    catch e
        e
    end
    @test err isa ArgumentError
    @test contains(err.msg, "init_field")
end

@testitem "alloc_field: an unpinned shape follows the float type of the system" begin
    import Peridynamics: PointScalar, PointVector, LocalPoints, HaloPoints, DofVector,
                         alloc_field, float_type, default_float_type

    # a stub system with a different float type: an unpinned shape has to follow it, while
    # a shape with an explicit element type keeps it for every simulation
    struct Float32System <: Peridynamics.AbstractSystem end
    Peridynamics.float_type(::Float32System) = Float32
    Peridynamics.get_n_loc_points(::Float32System) = 5
    Peridynamics.get_n_points(::Float32System) = 7
    system = Float32System()
    @test default_float_type() === Float64
    @test float_type(system) === Float32

    @test eltype(alloc_field(PointVector(), system, LocalPoints())) === Float32
    @test eltype(alloc_field(PointScalar(), system, LocalPoints())) === Float32
    @test eltype(alloc_field(DofVector(), system, LocalPoints())) === Float32
    @test eltype(alloc_field(PointVector{Float64}(), system, LocalPoints())) === Float64
    @test eltype(alloc_field(PointScalar{Bool}(), system, LocalPoints())) === Bool

    # dof fields have one entry per degree of freedom
    @test size(alloc_field(DofVector(), system, LocalPoints())) == (3 * 5,)
    @test size(alloc_field(DofVector(), system, HaloPoints())) == (3 * 7,)

end

@testitem "FullField / EmptyField: how a time solver claims a field" begin
    import Peridynamics: FullField, EmptyField, LocalPoints, HaloPoints, PointVector,
                         DofVector, alloc_solver_field, alloc_empty_field, alloc_empty_array

    struct MarkerSystem <: Peridynamics.AbstractSystem end
    Peridynamics.get_n_loc_points(::MarkerSystem) = 5
    Peridynamics.get_n_points(::MarkerSystem) = 7
    system = MarkerSystem()

    # `FullField()` allocates at the extent of the declaration
    @test size(alloc_solver_field(FullField(), PointVector(), system, LocalPoints(),
                                  nothing)) == (3, 5)
    @test size(alloc_solver_field(FullField(), PointVector(), system, HaloPoints(),
                                  nothing)) == (3, 7)
    @test size(alloc_solver_field(FullField(), DofVector(), system, LocalPoints(),
                                  nothing)) == (3 * 5,)

    # ... and `FullField(extent)` overrules it, which is what `NewtonKrylov` needs for
    # `b_int`
    @test size(alloc_solver_field(FullField(HaloPoints()), PointVector(), system,
                                  LocalPoints(), nothing)) == (3, 7)

    # the initial value of the declaration is kept
    @test alloc_solver_field(FullField(), PointVector(), system, LocalPoints(), 2.0) ==
          fill(2.0, 3, 5)

    # an empty field has the type the field has in the storage, so it follows the float type
    # of the simulation and works for a field declared with a container type as well
    S64 = Peridynamics.storage_type(BBMaterial())
    S32 = Peridynamics.storage_type(BBMaterial(), Float32)
    @test alloc_empty_field(S64, Val(:velocity)) isa Matrix{Float64}
    @test size(alloc_empty_field(S64, Val(:velocity))) == (0, 0)
    @test alloc_empty_field(S64, Val(:residual)) isa Vector{Float64}
    @test size(alloc_empty_field(S64, Val(:residual))) == (0,)
    @test alloc_empty_field(S32, Val(:velocity)) isa Matrix{Float32}
    @test alloc_empty_array(Array{Float32,3}) isa Array{Float32,3}
    @test size(alloc_empty_array(Array{Float32,3})) == (0, 0, 0)
    @test_throws ArgumentError alloc_empty_array(Float64)
end

@testitem "alloc_field: the shapes follow the spatial dimension of the system" begin
    using Peridynamics.LinearAlgebra
    import Peridynamics: PointScalar, PointVector, PointTensor, PointSymTensor, BondVector,
                         BondTensor, LocalPoints, HaloPoints, alloc_field, get_n_dim

    # a stub system that answers with a different spatial dimension: everything a field
    # shape allocates has to follow `get_n_dim`, so that a 2D system needs no new shapes
    struct Dim2System <: Peridynamics.AbstractSystem end
    Peridynamics.get_n_dim(::Dim2System) = 2
    Peridynamics.get_n_loc_points(::Dim2System) = 5
    Peridynamics.get_n_points(::Dim2System) = 7
    Peridynamics.get_n_bonds(::Dim2System) = 11
    system = Dim2System()
    @test get_n_dim(system) == 2

    @test size(alloc_field(PointScalar{Float64}(), system, LocalPoints())) == (5,)
    @test size(alloc_field(PointVector{Float64}(), system, LocalPoints())) == (2, 5)
    @test size(alloc_field(PointVector{Float64}(), system, HaloPoints())) == (2, 7)
    @test size(alloc_field(PointTensor{Float64}(), system, LocalPoints())) == (4, 5)
    @test size(alloc_field(PointSymTensor{Float64}(), system, LocalPoints())) == (3, 5)
    @test size(alloc_field(BondVector{Float64}(), system, LocalPoints())) == (2, 11)
    @test size(alloc_field(BondTensor{Float64}(), system, LocalPoints())) == (4, 11)

    # dof fields follow the dimension too, `PointField` deliberately does not
    @test size(alloc_field(Peridynamics.DofVector{Float64}(), system, LocalPoints())) ==
          (2 * 5,)
    @test size(alloc_field(Peridynamics.PointField{7,Float64}(), system, LocalPoints())) ==
          (7, 5)

    # the identity tensor has to follow the dimension as well
    R = alloc_field(BondTensor{Float64}(), system, LocalPoints(), I)
    @test all(isone, R[[1, 4], :])
    @test sum(R) == 2 * 11
end

@testitem "@storage_fields / @inherit: blocks are merged in declaration order" begin
    import Peridynamics: PointScalar, PointVector, storage_fields_expr, @storage_fields

    Peridynamics.@storage_fields FieldsA begin
        @lth position::PointVector
        displacement::PointVector
        b_int::PointVector
    end

    Peridynamics.@storage_fields FieldsB begin
        displacement::PointVector
        b_int::PointVector
        acceleration::PointVector
    end

    # the marker type is what `@inherit` resolves
    @test FieldsA <: Peridynamics.AbstractStorageFields
    @test [d.name for d in storage_fields_expr(FieldsA)] ==
          [:position, :displacement, :b_int]

    # two blocks may contribute the same field if they declare it identically
    Peridynamics.@storage_fields FieldsC begin
        @inherit FieldsA
        @inherit FieldsB
        damage::PointScalar
    end
    @test [d.name for d in storage_fields_expr(FieldsC)] ==
          [:position, :displacement, :b_int, :acceleration, :damage]

    # the declared type is resolved, the shape is stored with it
    decl = storage_fields_expr(FieldsC)[1]
    @test decl.name === :position
    @test decl.annotation === :lth
    @test decl.shape === PointVector()
    @test isnothing(decl.init)

    # a field of the definition itself overrides an inherited field in place
    Peridynamics.@storage_fields FieldsD begin
        @inherit FieldsA
        @htl b_int::PointVector
        damage::PointScalar
    end
    declsD = storage_fields_expr(FieldsD)
    @test [d.name for d in declsD] == [:position, :displacement, :b_int, :damage]
    @test declsD[3].annotation === :htl

    # chained inheritance keeps working
    Peridynamics.@storage_fields FieldsE begin
        @inherit FieldsD
        velocity::PointVector
    end
    @test [d.name for d in storage_fields_expr(FieldsE)] ==
          [:position, :displacement, :b_int, :damage, :velocity]
end

@testitem "@storage_fields / @inherit: conflicting and unresolvable blocks are errors" begin
    import Peridynamics: PointVector, storage_fields_expr

    Peridynamics.@storage_fields ConflictA begin
        b_int::PointVector
    end
    Peridynamics.@storage_fields ConflictB begin
        @htl b_int::PointVector
    end

    # two blocks that declare the same field differently is an error
    err = try
        @eval Peridynamics.@storage_fields ConflictC begin
            @inherit ConflictA
            @inherit ConflictB
        end
        nothing
    catch e
        e
    end
    @test !isnothing(err)
    @test contains(sprint(showerror, err), "conflicting declarations")

    # a name that cannot be resolved names the fix
    err = try
        @eval Peridynamics.@storage_fields BadInherit begin
            @inherit ThisBlockDoesNotExist
        end
        nothing
    catch e
        e
    end
    @test !isnothing(err)
    @test contains(sprint(showerror, err), "cannot resolve")

    # only storages and field blocks can be inherited from
    @test_throws ArgumentError storage_fields_expr(Int)

    # `@inherit` outside of a storage definition is an error
    @test_throws ArgumentError Peridynamics.@inherit ConflictA
end

@testitem "@storage: a storage of shaped fields allocates without init_field" begin
    import Peridynamics: BondScalar, BondTensor, PointScalar, PointVector, PointTensor,
                         AbstractBondSystemMaterial, NoCorrection, CriticalStretch,
                         StandardPointParameters, get_n_loc_points, get_n_points,
                         get_n_bonds, storage_fields_expr

    struct ShapeMat <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::CriticalStretch
    end
    ShapeMat() = ShapeMat(CriticalStretch())
    Peridynamics.@params ShapeMat StandardPointParameters

    Peridynamics.@storage ShapeMat struct ShapeStorage
        @inherit Peridynamics.VelocityVerletFields
        @inherit Peridynamics.BondFracFields
        @htl b_int::PointVector
        stress::PointTensor
        flag::PointScalar{Bool} = true
        bond_state::BondTensor
        bond_scalar::BondScalar
    end

    function Peridynamics.force_density_point!(storage::ShapeStorage, system, ::ShapeMat,
                                               params, t, Δt, i)
        return nothing
    end

    # the container types follow the shapes
    S = Peridynamics.storage_type(ShapeMat())
    @test S <: ShapeStorage
    @test fieldtype(S, :position) === Matrix{Float64}
    @test fieldtype(S, :b_int) === Matrix{Float64}
    @test fieldtype(S, :stress) === Matrix{Float64}
    @test fieldtype(S, :flag) === Vector{Bool}
    @test fieldtype(S, :bond_state) === Matrix{Float64}
    @test fieldtype(S, :bond_scalar) === Vector{Float64}
    @test fieldtype(S, :n_active_bonds) === Vector{Int}

    # the storage allocates without a single `init_field` method
    position = zeros(3, 10)
    position[1, :] = 0.0:9.0
    body = Body(ShapeMat(), position, ones(10))
    material!(body; horizon=1.5, rho=1.0, E=1.0, nu=0.25, Gc=1.0)
    pd = Peridynamics.PointDecomposition(body, 2)
    ps = Peridynamics.get_param_spec(body)
    chunk = Peridynamics.BodyChunk(body, VelocityVerlet(steps=1), pd, 1, ps)
    s, system = chunk.storage, chunk.system
    n_loc, n_pts, n_bonds = get_n_loc_points(system), get_n_points(system),
                            get_n_bonds(system)

    # the annotations survive the inheritance and the override
    @test Peridynamics.halo_to_loc_fields(s) == (:b_int,)
    @test Peridynamics.loc_to_halo_fields(s) == (:position,)
    @test Peridynamics.is_halo_field(s, Val(:b_int))
    @test !Peridynamics.is_halo_field(s, Val(:displacement))

    @test size(s.position) == (3, n_pts)
    @test s.position == system.position
    @test size(s.displacement) == (3, n_loc)
    @test size(s.b_int) == (3, n_pts) # the halo annotation wins over the time solver
    @test size(s.stress) == (9, n_loc)
    @test size(s.flag) == (n_loc,) && all(s.flag)
    @test size(s.bond_state) == (9, n_bonds)
    @test size(s.bond_scalar) == (n_bonds,)
    @test length(s.damage) == n_loc
    @test length(s.bond_active) == n_bonds && all(s.bond_active)
end

@testitem "@storage: a container-typed field keeps its type and needs init_field" begin
    import Peridynamics: AbstractBondSystemMaterial, AbstractTimeSolver, BondSystem,
                         NoCorrection, CriticalStretch, StandardPointParameters,
                         get_n_loc_points, get_n_bonds

    struct LegacyMat <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::CriticalStretch
    end
    LegacyMat() = LegacyMat(CriticalStretch())
    Peridynamics.@params LegacyMat StandardPointParameters

    # the syntax of every field declaration before field shapes existed
    Peridynamics.@storage LegacyMat struct LegacyStorage <: Peridynamics.AbstractStorage
        @lth position::PointVector{Float64}
        displacement::PointVector
        velocity::PointVector
        velocity_half::PointVector
        acceleration::PointVector
        b_int::PointVector
        b_ext::PointVector
        damage::PointScalar
        n_active_bonds::PointScalar{Int}
        bond_active::Vector{Bool}
        my_own_field::Vector{Float64}
    end

    function Peridynamics.init_field(::LegacyMat, ::AbstractTimeSolver, system::BondSystem,
                                     ::Val{:my_own_field})
        return zeros(get_n_bonds(system))
    end

    function Peridynamics.force_density_point!(storage::LegacyStorage, system, ::LegacyMat,
                                               params, t, Δt, i)
        return nothing
    end

    position = zeros(3, 10)
    position[1, :] = 0.0:9.0
    body = Body(LegacyMat(), position, ones(10))
    material!(body; horizon=1.5, rho=1.0, E=1.0, nu=0.25, Gc=1.0)
    pd = Peridynamics.PointDecomposition(body, 2)
    ps = Peridynamics.get_param_spec(body)
    chunk = Peridynamics.BodyChunk(body, VelocityVerlet(steps=1), pd, 1, ps)
    s, system = chunk.storage, chunk.system

    @test fieldtype(Peridynamics.storage_type(LegacyMat()), :position) === Matrix{Float64}
    @test size(s.displacement) == (3, get_n_loc_points(system))
    @test length(s.my_own_field) == get_n_bonds(system)
end

@testitem "@storage: a field a time solver sizes cannot be a container type" begin
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection, CriticalStretch,
                         StandardPointParameters

    struct UnshapedMat <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::CriticalStretch
    end
    UnshapedMat() = UnshapedMat(CriticalStretch())
    Peridynamics.@params UnshapedMat StandardPointParameters

    # `velocity` is a field of the Velocity Verlet solver, and the solver only says that it
    # needs it, so a container type leaves nobody who knows the size
    Peridynamics.@storage UnshapedMat struct UnshapedStorage <: Peridynamics.AbstractStorage
        @lth position::PointVector{Float64}
        displacement::PointVector
        velocity::Matrix{Float64}
        velocity_half::PointVector
        acceleration::PointVector
        b_int::PointVector
        b_ext::PointVector
        damage::PointScalar
        n_active_bonds::PointScalar{Int}
        bond_active::BondScalar{Bool}
    end

    function Peridynamics.force_density_point!(storage::UnshapedStorage, system,
                                               ::UnshapedMat, params, t, Δt, i)
        return nothing
    end

    position = zeros(3, 10)
    position[1, :] = 0.0:9.0
    body = Body(UnshapedMat(), position, ones(10))
    material!(body; horizon=1.5, rho=1.0, E=1.0, nu=0.25, Gc=1.0)
    pd = Peridynamics.PointDecomposition(body, 2)
    ps = Peridynamics.get_param_spec(body)
    err = try
        Peridynamics.BodyChunk(body, VelocityVerlet(steps=1), pd, 1, ps)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test contains(err.msg, "velocity")
    @test contains(err.msg, "@inherit VelocityVerletFields")
end

@testitem "get_field_decls: invalid declarations are rejected" begin
    import Peridynamics: PointScalar, get_field_decls

    # an initial value needs a field shape that knows how to apply it
    block = quote
        my_field::Vector{Float64} = 1.0
    end
    @test_throws ArgumentError get_field_decls(block.args, @__MODULE__)

    # untyped fields are still rejected
    block = quote
        my_field
    end
    @test_throws ArgumentError get_field_decls(block.args, @__MODULE__)
    block = quote
        @lth my_field
    end
    @test_throws ArgumentError get_field_decls(block.args, @__MODULE__)

    # `@inherit` needs at least one storage or field block
    block = quote
        @inherit
    end
    @test_throws ArgumentError get_field_decls(block.args, @__MODULE__)

    # ... and every one of them has to be resolvable
    block = quote
        @inherit NoSuchFieldBlock1 NoSuchFieldBlock2
    end
    @test_throws ArgumentError get_field_decls(block.args, @__MODULE__)

    # a field shape works and carries the initial value
    block = quote
        my_field::PointScalar{Bool} = true
    end
    decls = get_field_decls(block.args, @__MODULE__)
    @test length(decls) == 1
    @test decls[1].name === :my_field
    @test decls[1].annotation === :none
    @test decls[1].shape === PointScalar{Bool}()
    @test decls[1].init === true
    # only the halo exchange is annotated, point data follows from the shape
    @test Peridynamics.is_point_decl(decls[1])
    @test !Peridynamics.is_halo_decl(decls[1])
end

@testitem "derive_storage_type_params: one parameter per element type and dimension" begin
    import Peridynamics: PointScalar, PointVector, BondScalar, DofVector, SimFloat,
                         derive_storage_type_params, get_field_decls, storage_param_bound,
                         storage_param_default, storage_field_decl_type

    block = quote
        @lth position::PointVector{Float64}
        displacement::PointVector
        damage::PointScalar
        n_active_bonds::PointScalar{Int}
        bond_active::BondScalar{Bool}
        residual::DofVector
        legacy::Vector{Float64}
        static::Peridynamics.StaticArrays.MArray{Tuple{3,3},Float64,2,9}
    end
    decls = get_field_decls(block.args, @__MODULE__)
    (; params, field_params, uses_sim_float) = derive_storage_type_params(decls)

    # one parameter per distinct element type and number of dimensions, in the order in
    # which they first appear, and `FT` for everything that follows the simulation
    @test uses_sim_float
    @test [p.name for p in params] == [:M_F64, :M_FT, :V_FT, :V_Int, :V_Bool, :V_F64]
    @test [p.eltype for p in params] == [Float64, SimFloat, SimFloat, Int, Bool, Float64]
    @test [p.n_dims for p in params] == [2, 2, 1, 1, 1, 1]

    # `residual` is a `Vector` of the float type and shares its parameter with `damage`,
    # `legacy` is a `Vector{Float64}` and shares none of them, because a pinned element
    # type is not the float type of the simulation
    @test field_params[:residual] === :V_FT
    @test field_params[:damage] === :V_FT
    @test field_params[:legacy] === :V_F64

    # a field that is not an `Array` keeps the type it is declared with
    @test !haskey(field_params, :static)
    @test storage_field_decl_type(decls[end], field_params) ===
          Peridynamics.StaticArrays.MArray{Tuple{3,3},Float64,2,9}
    @test storage_field_decl_type(decls[2], field_params) === :M_FT

    # the bound is abstract so that any array backend fits, the default is the CPU array
    @test storage_param_bound(params[2]) == Expr(:curly, AbstractMatrix, :FT)
    @test storage_param_bound(params[4]) == Expr(:curly, AbstractVector, Int)
    @test storage_param_default(params[2]) == Expr(:curly, Matrix, :FT)
    @test storage_param_default(params[4]) == Expr(:curly, Vector, Int)
end

@testitem "@storage: the storage is generic in the float type and adaptable" begin
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection, CriticalStretch,
                         StandardPointParameters, storage_type, get_storage

    # a minimal stand-in for the array type of another backend
    struct WrappedArray{T,N} <: AbstractArray{T,N}
        a::Array{T,N}
    end
    Base.size(x::WrappedArray) = size(x.a)
    Base.getindex(x::WrappedArray, i...) = getindex(x.a, i...)
    Base.setindex!(x::WrappedArray, v, i...) = setindex!(x.a, v, i...)
    struct WrappedBackend end
    function Peridynamics.Adapt.adapt_storage(::WrappedBackend, a::Array{T,N}) where {T,N}
        return WrappedArray{T,N}(a)
    end

    struct AdaptMat <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::CriticalStretch
    end
    AdaptMat() = AdaptMat(CriticalStretch())
    Peridynamics.@params AdaptMat StandardPointParameters

    Peridynamics.@storage AdaptMat struct AdaptStorage
        @inherit Peridynamics.VelocityVerletFields
        @inherit Peridynamics.BondFracFields
    end

    function Peridynamics.force_density_point!(storage::AdaptStorage, system, ::AdaptMat,
                                               params, t, Δt, i)
        return nothing
    end

    # `position` is pinned to `Float64` for every float type, everything else follows it
    S64, S32 = storage_type(AdaptMat()), storage_type(AdaptMat(), Float32)
    @test S64 <: AdaptStorage && S32 <: AdaptStorage
    @test fieldtype(S64, :displacement) === Matrix{Float64}
    @test fieldtype(S32, :displacement) === Matrix{Float32}
    @test fieldtype(S32, :damage) === Vector{Float32}
    @test fieldtype(S32, :position) === Matrix{Float64}
    @test fieldtype(S32, :n_active_bonds) === Vector{Int}
    @test fieldtype(S32, :bond_active) === Vector{Bool}

    position = zeros(3, 10)
    position[1, :] = 0.0:9.0
    body = Body(AdaptMat(), position, ones(10))
    material!(body; horizon=1.5, rho=1.0, E=1.0, nu=0.25, Gc=1.0)
    pd = Peridynamics.PointDecomposition(body, 2)
    ps = Peridynamics.get_param_spec(body)
    chunk = Peridynamics.BodyChunk(body, VelocityVerlet(steps=1), pd, 1, ps)
    mat, solver, system = body.mat, VelocityVerlet(steps=1), chunk.system

    # one literal `Val` per field makes the whole construction type stable
    @test typeof(chunk.storage) === S64
    @test (@inferred get_storage(mat, solver, system)) isa AdaptStorage

    # ... and the parametric struct makes the storage movable to another backend
    adapted = Peridynamics.Adapt.adapt(WrappedBackend(), chunk.storage)
    @test adapted isa AdaptStorage
    @test adapted.position isa WrappedArray{Float64,2}
    @test adapted.damage isa WrappedArray{Float64,1}
    @test adapted.bond_active isa WrappedArray{Bool,1}
    @test adapted.position == chunk.storage.position
end
