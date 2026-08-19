# Boundary conditions (`src/physics/boundary_conditions.jl`): declaration on a body, conflicts,
# argument checks and their application to a chunk.
#
# The condition functions specialize `velocity_bc!` & co. on the function type, so the items
# use the named functions of `Fixtures` (`f_t`, `f_pt`, ...) instead of a fresh lambda per
# call.

@testmodule ConditionCase begin
    using Peridynamics
    "A two-point body with the sets `:a` (both points) and `:b` (the second point)."
    function two_points()
        body = Body(BBMaterial(), [0.0 1.0; 0.0 0.0; 0.0 0.0], [1.0, 1.0])
        material!(body, horizon=1.5, rho=1, E=210e9, Gc=1.0)
        point_set!(body, :a, 1:2)
        point_set!(body, :b, [2])
        return body
    end
end

@testitem "boundary conditions: applied to the chunk" setup=[Fixtures, ConditionCase] begin
    body = ConditionCase.two_points()
    velocity_bc!(Fixtures.f_t, body, :a, :x)      # f(t)
    velocity_bc!(Fixtures.f_pt, body, :a, :y)     # f(p, t)
    forcedensity_bc!(Fixtures.f_t, body, :a, :x)
    forcedensity_bc!(Fixtures.f_pt, body, :a, :y)
    velocity_ic!(body, :a, :x, 1.234)
    velocity_ic!(p -> p[1] * 3.456, body, :a, :y)
    chunk = Fixtures.handler(body, VelocityVerlet(steps=1); init=false).chunks[1]
    (; velocity_half, velocity, b_ext) = chunk.storage
    # the initial conditions are applied at construction
    @test velocity ≈ [1.234 1.234; 0.0 3.456; 0.0 0.0]
    # the boundary conditions at a time
    t = 2.345
    Peridynamics.apply_boundary_conditions!(chunk, t)
    @test velocity_half ≈ [t t; 0.0 1.0 * t; 0.0 0.0] # f_pt = p[1] * t with p[1] = 0, 1
    @test b_ext ≈ [t t; 0.0 1.0 * t; 0.0 0.0]
end

@testitem "displacement_bc!: applied incrementally" setup=[Fixtures, ConditionCase] begin
    body = ConditionCase.two_points()
    displacement_bc!(Fixtures.f_p, body, :a, :y) # f(p) = 2
    ts = DynamicRelaxation(; steps=10, stepsize=1.0)
    chunk = Fixtures.handler(body, ts; init=false).chunks[1]
    (; displacement) = chunk.storage
    @test iszero(displacement)
    # the prescribed displacement is reached at the end of the time, half of it at half the time
    Peridynamics.apply_incr_boundary_conditions!(chunk, 0.5)
    @test displacement ≈ [0.0 0.0; 1.0 1.0; 0.0 0.0]
    Peridynamics.apply_incr_boundary_conditions!(chunk, 1.0)
    @test displacement ≈ [0.0 0.0; 2.0 2.0; 0.0 0.0]
end

@testitem "boundary conditions: conflicts on the same points and dimension" setup=[Fixtures, ConditionCase] begin
    F = Fixtures
    body = ConditionCase.two_points()
    # a second condition on the same field, set and dimension, or on an overlapping set
    velocity_bc!(F.f_t, body, :a, :x)
    @test_throws ArgumentError velocity_bc!(F.f_2t, body, :a, :x)
    @test_throws ArgumentError velocity_bc!(F.f_pt, body, :a, :x)
    @test_throws ArgumentError velocity_bc!(F.f_pt, body, :b, :x)
    velocity_bc!(F.f_pt, body, :a, :y)
    @test_throws ArgumentError velocity_bc!(F.f_pt2, body, :a, :y)
    @test_throws ArgumentError velocity_bc!(F.f_t, body, :a, :y)
    @test_throws ArgumentError velocity_bc!(F.f_t, body, :b, :y)
    forcedensity_bc!(F.f_t, body, :a, :x)
    @test_throws ArgumentError forcedensity_bc!(F.f_2t, body, :a, :x)
    @test_throws ArgumentError forcedensity_bc!(F.f_pt, body, :a, :x)
    @test_throws ArgumentError forcedensity_bc!(F.f_pt, body, :b, :x)
    forcedensity_bc!(F.f_pt, body, :a, :y)
    @test_throws ArgumentError forcedensity_bc!(F.f_pt2, body, :a, :y)
    @test_throws ArgumentError forcedensity_bc!(F.f_t, body, :a, :y)
    @test_throws ArgumentError forcedensity_bc!(F.f_t, body, :b, :y)
    displacement_bc!(F.f_p, body, :a, :y)
    @test_throws ArgumentError displacement_bc!(F.f_p1, body, :a, :y)
    @test_throws ArgumentError displacement_bc!(F.f_p1, body, :b, :y)
    # the other dimensions and the other set are free
    velocity_bc!(F.f_t, body, :a, :z)
    forcedensity_bc!(F.f_t, body, :a, :z)
end

@testitem "boundary conditions: argument checks" setup=[Fixtures, ConditionCase] begin
    F = Fixtures
    body = ConditionCase.two_points()
    # the condition function must be f(t), f(p, t) or, for displacements, f(p)
    @test_throws ArgumentError velocity_bc!(F.f_bad_ab, body, :a, :z)
    @test_throws ArgumentError velocity_bc!(F.f_bad_ktu, body, :a, :z)
    @test_throws ArgumentError forcedensity_bc!(F.f_bad_ab, body, :a, :z)
    @test_throws ArgumentError forcedensity_bc!(F.f_bad_ktu, body, :a, :z)
    @test_throws ArgumentError displacement_bc!(F.f_bad_k, body, :a, :z)
    @test_throws ArgumentError displacement_bc!(F.f_bad_ab, body, :a, :z)
    @test_throws ArgumentError displacement_bc!(F.f_t, body, :a, :z)
    @test_throws ArgumentError displacement_bc!(F.f_pt, body, :a, :z)
    # unknown dimensions
    @test_throws ArgumentError velocity_bc!(F.f_one, body, :a, :k)
    @test_throws ArgumentError velocity_bc!(F.f_one, body, :a, 4)
    # unknown point set
    @test_throws ErrorException velocity_bc!(F.f_one, body, :c, :x)
end

@testitem "data boundary conditions: declaration, conflicts and checks" setup=[Fixtures, ConditionCase] begin
    F = Fixtures
    body = ConditionCase.two_points()
    # the data has one row per dimension and one column per point of the body
    @test_throws ArgumentError Peridynamics.velocity_databc!(body, zeros(2, 3), :a, [1, 2])
    @test_throws ArgumentError Peridynamics.velocity_databc!(body, zeros(3, 2), :a, [1, 2])
    Peridynamics.velocity_databc!(body, ones(1, 2), :a, [1])
    # conflicts with standard conditions
    @test_throws ArgumentError velocity_bc!(F.f_2t, body, :a, :x)
    Peridynamics.velocity_databc!(body, ones(1, 2), :a, [:y])
    @test_throws ArgumentError velocity_bc!(F.f_t, body, :a, :y)
    @test_throws ArgumentError velocity_bc!(F.f_t, body, :b, :y)
    Peridynamics.forcedensity_databc!(body, ones(1, 2), :a, [:x])
    @test_throws ArgumentError Peridynamics.forcedensity_databc!(body, 2 * ones(1, 2), :a, [:x])
    @test_throws ArgumentError forcedensity_bc!(F.f_pt, body, :a, :x)
    @test_throws ArgumentError forcedensity_bc!(F.f_pt, body, :b, :x)
    Peridynamics.forcedensity_databc!(body, ones(2, 2), :a, [:y, :z])
    for dim in (:y, :z), set in (:a, :b)
        @test_throws ArgumentError forcedensity_bc!(F.f_t, body, set, dim)
    end
    # unknown and malformed dimensions
    @test_throws ArgumentError Peridynamics.velocity_databc!(body, ones(1, 2), :a, [:k])
    @test_throws ArgumentError Peridynamics.velocity_databc!(body, ones(1, 2), :a, [4])
    @test_throws MethodError Peridynamics.velocity_databc!(body, ones(1, 2), :a, 3)
    @test_throws MethodError Peridynamics.velocity_databc!(body, ones(2, 2), :a, [:x, 2])
    @test_throws ArgumentError Peridynamics.velocity_databc!(body, ones(4, 2), :a, [1, 2, 3, 4])
    # conflicts when the standard conditions come first
    body = ConditionCase.two_points()
    velocity_bc!(F.f_2t, body, :a, :x)
    forcedensity_bc!(F.f_t, body, :a, :y)
    @test_throws ArgumentError Peridynamics.velocity_databc!(body, ones(1, 2), :a, [1])
    @test_throws ArgumentError Peridynamics.forcedensity_databc!(body, ones(2, 2), :a, [:y, :z])
end

@testitem "data boundary conditions: applied to the chunk" setup=[Fixtures, ConditionCase] begin
    body = ConditionCase.two_points()
    v_data = [1.0 2.0; 3.0 4.0]   # velocities in y and z per point
    f_data = [5.0 6.0]            # force density in x per point
    Peridynamics.velocity_databc!(body, v_data, :a, [:y, :z])
    Peridynamics.forcedensity_databc!(body, f_data, :b, [:x]) # only point 2
    msg = sprint(show, MIME("text/plain"), body; context=:compact => false)
    @test contains(msg, "2 boundary condition(s):")
    @test contains(msg, "Data BC on velocity: point_set=a, dims=UInt8[0x02, 0x03]")
    @test contains(msg, "Data BC on force density: point_set=b, dims=UInt8[0x01]")
    chunk = Fixtures.handler(body, VelocityVerlet(steps=1); init=false).chunks[1]
    Peridynamics.apply_boundary_conditions!(chunk, 0.0)
    @test chunk.storage.velocity_half ≈ [0.0 0.0; 1.0 2.0; 3.0 4.0]
    @test chunk.storage.b_ext ≈ [0.0 6.0; 0.0 0.0; 0.0 0.0]
end

@testitem "displacement_bc!: show of the condition and of the body" setup=[Fixtures] begin
    body = Fixtures.tetra4()
    displacement_bc!(Fixtures.f_p, body, :all_points, :x)
    bc = only(body.pos_single_dim_bcs)
    msg = sprint(show, bc)
    @test startswith(msg, "Pos. BC on displacement: ")
    @test contains(msg, "point_set=all_points") && contains(msg, "dim=1")
    msg = sprint(show, MIME("text/plain"), body)
    @test contains(msg, "Pos. BC on displacement")
end
