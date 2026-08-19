# `Body`: the point cloud with its point sets, material parameters, conditions and pre-cracks
# (`src/discretization/body.jl`). The bodies here are deterministic point clouds; a `Body` does
# not care where the points are until a system is built.

@testmodule BodyCase begin
    using Peridynamics
    "A body of `n` points on a line with unit volumes; no material yet."
    function line(mat=BBMaterial(), n=4)
        position = zeros(3, n)
        position[1, :] .= 0:(n - 1)
        return Body(mat, position, ones(n))
    end
    "The 4-point tetrahedron with unit volumes; no material yet."
    function tetra()
        position = [0.0 1.0 0.0 0.0
                    0.0 0.0 1.0 0.0
                    0.0 0.0 0.0 1.0]
        return Body(BBMaterial(), position, [1.0, 1.0, 1.0, 1.0])
    end
    """
    The 8-point body of the `log_msg_body` items: two parameter sets, two initial and two
    boundary conditions, a pre-crack, points without failure and a name. `matkwargs` are the
    `material!` keywords the material needs beyond `horizon`, `rho`, `E` and `Gc`.
    """
    function logged_body(mat; matkwargs...)
        position, volume = uniform_box(1, 1, 1, 0.5)
        body = Body(mat, position, volume)
        point_set!(body, :a, 1:2)
        material!(body; horizon=1, rho=1, E=1, Gc=1, matkwargs...)
        material!(body, :a; horizon=2, rho=2, E=2, Gc=2, matkwargs...)
        velocity_ic!(body, :a, :z, 1.0)
        velocity_ic!(p -> p[1] * 2.0, body, :a, :y)
        velocity_bc!(t -> t, body, :a, 1)
        forcedensity_bc!((p, t) -> p[1] + p[2] + p[3] + t, body, :a, 2)
        point_set!(body, :b, 3:4)
        precrack!(body, :a, :b)
        no_failure!(body, :a)
        Peridynamics.change_name!(body, :testbody)
        return body
    end
end

@testitem "Body: construction" setup=[BodyCase] begin
    for mat in (BBMaterial(), CKIMaterial())
        body = BodyCase.line(mat, 10)
        @test body.mat == mat
        @test body.n_points == n_points(body) == 10
        @test body.position[1, :] == 0:9
        @test body.volume == ones(10)
        @test body.fail_permit == fill(true, 10)
        @test isempty(body.single_dim_bcs) && isempty(body.posdep_single_dim_bcs)
        @test isempty(body.single_dim_ics) && isempty(body.posdep_single_dim_ics)
        @test isempty(body.point_sets_precracks)
        @test body.point_sets == Dict{Symbol,Vector{Int}}(:all_points => 1:10)
        @test eltype(body.point_params) === Peridynamics.point_param_type(mat)
        @test isempty(body.point_params)
        @test body.params_map == zeros(Int, 10)
    end
    # the material is the first argument, the point cloud must be a 3×n matrix with n volumes
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0] # only two rows
    @test_throws DimensionMismatch Body(BBMaterial(), position, ones(4))
    @test_throws DimensionMismatch Peridynamics.check_pos_and_vol(4, position, ones(4))
    @test_throws DimensionMismatch Peridynamics.check_pos_and_vol(3, position, ones(4))
    @test_throws ErrorException Body(BBMaterial(), zeros(3, 0), Float64[])
    @test_throws ErrorException Body(BBMaterial(), [NaN 1.0; 0.0 0.0; 0.0 0.0], ones(2))
    @test_throws ErrorException Body(BBMaterial(), [0.0 1.0; 0.0 0.0; 0.0 0.0], [1.0, NaN])
end

@testitem "point_set!: indices, predicates and errors" setup=[BodyCase] begin
    body = BodyCase.tetra()
    @test body.point_sets == Dict{Symbol,Vector{Int}}(:all_points => 1:4)
    point_set!(body, :a, 1:2)
    @test body.point_sets == Dict(:all_points => 1:4, :a => 1:2)
    # an existing name is an error
    @test_throws ErrorException point_set!(body, :a, 3:4)
    # a predicate on the position
    point_set!(x -> x > 0.5, body, :b)
    @test body.point_sets == Dict(:all_points => 1:4, :a => 1:2, :b => [2])
    # with do syntax
    point_set!(body, :c) do p
        p[3] > 0.0
    end
    @test body.point_sets == Dict(:all_points => 1:4, :a => 1:2, :b => [2], :c => [4])
    # indices outside the body
    @test_throws BoundsError point_set!(body, :d, 1:5)
    # the public accessor
    @test point_sets(body) === body.point_sets
    @test_throws ErrorException Peridynamics.check_if_set_is_defined(body.point_sets, :d)
end

@testitem "material!: parameters for the body and for point sets" setup=[BodyCase] begin
    body = BodyCase.tetra()
    point_set!(body, :a, 1:2)
    @test isempty(body.point_params) && body.params_map == zeros(Int, 4)
    # the whole body
    material!(body; horizon=1, E=1, rho=1, Gc=1)
    @test length(body.point_params) == 1
    p1 = body.point_params[1]
    @test (p1.δ, p1.E, p1.rho, p1.nu, p1.Gc) == (1.0, 1.0, 1.0, 0.25, 1.0)
    @test p1.G ≈ 0.4 && p1.K ≈ 2 / 3 # from E = 1 and nu = 1/4
    @test p1.εc ≈ sqrt(5 * p1.Gc / (9 * p1.K * p1.δ)) # critical stretch of a bond-based material
    @test body.params_map == [1, 1, 1, 1]
    # a point set gets its own parameters
    material!(body, :a; horizon=2, E=2, rho=2, Gc=2)
    @test length(body.point_params) == 2
    p2 = body.point_params[2]
    @test (p2.δ, p2.E, p2.rho, p2.Gc) == (2.0, 2.0, 2.0, 2.0)
    @test body.params_map == [2, 2, 1, 1]
    # the whole body again overwrites everything
    material!(body; horizon=3, E=3, rho=3, Gc=3)
    @test length(body.point_params) == 1
    @test body.point_params[1].δ == 3.0
    @test body.params_map == [1, 1, 1, 1]
end

@testitem "material!: with and without failure" setup=[BodyCase] begin
    body = BodyCase.tetra()
    # no critical stretch and no energy release rate: the points cannot fail
    for kwargs in ((; epsilon_c=0), (; Gc=0), (;))
        material!(body; horizon=1, E=1, rho=1, kwargs...)
        @test body.point_params[1].Gc == 0 && body.point_params[1].εc == 0
        @test body.fail_permit == zeros(Bool, 4)
    end
    # both given is ambiguous
    @test_throws ArgumentError material!(body; horizon=1, E=1, rho=1, epsilon_c=1, Gc=2)
    # overwriting toggles the failure permission
    material!(body; horizon=1, E=1, rho=1, Gc=1)
    @test all(body.fail_permit)
    material!(body; horizon=1, E=1, rho=1, Gc=0)
    @test !any(body.fail_permit)
    material!(body; horizon=1, E=1, rho=1, epsilon_c=1)
    @test all(body.fail_permit)
    # per point set
    point_set!(body, :a, 1:2)
    material!(body, :a; horizon=1, E=1, rho=1, Gc=0)
    @test body.fail_permit == [0, 0, 1, 1]
    material!(body; horizon=1, E=1, rho=1, Gc=0)
    @test body.fail_permit == [0, 0, 0, 0]
    material!(body, :a; horizon=1, E=1, rho=1, epsilon_c=1)
    @test body.fail_permit == [1, 1, 0, 0]
end

@testitem "no_failure!" setup=[BodyCase] begin
    body = BodyCase.tetra()
    point_set!(body, :a, 1:2)
    # needs material parameters for every point
    @test_throws ArgumentError no_failure!(body)
    material!(body, :a; horizon=1, E=1, rho=1, Gc=1)
    @test_throws ArgumentError no_failure!(body)
    # the whole body
    material!(body; horizon=1, E=1, rho=1, Gc=1)
    no_failure!(body)
    @test body.fail_permit == [0, 0, 0, 0]
    # a point set
    material!(body; horizon=1, E=1, rho=1, Gc=1)
    @test body.fail_permit == [1, 1, 1, 1]
    no_failure!(body, :a)
    @test body.fail_permit == [0, 0, 1, 1]
end

@testitem "velocity_bc! and forcedensity_bc!: declaration" setup=[Fixtures, BodyCase] begin
    # the dimension as an integer or a symbol, on different sets; the conditions are stored in
    # the order they are declared, with the field they act on
    for (bc!, field) in ((velocity_bc!, :velocity_half), (forcedensity_bc!, :b_ext))
        body = BodyCase.tetra()
        point_set!(body, :a, 1:2)
        point_set!(body, :b, 3:4)
        material!(body; horizon=1, E=1, rho=1, Gc=1)
        @test !Peridynamics.has_conditions(body)
        dims = [(:a, 1), (:a, 2), (:a, 3), (:b, :x), (:b, :y), (:b, :z)]
        for (i, (set, dim)) in enumerate(dims)
            bc!(Fixtures.f_one, body, set, dim)
            @test length(body.single_dim_bcs) == i
            bc = body.single_dim_bcs[i]
            @test all(bc.fun(t) == 1 for t in (-1, 0, 1, Inf, NaN))
            @test bc.field === field
            @test bc.point_set === set
            @test bc.dim == UInt8(dim isa Int ? dim : findfirst(==(dim), (:x, :y, :z)))
        end
        @test Peridynamics.has_conditions(body)
        # position dependent: f(p, t), on a fresh body (the dimensions above are taken)
        body = BodyCase.tetra()
        point_set!(body, :a, 1:2)
        material!(body; horizon=1, E=1, rho=1, Gc=1)
        bc!(Fixtures.f_pt, body, :a, 1)
        @test length(body.posdep_single_dim_bcs) == 1
        bc = body.posdep_single_dim_bcs[1]
        @test bc([2.0, 0.0, 0.0], 3.0) == 6.0 # p[1] * t
        @test bc.field === field && bc.point_set === :a && bc.dim == 0x01
    end
end

@testitem "Job: a body without conditions is rejected" setup=[BodyCase] begin
    body = BodyCase.tetra()
    material!(body; horizon=1, E=1, rho=1, Gc=1)
    @test !Peridynamics.has_conditions(body)
    @test_throws ErrorException Job(body, VelocityVerlet(steps=8))
end

@testitem "displacement bc with a solver that cannot apply it" begin
    # `displacement_bc!` adds a `PosSingleDimBC`, which only the `NewtonKrylov` solver applies.
    # With any other solver it used to be silently ignored.
    pos, vol = uniform_box(1, 1, 1, 0.25)
    make_body() = begin
        body = Body(BBMaterial(), pos, vol)
        material!(body; horizon=0.8, E=210e9, rho=8000, Gc=1)
        displacement_bc!(p -> 0.01 * p[1], body, :all_points, :x)
        body
    end

    @test_throws ArgumentError Job(make_body(), VelocityVerlet(steps=1))
    @test_throws ArgumentError Job(make_body(), DynamicRelaxation(steps=1))

    # the message has to say which solver is needed and which one was given
    err = try
        Job(make_body(), VelocityVerlet(steps=1))
    catch e
        e
    end
    @test contains(err.msg, "NewtonKrylov")
    @test contains(err.msg, "VelocityVerlet")

    # the same body with the solver that can apply the condition is fine
    @test Job(make_body(), NewtonKrylov(steps=1)) isa Job

    # a body whose load is prescribed with a condition the solver applies is unaffected
    body = Body(BBMaterial(), pos, vol)
    material!(body; horizon=0.8, E=210e9, rho=8000, Gc=1)
    velocity_bc!(t -> 1.0, body, :all_points, :x)
    @test Job(body, VelocityVerlet(steps=1)) isa Job
end

@testitem "show: Body" setup=[BodyCase] begin
    body = BodyCase.line(BBMaterial(), 10)
    compact(body) = sprint(show, MIME("text/plain"), body; context=:compact => true)
    full(body) = sprint(show, MIME("text/plain"), body; context=:compact => false)
    @test compact(body) == "10-point Body{BBMaterial}"
    @test contains(full(body), "1 point set(s):")
    @test contains(full(body), "10-point set `all_points`")
    point_set!(body, :a, 1:2)
    @test contains(full(body), "2 point set(s):") && contains(full(body), "2-point set `a`")
    material!(body, horizon=1, rho=1, E=1, Gc=1)
    @test contains(full(body), "1 point parameter(s):")
    @test contains(full(body), "δ=1.0, E=1.0, nu=0.25, rho=1.0, Gc=1.0")
    material!(body, :a, horizon=2, rho=2, E=2, Gc=2)
    @test !contains(full(body), "1 point parameter(s):")
    @test contains(full(body), "2 point parameter(s):")
    @test contains(full(body), "δ=2.0, E=2.0, nu=0.25, rho=2.0, Gc=2.0")
    velocity_ic!(body, :a, :z, 1.0)
    velocity_ic!(p -> p[1] * 2.0, body, :a, :y)
    velocity_bc!(t -> t, body, :a, 1)
    forcedensity_bc!((p, t) -> p[1] + p[2] + p[3] + t, body, :a, 2)
    msg = full(body)
    @test contains(msg, "2 boundary condition(s):")
    @test contains(msg, "BC on velocity: point_set=a, dim=1")
    @test contains(msg, "Pos.-dep. BC on force density: point_set=a, dim=2")
    @test contains(msg, "2 initial condition(s):")
    @test contains(msg, "IC on velocity: point_set=a, dim=3")
    @test contains(msg, "Pos.-dep. IC on velocity: point_set=a, dim=2")
    point_set!(body, :b, 3:4)
    precrack!(body, :a, :b)
    @test contains(full(body), "1 predefined crack(s)")
    no_failure!(body, :a)
    @test contains(full(body), "2 points with failure prohibited")
    Peridynamics.change_name!(body, :testbody)
    @test compact(body) == "10-point Body{BBMaterial} with name `testbody`"
    @test contains(full(body), "with name `testbody`")
    # the compact form never changes with the content
    @test compact(BodyCase.line(BBMaterial(), 10)) == "10-point Body{BBMaterial}"
end

# The body description written to the logfile, one item per material type because the
# material properties block differs. The bodies are built by `BodyCase.logged_body`.

@testitem "log_msg_body: BBMaterial" setup=[BodyCase] begin
    body = BodyCase.logged_body(BBMaterial())
    @test Peridynamics.log_msg_body(body) == """
        BODY `testbody`
          POINT CLOUD
            number of points ........................................................... 8
            min, max values x-direction ...................................... -0.25, 0.25
            min, max values y-direction ...................................... -0.25, 0.25
            min, max values z-direction ...................................... -0.25, 0.25
          POINT SETS
            number of points in set `a` ................................................ 2
            number of points in set `all_points` ....................................... 8
            number of points in set `b` ................................................ 2
          INITIAL CONDITIONS
            velocity condition ...................................... set `a`, dimension 3
            velocity condition ...................................... set `a`, dimension 2
          BOUNDARY CONDITIONS
            velocity condition ...................................... set `a`, dimension 1
            force density condition ................................. set `a`, dimension 2
          MATERIAL
            material type ..................................................... BBMaterial
            correction type .................................... Peridynamics.NoCorrection
            damage model type ............................... Peridynamics.CriticalStretch
            MATERIAL PROPERTIES #1
              horizon .................................................................. 1
              density .................................................................. 1
              Young's modulus .......................................................... 1
              Poisson's ratio ....................................................... 0.25
              shear modulus .......................................................... 0.4
              bulk modulus ..................................................... 0.6666667
              critical energy release rate ............................................. 1
              critical stretch ................................................. 0.9128709
            MATERIAL PROPERTIES #2
              horizon .................................................................. 2
              density .................................................................. 2
              Young's modulus .......................................................... 2
              Poisson's ratio ....................................................... 0.25
              shear modulus .......................................................... 0.8
              bulk modulus ...................................................... 1.333333
              critical energy release rate ............................................. 2
              critical stretch ................................................. 0.6454972
        """
end

@testitem "log_msg_body: OSBMaterial" setup=[BodyCase] begin
    body = BodyCase.logged_body(OSBMaterial(); nu=0.25)
    @test Peridynamics.log_msg_body(body) == """
        BODY `testbody`
          POINT CLOUD
            number of points ........................................................... 8
            min, max values x-direction ...................................... -0.25, 0.25
            min, max values y-direction ...................................... -0.25, 0.25
            min, max values z-direction ...................................... -0.25, 0.25
          POINT SETS
            number of points in set `a` ................................................ 2
            number of points in set `all_points` ....................................... 8
            number of points in set `b` ................................................ 2
          INITIAL CONDITIONS
            velocity condition ...................................... set `a`, dimension 3
            velocity condition ...................................... set `a`, dimension 2
          BOUNDARY CONDITIONS
            velocity condition ...................................... set `a`, dimension 1
            force density condition ................................. set `a`, dimension 2
          MATERIAL
            material type .................................................... OSBMaterial
            correction type .................................... Peridynamics.NoCorrection
            kernel function ................................................ linear_kernel
            damage model type ............................... Peridynamics.CriticalStretch
            MATERIAL PROPERTIES #1
              horizon .................................................................. 1
              density .................................................................. 1
              Young's modulus .......................................................... 1
              Poisson's ratio ....................................................... 0.25
              shear modulus .......................................................... 0.4
              bulk modulus ..................................................... 0.6666667
              critical energy release rate ............................................. 1
              critical stretch ................................................. 0.9128709
            MATERIAL PROPERTIES #2
              horizon .................................................................. 2
              density .................................................................. 2
              Young's modulus .......................................................... 2
              Poisson's ratio ....................................................... 0.25
              shear modulus .......................................................... 0.8
              bulk modulus ...................................................... 1.333333
              critical energy release rate ............................................. 2
              critical stretch ................................................. 0.6454972
        """
end

@testitem "log_msg_body: CMaterial" setup=[BodyCase] begin
    body = BodyCase.logged_body(CMaterial(); nu=0.25)
    @test Peridynamics.log_msg_body(body) == """
        BODY `testbody`
          POINT CLOUD
            number of points ........................................................... 8
            min, max values x-direction ...................................... -0.25, 0.25
            min, max values y-direction ...................................... -0.25, 0.25
            min, max values z-direction ...................................... -0.25, 0.25
          POINT SETS
            number of points in set `a` ................................................ 2
            number of points in set `all_points` ....................................... 8
            number of points in set `b` ................................................ 2
          INITIAL CONDITIONS
            velocity condition ...................................... set `a`, dimension 3
            velocity condition ...................................... set `a`, dimension 2
          BOUNDARY CONDITIONS
            velocity condition ...................................... set `a`, dimension 1
            force density condition ................................. set `a`, dimension 2
          MATERIAL
            material type ...................................................... CMaterial
            kernel function ................................................ linear_kernel
            constitutive model ....................... Peridynamics.SaintVenantKirchhoff()
            zero-energy mode stabilization .................. Peridynamics.ZEMSilling(0.5)
            damage model type ............................... Peridynamics.CriticalStretch
            maximum damage .......................................................... 0.85
            MATERIAL PROPERTIES #1
              horizon .................................................................. 1
              density .................................................................. 1
              Young's modulus .......................................................... 1
              Poisson's ratio ....................................................... 0.25
              shear modulus .......................................................... 0.4
              bulk modulus ..................................................... 0.6666667
              critical energy release rate ............................................. 1
              critical stretch ................................................. 0.9128709
            MATERIAL PROPERTIES #2
              horizon .................................................................. 2
              density .................................................................. 2
              Young's modulus .......................................................... 2
              Poisson's ratio ....................................................... 0.25
              shear modulus .......................................................... 0.8
              bulk modulus ...................................................... 1.333333
              critical energy release rate ............................................. 2
              critical stretch ................................................. 0.6454972
        """
end

@testitem "log_msg_body: BACMaterial" setup=[BodyCase] begin
    body = BodyCase.logged_body(BACMaterial(); nu=0.25)
    @test Peridynamics.log_msg_body(body) == """
        BODY `testbody`
          POINT CLOUD
            number of points ........................................................... 8
            min, max values x-direction ...................................... -0.25, 0.25
            min, max values y-direction ...................................... -0.25, 0.25
            min, max values z-direction ...................................... -0.25, 0.25
          POINT SETS
            number of points in set `a` ................................................ 2
            number of points in set `all_points` ....................................... 8
            number of points in set `b` ................................................ 2
          INITIAL CONDITIONS
            velocity condition ...................................... set `a`, dimension 3
            velocity condition ...................................... set `a`, dimension 2
          BOUNDARY CONDITIONS
            velocity condition ...................................... set `a`, dimension 1
            force density condition ................................. set `a`, dimension 2
          MATERIAL
            material type .................................................... BACMaterial
            kernel function ................................................ linear_kernel
            constitutive model ....................... Peridynamics.SaintVenantKirchhoff()
            damage model type ............................... Peridynamics.CriticalStretch
            maximum damage .......................................................... 0.85
            MATERIAL PROPERTIES #1
              horizon .................................................................. 1
              bond horizon ............................................................. 1
              density .................................................................. 1
              Young's modulus .......................................................... 1
              Poisson's ratio ....................................................... 0.25
              shear modulus .......................................................... 0.4
              bulk modulus ..................................................... 0.6666667
              critical energy release rate ............................................. 1
              critical stretch ................................................. 0.9128709
            MATERIAL PROPERTIES #2
              horizon .................................................................. 2
              bond horizon ............................................................. 2
              density .................................................................. 2
              Young's modulus .......................................................... 2
              Poisson's ratio ....................................................... 0.25
              shear modulus .......................................................... 0.8
              bulk modulus ...................................................... 1.333333
              critical energy release rate ............................................. 2
              critical stretch ................................................. 0.6454972
        """
end

@testitem "log_msg_body: CKIMaterial" setup=[BodyCase] begin
    body = BodyCase.logged_body(CKIMaterial(); nu=0.25)
    @test Peridynamics.log_msg_body(body) == """
        BODY `testbody`
          POINT CLOUD
            number of points ........................................................... 8
            min, max values x-direction ...................................... -0.25, 0.25
            min, max values y-direction ...................................... -0.25, 0.25
            min, max values z-direction ...................................... -0.25, 0.25
          POINT SETS
            number of points in set `a` ................................................ 2
            number of points in set `all_points` ....................................... 8
            number of points in set `b` ................................................ 2
          INITIAL CONDITIONS
            velocity condition ...................................... set `a`, dimension 3
            velocity condition ...................................... set `a`, dimension 2
          BOUNDARY CONDITIONS
            velocity condition ...................................... set `a`, dimension 1
            force density condition ................................. set `a`, dimension 2
          MATERIAL
            material type .................................................... CKIMaterial
            damage model type ............................... Peridynamics.CriticalStretch
            MATERIAL PROPERTIES #1
              horizon .................................................................. 1
              density .................................................................. 1
              Young's modulus .......................................................... 1
              Poisson's ratio ....................................................... 0.25
              shear modulus .......................................................... 0.4
              bulk modulus ..................................................... 0.6666667
              critical energy release rate ............................................. 1
              critical stretch ................................................. 0.9128709
              parameter one-neighbor interactions ............................... 3.819719
              parameter two-neighbor interactions ...................................... 0
              parameter three-neighbor interactions .................................... 0
            MATERIAL PROPERTIES #2
              horizon .................................................................. 2
              density .................................................................. 2
              Young's modulus .......................................................... 2
              Poisson's ratio ....................................................... 0.25
              shear modulus .......................................................... 0.8
              bulk modulus ...................................................... 1.333333
              critical energy release rate ............................................. 2
              critical stretch ................................................. 0.6454972
              parameter one-neighbor interactions .............................. 0.4774648
              parameter two-neighbor interactions ...................................... 0
              parameter three-neighbor interactions .................................... 0
        """
end

@testitem "Body: from an Abaqus inp file" begin
    file = joinpath(@__DIR__, "..", "AbaqusMeshConverter", "models", "CubeC3D8.inp")
    body = Body(BBMaterial(), file)
    @test size(body.position) == (3, 125)
    @test length(body.volume) == 125
    @test n_points(body) == 125
    @test body.volume ≈ fill(4^3, 125)
    sets = point_sets(body)
    @test sets[:l] == 101:125
    @test sets[:r] == 1:25
end
