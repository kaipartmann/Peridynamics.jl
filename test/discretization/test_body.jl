@testitem "Body BBMaterial" begin
    # setup
    n_points = 10
    mat, position, volume = BBMaterial(), rand(3, n_points), rand(n_points)
    body = Body(mat, position, volume)

    # test body creation
    @test body.mat == BBMaterial()
    @test body.n_points == n_points
    @test body.position == position
    @test body.volume == volume
    @test body.fail_permit == fill(true, n_points)
    @test body.single_dim_bcs == Vector{Peridynamics.SingleDimBC}()
    @test body.posdep_single_dim_bcs == Vector{Peridynamics.PosDepSingleDimBC}()
    @test body.single_dim_ics == Vector{Peridynamics.SingleDimIC}()
    @test body.posdep_single_dim_ics == Vector{Peridynamics.PosDepSingleDimIC}()
    @test body.point_sets_precracks == Vector{Peridynamics.PointSetsPreCrack}()
    @test body.point_sets == Dict{Symbol,Vector{Int}}(:all_points => 1:n_points)
    @test body.point_params == Vector{Peridynamics.BBPointParameters}()
    @test body.params_map == zeros(Int, n_points)
end

@testitem "Body CKIMaterial" begin
    # setup
    n_points = 10
    mat, position, volume = CKIMaterial(), rand(3, n_points), rand(n_points)
    body = Body(mat, position, volume)

    # test body creation
    @test body.mat == CKIMaterial()
    @test body.n_points == n_points
    @test body.position == position
    @test body.volume == volume
    @test body.fail_permit == fill(true, n_points)
    @test body.single_dim_bcs == Vector{Peridynamics.SingleDimBC}()
    @test body.single_dim_ics == Vector{Peridynamics.SingleDimIC}()
    @test body.point_sets_precracks == Vector{Peridynamics.PointSetsPreCrack}()
    @test body.point_sets == Dict{Symbol,Vector{Int}}(:all_points => 1:n_points)
    @test body.point_params == Vector{Peridynamics.BBPointParameters}()
    @test body.params_map == zeros(Int, n_points)
end

@testitem "point_set!" begin
    # setup
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1, 1, 1, 1]
    mat = BBMaterial()
    body = Body(mat, position, volume)

    # test body
    @test body.n_points == 4
    @test body.position == position
    @test body.point_sets == Dict{Symbol,Vector{Int}}(:all_points => 1:4)

    # add point set
    point_set!(body, :a, 1:2)
    @test body.point_sets == Dict(:all_points => 1:4, :a => 1:2)

    # add another point set via function definition
    point_set!(x -> x > 0.5, body, :b)
    @test body.point_sets == Dict(:all_points => 1:4, :a => 1:2, :b => [2])

    # add point set with do syntax
    point_set!(body, :c) do p
        p[3] > 0.0
    end
    @test body.point_sets == Dict(:all_points => 1:4, :a => 1:2, :b => [2], :c => [4])

    # point_set!
    @test_throws BoundsError point_set!(body, :d, 1:5)
end

@testitem "material!" begin
    # setup
    n_points = 4
    mat, position, volume = BBMaterial(), rand(3, n_points), rand(n_points)
    body = Body(mat, position, volume)

    # test body creation
    @test body.mat == BBMaterial()
    @test body.n_points == n_points
    @test body.point_params == Vector{Peridynamics.BBPointParameters}()
    @test body.params_map == zeros(Int, n_points)
    @test body.point_sets == Dict{Symbol,Vector{Int}}(:all_points => 1:n_points)

    # add point set
    point_set!(body, :a, 1:2)
    @test body.point_sets == Dict(:all_points => 1:n_points, :a => 1:2)

    # add material
    material!(body; horizon=1, E=1, rho=1, Gc=1)
    @test body.point_params == [
        Peridynamics.BBPointParameters(1.0, 1.0, 1.0, 0.25, 0.4, 0.6666666666666666, 0.4,
                                       0.4, 1.0, 0.9128709291752769, 3.819718634205488),
    ]
    @test body.params_map == [1, 1, 1, 1]

    # add material to set
    material!(body, :a; horizon=2, E=2, rho=2, Gc=2)
    @test body.point_params == [
        Peridynamics.BBPointParameters(1.0, 1.0, 1.0, 0.25, 0.4, 0.6666666666666666, 0.4,
                                       0.4, 1.0, 0.9128709291752769, 3.819718634205488),
        Peridynamics.BBPointParameters(2.0, 2.0, 2.0, 0.25, 0.8, 1.3333333333333333, 0.8,
                                       0.8, 2.0, 0.6454972243679028, 0.477464829275686),
    ]
    @test body.params_map == [2, 2, 1, 1]

    # add material to body -> overwriting everything!
    material!(body; horizon=3, E=3, rho=3, Gc=3)
    @test body.point_params == [
        Peridynamics.BBPointParameters(3.0, 3.0, 3.0, 0.25, 1.2, 2.0, 1.2, 1.2, 3.0,
                                       0.5270462766947299, 0.1414710605261292),
    ]
    @test body.params_map == [1, 1, 1, 1]
end

@testitem "material! without failure criteria" begin
    # setup
    n_points = 4
    mat, position, volume = BBMaterial(), rand(3, n_points), rand(n_points)
    body = Body(mat, position, volume)

    # if εc = 0
    material!(body; horizon=1, E=1, rho=1, epsilon_c=0)
    @test isapprox(body.point_params[1].Gc, 0; atol=eps())
    @test isapprox(body.point_params[1].εc, 0; atol=eps())
    @test isapprox(body.fail_permit, zeros(n_points); atol=eps())

    # if Gc = 0
    material!(body; horizon=1, E=1, rho=1, Gc=0)
    @test isapprox(body.point_params[1].Gc, 0; atol=eps())
    @test isapprox(body.point_params[1].εc, 0; atol=eps())
    @test isapprox(body.fail_permit, zeros(n_points); atol=eps())

    # if εc and Gc are both not defined
    material!(body; horizon=1, E=1, rho=1, )
    @test isapprox(body.point_params[1].Gc, 0; atol=eps())
    @test isapprox(body.point_params[1].εc, 0; atol=eps())
    @test isapprox(body.fail_permit, zeros(n_points); atol=eps())

    # if εc and Gc are both defined
    @test_throws ArgumentError material!(body; horizon=1, E=1, rho=1, epsilon_c=1, Gc=2)

    # if material without failure is overwritten by material with failure
    material!(body; horizon=1, E=1, rho=1, )
    @test isapprox(body.fail_permit, zeros(n_points); atol=eps())
    material!(body; horizon=1, E=1, rho=1, Gc=1)
    @test isapprox(body.fail_permit, ones(n_points); atol=eps())
    material!(body; horizon=1, E=1, rho=1, Gc=0)
    @test isapprox(body.fail_permit, zeros(n_points); atol=eps())
    material!(body; horizon=1, E=1, rho=1, epsilon_c=1)
    @test isapprox(body.fail_permit, ones(n_points); atol=eps())

    # pointsets:
    point_set!(body, :a, 1:2)
    material!(body, :a; horizon=1, E=1, rho=1, Gc=0)
    @test body.fail_permit == [0, 0, 1, 1]

    material!(body; horizon=1, E=1, rho=1, Gc=0)
    @test body.fail_permit == [0, 0, 0, 0]

    material!(body, :a; horizon=1, E=1, rho=1, epsilon_c=1)
    @test body.fail_permit == [1, 1, 0, 0]
end

@testitem "velocity_bc!" begin
    # setup
    n_points = 4
    mat, position, volume = BBMaterial(), rand(3, n_points), rand(n_points)
    body = Body(mat, position, volume)

    # test body creation
    @test body.n_points == n_points
    @test body.position == position
    @test body.volume == volume
    @test body.point_sets == Dict{Symbol,Vector{Int}}(:all_points => 1:n_points)

    # add point set
    point_set!(body, :a, 1:2)
    point_set!(body, :b, 3:4)
    @test body.point_sets == Dict(:all_points => 1:n_points, :a => 1:2, :b => 3:4)

    # add material
    material!(body; horizon=1, E=1, rho=1, Gc=1)
    @test body.point_params == [
        Peridynamics.BBPointParameters(1.0, 1.0, 1.0, 0.25, 0.4, 0.6666666666666666, 0.4,
                                       0.4, 1.0, 0.9128709291752769, 3.819718634205488),
    ]
    @test body.params_map == [1, 1, 1, 1]

    # velocity bc 1
    velocity_bc!(t -> 1, body, :a, 1)
    @test length(body.single_dim_bcs) == 1
    bc1 = body.single_dim_bcs[1]
    for t in [-1, 0, 1, Inf, NaN]
        @test bc1.fun(t) == 1
    end
    @test bc1.field === :velocity_half
    @test bc1.point_set === :a
    @test bc1.dim == 0x01

    # velocity bc 2
    velocity_bc!(t -> 2, body, :a, 2)
    @test length(body.single_dim_bcs) == 2
    bc2 = body.single_dim_bcs[2]
    for t in [-1, 0, 1, Inf, NaN]
        @test bc2.fun(t) == 2
    end
    @test bc2.field === :velocity_half
    @test bc2.point_set === :a
    @test bc2.dim == 0x02

    # velocity bc 3
    velocity_bc!(t -> 3, body, :a, 3)
    @test length(body.single_dim_bcs) == 3
    bc3 = body.single_dim_bcs[3]
    for t in [-1, 0, 1, Inf, NaN]
        @test bc3.fun(t) == 3
    end
    @test bc3.field === :velocity_half
    @test bc3.point_set === :a
    @test bc3.dim == 0x03

    # velocity bc 4
    velocity_bc!(t -> 4, body, :b, :x)
    @test length(body.single_dim_bcs) == 4
    bc4 = body.single_dim_bcs[4]
    for t in [-1, 0, 1, Inf, NaN]
        @test bc4.fun(t) == 4
    end
    @test bc4.field === :velocity_half
    @test bc4.point_set === :b
    @test bc4.dim == 0x01

    # velocity bc 5
    velocity_bc!(t -> 5, body, :b, :y)
    @test length(body.single_dim_bcs) == 5
    bc5 = body.single_dim_bcs[5]
    for t in [-1, 0, 1, Inf, NaN]
        @test bc5.fun(t) == 5
    end
    @test bc5.field === :velocity_half
    @test bc5.point_set === :b
    @test bc5.dim == 0x02

    # velocity bc 6
    velocity_bc!(t -> 6, body, :b, :z)
    @test length(body.single_dim_bcs) == 6
    bc6 = body.single_dim_bcs[6]
    for t in [-1, 0, 1, Inf, NaN]
        @test bc6.fun(t) == 6
    end
    @test bc6.field === :velocity_half
    @test bc6.point_set === :b
    @test bc6.dim == 0x03
end

@testitem "velocity_bc! position dependend" begin
    # setup
    n_points = 4
    mat, position, volume = BBMaterial(), rand(3, n_points), rand(n_points)
    body = Body(mat, position, volume)

    # test body creation
    @test body.n_points == n_points
    @test body.position == position
    @test body.volume == volume
    @test body.point_sets == Dict{Symbol,Vector{Int}}(:all_points => 1:n_points)

    # add point set
    point_set!(body, :a, 1:2)
    point_set!(body, :b, 3:4)
    @test body.point_sets == Dict(:all_points => 1:n_points, :a => 1:2, :b => 3:4)

    # add material
    material!(body; horizon=1, E=1, rho=1, Gc=1)
    @test body.point_params == [
        Peridynamics.BBPointParameters(1.0, 1.0, 1.0, 0.25, 0.4, 0.6666666666666666, 0.4,
                                       0.4, 1.0, 0.9128709291752769, 3.819718634205488),
    ]
    @test body.params_map == [1, 1, 1, 1]

    # velocity bc 1
    velocity_bc!((p, t) -> 1, body, :a, 1)
    @test length(body.posdep_single_dim_bcs) == 1
    bc1 = body.posdep_single_dim_bcs[1]
    _p = [0, 0, 0]
    for t in [-1, 0, 1, Inf, NaN]
        @test bc1(_p, t) == 1
    end
    @test bc1.field === :velocity_half
    @test bc1.point_set === :a
    @test bc1.dim == 0x01


    # velocity bc 2
    velocity_bc!((p, t) -> p[1] + p[2] + p[3] + t, body, :a, 2)
    @test length(body.posdep_single_dim_bcs) == 2
    bc2 = body.posdep_single_dim_bcs[2]
    _p = [0, 0, 0]
    for t in [-1, 0, 1, Inf, NaN]
        @test bc1(_p, t) == 1
    end
    @test bc2.field === :velocity_half
    @test bc2.point_set === :a
    @test bc2.dim == 0x02
    @test bc2([1,2,3], 1.0) == 7.0
end

@testitem "forcedensity_bc!" begin
    # setup
    n_points = 4
    mat, position, volume = BBMaterial(), rand(3, n_points), rand(n_points)
    body = Body(mat, position, volume)

    # test body creation
    @test body.n_points == n_points
    @test body.position == position
    @test body.volume == volume
    @test body.point_sets == Dict{Symbol,Vector{Int}}(:all_points => 1:n_points)

    # add point set
    point_set!(body, :a, 1:2)
    point_set!(body, :b, 3:4)
    @test body.point_sets == Dict(:all_points => 1:n_points, :a => 1:2, :b => 3:4)

    # velocity bc 1
    forcedensity_bc!(t -> 1, body, :a, 1)
    @test length(body.single_dim_bcs) == 1
    bc1 = body.single_dim_bcs[1]
    for t in [-1, 0, 1, Inf, NaN]
        @test bc1.fun(t) == 1
    end
    @test bc1.field === :b_ext
    @test bc1.point_set === :a
    @test bc1.dim == 0x01

    # velocity bc 2
    forcedensity_bc!(t -> 2, body, :a, 2)
    @test length(body.single_dim_bcs) == 2
    bc2 = body.single_dim_bcs[2]
    for t in [-1, 0, 1, Inf, NaN]
        @test bc2.fun(t) == 2
    end
    @test bc2.field === :b_ext
    @test bc2.point_set === :a
    @test bc2.dim == 0x02

    # velocity bc 3
    forcedensity_bc!(t -> 3, body, :a, 3)
    @test length(body.single_dim_bcs) == 3
    bc3 = body.single_dim_bcs[3]
    for t in [-1, 0, 1, Inf, NaN]
        @test bc3.fun(t) == 3
    end
    @test bc3.field === :b_ext
    @test bc3.point_set === :a
    @test bc3.dim == 0x03

    # velocity bc 4
    forcedensity_bc!(t -> 4, body, :b, :x)
    @test length(body.single_dim_bcs) == 4
    bc4 = body.single_dim_bcs[4]
    for t in [-1, 0, 1, Inf, NaN]
        @test bc4.fun(t) == 4
    end
    @test bc4.field === :b_ext
    @test bc4.point_set === :b
    @test bc4.dim == 0x01

    # velocity bc 5
    forcedensity_bc!(t -> 5, body, :b, :y)
    @test length(body.single_dim_bcs) == 5
    bc5 = body.single_dim_bcs[5]
    for t in [-1, 0, 1, Inf, NaN]
        @test bc5.fun(t) == 5
    end
    @test bc5.field === :b_ext
    @test bc5.point_set === :b
    @test bc5.dim == 0x02

    # velocity bc 6
    forcedensity_bc!(t -> 6, body, :b, :z)
    @test length(body.single_dim_bcs) == 6
    bc6 = body.single_dim_bcs[6]
    for t in [-1, 0, 1, Inf, NaN]
        @test bc6.fun(t) == 6
    end
    @test bc6.field === :b_ext
    @test bc6.point_set === :b
    @test bc6.dim == 0x03
end

@testitem "forcedensity_bc! position dependend" begin
    # setup
    n_points = 4
    mat, position, volume = BBMaterial(), rand(3, n_points), rand(n_points)
    body = Body(mat, position, volume)

    # test body creation
    @test body.n_points == n_points
    @test body.position == position
    @test body.volume == volume
    @test body.point_sets == Dict{Symbol,Vector{Int}}(:all_points => 1:n_points)

    # add point set
    point_set!(body, :a, 1:2)
    point_set!(body, :b, 3:4)
    @test body.point_sets == Dict(:all_points => 1:n_points, :a => 1:2, :b => 3:4)

    # velocity bc 1
    forcedensity_bc!((p,t) -> 1, body, :a, 1)
    @test length(body.posdep_single_dim_bcs) == 1
    bc1 = body.posdep_single_dim_bcs[1]
    _p = [0, 0, 0]
    for t in [-1, 0, 1, Inf, NaN]
        @test bc1(_p, t) == 1
    end
    @test bc1.field === :b_ext
    @test bc1.point_set === :a
    @test bc1.dim == 0x01

    # velocity bc 2
    forcedensity_bc!((p, t) -> p[1] + p[2] + p[3] + t, body, :a, 2)
    @test length(body.posdep_single_dim_bcs) == 2
    bc2 = body.posdep_single_dim_bcs[2]
    _p = [0, 0, 0]
    for t in [-1, 0, 1, Inf, NaN]
        @test bc1(_p, t) == 1
    end
    @test bc2.field === :b_ext
    @test bc2.point_set === :a
    @test bc2.dim == 0x02
    @test bc2([1,2,3], 1.0) == 7.0
end

@testitem "show Body" begin
    # setup
    n_points = 10
    mat, position, volume = BBMaterial(), rand(3, n_points), rand(n_points)
    body = Body(mat, position, volume)

    io = IOBuffer()

    show(IOContext(io, :compact=>true), MIME("text/plain"), body)
    msg = String(take!(io))
    @test msg == "10-point Body{BBMaterial{NoCorrection}}"

    show(IOContext(io, :compact=>false), MIME("text/plain"), body)
    msg = String(take!(io))
    @test contains(msg, "1 point set(s):")
    @test contains(msg, "10-point set `all_points`")

    point_set!(body, :a, 1:2)

    show(IOContext(io, :compact=>true), MIME("text/plain"), body)
    msg = String(take!(io))
    @test msg == "10-point Body{BBMaterial{NoCorrection}}"

    show(IOContext(io, :compact=>false), MIME("text/plain"), body)
    msg = String(take!(io))
    @test contains(msg, "2 point set(s):")
    @test contains(msg, "2-point set `a`")
    @test contains(msg, "10-point set `all_points`")

    material!(body, horizon=1, rho=1, E=1, Gc=1)

    show(IOContext(io, :compact=>true), MIME("text/plain"), body)
    msg = String(take!(io))
    @test msg == "10-point Body{BBMaterial{NoCorrection}}"

    show(IOContext(io, :compact=>false), MIME("text/plain"), body)
    msg = String(take!(io))
    @test contains(msg, "1 point parameter(s):")
    @test contains(msg, "δ=1.0, E=1.0, nu=0.25, rho=1.0, Gc=1.0")

    material!(body, :a, horizon=2, rho=2, E=2, Gc=2)

    show(IOContext(io, :compact=>true), MIME("text/plain"), body)
    msg = String(take!(io))
    @test msg == "10-point Body{BBMaterial{NoCorrection}}"

    show(IOContext(io, :compact=>false), MIME("text/plain"), body)
    msg = String(take!(io))
    @test !contains(msg, "1 point parameter(s):")
    @test contains(msg, "2 point parameter(s):")
    @test contains(msg, "δ=1.0, E=1.0, nu=0.25, rho=1.0, Gc=1.0")
    @test contains(msg, "δ=2.0, E=2.0, nu=0.25, rho=2.0, Gc=2.0")

    velocity_ic!(body, :a, :z, 1.0)
    velocity_ic!(p -> p[1] * 2.0, body, :a, :y)
    velocity_bc!(t -> t, body, :a, 1)
    forcedensity_bc!((p, t) -> p[1] + p[2] + p[3] + t, body, :a, 2)

    show(IOContext(io, :compact=>true), MIME("text/plain"), body)
    msg = String(take!(io))
    @test msg == "10-point Body{BBMaterial{NoCorrection}}"

    show(IOContext(io, :compact=>false), MIME("text/plain"), body)
    msg = String(take!(io))
    @test contains(msg, "2 point set(s):")
    @test contains(msg, "2-point set `a`")
    @test contains(msg, "10-point set `all_points`")
    @test contains(msg, "2 boundary condition(s):")
    @test contains(msg, "BC on velocity: point_set=a, dim=1")
    @test contains(msg, "Pos.-dep. BC on force density: point_set=a, dim=2")
    @test contains(msg, "2 initial condition(s):")
    @test contains(msg, "IC on velocity: point_set=a, dim=3")
    @test contains(msg, "Pos.-dep. IC on velocity: point_set=a, dim=2")

    point_set!(body, :b, 3:4)

    precrack!(body, :a, :b)

    show(IOContext(io, :compact=>true), MIME("text/plain"), body)
    msg = String(take!(io))
    @test msg == "10-point Body{BBMaterial{NoCorrection}}"

    show(IOContext(io, :compact=>false), MIME("text/plain"), body)
    msg = String(take!(io))
    @test contains(msg, "1 predefined crack(s)")

    failure_permit!(body, :a, false)

    show(IOContext(io, :compact=>true), MIME("text/plain"), body)
    msg = String(take!(io))
    @test msg == "10-point Body{BBMaterial{NoCorrection}}"

    show(IOContext(io, :compact=>false), MIME("text/plain"), body)
    msg = String(take!(io))
    @test contains(msg, "2 points with no failure permission")

    Peridynamics.change_name!(body, :testbody)

    show(IOContext(io, :compact=>true), MIME("text/plain"), body)
    msg = String(take!(io))
    @test msg == "10-point Body{BBMaterial{NoCorrection}} with name `testbody`"

    show(IOContext(io, :compact=>false), MIME("text/plain"), body)
    msg = String(take!(io))
    @test contains(msg, "with name `testbody`")
end

@testitem "log_msg_body" begin
    # setup
    n_points = 10
    position, volume = uniform_box(1, 1, 1, 0.5)
    body = Body(BBMaterial(), position, volume)
    point_set!(body, :a, 1:2)
    material!(body, horizon=1, rho=1, E=1, Gc=1)
    material!(body, :a, horizon=2, rho=2, E=2, Gc=2)
    velocity_ic!(body, :a, :z, 1.0)
    velocity_ic!(p -> p[1] * 2.0, body, :a, :y)
    velocity_bc!(t -> t, body, :a, 1)
    forcedensity_bc!((p, t) -> p[1] + p[2] + p[3] + t, body, :a, 2)
    point_set!(body, :b, 3:4)
    precrack!(body, :a, :b)
    failure_permit!(body, :a, false)
    Peridynamics.change_name!(body, :testbody)

    msg = Peridynamics.log_msg_body(body)

    @test msg == """
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
            material type ............. Peridynamics.BBMaterial{Peridynamics.NoCorrection}
            MATERIAL PROPERTIES #1
              horizon .................................................................. 1
              density .................................................................. 1
              Young's modulus .......................................................... 1
              Poisson's ratio ....................................................... 0.25
              shear modulus .......................................................... 0.4
              bulk modulus ..................................................... 0.6666667
            MATERIAL PROPERTIES #2
              horizon .................................................................. 2
              density .................................................................. 2
              Young's modulus .......................................................... 2
              Poisson's ratio ....................................................... 0.25
              shear modulus .......................................................... 0.8
              bulk modulus ...................................................... 1.333333
        """
end

@testitem "Body from inp file" begin
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
