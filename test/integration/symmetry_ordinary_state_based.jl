@testitem "symmetry OSBMaterial" begin
    # simulation
    Δx = 0.2
    width = 1
    grid⁺ = Δx/2:Δx:width/2
    grid⁻ = -reverse(grid⁺)
    grid = [grid⁻; grid⁺]
    pos = hcat(([x;y;z] for x in grid for y in grid for z in grid)...)
    n_points = size(pos, 2)
    vol = fill(Δx^3, n_points)
    body = Body(OSBMaterial(), pos, vol)
    failure_permit!(body, false)
    material!(body, horizon=3.015Δx, rho=7850, E=210e9, nu=0.25, Gc=1)
    point_set!(z -> z > width/2 - 0.6Δx, body, :set_a)
    point_set!(z -> z < -width/2 + 0.6Δx, body, :set_b)
    velocity_bc!(t -> 100, body, :set_a, :z)
    velocity_bc!(t -> -100, body, :set_b, :z)
    vv = VelocityVerlet(steps=100)
    temp_root = joinpath(@__DIR__, "temp_root_symmetry_test")
    rm(temp_root; recursive=true, force=true)
    job = Job(body, vv; path=temp_root, freq=10)
    dh = Peridynamics.submit_threads(job, 1)

    # check if the correct files were exported
    @test length(filter(x->endswith(x,".pvtu"), readdir(joinpath(temp_root,"vtk")))) == 11

    # find the points that should have symmetrical positions
    xp = findall(pos[1,:] .> 0)
    xm = findall(pos[1,:] .< 0)
    yp = findall(pos[2,:] .> 0)
    ym = findall(pos[2,:] .< 0)
    zp = findall(pos[3,:] .> 0)
    zm = findall(pos[3,:] .< 0)

    # check the symmetry
    end_position = first(dh.chunks).storage.position
    @test !(end_position ≈ pos)
    @test sort(end_position[1,xp]) ≈ sort(-end_position[1,xm])
    @test sort(end_position[2,xp]) ≈ sort(end_position[2,xm])
    @test sort(end_position[3,xp]) ≈ sort(end_position[3,xm])
    @test sort(end_position[1,yp]) ≈ sort(end_position[1,ym])
    @test sort(end_position[2,yp]) ≈ sort(-end_position[2,ym])
    @test sort(end_position[3,yp]) ≈ sort(end_position[3,ym])
    @test sort(end_position[1,zp]) ≈ sort(end_position[1,zm])
    @test sort(end_position[2,zp]) ≈ sort(end_position[2,zm])
    @test sort(end_position[3,zp]) ≈ sort(-end_position[3,zm])

    # check if points have moved and are not zero
    end_displacement = first(dh.chunks).storage.displacement
    for i in axes(end_displacement, 2)
        for d in axes(end_displacement, 1)
            @test abs(end_displacement[d, i]) > 0
        end
    end

    # delete all simulation files
    rm(temp_root; recursive=true, force=true)
end
