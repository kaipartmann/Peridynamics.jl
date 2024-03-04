@testitem "process_each_export" begin
    root = joinpath(@__DIR__, "temp_root_proc_each_export")
    root_post = joinpath(root, "post")
    isdir(root) && rm(root; recursive=true, force=true)

    l, Δx, δ, a = 1.0, 1 / 4, 3.015 / 4, 0.5
    pos, vol = uniform_box(l, l, l, Δx)
    b1 = Body(BBMaterial(), pos, vol)
    material!(b1; horizon=3.015Δx, E=2.1e5, rho=8e-6, Gc=2.7)
    point_set!(y -> y > l / 2 - Δx, b1, :set_top)
    point_set!(y -> y < -l / 2 + Δx, b1, :set_bottom)
    velocity_bc!(t -> 30, b1, :set_top, :y)
    velocity_bc!(t -> -30, b1, :set_bottom, :y)
    vv = VelocityVerlet(steps=2)
    job = Job(b1, vv; path=root, freq=1)

    @test_throws ArgumentError process_each_export((r, id) -> nothing, job)

    submit(job)

    file_1 = joinpath(root_post, "max_displacement_1.txt")
    file_2 = joinpath(root_post, "max_displacement_2.txt")
    file_3 = joinpath(root_post, "max_displacement_3.txt")

    mkpath(root_post)
    process_each_export(job) do result, file_id
        filename = string("max_displacement_", file_id, ".txt")
        open(joinpath(root_post, filename), "w+") do io
            maxdisp = maximum(result[:displacement][1,:])
            msg = string("maximum displacement x: ", maxdisp)
            write(io, msg)
        end
        return nothing
    end
    @test isfile(file_1)
    @test contains(read(file_1, String), "maximum displacement x: 0.0")
    @test isfile(file_2)
    @test contains(read(file_2, String), "maximum displacement x: 0.0")
    @test isfile(file_3)
    @test contains(read(file_3, String), "maximum displacement x: 2.4")

    rm(root_post; recursive=true, force=true)
    mkpath(root_post)
    process_each_export(job; serial=true) do result, file_id
        filename = string("max_displacement_", file_id, ".txt")
        open(joinpath(root_post, filename), "w+") do io
            msg = string("maximum displacement x: ", maximum(result[:displacement][1,:]))
            write(io, msg)
        end
        return nothing
    end
    @test isfile(file_1)
    @test contains(read(file_1, String), "maximum displacement x: 0.0")
    @test isfile(file_2)
    @test contains(read(file_2, String), "maximum displacement x: 0.0")
    @test isfile(file_3)
    @test contains(read(file_3, String), "maximum displacement x: 2.4")

    rm(root_post; recursive=true, force=true)
    mkpath(root_post)
    mpi_cmd = """
    using Peridynamics
    files = "$(joinpath(root, "vtk"))"
    process_each_export(files) do result, file_id
        filename = string("max_displacement_", file_id, ".txt")
        open(joinpath("$(root_post)", filename), "w+") do io
            maxdisp = maximum(result[:displacement][1,:])
            msg = string("maximum displacement x: ", maxdisp)
            write(io, msg)
        end
        return nothing
    end
    """
    run(`$(Peridynamics.MPI.mpiexec()) -n 3 $(Base.julia_cmd()) --project -e $(mpi_cmd)`)
    @test isfile(file_1)
    @test contains(read(file_1, String), "maximum displacement x: 0.0")
    @test isfile(file_2)
    @test contains(read(file_2, String), "maximum displacement x: 0.0")
    @test isfile(file_3)
    @test contains(read(file_3, String), "maximum displacement x: 2.4")

    rm(root; recursive=true, force=true)
end
