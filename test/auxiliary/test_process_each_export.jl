@testitem "process_each_export" begin
    root = joinpath(@__DIR__, "temp_root_proc_each_export")
    isdir(root) && rm(root; recursive=true, force=true)
    root_post_threads = joinpath(root, "post_threads")
    mkpath(root_post_threads)
    root_post_serial = joinpath(root, "post_serial")
    mkpath(root_post_serial)
    root_post_mpi = joinpath(root, "post_mpi")
    mkpath(root_post_mpi)

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

    process_each_export(job) do result, file_id
        filename = string("max_displacement_", file_id, ".txt")
        open(joinpath(root_post_threads, filename), "w+") do io
            maxdisp = maximum(result[:displacement][1,:])
            msg = string("maximum displacement x: ", maxdisp)
            write(io, msg)
        end
        return nothing
    end
    file_1_threads = joinpath(root_post_threads, "max_displacement_1.txt")
    @test isfile(file_1_threads)
    @test contains(read(file_1_threads, String), "maximum displacement x: 0.0")
    file_2_threads = joinpath(root_post_threads, "max_displacement_2.txt")
    @test isfile(file_2_threads)
    @test contains(read(file_2_threads, String), "maximum displacement x: 0.0")
    file_3_threads = joinpath(root_post_threads, "max_displacement_3.txt")
    @test isfile(file_3_threads)
    @test contains(read(file_3_threads, String), "maximum displacement x: 2.4")

    process_each_export(job; serial=true) do result, file_id
        filename = string("max_displacement_", file_id, ".txt")
        open(joinpath(root_post_serial, filename), "w+") do io
            msg = string("maximum displacement x: ", maximum(result[:displacement][1,:]))
            write(io, msg)
        end
        return nothing
    end
    file_1_serial = joinpath(root_post_serial, "max_displacement_1.txt")
    @test isfile(file_1_serial)
    @test contains(read(file_1_serial, String), "maximum displacement x: 0.0")
    file_2_serial = joinpath(root_post_serial, "max_displacement_2.txt")
    @test isfile(file_2_serial)
    @test contains(read(file_2_serial, String), "maximum displacement x: 0.0")
    file_3_serial = joinpath(root_post_serial, "max_displacement_3.txt")
    @test isfile(file_3_serial)
    @test contains(read(file_3_serial, String), "maximum displacement x: 2.4")

    mpi_cmd = """
    using Peridynamics
    files = "$(joinpath(root, "vtk"))"
    process_each_export(files) do result, file_id
        filename = string("max_displacement_", file_id, ".txt")
        open(joinpath("$(root_post_mpi)", filename), "w+") do io
            maxdisp = maximum(result[:displacement][1,:])
            msg = string("maximum displacement x: ", maxdisp)
            write(io, msg)
        end
        return nothing
    end
    """
    run(`$(Peridynamics.MPI.mpiexec()) -n 3 $(Base.julia_cmd()) --project -e $(mpi_cmd)`)
    file_1_mpi = joinpath(root_post_mpi, "max_displacement_1.txt")
    @test isfile(file_1_mpi)
    @test contains(read(file_1_mpi, String), "maximum displacement x: 0.0")
    file_2_mpi = joinpath(root_post_mpi, "max_displacement_2.txt")
    @test isfile(file_2_mpi)
    @test contains(read(file_2_mpi, String), "maximum displacement x: 0.0")
    file_3_mpi = joinpath(root_post_mpi, "max_displacement_3.txt")
    @test isfile(file_3_mpi)
    @test contains(read(file_3_mpi, String), "maximum displacement x: 2.4")

    rm(root; recursive=true, force=true)
end
