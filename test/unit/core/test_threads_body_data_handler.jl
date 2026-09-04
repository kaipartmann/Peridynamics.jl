@testitem "ThreadsBodyDataHandler BBMaterial VelocityVerlet" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    point_set!(body, :a, 1:2)
    point_set!(body, :b, 3:4)
    velocity_ic!(body, :a, :x, 1.0)
    velocity_bc!(t->t, body, :a, :x)
    forcedensity_bc!(t->t, body, :a, :x)
    precrack!(body, :a, :b)
    ts = VelocityVerlet(steps=10)
    dh = Peridynamics.threads_data_handler(body, ts, 2)

    @test dh.n_chunks == 2
    @test length(dh.chunks) == 2
    @test length(dh.lth_exs[1]) == 1
    @test length(dh.lth_exs[2]) == 1
    @test length(dh.htl_exs[1]) == 1
    @test length(dh.htl_exs[2]) == 1

    io = IOBuffer()

    show(IOContext(io, :compact=>true), MIME("text/plain"), dh)
    msg = String(take!(io))
    @test contains(msg, "DataHandler()")

    show(IOContext(io, :compact=>false), MIME("text/plain"), dh)
    msg = String(take!(io))
    @test contains(msg, "DataHandler()")
end

@testitem "threads_data_handler: max_n_chunks clamps the decomposition, first_chunk" begin
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection, CriticalStretch

    # a material whose system exists once per body, e.g. one that transforms the whole body
    # at once, so it may never be decomposed
    struct OneChunkMat <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::CriticalStretch
    end
    OneChunkMat() = OneChunkMat(CriticalStretch())
    Peridynamics.@params OneChunkMat struct OneChunkParams
        @inherit StandardParameters
    end
    Peridynamics.@storage OneChunkMat struct OneChunkStorage
        @inherit Peridynamics.VelocityVerletFields
        @inherit Peridynamics.BondFracFields
    end
    function Peridynamics.force_density_point!(::OneChunkStorage, system, ::OneChunkMat,
                                               params, t, Δt, i)
        return nothing
    end
    Peridynamics.max_n_chunks(::OneChunkMat) = 1

    # the default is no limit at all
    @test Peridynamics.max_n_chunks(BBMaterial()) == typemax(Int)
    @test Peridynamics.max_n_chunks(OneChunkMat()) == 1

    position = zeros(3, 10)
    position[1, :] .= 0.0:9.0
    solver = VelocityVerlet(steps=1)

    # a material without a limit gets the chunks it asked for
    body = Body(BBMaterial(), position, ones(10))
    material!(body; horizon=1.5, rho=1, E=1, nu=0.25, Gc=1)
    dh = Peridynamics.threads_data_handler(body, solver, 4)
    @test dh.n_chunks == 4
    @test Peridynamics.first_chunk(dh) === dh.chunks[1]
    @test Peridynamics.get_body_name(dh) === dh.chunks[1].body_name

    # a material with a limit is clamped to it, and the one chunk owns every point
    limited = Body(OneChunkMat(), position, ones(10))
    material!(limited; horizon=1.5, rho=1, E=1, nu=0.25, Gc=1)
    dh_limited = Peridynamics.threads_data_handler(limited, solver, 4)
    @test dh_limited.n_chunks == 1
    @test Peridynamics.get_n_loc_points(Peridynamics.first_chunk(dh_limited)) == 10

    # the number of chunks is still never larger than the number of points
    small = Body(BBMaterial(), position[:, 1:2], ones(2))
    material!(small; horizon=1.5, rho=1, E=1, nu=0.25, Gc=1)
    @test Peridynamics.threads_data_handler(small, solver, 8).n_chunks == 2

    # an MPI run cannot clamp, because the ranks come from the outside
    @test Peridynamics.check_max_n_chunks(BBMaterial(), 8) === nothing
    @test Peridynamics.check_max_n_chunks(OneChunkMat(), 1) === nothing
    msg = try
        Peridynamics.check_max_n_chunks(OneChunkMat(), 4)
    catch e
        sprint(showerror, e)
    end
    @test contains(msg, "OneChunkMat")
    @test contains(msg, "at most 1 chunk")
    @test contains(msg, "4 MPI ranks")
end
