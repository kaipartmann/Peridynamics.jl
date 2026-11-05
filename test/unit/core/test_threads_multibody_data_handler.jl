# `ThreadsMultibodyDataHandler`: one threads data handler per body of a `MultibodySetup`, the
# contacts and the position caches the contact search needs (`src/core/threads_multibody_data_handler.jl`).

@testitem "ThreadsMultibodyDataHandler: structure and accessors" setup=[Fixtures] begin
    b1 = Fixtures.line10(BBMaterial(); Gc=1)
    velocity_ic!(b1, :all_points, :x, 1.0)
    b2 = Fixtures.line10(OSBMaterial(); Gc=1)
    b2.position[2, :] .+= 1.2 # next to the first body, within the contact radius
    ms = MultibodySetup(:first => b1, :second => b2)
    contact!(ms, :first, :second; radius=1.5)
    solver = VelocityVerlet(steps=1)
    dh = Peridynamics.threads_data_handler(ms, solver, 2)
    @test dh isa Peridynamics.ThreadsMultibodyDataHandler
    @test dh.n_bodies == 2
    @test dh.body_names == [:first, :second]
    @test dh.body_idxs == Dict(:first => 1, :second => 2)
    @test length(dh.srf_contacts) == 2 # one per direction of the contact pair
    @test Peridynamics.get_body_dh(dh, :first) === Peridynamics.get_body_dh(dh, 1)
    @test Peridynamics.get_body_dh(dh, 1) isa Peridynamics.ThreadsBodyDataHandler
    @test Peridynamics.get_body_dh(dh, 1).n_chunks == 2
    @test Peridynamics.get_body_name(dh, 2) == "second"
    @test collect(Peridynamics.each_body_idx(dh)) == [1, 2]
    @test collect(Peridynamics.each_body_name(dh)) == [:first, :second]
    @test length(collect(Peridynamics.each_body_dh(dh))) == 2
    # the caches hold the positions and volumes of the whole bodies
    @test size(dh.position_caches[1]) == (3, 10) && size(dh.position_caches[2]) == (3, 10)
    @test dh.volume_caches[1] == ones(10)
    # after initialization the caches are filled with the current positions
    Peridynamics.init_time_solver!(solver, dh)
    Peridynamics.initialize!(dh, solver)
    Peridynamics.update_caches!(dh)
    @test dh.position_caches[1] ≈ b1.position
    @test dh.position_caches[2] ≈ b2.position
    # a force calculation with contact leaves finite forces on both bodies
    Peridynamics.calc_force_density!(dh, solver.Δt, solver.Δt)
    Peridynamics.update_caches!(dh)
    Peridynamics.calc_contact_force_densities!(dh)
    for body_dh in Peridynamics.each_body_dh(dh)
        for chunk in body_dh.chunks
            @test all(isfinite, chunk.storage.b_int)
        end
    end
    # unknown body names are an error
    @test_throws KeyError Peridynamics.get_body_dh(dh, :third)
end
