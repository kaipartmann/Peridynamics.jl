@testitem "Contact example with 4 points" begin
    pos1 = [
        0.0 1.0
        0.0 0.0
        0.0 0.0
    ]
    pos2 = [
        1.5 2.5
        0.0 0.0
        0.0 0.0
    ]
    point_spacing = 1.0
    δ = 1.5 * point_spacing
    E = 210e9
    rho = 7850.0
    n_points = 2
    vol = fill(point_spacing^3, n_points)
    body1 = Body(BBMaterial(), pos1, vol)
    material!(body1, horizon=δ, rho=rho, E=E, Gc=1.0)
    body2 = Body(BBMaterial(), pos2, vol)
    material!(body2, horizon=δ, rho=rho, E=E, Gc=1.0)
    ms = MultibodySetup(:body1 => body1, :body2 => body2)
    contact!(ms, :body1, :body2; radius=point_spacing)
    job = Job(ms, VelocityVerlet(steps=1))
    dh = Peridynamics.submit_threads(job, 1)

    @test ms.srf_contacts[1].penalty_factor == 1e12
    @test ms.srf_contacts[2].penalty_factor == 1e12
    b_int_contact = 0.5 * 9 * 1e12 / (π * δ^4) * (1-0.5) / 0.5
    s1 = dh.body_dhs[1].chunks[1].storage
    s2 = dh.body_dhs[2].chunks[1].storage
    @test s1.b_int[1,2] ≈ -b_int_contact
    @test s2.b_int[1,1] ≈ b_int_contact
    acc_contact = b_int_contact / rho
    @test s1.acceleration[1,2] ≈ -acc_contact
    @test s2.acceleration[1,1] ≈ acc_contact
    vel_contact = acc_contact * job.time_solver.Δt * 0.5
    @test s1.velocity[1,2] ≈ -vel_contact
    @test s2.velocity[1,1] ≈ vel_contact
end
