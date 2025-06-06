@testitem "uniform displacement CRMaterial" begin
    using Peridynamics.StaticArrays
    using Peridynamics.LinearAlgebra

    # Create a uniform box with the specified dimensions and discretization
    l, w, h = 1.0, 1.0, 1.0
    ΔX = w/5
    δ = 3.1 * ΔX
    ref_position, volume = uniform_box(l, w, h, ΔX)
    mat = CRMaterial()
    body = Body(mat, ref_position, volume)
    material!(body; horizon=δ, rho=8000, E=210e9, nu=0.25)
    no_failure!(body)
    ts = VelocityVerlet(steps=1)

    # create the internal data structure for the test
    dh = Peridynamics.threads_data_handler(body, ts, 1) # just one thread
    Peridynamics.initialize!(dh, ts)
    chunk = dh.chunks[1]
    (; storage, system, paramsetup) = chunk
    (; bonds, volume) = system
    (; position) = chunk.storage

    # define a uniform deformation gradient for the whole body
    F_a = @SMatrix [1.5 0.0 0.0
                    0.0 1.5 0.0
                    0.0 0.0 1.5]
    J = det(F_a)

    # apply a uniform displacement to the body
    for i in Peridynamics.each_point_idx(chunk)
        X = Peridynamics.get_vector(position, i)
        x = F_a * X
        for d in 1:3
            position[d, i] = x[d]
        end
    end

    # run the force calculation step to compute the deformation gradient
    Peridynamics.calc_force_density!(dh, 0, 0)

    # check if the deformation gradient is equal to the expected value
    for i in Peridynamics.each_point_idx(chunk)
        (; F) = Peridynamics.calc_deformation_gradient(storage, system, mat, paramsetup, i)
        @test F ≈ F_a
    end
end

@testitem "uniform displacement RKCRMaterial" begin
    using Peridynamics.StaticArrays
    using Peridynamics.LinearAlgebra

    # Create a uniform box with the specified dimensions and discretization
    l, w, h = 1.0, 1.0, 1.0
    ΔX = w/5
    δ = 3.1 * ΔX
    ref_position, volume = uniform_box(l, w, h, ΔX)
    mat = RKCRMaterial()
    body = Body(mat, ref_position, volume)
    material!(body; horizon=δ, rho=8000, E=210e9, nu=0.25)
    no_failure!(body)
    ts = VelocityVerlet(steps=1)

    # create the internal data structure for the test
    dh = Peridynamics.threads_data_handler(body, ts, 1) # just one thread
    Peridynamics.initialize!(dh, ts)
    chunk = dh.chunks[1]
    (; storage, system, paramsetup) = chunk
    (; bonds, volume) = system
    (; position, defgrad) = chunk.storage

    # define a uniform deformation gradient for the whole body
    F_a = @SMatrix [1.5 0.0 0.0
                    0.0 1.5 0.0
                    0.0 0.0 1.5]
    J = det(F_a)

    # apply a uniform displacement to the body
    for i in Peridynamics.each_point_idx(chunk)
        X = Peridynamics.get_vector(position, i)
        x = F_a * X
        for d in 1:3
            position[d, i] = x[d]
        end
    end

    # run the force calculation step to compute the deformation gradient
    Peridynamics.calc_force_density!(dh, 0, 0)

    # check if the deformation gradient is equal to the expected value
    for i in Peridynamics.each_point_idx(chunk)
        @test Peridynamics.get_tensor(defgrad, i) ≈ F_a
    end
end
