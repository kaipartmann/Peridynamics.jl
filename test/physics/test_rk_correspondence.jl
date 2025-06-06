@testitem "damage changed flag" begin
    pos, vol = uniform_box(1, 1, 1, 0.4)
    body = Body(RKCMaterial(), pos, vol)
    material!(body; horizon=1.5, rho=1, E=210e9, nu=0.25, Gc=1.0)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; n_neighbors) = chunk.system
    (; update_gradients, damage, n_active_bonds, bond_active, gradient_weight) = chunk.storage

    @test n_neighbors == fill(7, 8)
    @test n_active_bonds == fill(7, 8)
    @test bond_active == fill(true, 56)
    @test gradient_weight[:,1] ≈ [1.1548446667165417, -0.3316495750857463, -0.3316495750857462]
    @test gradient_weight[:,2] ≈ [-0.33164957508574605, 1.1548446667165417, -0.3316495750857461]
    @test update_gradients == fill(false, 8)

    Peridynamics.calc_weights_and_defgrad!(chunk, 0.0, 0.0)

    # everything should be the same, the gradients are initialized and damage did not change
    @test n_active_bonds == fill(7, 8)
    @test bond_active == fill(true, 56)
    @test gradient_weight[:,1] ≈ [1.1548446667165417, -0.3316495750857463, -0.3316495750857462]
    @test gradient_weight[:,2] ≈ [-0.33164957508574605, 1.1548446667165417, -0.3316495750857461]
    @test update_gradients == fill(false, 8)

    # bond 1 failed somehow
    bond_active[1] = false
    Peridynamics.calc_weights_and_defgrad!(chunk, 0.0, 0.0)

    # now the changes should be reflected in the chunk
    @test n_active_bonds == [6, 7, 7, 7, 7, 7, 7, 7]
    @test bond_active == [false; fill(true, 55)]
    @test iszero(gradient_weight[:,1])
    @test gradient_weight[:,2] ≈ [-0.6163778391975847, 1.2366132461014017, -0.24988099570088682]
    @test update_gradients == fill(false, 8)
end


# @testitem "gradient weights" begin
#     pos, vol = uniform_box(1, 1, 1, 0.4)
#     body1 = Body(RKCMaterial(), pos, vol)
#     material!(body1; horizon=2, rho=1, E=210e9, nu=0.25, Gc=1.0)

#     dh1 = Peridynamics.threads_data_handler(body1, VelocityVerlet(steps=1), 1)
#     # (; system, mat, storage, paramsetup) = dh.chunks[1]
#     # (; gradient_weight) = storage
#     # params = paramsetup
#     # t = 0.0
#     # Δt = 0.0

#     body2 = deepcopy(Body(RKCMaterial(), pos, vol))
#     material!(body2; horizon=2, rho=1, E=210e9, nu=0.25, Gc=1.0)

#     dh2 = Peridynamics.threads_data_handler(body2, VelocityVerlet(steps=1), 1)
#     # (; system, mat, storage, paramsetup) = dh.chunks[1]
#     # (; gradient_weight) = storage

#     @test dh1.chunks[1].storage.gradient_weight ≈ dh2.chunks[1].storage.gradient_weight

#     dispvar = 1e-8 * randn(3, 8)
#     posvar = deepcopy(dh1.chunks[1].storage.position) .+ dispvar
#     dh1.chunks[1].storage.displacement .= copy(dispvar)
#     dh1.chunks[1].storage.position .= copy(posvar)
#     dh2.chunks[1].storage.displacement .= copy(dispvar)
#     dh2.chunks[1].storage.position .= copy(posvar)

#     Peridynamics.calc_weights_and_defgrad!(dh1.chunks[1], 0.0, 0.0)
#     Peridynamics.calc_weights_and_defgrad!(dh2.chunks[1], 0.0, 0.0)

#     @test dh1.chunks[1].storage.defgrad ≈ dh2.chunks[1].storage.defgrad

#     # # @test gradient_weight[:,1] ≈ [18.04444791744596, -5.182024610714787, -5.182024610714787]

#     # # @test gradient_weight[:,1] ≈ [1.1548446667165413, -0.33164957508574583, -0.331649575085746]
#     # # # @test gradient_weight[:,2] ≈ [-0.3316495750857459, 1.154844666716541, -0.3316495750857463]
#     # # # @test gradient_weight[:,3] ≈ [0.5612621576497658, 0.5612621576497656, -0.45224360054794666]

#     # gradient_weight .= 0.0

#     # for i in 1:3
#     #     Peridynamics.rkc_weights!(storage, system, mat, params, t, Δt, i)
#     # end

#     # @test gradient_weight[:,1] ≈ [18.04444791744596, -5.182024610714787, -5.182024610714787]

#     # # @test gradient_weight[:,1] ≈ [1.1548446667165413, -0.33164957508574583, -0.331649575085746]
#     # # # @test gradient_weight[:,2] ≈ [-0.3316495750857459, 1.154844666716541, -0.3316495750857463]
#     # # # @test gradient_weight[:,3] ≈ [0.5612621576497658, 0.5612621576497656, -0.45224360054794666]


#     # gradient_weight .= 0.0

#     # for i in 1:3
#     #     Peridynamics.my_rkc_weights!(storage, system, mat, params, t, Δt, i)
#     # end

#     # @test gradient_weight[:,1] ≈ [18.04444791744596, -5.182024610714787, -5.182024610714787]

#     # # @test gradient_weight[:,1] ≈ [1.1548446667165413, -0.33164957508574583, -0.331649575085746]
#     # # # @test gradient_weight[:,2] ≈ [-0.3316495750857459, 1.154844666716541, -0.3316495750857463]
#     # # # @test gradient_weight[:,3] ≈ [0.5612621576497658, 0.5612621576497656, -0.45224360054794666]

# end
