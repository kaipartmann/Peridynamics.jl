# Under a homogeneous deformation `x = F X`, the deformation gradient of every correspondence
# formulation is exact: the nonlocal average of an affine map is the map itself. This holds for
# every point including the surface points with truncated families.

@testitem "deformation gradient: exact for homogeneous deformations" setup=[Fixtures] begin
    using Peridynamics.StaticArrays
    F_a = @SMatrix [1.1 0.01 0.0; 0.0 1.2 0.02; 0.03 0.0 1.3]
    for (name, mat) in Fixtures.CORE_MATERIALS
        Fixtures.is_correspondence(mat) || continue
        @testset "$name" begin
            body = Fixtures.cube(mat; n=5, m=3.1)
            ts = VelocityVerlet(steps=1)
            dh = Fixtures.handler(body, ts)
            chunk = dh.chunks[1]
            Fixtures.apply_deformation!(chunk, F_a, ts.Δt)
            Peridynamics.calc_force_density!(dh, ts.Δt, ts.Δt)
            (; storage, system, paramsetup) = chunk
            if mat isa Peridynamics.AbstractCorrespondenceMaterial
                # the classic formulations compute the deformation gradient on the fly
                params = Peridynamics.get_params(paramsetup, 1)
                @test all(Peridynamics.each_point_idx(chunk)) do i
                    (; F) = Peridynamics.calc_deformation_gradient!(storage, system, mat,
                                                                    params, i)
                    F ≈ F_a
                end
            elseif mat isa Peridynamics.AbstractRKCMaterial
                @test all(Peridynamics.get_tensor(storage.defgrad, i) ≈ F_a
                          for i in Peridynamics.each_point_idx(chunk))
            else # bond-associated: one deformation gradient per bond
                params = Peridynamics.get_params(paramsetup, 1)
                @test all(Peridynamics.each_point_idx(chunk)) do i
                    all(Peridynamics.each_bond_idx(system, i)) do bond_id
                        (; F) = Peridynamics.calc_deformation_gradient(storage, system, mat,
                                                                       params, i, bond_id)
                        F ≈ F_a
                    end
                end
            end
        end
    end
end
