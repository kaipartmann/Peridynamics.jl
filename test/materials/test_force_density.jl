# Invariants of the internal force density that every material must satisfy on a deformed body,
# without running a time loop. These items loop `Fixtures.MATERIALS`; a new material is covered
# by adding it there.

@testitem "force density: zero at the reference configuration" setup=[Fixtures] begin
    # zero up to rounding: the current and the reference bond length are computed in different
    # orders of operations, and the difference of order `eps()` is multiplied by the stiffness,
    # so the forces are compared with the force density scale `E / Δx` of a unit strain
    for (name, mat) in Fixtures.MATERIALS
        @testset "$name" begin
            n = Fixtures.default_n(mat)
            body = Fixtures.cube(mat; n)
            solver = VelocityVerlet(steps=1)
            dh = Fixtures.handler(body, solver)
            chunk = dh.chunks[1]
            Peridynamics.calc_force_density!(dh, solver.Δt, solver.Δt)
            scale = body.point_params[1].E * n
            @test maximum(abs, chunk.storage.b_int) < 1e-12 * scale
        end
    end
end

@testitem "force density: conserves linear momentum" setup=[Fixtures] begin
    # The internal forces are interaction forces, so their sum over the body vanishes for any
    # deformation — including the inhomogeneous one of `material_fixture`, which activates the
    # zero-energy mode stabilizations. Checked relative to the total force magnitude.
    for (name, mat) in Fixtures.MATERIALS
        @testset "$name" begin
            fixture = Fixtures.material_fixture(mat)
            Fixtures.force_density!(fixture)
            (; b_int) = fixture.chunk.storage
            (; volume) = fixture.chunk.system
            n_loc = Peridynamics.get_n_loc_points(fixture.chunk.system)
            total = sum(b_int[:, i] * volume[i] for i in 1:n_loc)
            scale = sum(abs.(b_int[:, i]) * volume[i] for i in 1:n_loc)
            @test !iszero(scale) # the deformation produced forces at all
            @test all(abs.(total) .< 1e-10 * scale)
        end
    end
end

@testitem "force density: mirror symmetric for a mirror symmetric deformation" setup=[Fixtures] begin
    # A body that is mirror-symmetric in x, y and z under a diagonal deformation gradient has a
    # mirror-symmetric force density: the component normal to a mirror plane flips its sign, the
    # others are equal. This is the static core of what the symmetry simulations in
    # `test/simulation/` check with a time loop.
    using Peridynamics.StaticArrays
    F = @SMatrix [1.01 0.0 0.0; 0.0 0.995 0.0; 0.0 0.0 1.02]
    for (name, mat) in Fixtures.MATERIALS
        @testset "$name" begin
            # `cube` from `uniform_box` is symmetric about the origin
            n = Fixtures.default_n(mat)
            if mat isa CKIMaterial
                # explicit interaction parameters, otherwise no triplets exist
                body = @test_logs (:warn, r"specified manually") Fixtures.cube(mat; n,
                                                                               Fixtures.cki_kwargs()...)
            else
                body = Fixtures.cube(mat; n)
            end
            solver = VelocityVerlet(steps=1)
            dh = Fixtures.handler(body, solver)
            chunk = dh.chunks[1]
            Fixtures.apply_deformation!(chunk, F, solver.Δt)
            Peridynamics.calc_force_density!(dh, solver.Δt, solver.Δt)
            X = chunk.system.position
            b = chunk.storage.b_int
            n_loc = Peridynamics.get_n_loc_points(chunk.system)
            @test !iszero(b)
            # for each point find its mirror image in direction d and compare the forces
            for d in 1:3
                mirror = [i for i in 1:n_loc] # placeholder, filled below
                for i in 1:n_loc
                    Xm = X[:, i] .* (d == 1 ? [-1, 1, 1] : d == 2 ? [1, -1, 1] : [1, 1, -1])
                    mirror[i] = findfirst(j -> isapprox(X[:, j], Xm; atol=1e-12), 1:n_loc)
                end
                for e in 1:3
                    sign = d == e ? -1 : 1
                    @test all(isapprox(b[e, i], sign * b[e, mirror[i]]; atol=1e-10 * maximum(abs, b))
                              for i in 1:n_loc)
                end
            end
        end
    end
end

@testitem "force density: bonds break beyond the critical stretch" setup=[Fixtures] begin
    # A corner point of a 3³ cube is pulled far away; with a tiny energy release rate the
    # critical stretch is far below that, so all its bonds break in the first force
    # calculation: its damage is 1, and every former neighbor lost exactly one bond.
    for (name, mat) in Fixtures.MATERIALS
        @testset "$name" begin
            pos, vol = uniform_box(1, 1, 1, 1 / 3)
            body = Body(mat, pos, vol)
            kwargs = (; horizon=3.015 / 3, rho=7850.0, E=210e9, nu=0.25, Gc=1e-6)
            if mat isa CKIMaterial
                @test_logs (:warn, r"specified manually") material!(body; kwargs...,
                                                                      Fixtures.cki_kwargs()...)
            else
                material!(body; kwargs...)
            end
            @test all(body.fail_permit)
            solver = VelocityVerlet(steps=1)
            dh = Fixtures.handler(body, solver)
            chunk = dh.chunks[1]
            (; storage, system) = chunk
            pulled = 1
            storage.position[:, pulled] .+= 2.0 # far outside the cube
            Fixtures.force_density!(chunk, solver.Δt, solver.Δt)
            Peridynamics.calc_damage!(chunk)
            # the bonds (or one-neighbor interactions) of every point, as neighbor lists
            interactions = system isa Peridynamics.InteractionSystem ? system.one_nis : system.bonds
            n_active = system isa Peridynamics.InteractionSystem ? storage.n_active_one_nis :
                       storage.n_active_bonds
            @test storage.damage[pulled] ≈ 1.0
            @test n_active[pulled] == 0
            for i in 2:n_points(body)
                neighbors = [interactions[b].neighbor for b in Peridynamics.each_bond_idx(system, i)]
                lost = count(==(pulled), neighbors)
                @test n_active[i] == length(neighbors) - lost
                @test storage.damage[i] ≈ lost / length(neighbors)
            end
        end
    end
end
