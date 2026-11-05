# A body that is mirror-symmetric in all three directions and loaded symmetrically must stay
# symmetric, for every material and every time solver. These are the per-commit smoke runs of
# the solvers: a few steps on a small cube, no accuracy claim. Each item loops over all materials,
# so adding a material to `Fixtures.MATERIALS` covers it here automatically.

@testmodule SymmetryCase begin
    using Peridynamics, Test

    # a cube of width 1 whose points are mirror images of each other in x, y and z
    function symmetric_cube(mat, Δx)
        grid⁺ = Δx/2:Δx:0.5
        grid = [-reverse(grid⁺); grid⁺]
        pos = hcat(([x; y; z] for x in grid for y in grid for z in grid)...)
        vol = fill(Δx^3, size(pos, 2))
        body = Body(mat, pos, vol)
        # the top and bottom layer of points, loaded in ±z
        point_set!(z -> z > 0.5 - 0.6Δx, body, :set_a)
        point_set!(z -> z < -0.5 + 0.6Δx, body, :set_b)
        return body
    end

    # `material!` with the parameters every symmetry run uses, including the CKI interaction set
    function symmetric_material!(body, mat, Δx)
        horizon = horizon_ratio(mat) * Δx
        if mat isa CKIMaterial
            # explicit interaction parameters, otherwise no triplets are built (see `cki_kwargs`)
            @test_logs (:warn, r"specified manually") material!(body; horizon, rho=7850,
                                                                E=210e9, nu=0.25, C1=1e11,
                                                                C2=1e11, C3=1e11)
        else
            material!(body; horizon, rho=7850, E=210e9, nu=0.25)
        end
        return body
    end
    # same as `Fixtures.horizon_ratio`; test modules cannot depend on each other
    horizon_ratio(::Peridynamics.AbstractMaterial) = 3.015
    horizon_ratio(::Peridynamics.AbstractInteractionSystemMaterial) = 2.015

    # run `body` with `solver`, exporting every `freq` steps, and check that the final positions
    # are mirror-symmetric and that every point has moved in every direction. `n_files` is the
    # expected number of exported time steps. Returns nothing; the assertions are made here
    function check_symmetry(body, solver, freq, n_files)
        root = mktempdir()
        job = Job(body, solver; path=root, freq)
        Peridynamics.set_quiet!(true)
        Peridynamics.submit_threads(job, 1)
        files = filter(endswith(".pvtu"), readdir(joinpath(root, "vtk"); join=true))
        @test length(files) == n_files
        r0 = read_vtk(first(files))
        r = read_vtk(last(files))
        pos = body.position
        @test r0[:position] ≈ pos
        @test iszero(r0[:displacement])
        x = r[:position]
        @test !(x ≈ pos) # something happened
        @test all(abs.(r[:displacement]) .> 0) # every point moved in every direction
        # mirror images in x: same |x| (as a multiset), same y and z
        for d in 1:3
            p = findall(pos[d, :] .> 0)
            m = findall(pos[d, :] .< 0)
            for e in 1:3
                if d == e
                    @test sort(x[e, p]) ≈ sort(-x[e, m])
                else
                    @test sort(x[e, p]) ≈ sort(x[e, m])
                end
            end
        end
        return nothing
    end

    # every material in `materials` pulled apart with velocity conditions, 10 explicit steps
    function velocity_verlet(materials)
        for (name, mat) in materials
            @testset "$name" begin
                Δx = mat isa CKIMaterial ? 0.25 : 0.2
                body = symmetric_cube(mat, Δx)
                symmetric_material!(body, mat, Δx)
                velocity_bc!(t -> 10, body, :set_a, :z)
                velocity_bc!(t -> -10, body, :set_b, :z)
                check_symmetry(body, VelocityVerlet(steps=10), 5, 3)
            end
        end
        return nothing
    end

    # every material in `materials` pulled apart with force density conditions, 10 relaxation steps
    function dynamic_relaxation(materials)
        for (name, mat) in materials
            @testset "$name" begin
                Δx = mat isa CKIMaterial ? 0.25 : 0.2
                body = symmetric_cube(mat, Δx)
                symmetric_material!(body, mat, Δx)
                forcedensity_bc!(t -> 1e10, body, :set_a, :z)
                forcedensity_bc!(t -> -1e10, body, :set_b, :z)
                check_symmetry(body, DynamicRelaxation(steps=10), 5, 3)
            end
        end
        return nothing
    end

    # every supported material in `materials` pulled apart quasi-statically, 2 load steps
    function newton_krylov(materials)
        solver = NewtonKrylov(steps=2, stepsize=1.0, maxiter=25, tol=1e-3)
        for (name, mat) in materials
            mat isa Union{CRMaterial,RKCRMaterial} && continue # not supported by the solver
            @testset "$name" begin
                Δx = 0.25 # 64 points; the Newton-Krylov solve is the expensive part
                body = symmetric_cube(mat, Δx)
                symmetric_material!(body, mat, Δx)
                forcedensity_bc!(p -> 1e10, body, :set_a, :z)
                forcedensity_bc!(p -> -1e10, body, :set_b, :z)
                check_symmetry(body, solver, 1, 3)
            end
        end
        return nothing
    end
end

# Per commit, the solvers run the core materials; the variants (corrections, stabilizations,
# kernels) share the solver interaction with one of them and run in the `extras` job.

@testitem "symmetry: VelocityVerlet" tags=[:simulation] setup=[Fixtures, SymmetryCase] begin
    SymmetryCase.velocity_verlet(Fixtures.CORE_MATERIALS)
end

@testitem "symmetry: VelocityVerlet, material variants" tags=[:simulation, :slow] setup=[Fixtures, SymmetryCase] begin
    SymmetryCase.velocity_verlet(Fixtures.VARIANT_MATERIALS)
end

@testitem "symmetry: DynamicRelaxation" tags=[:simulation] setup=[Fixtures, SymmetryCase] begin
    SymmetryCase.dynamic_relaxation(Fixtures.CORE_MATERIALS)
end

@testitem "symmetry: DynamicRelaxation, material variants" tags=[:simulation, :slow] setup=[Fixtures, SymmetryCase] begin
    SymmetryCase.dynamic_relaxation(Fixtures.VARIANT_MATERIALS)
end

@testitem "symmetry: NewtonKrylov" tags=[:simulation] setup=[Fixtures, SymmetryCase] begin
    # the Newton-Krylov solve is compiled per material and is the expensive part; per commit
    # the bond-based and the correspondence material run here, the ones with a system of their
    # own in the `extras` job
    SymmetryCase.newton_krylov(filter(p -> first(p) in ("BB", "C"), Fixtures.CORE_MATERIALS))
end

@testitem "symmetry: NewtonKrylov, other materials" tags=[:simulation, :slow] setup=[Fixtures, SymmetryCase] begin
    # The solver itself is material-agnostic and the force densities are covered by the
    # VelocityVerlet items, so only the materials with a system of their own run here: the
    # bond-associated (BAC), the interaction (CKI) and the gradient-weight (RKC) system.
    SymmetryCase.newton_krylov(filter(p -> first(p) in ("BAC", "CKI", "RKC-C1"), Fixtures.MATERIALS))
end
