# Quasi-static problems with a closed-form solution: a bar under uniform tension with its left
# end fixed, solved with the dynamic relaxation and the Newton-Krylov solver and with the data
# boundary conditions. These per-commit runs are short and coarse and assert a loose band; the
# accuracy of the converged solutions is asserted in the verification suite (see
# `test/verification/verification.jl`, "uniform tension accuracy"), which shares this module.

@testmodule UniformTension begin
    using Peridynamics

    "Length and cross section of the bar, the pulling force and the Young's modulus."
    const l, w, h = 1.0, 0.1, 0.1
    const F = 2e6
    const E = 200e9

    "The exact elongation of the bar: `F l / (E A)`."
    analytic_elongation() = F / (E * w * h) * l

    """
        bar(mat, Δx; kwargs...)

    The bar of `mat` with spacing `Δx`, three extra layers on the left for the fixed end, the
    point sets `:left` (fixed) and `:right` (loaded), and the force density that distributes `F`
    over the right layer. `kwargs` are passed to `material!`.
    """
    function bar(mat, Δx; kwargs...)
        pos, vol = uniform_box(l + 3Δx, w, h, Δx; center=(-1.5Δx, 0, 0))
        body = Body(mat, pos, vol)
        material!(body; horizon=3.015Δx, rho=7850.0, E, nu=0.25, kwargs...)
        point_set!(x -> x < -l / 2, body, :left)
        point_set!(x -> x > l / 2 - Δx, body, :right)
        b_right = F / sum(vol[body.point_sets[:right]])
        return body, b_right
    end

    "Run `body` with `solver`, export the last step, and return the mean x-displacement of `:right`."
    function elongation(body, solver, steps)
        path = mktempdir()
        job = Job(body, solver; freq=steps, path, fields=(:displacement,))
        submit(job; quiet=true)
        results = read_vtk(last(Peridynamics.find_vtk_files(joinpath(path, "vtk"))))
        right = body.point_sets[:right]
        return sum(results[:displacement][1, right]) / length(right)
    end

    "Relative error of the elongation of `body` solved with `solver`."
    function elongation_error(body, solver, steps)
        Δl = analytic_elongation()
        return abs(elongation(body, solver, steps) - Δl) / Δl
    end

    "The bar with velocity conditions fixing the left end: for the dynamic relaxation."
    function relaxation_bar(Δx)
        body, b_right = bar(CMaterial(), Δx; epsilon_c=1.0)
        for dim in (:x, :y, :z)
            velocity_bc!(t -> 0.0, body, :left, dim)
        end
        forcedensity_bc!(t -> b_right, body, :right, :x)
        return body
    end

    "The bar with displacement conditions fixing the left end: for the Newton-Krylov solver."
    function newton_bar(Δx)
        body, b_right = bar(CMaterial(zem=ZEMSilling(Cs=0.5)), Δx)
        for dim in (:x, :y, :z)
            displacement_bc!(p -> 0.0, body, :left, dim)
        end
        forcedensity_bc!(p -> b_right, body, :right, :x)
        return body
    end

    "The bar with data boundary conditions (one value per point) for both ends."
    function databc_bar(Δx)
        body, b_right = bar(BBMaterial{EnergySurfaceCorrection}(), Δx; epsilon_c=1.0)
        # one row per dimension, one column per point; only the points of the set are read
        f_matrix = zeros(1, body.n_points)
        v_matrix = zeros(3, body.n_points)
        for i in body.point_sets[:right]
            f_matrix[1, i] = b_right
        end
        Peridynamics.velocity_databc!(body, v_matrix, :left, [1, 2, 3])
        Peridynamics.forcedensity_databc!(body, f_matrix, :right, [:x])
        return body
    end
end

@testitem "uniform tension: DynamicRelaxation" tags=[:simulation] setup=[UniformTension] begin
    # 1000 relaxation steps on 297 points: the elongation is within a few percent of the exact
    # one (measured 2.2 %; 2000 steps give 0.5 %, see the verification suite)
    body = UniformTension.relaxation_bar(1 / 30)
    @test UniformTension.elongation_error(body, DynamicRelaxation(steps=1000), 1000) < 0.05
end

@testitem "uniform tension: NewtonKrylov" tags=[:simulation] setup=[UniformTension] begin
    # five load steps on 297 points converge to the exact elongation (measured 0.5 %)
    body = UniformTension.newton_bar(1 / 30)
    solver = NewtonKrylov(steps=5, tol=1e-3, maxiter=50)
    @test UniformTension.elongation_error(body, solver, 5) < 0.02
end

@testitem "uniform tension: data boundary conditions" tags=[:simulation] setup=[UniformTension] begin
    # the same bar with the conditions given as data, one value per point and dimension; the
    # bond-based material with surface correction converges to 13 % at this spacing (8.6 % at
    # Δx = 1/40, see the verification suite), so this is a band on a converged coarse solution
    body = UniformTension.databc_bar(1 / 30)
    # the data conditions are shown with the body
    msg = sprint(show, MIME("text/plain"), body; context=:compact => false)
    @test contains(msg, "2 boundary condition(s):")
    @test contains(msg, "Data BC on velocity: point_set=left, dims=UInt8[0x01, 0x02, 0x03]")
    @test contains(msg, "Data BC on force density: point_set=right, dims=UInt8[0x01]")
    @test UniformTension.elongation_error(body, DynamicRelaxation(steps=1000), 1000) < 0.2
end
