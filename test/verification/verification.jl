@testsnippet Verif begin
    # a snippet does not get the `using Peridynamics` that a test item gets by default
    using Peridynamics
    include(joinpath(pkgdir(Peridynamics), "jobs", "jobs.jl"))

    "Convergence order observed between two grids with spacings `h1`, `h2` and errors `e1`, `e2`."
    observed_order(h1, e1, h2, e2) = log(e1 / e2) / log(h1 / h2)

    """
        interior(chunk, δ)

    Indices of the points further than the horizon `δ` from every surface of the unit cube.

    A point near a surface has an incomplete family and therefore an effective stiffness below
    the bulk value. That surface effect is a property of the formulation and not an error, so a
    comparison against a classical closed form has to be made away from it.
    """
    function interior(chunk, δ)
        X = chunk.system.position
        return [i for i in axes(X, 2) if all(abs(X[d, i]) < 0.5 - δ for d in 1:3)]
    end

    """
        checkerboard!(chunk, amplitude, Δx)

    Displace every point of `chunk` by `±amplitude` in the alternating pattern of the lattice.

    This is the classical zero-energy mode, a deformation a correspondence formulation cannot
    see. A formulation with working stabilization answers it with a restoring force, one without
    answers it with almost nothing.
    """
    function checkerboard!(chunk, amplitude, Δx)
        ref = chunk.system.position
        (; position, displacement) = chunk.storage
        x0 = minimum(@view ref[1, :])
        y0 = minimum(@view ref[2, :])
        z0 = minimum(@view ref[3, :])
        for i in axes(ref, 2)
            # `uniform_box` is cell centered, so `ref[d, i] / Δx` is half-integral and
            # tie-to-even rounding would not alternate. Index from the first lattice point.
            ix = round(Int, (ref[1, i] - x0) / Δx)
            iy = round(Int, (ref[2, i] - y0) / Δx)
            iz = round(Int, (ref[3, i] - z0) / Δx)
            sign = iseven(ix + iy + iz) ? 1.0 : -1.0
            for d in 1:3
                u = amplitude * sign
                position[d, i] = ref[d, i] + u
                displacement[d, i] = u
            end
        end
        return chunk
    end

    "Root mean square of the internal force density over the given points."
    function rms_force(chunk, ids)
        b = chunk.storage.b_int
        return sqrt(sum(sum(abs2, @view b[:, i]) for i in ids) / length(ids))
    end

    """
        front_speed(; mat, npyz, m, t1, t2)

    Speed of the wave front in [`wave_in_bar`](@ref), from its position at two instants.

    Taking a difference of two positions removes the time the pulse needs to be injected, which
    would otherwise show up as a wave that is too slow.
    """
    function front_speed(; mat, npyz=4, m=3.015, t1=1.5e-5, t2=3.0e-5, rel_threshold=0.02)
        x = map((t1, t2)) do time
            dh = submit(wave_in_bar(; mat, npyz, m, time); quiet=true)
            X = dh.chunks[1].system.position
            v = abs.(@view dh.chunks[1].storage.velocity[1, :])
            return maximum(X[1, i] for i in findall(v .> rel_threshold * maximum(v)))
        end
        return (x[2] - x[1]) / (t2 - t1)
    end
end

@testitem "crack tip speed" setup=[Verif] begin
    # A crack in a brittle solid runs at a well defined fraction of the Rayleigh wave speed and
    # cannot exceed it, which is what makes it a verification case: unlike the amount of damage,
    # it is a material property with a physical bound.
    #
    # The last two assertions matter as much as the speed, because a damage model can produce a
    # plausible speed while sending the crack off the symmetry line or letting the tip jump
    # backwards.
    #
    # Measured here: 0.541 (BB), 0.571 (BB-ESC), 0.658 (OSB), 0.654 (C), with the tip within
    # 1.5 Δx of the symmetry line in every case.
    npxy, l = 40, 0.1
    Δx = l / npxy
    c_R = rayleigh_wave_speed(32e9, 0.25, 2500.0)
    # both gauges sit inside the running phase, see `crack_speed`
    gauge_1, gauge_2 = 0.10l, 0.35l
    for (name, mat) in ("BB" => BBMaterial(), "C" => CMaterial())
        times, tip_x, tip_y = mktempdir() do path
            crack_tip_history(mode_i(; mat, npxy, l, path, freq=10); cluster_radius=2Δx)
        end
        v = crack_speed(times, tip_x, gauge_1, gauge_2)
        # Monotonicity over the running phase only: once the crack has arrested the reported tip
        # jitters by a fraction of Δx, because points in the wake cross out of the damage band.
        running = filter(x -> isfinite(x) && x ≤ gauge_2, tip_x)
        off_axis = maximum(filter(isfinite, tip_y))
        @info "crack tip speed" material=name speed=v ratio_to_c_R=v/c_R off_axis
        @test 0.4 < v / c_R < 0.8
        @test off_axis ≤ 2Δx
        @test issorted(running)
    end
end

@testitem "crack tip speed convergence" tags=[:skipci] setup=[Verif] begin
    # The measured speed has to be a property of the material and not of the grid, so refining
    # the plate must not move it much: 0.525 c_R at npxy = 40 against 0.496 c_R at npxy = 60.
    # Tagged `:skipci` because the finer grid takes about fifteen seconds.
    l = 0.1
    c_R = rayleigh_wave_speed(32e9, 0.25, 2500.0)
    ratios = map((40, 60)) do npxy
        Δx = l / npxy
        times, tip_x, _ = mktempdir() do path
            crack_tip_history(mode_i(; npxy, l, path, freq=10); cluster_radius=2Δx)
        end
        r = crack_speed(times, tip_x, 0.10l, 0.35l) / c_R
        @info "crack tip speed convergence" npxy ratio_to_c_R=r
        return r
    end
    @test all(r -> 0.4 < r < 0.8, ratios)
    @test abs(ratios[2] - ratios[1]) / ratios[1] < 0.15
end

@testitem "zero-energy mode stabilization" setup=[Verif] begin
    # Switching the stabilization off is a controlled way to produce the defect it exists to
    # prevent, so the restoring force must collapse when it is off. Measured here: 5.5e9 without
    # against 6.3e11 with the default `ZEMSilling(Cs=0.5)`, a factor of 115.
    #
    # The strain energy density cannot be used instead: for the correspondence materials it is
    # computed from the deformation gradient alone and does not include the stabilization
    # energy, so it is identical whether the stabilization is on or off.
    npyz = 12
    Δx = 1.0 / npyz
    amplitude = 0.05Δx
    function checkerboard_force(mat)
        fixture = chunk_fixture(tension(; mat, npyz))
        chunk = fixture.chunk
        checkerboard!(chunk, amplitude, Δx)
        force_density!(fixture)
        return rms_force(chunk, interior(chunk, 3.015Δx))
    end

    pattern_fixture = chunk_fixture(tension(; mat=CMaterial(), npyz))
    pattern_chunk = checkerboard!(pattern_fixture.chunk, amplitude, Δx)
    signs = reshape(sign.(@view(pattern_chunk.storage.displacement[1, :])), npyz, npyz,
                    npyz)
    @test all(signs[i, j, k] == -signs[i + 1, j, k] for i in 1:npyz-1,
              j in 1:npyz, k in 1:npyz) &&
          all(signs[i, j, k] == -signs[i, j + 1, k] for i in 1:npyz,
              j in 1:npyz-1, k in 1:npyz) &&
          all(signs[i, j, k] == -signs[i, j, k + 1] for i in 1:npyz,
              j in 1:npyz, k in 1:npyz-1)

    unstabilized = checkerboard_force(CMaterial(zem=ZEMSilling(Cs=0.0)))
    for (name, mat) in ("ZEMSilling" => CMaterial(), "ZEMWan" => CMaterial(zem=ZEMWan()),
                        "CRMaterial" => CRMaterial())
        ratio = checkerboard_force(mat) / unstabilized
        @info "zero-energy mode stabilization" stabilization=name ratio
        @test ratio > 10
    end
end

@testitem "wave speed" setup=[Verif] begin
    # A pulse in a slender bar travels at `c₀ = sqrt(E / rho)`. This needs no field beyond the
    # velocity and therefore works for every material.
    #
    # The resolution is coarse so that this can run on every commit, and at that resolution the
    # discretization error is large: the tolerances are the measured errors with headroom, not a
    # statement about the quality of the formulations. That they are discretization errors and
    # not defects is what the convergence item below shows.
    c₀ = sqrt(210e9 / 7850.0)
    TOL = Dict("BB" => 0.20, "OSB" => 0.20, "C" => 0.05)  # measured: -16.9%, +15.4%, +0.5%
    for (name, mat) in ("BB" => BBMaterial(), "OSB" => OSBMaterial(), "C" => CMaterial())
        rel = front_speed(; mat) / c₀ - 1
        @info "wave speed" material=name relative_error=rel
        @test abs(rel) < TOL[name]
    end
end

@testitem "wave speed convergence" tags=[:skipci] setup=[Verif] begin
    # Refining the discretization has to move the wave speed towards `c₀`, which is a stronger
    # statement than any single tolerance. Tagged `:skipci` because the finest grid takes about
    # half a minute.
    c₀ = sqrt(210e9 / 7850.0)
    errors = map((4, 6, 8)) do npyz
        err = abs(front_speed(; mat=BBMaterial(), npyz) / c₀ - 1)
        @info "wave speed convergence" npyz relative_error=err
        return err
    end
    @test issorted(errors; rev=true)
    @test observed_order(1 / 4, errors[1], 1 / 8, errors[3]) > 0
end
