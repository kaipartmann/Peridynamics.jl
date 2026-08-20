@testsnippet Verif begin
    # a snippet does not get the `using Peridynamics` that a test item gets by default
    using Peridynamics
    include(joinpath(pkgdir(Peridynamics), "jobs", "jobs.jl"))

    # add a material here and every item below covers it
    MATERIALS = [
        "BB" => BBMaterial(),
        "BB-ESC" => BBMaterial{EnergySurfaceCorrection}(),
        "DHBB" => DHBBMaterial(),
        "GBB" => GBBMaterial(),
        "OSB" => OSBMaterial(),
        "OSB-ESC" => OSBMaterial{EnergySurfaceCorrection}(),
        "C" => CMaterial(),
        "C-ZEMWan" => CMaterial(zem=ZEMWan()),
        "CR" => CRMaterial(),
        "BAC" => BACMaterial(),
        "RKC-C1" => RKCMaterial(monomial=:C1),
        "RKC-RK1" => RKCMaterial(monomial=:RK1),
        "RKCR" => RKCRMaterial(),
        "CKI" => CKIMaterial(),
    ]

    # materials whose crack in `mode_i` can be measured
    FRACTURE_MATERIALS = MATERIALS

    # materials whose bar in `wave_in_bar` can be refined
    WAVE_REFINEMENT_MATERIALS = MATERIALS

    # Materials whose plate in `mode_i` can be refined. `BACMaterial` is out: the plate is only
    # ~4 point spacings thick, so refining changes how much each bond-associated family is
    # truncated. That is the plate, not the formulation, and at `npxy = 75` no material crosses
    # the gauge window anymore.
    FRACTURE_REFINEMENT_MATERIALS = filter(p -> !(last(p) isa BACMaterial), MATERIALS)

    # Materials that average the deformation gradient over the family. The checkerboard is a
    # zero-energy mode of exactly that average, so only these have to answer it.
    correspondence_class(mat) = mat isa Union{Peridynamics.AbstractCorrespondenceMaterial,
                                              Peridynamics.AbstractRKCMaterial,
                                              Peridynamics.AbstractBondAssociatedSystemMaterial}
    CORRESPONDENCE_MATERIALS = filter(p -> correspondence_class(last(p)), MATERIALS)

    # material and the same material without stabilization. Only correspondence materials have
    # a `zem` knob, the others handle zero-energy modes through their quadrature.
    ZEM_MATERIALS = [
        "C-ZEMSilling" => (CMaterial(), CMaterial(zem=ZEMSilling(Cs=0.0))),
        "C-ZEMWan" => (CMaterial(zem=ZEMWan()), CMaterial(zem=ZEMSilling(Cs=0.0))),
        "CR-ZEMSilling" => (CRMaterial(), CRMaterial(zem=ZEMSilling(Cs=0.0))),
        "CR-ZEMWan" => (CRMaterial(zem=ZEMWan()), CRMaterial(zem=ZEMSilling(Cs=0.0))),
    ]

    # interaction systems evaluate triplets, so their interaction count grows with the cube
    horizon_ratio(::Peridynamics.AbstractMaterial) = 3.015
    horizon_ratio(::Peridynamics.AbstractInteractionSystemMaterial) = 2.015

    # convergence order between two grids with spacings `h1`, `h2` and errors `e1`, `e2`
    observed_order(h1, e1, h2, e2) = log(e1 / e2) / log(h1 / h2)

    # Points further than `δ` from every surface of the unit cube. Surface points have incomplete
    # families and thus lower stiffness, a feature and not an error, so comparisons with closed
    # forms have to stay away from them.
    function interior(chunk, δ)
        X = chunk.system.position
        return [i for i in axes(X, 2) if all(abs(X[d, i]) < 0.5 - δ for d in 1:3)]
    end

    # Displaces every point by `±amplitude` in the alternating lattice pattern, the classical
    # zero-energy mode. Imposed as an increment over `Δt` with matching velocities, or the
    # rotated formulations report exactly zero and pass for the wrong reason: they build their
    # stress from the velocity gradient.
    function checkerboard!(chunk, amplitude, Δx, Δt)
        ref = chunk.system.position
        (; position, displacement, velocity, velocity_half) = chunk.storage
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
                velocity[d, i] = u / Δt
                velocity_half[d, i] = u / Δt
            end
        end
        return chunk
    end

    # root mean square of the internal force density over the given points
    function rms_force(chunk, ids)
        b = chunk.storage.b_int
        return sqrt(sum(sum(abs2, @view b[:, i]) for i in ids) / length(ids))
    end

    # Restoring force density against the zero-energy mode, normalized with the bulk scale
    # `E * ε / δ`, `ε = 2 * amplitude / Δx`. Dimensionless and of order one, so all materials
    # are comparable.
    function checkerboard_resistance(mat; npyz::Int=12, E::Real=210e9)
        m = horizon_ratio(mat)
        Δx = 1.0 / npyz
        δ = m * Δx
        amplitude = 0.05Δx
        fixture = chunk_fixture(tension(; mat, npyz, m, E))
        chunk = fixture.chunk
        checkerboard!(chunk, amplitude, Δx, fixture.Δt)
        force_density!(fixture)
        return rms_force(chunk, interior(chunk, δ)) / (E * 2amplitude / Δx / δ)
    end

    # Largest `|v_x|` of every lattice plane of the bar, sorted by position. All chunks are
    # needed, a single chunk does not necessarily span the whole bar.
    function velocity_profile(dh)
        peak = Dict{Float64,Float64}()
        for chunk in dh.chunks
            X = chunk.system.position
            v = chunk.storage.velocity
            for i in 1:Peridynamics.get_n_loc_points(chunk)
                peak[X[1, i]] = max(get(peak, X[1, i], 0.0), abs(v[1, i]))
            end
        end
        positions = sort!(collect(keys(peak)))
        return positions, [peak[x] for x in positions]
    end

    # Position of the wave front, interpolated between the two lattice planes it lies between.
    # Without that it is quantized to the point spacing, which hides the error of the accurate
    # formulations.
    function front_position(dh; rel_threshold=0.02)
        positions, profile = velocity_profile(dh)
        threshold = rel_threshold * maximum(profile)
        k = findlast(≥(threshold), profile)
        (isnothing(k) || k == length(positions)) && return positions[end]
        # the profile decays across the front, so the crossing lies between plane k and k + 1
        w = (profile[k] - threshold) / (profile[k] - profile[k + 1])
        return positions[k] + w * (positions[k + 1] - positions[k])
    end

    # Speed of the wave front from its position at two instants. The difference removes the time
    # the pulse needs to be injected, which would look like a too slow wave.
    function front_speed(; mat, npyz=4, m=horizon_ratio(mat), t1=1.5e-5, t2=3.0e-5)
        x = map((t1, t2)) do time
            return front_position(submit(wave_in_bar(; mat, npyz, m, time); quiet=true))
        end
        return (x[2] - x[1]) / (t2 - t1)
    end

    # Crack speed in `mode_i` over the Rayleigh wave speed, plus how far the tip left the symmetry
    # line, its largest backward step and the number of samples, all only between the two gauges.
    # `freq` is part of the measurement: the speed inherits the sample spacing as an error and a
    # too coarse export can miss the window completely.
    function crack_tip_measurements(mat; npxy::Int=45, l::Real=0.1, freq::Int=4,
                                    gauge_1::Real=0.10l, gauge_2::Real=0.35l)
        Δx = l / npxy
        m = horizon_ratio(mat)
        times, tip_x, tip_y = mktempdir() do path
            job = mode_i(; mat, npxy, l, m, path, freq)
            crack_tip_history(job; cluster_radius=2Δx)
        end
        v = crack_speed(times, tip_x, gauge_1, gauge_2)
        c_R = rayleigh_wave_speed(32e9, 0.25, 2500.0)
        window = [i for i in eachindex(tip_x)
                  if isfinite(tip_x[i]) && gauge_1 ≤ tip_x[i] ≤ gauge_2]
        return (; ratio=v / c_R, off_axis=maximum(tip_y[window]; init=0.0),
                backstep=maximum(-diff(tip_x[window]); init=0.0), n_window=length(window), Δx)
    end

    # Torsional vibration of the free end of `bartwist`: measured angular frequency over the
    # analytical one, the peak angle, and that peak over `Ω_0 / ω`, the amplitude a harmonic
    # oscillator released from `Ω_0` would reach.
    function twist_response(mat; npxy::Int=10, lz::Real=0.5, E::Real=17e6, nu::Real=0.25,
                            rho::Real=1100.0, Ω_0::Real=300.0, freq::Int=5)
        m = horizon_ratio(mat)
        times, angles = mktempdir() do path
            twist_history(bartwist(; mat, npxy, lz, m, E, nu, rho, Ω_0, freq, path))
        end
        r = twist_measurements(times, angles)
        return (; ratio=r.omega / torsional_frequency(E, nu, rho, lz), r.peak,
                peak_ratio=r.peak * r.omega / Ω_0)
    end

    # Runs `job` with `n_chunks` chunks and returns displacement and damage indexed as the body
    # indexes them, so different decompositions can be compared point by point.
    function global_fields(job, n_chunks::Int)
        dh = Peridynamics.threads_data_handler(job.spatial_setup, job.time_solver, n_chunks)
        Peridynamics.init_time_solver!(job.time_solver, dh)
        Peridynamics.initialize!(dh, job.time_solver)
        Peridynamics.solve!(dh, job.time_solver, job.options)
        n_points = length(job.spatial_setup.volume)
        displacement = fill(NaN, 3, n_points)
        damage = fill(NaN, n_points)
        for chunk in dh.chunks
            for (li, i) in enumerate(chunk.system.chunk_handler.loc_points)
                for d in 1:3
                    displacement[d, i] = chunk.storage.displacement[d, li]
                end
                damage[i] = chunk.storage.damage[li]
            end
        end
        # every point of the body has to be the local point of exactly one chunk
        @assert !any(isnan, displacement) && !any(isnan, damage)
        return (; displacement, damage)
    end
end

@testitem "wave speed" tags=[:verification, :skipci] setup=[Verif] begin
    # A pulse in a slender bar travels at `c₀ = sqrt(E / rho)`, which needs nothing but the
    # velocity and works for every material. The grid is coarse, so `TOL` is one loose bound for
    # all materials instead of a table of measured errors; the convergence item shows the rest is
    # discretization error.
    TOL = 0.25
    c₀ = sqrt(210e9 / 7850.0)
    for (name, mat) in MATERIALS
        rel = front_speed(; mat) / c₀ - 1
        @info "wave speed" material=name relative_error=rel
        @test abs(rel) < TOL
    end
end

@testitem "wave speed convergence" tags=[:verification, :skipci] setup=[Verif] begin
    # Refining has to move the wave speed towards `c₀`, which is stronger than any single
    # tolerance. The finest grid alone takes minutes for the RK materials, which makes this and
    # the crack tip convergence item the two long ones of the file.
    #
    # An order needs the discretization error to dominate. Materials already at the noise floor
    # on the coarsest grid cancelled their own error there and get slightly worse when refined,
    # so below `NOISE_FLOOR` only accuracy is asserted.
    NOISE_FLOOR = 0.05
    c₀ = sqrt(210e9 / 7850.0)
    grids = (4, 6, 8)
    for (name, mat) in WAVE_REFINEMENT_MATERIALS
        errors = map(npyz -> abs(front_speed(; mat, npyz) / c₀ - 1), grids)
        order = observed_order(1 / grids[1], errors[1], 1 / grids[end], errors[end])
        @info "wave speed convergence" material=name errors order
        if errors[1] > NOISE_FLOOR
            @test issorted(errors; rev=true)
            @test order > 0
        else
            # already at the floor, so only ask that refining does not take it back out
            @test errors[end] < NOISE_FLOOR
        end
    end
end

@testitem "crack tip speed" tags=[:verification, :skipci] setup=[Verif] begin
    # A crack runs at a well defined fraction of `c_R` and cannot exceed it. Unlike the amount of
    # damage that is a physical bound, which is what makes this a verification case, so the upper
    # bound is `c_R` itself and not the tighter value the measurements would allow.
    #
    # The shape matters as much as the speed: a damage model can produce a plausible speed while
    # sending the crack off the symmetry line. Backward jitter of the tip is fine though, the
    # cluster it is averaged over gains and loses points.
    for (name, mat) in FRACTURE_MATERIALS
        r = crack_tip_measurements(mat)
        @info("crack tip speed", material=name, ratio_to_c_R=r.ratio, n_window=r.n_window,
              off_axis=r.off_axis / r.Δx, backstep=r.backstep / r.Δx)
        # a crack that never crossed both gauges leaves `ratio` at `NaN`, which fails this
        @test 0.4 < r.ratio < 1.0
        @test r.off_axis ≤ 2r.Δx
        @test r.backstep ≤ r.Δx
        # a handful of samples in the window separates a measured speed from a sampling artifact
        @test r.n_window ≥ 5
    end
end

@testitem "crack tip speed convergence" tags=[:verification, :skipci] setup=[Verif] begin
    # The speed has to be a property of the material and not of the grid, so refining the plate
    # must not move it much. The finer plate alone takes a minute and a half for the RK
    # materials.
    #
    # Coarsening instead is no option: at `npxy = 30` the crack of the bond-based materials
    # crosses the whole gauge window inside one export interval, so nothing is sampled in it.
    # `BACMaterial` is not in this list, see `FRACTURE_REFINEMENT_MATERIALS` for why.
    for (name, mat) in FRACTURE_REFINEMENT_MATERIALS
        ratios = map(npxy -> crack_tip_measurements(mat; npxy).ratio, (45, 60))
        change = abs(ratios[2] - ratios[1]) / ratios[1]
        @info "crack tip speed convergence" material=name ratios change
        @test all(r -> 0.4 < r < 1.0, ratios)
        @test change < 0.20
    end
end

@testitem "zero-energy mode resistance" tags=[:verification, :skipci] setup=[Verif] begin
    # Every correspondence-class material has to answer the checkerboard with a restoring force
    # of the order of the bulk stiffness, through a stabilization term, a bond-associated
    # quadrature or a reproducing kernel alike. See `checkerboard_resistance` for the
    # normalization and `CORRESPONDENCE_MATERIALS` for who is asked.
    npyz = 12
    Δx = 1.0 / npyz

    # the pattern has to alternate along every axis, or some other deformation is measured
    pattern_fixture = chunk_fixture(tension(; mat=CMaterial(), npyz))
    pattern_chunk = checkerboard!(pattern_fixture.chunk, 0.05Δx, Δx, pattern_fixture.Δt)
    signs = reshape(sign.(@view(pattern_chunk.storage.displacement[1, :])), npyz, npyz, npyz)
    @test all(signs[i, j, k] == -signs[i + 1, j, k] for i in 1:npyz-1,
              j in 1:npyz, k in 1:npyz) &&
          all(signs[i, j, k] == -signs[i, j + 1, k] for i in 1:npyz,
              j in 1:npyz-1, k in 1:npyz) &&
          all(signs[i, j, k] == -signs[i, j, k + 1] for i in 1:npyz,
              j in 1:npyz, k in 1:npyz-1)

    for (name, mat) in CORRESPONDENCE_MATERIALS
        resistance = checkerboard_resistance(mat)
        @info "zero-energy mode resistance" material=name resistance
        @test resistance > 0.5
    end
end

@testitem "zero-energy mode stabilization" tags=[:verification, :skipci] setup=[Verif] begin
    # Sharper form of the item above for the materials with a `zem` knob: a formulation against
    # itself instead of against a threshold. The strain energy density cannot do this, it comes
    # from the deformation gradient alone and is the same with the stabilization off.
    for (name, (mat, unstabilized)) in ZEM_MATERIALS
        ratio = checkerboard_resistance(mat) / checkerboard_resistance(unstabilized)
        @info "zero-energy mode stabilization" stabilization=name ratio
        @test ratio > 10
    end
end

@testitem "torsional vibration" tags=[:verification, :skipci] setup=[Verif] begin
    # A clamped rod released from the velocity field of its first torsional mode swings at
    # `torsional_frequency`. The free end turns by about 70 degrees, so this is the large rotation
    # case: a formulation that is not objective takes part of that rotation for deformation and
    # the frequency is where it shows up.
    #
    # `peak_ratio` is the sharper check. A conservative vibration released from `Ω_0` reaches the
    # angle `Ω_0 / ω` wherever the discretization puts `ω`, and spurious energy from rotation
    # mistaken for strain breaks that relation. A frequency alone can be blamed on the grid, which
    # is also why its window is wide: the rod is ten points across, so the whole cross section
    # lies within a horizon of a free surface.
    for (name, mat) in MATERIALS
        r = twist_response(mat)
        @info("torsional vibration", material=name, frequency_ratio=r.ratio, peak=r.peak,
              peak_ratio=r.peak_ratio)
        # a rod that never swings back through zero leaves `ratio` at `NaN`, which fails this
        @test 0.7 < r.ratio < 1.3
        @test 0.95 < r.peak_ratio < 1.05
        # the case is only a large rotation case if the end really turns that far
        @test r.peak > 0.9
    end
end

@testitem "decomposition invariance" tags=[:verification, :skipci] setup=[Verif] begin
    # Results must not depend on how the body is cut into chunks. Every chunk builds its own
    # system from local points and a halo, which is where indices that are local in one place and
    # global in another hide from a single-chunk run. `mode_i` and not a uniform body, so the
    # precrack covers the bond and damage bookkeeping across chunk boundaries too.
    #
    # This caught a real bug: a bond-associated system read reference positions from the body with
    # localized indices, which made `BACMaterial` wrong on everything but the first chunk. The MPI
    # comparison cannot see it, it uses one rank per chunk, so both sides are wrong alike.
    npxy, steps = 20, 50
    for (name, mat) in MATERIALS
        reference = global_fields(mode_i(; mat, npxy, m=horizon_ratio(mat), steps), 1)
        scale = maximum(abs, reference.displacement)
        for n_chunks in (2, 4, 8)
            fields = global_fields(mode_i(; mat, npxy, m=horizon_ratio(mat), steps), n_chunks)
            u_diff = maximum(abs, fields.displacement .- reference.displacement) / scale
            d_diff = maximum(abs, fields.damage .- reference.damage)
            @info "decomposition invariance" material=name n_chunks u_diff d_diff
            @test u_diff < 1e-9
            @test d_diff < 1e-12
        end
    end
end

@testitem "uniform tension accuracy" tags=[:verification, :skipci, :simulation] setup=[UniformTension] begin
    # The converged solutions of the bars of `test/simulation/test_analytic_solutions.jl`: the
    # correspondence material reaches the exact elongation within 2 % with both solvers, the
    # bond-based material with surface correction within 10 % at Δx = 1/40. The per-commit
    # items only assert a loose band on shorter runs.
    @test UniformTension.elongation_error(UniformTension.relaxation_bar(1 / 30),
                                          DynamicRelaxation(steps=2000), 2000) < 0.02
    @test UniformTension.elongation_error(UniformTension.newton_bar(1 / 30),
                                          NewtonKrylov(steps=5, tol=1e-3, maxiter=50), 5) < 0.02
    @test UniformTension.elongation_error(UniformTension.databc_bar(1 / 40),
                                          DynamicRelaxation(steps=2000), 2000) < 0.1
end
