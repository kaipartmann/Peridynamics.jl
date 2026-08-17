@testsnippet Perf begin
    # a snippet does not get the `using Peridynamics` that a test item gets by default
    using Peridynamics

    # One entry per material. A new material only needs a line here to be checked.
    MATERIALS = [
        "BB" => BBMaterial(),
        "BB-ESC" => BBMaterial{EnergySurfaceCorrection}(),
        "DHBB" => DHBBMaterial(),
        "GBB" => GBBMaterial(),
        "OSB" => OSBMaterial(),
        "C-ZEMSilling" => CMaterial(model=SaintVenantKirchhoff(), zem=ZEMSilling()),
        "C-ZEMWan" => CMaterial(model=SaintVenantKirchhoff(), zem=ZEMWan()),
        "CR" => CRMaterial(),
        "BAC" => BACMaterial(),
        "RKC" => RKCMaterial(),
        "RKCR" => RKCRMaterial(),
        "CKI" => CKIMaterial(),
    ]

    # Resolution of the cube a material is checked on. Nothing here depends on the size, so it
    # is only large enough to give the interior points a full neighborhood. The interaction
    # system materials evaluate triplets of neighbors instead of pairs, so they get a smaller
    # body and a smaller horizon to keep the runtime in the same order of magnitude.
    body_size(::Peridynamics.AbstractMaterial) = (npyz=6, m=3.015)
    body_size(::Peridynamics.AbstractInteractionSystemMaterial) = (npyz=5, m=2.015)

    """
        material_fixture(mat)

    A deformed body chunk for `mat` together with the time and the step size at which its force
    density is evaluated, exactly as a time solver would see them, but with no threading, MPI,
    file output or logging around it.

    The step size is part of the fixture because the rotated formulations update their stress
    history incrementally, so an arbitrary value would change both the state and the code path.
    """
    function material_fixture(mat)
        npyz, m = body_size(mat)
        l = 1.0
        Δx = l / npyz
        pos, vol = uniform_box(l, l, l, Δx)
        body = Body(mat, pos, vol)
        material!(body; horizon=m*Δx, rho=7850.0, E=210e9, nu=0.25, Gc=2.7)
        # the critical stretch is of the order of 1e-5 here, so any deformation worth measuring
        # would break every bond at once. `calc_failure!` still visits every bond either way,
        # so the code path stays representative.
        no_failure!(body)
        solver = VelocityVerlet(steps=1)
        dh = Peridynamics.threads_data_handler(body, solver, 1)
        Peridynamics.init_time_solver!(solver, dh)
        Peridynamics.initialize!(dh, solver)
        return deform!((chunk=dh.chunks[1], t=solver.Δt, Δt=solver.Δt))
    end

    """
        deform!(fixture; F, amplitude)

    Impose a small non-identity `F` plus an inhomogeneous perturbation on the chunk of `fixture`
    and return it.

    Both parts matter: an undeformed body has zero bond strain and short-circuits, and a purely
    homogeneous deformation makes the zero-energy mode stabilization of the correspondence
    formulations vanish, which is the expensive part of exactly those materials. The
    displacement is treated as an increment over the step size, so the velocity fields are
    populated too, which `CRMaterial` and `RKCRMaterial` need for the velocity gradient.
    """
    function deform!(fixture; F=[1.01 0.003 0.0; 0.0 0.995 0.0; 0.0 0.0 1.0], amplitude=1e-4)
        (; chunk, Δt) = fixture
        ref = chunk.system.position
        (; position, displacement, velocity, velocity_half) = chunk.storage
        for i in axes(ref, 2)
            X1, X2, X3 = ref[1, i], ref[2, i], ref[3, i]
            x1 = F[1, 1] * X1 + F[1, 2] * X2 + F[1, 3] * X3 + amplitude * sinpi(3X2)
            x2 = F[2, 1] * X1 + F[2, 2] * X2 + F[2, 3] * X3 + amplitude * sinpi(3X3)
            x3 = F[3, 1] * X1 + F[3, 2] * X2 + F[3, 3] * X3 + amplitude * sinpi(3X1)
            for (d, x, X) in ((1, x1, X1), (2, x2, X2), (3, x3, X3))
                u = x - X
                position[d, i] = x
                displacement[d, i] = u
                velocity[d, i] = u / Δt
                velocity_half[d, i] = u / Δt
            end
        end
        return fixture
    end

    """
        force_density!(fixture)

    Run the complete force density calculation of `fixture`, the hot loop of every simulation.

    The reproducing kernel materials override `calc_force_density!` for the data handler and not
    for the chunk, because they need a halo exchange between the gradient weights and the
    forces. Calling the chunk method alone would silently skip `calc_weights_and_defgrad!` and
    check half the work, so those materials get their own method here.
    """
    force_density!(fixture) = force_density!(fixture.chunk, fixture.t, fixture.Δt)

    force_density!(chunk, t, Δt) = Peridynamics.calc_force_density!(chunk, t, Δt)

    function force_density!(chunk::Peridynamics.BodyChunk{<:Peridynamics.AbstractBondSystem,
                                                          <:Peridynamics.AbstractRKCMaterial},
                            t, Δt)
        Peridynamics.calc_weights_and_defgrad!(chunk, t, Δt)
        Peridynamics.calc_force_density!(chunk, t, Δt)
        return nothing
    end

    """
        gradient_weights!(fixture)

    Recompute the gradient weights of a reproducing kernel material for every point of `fixture`.

    [`force_density!`](@ref) does not reach this. The weights are only recomputed where damage
    has just grown, and `calc_damage!` rewrites the `update_gradients` flag at the start of
    every force calculation, so an undamaged body never enters it and the cost stays invisible.
    """
    gradient_weights!(fixture) = Peridynamics.initialize!(fixture.chunk)

    "Whether [`gradient_weights!`](@ref) does anything for `mat`."
    has_gradient_weights(mat) = mat isa Peridynamics.AbstractRKCMaterial
end

@testitem "force density allocations" setup=[Perf] begin
    # The force density calculation runs once per point per time step, so anything it allocates
    # is paid for millions of times over a simulation. The target is zero for every material.
    for (name, mat) in MATERIALS
        fixture = material_fixture(mat)
        # twice: once to compile, once more because the first real call also settles the surface
        # corrections and the gradient weights
        force_density!(fixture)
        force_density!(fixture)
        # a fresh fixture, because the rotated materials update their stress history
        fixture = material_fixture(mat)
        bytes = @allocated force_density!(fixture)
        @info "force density allocations" material=name bytes
        @test bytes == 0
    end
end

@testitem "gradient weight allocations" setup=[Perf] begin
    # `force_density!` never reaches the gradient weight calculation, so it is measured
    # separately here. See `gradient_weights!` in the `Perf` snippet.
    #
    # These allocate for a known reason, see the type stability item below. `@test_broken` keeps
    # them visible in every run, and once the cause is fixed Julia reports an unexpected pass.
    for (name, mat) in MATERIALS
        has_gradient_weights(mat) || continue
        fixture = material_fixture(mat)
        gradient_weights!(fixture)
        gradient_weights!(fixture)
        bytes = @allocated gradient_weights!(fixture)
        @info "gradient weight allocations" material=name bytes
        @test_broken bytes == 0
    end
end

@testitem "force density type stability" tags=[:skipci] setup=[Perf] begin
    # JET is tied closely to the compiler internals of a given Julia version while CI runs 1.10,
    # 1.11 and 1.12, so this is tagged `:skipci` and run deliberately:
    #
    #     julia -t 6 test/runtestitems.jl "type stability"
    #
    # `Test.@inferred` is useless on `force_density!`: it compares the inferred against the
    # actual *return* type, and that is `Nothing` no matter how unstable the body is.
    using JET
    for (name, mat) in MATERIALS
        # `RKCMaterial` carries `monomial` as a `Symbol` field rather than as a type parameter,
        # so `q_dim = get_q_dim(monomial)` in `rkc_weights!` is a runtime value and the moment
        # matrix built from it cannot be stack allocated. JET names it exactly:
        # `get_q_dim(%2::Val)::Any` and `zero(Type{SMatrix{_A,_B,Float64,_C}})::Any`.
        broken = mat isa Peridynamics.AbstractRKCMaterial
        fixture = material_fixture(mat)
        @testset "$name" begin
            @test_opt broken=broken target_modules=(Peridynamics,) force_density!(fixture)
            if has_gradient_weights(mat)
                @test_opt broken=broken target_modules=(Peridynamics,) gradient_weights!(fixture)
            end
        end
    end
end
