# Performance guarantees of the hot loops: no allocations and type stability for every material.
# These items are `:perf` and run in the `extras` CI job. The fixtures live in `test/setup/`.

@testitem "force density allocations" tags=[:perf] setup=[Fixtures] begin
    # The force density calculation runs once per point per time step, so anything it allocates
    # is paid for millions of times over a simulation. The target is zero for every material.
    for (name, mat) in Fixtures.MATERIALS
        fixture = Fixtures.material_fixture(mat)
        # twice: once to compile, once more because the first real call also settles the surface
        # corrections and the gradient weights
        Fixtures.force_density!(fixture)
        Fixtures.force_density!(fixture)
        # a fresh fixture, because the rotated materials update their stress history
        fixture = Fixtures.material_fixture(mat)
        bytes = @allocated Fixtures.force_density!(fixture)
        @debug "force density allocations" material=name bytes
        @testset "$name" begin
            @test bytes == 0
        end
    end
end

@testitem "gradient weight allocations" tags=[:perf] setup=[Fixtures] begin
    # `force_density!` never reaches the gradient weight calculation, so it is measured
    # separately here. See `gradient_weights!` in `test/setup/fixtures.jl`.
    #
    # These allocate for a known reason, see the type stability item below. `@test_broken` keeps
    # them visible in every run, and once the cause is fixed Julia reports an unexpected pass.
    for (name, mat) in Fixtures.MATERIALS
        Fixtures.has_gradient_weights(mat) || continue
        fixture = Fixtures.material_fixture(mat)
        Fixtures.gradient_weights!(fixture)
        Fixtures.gradient_weights!(fixture)
        bytes = @allocated Fixtures.gradient_weights!(fixture)
        @debug "gradient weight allocations" material=name bytes
        @testset "$name" begin
            @test_broken bytes == 0
        end
    end
end

@testitem "force density type stability" tags=[:perf, :skipci] setup=[Fixtures] begin
    # JET is tied closely to the compiler internals of a given Julia version while CI runs 1.10,
    # 1.11 and 1.12, so this is tagged `:skipci` and run deliberately:
    #
    #     julia -t 6 test/runtestitems.jl "type stability"
    #
    # `Test.@inferred` is useless on `force_density!`: it compares the inferred against the
    # actual *return* type, and that is `Nothing` no matter how unstable the body is.
    using JET
    for (name, mat) in Fixtures.MATERIALS
        # `RKCMaterial` carries `monomial` as a `Symbol` field rather than as a type parameter,
        # so `q_dim = get_q_dim(monomial)` in `rkc_weights!` is a runtime value and the moment
        # matrix built from it cannot be stack allocated. JET names it exactly:
        # `get_q_dim(%2::Val)::Any` and `zero(Type{SMatrix{_A,_B,Float64,_C}})::Any`.
        broken = mat isa Peridynamics.AbstractRKCMaterial
        fixture = Fixtures.material_fixture(mat)
        @testset "$name" begin
            @test_opt broken=broken target_modules=(Peridynamics,) Fixtures.force_density!(fixture)
            if Fixtures.has_gradient_weights(mat)
                @test_opt broken=broken target_modules=(Peridynamics,) Fixtures.gradient_weights!(fixture)
            end
        end
    end
end
