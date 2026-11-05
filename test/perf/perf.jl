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
    # The monomial of the RKC family is a type parameter, so the size of the moment matrix
    # in `rkc_weights!` is a compile-time constant and the matrix lives on the stack.
    for (name, mat) in Fixtures.MATERIALS
        Fixtures.has_gradient_weights(mat) || continue
        fixture = Fixtures.material_fixture(mat)
        Fixtures.gradient_weights!(fixture)
        Fixtures.gradient_weights!(fixture)
        bytes = @allocated Fixtures.gradient_weights!(fixture)
        @debug "gradient weight allocations" material=name bytes
        @testset "$name" begin
            @test bytes == 0
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
        fixture = Fixtures.material_fixture(mat)
        @testset "$name" begin
            @test_opt target_modules=(Peridynamics,) Fixtures.force_density!(fixture)
            if Fixtures.has_gradient_weights(mat)
                @test_opt target_modules=(Peridynamics,) Fixtures.gradient_weights!(fixture)
            end
        end
    end
end

@testitem "parameter property forwarding" tags=[:perf] setup=[Fixtures] begin
    # `params.Gc` reads a value that lives in the parameters of the damage model; the
    # generated `getproperty` has to compile to a direct load, so a loop of flat reads
    # allocates nothing and costs the same as reading a field
    body = Fixtures.cube(BBMaterial(); n=4)
    params = only(body.point_params)
    function read_params(params, n)
        s = 0.0
        for _ in 1:n
            s += params.Gc + params.εc + params.δ + params.E
        end
        return s
    end
    read_params(params, 2)
    bytes = @allocated read_params(params, 1000)
    if VERSION ≥ v"1.12"
        @test bytes == 0 # allocates in v1.10
    end
    @test read_params(params, 1) ≈ params.Gc + params.εc + params.δ + params.E
    # the nested access compiles away as well
    nested(params, n) = sum(_ -> params.dmg_params.Gc, 1:n)
    nested(params, 2)
    if VERSION ≥ v"1.12"
        @test (@allocated nested(params, 1000)) == 0 # allocates in v1.10
    end
end
