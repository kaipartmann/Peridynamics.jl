# Shared fixtures for the whole test suite.
#
# A `@testmodule` is evaluated once per test process, unlike a `@testsnippet`, whose source is
# included into every single test item that uses it and therefore recompiled every time. Test
# items opt in with `setup=[Fixtures]` and access everything qualified, e.g.
#
#     @testitem "velocity_bc!: ..." setup=[Fixtures] begin
#         body = Fixtures.tetra4(BBMaterial())
#         for (name, mat) in Fixtures.MATERIALS ... end
#     end
#
# Nothing is exported on purpose, so short helper names like `chunk` can never shadow or collide
# with variables of a test item. The module body only defines things: no `@test`, no I/O.
#
# `@testmodule`s cannot depend on each other (the macro takes no `setup` keyword), so this is the
# single module for fixtures. Its parts live in plain files next to it:
#   materials.jl        the material tables every material-generic test loops over
#   fixtures.jl         bodies, chunks, data handlers and deformations
#   study_fixtures.jl   job creators for the `Study` tests
@testmodule Fixtures begin
    using Peridynamics, Test, Random
    using Peridynamics.StaticArrays, Peridynamics.LinearAlgebra

    include(joinpath(@__DIR__, "materials.jl"))
    include(joinpath(@__DIR__, "fixtures.jl"))
    include(joinpath(@__DIR__, "study_fixtures.jl"))
end
