# `ConditionHandler`: the conditions of a body localized to a chunk, with the degrees of freedom
# the position-dependent (displacement) conditions constrain (`src/discretization/condition_handler.jl`).

@testitem "ConditionHandler: localized point sets and constrained dofs" setup=[Fixtures, BodyCase] begin
    body = BodyCase.tetra()
    material!(body; horizon=2, E=1, rho=1)
    point_set!(body, :a, [1, 2])
    point_set!(body, :b, [3, 4])
    velocity_bc!(Fixtures.f_t, body, :a, :x)
    displacement_bc!(Fixtures.f_p, body, :b, :y) # constrains the y-dof of points 3 and 4
    # a single chunk holds the conditions of the body
    chunk = Fixtures.chunk(body, NewtonKrylov(steps=1))
    ch = chunk.condhandler
    @test ch.loc_point_sets[:a] == [1, 2] && ch.loc_point_sets[:b] == [3, 4]
    @test ch.single_dim_bcs === body.single_dim_bcs
    @test ch.pos_single_dim_bcs === body.pos_single_dim_bcs
    @test isempty(ch.data_bcs)
    # dof = 3(i-1) + d: the y-dofs of points 3 and 4 are 8 and 11
    @test Peridynamics.constrained_dofs(ch) == [8, 11]
    @test sort(Peridynamics.free_dofs(ch)) == setdiff(1:12, [8, 11])
    # two chunks: the point sets and the dofs are local to each chunk
    dh = Fixtures.handler(body, NewtonKrylov(steps=1); n_chunks=2, init=false)
    for chunk in dh.chunks
        local_ch = chunk.condhandler
        n_loc = Peridynamics.get_n_loc_points(chunk.system)
        @test all(all(1 .<= ids .<= n_loc) for ids in values(local_ch.loc_point_sets))
        @test sort(vcat(Peridynamics.free_dofs(local_ch),
                        Peridynamics.constrained_dofs(local_ch))) == 1:(3 * n_loc)
    end
end
