@testitem "time solver interface" begin
    struct MyCustomDataHandler <: Peridynamics.AbstractDataHandler end
    struct MyCustomTimeSolver <: Peridynamics.AbstractTimeSolver end
    struct MyCustomSystem <: Peridynamics.AbstractSystem end
    dh = MyCustomDataHandler()
    solver = MyCustomTimeSolver()
    system = MyCustomSystem()

    @test_throws MethodError Peridynamics.init_time_solver!(solver, dh)

    @test_throws ErrorException Peridynamics.solve!(dh, solver, Dict())

    @test_throws Peridynamics.InterfaceError begin
        Peridynamics.req_point_data_fields_timesolver(Peridynamics.AbstractTimeSolver)
    end

    @test_throws Peridynamics.InterfaceError begin
        Peridynamics.req_bond_data_fields_timesolver(Peridynamics.AbstractTimeSolver)
    end

    @test_throws Peridynamics.InterfaceError begin
        Peridynamics.req_data_fields_timesolver(Peridynamics.AbstractTimeSolver)
    end

    @test isnothing(Peridynamics.init_field_solver(solver, system, Val(:myrandomfield)))
end

@testitem "required fields all timesolvers" begin
    req_fields_solvers = (:position, :displacement, :velocity, :velocity_half,
                          :acceleration, :b_int, :b_ext, :velocity_half_old, :b_int_old,
                          :density_matrix)
    @test Peridynamics.required_fields_timesolvers() === req_fields_solvers
end
