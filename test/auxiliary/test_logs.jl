@testitem "quiet" begin
    Peridynamics.set_quiet!(true)
    @test Peridynamics.quiet() == true
    Peridynamics.set_quiet!(false)
    @test Peridynamics.quiet() == false
end

@testitem "set_progress_bars!" begin
    Peridynamics.set_quiet!(false)
    @test Peridynamics.quiet() == false
    force_mpi_run!()
    @test Peridynamics.mpi_run() == true
    enable_mpi_progress_bars!()
    @test Peridynamics.mpi_progress_bars() == true
    @test Peridynamics.progress_bars() == false
    @test_logs (:warn,) (Peridynamics.set_progress_bars!())
    @test Peridynamics.progress_bars() == true
    reset_mpi_progress_bars!()
    @test Peridynamics.mpi_progress_bars() == false
    Peridynamics.set_progress_bars!()
    @test Peridynamics.progress_bars() == false
end
