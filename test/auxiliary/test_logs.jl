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
    @test_logs (:warn,) (Peridynamics.set_progress_bars!())
    @test Peridynamics.progress_bars() == true
    reset_mpi_progress_bars!()
    @test Peridynamics.mpi_progress_bars() == false
    Peridynamics.set_progress_bars!()
    @test Peridynamics.progress_bars() == false
end

@testitem "log_create_data_handler" begin
    # setup
    Peridynamics.set_quiet!(false)
    @test Peridynamics.quiet() == false
    force_mpi_run!()
    @test Peridynamics.mpi_run() == true
    enable_mpi_progress_bars!()
    @test Peridynamics.mpi_progress_bars() == true
    @test_logs (:warn,) (Peridynamics.set_progress_bars!())
    @test Peridynamics.progress_bars() == true

    io = IOBuffer()
    Peridynamics.log_create_data_handler_start(io)
    msg_start = String(take!(io))
    @test msg_start == "DATA HANDLER CREATION ... ⏳"
    Peridynamics.log_create_data_handler_end(io)
    msg_end = String(take!(io))
    @test msg_end == "\rDATA HANDLER CREATION COMPLETED ✓\n"

    # reset
    reset_mpi_progress_bars!()
    @test Peridynamics.mpi_progress_bars() == false
    Peridynamics.set_progress_bars!()
    @test Peridynamics.progress_bars() == false
end
