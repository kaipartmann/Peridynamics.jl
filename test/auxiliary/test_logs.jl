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
    @test msg_end == "\rDATA HANDLER CREATION COMPLETED ✔\n"

    # reset
    reset_mpi_progress_bars!()
    @test Peridynamics.mpi_progress_bars() == false
    Peridynamics.set_progress_bars!()
    @test Peridynamics.progress_bars() == false
end

@testitem "msg_qty" begin
    msg = Peridynamics.msg_qty("test1", "test2"; linewidth=1)
    @test msg == "  test1 test2\n"
    msg = Peridynamics.msg_qty("   test1   ", "   test2   "; linewidth=1)
    @test msg == "  test1 test2\n"
    msg = Peridynamics.msg_qty("test1", "test2"; linewidth=1, indentation=0)
    @test msg == "test1 test2\n"
    msg = Peridynamics.msg_qty("test1", "test2"; linewidth=1, indentation=5)
    @test msg == "     test1 test2\n"
    msg = Peridynamics.msg_qty("test1", "test2"; linewidth=16, filler='a')
    @test msg == "  test1 aa test2\n"
    msg = Peridynamics.msg_qty("test1", "test2"; linewidth=16)
    @test msg == "  test1 .. test2\n"
    msg = Peridynamics.msg_qty("test1", "test2"; linewidth=16, separator=":")
    @test msg == "  test1: . test2\n"
    msg = Peridynamics.msg_qty("test1", "test2"; linewidth=16, separator="=")
    @test msg == "  test1= . test2\n"
    msg = Peridynamics.msg_qty("test", 1; linewidth=10)
    @test msg == "  test . 1\n"
    msg = Peridynamics.msg_qty("   test  ", 1; linewidth=10)
    @test msg == "  test . 1\n"
    msg = Peridynamics.msg_qty("test", 1; linewidth=9)
    @test msg == "  test  1\n"
    msg = Peridynamics.msg_qty("test", 1; linewidth=8)
    @test msg == "  test 1\n"
    msg = Peridynamics.msg_qty("test", 1; linewidth=1)
    @test msg == "  test 1\n"
    msg = Peridynamics.msg_qty("test", 123456789101112; linewidth=23)
    @test msg == "  test ... 1.234568e+14\n"
    msg = Peridynamics.msg_qty("test", 123456789101112; linewidth=19)
    @test msg == "  test 1.234568e+14\n"
    msg = Peridynamics.msg_qty("test", 123456789101112; linewidth=1)
    @test msg == "  test 1.234568e+14\n"
    msg = Peridynamics.msg_qty("test", 1.123456789101112; linewidth=20)
    @test msg == "  test .... 1.123457\n"
    msg = Peridynamics.msg_qty("test", 1.123456789101112; linewidth=15)
    @test msg == "  test 1.123457\n"
    msg = Peridynamics.msg_qty("test", 1.123456789101112; linewidth=1)
    @test msg == "  test 1.123457\n"
    msg = Peridynamics.msg_qty("test", 1; separator="=", linewidth=15, indentation=4)
    @test msg == "    test= ... 1\n"

    msg = Peridynamics.msg_qty("test1", "test2"; leftwidth=1)
    @test msg == "  test1 test2\n"
    msg = Peridynamics.msg_qty("   test1   ", "   test2   "; leftwidth=1)
    @test msg == "  test1 test2\n"
    msg = Peridynamics.msg_qty("test1", "test2"; leftwidth=1, indentation=0)
    @test msg == "test1 test2\n"
    msg = Peridynamics.msg_qty("test1", "test2"; leftwidth=1, indentation=5)
    @test msg == "     test1 test2\n"
    msg = Peridynamics.msg_qty("test1", "test2"; leftwidth=16, filler='a')
    @test msg == "  test1 aaaaaaa test2\n"
    msg = Peridynamics.msg_qty("test1", "test2"; leftwidth=16)
    @test msg == "  test1 ....... test2\n"
    msg = Peridynamics.msg_qty("test1", "test2"; leftwidth=13, separator=":")
    @test msg == "  test1: ... test2\n"
    msg = Peridynamics.msg_qty("test1", "test2"; leftwidth=11, separator="=")
    @test msg == "  test1= . test2\n"
    msg = Peridynamics.msg_qty("test1", "test2"; leftwidth=10, separator="=")
    @test msg == "  test1=  test2\n"
    msg = Peridynamics.msg_qty("test1", "test2"; leftwidth=9, separator="=")
    @test msg == "  test1= test2\n"

    msg = Peridynamics.msg_qty(:a, :b; linewidth=1)
    @test msg == "  a b\n"
    msg = Peridynamics.msg_qty(1.2345, 2.3456; linewidth=1)
    @test msg == "  1.2345 2.3456\n"

    @test_throws ArgumentError Peridynamics.msg_qty("a", "b"; linewidth=10, leftwidth=9)
end
