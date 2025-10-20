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

@testitem "msg_path" begin
    # Test 1: Short path - should behave exactly like msg_qty
    short_path = "/Users/test"
    msg1 = Peridynamics.msg_path("root path", short_path)
    msg2 = Peridynamics.msg_qty("root path", short_path)
    @test msg1 == msg2

    # Test 2: Path that fits within default linewidth
    medium_path = "/Users/kaipartmann/Code/Peridynamics.jl"
    msg = Peridynamics.msg_path("root path", medium_path)
    @test occursin("root path", msg)
    @test occursin(medium_path, msg)
    @test count('\n', msg) == 1  # Single line

    # Test 3: Long path requiring line break
    long_path = "/Users/kaipartmann/Code/PeridynamicSimulations/out/btt-cube_GBBMaterial_sigma017/post/ctpx"
    msg = Peridynamics.msg_path("root path", long_path)
    @test occursin("root path", msg)
    # Path is split, so check that parts are present
    @test occursin("/Users/kaipartmann/Code/PeridynamicSimulations/out/", msg)
    @test occursin("btt-cube_GBBMaterial_sigma017/post/ctpx", msg)
    @test count('\n', msg) == 2  # Two lines
    @test occursin("(continued)", msg)

    # Test 4: Verify line splitting at directory separators
    lines = split(msg, '\n'; keepempty=false)
    @test length(lines) == 2
    @test occursin("root path", lines[1])
    @test occursin("(continued)", lines[2])

    # Test 5: Custom linewidth forcing split
    msg = Peridynamics.msg_path("path", long_path; linewidth=50)
    @test occursin("(continued)", msg)
    @test count('\n', msg) >= 2

    # Test 6: Custom indentation
    msg = Peridynamics.msg_path("path", short_path; indentation=5)
    @test startswith(msg, "     ")

    # Test 7: With separator
    msg = Peridynamics.msg_path("path", short_path; separator=":")
    @test occursin("path:", msg)

    # Test 8: With leftwidth parameter
    msg = Peridynamics.msg_path("path", medium_path; leftwidth=20)
    @test occursin("path", msg)
    @test occursin(medium_path, msg)

    # Test 9: Very long path requiring multiple continuation lines
    very_long_path = "/very/long/path/that/keeps/going/and/going/with/many/directories/nested/deeply/in/the/filesystem/structure/with/even/more/subdirectories"
    msg = Peridynamics.msg_path("test", very_long_path; linewidth=60)
    lines = split(msg, '\n'; keepempty=false)
    @test length(lines) >= 2
    # Verify all lines respect linewidth
    for line in lines
        @test length(line) <= 60
    end

    # Test 10: Verify default linewidth respected
    msg = Peridynamics.msg_path("description", long_path)
    lines = split(msg, '\n'; keepempty=false)
    for line in lines
        @test length(line) <= 82  # default_linewidth
    end

    # Test 11: Custom filler character
    msg = Peridynamics.msg_path("path", short_path; filler='-')
    @test occursin("-", msg) || !occursin(".", msg)

    # Test 12: Empty path
    msg = Peridynamics.msg_path("path", "")
    @test occursin("path", msg)

    # Test 13: Path with spaces
    path_with_spaces = "/Users/test user/my documents/file"
    msg = Peridynamics.msg_path("path", path_with_spaces)
    @test occursin(path_with_spaces, msg)

    # Test 14: Windows-style path separator
    windows_path = "C:\\Users\\test\\Documents\\very\\long\\path\\that\\needs\\breaking\\into\\multiple\\lines\\for\\display"
    msg = Peridynamics.msg_path("path", windows_path; linewidth=60)
    # Check that path components are present (may be split)
    @test occursin("Users", msg)
    @test occursin("Documents", msg)
    if count('\n', msg) > 1
        @test occursin("(continued)", msg)
    end

    # Test 15: Path with no directory separators (edge case)
    no_sep_path = "verylongpathwithnoseparatorsatallthatshouldstillbesplitcorrectly"
    msg = Peridynamics.msg_path("test", no_sep_path; linewidth=40)
    # Check that all parts of the path are present (may be split)
    @test occursin("verylongpath", msg)
    @test occursin("rrectly", msg)  # Last part of "correctly"

    # Test 16: Custom indentation with long path
    msg = Peridynamics.msg_path("root path", long_path; indentation=4)
    lines = split(msg, '\n'; keepempty=false)
    @test startswith(lines[1], "    ")
    if length(lines) > 1
        @test startswith(lines[2], "    ")
    end

    # Test 17: Path exactly at boundary
    boundary_path = "/Users/kaipartmann/Code/Project"
    msg = Peridynamics.msg_path("path", boundary_path; linewidth=50)
    @test occursin(boundary_path, msg)

    # Test 18: Custom continuation label
    msg = Peridynamics.msg_path("root path", long_path; continuation_label=">>>")
    @test occursin(">>>", msg)
    @test !occursin("(continued)", msg)
    lines = split(msg, '\n'; keepempty=false)
    @test length(lines) == 2
    @test occursin(">>>", lines[2])

    # Test 19: Custom continuation label with custom indentation
    msg = Peridynamics.msg_path("test", long_path; continuation_label="...", indentation=5)
    @test occursin("...", msg)
    lines = split(msg, '\n'; keepempty=false)
    @test startswith(lines[1], "     ")
    if length(lines) > 1
        @test startswith(lines[2], "     ")
        @test occursin("...", lines[2])
    end

    # Test 20: Empty continuation label
    msg = Peridynamics.msg_path("path", long_path; continuation_label="")
    lines = split(msg, '\n'; keepempty=false)
    @test length(lines) == 2
    # Second line should not have "(continued)" but still have the path

    # Test 21: Long continuation label
    msg = Peridynamics.msg_path("path", long_path; continuation_label="[continuation]")
    @test occursin("[continuation]", msg)

    # Test 22: Test _find_path_break helper with path shorter than budget
    result = Peridynamics._find_path_break("/short", 100)
    @test result == 6

    # Test 23: Test _find_path_break with path longer than budget
    result = Peridynamics._find_path_break("/very/long/path", 8)
    @test result == 6  # Should break at /long/ (the last separator within budget)

    # Test 24: Test _find_path_break with zero or negative budget
    result = Peridynamics._find_path_break("/path", 0)
    @test result == 0

    result = Peridynamics._find_path_break("/path", -5)
    @test result == 0

    # Test 25: Test _find_path_break with no separator in budget
    result = Peridynamics._find_path_break("verylongword", 5)
    @test result == 5

    # Test 26: Whitespace handling in description and path
    msg = Peridynamics.msg_path("  path  ", "  /test/path  ")
    @test occursin("path", msg)
    @test occursin("/test/path", msg)
    @test !occursin("  path  ", msg)  # Should be stripped
end

@testitem "msg_vec" begin
    # Test 1: Short vector - should fit on one line
    short_vec = [1.0, 2.0, 3.0]
    msg = Peridynamics.msg_vec("values", short_vec)
    @test occursin("values", msg)
    @test occursin("[1, 2, 3]", msg)
    @test count('\n', msg) == 1  # Single line

    # Test 2: Empty vector
    empty_vec = Float64[]
    msg = Peridynamics.msg_vec("empty", empty_vec)
    @test occursin("empty", msg)
    @test occursin("[]", msg)
    @test count('\n', msg) == 1

    # Test 3: Single element vector
    single_vec = [42.0]
    msg = Peridynamics.msg_vec("single", single_vec)
    @test occursin("single", msg)
    @test occursin("[42]", msg)

    # Test 4: Long vector requiring line break
    long_vec = collect(1.0:20.0)
    msg = Peridynamics.msg_vec("long values", long_vec)
    @test occursin("long values", msg)
    @test count('\n', msg) >= 2  # Should span multiple lines
    @test occursin("(continued)", msg)

    # Test 5: Verify vector elements are present in long vector
    lines = split(msg, '\n'; keepempty=false)
    @test length(lines) >= 2
    full_msg = join(lines, "")
    @test occursin("1", full_msg)
    @test occursin("20", full_msg)

    # Test 6: Custom linewidth forcing split
    medium_vec = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    msg = Peridynamics.msg_vec("data", medium_vec; linewidth=25)
    # With linewidth=25, this should split
    if count('\n', msg) >= 2
        @test occursin("(continued)", msg)
    end

    # Test 7: Custom indentation
    msg = Peridynamics.msg_vec("vec", short_vec; indentation=5)
    @test startswith(msg, "     ")

    # Test 8: Custom separator
    msg = Peridynamics.msg_vec("vec", short_vec; separator=":")
    @test occursin("vec:", msg)

    # Test 9: Custom filler
    msg = Peridynamics.msg_vec("vec", short_vec; filler='-')
    @test occursin("-", msg) || !occursin(".", msg)

    # Test 10: Custom continuation label
    msg = Peridynamics.msg_vec("data", long_vec; continuation_label=">>>", linewidth=50)
    # Only check if it actually wraps
    if count('\n', msg) >= 2
        @test occursin(">>>", msg)
        @test !occursin("(continued)", msg)
    else
        @test !occursin("(continued)", msg)  # Should not have default label
    end

    # Test 11: Empty continuation label
    msg = Peridynamics.msg_vec("data", long_vec; continuation_label="", linewidth=50)
    lines = split(msg, '\n'; keepempty=false)
    # Should wrap with narrow linewidth
    @test length(lines) >= 1

    # Test 12: Custom vec_delimiter
    msg = Peridynamics.msg_vec("vec", short_vec; vec_delimiter="; ")
    @test occursin("[1; 2; 3]", msg)

    # Test 13: Custom vec_brackets
    msg = Peridynamics.msg_vec("vec", short_vec; vec_brackets=("{", "}"))
    @test occursin("{1, 2, 3}", msg)

    # Test 14: Custom vec_brackets and delimiter
    msg = Peridynamics.msg_vec("vec", short_vec; vec_delimiter=" | ", vec_brackets=("(", ")"))
    @test occursin("(1 | 2 | 3)", msg)

    # Test 15: Very long vector with multiple continuation lines
    very_long_vec = collect(1.0:50.0)
    msg = Peridynamics.msg_vec("array", very_long_vec; linewidth=60)
    lines = split(msg, '\n'; keepempty=false)
    @test length(lines) >= 3
    # Verify all lines respect linewidth
    for line in lines
        @test length(line) <= 60
    end

    # Test 16: Default linewidth respected
    msg = Peridynamics.msg_vec("values", long_vec)
    lines = split(msg, '\n'; keepempty=false)
    for line in lines
        @test length(line) <= 82  # default_linewidth
    end

    # Test 17: Custom leftwidth parameter
    msg = Peridynamics.msg_vec("data", short_vec; leftwidth=20)
    @test occursin("data", msg)
    @test occursin("[1, 2, 3]", msg)

    # Test 18: Vector with floating point numbers
    float_vec = [1.234567, 2.345678, 3.456789]
    msg = Peridynamics.msg_vec("floats", float_vec)
    @test occursin("floats", msg)
    # Should be formatted with %.7g
    @test occursin("1.234567", msg)

    # Test 19: Vector with very small numbers (scientific notation)
    small_vec = [1e-10, 2e-10, 3e-10]
    msg = Peridynamics.msg_vec("small", small_vec)
    @test occursin("small", msg)
    @test occursin("e-", msg)  # Should have scientific notation

    # Test 20: Vector with large numbers (scientific notation)
    large_vec = [1e15, 2e15, 3e15]
    msg = Peridynamics.msg_vec("large", large_vec)
    @test occursin("large", msg)
    @test occursin("e+", msg)  # Should have scientific notation

    # Test 21: Custom indentation with long vector
    msg = Peridynamics.msg_vec("data", long_vec; indentation=4)
    lines = split(msg, '\n'; keepempty=false)
    @test startswith(lines[1], "    ")
    if length(lines) > 1
        @test startswith(lines[2], "    ")
    end

    # Test 22: Integer vector
    int_vec = [1, 2, 3, 4, 5]
    msg = Peridynamics.msg_vec("ints", int_vec)
    @test occursin("ints", msg)
    @test occursin("[1, 2, 3, 4, 5]", msg)

    # Test 23: Test _find_vec_break with short string
    result = Peridynamics._find_vec_break("[1, 2, 3]", 100)
    @test result == 9

    # Test 24: Test _find_vec_break with comma within budget
    result = Peridynamics._find_vec_break("[1, 2, 3, 4]", 7)
    @test result == 7  # Should break after ", " following the comma at position 3

    # Test 25: Test _find_vec_break with zero or negative budget
    result = Peridynamics._find_vec_break("[1, 2]", 0)
    @test result == 0

    result = Peridynamics._find_vec_break("[1, 2]", -5)
    @test result == 0

    # Test 26: Test _find_vec_break with no comma in budget
    result = Peridynamics._find_vec_break("[123456789]", 5)
    @test result == 5

    # Test 27: Whitespace handling in description
    msg = Peridynamics.msg_vec("  data  ", short_vec)
    @test occursin("data", msg)
    @test !occursin("  data  ", msg)  # Should be stripped

    # Test 28: Vector that fits exactly at boundary
    boundary_vec = [1.0, 2.0, 3.0, 4.0]
    msg = Peridynamics.msg_vec("data", boundary_vec; linewidth=35)
    @test occursin("data", msg)

    # Test 29: Test behavior with different numeric types
    mixed_vec = [1.0, 2.5, 3.333333, 4.0]
    msg = Peridynamics.msg_vec("mixed", mixed_vec)
    @test occursin("mixed", msg)
    @test occursin("[", msg)
    @test occursin("]", msg)

    # Test 30: Very long continuation label
    msg = Peridynamics.msg_vec("data", long_vec; continuation_label="[continuation of data]", linewidth=50)
    if count('\n', msg) >= 2
        @test occursin("[continuation of data]", msg)
    end

    # Test 31: Vector with Symbols
    symbol_vec = [:a, :b, :c, :d]
    msg = Peridynamics.msg_vec("symbols", symbol_vec; linewidth=20)
    @test occursin("symbols", msg)
    @test occursin("[:a, :b, :c, :d]", msg)
end
