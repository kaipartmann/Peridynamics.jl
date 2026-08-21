@testitem "get_method_of_function" begin
    func_method = Peridynamics.get_method_of_function(x -> 2x)
    @test isa(func_method, Method)

    func_method = Peridynamics.get_method_of_function(x -> 2x)
    @test isa(func_method, Method)

    f(x::Int) = 2x
    f(x::Float64) = 2.0x
    @test_throws ArgumentError Peridynamics.get_method_of_function(f)
end

@testitem "get_argument_names_of_function" begin # lint-ok: rand is only the body of an example function
    func_method = Peridynamics.get_method_of_function(x -> 2x)
    argnames = Peridynamics.get_argument_names_of_function(func_method)
    @test argnames == [:x]

    func_method = Peridynamics.get_method_of_function(y -> 2y)
    argnames = Peridynamics.get_argument_names_of_function(func_method)
    @test argnames == [:y]

    func_method = Peridynamics.get_method_of_function((k,t) -> 2t + k)
    argnames = Peridynamics.get_argument_names_of_function(func_method)
    @test argnames == [:k,:t]

    func_method = Peridynamics.get_method_of_function(() -> rand())
    argnames = Peridynamics.get_argument_names_of_function(func_method)
    @test argnames == Vector{Symbol}()

    f(a, b, c, d) = a + b + c + d
    func_method = Peridynamics.get_method_of_function(f)
    argnames = Peridynamics.get_argument_names_of_function(func_method)
    @test argnames == [:a, :b, :c, :d]

    h() = rand()
    func_method = Peridynamics.get_method_of_function(h)
    argnames = Peridynamics.get_argument_names_of_function(func_method)
    @test argnames == Vector{Symbol}()
end

@testitem "check_kwargs" begin
    @test Peridynamics.check_kwargs(Dict{Symbol,Any}(:a => 1), (:a, :b)) === nothing
    @test_throws ArgumentError Peridynamics.check_kwargs(Dict{Symbol,Any}(:c => 1), (:a, :b))
end

@testitem "check_kwargs: an unknown keyword names the closest allowed one" begin
    allowed = (:horizon, :rho, :E, :nu, :G, :K, :lambda, :mu, :Gc, :epsilon_c)

    # Levenshtein distance: insertions, deletions and substitutions count one each
    @test Peridynamics.kwarg_distance("epsilon_c", "epsilon_c") == 0
    @test Peridynamics.kwarg_distance("epsilon_C", "epsilon_c") == 1
    @test Peridynamics.kwarg_distance("bta", "beta") == 1
    @test Peridynamics.kwarg_distance("horizn", "horizon") == 1
    @test Peridynamics.kwarg_distance("", "rho") == 3

    # a close miss is suggested, a far one is not (no suggestion beats a wrong one)
    @test Peridynamics.closest_kwarg(:epsilon_C, allowed) === :epsilon_c
    @test Peridynamics.closest_kwarg(:horizn, allowed) === :horizon
    @test Peridynamics.closest_kwarg(:Gcc, allowed) === :Gc
    # a tie (`Ec` is one edit from both `E` and `Gc`) goes to the first allowed keyword
    @test Peridynamics.closest_kwarg(:Ec, allowed) === :E
    @test Peridynamics.closest_kwarg(:density, allowed) === nothing
    @test Peridynamics.closest_kwarg(:rho, ()) === nothing

    msg = Peridynamics.unknown_kwarg_msg(:epsilon_C, allowed)
    @test contains(msg, "keyword `epsilon_C` not allowed")
    @test contains(msg, "Did you mean `epsilon_c`?")
    @test contains(msg, "allowed keywords: horizon, rho, E")
    msg = Peridynamics.unknown_kwarg_msg(:density, allowed)
    @test !contains(msg, "Did you mean")
    @test contains(msg, "allowed keywords:")

    # and `check_kwargs` throws exactly that message
    err = try
        Peridynamics.check_kwargs(Dict{Symbol,Any}(:epsilon_C => 1), allowed)
    catch e
        e
    end
    @test err isa ArgumentError
    @test contains(err.msg, "Did you mean `epsilon_c`?")
end
