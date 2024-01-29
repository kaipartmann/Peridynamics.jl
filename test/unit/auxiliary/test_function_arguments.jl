using Peridynamics, Test

let
    func_method = Peridynamics.get_method_of_function(x -> 2x)
    @test isa(func_method, Method)
end

let
    func_method = Peridynamics.get_method_of_function(x -> 2x)
    @test isa(func_method, Method)
end

let
    f(x::Int) = 2x
    f(x::Float64) = 2.0x
    @test_throws ArgumentError Peridynamics.get_method_of_function(f)
end

let
    func_method = Peridynamics.get_method_of_function(x -> 2x)
    argnames = Peridynamics.get_argument_names_of_function(func_method)
    @test argnames == [:x]
end

let
    func_method = Peridynamics.get_method_of_function(y -> 2y)
    argnames = Peridynamics.get_argument_names_of_function(func_method)
    @test argnames == [:y]
end

let
    func_method = Peridynamics.get_method_of_function((k,t) -> 2t + k)
    argnames = Peridynamics.get_argument_names_of_function(func_method)
    @test argnames == [:k,:t]
end

let
    func_method = Peridynamics.get_method_of_function(() -> rand())
    argnames = Peridynamics.get_argument_names_of_function(func_method)
    @test argnames == Vector{Symbol}()
end

let
    f(a, b, c, d) = a + b + c + d
    func_method = Peridynamics.get_method_of_function(f)
    argnames = Peridynamics.get_argument_names_of_function(func_method)
    @test argnames == [:a, :b, :c, :d]
end

let
    f() = rand()
    func_method = Peridynamics.get_method_of_function(f)
    argnames = Peridynamics.get_argument_names_of_function(func_method)
    @test argnames == Vector{Symbol}()
end
