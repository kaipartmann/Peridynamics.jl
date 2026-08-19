# The `@params` macro and the point parameter interface of `src/core/parameters.jl`.

@testitem "material declaration: required parameters and allowed kwargs" begin
    import Peridynamics: NoCorrection, InterfaceError

    struct TestMaterial1 <: Peridynamics.AbstractBondSystemMaterial{NoCorrection} end
    @test isnothing(Peridynamics.typecheck_material(TestMaterial1))
    @test Peridynamics.required_point_parameters(TestMaterial1) === (:δ, :rho, :E, :nu, :G,
           :K, :λ, :μ)
    @test Peridynamics.allowed_material_kwargs(TestMaterial1()) === (:horizon, :rho, :E,
           :nu, :G, :K, :lambda, :mu, :Gc, :epsilon_c)

    struct WrongTestMaterial end
    @test_throws ArgumentError Peridynamics.typecheck_material(WrongTestMaterial)

    struct WrongTestMaterial2 <: Peridynamics.AbstractMaterial end
    @test isnothing(Peridynamics.typecheck_material(WrongTestMaterial2))
    @test_throws InterfaceError Peridynamics.required_point_parameters(WrongTestMaterial2)
    @test_throws InterfaceError Peridynamics.allowed_material_kwargs(WrongTestMaterial2())
end

@testitem "@params: point parameter declaration" begin
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection, AbstractPointParameters,
                         InterfaceError, typecheck_params, constructor_check,
                         point_param_type, get_point_params, macrocheck_input_material,
                         macrocheck_input_params
    struct TestMaterial2 <: AbstractBondSystemMaterial{NoCorrection} end
    struct TestPointParameters2 <: AbstractPointParameters
        δ::Float64
        rho::Float64
        E::Float64
        nu::Float64
        G::Float64
        K::Float64
        λ::Float64
        μ::Float64
        Gc::Float64
        εc::Float64
    end
    tpp2 = TestPointParameters2(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    TestPointParameters2(::TestMaterial2, ::Dict{Symbol,Any}) = nothing

    @test isnothing(typecheck_params(TestMaterial2, TestPointParameters2))
    ie1 = InterfaceError(TestMaterial2, "point_param_type")
    @test_throws ie1 point_param_type(TestMaterial2())
    ie3 = InterfaceError(TestMaterial2, "get_point_params")
    @test_throws ie3 get_point_params(TestMaterial2(), Dict{Symbol,Any}())

    struct PointParametersNoSubtype
        δ::Float64
        rho::Float64
        E::Float64
        nu::Float64
        G::Float64
        K::Float64
        λ::Float64
        μ::Float64
        Gc::Float64
        εc::Float64
    end
    @test_throws ArgumentError typecheck_params(TestMaterial2, PointParametersNoSubtype)
    @test_throws InterfaceError constructor_check(TestMaterial2, PointParametersNoSubtype)

    struct PointParametersMissingHorizon <: AbstractPointParameters
        rho::Float64
        E::Float64
        nu::Float64
        G::Float64
        K::Float64
        λ::Float64
        μ::Float64
        Gc::Float64
        εc::Float64
    end
    @test_throws ErrorException typecheck_params(TestMaterial2,
                                                 PointParametersMissingHorizon)

    @test_throws InterfaceError point_param_type(TestMaterial2())

    Peridynamics.@params TestMaterial2 TestPointParameters2
    @test hasmethod(point_param_type, Tuple{TestMaterial2})
    @test Peridynamics.point_param_type(TestMaterial2()) == TestPointParameters2
    @test hasmethod(get_point_params, Tuple{TestMaterial2,Dict{Symbol,Any}})

    @test isnothing(macrocheck_input_material(:MyMaterial))
    @test isnothing(macrocheck_input_material(:(MyModule.MyMaterial)))
    @test_throws ArgumentError macrocheck_input_material(:(1 + 1))
    @test isnothing(macrocheck_input_params(:MyParams))
    @test isnothing(macrocheck_input_params(:(MyModule.MyParams)))
    @test_throws ArgumentError macrocheck_input_params(:(1 + 1))
end
