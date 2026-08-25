# The `@params` macro and the point parameter interface of `src/core/parameters.jl`. The
# declaration language itself is covered in `test_param_fields.jl`.

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

@testitem "@params: linking a hand-written point parameter type" begin
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
    # the error names the type and the method, and its hint shows how to fill the hole
    for (f, name) in ((() -> point_param_type(TestMaterial2()), "point_param_type"),
                      (() -> get_point_params(TestMaterial2(), Dict{Symbol,Any}()),
                       "get_point_params"))
        err = try
            f()
        catch e
            e
        end
        @test err isa InterfaceError
        @test err.type === TestMaterial2
        @test err.func == name
        @test occursin("@params", err.hint)
    end

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

    Peridynamics.@params TestMaterial2 TestPointParameters2
    @test hasmethod(point_param_type, Tuple{TestMaterial2})
    @test Peridynamics.point_param_type(TestMaterial2()) == TestPointParameters2
    # a hand-written type is not generic in the float type and ignores the request
    @test Peridynamics.point_param_type(TestMaterial2(), Float32) == TestPointParameters2
    @test hasmethod(get_point_params, Tuple{TestMaterial2,Dict{Symbol,Any}})

    @test isnothing(macrocheck_input_material(:MyMaterial))
    @test isnothing(macrocheck_input_material(:(MyModule.MyMaterial)))
    @test_throws ArgumentError macrocheck_input_material(:(1 + 1))
    @test isnothing(macrocheck_input_params(:MyParams))
    @test isnothing(macrocheck_input_params(:(MyModule.MyParams)))
    @test_throws ArgumentError macrocheck_input_params(:(1 + 1))
end

@testitem "@params: a material family shares a constructor, not a binding" begin
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection, point_param_type,
                         get_point_params

    # binding a family to a struct would answer `point_param_type` for every material of it
    err = try
        @eval Peridynamics.@params AbstractBondSystemMaterial struct PFFamilyParams
            @inherit StandardParameters
        end
    catch e
        e
    end
    @test err isa LoadError
    @test err.error isa ArgumentError
    @test occursin("material family", err.error.msg)

    # a family constructor serves every material that links the shared type
    struct PFFamilyMat <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::Peridynamics.CriticalStretch
    end
    PFFamilyMat() = PFFamilyMat(Peridynamics.CriticalStretch())
    Peridynamics.@params PFFamilyMat Peridynamics.StandardPointParameters
    @test point_param_type(PFFamilyMat()) === Peridynamics.StandardPointParameters{Float64}
    p = Dict{Symbol,Any}(:horizon => 1.0, :rho => 1.0, :E => 1.0, :nu => 0.25, :Gc => 1.0)
    par = get_point_params(PFFamilyMat(), p)
    @test par isa Peridynamics.StandardPointParameters{Float64}
    @test par.bc ≈ 18 * par.K / (π * par.δ^4)

    # the constructor-only form has to declare exactly the fields of the existing type
    struct PFFamilyMat2 <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::Peridynamics.CriticalStretch
    end
    @test_throws ArgumentError @eval Peridynamics.@params PFFamilyMat2 Peridynamics.StandardPointParameters begin
        @inherit DiscretizationParameters
    end
end

@testitem "instantiate_point_params: generated and hand-written types" begin
    import Peridynamics: instantiate_point_params, StandardPointParameters,
                         AbstractPointParameters

    # a `@params`-generated type is a `UnionAll` and takes the requested float type
    @test instantiate_point_params(StandardPointParameters, Float64) ===
          StandardPointParameters{Float64}
    @test instantiate_point_params(StandardPointParameters, Float32) ===
          StandardPointParameters{Float32}

    # a hand-written concrete type ignores it
    struct PFHandwritten <: AbstractPointParameters
        δ::Float64
    end
    @test instantiate_point_params(PFHandwritten, Float32) === PFHandwritten
end

@testitem "StandardPointParameters: fields, family constructor and float conversion" begin
    import Peridynamics: StandardPointParameters

    @test StandardPointParameters isa UnionAll
    @test fieldnames(StandardPointParameters) == (:δ, :rho, :E, :nu, :G, :K, :λ, :μ, :Gc,
                                                  :εc, :bc)
    @test isbitstype(StandardPointParameters{Float64})

    # the family constructor accepts any material and resolves the standard parameters
    p = Dict{Symbol,Any}(:horizon => 2.0, :rho => 3.0, :E => 1.0, :nu => 0.25,
                         :epsilon_c => 0.01)
    par = StandardPointParameters(OSBMaterial(), p)
    @test par isa StandardPointParameters{Float64}
    @test par.δ == 2.0
    @test par.εc == 0.01
    @test par.Gc ≈ 9.0 / 5.0 * par.K * par.δ * par.εc^2

    par32 = StandardPointParameters{Float32}(par)
    @test par32 isa StandardPointParameters{Float32}
    @test par32.δ === 2.0f0
    @test par32.εc ≈ 0.01f0
end
