# The `@params` macro and the point parameter interface of `src/core/parameters.jl`. The
# declaration language itself is covered in `test_param_fields.jl`.

@testitem "material declaration: required parameters and allowed kwargs" begin
    import Peridynamics: NoCorrection, InterfaceError

    struct TestMaterial1 <: Peridynamics.AbstractBondSystemMaterial{NoCorrection} end
    @test isnothing(Peridynamics.typecheck_material(TestMaterial1))
    @test Peridynamics.required_point_parameters(TestMaterial1) === (:δ, :rho, :E, :nu, :G,
           :K, :λ, :μ)
    @test Peridynamics.allowed_material_kwargs(TestMaterial1()) === (:horizon, :rho, :E,
           :nu, :G, :K, :lambda, :mu)

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
    @test point_param_type(PFFamilyMat()) ===
          Peridynamics.StandardPointParameters{Float64,
                                               Peridynamics.CriticalStretchParameters{Float64}}
    p = Dict{Symbol,Any}(:horizon => 1.0, :rho => 1.0, :E => 1.0, :nu => 0.25, :Gc => 1.0)
    par = get_point_params(PFFamilyMat(), p)
    @test par isa Peridynamics.StandardPointParameters{Float64}   # partial `isa` still holds
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
    @test fieldnames(StandardPointParameters) == (:δ, :rho, :E, :nu, :G, :K, :λ, :μ, :bc,
                                                  :dmg_params)
    CSP = Peridynamics.CriticalStretchParameters
    @test isbitstype(StandardPointParameters{Float64,CSP{Float64}})

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

@testsnippet ParamModels begin
    # A constitutive model with parameters of its own: the generated stress is the
    # Saint-Venant-Kirchhoff stress scaled by `stiffness_scale`, so every read of the
    # model-owned parameter shows up as an exact scaling. `H_rel` reads the shear modulus
    # `μ` of the material parameters declared above the marker.
    struct PMScaledSVK <: Peridynamics.AbstractConstitutiveModel end
    Peridynamics.@cm_params PMScaledSVK struct PMScaledSVKParameters
        @log "yield stress" sigma_y
        @kwarg hardening_modulus H = 0.0
        @derived H_rel = H / μ
        @log "stiffness scale" stiffness_scale = 1.0
    end
    function Peridynamics.first_piola_kirchhoff(::PMScaledSVK,
                                                storage::Peridynamics.AbstractStorage,
                                                params, F)
        svk = Peridynamics.first_piola_kirchhoff(Peridynamics.SaintVenantKirchhoff(),
                                                 storage, params, F)
        return params.stiffness_scale * svk
    end

    # a damage model with the standard fracture keywords plus one of its own
    struct PMDamage <: Peridynamics.AbstractDamageModel end
    Peridynamics.@dmg_params PMDamage struct PMDamageParameters
        @inherit FractureParameters
        @log "stretch scale" stretch_scale = 1.0
    end
    function Peridynamics.get_frac_params(::PMDamage, δ, K; kwargs...)
        return Peridynamics.get_frac_params(CriticalStretch(), δ, K; kwargs...)
    end
    function Peridynamics.has_fracture(::PMDamage, params)
        return Peridynamics.has_fracture(CriticalStretch(), params)
    end
    function Peridynamics.calc_failure!(storage, system, mat, ::PMDamage, paramsetup, i)
        return Peridynamics.calc_failure!(storage, system, mat, CriticalStretch(),
                                          paramsetup, i)
    end
end

@testitem "@cm_params: the parameters a constitutive model owns" setup=[ParamModels] begin
    using Peridynamics: constitutive_param_type, get_cm_params, constitutive_param_kwargs,
                        convert_nested_params, log_material_parameters

    @test PMScaledSVKParameters isa UnionAll
    @test PMScaledSVKParameters{Float64} <: Peridynamics.AbstractConstitutiveParameters
    @test isbitstype(PMScaledSVKParameters{Float64})
    @test constitutive_param_type(PMScaledSVK(), Float64) === PMScaledSVKParameters{Float64}
    @test constitutive_param_kwargs(PMScaledSVK()) == (:sigma_y, :hardening_modulus,
                                                       :stiffness_scale)

    # the constructor sees the material parameters declared above the marker
    p = Dict{Symbol,Any}(:sigma_y => 300.0, :hardening_modulus => 1.0)
    mp = get_cm_params(PMScaledSVK(), Float64, (; μ=2.0), p)
    @test mp isa PMScaledSVKParameters{Float64}
    @test mp.sigma_y == 300.0
    @test mp.H == 1.0
    @test mp.H_rel == 0.5
    @test mp.stiffness_scale == 1.0

    # a keyword without default is required
    @test_throws UndefKeywordError get_cm_params(PMScaledSVK(), Float64, (; μ=2.0),
                                                 Dict{Symbol,Any}())
    # a missing material-level parameter names the model and the parameter
    err = try
        get_cm_params(PMScaledSVK(), Float64, (;), Dict{Symbol,Any}(:sigma_y => 1.0))
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("PMScaledSVK", err.msg)
    @test occursin("μ", err.msg)

    # float conversion, directly and through the nested-parameter hook
    @test PMScaledSVKParameters{Float32}(mp) isa PMScaledSVKParameters{Float32}
    @test convert_nested_params(Float32, mp) isa PMScaledSVKParameters{Float32}
    @test convert_nested_params(Float32, nothing) === nothing

    # the `@log` labels travel with the model parameters
    msg = log_material_parameters(mp)
    @test occursin("yield stress", msg)
    @test occursin("stiffness scale", msg)
    @test !occursin("H_rel", msg)
end

@testitem "@dmg_params: a damage model registers its own keywords" setup=[ParamModels] begin
    using Peridynamics: point_param_type, damage_param_kwargs, all_material_kwargs

    @test damage_param_kwargs(PMDamage()) == (:Gc, :epsilon_c, :stretch_scale)

    pos, vol = uniform_box(1.0, 1.0, 1.0, 0.5)
    mat = BBMaterial(; dmgmodel=PMDamage())
    @test :stretch_scale in all_material_kwargs(mat)

    body = Body(mat, pos, vol)
    material!(body; horizon=1.5, rho=8e-6, E=2.1e5, Gc=2.7, stretch_scale=2.0)
    par = only(body.point_params)
    @test typeof(par) === point_param_type(mat)
    @test par.dmg_params isa PMDamageParameters{Float64}
    @test par.Gc == 2.7
    @test par.stretch_scale == 2.0
    @test all(body.fail_permit)

    # an unknown keyword is still rejected
    @test_throws ArgumentError material!(body; horizon=1.5, rho=8e-6, E=2.1e5, nope=1)
end

@testitem "@cm_params: the model parameters feed the force density" setup=[ParamModels,
                                                                           Fixtures] begin
    # the scaled Saint-Venant-Kirchhoff model: doubling `stiffness_scale` doubles the
    # forces, which proves the flat read of a model-owned parameter in the force path
    function force_norm(scale)
        body = Fixtures.cube(CMaterial(; model=PMScaledSVK()); n=4, sigma_y=300.0,
                             stiffness_scale=scale)
        dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
        chunk = dh.chunks[1]
        chunk.storage.position .*= 1.001
        Peridynamics.calc_force_density!(chunk, 0.0, 1e-7)
        return copy(chunk.storage.b_int)
    end
    b_ref = force_norm(1.0)
    b_doubled = force_norm(2.0)
    @test !iszero(b_ref)
    @test b_doubled ≈ 2.0 .* b_ref
end

@testitem "@cm_params/@dmg_params: the macro input checks" setup=[ParamModels] begin
    # a model type that is no model is rejected when the definition is evaluated
    struct PMNotAModel end
    @test_throws ArgumentError @eval Peridynamics.@dmg_params PMNotAModel struct PMNMParams
        a = 1.0
    end

    # a model cannot carry the parameters of another model
    @test_throws LoadError @eval Peridynamics.@dmg_params PMDamage struct PMMarkerInside
        dmg_params::Peridynamics.DamageParameters
    end

    # `mat` is not available inside a model parameter block
    @test_throws LoadError @eval Peridynamics.@dmg_params PMDamage struct PMMatRead
        @derived x = float(mat.n_cycles)
    end

    # the struct takes no supertype and no type parameters, and needs declarations
    @test_throws LoadError @eval Peridynamics.@cm_params PMScaledSVK struct PMSuper <:
                                                                            Peridynamics.AbstractConstitutiveParameters
        a = 1.0
    end
    @test_throws LoadError @eval Peridynamics.@cm_params PMScaledSVK struct PMEmpty end

    # a model with parameters needs the marker in the point parameters of the material
    struct PMNoMarkerMat <: Peridynamics.AbstractBondSystemMaterial{Peridynamics.NoCorrection}
        dmgmodel::PMDamage
    end
    PMNoMarkerMat() = PMNoMarkerMat(PMDamage())
    Peridynamics.@params PMNoMarkerMat struct PMNoMarkerParams
        @inherit DiscretizationParameters ElasticParameters
    end
    err = try
        Peridynamics.check_model_params(PMNoMarkerMat())
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("PMDamage", err.msg)
    @test occursin("dmg_params::DamageParameters", err.msg)
end

@testitem "@params: header forms, empty bodies and input checks" begin
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection, point_param_type,
                         get_point_params, allowed_material_kwargs

    # an explicit supertype in the struct header is kept
    Peridynamics.@params struct PFSuperParams <: Peridynamics.AbstractPointParameters
        sp_a = 1.0
    end
    @test PFSuperParams <: Peridynamics.AbstractPointParameters

    # a definition without declarations is refused
    struct PFEmptyMat <: AbstractBondSystemMaterial{NoCorrection} end
    @test_throws LoadError @eval Peridynamics.@params PFEmptyMat struct PFEmptyParams end

    # direct macro input checks
    @test_throws ArgumentError Peridynamics.macrocheck_input_params_block(:(1 + 1))
    @test_throws ArgumentError Peridynamics.macrocheck_input_params_struct(:(1 + 1))

    # a hand-written `point_param_type` without float-type method ignores the request
    struct PFOneArgMat <: AbstractBondSystemMaterial{NoCorrection} end
    struct PFOneArgParams <: Peridynamics.AbstractPointParameters
        δ::Float64
    end
    Peridynamics.point_param_type(::PFOneArgMat) = PFOneArgParams
    @test point_param_type(PFOneArgMat(), Float32) === PFOneArgParams

    # a definition whose parameters all pin their type is not generic in the float type
    struct PFPinnedMat <: Peridynamics.AbstractMaterial end
    Peridynamics.required_point_parameters(::Type{PFPinnedMat}) = ()
    Peridynamics.@params PFPinnedMat struct PFPinnedParams
        n_substeps::Int = 2
    end
    @test !(PFPinnedParams isa UnionAll)
    @test point_param_type(PFPinnedMat()) === PFPinnedParams
    @test point_param_type(PFPinnedMat(), Float32) === PFPinnedParams
    @test allowed_material_kwargs(PFPinnedMat()) == (:n_substeps,)
    par = get_point_params(PFPinnedMat(), Dict{Symbol,Any}())
    @test par === PFPinnedParams(2)

    # a custom material family shares a constructor through the block form
    abstract type PFTestFam <: AbstractBondSystemMaterial{NoCorrection} end
    Peridynamics.@params PFTestFam Peridynamics.StandardPointParameters begin
        @inherit StandardParameters
    end
    struct PFTestFamMat <: PFTestFam
        dmgmodel::Peridynamics.CriticalStretch
    end
    PFTestFamMat() = PFTestFamMat(Peridynamics.CriticalStretch())
    Peridynamics.@params PFTestFamMat Peridynamics.StandardPointParameters
    p = Dict{Symbol,Any}(:horizon => 1.0, :rho => 1.0, :E => 1.0, :nu => 0.25)
    @test get_point_params(PFTestFamMat(), p) isa Peridynamics.StandardPointParameters
end
