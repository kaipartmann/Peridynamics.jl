# The `@params` macro and the point parameter interface of `src/core/parameters.jl`. The
# declaration language itself is covered in `test_param_fields.jl`.

@testitem "material declaration: required parameters and the interface without @params" begin
    import Peridynamics: NoCorrection, InterfaceError

    struct TestMaterial1 <: Peridynamics.AbstractBondSystemMaterial{NoCorrection} end
    @test isnothing(Peridynamics.typecheck_material(TestMaterial1))
    @test Peridynamics.required_point_parameters(TestMaterial1) === (:δ, :rho, :E, :nu, :G,
           :K, :λ, :μ)

    # without a `@params` declaration every method of the interface names the macro
    for f in (() -> Peridynamics.point_param_type(TestMaterial1()),
              () -> Peridynamics.get_point_params(TestMaterial1(), Dict{Symbol,Any}()),
              () -> Peridynamics.allowed_material_kwargs(TestMaterial1()))
        err = try
            f()
        catch e
            e
        end
        @test err isa InterfaceError
        @test err.type === TestMaterial1
        @test occursin("@params TestMaterial1 struct", err.hint)
    end

    struct WrongTestMaterial end
    @test_throws ArgumentError Peridynamics.typecheck_material(WrongTestMaterial)

    struct WrongTestMaterial2 <: Peridynamics.AbstractMaterial end
    @test isnothing(Peridynamics.typecheck_material(WrongTestMaterial2))
    @test_throws InterfaceError Peridynamics.required_point_parameters(WrongTestMaterial2)
end

@testitem "typecheck_params: the required parameters of the material family" begin
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection, AbstractPointParameters,
                         typecheck_params

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
    end
    @test isnothing(typecheck_params(TestMaterial2, TestPointParameters2))

    struct PointParametersNoSubtype
        δ::Float64
        rho::Float64
    end
    @test_throws ArgumentError typecheck_params(TestMaterial2, PointParametersNoSubtype)

    struct PointParametersMissingHorizon <: AbstractPointParameters
        rho::Float64
        E::Float64
        nu::Float64
        G::Float64
        K::Float64
        λ::Float64
        μ::Float64
    end
    @test_throws ErrorException typecheck_params(TestMaterial2,
                                                 PointParametersMissingHorizon)
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

@testitem "@params: using the point parameters of another material" setup=[ParamModels] begin
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection, AbstractPointParameters,
                         point_param_type, get_point_params, allowed_material_kwargs

    # the definition binds one material ...
    struct PLBaseMat{D} <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::D
    end
    PLBaseMat() = PLBaseMat(CriticalStretch())
    Peridynamics.@params PLBaseMat struct PLParams
        @inherit StandardParameters
        @log "scale" scale = 1.0
    end

    # ... and the second form binds another one to the same type, with the same keywords
    struct PLOtherMat{D} <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::D
    end
    PLOtherMat() = PLOtherMat(CriticalStretch())
    Peridynamics.@params PLOtherMat PLParams
    @test allowed_material_kwargs(PLOtherMat()) == allowed_material_kwargs(PLBaseMat())
    CSP = Peridynamics.CriticalStretchParameters
    @test point_param_type(PLOtherMat()) === PLParams{Float64,CSP{Float64}}
    @test point_param_type(PLOtherMat(), Float32) === PLParams{Float32,CSP{Float32}}

    p = Dict{Symbol,Any}(:horizon => 1.0, :rho => 1.0, :E => 1.0, :nu => 0.25, :Gc => 1.0,
                         :scale => 2.0)
    par = get_point_params(PLOtherMat(), p)
    @test typeof(par) === point_param_type(PLOtherMat())
    @test par.scale == 2.0
    @test par.bc ≈ 18 * par.K / (π * par.δ^4)

    # the marker fields resolve with the models of the material that uses the type
    mat = PLOtherMat(PMDamage())
    @test point_param_type(mat) === PLParams{Float64,PMDamageParameters{Float64}}
    @test :stretch_scale in Peridynamics.all_material_kwargs(mat)
    par = get_point_params(mat, merge(p, Dict{Symbol,Any}(:stretch_scale => 3.0)))
    @test par.stretch_scale == 3.0

    # the qualified name of a type of the package works as well
    struct PLPackageMat{D} <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::D
    end
    PLPackageMat() = PLPackageMat(CriticalStretch())
    Peridynamics.@params PLPackageMat Peridynamics.OSBPointParameters
    @test point_param_type(PLPackageMat()) === point_param_type(OSBMaterial())

    # only point parameters defined with `@params` can be used for another material
    struct PLHandwritten <: AbstractPointParameters
        δ::Float64
        rho::Float64
        E::Float64
        nu::Float64
        G::Float64
        K::Float64
        λ::Float64
        μ::Float64
    end
    struct PLHandMat <: AbstractBondSystemMaterial{NoCorrection} end
    err = try
        @eval Peridynamics.@params PLHandMat PLHandwritten
    catch e
        e
    end
    @test err isa LoadError
    @test err.error isa ArgumentError
    @test occursin("not defined with `@params`", err.error.msg)
    @test occursin("allowed_material_kwargs", err.error.msg)

    # a name that does not resolve, or resolves to something else
    err = try
        @eval Peridynamics.@params PLHandMat PLNotDefinedAnywhere
    catch e
        e
    end
    @test err isa LoadError
    @test err.error isa ArgumentError
    @test occursin("cannot resolve", err.error.msg)
    err = try
        @eval Peridynamics.@params PLHandMat Float64
    catch e
        e
    end
    @test err isa LoadError
    @test err.error isa ArgumentError
    @test occursin("not a point parameter type", err.error.msg)

    # the material still has to provide what its family requires
    struct PLShortMat <: AbstractBondSystemMaterial{NoCorrection} end
    Peridynamics.required_point_parameters(::Type{PLShortMat}) = (:δ, :rho, :nope)
    @test_throws ErrorException @eval Peridynamics.@params PLShortMat PLParams
end

@testitem "@params: the generated constructor accepts any material" begin
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection, point_param_type,
                         get_point_params

    # the point parameters of one material can be constructed for another, which is what
    # the second form of `@params` relies on; only the declarations read the material
    p = Dict{Symbol,Any}(:horizon => 2.0, :rho => 3.0, :E => 1.0, :nu => 0.25,
                         :epsilon_c => 0.01)
    par = Peridynamics.OSBPointParameters(BBMaterial(), p)
    @test par isa Peridynamics.OSBPointParameters{Float64}
    @test par.δ == 2.0
    @test par.εc == 0.01
    @test par.Gc ≈ 9.0 / 5.0 * par.K * par.δ * par.εc^2
    @test Peridynamics.OSBPointParameters{Float32}(BBMaterial(), p) isa
          Peridynamics.OSBPointParameters{Float32}

    # binding a material family works like binding a family with `@storage`: every
    # material of the family answers with the type unless it declares its own
    abstract type PFTestFam <: AbstractBondSystemMaterial{NoCorrection} end
    Peridynamics.@params PFTestFam struct PFFamParams
        @inherit StandardParameters
    end
    struct PFTestFamMat <: PFTestFam
        dmgmodel::Peridynamics.CriticalStretch
    end
    PFTestFamMat() = PFTestFamMat(Peridynamics.CriticalStretch())
    @test point_param_type(PFTestFamMat()) ===
          PFFamParams{Float64,Peridynamics.CriticalStretchParameters{Float64}}
    @test get_point_params(PFTestFamMat(), p) isa PFFamParams{Float64}
    struct PFTestFamOwnMat <: PFTestFam
        dmgmodel::Peridynamics.CriticalStretch
    end
    PFTestFamOwnMat() = PFTestFamOwnMat(Peridynamics.CriticalStretch())
    Peridynamics.@params PFTestFamOwnMat struct PFFamOwnParams
        @inherit StandardParameters
        own = 1.0
    end
    @test point_param_type(PFTestFamOwnMat()) <: PFFamOwnParams
end

@testitem "BBPointParameters and OSBPointParameters: fields, models and float conversion" begin
    import Peridynamics: BBPointParameters, DHBBPointParameters, OSBPointParameters

    CSP = Peridynamics.CriticalStretchParameters
    @test BBPointParameters isa UnionAll
    @test fieldnames(BBPointParameters) == (:δ, :rho, :E, :nu, :G, :K, :λ, :μ, :bc,
                                            :dmg_params)
    @test fieldnames(OSBPointParameters) == fieldnames(BBPointParameters)
    @test fieldnames(DHBBPointParameters) == fieldnames(BBPointParameters)
    @test isbitstype(BBPointParameters{Float64,CSP{Float64}})
    @test isbitstype(OSBPointParameters{Float64,CSP{Float64}})

    # the bond-based materials fix the Poisson's ratio, the state-based one does not
    p = Dict{Symbol,Any}(:horizon => 2.0, :rho => 3.0, :E => 1.0, :nu => 0.3, :Gc => 1.0)
    @test_throws ArgumentError BBPointParameters(BBMaterial(), p)
    par = OSBPointParameters(OSBMaterial(), p)
    @test par.nu == 0.3
    @test par.Gc == 1.0
    @test par.dmg_params isa CSP{Float64}

    par32 = OSBPointParameters{Float32}(par)
    @test par32 isa OSBPointParameters{Float32,CSP{Float32}}
    @test par32.δ === 2.0f0
    @test par32.Gc === 1.0f0
end

@testitem "@params: FT in the type of a parameter follows the float type" begin
    import Peridynamics: AbstractBondSystemMaterial, NoCorrection, point_param_type,
                         get_point_params
    using Peridynamics.StaticArrays

    struct PFFTMat{D} <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::D
    end
    PFFTMat() = PFFTMat(CriticalStretch())
    Peridynamics.@params PFFTMat struct PFFTParams
        @inherit StandardParameters
        @derived C::SArray{NTuple{4,3},FT,4,81} = Peridynamics.get_hooke_matrix(nu, λ, μ)
        @derived n::SVector{3,FT} = SVector{3,Float64}(δ, rho, E)
    end

    # the field types are built from the float type parameter of the struct
    CSP = Peridynamics.CriticalStretchParameters
    @test fieldtype(PFFTParams{Float64,CSP{Float64}}, :C) ===
          SArray{NTuple{4,3},Float64,4,81}
    @test fieldtype(PFFTParams{Float32,CSP{Float32}}, :C) ===
          SArray{NTuple{4,3},Float32,4,81}
    @test fieldtype(PFFTParams{Float32,CSP{Float32}}, :n) === SVector{3,Float32}
    @test point_param_type(PFFTMat(), Float32) === PFFTParams{Float32,CSP{Float32}}
    @test isbitstype(point_param_type(PFFTMat()))

    p = Dict{Symbol,Any}(:horizon => 2.0, :rho => 3.0, :E => 1.0, :nu => 0.25)
    par = get_point_params(PFFTMat(), p)
    @test par.C isa SArray{NTuple{4,3},Float64,4,81}
    @test par.C ≈ Peridynamics.get_hooke_matrix(par.nu, par.λ, par.μ)
    @test par.n == SVector{3,Float64}(2.0, 3.0, 1.0)

    # the converting constructor converts the element type of such a parameter too
    par32 = PFFTParams{Float32}(par)
    @test par32.C isa SArray{NTuple{4,3},Float32,4,81}
    @test par32.n === SVector{3,Float32}(2.0f0, 3.0f0, 1.0f0)
    @test par32.C ≈ par.C

    # the table shows the type as written
    @test occursin("`SArray{NTuple{4, 3}, FT, 4, 81}`", Peridynamics.block_table(PFFTMat()))

    # the shipped correspondence parameters are declared this way
    @test fieldtype(point_param_type(CMaterial(), Float32), :C) ===
          SArray{NTuple{4,3},Float32,4,81}

    # such a declaration can be inherited from another module: the names of the type
    # expression are resolved where it was written
    struct PFFTInheritMat{D} <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::D
    end
    PFFTInheritMat() = PFFTInheritMat(CriticalStretch())
    Peridynamics.@params PFFTInheritMat struct PFFTInheritParams
        @inherit CPointParameters
    end
    @test fieldtype(point_param_type(PFFTInheritMat(), Float32), :C) ===
          SArray{NTuple{4,3},Float32,4,81}

    # `FT` cannot name a parameter, and every other name in the type has to resolve
    struct PFFTBadMat <: AbstractBondSystemMaterial{NoCorrection} end
    @test_throws LoadError @eval Peridynamics.@params_fields PFFTNameBlock begin
        FT = 1.0
    end
    @test_throws LoadError @eval Peridynamics.@params_fields PFFTGroupNameBlock begin
        @derived (; FT, x) = f()
    end
    err = try
        @eval Peridynamics.@params_fields PFFTUnknownBlock begin
            c::NoSuchArray{FT}
        end
    catch e
        e
    end
    @test err isa LoadError
    @test err.error isa ArgumentError
    @test occursin("NoSuchArray", err.error.msg)
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
    struct PFSuperMat <: Peridynamics.AbstractMaterial end
    Peridynamics.required_point_parameters(::Type{PFSuperMat}) = ()
    Peridynamics.@params PFSuperMat struct PFSuperParams <:
                                            Peridynamics.AbstractPointParameters
        sp_a = 1.0
    end
    @test PFSuperParams <: Peridynamics.AbstractPointParameters
    @test point_param_type(PFSuperMat()) === PFSuperParams{Float64}

    # a definition without declarations is refused
    struct PFEmptyMat <: AbstractBondSystemMaterial{NoCorrection} end
    @test_throws LoadError @eval Peridynamics.@params PFEmptyMat struct PFEmptyParams end

    # the macro has exactly two forms; everything else names them
    for expr in (:(Peridynamics.@params struct PFBare
                       a = 1.0
                   end),
                 :(Peridynamics.@params PFEmptyMat PFBare begin
                       a = 1.0
                   end),
                 :(Peridynamics.@params PFEmptyMat),
                 :(Peridynamics.@params PFEmptyMat 1 + 1))
        err = try
            @eval $(expr)
        catch e
            e
        end
        @test err isa LoadError
        @test err.error isa ArgumentError
        @test occursin("struct MyPointParameters", err.error.msg)
        @test occursin("MyOtherMaterial MyPointParameters", err.error.msg)
    end
    @test_throws ArgumentError Peridynamics.macrocheck_input_material(:(1 + 1))
    @test isnothing(Peridynamics.macrocheck_input_material(:(MyModule.MyMaterial)))

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

    # ... and a material that uses such a type is not either
    struct PFPinnedMat2 <: Peridynamics.AbstractMaterial end
    Peridynamics.required_point_parameters(::Type{PFPinnedMat2}) = ()
    Peridynamics.@params PFPinnedMat2 PFPinnedParams
    @test point_param_type(PFPinnedMat2(), Float32) === PFPinnedParams
    @test get_point_params(PFPinnedMat2(), Dict{Symbol,Any}(:n_substeps => 3)) ===
          PFPinnedParams(3)
end
