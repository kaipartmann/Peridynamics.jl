# The parameter declaration framework of `src/core/param_fields.jl`: `@params_fields`,
# the `@kwarg`/`@derived`/`@log` annotations and the code `@params` generates from the
# declarations.

@testitem "@params: the generated parametric type, constructor and keywords" begin
    using Peridynamics: NoCorrection, AbstractBondSystemMaterial, AbstractDamageModel,
                        CriticalStretch, point_param_type, get_point_params,
                        allowed_material_kwargs

    struct PFMat1{D} <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::D
    end
    PFMat1(; dmgmodel::AbstractDamageModel=CriticalStretch()) = PFMat1(dmgmodel)

    Peridynamics.@params PFMat1 struct PFParams1
        @inherit StandardParameters
        sigma_y = Inf
        hardening = 0.0
        n_substeps::Int = 1
    end

    # the struct is parametric in the float type and in the model parameter type of the
    # damage-model marker of `StandardParameters`, and a pinned type stays what it is
    @test PFParams1 isa UnionAll
    @test fieldnames(PFParams1) == (:δ, :rho, :E, :nu, :G, :K, :λ, :μ, :bc,
                                    :dmg_params, :sigma_y, :hardening, :n_substeps)
    @test fieldtype(PFParams1{Float64}, :δ) === Float64
    @test fieldtype(PFParams1{Float32}, :δ) === Float32
    @test fieldtype(PFParams1{Float32}, :n_substeps) === Int

    # `point_param_type` has to answer with a concrete type, `Body` is parameterized with
    # it; the marker parameters are answered by the models of the material instance
    CSP = Peridynamics.CriticalStretchParameters
    @test point_param_type(PFMat1()) === PFParams1{Float64,CSP{Float64}}
    @test point_param_type(PFMat1(), Float32) === PFParams1{Float32,CSP{Float32}}
    @test isbitstype(point_param_type(PFMat1()))

    # the allowed keywords follow from the declarations; the fracture keywords belong to
    # the damage model now and arrive through `all_material_kwargs`
    @test allowed_material_kwargs(PFMat1()) == (:horizon, :rho, :E, :nu, :G, :K, :lambda,
                                                :mu, :sigma_y, :hardening, :n_substeps)
    @test Peridynamics.all_material_kwargs(PFMat1()) == (:horizon, :rho, :E, :nu, :G, :K,
                                                         :lambda, :mu, :sigma_y,
                                                         :hardening, :n_substeps, :Gc,
                                                         :epsilon_c)

    p = Dict{Symbol,Any}(:horizon => 1.5, :rho => 8e-6, :E => 2.1e5, :nu => 0.25,
                         :Gc => 2.7, :sigma_y => 300.0)
    par = get_point_params(PFMat1(), p)
    @test par isa PFParams1{Float64}
    @test par.δ == 1.5
    @test par.rho == 8e-6
    @test par.E == 2.1e5
    @test par.nu == 0.25
    @test par.K ≈ 2.1e5 / (3 * (1 - 2 * 0.25))
    @test par.bc ≈ 18 * par.K / (π * par.δ^4)   # `@derived` sees the parameters before it
    @test par.sigma_y == 300.0
    @test par.hardening == 0.0                  # default of an unspecified keyword
    @test par.n_substeps === 1

    # converting the whole set of parameters to another float type
    par32 = PFParams1{Float32}(par)
    @test par32 isa PFParams1{Float32}
    @test par32.δ === 1.5f0
    @test par32.n_substeps === 1
end

@testitem "@params: declaration and keyword errors" begin
    using Peridynamics: NoCorrection, AbstractBondSystemMaterial, CriticalStretch,
                        get_point_params, check_material_kwargs

    struct PFMat2 <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::Peridynamics.CriticalStretch
    end
    PFMat2() = PFMat2(CriticalStretch())

    Peridynamics.@params PFMat2 struct PFParams2
        @inherit StandardParameters
    end

    # a required keyword that is not given names itself
    p = Dict{Symbol,Any}(:rho => 1.0, :E => 1.0, :nu => 0.25)
    @test_throws UndefKeywordError get_point_params(PFMat2(), p)

    # a keyword that no declaration consumes is rejected
    @test_throws ArgumentError check_material_kwargs(PFMat2(), Dict{Symbol,Any}(:nope => 1))

    # point parameters may not declare their own type parameters
    @test_throws LoadError @eval Peridynamics.@params PFMat2 struct PFParamsBad{T}
        @inherit StandardParameters
    end

    # a group of parameters is computed and therefore has to say `@derived`
    @test_throws LoadError @eval Peridynamics.@params_fields PFFieldsBareGroup begin
        (; E, nu) = Peridynamics.get_elastic_params(; E, nu)
    end

    # only shorthand keyword arguments, so that "after `;` is a `material!` keyword" holds
    @test_throws LoadError @eval Peridynamics.@params_fields PFFieldsExplicitKwarg begin
        @derived (; E, nu) = Peridynamics.get_elastic_params(; E=1.0, nu)
    end

    # `p`, the keyword dictionary, is not part of the declaration language any more
    @test_throws LoadError @eval Peridynamics.@params_fields PFFieldsKeywordDict begin
        @derived (; E, nu) = Peridynamics.get_elastic_params(p)
    end

    # a group cannot carry an annotation that applies to a single parameter
    @test_throws LoadError @eval Peridynamics.@params_fields PFFieldsGroupLabel begin
        @log "elastic" @derived (; E, nu) = Peridynamics.get_elastic_params(; E, nu)
    end

    # a right-hand side sees only what is declared above it
    @test_throws LoadError @eval Peridynamics.@params_fields PFFieldsForwardRef begin
        @derived twice = 2 * later
        later = 1.0
    end
end

@testitem "@derived: parameter groups supplied by one call" begin
    using Peridynamics: param_fields_expr, get_point_params, allowed_material_kwargs,
                        NoCorrection, AbstractBondSystemMaterial, CriticalStretch

    two_of(; first_kw, second_kw=nothing) = (; a=float(first_kw),
                                             b=isnothing(second_kw) ? 0.0 : float(second_kw))
    one_of(; third_kw=nothing) = (; c=isnothing(third_kw) ? -1.0 : float(third_kw))

    # several groups in one body, in any order, mixed with ordinary declarations
    Peridynamics.@params_fields PFGroups begin
        @derived (; a, b) = two_of(; first_kw, second_kw)
        d = 4.0
        @derived (; c) = one_of(; third_kw)
        @derived s = a + b + c + d
    end
    # the fields follow the declaration order, a group in the position of its `@derived`
    spec = param_fields_expr(PFGroups)
    @test [decl.name for decl in spec.decls] == [:a, :b, :d, :c, :s]

    # the allowed keywords are exactly the ones written after the `;` of the calls, plus the
    # keywords of the ordinary declarations
    @test spec.kwargs == [:first_kw, :second_kw, :d, :third_kw]

    # a group member is not a keyword of its own
    @test all(decl.kwarg === :none for decl in spec.decls if decl.name in (:a, :b, :c, :s))

    struct PFMatGroups <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::CriticalStretch
    end
    PFMatGroups() = PFMatGroups(CriticalStretch())

    Peridynamics.@params PFMatGroups struct PFParamsGroups
        @inherit StandardParameters
        @inherit PFGroups
    end

    @test :first_kw in allowed_material_kwargs(PFMatGroups())
    p = Dict{Symbol,Any}(:horizon => 1.0, :rho => 1.0, :E => 1.0, :nu => 0.25, :Gc => 1.0,
                         :first_kw => 2.0)
    par = get_point_params(PFMatGroups(), p)
    @test par.a == 2.0
    @test par.b == 0.0    # a keyword that was not given is not forwarded to the call
    @test par.c == -1.0
    @test par.s ≈ 2.0 + 0.0 + (-1.0) + 4.0

    # a member of a group is labelled by naming it, because the group declares no label
    Peridynamics.@params_fields PFGroupLabels begin
        @derived (; a, b) = two_of(; first_kw)
        @log "the first one" a
    end
    labels = Dict(decl.name => decl.label for decl in param_fields_expr(PFGroupLabels).decls)
    @test labels[:a] == "the first one"
    @test labels[:b] == ""
end

@testitem "@params_fields: inheritance and merge rules" begin
    using Peridynamics: ParamFieldDecl, param_fields_expr

    Peridynamics.@params_fields PFBlockA begin
        a = 1.0
        b = 2.0
    end

    Peridynamics.@params_fields PFBlockB begin
        @inherit PFBlockA
        c = 3.0
    end

    @test [d.name for d in param_fields_expr(PFBlockB).decls] == [:a, :b, :c]
    @test param_fields_expr(PFBlockB).kwargs == [:a, :b, :c]

    # a declaration in the body overrides an inherited one and keeps its position
    Peridynamics.@params_fields PFBlockC begin
        @inherit PFBlockB
        @derived b = 42.0
    end
    decls = param_fields_expr(PFBlockC).decls
    @test [d.name for d in decls] == [:a, :b, :c]
    @test decls[2].kwarg === :none                      # `@derived` is not a keyword
    @test param_fields_expr(PFBlockC).kwargs == [:a, :b, :c]

    # two blocks that declare the same parameter differently cannot both be inherited
    Peridynamics.@params_fields PFBlockD begin
        b = 99.0
    end
    @test_throws LoadError @eval Peridynamics.@params_fields PFBlockConflict begin
        @inherit PFBlockA PFBlockD
    end

    # the point parameters of a material are a block as well, and a declaration in the
    # body overrides an inherited one in place: the dual-horizon parameters are the
    # bond-based ones with another bond constant, in the same order
    bb = param_fields_expr(Peridynamics.BBPointParameters).decls
    dhbb = param_fields_expr(Peridynamics.DHBBPointParameters).decls
    @test [d.name for d in dhbb] == [d.name for d in bb]
    @test count(a != b for (a, b) in zip(bb, dhbb)) == 1
    @test only(d for d in dhbb if d.name === :bc).source == "(0.5 * 18 * K) / (π * δ ^ 4)"

    # only `@params` and `@params_fields` definitions can be inherited from
    @test_throws ArgumentError param_fields_expr(Float64)
end

@testitem "@params: the annotations @kwarg, @derived and @log" begin
    using Peridynamics: NoCorrection, AbstractBondSystemMaterial, CriticalStretch,
                        get_point_params, allowed_material_kwargs, log_material_parameters

    struct PFMat3 <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::Peridynamics.CriticalStretch
    end
    PFMat3() = PFMat3(CriticalStretch())

    Peridynamics.@params PFMat3 struct PFParams3
        @inherit StandardParameters
        @kwarg gamma_c gammac = 1e-10
        @log "yield stress" sigma_y = Inf
        @derived twice_bc = 2 * bc
    end

    # `@kwarg` renames the keyword, `@derived` removes it
    kwargs = allowed_material_kwargs(PFMat3())
    @test :gamma_c in kwargs
    @test !(:gammac in kwargs)
    @test :sigma_y in kwargs
    @test !(:twice_bc in kwargs)

    p = Dict{Symbol,Any}(:horizon => 1.0, :rho => 1.0, :E => 1.0, :nu => 0.25, :Gc => 1.0,
                         :gamma_c => 0.5)
    par = get_point_params(PFMat3(), p)
    @test par.gammac == 0.5
    @test par.sigma_y == Inf
    @test par.twice_bc ≈ 2 * par.bc

    # `@log` puts the parameter into the simulation log, an unlabelled one stays out
    msg = log_material_parameters(par)
    @test occursin("yield stress", msg)
    @test !occursin("twice_bc", msg)
end

@testitem "point_param_type: shipped parameters are concrete, isbits and float-generic" begin
    using Peridynamics: point_param_type, BBPointParameters, DHBBPointParameters,
                        OSBPointParameters, CPointParameters, BACPointParameters,
                        CKIPointParameters, RKCPointParameters
    using Peridynamics.StaticArrays: SArray

    materials = (BBMaterial(), DHBBMaterial(), GBBMaterial(), OSBMaterial(), CMaterial(),
                 CRMaterial(), BACMaterial(), CKIMaterial(), RKCMaterial(), RKCRMaterial())
    for mat in materials
        P = point_param_type(mat)
        # concrete, so that `Body` and `BodyChunk` stay concrete, and isbits, so that the
        # parameters can be captured by value in a kernel
        @test isconcretetype(P)
        @test isbitstype(P)
        @test P === point_param_type(mat, Float64)
        @test point_param_type(mat, Float32) !== P
    end

    # every shipped point parameter type is generic in the float type and carries the
    # marker for the damage model parameters; the marker for the constitutive model exists
    # exactly where the material carries one — the correspondence family
    for P in (BBPointParameters, DHBBPointParameters, OSBPointParameters, CPointParameters,
              BACPointParameters, CKIPointParameters, RKCPointParameters)
        @test P isa UnionAll
        @test Peridynamics.has_dmg_param_marker(P)
    end
    for P in (CPointParameters, BACPointParameters, RKCPointParameters)
        @test Peridynamics.has_cm_param_marker(P)
    end
    for P in (BBPointParameters, DHBBPointParameters, OSBPointParameters, CKIPointParameters)
        @test !Peridynamics.has_cm_param_marker(P)
    end

    # materials with identical parameters share the type, the others have their own
    @test point_param_type(GBBMaterial()) === point_param_type(BBMaterial())
    @test point_param_type(CRMaterial()) === point_param_type(CMaterial())
    @test point_param_type(RKCRMaterial()) === point_param_type(RKCMaterial())
    @test point_param_type(DHBBMaterial()) !== point_param_type(BBMaterial())

    # the stiffness tensor of the correspondence parameters follows the float type
    @test fieldtype(point_param_type(CMaterial(), Float32), :C) ===
          SArray{NTuple{4,3},Float32,4,81}
end

@testitem "material!: point parameters of every material" begin
    using Peridynamics: get_point_params

    pos, vol = uniform_box(1.0, 1.0, 1.0, 0.5)

    # `nu` is not given for the bond-based materials, which fix it at 1/4
    cases = ((BBMaterial(), (;)), (DHBBMaterial(), (;)), (GBBMaterial(), (;)),
             (OSBMaterial(), (; nu=0.25)), (CMaterial(), (; nu=0.25)),
             (CRMaterial(), (; nu=0.25)), (BACMaterial(), (; nu=0.25)),
             (CKIMaterial(), (; nu=0.25)), (RKCMaterial(), (; nu=0.25)),
             (RKCRMaterial(), (; nu=0.25)))
    for (mat, extra) in cases
        body = Body(mat, pos, vol)
        material!(body; horizon=1.5, rho=8e-6, E=2.1e5, Gc=2.7, extra...)
        par = only(body.point_params)
        @test par.δ == 1.5
        @test par.rho == 8e-6
        @test par.E == 2.1e5
        @test par.nu ≈ 0.25
        @test par.Gc == 2.7
        @test par.K ≈ 2.1e5 / (3 * (1 - 2 * 0.25))
        @test typeof(par) === Peridynamics.point_param_type(mat)
    end

    # the bond constant of the dual-horizon model is half of the standard one
    bc(mat) = (b = Body(mat, pos, vol);
               material!(b; horizon=1.5, rho=8e-6, E=2.1e5, Gc=2.7);
               only(b.point_params).bc)
    @test bc(DHBBMaterial()) ≈ 0.5 * bc(BBMaterial())

    # bond-based peridynamics rejects any Poisson's ratio other than 1/4
    for mat in (BBMaterial(), DHBBMaterial(), GBBMaterial())
        body = Body(mat, pos, vol)
        @test_throws ArgumentError material!(body; horizon=1.5, rho=8e-6, E=2.1e5, nu=0.3,
                                             Gc=2.7)
    end

    # the bond horizon of a bond-associated material defaults to the horizon
    body = Body(BACMaterial(), pos, vol)
    material!(body; horizon=1.5, rho=8e-6, E=2.1e5, nu=0.25, Gc=2.7)
    @test only(body.point_params).δb == 1.5
    body = Body(BACMaterial(), pos, vol)
    material!(body; horizon=1.5, bond_horizon=2.0, rho=8e-6, E=2.1e5, nu=0.25, Gc=2.7)
    @test only(body.point_params).δb == 2.0
end

@testitem "@params: the marker fields of the model parameters" begin
    using Peridynamics: param_fields_expr, is_cm_param_decl, is_dmg_param_decl, block_table

    # a marker is one bare field of a block or a definition and registers no keyword
    Peridynamics.@params_fields PFMarkerBlock begin
        a = 1.0
        cm_params::Peridynamics.ConstitutiveParameters
        dmg_params::Peridynamics.DamageParameters
    end
    spec = param_fields_expr(PFMarkerBlock)
    @test [d.name for d in spec.decls] == [:a, :cm_params, :dmg_params]
    @test spec.kwargs == [:a]
    @test is_cm_param_decl(spec.decls[2])
    @test is_dmg_param_decl(spec.decls[3])

    # the table says who owns the declarations behind the marker
    table = block_table(PFMarkerBlock)
    @test occursin("owned by the constitutive model", table)
    @test occursin("owned by the damage model", table)

    # one marker per kind, and nothing of the declaration language applies to a marker
    @test_throws LoadError @eval Peridynamics.@params_fields PFMarkerTwice begin
        one::Peridynamics.DamageParameters
        two::Peridynamics.DamageParameters
    end
    @test_throws LoadError @eval Peridynamics.@params_fields PFMarkerDefault begin
        dmg_params::Peridynamics.DamageParameters = 1.0
    end
    @test_throws LoadError @eval Peridynamics.@params_fields PFMarkerDerived begin
        @derived dmg_params::Peridynamics.DamageParameters = 1.0
    end
end

@testitem "getproperty: flat reads over the marker fields" begin
    pos, vol = uniform_box(1.0, 1.0, 1.0, 0.5)
    body = Body(BBMaterial(), pos, vol)
    material!(body; horizon=1.5, rho=8e-6, E=2.1e5, Gc=2.7)
    par = only(body.point_params)

    # `Gc` physically lives in the parameters of the damage model
    @test par.Gc == 2.7
    @test par.Gc === Base.getfield(Base.getfield(par, :dmg_params), :Gc)
    @test par.dmg_params isa Peridynamics.CriticalStretchParameters
    # a bond-based material has no constitutive model, so its parameters have no slot
    @test !hasproperty(par, :cm_params)
    # a correspondence material has one, empty for the parameterless standard model
    rkc_body = Body(RKCMaterial(), pos, vol)
    material!(rkc_body; horizon=1.5, rho=8e-6, E=2.1e5, nu=0.25, Gc=2.7)
    @test only(rkc_body.point_params).cm_params === nothing
    @test :Gc in propertynames(par)
    @test hasproperty(par, :εc)
    @test hasproperty(par, :dmg_params)
    @test !hasproperty(par, :notaparameter)
    @test_throws Exception par.notaparameter

    # a parameter name that exists in the material and in a model is rejected when the
    # material is set up
    struct PFCollide <: Peridynamics.AbstractDamageModel end
    Peridynamics.@dmg_params PFCollide struct PFCollideParameters
        δ = 1.0
    end
    body2 = Body(BBMaterial(; dmgmodel=PFCollide()), pos, vol)
    err = try
        material!(body2; horizon=1.5, rho=8e-6, E=2.1e5)
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("δ", err.msg)
    @test occursin("more than once", err.msg)
end

@testitem "@params_fields: the error paths of the declaration parser" begin
    import Peridynamics: param_fields_expr

    # the annotations exist only inside a definition
    @test_throws ArgumentError @eval Peridynamics.@derived x = 1.0
    @test_throws ArgumentError @eval Peridynamics.@kwarg kw x
    @test_throws ArgumentError @eval Peridynamics.@log "label" x

    # invalid macro inputs
    @test_throws LoadError @eval Peridynamics.@params_fields (1 + 1) begin
        a = 1.0
    end
    @test_throws LoadError @eval Peridynamics.@params_fields PFBadBlock 42

    # things that are no parameter declaration
    @test_throws LoadError @eval Peridynamics.@params_fields PFLiteral begin
        1.0
    end
    @test_throws LoadError @eval Peridynamics.@params_fields PFCallDecl begin
        some_call(1.0)
    end

    # `@inherit` needs a resolvable type that registered declarations
    @test_throws LoadError @eval Peridynamics.@params_fields PFInheritUnknown begin
        @inherit NoSuchBlockName123
    end
    @test_throws LoadError @eval Peridynamics.@params_fields PFInheritValue begin
        @inherit π
    end

    # a group needs a call on the right-hand side, and members are `name` or `name::Type`
    @test_throws LoadError @eval Peridynamics.@params_fields PFGroupNoCall begin
        @derived (; a, b) = 1.0
    end
    @test_throws LoadError @eval Peridynamics.@params_fields PFGroupBadMember begin
        @derived (; a, f(b)) = Peridynamics.get_horizon(; horizon)
    end
    # a typed member pins its type
    Peridynamics.@params_fields PFGroupTypedMember begin
        @derived (; gtm_a::Int, gtm_b) = pf_typed_member_provider(; gtm_kw)
    end
    spec = param_fields_expr(PFGroupTypedMember)
    @test spec.decls[1].type === Int

    # explicit keyword values are rejected in the single-declaration path as well
    @test_throws LoadError @eval Peridynamics.@params_fields PFExplicitSingle begin
        @derived x = Peridynamics.get_horizon(; horizon=1.0)
    end

    # unknown and malformed annotations
    @test_throws LoadError @eval Peridynamics.@params_fields PFUnknownAnnotation begin
        @inbounds x = 1.0
    end
    @test_throws LoadError @eval Peridynamics.@params_fields PFBadDerived begin
        @derived x y
    end
    @test_throws LoadError @eval Peridynamics.@params_fields PFBadKwarg begin
        @kwarg onlyone
    end
    @test_throws LoadError @eval Peridynamics.@params_fields PFBadLog begin
        @log missing_label x
    end

    # nothing of the declaration language applies to a marker field
    @test_throws LoadError @eval Peridynamics.@params_fields PFMarkerDerivedBare begin
        @derived dmg_params::Peridynamics.DamageParameters
    end
    @test_throws LoadError @eval Peridynamics.@params_fields PFMarkerKwarg begin
        @kwarg kw dmg_params::Peridynamics.DamageParameters
    end
    @test_throws LoadError @eval Peridynamics.@params_fields PFMarkerLabel begin
        @log "nope" dmg_params::Peridynamics.DamageParameters
    end
    # one marker per kind, also for the constitutive model
    @test_throws LoadError @eval Peridynamics.@params_fields PFCmTwice begin
        one::Peridynamics.ConstitutiveParameters
        two::Peridynamics.ConstitutiveParameters
    end

    # a pinned type has to be resolvable
    @test_throws LoadError @eval Peridynamics.@params_fields PFBadType begin
        x::NoSuchType123 = 1.0
    end

    # a qualified callee that does not resolve is kept as written
    Peridynamics.@params_fields PFUnresolvedQualified begin
        @derived uq = Base.no_such_function_123(1.0)
    end
    @test occursin("no_such_function_123", param_fields_expr(PFUnresolvedQualified).decls[1].source)

    # direct input checks of the model parameter macros
    @test isnothing(Peridynamics.macrocheck_input_model(:(SomeModule.SomeModel)))
    @test_throws ArgumentError Peridynamics.macrocheck_input_model(:(1 + 1))
    @test_throws ArgumentError Peridynamics.macrocheck_input_model_params(:NotAStruct)
    @test_throws ArgumentError Peridynamics.typecheck_model_params(:cm, 1)
end

@testitem "getproperty: an ambiguous parameter name names both owners" begin
    using Peridynamics: NoCorrection, AbstractBondSystemMaterial, get_point_params

    struct PFAmbCM <: Peridynamics.AbstractConstitutiveModel end
    Peridynamics.@cm_params PFAmbCM struct PFAmbCMParams
        q = 1.0
    end
    struct PFAmbDM <: Peridynamics.AbstractDamageModel end
    Peridynamics.@dmg_params PFAmbDM struct PFAmbDMParams
        q = 2.0
    end
    struct PFAmbMat <: AbstractBondSystemMaterial{NoCorrection}
        dmgmodel::PFAmbDM
    end
    PFAmbMat() = PFAmbMat(PFAmbDM())
    Peridynamics.get_constitutive_model(::PFAmbMat) = PFAmbCM()
    Peridynamics.@params PFAmbMat struct PFAmbParams
        @inherit StandardParameters
        cm_params::Peridynamics.ConstitutiveParameters
    end

    # constructing directly bypasses `material!` and its collision check
    p = Dict{Symbol,Any}(:horizon => 1.0, :rho => 1.0, :E => 1.0, :nu => 0.25)
    par = get_point_params(PFAmbMat(), p)
    @test par.cm_params.q == 1.0
    @test par.dmg_params.q == 2.0
    err = try
        par.q
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("cm_params", err.msg)
    @test occursin("dmg_params", err.msg)

    # `material!` refuses the setup up front
    pos, vol = uniform_box(1.0, 1.0, 1.0, 0.5)
    body = Body(PFAmbMat(), pos, vol)
    @test_throws ArgumentError material!(body; horizon=1.0, rho=1.0, E=1.0, nu=0.25)
end

@testitem "@dmg_params: pinned-type parameters and the paramless defaults" begin
    import Peridynamics: damage_param_type, damage_param_kwargs, get_dmg_params,
                         convert_nested_params

    # a model whose parameters all pin their type is not generic in the float type
    struct PFCountDamage <: Peridynamics.AbstractDamageModel end
    Peridynamics.@dmg_params PFCountDamage struct PFCountParams
        n_max::Int = 3
    end
    @test !(PFCountParams isa UnionAll)
    @test damage_param_type(PFCountDamage(), Float64) === PFCountParams
    @test damage_param_type(PFCountDamage(), Float32) === PFCountParams
    mp = get_dmg_params(PFCountDamage(), Float64, BBMaterial(), (;),
                        Dict{Symbol,Any}(:n_max => 5))
    @test mp === PFCountParams(5)
    @test convert_nested_params(Float32, mp) === mp
    @test damage_param_kwargs(PFCountDamage()) == (:n_max,)

    # a damage model without parameters resolves the marker to `nothing` and
    # prohibits failure by default
    struct PFPlainDamage <: Peridynamics.AbstractDamageModel end
    @test damage_param_type(PFPlainDamage(), Float64) === Nothing
    @test damage_param_kwargs(PFPlainDamage()) == ()
    pos, vol = uniform_box(1.0, 1.0, 1.0, 0.5)
    body = Body(BBMaterial(; dmgmodel=PFPlainDamage()), pos, vol)
    material!(body; horizon=1.5, rho=8e-6, E=2.1e5)
    par = only(body.point_params)
    @test par.dmg_params === nothing
    @test !any(body.fail_permit)
end
