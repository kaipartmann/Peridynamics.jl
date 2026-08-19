# Material tables. Every material-generic test loops over one of these, so a new material only
# needs a line here to be covered by all of them.

"One entry per material configuration, with its default constitutive model where it has one."
const MATERIALS = [
    "BB" => BBMaterial(),
    "BB-ESC" => BBMaterial{EnergySurfaceCorrection}(),
    "DHBB" => DHBBMaterial(),
    "GBB" => GBBMaterial(),
    "OSB" => OSBMaterial(),
    "OSB-ESC" => OSBMaterial{EnergySurfaceCorrection}(),
    "C" => CMaterial(),
    "C-ZEMWan" => CMaterial(zem=ZEMWan()),
    "CR" => CRMaterial(),
    "BAC" => BACMaterial(),
    "RKC-C1" => RKCMaterial(monomial=:C1),
    "RKC-RK1" => RKCMaterial(monomial=:RK1),
    "RKCR" => RKCRMaterial(),
    "CKI" => CKIMaterial(),
]

"""
The materials with a storage and force density implementation of their own. Every other entry
of `MATERIALS` is a variant (a correction, a stabilization, a kernel) that shares the solver
interaction with one of these. Per-commit solver runs loop over these; the variants run in the
`extras` job (`:slow`), see `test/simulation/`.
"""
const CORE_MATERIALS = filter(p -> first(p) in ("BB", "OSB", "C", "CR", "BAC", "RKC-C1", "RKCR", "CKI"),
                              MATERIALS)
const VARIANT_MATERIALS = filter(p -> !(p in CORE_MATERIALS), MATERIALS)

"The constitutive models and the materials that accept each of them."
const CONSTITUTIVE_MODELS = [
    "SVK" => SaintVenantKirchhoff(),
    "LE" => LinearElastic(),
    "NH" => NeoHooke(),
    "NHP" => NeoHookePenalty(),
]

"""
Every (material, constitutive model) combination that the package supports. The rotated
formulations only accept the Saint-Venant-Kirchhoff and the linear elastic model.
"""
const MODEL_MATERIALS = let pairs = Pair{String,Any}[]
    for (mname, model) in CONSTITUTIVE_MODELS
        push!(pairs, "C-$mname" => CMaterial(model=model))
        push!(pairs, "RKC-$mname" => RKCMaterial(model=model))
        push!(pairs, "BAC-$mname" => BACMaterial(model=model))
        if model isa Union{SaintVenantKirchhoff,LinearElastic}
            push!(pairs, "CR-$mname" => CRMaterial(model=model))
            push!(pairs, "RKCR-$mname" => RKCRMaterial(model=model))
        end
    end
    pairs
end

"Materials that average the deformation gradient over the family (zero-energy mode relevant)."
is_correspondence(mat) = mat isa Union{Peridynamics.AbstractCorrespondenceMaterial,
                                      Peridynamics.AbstractRKCMaterial,
                                      Peridynamics.AbstractBondAssociatedSystemMaterial}
const CORRESPONDENCE_MATERIALS = filter(p -> is_correspondence(last(p)), MATERIALS)

"""
Horizon-to-spacing ratio used for `mat`. The interaction system materials evaluate triplets of
neighbors, and the triplet search scales with the ninth power of this ratio, so they get a
smaller one to keep the runtime in the same order of magnitude as the bond systems.
"""
horizon_ratio(::Peridynamics.AbstractMaterial) = 3.015
horizon_ratio(::Peridynamics.AbstractInteractionSystemMaterial) = 2.015

"""
    supports(solver, mat)

Whether `solver` can run a body of material `mat`. `NewtonKrylov` rejects the rotated
formulations, which update a stress history per time step, see
`validate_chunk_for_newton_krylov!`. Loops over solvers × materials filter with this.
"""
supports(::Peridynamics.AbstractTimeSolver, mat) = true
supports(::NewtonKrylov, mat) = !(mat isa Union{CRMaterial,RKCRMaterial})

"""
    cki_kwargs()

Explicit interaction parameters for `CKIMaterial`. With only `E` and `nu` given, `C3` is derived
as `32/π⁴ (λ - μ)/δ¹²`, which is exactly zero for `nu = 0.25`, so no three-neighbor interactions
are built at all and the corresponding force loops never run. Pass these to `material!` when a
test must reach the two- and three-neighbor code paths. Setting them manually emits a warning
that the test has to expect with `@test_logs (:warn, r"specified manually")`.
"""
cki_kwargs() = (; C1=1e11, C2=1e11, C3=1e11)
