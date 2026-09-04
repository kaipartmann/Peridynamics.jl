# =========================================================================================
#  THE EXTENSION API
# =========================================================================================
#
#  Peridynamics.jl has three API tiers. This file declares the second one.
#
#    1. The user API. Everything `export`ed in `src/Peridynamics.jl`, which is what a
#       simulation script is written with. Stable. A change is breaking and is listed in
#       `NEWS.md`.
#
#    2. The extension API. The names declared `public` below, which is what you need to add
#       a material, a constitutive model, a damage model or a set of point parameters. Not
#       exported, so it is always written as `Peridynamics.foo` or pulled in explicitly
#       with `using Peridynamics: foo`. Stable within a minor release. A rename is breaking
#       and is listed in `NEWS.md`. While the package is at 0.x this tier may still change
#       in a minor version bump, see `docs/src/api_stability.md`.
#
#    3. Everything else is internal and may change in any release without notice, even
#       though much of it is documented.
#
#  A name belongs in tier 2 if the custom material tutorial or the manual names it, and if
#  its spelling and signature are settled. Deliberately not declared here, so that they stay
#  free to change: the family-level stress hooks (`calc_first_piola_kirchhoff!`,
#  `rkc_stress_integral!`, `monomial`, ...), everything that `@storage` and `@params`
#  generate for you (`allowed_material_kwargs`, `req_storage_fields`, `point_param_type`,
#  ...), all macro expansion machinery, the whole parallelization layer, and the solver and
#  system extension paths.
#
#  Adding a name here is not breaking, removing one is. Start narrow, grow on demand.
#  `test/quality/test_public_api.jl` holds a snapshot of this list, so nothing is promoted
#  or demoted by accident.
# =========================================================================================

# The declaration macros: the storage layer, the state of a constitutive or a damage model,
# the point parameter layer, the parameters a model owns, and the annotations that are
# recognized inside their bodies.
public @storage, @storage_fields, @cm_storage, @dmg_storage, @inherit, @lth, @htl
public @params, @params_fields, @cm_params, @dmg_params, @kwarg, @derived, @log

# Field shapes. What a storage field is declared with. They determine its size, its element
# type and how it is allocated.
public PointScalar, PointVector, PointTensor, PointSymTensor, PointField
public BondScalar, BondVector, BondTensor, BondSymTensor, BondField
public DofVector

# Field shape modifiers and the markers of a storage and a point parameter definition.
public SimFloat, LocalPoints, HaloPoints, FullField, EmptyField
public ConstitutiveState, DamageState
public ConstitutiveParameters, DamageParameters

# Storage field blocks, to be `@inherit`ed by a storage or a damage state. Every storage
# needs the block of the time solver it is used with. The fracture bookkeeping blocks are
# inherited inside a `@dmg_storage` declaration, the others follow from the system and the
# material family.
public VelocityVerletFields, DynamicRelaxationFields, NewtonKrylovFields
public BondLengthCache, BondFracFields, InteractionFracFields, RKCFields
public BondFracState, InteractionFracState

# Point parameter blocks, to be `@inherit`ed by a set of point parameters.
public DiscretizationParameters, ElasticParameters, BBElasticParameters
public FractureParameters, BondHorizonParameters, InteractionParameters
public StandardParameters

# The point parameters and the storages of the shipped materials. A material of your own
# inherits from them, or uses the point parameters as they are with `@params MyMat
# BBPointParameters`. What they expose is generated from their declarations.
public BBPointParameters, DHBBPointParameters, OSBPointParameters, CPointParameters
public RKCPointParameters, BACPointParameters, CKIPointParameters
public BBStorage, DHBBStorage, GBBStorage, OSBStorage, CStorage, CRStorage
public RKCStorage, RKCRStorage, BACStorage, CKIStorage

# What a block exposes, generated from its declarations.
public block_table

# The abstract types a new material, model or storage is a subtype of.
public AbstractMaterial, AbstractBondSystemMaterial, AbstractBondBasedMaterial
public AbstractCorrespondenceMaterial, AbstractRKCMaterial
public AbstractBondAssociatedSystemMaterial, AbstractInteractionSystemMaterial
public AbstractConstitutiveModel, AbstractConstitutiveState
public AbstractDamageModel, AbstractDamageState
public AbstractStorage, AbstractPointParameters, AbstractParameterSetup
public AbstractSystem, AbstractBondSystem, AbstractTimeSolver

# The systems a material is dispatched on, and what it takes to write one: the declaration
# macros, the sizes a constructor allocates against and the two functions that say which
# system a material gets and which materials a system accepts.
public BondSystem, InteractionSystem
public @system
public SystemSizes, system_type, host_system_type, check_system_compat
public max_n_chunks, first_chunk

# The errors the interfaces throw. Catch them in tests, or throw them from your own
# interface.
public InterfaceError, StorageContractError, HistoryDependenceError

# The material interface.
public storage_type, init_field, force_density_point!

# The constitutive model interface.
public get_constitutive_model, first_piola_kirchhoff, strain_energy_density
public constitutive_state, constitutive_storage_type
public is_history_dependent, supports_history_dependence

# The damage model interface. The relation between `Gc` and `εc` depends on the
# micro-modulus, so a material with a non-constant one defines the two hooks.
public get_dmgmodel, get_frac_params, has_fracture, calc_failure!, calc_damage!
public damage_state, damage_storage_type
public critical_stretch, energy_release_rate
# The damage model is a plug-in box: everything outside of it reads the fracture
# bookkeeping through these functions and never by field name, and a material that kills
# bonds writes through them.
public bond_is_active, get_damage, break_bond!, break_bonds!
# A damage model may soften a bond instead of deleting it. `bond_integrity` scales the load
# a bond still carries, `kinematic_weight` scales what it contributes to the deformation
# gradient. Both default to one. A material says whether its force path honors them.
public bond_integrity, kinematic_weight
public supports_bond_integrity, supports_kinematic_weight

# Accessing a system, its points and its bonds from inside a force density calculation.
public get_params, each_point_idx, each_bond_idx
public get_n_points, get_n_loc_points, get_n_bonds, get_n_dim
public kernel, surface_correction_factor, float_type
# The kinematics of a bond. Read the current length and the stretch of a bond with these and
# never by gathering the two positions, then a material with a bond length cache and one
# without it both run at full speed.
public current_bond_length, bond_stretch, update_bond_lengths!
public get_neighbor, reference_bond_length, bond_may_fail

# Reading and writing the columns of a storage field as static vectors and tensors.
public dims, get_vector, get_vector_diff, update_vector!, update_add_vector!
public get_tensor, update_tensor!, get_sym_tensor, update_sym_tensor!

# Strain measures a finite strain constitutive model needs, evaluated in closed form.
public sym_eigvals, hencky_and_invstretch

# Exporting a field of your own storage to VTK.
public export_field, custom_field
