#=
Every abstract type of the package, declared in one place so that the type hierarchy can be
read at a glance. The types of the extension API carry `extension_api_note()`, everything
else is internal. Supertypes are declared before their subtypes, so the order of the
sections matters.
=#

# --------------------------------------------------------------------------------------
# materials, models and their families
# --------------------------------------------------------------------------------------

"""
    AbstractMaterial

$(extension_api_note())

Supertype of every material. A new material is usually not a direct subtype of this one but
of the family that matches the system it is discretized on, which is
[`AbstractBondSystemMaterial`](@ref) for bonds and [`AbstractInteractionSystemMaterial`](@ref)
for interactions.

Every material declares its point parameters with [`@params`](@ref), its storage with
[`@storage`](@ref) and defines [`force_density_point!`](@ref). The tutorial on custom
materials walks through all three.
"""
abstract type AbstractMaterial end

"""
    AbstractBondSystemMaterial{Correction}

$(extension_api_note())

Supertype of every material that is discretized on a [`BondSystem`](@ref), that is on the
bonds between a point and its neighbors. `Correction` is the surface correction of the
material, e.g. `NoCorrection` or `EnergySurfaceCorrection`.

A subtype defines [`force_density_point!`](@ref) with a [`BondSystem`](@ref) as the second
argument, and it needs a field `dmgmodel`, because every bond system material is asked for
its damage model before the force density is evaluated.
"""
abstract type AbstractBondSystemMaterial{Correction} <: AbstractMaterial end

"""
    AbstractBondBasedMaterial{Correction}

$(extension_api_note())

Supertype of the bond-based materials, in which the force density of a bond depends only on
that bond. `Correction` is the surface correction of the material.

This is a subtype of [`AbstractBondSystemMaterial`](@ref). [`BBMaterial`](@ref),
[`DHBBMaterial`](@ref) and [`GBBMaterial`](@ref) are the materials of this family.
"""
abstract type AbstractBondBasedMaterial{CM} <: AbstractBondSystemMaterial{CM} end

"""
    AbstractCorrespondenceMaterial{CM,ZEM}

$(extension_api_note())

Supertype of the correspondence materials that evaluate a constitutive model once per point.
`CM` is the constitutive model and `ZEM` the zero-energy mode stabilization.

The constitutive model is reached with [`get_constitutive_model`](@ref), so a model written
against [`first_piola_kirchhoff`](@ref) works with every material of this family.
[`CMaterial`](@ref) and [`CRMaterial`](@ref) are the materials of this family.
"""
abstract type AbstractCorrespondenceMaterial{CM,ZEM} <: AbstractBondSystemMaterial{ZEM} end

"""
    AbstractRKCMaterial{CM,Correction,Monomial}

$(extension_api_note())

Supertype of the reproducing kernel correspondence materials, which evaluate a constitutive
model at a bond-associated quadrature point. `CM` is the constitutive model, `Correction`
the surface correction, and `Monomial` the symbol of the monomial basis of the reproducing
kernel (`:C1`, `:RK1`, `:RK2` or `:PD2`).

The index passed to [`first_piola_kirchhoff`](@ref) is a bond index for this family. The
materials of this family honor [`bond_integrity`](@ref) and [`kinematic_weight`](@ref), so a
damage model that softens bonds works with them. [`RKCMaterial`](@ref) and
[`RKCRMaterial`](@ref) are the materials of this family.
"""
abstract type AbstractRKCMaterial{CM,C,M} <: AbstractBondSystemMaterial{C} end

"""
    AbstractBondAssociatedSystemMaterial

$(extension_api_note())

Supertype of the materials that are discretized on a bond-associated system, which evaluate
their constitutive model on the bond-associated neighborhood of every bond.
[`BACMaterial`](@ref) is the material of this family, and the index passed to
[`first_piola_kirchhoff`](@ref) is a bond index.
"""
abstract type AbstractBondAssociatedSystemMaterial <: AbstractBondSystemMaterial{Nothing} end

"""
    AbstractInteractionSystemMaterial

$(extension_api_note())

Supertype of every material that is discretized on an [`InteractionSystem`](@ref), that is
on one-, two- and three-neighbor interactions instead of bonds. [`CKIMaterial`](@ref) is the
material of this family.

A subtype defines [`force_density_point!`](@ref) with an [`InteractionSystem`](@ref) as the
second argument.
"""
abstract type AbstractInteractionSystemMaterial <: AbstractMaterial end

"""
    AbstractConstitutiveModel

$(extension_api_note())

Supertype of every constitutive model. A constitutive model turns a deformation gradient
into a stress. It does not depend on the material family that evaluates it, so the same
model runs on [`CMaterial`](@ref), [`RKCMaterial`](@ref) and [`BACMaterial`](@ref) alike.

A model defines [`first_piola_kirchhoff`](@ref). A model with parameters of its own declares
them with [`@cm_params`](@ref), and a model that carries a history declares its state with
[`@cm_storage`](@ref), which makes [`is_history_dependent`](@ref) true.
"""
abstract type AbstractConstitutiveModel end

"""
    AbstractConstitutiveState

$(extension_api_note())

Supertype of the state of a history-dependent constitutive model. Types of this kind are
generated by [`@cm_storage`](@ref) and are never written by hand.

The state that a storage carries is reached with [`constitutive_state`](@ref).
"""
abstract type AbstractConstitutiveState end

"""
    AbstractDamageModel

$(extension_api_note())

Supertype of every damage model. A damage model decides when a bond or an interaction fails.

A subtype defines [`calc_failure!`](@ref). A model that reads the fracture keywords of
[`material!`](@ref) declares its parameters with [`@dmg_params`](@ref), converts them with
[`get_frac_params`](@ref) and says with [`has_fracture`](@ref) whether fracture is enabled.
[`CriticalStretch`](@ref) is the damage model of this package.

A model that needs per-bond variables of its own declares them with
[`@dmg_storage`](@ref) and reaches them with [`damage_state`](@ref). It then works with every
material that carries a `dmg_state::DamageState` field, without either knowing the other. A
model that softens a bond instead of deleting it also defines [`bond_integrity`](@ref),
[`kinematic_weight`](@ref) and [`calc_damage!`](@ref).
"""
abstract type AbstractDamageModel end

"""
    AbstractDamageState

$(extension_api_note())

Supertype of the state of a damage model, that is of the per-bond or per-point variables it
integrates, such as an accumulated damage. Types of this kind are generated by
[`@dmg_storage`](@ref) and are never written by hand.

The state that a storage carries is reached with [`damage_state`](@ref). Unlike a
constitutive state, a damage state does not make anything history dependent. A damage model
advances its state in [`calc_failure!`](@ref), which every time solver calls exactly once per
step.
"""
abstract type AbstractDamageState end

"""
    AbstractConstitutiveParameters

$(internal_api_warning())

Supertype of the parameters a constitutive model owns. Types of this kind are generated by
[`@cm_params`](@ref) and fill the `cm_params::ConstitutiveParameters` marker field of a
point parameter type.
"""
abstract type AbstractConstitutiveParameters end

"""
    AbstractDamageParameters

$(internal_api_warning())

Supertype of the parameters a damage model owns. Types of this kind are generated by
[`@dmg_params`](@ref) and fill the `dmg_params::DamageParameters` marker field of a point
parameter type.
"""
abstract type AbstractDamageParameters end

"""
    AbstractCorrection

$(internal_api_warning())

Supertype of the surface corrections of a bond system, e.g. `NoCorrection` and
`EnergySurfaceCorrection`, and of the zero-energy mode stabilizations.
"""
abstract type AbstractCorrection end

"""
    AbstractZEMStabilization

$(internal_api_warning())

Supertype of the zero-energy mode stabilizations of the correspondence materials, e.g.
`ZEMSilling` and `ZEMWan`.
"""
abstract type AbstractZEMStabilization <: AbstractCorrection end

# --------------------------------------------------------------------------------------
# point parameters and storages
# --------------------------------------------------------------------------------------

"""
    AbstractParameterSetup

$(extension_api_note())

Supertype of what a body chunk carries as its point parameters: either one set of parameters
for the whole chunk, or a handler that resolves them per point when [`material!`](@ref) was
called more than once. Either way, [`get_params`](@ref) reads the parameters of a point from
it. It is what [`force_density_point!`](@ref) receives as its parameters argument, and the
type to annotate that argument with when a method has to, e.g. when [`calc_failure!`](@ref)
is defined for one damage model.
"""
abstract type AbstractParameterSetup end

"""
    AbstractPointParameters

$(extension_api_note())

Supertype of every set of point parameters. Point parameters are the material properties of
a single point. They are `isbits`, so that they can be captured by value.

Sets of point parameters are generated by [`@params`](@ref) and are never written by hand.
The parameters of point `i` are reached with [`get_params`](@ref).
"""
abstract type AbstractPointParameters <: AbstractParameterSetup end

"""
    AbstractParameterHandler

$(internal_api_warning())

Supertype of the handlers that resolve the point parameters per point when a body has more
than one parameter set, see `ParameterHandler`.
"""
abstract type AbstractParameterHandler <: AbstractParameterSetup end

"""
    AbstractPointParameterFields

$(internal_api_warning())

Supertype of the parameter blocks defined with [`@params_fields`](@ref). A block is a marker
type that carries its declarations, so that [`@inherit`](@ref) can include them.
"""
abstract type AbstractPointParameterFields end

"""
    AbstractParamSpec

$(internal_api_warning())

Supertype of the specifications that say how a body chunk carries its point parameters,
see `SingleParamChunk` and `MultiParamChunk`.
"""
abstract type AbstractParamSpec end

"""
    AbstractStorage

$(extension_api_note())

Supertype of every storage. A storage holds all point and bond fields of a body chunk, one
array per quantity, and it is the only thing that changes during a simulation.

Storages are generated by [`@storage`](@ref) and are never written by hand. The storage type
that belongs to a material is [`storage_type`](@ref).
"""
abstract type AbstractStorage end

"""
    AbstractStorageFields

$(internal_api_warning())

Supertype of the field blocks defined with [`@storage_fields`](@ref). A block is a marker
type that carries its declarations, so that [`@inherit`](@ref) can include them.
"""
abstract type AbstractStorageFields end

"""
    AbstractFieldShape{T}

$(internal_api_warning())

Supertype of every field shape with element type `T`, e.g. [`PointVector`](@ref). A shape
says what a storage field means, and the container type and the allocation of the field
follow from it.
"""
abstract type AbstractFieldShape{T} end

"""
    AbstractSolverField

$(internal_api_warning())

Supertype of the answers a time solver gives when asked whether it needs a storage field,
see [`FullField`](@ref) and [`EmptyField`](@ref).
"""
abstract type AbstractSolverField end

# --------------------------------------------------------------------------------------
# discretization
# --------------------------------------------------------------------------------------

"""
    AbstractSpatialSetup

$(internal_api_warning())

Supertype of everything a [`Job`](@ref) can be created for, that is a [`Body`](@ref) or a
[`MultibodySetup`](@ref).
"""
abstract type AbstractSpatialSetup end

"""
    AbstractBody{T<:AbstractMaterial}

$(internal_api_warning())

Supertype of a body with the material `T`, see [`Body`](@ref).
"""
abstract type AbstractBody{T<:AbstractMaterial} <: AbstractSpatialSetup end

"""
    AbstractMultibodySetup

$(internal_api_warning())

Supertype of a setup of several bodies that interact through contact, see
[`MultibodySetup`](@ref).
"""
abstract type AbstractMultibodySetup <: AbstractSpatialSetup end

"""
    AbstractPredefinedCrack

$(internal_api_warning())

Supertype of the cracks that are cut into a body before the simulation starts, see
[`precrack!`](@ref).
"""
abstract type AbstractPredefinedCrack end

"""
    AbstractCondition

$(internal_api_warning())

Supertype of the boundary and initial conditions of a body, e.g. what
[`velocity_bc!`](@ref) and [`velocity_ic!`](@ref) add.
"""
abstract type AbstractCondition end

"""
    AbstractSystem

$(extension_api_note())

Supertype of every system. A system is the discretization of a body chunk: the points, their
volumes and their neighborhood relations. It is created once during setup and never changes.

[`BondSystem`](@ref) and [`InteractionSystem`](@ref) are the systems of this package.
"""
abstract type AbstractSystem end

"""
    AbstractBondSystem

$(extension_api_note())

Supertype of the systems whose neighborhood relation is a bond, that is
[`BondSystem`](@ref) and the bond-associated system.

The points of such a system are iterated with [`each_point_idx`](@ref) and the bonds of
point `i` with [`each_bond_idx`](@ref).
"""
abstract type AbstractBondSystem <: AbstractSystem end

"""
    AbstractChunkHandler

$(internal_api_warning())

Supertype of the handlers that know which points of a body a chunk owns and which ones it
reads from other chunks, see `ChunkHandler`.
"""
abstract type AbstractChunkHandler end

"""
    AbstractBodyChunk{S<:AbstractSystem,T<:AbstractMaterial}

$(internal_api_warning())

Supertype of a chunk of a body: the part of a body one thread or one MPI rank is
responsible for, with its system `S`, its material `T`, its storage and its point
parameters. See `BodyChunk`.
"""
abstract type AbstractBodyChunk{S<:AbstractSystem,T<:AbstractMaterial} end

# --------------------------------------------------------------------------------------
# running a simulation
# --------------------------------------------------------------------------------------

"""
    AbstractTimeSolver

$(extension_api_note())

Supertype of every time solver. A time solver decides how the simulation advances in time
and which fields a storage needs, e.g. through [`VelocityVerletFields`](@ref).

[`VelocityVerlet`](@ref), [`DynamicRelaxation`](@ref) and [`NewtonKrylov`](@ref) are the
solvers of this package. Writing a new one is an internal interface at the moment and is not
covered by the extension API.
"""
abstract type AbstractTimeSolver end

"""
    AbstractJob

$(internal_api_warning())

Supertype of a simulation job, see [`Job`](@ref).
"""
abstract type AbstractJob end

"""
    AbstractJobOptions

$(internal_api_warning())

Supertype of the options of a job, e.g. the export path and the exported fields, see
`JobOptions`.
"""
abstract type AbstractJobOptions end

"""
    AbstractDataHandler

$(internal_api_warning())

Supertype of the data handlers, which own the body chunks of a simulation and run the halo
exchange between them. There is one handler per parallelization backend and per kind of
spatial setup.
"""
abstract type AbstractDataHandler end

"""
    AbstractThreadsDataHandler

$(internal_api_warning())

Supertype of the data handlers of a multithreaded simulation, see
`ThreadsBodyDataHandler` and `ThreadsMultibodyDataHandler`.
"""
abstract type AbstractThreadsDataHandler <: AbstractDataHandler end

"""
    AbstractMPIDataHandler

$(internal_api_warning())

Supertype of the data handlers of an MPI simulation, see `MPIBodyDataHandler`.
"""
abstract type AbstractMPIDataHandler <: AbstractDataHandler end

"""
    AbstractThreadsBodyDataHandler{Sys,M,P,S}

$(internal_api_warning())

Supertype of the multithreaded data handler of a single body with the system `Sys`, the
material `M`, the parameter setup `P` and the storage `S`.
"""
abstract type AbstractThreadsBodyDataHandler{Sys,M,P,S} <: AbstractThreadsDataHandler end

"""
    AbstractThreadsMultibodyDataHandler

$(internal_api_warning())

Supertype of the multithreaded data handler of a [`MultibodySetup`](@ref).
"""
abstract type AbstractThreadsMultibodyDataHandler <: AbstractThreadsDataHandler end

"""
    AbstractMPIBodyDataHandler{Sys,M,P,S}

$(internal_api_warning())

Supertype of the MPI data handler of a single body with the system `Sys`, the material `M`,
the parameter setup `P` and the storage `S`.
"""
abstract type AbstractMPIBodyDataHandler{Sys,M,P,S} <: AbstractMPIDataHandler end

"""
    AbstractMPIMultibodyDataHandler

$(internal_api_warning())

Supertype of the MPI data handler of a [`MultibodySetup`](@ref).
"""
abstract type AbstractMPIMultibodyDataHandler <: AbstractMPIDataHandler end
