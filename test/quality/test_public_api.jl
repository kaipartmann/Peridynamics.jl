# A snapshot of the API surface of Peridynamics.jl. Its whole point is to be annoying: any
# name that is exported or declared `public` in `src/public_api.jl` shows up here, so a name
# can never be promoted or demoted by accident, only by a deliberate edit that a reviewer
# sees. Adding a name to the extension API is not breaking, removing one is.

@testmodule APISnapshot begin
    using Peridynamics

    # every name this package exports. `names` also reports the module itself.
    const EXPORTED = [
        "@mpiroot", "@mpitime", "BACMaterial", "BBMaterial", "Body", "CKIMaterial",
        "CMaterial", "CRMaterial", "CriticalStretch", "DHBBMaterial", "DynamicRelaxation",
        "EnergySurfaceCorrection", "GBBMaterial", "Job", "LinearElastic", "MultibodySetup",
        "NeoHooke", "NeoHookePenalty", "NewtonKrylov", "NoCorrection", "OSBMaterial",
        "Peridynamics", "RKCMaterial", "RKCRMaterial", "SaintVenantKirchhoff", "Study",
        "VelocityVerlet", "ZEMSilling", "ZEMWan", "const_one_kernel", "contact!",
        "cubic_b_spline_kernel", "cubic_b_spline_kernel_norm", "disable_mpi_timers!",
        "displacement_bc!", "enable_mpi_progress_bars!", "enable_mpi_timers!",
        "force_mpi_run!", "force_threads_run!", "forcedensity_bc!", "linear_kernel",
        "material!", "mpi_barrier", "mpi_isroot", "n_points", "no_failure!", "point_set!",
        "point_sets", "precrack!", "process_each_export", "process_each_job", "read_inp",
        "read_vtk", "reset_mpi_progress_bars!", "rotate!", "round_cylinder", "round_sphere",
        "submit", "submit!", "trunc_pyramid", "uniform_box", "uniform_cylinder",
        "uniform_sphere", "velocity_bc!", "velocity_ic!",
    ]

    # not exported, but declared `public` in `src/public_api.jl`: the extension API
    const PUBLIC = [
        "@cm_params", "@cm_storage", "@derived", "@dmg_params", "@dmg_storage", "@htl",
        "@inherit", "@kwarg", "@log", "@lth", "@params", "@params_fields", "@storage",
        "@storage_fields",
        "AbstractBondAssociatedSystemMaterial", "AbstractBondBasedMaterial",
        "AbstractBondSystem", "AbstractBondSystemMaterial", "AbstractConstitutiveModel",
        "AbstractConstitutiveState", "AbstractCorrespondenceMaterial", "AbstractDamageModel",
        "AbstractDamageState", "AbstractInteractionSystemMaterial", "AbstractMaterial",
        "AbstractParameterSetup", "AbstractPointParameters", "AbstractRKCMaterial",
        "AbstractStorage", "AbstractSystem", "AbstractTimeSolver",
        "BACPointParameters", "BACStorage", "BBElasticParameters", "BBPointParameters",
        "BBStorage", "Bond", "BondField", "BondFracFields", "BondHorizonParameters",
        "BondLengthCache",
        "BondScalar", "BondSymTensor", "BondSystem", "BondTensor", "BondVector",
        "CKIPointParameters",
        "CKIStorage", "CPointParameters", "CRStorage", "CStorage", "ConstitutiveParameters",
        "ConstitutiveState", "DHBBPointParameters", "DHBBStorage", "DamageParameters",
        "DamageState", "DiscretizationParameters", "DofVector", "DynamicRelaxationFields",
        "ElasticParameters", "EmptyField", "FractureParameters", "FullField", "GBBStorage",
        "HaloPoints", "HistoryDependenceError", "InteractionFracFields",
        "InteractionParameters", "InteractionSystem", "InterfaceError", "LocalPoints",
        "NewtonKrylovFields", "OSBPointParameters", "OSBStorage", "PointField",
        "PointScalar", "PointSymTensor", "PointTensor", "PointVector", "RKCFields",
        "RKCPointParameters", "RKCRStorage", "RKCStorage", "SimFloat", "StandardParameters",
        "StorageContractError", "VelocityVerletFields",
        "block_table", "bond_integrity", "bond_stretch", "calc_damage!", "calc_failure!",
        "constitutive_state", "constitutive_storage_type", "critical_stretch",
        "current_bond_length",
        "custom_field", "damage_state", "damage_storage_type", "each_bond_idx",
        "each_point_idx", "energy_release_rate", "export_field", "first_piola_kirchhoff",
        "float_type", "force_density_point!", "get_constitutive_model",
        "get_dmgmodel", "get_frac_params", "get_n_bonds", "get_n_loc_points",
        "get_n_points", "get_params", "get_sym_tensor",
        "get_tensor", "get_vector", "get_vector_diff", "has_fracture",
        "hencky_and_invstretch",
        "init_field", "is_history_dependent", "kernel", "kinematic_weight", "storage_type",
        "strain_energy_density", "supports_bond_integrity", "supports_history_dependence",
        "supports_kinematic_weight", "surface_correction_factor", "sym_eigvals",
        "update_add_vector!", "update_bond_lengths!", "update_sym_tensor!", "update_tensor!",
        "update_vector!",
    ]

    all_names() = names(Peridynamics)
    exported_names() = sort!(String.(filter(n -> Base.isexported(Peridynamics, n),
                                            all_names())))
    public_names() = sort!(String.(filter(n -> !Base.isexported(Peridynamics, n),
                                          all_names())))

    const LITERATE = normpath(@__DIR__, "..", "..", "docs", "src", "literate")
    const TUTORIALS = [joinpath(LITERATE, "tutorial_custom_material.jl"),
                       joinpath(LITERATE, "tutorial_custom_damage_model.jl"),
                       joinpath(LITERATE, "tutorial_custom_constitutive_model.jl")]

    # every name a tutorial reaches through the package: `Peridynamics.<name>` and the
    # names of `using Peridynamics: a, b` lines, which may span several lines
    function tutorial_names(path)
        src = read(path, String)
        found = Set{String}()
        for m in eachmatch(r"Peridynamics\.(@?[A-Za-z_][A-Za-z0-9_!]*)", src)
            # "Peridynamics.jl" is the package, not a name of it
            m.captures[1] == "jl" || push!(found, m.captures[1])
        end
        for m in eachmatch(r"using Peridynamics:((?:[^\n]*,\s*\n)*[^\n]*)", src)
            for name in split(m.captures[1], ",")
                isempty(strip(name)) || push!(found, strip(name))
            end
        end
        return sort!(collect(found))
    end
end

@testitem "public API: the exported names are the snapshot" tags=[:lint] setup=[APISnapshot] begin
    # the user API. A change here is breaking and belongs in `NEWS.md`.
    @test APISnapshot.exported_names() == APISnapshot.EXPORTED
end

@testitem "public API: the public names are the snapshot" tags=[:lint] setup=[APISnapshot] begin
    if VERSION >= v"1.11"
        # `public` only exists from 1.11 on, and `src/public_api.jl` is only included there
        @test APISnapshot.public_names() == APISnapshot.PUBLIC
    else
        # on the LTS the file is never parsed, so nothing but the exports is visible
        @test isempty(APISnapshot.public_names())
    end
end

@testitem "public API: every public name is defined and documented" tags=[:lint] setup=[APISnapshot] begin
    meta = Base.Docs.meta(Peridynamics)
    for name in Symbol.(APISnapshot.PUBLIC)
        # a `public` declaration for a name that does not exist is silently accepted by
        # Julia, so check it here rather than finding out in someone else's package
        @test isdefined(Peridynamics, name)
        # the extension API is the documented API: `checkdocs = :public` in `docs/make.jl`
        # then also enforces that each of these appears on a reference page
        @test haskey(meta, Base.Docs.Binding(Peridynamics, name))
    end
end

@testitem "public API: the tiers do not overlap" tags=[:lint] setup=[APISnapshot] begin
    @test isempty(intersect(APISnapshot.EXPORTED, APISnapshot.PUBLIC))
    # a public name is reachable as `Peridynamics.<name>` but must not leak into scope
    if VERSION >= v"1.11"
        @test Base.ispublic(Peridynamics, :storage_type)
        @test !Base.isexported(Peridynamics, :storage_type)
        @test !isdefined(@__MODULE__, :storage_type)
    end
end

@testitem "public API: the tutorials use only exported or public names" tags=[:lint] setup=[APISnapshot] begin
    # the tutorials are the promise that an extension needs nothing internal, so every
    # `Peridynamics.<name>` in them has to be in one of the two tiers. `LinearAlgebra` and
    # `StaticArrays` are reached through the package and are not names of it.
    reexported = ("LinearAlgebra", "StaticArrays")
    tiers = union(Set(APISnapshot.EXPORTED), Set(APISnapshot.PUBLIC))
    for path in APISnapshot.TUTORIALS
        found = filter(n -> !(n in reexported), APISnapshot.tutorial_names(path))
        @test !isempty(found)
        internal = filter(n -> !(n in tiers), found)
        @test isempty(internal)
    end
end
