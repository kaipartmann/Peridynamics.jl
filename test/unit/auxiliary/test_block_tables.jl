@testitem "block_table: the table of a field block and of a storage" begin
    import Peridynamics: block_table, StorageFieldDecl, ConstitutiveState, VelocityVerletFields,
                         NewtonKrylovFields, storage_type

    # a field block: one row per field with its shape, its entries and its halo exchange
    table = block_table(VelocityVerletFields)
    @test startswith(table, "| field | shape | entries | halo exchange |\n|:---|:---|:---|:---|\n")
    @test contains(table, "| `position` | `PointVector{Float64}` | points | local → halo |")
    @test contains(table, "| `b_int` | `PointVector` | points | – |")
    @test !contains(table, "Peridynamics.")
    @test contains(block_table(NewtonKrylovFields), "| `residual` | `DofVector` | degrees of freedom | – |")

    # a storage type, also instantiated, has the table of all its fields; the fracture
    # bookkeeping is in the state of the damage model and has a table of its own
    table = block_table(storage_type(CMaterial()))
    @test contains(table, "| `b_int` | `PointVector` | points | halo → local |")
    @test contains(table, "| `dmg_state` | `DamageState` | state of the damage model | – |")
    table = block_table(Peridynamics.BondFracState)
    @test contains(table, "| `bond_active` | `BondScalar{Bool}` | bonds | – |")
    table = block_table(storage_type(BACMaterial()))
    @test contains(table, "| `bond_stress` | `Matrix{Float64}` | – | – |")

    # the nested state of the constitutive model is no array
    decls = [StorageFieldDecl(:cm_state, :none, ConstitutiveState, nothing, nothing)]
    @test contains(block_table(decls), "| `cm_state` | `ConstitutiveState` | state of the constitutive model | – |")

    # only what `@storage` and `@storage_fields` define has a table
    err = try
        block_table(Int)
    catch e
        e
    end
    @test err isa ArgumentError
    @test contains(err.msg, "does not declare any parameters or storage fields")
    @test_throws ArgumentError block_table(BBMaterial)
end

@testitem "show: a field block prints its table at the REPL and in the docs" begin
    import Peridynamics: VelocityVerletFields

    plain = sprint(show, MIME"text/plain"(), VelocityVerletFields)
    @test contains(plain, "VelocityVerletFields")
    @test contains(plain, "@inherit")
    @test contains(plain, "position")
    @test contains(plain, "local → halo")
    html = sprint(show, MIME"text/html"(), VelocityVerletFields)
    @test contains(html, "<table>")
    @test contains(html, "<code>position</code>")

    # a block that registers no fields falls back to its name
    struct LonelyBlock <: Peridynamics.AbstractStorageFields end
    plain = sprint(show, MIME"text/plain"(), LonelyBlock)
    @test contains(plain, "LonelyBlock")
    @test !contains(plain, "@inherit")
    @test !contains(sprint(show, MIME"text/html"(), LonelyBlock), "<table>")
end

@testitem "block_table: the table of a parameter block and of point parameters" begin
    import Peridynamics: block_table, DiscretizationParameters, point_param_type

    params = block_table(DiscretizationParameters)
    @test startswith(params,
                     "| parameter | type | `material!` keyword | value | simulation log |")
    @test occursin("| `δ` |", params)
    @test occursin("| `rho` |", params)
    # the keywords appear inside the call that reads them, which is the connection between
    # what `material!` takes and what the parameters are
    @test occursin("`(; δ, rho) = get_discretization_params(; horizon, rho)`", params)
    @test occursin("Keywords of `material!`: `horizon`, `rho`.", params)
    @test occursin("horizon", params)   # the simulation log label

    # a material instance answers with the table of its point parameters
    @test block_table(BBMaterial()) == block_table(point_param_type(BBMaterial()))
    @test occursin("Keywords of `material!`: `horizon`, `rho`, `E`",
                   block_table(BBMaterial()))

    # typing the name of a parameter block at the REPL shows the same table, rendered
    shown = sprint(show, MIME"text/plain"(), DiscretizationParameters)
    @test occursin("a block you can @inherit", shown)
    @test occursin("get_discretization_params(; horizon, rho)", shown)
    @test occursin("Keywords of material!: horizon, rho.", shown)
    # and a real table in the documentation
    @test occursin("<table>", sprint(show, MIME"text/html"(), DiscretizationParameters))
end
