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

    # a storage type, also instantiated, has the table of all its fields
    table = block_table(storage_type(CMaterial()))
    @test contains(table, "| `b_int` | `PointVector` | points | halo → local |")
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
