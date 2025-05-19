@testitem "required_fields" begin
    rf_abstractmat = (:position, :displacement, :velocity, :velocity_half, :acceleration,
                      :b_int, :b_ext, :velocity_half_old, :b_int_old, :density_matrix)
    @test Peridynamics.required_fields(Peridynamics.AbstractMaterial) === rf_abstractmat

    rf_bb = (:position, :displacement, :velocity, :velocity_half, :acceleration, :b_int,
             :b_ext, :velocity_half_old, :b_int_old, :density_matrix, :damage,
             :n_active_bonds, :bond_active)
    @test Peridynamics.required_fields(BBMaterial) === rf_bb

    rf_cki = (:position, :displacement, :velocity, :velocity_half, :acceleration, :b_int,
              :b_ext, :velocity_half_old, :b_int_old, :density_matrix, :damage,
              :n_active_one_nis, :one_ni_active)
    @test Peridynamics.required_fields(CKIMaterial) === rf_cki
end

@testitem "get_storage_header" begin
    import Peridynamics: get_storage_header

    # Test case 1: Simple struct declaration
    expr1 = :(struct MyStorage end)
    header1, type1 = get_storage_header(expr1)
    @test header1 == Expr(:(<:), :MyStorage, :(Peridynamics.AbstractStorage))
    @test type1 == :MyStorage

    # Test case 2: Struct declaration with a subtype
    expr2 = :(struct MyStorage <: Peridynamics.AbstractStorage end)
    header2, type2 = get_storage_header(expr2)
    @test header2 == Expr(:(<:), :MyStorage, :(Peridynamics.AbstractStorage))
    @test type2 == :MyStorage

    # Test case 3: Parametric struct declaration with a subtype
    expr3 = :(struct MyStorage{A,B,C} <: Peridynamics.AbstractStorage end)
    header3, type3 = get_storage_header(expr3)
    @test header3 == Expr(:(<:), :(MyStorage{A,B,C}), :(Peridynamics.AbstractStorage))
    @test type3 == :MyStorage

    # Test case 4: Parametric struct declaration without a subtype
    expr4 = :(struct MyStorage{A,B,C} end)
    header4, type4 = get_storage_header(expr4)
    @test header4 == Expr(:(<:), :(MyStorage{A,B,C}), :(Peridynamics.AbstractStorage))
    @test type4 == :MyStorage
end

@testitem "macrochecks @storage macro" begin

    input = :(Peridynamics.BBStorage)
    @test isnothing(Peridynamics.macrocheck_input_storage_type(input))

    input = :(MyNonexistingStorage)
    @test isnothing(Peridynamics.macrocheck_input_storage_type(input))

    input = :(Base.MyNonexistingStorage)
    @test isnothing(Peridynamics.macrocheck_input_storage_type(input))

    input = :(:MyNonexistingStorageAsSymbol)
    @test_throws ArgumentError Peridynamics.macrocheck_input_storage_type(input)

    input = :(
        struct MyTestStorage <: Peridynamics.AbstractStorage
            a::Int
        end
    )
    @test isnothing(Peridynamics.macrocheck_input_storage_struct(input))

    input = :(
        struct MyTestStorage{A,B,C} <: Peridynamics.AbstractStorage
            c::Vector{Float64}
        end
    )
    @test isnothing(Peridynamics.macrocheck_input_storage_struct(input))

    input = :(
        @kwdef struct MyTestStorage <: Peridynamics.AbstractStorage
            c::Vector{Float64}
        end
    )
    @test_throws ArgumentError Peridynamics.macrocheck_input_storage_struct(input)

    input = :(Peridynamics.VelocityVerlet)
    @test isnothing(Peridynamics.macrocheck_input_timesolver(input))

    input = :(MyNonexistingSolver)
    @test isnothing(Peridynamics.macrocheck_input_timesolver(input))

    input = :(Base.MyNonexistingSolver)
    @test isnothing(Peridynamics.macrocheck_input_timesolver(input))

    input = :(:MyNonexistingSolverAsSymbol)
    @test_throws ArgumentError Peridynamics.macrocheck_input_timesolver(input)

    input = :(:testfield)
    @test isnothing(Peridynamics.macrocheck_input_field(input))

    input = :testfield
    @test_throws ArgumentError Peridynamics.macrocheck_input_field(input)
end
