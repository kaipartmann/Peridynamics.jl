@testitem "internal_api_warning" begin
    if VERSION < v"1.11"
        @test !isempty(Peridynamics.internal_api_warning())
    else
        @test isempty(Peridynamics.internal_api_warning())
    end
end

@testitem "experimental_api_warning" begin
    msg = Peridynamics.experimental_api_warning()
    @test contains(msg, "Experimental feature")
    @test contains(msg, "not") && contains(msg, "public API")
end

@testitem "extension_api_note: marks the extension tier" begin
    msg = Peridynamics.extension_api_note()
    @test contains(msg, "Extension API")
    @test contains(msg, "Peridynamics.<name>")
    # the tiers must not be confusable with each other
    @test !contains(msg, "Internal use only")
    @test !contains(Peridynamics.internal_api_warning(), "Extension API")
end
