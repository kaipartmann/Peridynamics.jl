@testitem "api tier markers: every tier has its own marker" begin
    # every docstring carries the marker of its tier, on every supported Julia version
    @test occursin("Internal use only", Peridynamics.internal_api_warning())
    @test occursin("Extension API", Peridynamics.extension_api_note())
    @test occursin("Peridynamics.<name>", Peridynamics.extension_api_note())
    @test occursin("Experimental feature", Peridynamics.experimental_api_warning())
end

@testitem "api tier markers: the tiers are not confusable" begin
    @test !occursin("Extension API", Peridynamics.internal_api_warning())
    @test !occursin("Internal use only", Peridynamics.extension_api_note())
    @test !occursin("Experimental feature", Peridynamics.extension_api_note())
end
