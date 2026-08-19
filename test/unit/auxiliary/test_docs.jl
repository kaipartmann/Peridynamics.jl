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
