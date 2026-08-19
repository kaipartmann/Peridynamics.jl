@testitem "internal_api_warning" begin
    if VERSION < v"1.11"
        @test !isempty(Peridynamics.internal_api_warning())
    else
        @test isempty(Peridynamics.internal_api_warning())
    end
end
