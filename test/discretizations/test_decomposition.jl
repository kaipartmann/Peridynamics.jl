@testitem "defaultdist" begin
    @test Peridynamics.defaultdist(4, 2) == [1:2, 3:4]
    @test Peridynamics.defaultdist(3, 2) == [1:2, 3:3]
    @test Peridynamics.defaultdist(2, 3) == [1:1, 2:2, 0:-1]
end

@testitem "find_localizer" begin
    point_ids = [101, 102, 103]
    localizer = Peridynamics.find_localizer(point_ids)
    @test localizer == Dict(101 => 1, 102 => 2, 103 => 3)
end
