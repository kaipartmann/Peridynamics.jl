@testitem "get_submit_options" begin
    o = Dict{Symbol,Any}(:quiet => true)
    quiet = Peridynamics.get_submit_options(o)
    @test quiet == true

    o = Dict{Symbol,Any}()
    quiet = Peridynamics.get_submit_options(o)
    @test quiet == false
end
