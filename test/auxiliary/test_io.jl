@testitem "ExportOptions BBVerletStorage" begin
    o = Dict{Symbol,Any}(:path => "rootpath", :freq => 10)
    eo = Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    @test eo.exportflag == true
    @test eo.root == "rootpath"
    @test eo.vtk == joinpath("rootpath", "vtk")
    @test eo.logfile == joinpath("rootpath", "logfile.log")
    @test eo.freq == 10
    @test eo.fields == ((:displacement, Matrix{Float64}), (:damage, Vector{Float64}))

    o = Dict{Symbol,Any}(:path => "rootpath")
    eo = Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    @test eo.exportflag == true
    @test eo.root == "rootpath"
    @test eo.vtk == joinpath("rootpath", "vtk")
    @test eo.logfile == joinpath("rootpath", "logfile.log")
    @test eo.freq == 10
    @test eo.fields == ((:displacement, Matrix{Float64}), (:damage, Vector{Float64}))

    o = Dict{Symbol,Any}(:freq => 10)
    msg = "if `freq` is spedified, the keyword `path` is also needed!\n"
    @test_throws ArgumentError(msg) begin
        Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    end

    o = Dict{Symbol,Any}()
    eo = Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    @test eo.exportflag == false
    @test eo.root == ""
    @test eo.vtk == ""
    @test eo.logfile == ""
    @test eo.freq == 0
    @test eo.fields == NTuple{0,Tuple{Symbol,DataType}}()

    o = Dict{Symbol,Any}(:path => "rootpath", :freq => -10)
    msg = "`freq` should be larger than zero!\n"
    @test_throws ArgumentError(msg) begin
        Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    end

    o = Dict{Symbol,Any}(:path => "rootpath", :write => (:displacement, :damage, :b_int))
    eo = Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    @test eo.exportflag == true
    @test eo.root == "rootpath"
    @test eo.vtk == joinpath("rootpath", "vtk")
    @test eo.logfile == joinpath("rootpath", "logfile.log")
    @test eo.freq == 10
    @test eo.fields == ((:displacement, Matrix{Float64}),
                        (:damage, Vector{Float64}),
                        (:b_int, Matrix{Float64}))
end

@testitem "export_results" begin

end
