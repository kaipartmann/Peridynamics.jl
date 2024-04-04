@testitem "ExportOptions BBVerletStorage" begin
    o = Dict{Symbol,Any}(:path => "rootpath", :freq => 10)
    eo = Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    @test eo.exportflag == true
    @test eo.root == "rootpath"
    @test eo.vtk == joinpath("rootpath", "vtk")
    @test eo.logfile == joinpath("rootpath", "logfile.log")
    @test eo.freq == 10
    @test eo.fields == [field for field in Peridynamics.DEFAULT_EXPORT_FIELDS]

    o = Dict{Symbol,Any}(:path => "rootpath")
    eo = Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    @test eo.exportflag == true
    @test eo.root == "rootpath"
    @test eo.vtk == joinpath("rootpath", "vtk")
    @test eo.logfile == joinpath("rootpath", "logfile.log")
    @test eo.freq == 10
    @test eo.fields == [field for field in Peridynamics.DEFAULT_EXPORT_FIELDS]

    o = Dict{Symbol,Any}(:freq => 10)
    msg = "if `freq` is specified, the keyword `path` is also needed!\n"
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
    @test eo.fields == Symbol[]

    o = Dict{Symbol,Any}(:path => "rootpath", :freq => -10)
    msg = "`freq` should be larger than zero!\n"
    @test_throws ArgumentError(msg) begin
        Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    end

    o = Dict{Symbol,Any}(:path => "rootpath", :fields => (:displacement, :b_ext))
    eo = Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    @test eo.exportflag == true
    @test eo.root == "rootpath"
    @test eo.vtk == joinpath("rootpath", "vtk")
    @test eo.logfile == joinpath("rootpath", "logfile.log")
    @test eo.freq == 10
    @test eo.fields == [:displacement, :b_ext]
end
