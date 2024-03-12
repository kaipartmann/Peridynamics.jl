@testitem "ExportOptions BBVerletStorage" begin
    o = Dict{Symbol,Any}(:path => "rootpath", :freq => 10)
    eo = Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    @test eo.exportflag == true
    @test eo.root == abspath("rootpath")
    @test eo.vtk == abspath(joinpath("rootpath", "vtk"))
    @test eo.logfile == abspath(joinpath("rootpath", "logfile.log"))
    @test eo.freq == 10
    @test eo.eff == Peridynamics.export_disp_and_dmg

    o = Dict{Symbol,Any}(:path => "rootpath")
    eo = Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    @test eo.exportflag == true
    @test eo.root == abspath("rootpath")
    @test eo.vtk == abspath(joinpath("rootpath", "vtk"))
    @test eo.logfile == abspath(joinpath("rootpath", "logfile.log"))
    @test eo.freq == 10
    @test eo.eff == Peridynamics.export_disp_and_dmg

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
    @test eo.eff isa Function
    @test eo.eff() == ()

    o = Dict{Symbol,Any}(:path => "rootpath", :freq => -10)
    msg = "`freq` should be larger than zero!\n"
    @test_throws ArgumentError(msg) begin
        Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    end

    eff = x -> (("disp", x.displacement), ("dmg", x.damage), ("b_int", x.b_int))
    o = Dict{Symbol,Any}(:path => "rootpath", :eff => eff)
    eo = Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    @test eo.exportflag == true
    @test eo.root == abspath("rootpath")
    @test eo.vtk == abspath(joinpath("rootpath", "vtk"))
    @test eo.logfile == abspath(joinpath("rootpath", "logfile.log"))
    @test eo.freq == 10
    @test eo.eff == eff
end

@testitem "exported fields function" begin
    position = rand(3,10)
    displacement = rand(3,10)
    velocity = rand(3,10)
    velocity_half = rand(3,10)
    acceleration = rand(3,10)
    b_int = rand(3,10)
    b_ext = rand(3,10)
    damage = rand(10)
    bond_active = rand(Bool,20)
    n_active_bonds = rand(Int,20)
    s = Peridynamics.BBVerletStorage(position, displacement, velocity, velocity_half,
                                     acceleration, b_int, b_ext, damage, bond_active,
                                     n_active_bonds)

    o = Dict{Symbol,Any}(:path => "rootpath")
    eo = Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    tpl = eo(s)
    @test length(tpl) == 2
    @test tpl[1] == ("displacement", displacement)
    @test tpl[2] == ("damage", damage)

    _eff = x -> (("disp", x.displacement), ("dmg", x.damage), ("b_int", x.b_int))
    o = Dict{Symbol,Any}(:path => "rootpath", :eff => _eff)
    eo = Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    tpl = eo(s)
    @test length(tpl) == 3
    @test tpl[1] == ("disp", displacement)
    @test tpl[2] == ("dmg", damage)
    @test tpl[3] == ("b_int", b_int)

    _eff = @eff :displacement :acceleration :damage :b_int
    o = Dict{Symbol,Any}(:path => "rootpath", :eff => _eff)
    eo = Peridynamics.get_export_options(Peridynamics.BBVerletStorage, o)
    tpl = eo(s)
    @test length(tpl) == 4
    @test tpl[1] == ("displacement", displacement)
    @test tpl[2] == ("acceleration", acceleration)
    @test tpl[3] == ("damage", damage)
    @test tpl[4] == ("b_int", b_int)
end
