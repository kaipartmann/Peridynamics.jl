@testitem "JobOptions" begin
    function tests_body_specific(body::Body, solver)
        default_fields_body = [field for field in Peridynamics.default_export_fields()]
        body_name = string(Peridynamics.get_name(body))
        vtk_root = joinpath("rootpath", "vtk")
        vtk_filebase = joinpath(vtk_root,
                                isempty(body_name) ? "timestep" : body_name * "_timestep")

        o = Dict{Symbol,Any}(:path => "rootpath", :freq => 10)
        jo = Peridynamics.get_job_options(body, solver, o)
        @test jo.export_allowed == true
        @test jo.root == "rootpath"
        @test jo.vtk == joinpath("rootpath", "vtk")
        @test jo.logfile == joinpath("rootpath", "logfile.log")
        @test jo.freq == 10
        @test jo.fields == default_fields_body
        @test jo.vtk_filebase == vtk_filebase

        o = Dict{Symbol,Any}(:path => "rootpath")
        jo = Peridynamics.get_job_options(body, solver, o)
        @test jo.export_allowed == true
        @test jo.root == "rootpath"
        @test jo.vtk == joinpath("rootpath", "vtk")
        @test jo.logfile == joinpath("rootpath", "logfile.log")
        @test jo.freq == 10
        @test jo.fields == default_fields_body
        @test jo.vtk_filebase == vtk_filebase

        o = Dict{Symbol,Any}(:freq => 10)
        msg = "if `freq` is specified, the keyword `path` is also needed!\n"
        @test_throws ArgumentError(msg) begin
            Peridynamics.get_job_options(body, solver, o)
        end

        o = Dict{Symbol,Any}()
        jo = Peridynamics.get_job_options(body, solver, o)
        @test jo.export_allowed == false
        @test jo.root == ""
        @test jo.vtk == ""
        @test jo.logfile == ""
        @test jo.freq == 0
        @test jo.fields == Vector{Symbol}()
        @test jo.vtk_filebase == ""

        o = Dict{Symbol,Any}(:path => "rootpath", :freq => -10)
        msg = "`freq` should be larger than zero!\n"
        @test_throws ArgumentError(msg) begin
            Peridynamics.get_job_options(body, solver, o)
        end

        o = Dict{Symbol,Any}(:path => "rootpath", :fields => (:displacement, :b_ext))
        jo = Peridynamics.get_job_options(body, solver, o)
        @test jo.export_allowed == true
        @test jo.root == "rootpath"
        @test jo.vtk == joinpath("rootpath", "vtk")
        @test jo.logfile == joinpath("rootpath", "logfile.log")
        @test jo.freq == 10
        @test jo.fields == [:displacement, :b_ext]
        @test jo.vtk_filebase == vtk_filebase

        o = Dict{Symbol,Any}(:path => "rootpath", :fields => :damage)
        jo = Peridynamics.get_job_options(body, solver, o)
        @test jo.export_allowed == true
        @test jo.root == "rootpath"
        @test jo.vtk == joinpath("rootpath", "vtk")
        @test jo.logfile == joinpath("rootpath", "logfile.log")
        @test jo.freq == 10
        @test jo.fields == [:damage]
        @test jo.vtk_filebase == vtk_filebase

        o = Dict{Symbol,Any}(:path => "rootpath", :fields => [:b_ext])
        jo = Peridynamics.get_job_options(body, solver, o)
        @test jo.export_allowed == true
        @test jo.root == "rootpath"
        @test jo.vtk == joinpath("rootpath", "vtk")
        @test jo.logfile == joinpath("rootpath", "logfile.log")
        @test jo.freq == 10
        @test jo.fields == [:b_ext]
        @test jo.vtk_filebase == vtk_filebase

        o = Dict{Symbol,Any}(:path => "rootpath", :fields => (:displacement, :b_ext, :abcd))
        @test_throws ArgumentError begin
            Peridynamics.get_job_options(body, solver, o)
        end

        return nothing
    end

    function tests_multibody_specific(ms::MultibodySetup, solver)
        default_fields_body = [field for field in Peridynamics.default_export_fields()]
        default_fields_multibody = Dict{Symbol,Vector{Symbol}}()
        vtk_filebase = Dict{Symbol,String}()
        vtk_root = joinpath("rootpath", "vtk")
        for body_name in Peridynamics.each_body_name(ms)
            default_fields_multibody[body_name] = default_fields_body
            vtk_filebase[body_name] = joinpath(vtk_root, string(body_name) * "_timestep")
        end

        o = Dict{Symbol,Any}(:path => "rootpath", :freq => 10)
        jo = Peridynamics.get_job_options(ms, solver, o)
        @test jo.export_allowed == true
        @test jo.root == "rootpath"
        @test jo.vtk == joinpath("rootpath", "vtk")
        @test jo.logfile == joinpath("rootpath", "logfile.log")
        @test jo.freq == 10
        @test jo.fields == default_fields_multibody
        @test jo.vtk_filebase == vtk_filebase

        o = Dict{Symbol,Any}(:path => "rootpath")
        jo = Peridynamics.get_job_options(ms, solver, o)
        @test jo.export_allowed == true
        @test jo.root == "rootpath"
        @test jo.vtk == joinpath("rootpath", "vtk")
        @test jo.logfile == joinpath("rootpath", "logfile.log")
        @test jo.freq == 10
        @test jo.fields == default_fields_multibody
        @test jo.vtk_filebase == vtk_filebase

        o = Dict{Symbol,Any}(:freq => 10)
        msg = "if `freq` is specified, the keyword `path` is also needed!\n"
        @test_throws ArgumentError(msg) begin
            Peridynamics.get_job_options(ms, solver, o)
        end

        o = Dict{Symbol,Any}()
        jo = Peridynamics.get_job_options(ms, solver, o)
        @test jo.export_allowed == false
        @test jo.root == ""
        @test jo.vtk == ""
        @test jo.logfile == ""
        @test jo.freq == 0
        @test jo.fields == Dict{Symbol,Vector{Symbol}}()
        @test jo.vtk_filebase == Dict{Symbol,String}()

        o = Dict{Symbol,Any}(:path => "rootpath", :freq => -10)
        msg = "`freq` should be larger than zero!\n"
        @test_throws ArgumentError(msg) begin
            Peridynamics.get_job_options(ms, solver, o)
        end

        o = Dict{Symbol,Any}(:path => "rootpath", :fields => (:displacement, :b_ext))
        jo = Peridynamics.get_job_options(ms, solver, o)
        @test jo.export_allowed == true
        @test jo.root == "rootpath"
        @test jo.vtk == joinpath("rootpath", "vtk")
        @test jo.logfile == joinpath("rootpath", "logfile.log")
        @test jo.freq == 10
        fields = [:displacement, :b_ext]
        fields_spec = Dict{Symbol,Vector{Symbol}}()
        for body_name in Peridynamics.each_body_name(ms)
            fields_spec[body_name] = fields
        end
        @test jo.fields == fields_spec
        @test jo.vtk_filebase == vtk_filebase

        o = Dict{Symbol,Any}(:path => "rootpath", :fields => :damage)
        jo = Peridynamics.get_job_options(ms, solver, o)
        @test jo.export_allowed == true
        @test jo.root == "rootpath"
        @test jo.vtk == joinpath("rootpath", "vtk")
        @test jo.logfile == joinpath("rootpath", "logfile.log")
        @test jo.freq == 10
        fields = [:damage]
        fields_spec = Dict{Symbol,Vector{Symbol}}()
        for body_name in Peridynamics.each_body_name(ms)
            fields_spec[body_name] = fields
        end
        @test jo.fields == fields_spec
        @test jo.vtk_filebase == vtk_filebase

        o = Dict{Symbol,Any}(:path => "rootpath", :fields => [:b_int])
        jo = Peridynamics.get_job_options(ms, solver, o)
        @test jo.export_allowed == true
        @test jo.root == "rootpath"
        @test jo.vtk == joinpath("rootpath", "vtk")
        @test jo.logfile == joinpath("rootpath", "logfile.log")
        @test jo.freq == 10
        fields = [:b_int]
        fields_spec = Dict{Symbol,Vector{Symbol}}()
        for body_name in Peridynamics.each_body_name(ms)
            fields_spec[body_name] = fields
        end
        @test jo.fields == fields_spec
        @test jo.vtk_filebase == vtk_filebase

        o = Dict{Symbol,Any}(:path => "rootpath", :fields => (:displacement, :b_ext, :abcd))
        @test_throws ArgumentError begin
            Peridynamics.get_job_options(ms, solver, o)
        end

        o = Dict{Symbol,Any}(:path => "rootpath", :fields => ":displacement")
        @test_throws ArgumentError begin
            Peridynamics.get_job_options(ms, solver, o)
        end

        o = Dict{Symbol,Any}(:path => "rootpath", :fields => (":displacement", ":damage"))
        @test_throws ArgumentError begin
            Peridynamics.get_job_options(ms, solver, o)
        end

        o = Dict{Symbol,Any}(:path => "rootpath", :fields => [":displacement", ":damage"])
        @test_throws ArgumentError begin
            Peridynamics.get_job_options(ms, solver, o)
        end

        return nothing
    end

    # setup
    bbb = Body(BBMaterial(), rand(3, 10), rand(10))
    bosb = Body(OSBMaterial(), rand(3, 10), rand(10))
    bnosb = Body(NOSBMaterial(), rand(3, 10), rand(10))
    ms = MultibodySetup(:a => bbb, :b => bosb, :c => bnosb)
    vv = VelocityVerlet(steps=1)
    dr = DynamicRelaxation(steps=1)

    tests_body_specific(bbb, vv)
    tests_body_specific(bosb, vv)
    tests_body_specific(bnosb, vv)
    tests_body_specific(bbb, dr)
    tests_body_specific(bosb, dr)
    tests_body_specific(bnosb, dr)

    tests_multibody_specific(ms, vv)

    default_fields_body = [field for field in Peridynamics.default_export_fields()]
    default_fields_multibody = Dict(:a => default_fields_body, :b => default_fields_body,
                                    :c => default_fields_body)

    o = Dict{Symbol,Any}(:path => "rootpath", :fields => Dict(:a => :damage))
    jo = Peridynamics.get_job_options(ms, vv, o)
    @test jo.export_allowed == true
    @test jo.root == "rootpath"
    @test jo.vtk == joinpath("rootpath", "vtk")
    @test jo.logfile == joinpath("rootpath", "logfile.log")
    @test jo.freq == 10
    @test jo.fields == Dict(:a => [:damage], :b => default_fields_body,
               :c => default_fields_body)

    o = Dict{Symbol,Any}(:path => "rootpath", :fields => Dict(:b => ("test", :damage)))
    @test_throws ArgumentError begin
        Peridynamics.get_job_options(ms, vv, o)
    end
end
