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

@testitem "export_results Body" begin
    temp_root = joinpath(@__DIR__, "temp_root_export_results_test")
    rm(temp_root; recursive=true, force=true)

    pos, vol = uniform_box(1, 1, 1, 0.5)
    body = Body(BBMaterial(), pos, vol)
    material!(body, horizon=1.5, E=1, rho=1, Gc=1)
    vv = VelocityVerlet(steps=1)
    dh = Peridynamics.threads_data_handler(body, vv, 1)
    chunk = dh.chunks[1]

    o = Dict{Symbol,Any}(:path => temp_root, :freq => 1)
    options = Peridynamics.get_job_options(body, vv, o)

    u_rand = rand(3, 8)
    chunk.storage.displacement .= u_rand

    n = 1
    t = 0.0001
    chunk_id = 1
    n_chunks = 1

    Peridynamics.export_results(dh, options, chunk_id, n, t)

    pvtu_file = joinpath(temp_root, "vtk", "timestep_000001.pvtu")
    @test isfile(pvtu_file)

    r = read_vtk(pvtu_file)
    @test r[:position] ≈ pos
    @test r[:displacement] ≈ u_rand
    @test r[:damage] == zeros(8)
    @test r[:time] == [t]

    rm(temp_root; recursive=true, force=true)
end

@testitem "export_results MultibodySetup" begin
    temp_root = joinpath(@__DIR__, "temp_root_export_results_test")
    rm(temp_root; recursive=true, force=true)

    pos_a, vol_a = uniform_box(1, 1, 1, 0.5)
    body_a = Body(BBMaterial(), pos_a, vol_a)
    material!(body_a, horizon=1.5, E=1, rho=1, Gc=1)

    pos_b, vol_b = uniform_box(1, 1, 1, 0.5, center_z=1.5)
    body_b = Body(BBMaterial(), pos_b, vol_b)
    material!(body_b, horizon=1.5, E=1, rho=1, Gc=1)

    ms = MultibodySetup(:a => body_a, :b => body_b)

    vv = VelocityVerlet(steps=1)

    dh = Peridynamics.threads_data_handler(ms, vv, 1)
    dh_a = dh.body_dhs[1]
    chunk_a = dh_a.chunks[1]
    @test chunk_a.body_name === :a
    dh_b = dh.body_dhs[2]
    chunk_b = dh_b.chunks[1]
    @test chunk_b.body_name === :b

    o = Dict{Symbol,Any}(:path => temp_root, :freq => 1)
    options = Peridynamics.get_job_options(ms, vv, o)

    u_rand_a = rand(3, 8)
    u_rand_b = rand(3, 8)
    chunk_a.storage.displacement .= u_rand_a
    chunk_b.storage.displacement .= u_rand_b

    n = 1
    t = 0.0001
    chunk_id = 1
    n_chunks = 1

    Peridynamics.export_results(dh_a, options, chunk_id, n, t)
    Peridynamics.export_results(dh_b, options, chunk_id, n, t)

    pvtu_file_a = joinpath(temp_root, "vtk", "a_timestep_000001.pvtu")
    @test isfile(pvtu_file_a)
    pvtu_file_b = joinpath(temp_root, "vtk", "b_timestep_000001.pvtu")
    @test isfile(pvtu_file_b)

    r_a = read_vtk(pvtu_file_a)
    @test r_a[:position] ≈ pos_a
    @test r_a[:displacement] ≈ u_rand_a
    @test r_a[:damage] == zeros(8)
    @test r_a[:time] == [t]

    r_b = read_vtk(pvtu_file_b)
    @test r_b[:position] ≈ pos_b
    @test r_b[:displacement] ≈ u_rand_b
    @test r_b[:damage] == zeros(8)
    @test r_b[:time] == [t]

    rm(temp_root; recursive=true, force=true)
end
