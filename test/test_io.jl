using Peridynamics
using Test

##------------------------------------------------------------------------------------------
# export_vtk

@testset "export_vtk" begin
    if Threads.nthreads() <= 2
        positions = [
            0.0 1.0
            0.0 0.0
            0.0 0.0
        ]
        point_spacing = 1.0
        δ = 1.5 * point_spacing
        E = 210e9
        n_points = 2
        volumes = fill(point_spacing^3, n_points)
        pc = PointCloud(positions, volumes)
        mat = BondBasedMaterial(horizon=δ, rho=7850.0, E=E, Gc=1.0)
        body = Peridynamics.create_simmodel(mat, pc)

        # export vtk
        rm.(filter(x->endswith(x,".vtu"), readdir(@__DIR__; join=true)), force=true)
        Peridynamics.export_vtk(body, "testfile", 0, 0.0)
        @test isfile("testfile_t0.vtu")
        rm.(filter(x->endswith(x,".vtu"), readdir(@__DIR__; join=true)), force=true)
    else
        @warn "Test omitted! Threads.nthreads() should be <= 2"
    end
end

##------------------------------------------------------------------------------------------
# ExportSettings

@testset "ExportSettings" begin
    es2 = ExportSettings("test/path", 10)
    @test es2.path == "test/path"
    @test es2.exportfreq == 10
    @test es2.resultfile_prefix == ""
    @test es2.logfile == ""
    @test es2.exportflag == true
end

@testset "ExportSettings empty" begin
    es1 = ExportSettings()
    @test es1.path == ""
    @test es1.exportfreq == 0
    @test es1.resultfile_prefix == ""
    @test es1.logfile == ""
    @test es1.exportflag == false
end
