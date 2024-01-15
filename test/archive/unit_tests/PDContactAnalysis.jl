using Peridynamics, Test

bs1 = BodySetup(PointCloud(1, 1, 1, 0.5), BBMaterial(; horizon=1, rho=1, E=1, Gc=1))
c1 = Contact((1,2), 0.1)

pdca1 = PDContactAnalysis(;
    name="pdca1",
    body_setup=[bs1, bs1],
    contact=[c1],
    td=VelocityVerlet(10),
    es=ExportSettings(),
)
@test pdca1.es.resultfile_prefix == "pdca1"
@test pdca1.es.logfile == "pdca1.log"
@test pdca1.es.exportfreq == typemax(Int)

pdca2 = PDContactAnalysis(;
    name="pdca2",
    body_setup=[bs1, bs1],
    contact=[c1],
    td=VelocityVerlet(10),
    es=ExportSettings("test",2),
)
@test pdca2.es.resultfile_prefix == joinpath("test", "pdca2")
@test pdca2.es.logfile == joinpath("test", "pdca2.log")
@test pdca2.es.exportfreq == 2

@test_throws TypeError PDContactAnalysis(;
    name="pdca2",
    body_setup=[bs1, bs1],
    contact=[c1],
    td=DynamicRelaxation(10),
    es=ExportSettings("test",2),
)
