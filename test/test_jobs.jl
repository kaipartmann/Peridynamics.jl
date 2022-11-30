using Peridynamics
using Test

pdsba1 = PDSingleBodyAnalysis(;
    name="pdsba1",
    pc=PointCloud(1, 1, 1, 0.5),
    mat=BondBasedMaterial(; horizon=1, rho=1, E=1, Gc=1),
    td=TimeDiscretization(10),
    es=ExportSettings(),
)
@test pdsba1.es.resultfile_prefix == "pdsba1"
@test pdsba1.es.logfile =="pdsba1.log"
@test pdsba1.es.exportfreq == 11
@test isempty(pdsba1.precracks)
@test isempty(pdsba1.bcs)
@test isempty(pdsba1.ics)

pdsba2 = PDSingleBodyAnalysis(;
    name="pdsba2",
    pc=PointCloud(1, 1, 1, 0.5),
    mat=BondBasedMaterial(; horizon=1, rho=1, E=1, Gc=1),
    td=TimeDiscretization(10),
    es=ExportSettings("test",2),
)
@test pdsba2.es.resultfile_prefix == joinpath("test", "pdsba2")
@test pdsba2.es.logfile == joinpath("test", "pdsba2.log")
@test pdsba2.es.exportfreq == 2
@test isempty(pdsba2.precracks)
@test isempty(pdsba2.bcs)
@test isempty(pdsba2.ics)
