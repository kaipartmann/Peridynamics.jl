using Peridynamics
pointcloud = read_inp("examples/models/TensileTestMesh.inp")
# pointcloud.failure_flag[pointcloud.point_sets["bottom"]] .= false
# pointcloud.failure_flag[pointcloud.point_sets["top"]] .= false
material = BondBasedMaterial(; horizon=0.008, rho=2700, E=70e9, Gc=100)
boundary_conditions = [
    VelocityBC(t -> -0.8, pointcloud.point_sets["bottom"], 1),
    VelocityBC(t -> 0.8, pointcloud.point_sets["top"], 1),
]
job = PDSingleBodyAnalysis(;
    name="TensileTest",
    pc=pointcloud,
    mat=material,
    bcs=boundary_conditions,
    td=TimeDiscretization(500),
    es=ExportSettings("examples/results/TensileTest", 10),
)
submit(job)
