"""
    init_body(mat::PDMaterial, pc::PointCloud)

Material interface function for initialization of the `body <: AbstractPDBody` used for
the simulation.
"""
function init_body end

function compute_forcedensity! end
