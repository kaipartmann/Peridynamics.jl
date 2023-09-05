"""
All calculated or defined parameters that should not change ever inside the time loop!
--> independent of multithreading!
--> preprocessing of the problem
"""
struct SimParameters{M <: AbstractPDMaterial, T <: AbstractTimeDiscretization}
    name::String
    pc::PointCloud
    unique_bonds::Bool
    n_bonds::Int
    bond_data::Vector{Tuple{Int, Int, Float64, Bool}}
    n_family_members::Vector{Int}
    owned_points::Vector{UnitRange{Int}}
    owned_bonds::Vector{UnitRange{Int}}
    single_tids::Vector{Tuple{Int, Int}}
    multi_tids::Vector{Tuple{Int, Vector{Int}}}
    mat::M
    cells::Vector{MeshCell{VTKCellType, Tuple{Int64}}}
    precracks::Vector{PreCrack}
    bcs::Vector{<:AbstractBC}
    ics::Vector{<:AbstractIC}
    td::T
    es::ExportSettings
end

"""
All things that are calculated and change globally ???
"""
struct GlobalStorage
    owned_points::Vector{UnitRange{Int}}
    owned_bonds::Vector{UnitRange{Int}}
    single_tids::Vector{Tuple{Int, Int}}
    multi_tids::Vector{Tuple{Int, Vector{Int}}}
    position::Matrix{Float64}
    displacement::Matrix{Float64}
    velocity::Matrix{Float64}
    velocity_half::Matrix{Float64}
    acceleration::Matrix{Float64}
    damage::Matrix{Float64}
    bond_failure::Vector{Int}
end

"""
All things that are thread locally
"""
struct ThreadLocalStorage
    pointidmap::Dict{Int, Int}
    bondidmap::Dict{Int, Int}
    b_int::Matrix{Float64}
    n_active_family_members::Vector{Int}
end

"""

"""
struct PDProblem
    sp::SimParameters
    gs::GlobalStorage
    tls::Vector{ThreadLocalStorage}
end

function init_pdproblem(sim::AbstractPDAnalysis)
    n_threads = nthreads()
    name = sim.name
    pc = sim.pc
    mat = sim.mat
    td = sim.td
    bcs = sim.bcs
    ics = sim.ics
    es = sim.es
    precracks = sim.precracks
    unique_bonds = true
    owned_points = defaultdist(pc.n_points, n_threads)
    point_mappings = find_mappings(owned_points)

    # find the bonds
    bond_data, n_family_members = find_unique_bonds(pc, mat, owned_points)
    n_bonds = length(bond_data)
    owned_bonds = defaultdist(n_bonds, n_threads)
    bond_mappings = find_mappings(owned_bonds)

    n_active_family_members = zeros(Int, (pc.n_points, n_threads))
    _sum_tids = zeros(Bool, (n_points, n_threads))
    _sum_tids .= false
    @threads for tid in 1:n_threads
        for current_bond in owned_bonds[tid]
            (i, j, _, _) = bond_data[current_bond]
            n_active_family_members[i, tid] += 1
            n_active_family_members[j, tid] += 1
            _sum_tids[i, tid] = true
            _sum_tids[j, tid] = true
        end
    end
    sum_tids = [findall(row) for row in eachrow(_sum_tids)]
    single_tids, multi_tids = find_tids(sum_tids)

    cells = get_cells(n_points)
    position = pc.position
    displacement = zeros(Float64, (3, n_points))
    velocity = zeros(Float64, (3, n_points))
    velocity_half = zeros(Float64, (3, n_points))
    acceleration = zeros(Float64, (3, n_points))
    b_ext = zeros(Float64, (3, n_points))
    b_int = zeros(Float64, (3, n_points))
    damage = zeros(Int, n_points)

    pdparams = SimParameters(name, pc, unique_bonds, n_bonds, bond_data, n_family_members,
                             owned_points, owned_bonds, single_tids, multi_tids, mat, cells,
                             precracks, bcs, ics, td, es)

    return pdparams
end
