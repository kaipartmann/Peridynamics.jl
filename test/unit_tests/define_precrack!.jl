using Peridynamics, Test

precrack = PreCrack([1], [2])

# BondBasedBody
n_points = 2
n_bonds = 1
n_threads = 1
unique_bonds = true
owned_points = UnitRange{Int64}[1:2]
owned_bonds = UnitRange{Int64}[1:1]
single_tids = Tuple{Int64, Int64}[]
multi_tids = Tuple{Int64, Vector{Int64}}[]
bond_data = Tuple{Int64, Int64, Float64, Bool}[(1, 2, 1.0, 1)]
volume = [1.0, 1.0]
cells = Peridynamics.get_cells(n_points)
n_family_members = [1, 1]
n_active_family_members = [1; 1;;]
position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
displacement = [0.0 0.0; 0.0 0.0; 0.0 0.0]
velocity = [0.0 0.0; 0.0 0.0; 0.0 0.0]
velocity_half = [0.0 0.0; 0.0 0.0; 0.0 0.0]
acceleration = [0.0 0.0; 0.0 0.0; 0.0 0.0]
b_int = [0.0 0.0; 0.0 0.0; 0.0 0.0;;;]
b_ext = [0.0 0.0; 0.0 0.0; 0.0 0.0]
damage = [0.0, 0.0]
bond_failure = [1]
body = Peridynamics.BondBasedBody(n_points, n_bonds, n_threads, unique_bonds, owned_points,
                           owned_bonds, single_tids, multi_tids, bond_data, volume,
                           cells, n_family_members, n_active_family_members, position,
                           displacement, velocity, velocity_half, acceleration, b_int,
                           b_ext, damage, bond_failure)

@test body.n_active_family_members == [1; 1;;]
@test body.bond_failure == [1]

Peridynamics.define_precrack!(body, precrack)

@test body.n_active_family_members == [0; 0;;]
@test body.bond_failure == [0]
