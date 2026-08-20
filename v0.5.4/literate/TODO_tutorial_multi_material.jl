# TODO

using Peridynamics
using GLMakie

lxy = 1
lz = 0.1lxy
nz = 5
dx = lz/nz
pc = PointCloud(lxy, lxy, lz, dx)

E_variation = 60e3:250e3
rho_variation = 2e-6:8e-6
epsilon_c_variation = 0.002:0.00001:0.01

pore_radius_variation = 0.05lxy:0.0001lxy:0.1lxy

n_pores = 200

mats = Tuple(BBMaterial(horizon=3dx, rho=rand(rho_variation), E=rand(E_variation), epsilon_c=rand(epsilon_c_variation)) for _ in 1:n_pores)

matofpoint = ones(Int, pc.n_points)
radius_of_pores = zeros(n_pores)

for i in 1:n_pores
    # get random pore position inside model
	pos = [rand(-lxy/2:0.1dx:lxy/2), rand(-lxy/2:0.1dx:lxy/2), rand(-lz/2:0.1dx:lz/2)]
    # get random pore radius
	rad = rand(pore_radius_variation)
	# add pore radius to radius vector
    radius_of_pores[i] = rad
	# get points inside of pore
    for (point_id, point) in enumerate(eachcol(pc.position))
        dist = sqrt((point[1]-pos[1])^2 + (point[2]-pos[2])^2 + (point[3]-pos[3])^2)
        if dist ≤ rad
            # assign material of pore to these points
            matofpoint[point_id] = i
        end
    end
end

# mat = MultiMaterial(mats, matofpoint)

pore_size = radius_of_pores[matofpoint]

#-
fig = Figure(resolution = (1000,1000), figure_padding=0) #src
ax = Axis3(fig[1,1]; aspect = :data, elevation=0.37π, azimuth=0.45π) #src
hidespines!(ax) #src
hidedecorations!(ax) #src
meshscatter!(ax, pc.position; markersize=0.5dx, color=pore_size, colormap=:lightrainbow) #src
save(joinpath(@__DIR__, "..", "assets", "multi_material_test.png"), fig; px_per_unit=3)
fig
