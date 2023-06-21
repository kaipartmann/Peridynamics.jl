struct ContinuumBasedMaterial <: AbstractPDMaterial
    δ::Float64
    rho::Float64
    E::Float64
    nu::Float64
    G::Float64
    K::Float64
    C1::Float64
    C2::Float64
    C3::Float64
    Gc::Float64
    εc::Float64
end

function ContinuumBasedMaterial(;
    horizon::Real,
    rho::Real,
    E::Real,
    nu::Real,
    Gc::Real=-1,
    epsilon_c::Real=-1,
    C1::Real=0,
    C2::Real=0,
    C3::Real=0,
)
    K = E / (3 * (1 - 2 * nu))
    G = E / (2 * (1 + nu))
    λ = E * nu /((1 + nu) * (1 - 2nu))
    μ = G
    if (Gc !== -1) && (epsilon_c == -1)
        epsilon_c = sqrt(5.0 * Gc / (9.0 * K * horizon))
    elseif (Gc == -1) && (epsilon_c !== -1)
        Gc = 9.0/5.0 * K * horizon * epsilon_c^2
    elseif (Gc !== -1) && (epsilon_c !== -1)
        msg = "Duplicate definition: define either Gc or epsilon_c, not both!"
        throw(ArgumentError(msg))
    elseif (Gc == -1) && (epsilon_c == -1)
        throw(ArgumentError("Either Gc or epsilon_c have to be defined!"))
    end
    if C1 == 0 && C2 == 0 && C3 == 0
        C1 = 30/π * μ/horizon^4
        C2 = 0.0
        C3 = 32/π^4 * (λ - μ)/horizon^12
    end
    return ContinuumBasedMaterial(
        horizon,
        rho,
        E,
        nu,
        G,
        K,
        C1,
        C2,
        C3,
        Gc,
        epsilon_c,
    )
end

function Base.show(io::IO, ::MIME"text/plain", mat::ContinuumBasedMaterial)
    println(io, typeof(mat), ":")
    for field in fieldnames(typeof(mat))
        @printf(io, "  %s: %g\n", string(field), getfield(mat, field))
    end
    return nothing
end

struct ContinuumBasedSimBody <: AbstractPDBody
    n_points::Int
    n_bonds::Int
    n_two_ni::Int
    n_three_ni::Int
    n_threads::Int
    unique_bonds::Bool
    has_only_one_ni::Bool
    has_one_and_two_ni::Bool
    has_one_and_three_ni::Bool
    has_one_two_three_ni::Bool
    owned_points::Vector{UnitRange{Int}}
    owned_bonds::Vector{UnitRange{Int}}
    owned_two_ni::Vector{UnitRange{Int}}
    owned_three_ni::Vector{UnitRange{Int}}
    single_tids::Vector{Tuple{Int,Int}}
    multi_tids::Vector{Tuple{Int,Vector{Int}}}
    bond_data::Vector{Tuple{Int,Int,Float64,Bool}}
    two_ni_data::Vector{Tuple{Int,Int,Float64}}
    three_ni_data::Vector{Tuple{Int,Int,Int,Float64}}
    volume::Vector{Float64}
    cells::Vector{MeshCell{VTKCellType,Tuple{Int64}}}
    n_family_members::Vector{Int}
    n_active_family_members::Matrix{Int}
    position::Matrix{Float64}
    displacement::Matrix{Float64}
    velocity::Matrix{Float64}
    velocity_half::Matrix{Float64}
    acceleration::Matrix{Float64}
    b_int::Array{Float64,3}
    b_ext::Matrix{Float64}
    damage::Vector{Float64}
    bond_failure::Vector{Int}
    volume_hood::Vector{Float64}
    volume_one_ni::Vector{Float64}
    volume_two_ni::Vector{Float64}
    volume_three_ni::Vector{Float64}
    energy::Vector{Float64}
end

function ContinuumBasedSimBody(
    mat::PDMaterial{ContinuumBasedMaterial},
    pc::PointCloud,
)
    n_threads = nthreads()
    n_points = pc.n_points
    @assert n_points >= n_threads "n_points < n_threads"
    owned_points = defaultdist(n_points, n_threads)
    volume = pc.volume
    cells = get_cells(n_points)
    position = pc.position
    displacement = zeros(Float64, (3, n_points))
    velocity = zeros(Float64, (3, n_points))
    velocity_half = zeros(Float64, (3, n_points))
    acceleration = zeros(Float64, (3, n_points))
    b_ext = zeros(Float64, (3, n_points))
    b_int = zeros(Float64, (3, n_points, n_threads))
    damage = zeros(Int, n_points)
    energy = Vector{Float64}(undef, 0) # TODO
    unique_bonds = false
    bond_data, n_family_members = find_bonds(pc, mat, owned_points)
    n_bonds = length(bond_data)
    owned_bonds = defaultdist(n_bonds, n_threads)
    bond_failure = ones(Int, n_bonds)
    n_active_family_members = zeros(Int, (n_points, n_threads))

    hood_range = find_hood_range(n_family_members, n_points)
    volume_hood = calc_volume_of_neighborhood(
        pc,
        bond_data,
        n_family_members,
        hood_range,
        mat,
    )
    volume_one_ni = [volume_hood[i]/n for (i, n) in enumerate(n_family_members)]

    # two ni
    has_two_ni = false
    if mat_has_twoni(mat)
        two_ni_data, n_two_ni, n_two_ni_per_point = find_two_ni(
            pc,
            mat,
            bond_data,
            hood_range,
            owned_points,
        )
        if !isempty(two_ni_data)
            has_two_ni = true
            volume_two_ni = [volume_hood[i]/n for (i, n) in enumerate(n_two_ni_per_point)]
        end
    end
    if !has_two_ni
        two_ni_data = Vector{Tuple{Int,Int,Float64}}(undef, 0)
        n_two_ni = 0
        n_two_ni_per_point = Vector{Int}(undef, 0)
        volume_two_ni = Vector{Float64}(undef, 0)
    end
    owned_two_ni = defaultdist(n_two_ni, n_threads)

    # three ni
    has_three_ni = false
    if mat_has_threeni(mat)
        three_ni_data, n_three_ni, n_three_ni_per_point = find_three_ni(
            pc,
            mat,
            bond_data,
            hood_range,
            owned_points,
        )
        if !isempty(three_ni_data)
            has_three_ni = true
            volume_three_ni =
                [volume_hood[i]/n for (i, n) in enumerate(n_three_ni_per_point)]
        end
    end
    if !has_three_ni
        three_ni_data = Vector{Tuple{Int,Int,Int,Float64}}(undef, 0)
        n_three_ni = 0
        n_three_ni_per_point = Vector{Int}(undef, 0)
        volume_three_ni = Vector{Float64}(undef, 0)
    end
    owned_three_ni = defaultdist(n_three_ni, n_threads)

    _sum_tids = zeros(Bool, (n_points, n_threads))
    _sum_tids .= false
    @threads for tid in 1:n_threads
        for current_bond in owned_bonds[tid]
            (i, _, _, _) = bond_data[current_bond]
            n_active_family_members[i, tid] += 1
            _sum_tids[i, tid] = true
        end
        for two_ni in owned_two_ni[tid]
            j_int, _, _ = two_ni_data[two_ni]
            i, _, _ = bond_data[j_int]
            _sum_tids[i, tid] = true
        end
        for three_ni in owned_three_ni[tid]
            j_int, _, _, _ = three_ni_data[three_ni]
            i, _, _ = bond_data[j_int]
            _sum_tids[i, tid] = true
        end
    end
    sum_tids = [findall(row) for row in eachrow(_sum_tids)]
    single_tids, multi_tids = find_tids(sum_tids)

    has_only_one_ni = false
    has_one_and_two_ni = false
    has_one_and_three_ni = false
    has_one_two_three_ni = false
    if has_two_ni && has_three_ni
        has_one_two_three_ni = true
    elseif has_two_ni && !has_three_ni
        has_one_and_two_ni = true
    elseif has_three_ni && !has_two_ni
        has_one_and_three_ni = true
    elseif !has_two_ni && !has_three_ni
        has_only_one_ni = true
    end

    return ContinuumBasedSimBody(
        n_points,
        n_bonds,
        n_two_ni,
        n_three_ni,
        n_threads,
        unique_bonds,
        has_only_one_ni,
        has_one_and_two_ni,
        has_one_and_three_ni,
        has_one_two_three_ni,
        owned_points,
        owned_bonds,
        owned_two_ni,
        owned_three_ni,
        single_tids,
        multi_tids,
        bond_data,
        two_ni_data,
        three_ni_data,
        volume,
        cells,
        n_family_members,
        n_active_family_members,
        position,
        displacement,
        velocity,
        velocity_half,
        acceleration,
        b_int,
        b_ext,
        damage,
        bond_failure,
        volume_hood,
        volume_one_ni,
        volume_two_ni,
        volume_three_ni,
        energy,
    )
end

function Base.show(io::IO, ::MIME"text/plain", body::ContinuumBasedSimBody)
    println(io, body.n_points, "-points ", typeof(body), ":")
    println(io, "  ", body.n_bonds, " bonds")
    println(io, "  ", body.n_two_ni, " two-neighbor interactions")
    println(io, "  ", body.n_three_ni, " three-neighbor interactions")
    return nothing
end

function mat_has_twoni(mat::MultiMaterial{ContinuumBasedMaterial})
    return sum([isapprox(m.C2, 0; atol=20eps()) ? 0 : 1 for m in mat[:]]) > 0 ? true : false
end

function mat_has_twoni(mat::ContinuumBasedMaterial)
    return isapprox(mat.C2, 0; atol=20eps()) ? false : true
end

function mat_has_threeni(mat::MultiMaterial{ContinuumBasedMaterial})
    return sum([isapprox(m.C3, 0; atol=20eps()) ? 0 : 1 for m in mat[:]]) > 0 ? true : false
end

function mat_has_threeni(mat::ContinuumBasedMaterial)
    return isapprox(mat.C3, 0; atol=20eps()) ? false : true
end

function find_hood_range(n_family_members::Vector{Int}, n_points::Int)
    hood_range = fill(0:0, n_points)
    current_one_ni = 1
    for i in 1:n_points
        hood_range[i] = current_one_ni:(current_one_ni + n_family_members[i] - 1)
        current_one_ni += n_family_members[i]
    end
    return hood_range
end

function calc_volume_of_neighborhood(
    pc::PointCloud,
    bond_data::Vector{Tuple{Int,Int,Float64,Bool}},
    n_family_members::Vector{Int},
    hood_range::Vector{UnitRange{Int}},
    mat::PDMaterial{ContinuumBasedMaterial}
)
    Δx_mean = calc_mean_point_spacing(pc, bond_data, n_family_members, hood_range)
    n_members_full_neighborhood = zeros(Int, pc.n_points)
    δ = [mat[i].δ for i in 1:pc.n_points]
    for i in 1:pc.n_points
        max_val = ceil(δ[i] / Δx_mean[i]) * Δx_mean[i]
        grid = range(start=-max_val, stop=max_val, step=Δx_mean[i])
        _position = hcat(([x; y; z] for x in grid for y in grid for z in grid)...)
        idx = findall(
            sqrt.(_position[1, :].^2 + _position[2, :].^2 + _position[3, :].^2) .<= δ[i]
        )
        n_members_full_neighborhood[i] = length(idx)
    end
    β = [n_family_members[i] / n_members_full_neighborhood[i] for i in 1:pc.n_points]
    volume_neighborhood = @. 4 / 3 * π * δ^3 * β
    return volume_neighborhood
end

function calc_mean_point_spacing(
    pc::PointCloud,
    bond_data::Vector{Tuple{Int,Int,Float64,Bool}},
    n_family_members::Vector{Int},
    hood_range::Vector{UnitRange{Int}},
)
    Δx_mean = zeros(pc.n_points)
    point_distance = [Xi for (_, _, Xi, _) in bond_data]
    Xi_min = [minimum(point_distance[hood_range[i]]) for i in 1:pc.n_points]
    for (i, j, _, _) in bond_data
        Δx_mean[i] += Xi_min[j]
    end
    for i in eachindex(Δx_mean)
        Δx_mean[i] /= n_family_members[i]
    end
    return Δx_mean
end

function find_two_ni(
    pc::PointCloud,
    mat::PDMaterial{ContinuumBasedMaterial},
    one_ni_data::Vector{Tuple{Int,Int,Float64,Bool}},
    hood_range::Vector{UnitRange{Int}},
    owned_points::Vector{UnitRange{Int}},
)
    n_threads = nthreads()
    _two_ni_data = fill([(0, 0, 0.0)], n_threads)
    number_of_twoni = zeros(Int, pc.n_points)
    p = Progress(pc.n_points;
        dt=1,
        desc="Two-NI search...    ",
        barlen=30,
        color=:normal,
        enabled=!is_logging(stderr),
    )
    owned_points_with_twoni = points_with_twoni(mat, owned_points)
    @threads for tid in 1:n_threads
        local_two_ni_data = Vector{Tuple{Int,Int,Float64}}(undef, 0)
        sizehint!(local_two_ni_data, pc.n_points * 1000)
        for i in owned_points_with_twoni[tid]
            num = 0
            δ = mat[i].δ
            ijseen = Set{Tuple{Int,Int}}()
            for j_int in hood_range[i], k_int in hood_range[i]
                _, j, _, _ = one_ni_data[j_int]
                _, k, _, _ = one_ni_data[k_int]
                if k !== j && !in((j, k), ijseen)
                    Ξijx = pc.position[1, j] - pc.position[1, i]
                    Ξijy = pc.position[2, j] - pc.position[2, i]
                    Ξijz = pc.position[3, j] - pc.position[3, i]
                    Ξikx = pc.position[1, k] - pc.position[1, i]
                    Ξiky = pc.position[2, k] - pc.position[2, i]
                    Ξikz = pc.position[3, k] - pc.position[3, i]
                    Ξjkx = pc.position[1, k] - pc.position[1, j]
                    Ξjky = pc.position[2, k] - pc.position[2, j]
                    Ξjkz = pc.position[3, k] - pc.position[3, j]
                    Ξjk = sqrt(Ξjkx * Ξjkx + Ξjky * Ξjky + Ξjkz * Ξjkz)
                    surface = surf_two_neigh(Ξijx, Ξijy, Ξijz, Ξikx, Ξiky, Ξikz)
                    if surface > eps() && Ξjk <= δ
                        num += 1
                        push!(local_two_ni_data, (j_int, k_int, surface))
                        push!(ijseen, (k, j))
                    end
                end
            end
            number_of_twoni[i] = num
            next!(p)
        end
        _two_ni_data[tid] = local_two_ni_data
    end
    finish!(p)
    two_ni_data = reduce(append!, _two_ni_data)
    n_two_ni = length(two_ni_data)

    return two_ni_data, n_two_ni, number_of_twoni
end

function points_with_twoni(::ContinuumBasedMaterial, owned_points::Vector{UnitRange{Int}})
    return owned_points
end

function get_points_with_interaction(
    mat::MultiMaterial{ContinuumBasedMaterial},
    owned_points::Vector{UnitRange{Int}},
    interaction::Symbol,
)
    nt = nthreads()
    _point_ids = Vector{Vector{Int}}(undef, nt)
    n_points_with_twoni = 0
    @threads for tid in 1:nt
        local_point_ids = Vector{Int}()
        for id in owned_points[tid]
            if isapprox(getfield(mat[id], interaction), 0)
                n_points_with_twoni += 1
                push!(local_point_ids, id)
            end
        end
        _point_ids[tid] = local_point_ids
    end
    point_ids = reduce(append!, _point_ids)
    owned_point_ids = defaultdist(n_points_with_twoni, nt)

    return [point_ids[opid] for opid in owned_point_ids]
end

function points_with_twoni(
    mat::MultiMaterial{ContinuumBasedMaterial},
    owned_points::Vector{UnitRange{Int}},
)
    return get_points_with_interaction(mat, owned_points, :C2)
end


function find_three_ni(
    pc::PointCloud,
    mat::PDMaterial{ContinuumBasedMaterial},
    one_ni_data::Vector{Tuple{Int,Int,Float64,Bool}},
    hood_range::Vector{UnitRange{Int}},
    owned_points::Vector{UnitRange{Int}},
)
    n_threads = nthreads()
    _three_ni_data = fill([(0, 0, 0, 0.0)], n_threads)
    n_three_ni_per_point = zeros(Int, pc.n_points)
    p = Progress(pc.n_points;
        dt=1,
        desc="Three-NI search...  ",
        barlen=30,
        color=:normal,
        enabled=!is_logging(stderr),
    )
    owned_points_with_threeni = points_with_threeni(mat, owned_points)
    @threads for tid in 1:n_threads
        local_three_ni_data = Vector{Tuple{Int,Int,Int,Float64}}(undef, 0)
        sizehint!(local_three_ni_data, pc.n_points * 1000)
        for i in owned_points_with_threeni[tid]
            num = 0
            δ = mat[i].δ
            jklseen = Set{Tuple{Int,Int,Int}}()
            for j_int in hood_range[i], k_int in hood_range[i], l_int in hood_range[i]
                _, j, _, _ = one_ni_data[j_int]
                _, k, _, _ = one_ni_data[k_int]
                _, l, _, _ = one_ni_data[l_int]
                if k !== j && l !== j && l !== k && !in((j, k, l), jklseen)
                    Ξijx = pc.position[1, j] - pc.position[1, i]
                    Ξijy = pc.position[2, j] - pc.position[2, i]
                    Ξijz = pc.position[3, j] - pc.position[3, i]
                    Ξikx = pc.position[1, k] - pc.position[1, i]
                    Ξiky = pc.position[2, k] - pc.position[2, i]
                    Ξikz = pc.position[3, k] - pc.position[3, i]
                    Ξilx = pc.position[1, l] - pc.position[1, i]
                    Ξily = pc.position[2, l] - pc.position[2, i]
                    Ξilz = pc.position[3, l] - pc.position[3, i]
                    Ξjkx = pc.position[1, k] - pc.position[1, j]
                    Ξjky = pc.position[2, k] - pc.position[2, j]
                    Ξjkz = pc.position[3, k] - pc.position[3, j]
                    Ξjlx = pc.position[1, l] - pc.position[1, j]
                    Ξjly = pc.position[2, l] - pc.position[2, j]
                    Ξjlz = pc.position[3, l] - pc.position[3, j]
                    Ξlkx = pc.position[1, k] - pc.position[1, l]
                    Ξlky = pc.position[2, k] - pc.position[2, l]
                    Ξlkz = pc.position[3, k] - pc.position[3, l]
                    _Ξjk = sqrt(Ξjkx * Ξjkx + Ξjky * Ξjky + Ξjkz * Ξjkz)
                    _Ξjl = sqrt(Ξjlx * Ξjlx + Ξjly * Ξjly + Ξjlz * Ξjlz)
                    _Ξlk = sqrt(Ξlkx * Ξlkx + Ξlky * Ξlky + Ξlkz * Ξlkz)
                    Aijkx = Ξijy * Ξikz - Ξijz * Ξiky
                    Aijky = Ξijz * Ξikx - Ξijx * Ξikz
                    Aijkz = Ξijx * Ξiky - Ξijy * Ξikx
                    volume = abs(Aijkx * Ξilx + Aijky * Ξily + Aijkz * Ξilz)
                    if volume > eps() && _Ξjk <= δ && _Ξjl <= δ && _Ξlk <= δ
                        num += 1
                        push!(local_three_ni_data, (j_int, k_int, l_int, volume))
                        push!(jklseen, (l, j, k))
                        push!(jklseen, (k, l, j))
                        push!(jklseen, (j, l, k))
                        push!(jklseen, (l, k, j))
                        push!(jklseen, (k, j, l))
                    end
                end
            end
            n_three_ni_per_point[i] = num
            next!(p)
        end
        _three_ni_data[tid] = local_three_ni_data
    end
    finish!(p)
    three_ni_data = reduce(append!, _three_ni_data)
    n_three_ni = length(three_ni_data)

    return three_ni_data, n_three_ni, n_three_ni_per_point
end

function points_with_threeni(::ContinuumBasedMaterial, owned_points::Vector{UnitRange{Int}})
    return owned_points
end

function points_with_threeni(
    mat::MultiMaterial{ContinuumBasedMaterial},
    owned_points::Vector{UnitRange{Int}},
)
    return get_points_with_interaction(mat, owned_points, :C3)
end

function surf_two_neigh(
    ξijx::Float64,
    ξijy::Float64,
    ξijz::Float64,
    ξikx::Float64,
    ξiky::Float64,
    ξikz::Float64,
)
    return sqrt(
        (ξijy * ξikz - ξijz * ξiky)^2 +
        (ξijz * ξikx - ξijx * ξikz)^2 +
        (ξijx * ξiky - ξijy * ξikx)^2
    )
end

create_simmodel(mat::PDMaterial{ContinuumBasedMaterial}, pc::PointCloud) =
    ContinuumBasedSimBody(mat, pc)


function compute_forcedensity!(
    body::ContinuumBasedSimBody,
    mat::PDMaterial{ContinuumBasedMaterial},
)
    body.b_int .= 0.0
    body.n_active_family_members .= 0
    if body.has_only_one_ni
        @threads for tid in 1:body.n_threads
            compute_forcedensity_one!(body, mat, tid)
        end
    elseif body.has_one_and_two_ni
        @threads for tid in 1:body.n_threads
            compute_forcedensity_one!(body, mat, tid)
            compute_forcedensity_two!(body, mat, tid)
        end
    elseif body.has_one_and_three_ni
        @threads for tid in 1:body.n_threads
            compute_forcedensity_one!(body, mat, tid)
            compute_forcedensity_three!(body, mat, tid)
        end
    elseif body.has_one_two_three_ni
        @threads for tid in 1:body.n_threads
            compute_forcedensity_one!(body, mat, tid)
            compute_forcedensity_two!(body, mat, tid)
            compute_forcedensity_three!(body, mat, tid)
        end
    end
    return nothing
end

function compute_forcedensity_one!(
    body::ContinuumBasedSimBody,
    mat::PDMaterial{ContinuumBasedMaterial},
    tid::Int,
)
    for one_ni in body.owned_bonds[tid]
        i, j, length_ref, failure_flag = body.bond_data[one_ni]
        ξijx = body.position[1, j] - body.position[1, i]
        ξijy = body.position[2, j] - body.position[2, i]
        ξijz = body.position[3, j] - body.position[3, i]
        length = sqrt(ξijx * ξijx + ξijy * ξijy + ξijz * ξijz)
        ε = (length - length_ref) / length_ref
        if ε > mat[i].εc || ε > mat[j].εc
            if failure_flag
                body.bond_failure[one_ni] = 0
            end
        end
        body.n_active_family_members[i, tid] += body.bond_failure[one_ni]
        temp = body.bond_failure[one_ni] * mat[i].C1 * (1 / length_ref - 1 / length) *
            body.volume_one_ni[i]
        body.b_int[1, i, tid] += temp * ξijx
        body.b_int[2, i, tid] += temp * ξijy
        body.b_int[3, i, tid] += temp * ξijz
    end
    return nothing
end

function compute_forcedensity_two!(
    body::ContinuumBasedSimBody,
    mat::PDMaterial{ContinuumBasedMaterial},
    tid::Int,
)
    for two_ni in body.owned_two_ni[tid]
        j_int, k_int, surface_ref = body.two_ni_data[two_ni]
        i, j, _ = body.bond_data[j_int]
        _, k, _ = body.bond_data[k_int]
        ξijx = body.position[1, j] - body.position[1, i]
        ξijy = body.position[2, j] - body.position[2, i]
        ξijz = body.position[3, j] - body.position[3, i]
        ξikx = body.position[1, k] - body.position[1, i]
        ξiky = body.position[2, k] - body.position[2, i]
        ξikz = body.position[3, k] - body.position[3, i]
        aijkx = ξijy * ξikz - ξijz * ξiky
        aijky = ξijz * ξikx - ξijx * ξikz
        aijkz = ξijx * ξiky - ξijy * ξikx
        surface = sqrt(aijkx * aijkx + aijky * aijky + aijkz * aijkz)
        if surface == 0 # to avoid divide by zero error for failed interactions
            surface = 1e-40
        end
        temp = body.bond_failure[j_int] * body.bond_failure[k_int] * 2 * mat[i].C2 *
            (1 / surface_ref - 1 / surface) * body.volume_two_ni[i]
        body.b_int[1, i, tid] += temp * (ξiky * aijkz - ξikz * aijky)
        body.b_int[2, i, tid] += temp * (ξikz * aijkx - ξikx * aijkz)
        body.b_int[3, i, tid] += temp * (ξikx * aijky - ξiky * aijkx)
        # variable switch: i,j = j,i
        body.b_int[1, i, tid] += temp * (ξijy * (ξikx * ξijy - ξiky * ξijx) -
            ξijz * (ξikz * ξijx - ξikx * ξijz))
        body.b_int[2, i, tid] += temp * (ξijz * (ξiky * ξijz - ξikz * ξijy) -
            ξijx * (ξikx * ξijy - ξiky * ξijx))
        body.b_int[3, i, tid] += temp * (ξijx * (ξikz * ξijx - ξikx * ξijz) -
            ξijy * (ξiky * ξijz - ξikz * ξijy))
    end
    return nothing
end

function compute_forcedensity_three!(
    body::ContinuumBasedSimBody,
    mat::PDMaterial{ContinuumBasedMaterial},
    tid::Int,
)
    for three_ni in body.owned_three_ni[tid]
        j_int, k_int, l_int, volume_ref = body.three_ni_data[three_ni]
        i, j, _ = body.bond_data[j_int]
        _, k, _ = body.bond_data[k_int]
        _, l, _ = body.bond_data[l_int]
        ξijx = body.position[1, j] - body.position[1, i]
        ξijy = body.position[2, j] - body.position[2, i]
        ξijz = body.position[3, j] - body.position[3, i]
        ξikx = body.position[1, k] - body.position[1, i]
        ξiky = body.position[2, k] - body.position[2, i]
        ξikz = body.position[3, k] - body.position[3, i]
        ξilx = body.position[1, l] - body.position[1, i]
        ξily = body.position[2, l] - body.position[2, i]
        ξilz = body.position[3, l] - body.position[3, i]
        # ijk
        aijkx = ξijy * ξikz - ξijz * ξiky
        aijky = ξijz * ξikx - ξijx * ξikz
        aijkz = ξijx * ξiky - ξijy * ξikx
        volume = aijkx * ξilx + aijky * ξily + aijkz * ξilz
        if volume == 0 # avoid to divide by zero error for failed interactions
            volume = 1e-40
        end
        abs_volume = abs(volume)
        temp = body.bond_failure[j_int] * body.bond_failure[k_int] *
            body.bond_failure[l_int] * 3 * mat[i].C3 * (1 / volume_ref - 1 / abs_volume) *
            volume * body.volume_three_ni[i]
        body.b_int[1, i, tid] += temp * (ξiky * ξilz - ξikz * ξily)
        body.b_int[2, i, tid] += temp * (ξikz * ξilx - ξikx * ξilz)
        body.b_int[3, i, tid] += temp * (ξikx * ξily - ξiky * ξilx)
        # kij  |  i->k, j->i, k->j  |  vkij == v
        body.b_int[1, i, tid] += temp * (ξijy * ξikz - ξijz * ξiky)
        body.b_int[2, i, tid] += temp * (ξijz * ξikx - ξijx * ξikz)
        body.b_int[3, i, tid] += temp * (ξijx * ξiky - ξijy * ξikx)
        # jki  |  i->j, j->k, k->i  |  vjki == v
        body.b_int[1, i, tid] += temp * (ξily * ξijz - ξilz * ξijy)
        body.b_int[2, i, tid] += temp * (ξilz * ξijx - ξilx * ξijz)
        body.b_int[3, i, tid] += temp * (ξilx * ξijy - ξily * ξijx)
        # ikj  |  i->i, j->k, k->j  |  vikj == -v
        body.b_int[1, i, tid] += temp * (ξily * ξikz - ξilz * ξiky)
        body.b_int[2, i, tid] += temp * (ξilz * ξikx - ξilx * ξikz)
        body.b_int[3, i, tid] += temp * (ξilx * ξiky - ξily * ξikx)
        # kji  |  i->k, j->j, k->i  |  vkji == -v
        body.b_int[1, i, tid] += temp * (ξiky * ξijz - ξikz * ξijy)
        body.b_int[2, i, tid] += temp * (ξikz * ξijx - ξikx * ξijz)
        body.b_int[3, i, tid] += temp * (ξikx * ξijy - ξiky * ξijx)
        # jik  |  i->j, j->i, k->k  |  vjik == -v
        body.b_int[1, i, tid] += temp * (ξijy * ξilz - ξijz * ξily)
        body.b_int[2, i, tid] += temp * (ξijz * ξilx - ξijx * ξilz)
        body.b_int[3, i, tid] += temp * (ξijx * ξily - ξijy * ξilx)
    end
    return nothing
end

function describe_mat(mat::ContinuumBasedMaterial)
    msg = @sprintf "  - Horizon δ [m]:                      %30g\n" mat.δ
    msg *= @sprintf "  - Density ρ [kg/m³]:                  %30g\n" mat.rho
    msg *= @sprintf "  - Young's modulus E [N/m²]:           %30g\n" mat.E
    msg *= @sprintf "  - Poisson ratio ν [-]:                %30g\n" mat.nu
    msg *= @sprintf "  - One-NI constant C₁ [N/m⁶]:          %30g\n" mat.C1
    msg *= @sprintf "  - Two-NI constant C₂ [N/m¹⁰]:         %30g\n" mat.C2
    msg *= @sprintf "  - Three-NI constant C₃ [N/m¹⁴]:       %30g\n" mat.C3
    msg *= @sprintf "  - Griffith's Parameter Gc [N/m]:      %30g\n" mat.Gc
    msg *= @sprintf "  - Critical bond stretch [-]:          %30g\n" mat.εc
    return msg
end

function describe_interactions(body::ContinuumBasedSimBody)
    msg = @sprintf "  - Number of bonds [-]:                %30d\n" body.n_bonds
    msg *= @sprintf "  - Number of two-neighbor int. [-]:    %30d\n" body.n_two_ni
    msg *= @sprintf "  - Number of three-neighbor int. [-]:  %30d\n" body.n_three_ni
    return msg
end
