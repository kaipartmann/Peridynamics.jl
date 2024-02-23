
function uniform_box(lx::Real, ly::Real, lz::Real, Δx::Real;
                     center_x::Real=0, center_y::Real=0, center_z::Real=0)
    _gridx = range((-lx + Δx) / 2, (lx - Δx) / 2; step=Δx)
    gridx = _gridx .- sum(_gridx) / length(_gridx)
    _gridy = range((-ly + Δx) / 2, (ly - Δx) / 2; step=Δx)
    gridy = _gridy .- sum(_gridy) / length(_gridy)
    _gridz = range((-lz + Δx) / 2, (lz - Δx) / 2; step=Δx)
    gridz = _gridz .- sum(_gridz) / length(_gridz)
    _position = vec(collect(Iterators.product(gridx, gridy, gridz)))
    position = copy(reinterpret(reshape, eltype(eltype(_position)), _position))
    volume = fill(Δx^3, size(position, 2))
    center_x != 0 && (position[1, :] .+= center_x)
    center_y != 0 && (position[2, :] .+= center_y)
    center_z != 0 && (position[3, :] .+= center_z)
    return position, volume
end
