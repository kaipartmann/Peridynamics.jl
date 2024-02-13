function get_frac_params(p::Dict{Symbol,Any}, δ::Float64, K::Float64)
    local Gc::Float64
    local εc::Float64

    if haskey(p, :Gc) && !haskey(p, :epsilon_c)
        Gc = float(p[:Gc])
        εc = sqrt(5.0 * Gc / (9.0 * K * δ))
    elseif !haskey(p, :Gc) && haskey(p, :epsilon_c)
        εc = float(p[:epsilon_c])
        Gc = 9.0 / 5.0 * K * δ * εc^2
    elseif haskey(p, :Gc) && haskey(p, :epsilon_c)
        msg = "insufficient keywords for calculation of fracture parameters!\n"
        msg *= "Define either Gc or epsilon_c, not both!\n"
        throw(ArgumentError(msg))
    else
        msg = "insufficient keywords for calculation of fracture parameters!\n"
        msg *= "Define either Gc or epsilon_c!\n"
        throw(ArgumentError(msg))
    end

    return Gc, εc
end
