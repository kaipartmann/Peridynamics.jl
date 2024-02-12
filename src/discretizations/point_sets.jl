
struct PointSetHandler{M<:AbstractMaterial}
    point_sets::Dict{Symbol,AbstractVector{<:Integer}}
    materials::Dict{Symbol,M}
end

function PointSetHandler(::Type{M}) where {M<:AbstractMaterial}
    point_sets = Dict{Symbol,AbstractVector{<:Integer}}()
    materials = Dict{Symbol,M}()
    return PointSetHandler{M}(point_sets, materials)
end

function PointSetHandler(::Type{M}, psh::PointSetHandler) where {M<:AbstractMaterial}
    return PointSetHandler{M}(psh.point_sets, Dict{Symbol,M}())
end

function material_is_undefined(materials::Dict{Symbol,M}) where {M<:AbstractMaterial}
    isempty(materials) && !isconcretetype(M) && return true
    return false
end

function check_if_set_is_defined(point_sets::Dict{Symbol,V}, name::Symbol) where {V}
    if !haskey(point_sets, name)
        error("there is no point set with name $(name)!")
    end
    return nothing
end

function _point_set!(point_sets::Dict{Symbol,V}, name::Symbol,
                     points::V) where {V<:AbstractVector{<:Integer}}
    if haskey(point_sets, name)
        error("there is already a point set with name $(name)!")
    end
    point_sets[name] = points
    return nothing
end

function _material!(materials::Dict{Symbol,M1}, name::Symbol,
                    mat::M2) where {M1,M2<:AbstractMaterial}
    if M1 != M2
        msg = "body is already assigned with material type $(M1))), "
        msg *= "cannot add material with type $(M2)!\n"
        throw(ArgumentError(msg))
    end
    if haskey(materials, name)
        msg = "material for set $(name) already defined!\n"
        throw(ArgumentError(msg))
    end
    materials[name] = mat
    return nothing
end
