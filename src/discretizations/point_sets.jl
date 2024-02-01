
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

function material_is_undefined(psh::PointSetHandler{M}) where {M<:AbstractMaterial}
    isempty(psh.materials) && !isconcretetype(M) && return true
    return false
end

@inline material_type(::PointSetHandler{M}) where {M<:AbstractMaterial} = M

function point_set_handler_is_undefined(psh::PointSetHandler{M}) where {M<:AbstractMaterial}
    isempty(psh.point_sets) && isempty(psh.materials) && return true
    return false
end

function check_if_set_is_defined(psh::PointSetHandler{M}, name::Symbol) where {M}
    if !haskey(psh.point_sets, name)
        error("there is no point set with name $(name)!")
    end
    return nothing
end

function _point_set!(psh::PointSetHandler{M}, name::Symbol,
                     points::V) where {V<:AbstractVector,M}
    if haskey(psh.point_sets, name)
        error("there is already a point set with name $(name)!")
    end
    psh.point_sets[name] = points
    return nothing
end

function _material!(psh::PointSetHandler{M1}, name::Symbol,
                    mat::M2) where {M1,M2<:AbstractMaterial}
    if M1 != M2
        msg = "body is already assigned with material type $(M1))), "
        msg *= "cannot add material with type $(M2)!\n"
        throw(ArgumentError(msg))
    end
    if haskey(psh.materials, name)
        msg = "material for set $(name) already defined!\n"
        throw(ArgumentError(msg))
    end
    psh.materials[name] = mat
    return nothing
end
