struct UndefPointSetHandler <: AbstractPointSetHandler end

struct NoMatPointSetHandler <: AbstractPointSetHandler
    point_sets::Dict{Symbol,AbstractVector}
    conditions::Dict{Symbol,AbstractCondition}
    precracks::Dict{Symbol,AbstractPredefinedCrack}
end

function NoMatPointSetHandler()
    point_sets = Dict{Symbol,AbstractVector}()
    conditions = Dict{Symbol,AbstractCondition}()
    precracks = Dict{Symbol,AbstractPredefinedCrack}()
    NoMatPointSetHandler(point_sets, conditions, precracks)
end

struct PointSetHandler{M<:AbstractMaterial} <: AbstractPointSetHandler
    point_sets::Dict{Symbol,AbstractVector}
    materials::Dict{Symbol,M}
    conditions::Dict{Symbol,AbstractCondition}
    precracks::Dict{Symbol,AbstractPredefinedCrack}

    function PointSetHandler(::Type{M},
                             nmpsh::NoMatPointSetHandler) where {M<:AbstractMaterial}
        materials = Dict{Symbol,M}()
        return new{M}(nmpsh.point_sets, materials, nmpsh.conditions, nmpsh.precracks)
    end
end

function PointSetHandler(::Type{M}, ::UndefPointSetHandler) where {M<:AbstractMaterial}
    return PointSetHandler(M, NoMatPointSetHandler())
end

function PointSetHandler(::Type{M}) where {M<:AbstractMaterial}
    return PointSetHandler(M, NoMatPointSetHandler())
end

function _point_set!(::UndefPointSetHandler, ::Symbol, ::AbstractVector)
    error("trying to assign to a empty point set handler! Something is wrong!")
end

function _point_set!(psh::AbstractPointSetHandler, name::Symbol,
                     points::V) where {V<:AbstractVector}
    if haskey(psh.point_sets, name)
        error("there is already a point set with name $(name)!")
    end
    psh.point_sets[name] = points
    return nothing
end

function _material!(psh::PointSetHandler{M1}, name::Symbol,
                    mat::M2) where {M1,M2<:AbstractMaterial}
    if M1 != M2
        error("wrong material type $(M2)! This body should be used with $(M1)")
    end
    if !haskey(psh.point_sets, name)
        error("there is no point set with name $(name)!")
    end
    if haskey(psh.materials, name)
        error("material for set $(name) already defined!")
    end
    psh.materials[name] = mat
    return nothing
end
