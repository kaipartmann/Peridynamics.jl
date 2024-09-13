function system_type(mat::AbstractMaterial)
    throw(MethodError(system_type, mat))
end

# function get_system_type(mat::AbstractMaterial)
#     Sys = system_type(mat)
#     FullSys = full_system_type(Sys, mat)
#     return FullSys
# end

# function full_system_type(::Type{System}, mat::AbstractMaterial) where {System}
#     return system_type(mat)
# end

function get_system(::AbstractBody{M}, ::PointDecomposition, ::Int) where {M}
    msg = "system for material $M not specified!\n"
    return error(msg)
end

function log_system(options::AbstractJobOptions, dh::AbstractDataHandler)
    log_system(system_type(dh), options, dh)
    return nothing
end

# macro system(material, system)
#     macrocheck_input_material(material)
#     macrocheck_input_system(system)
#     local _checks = quote
#         Peridynamics.typecheck_material($(esc(material)))
#         Peridynamics.typecheck_system($(esc(system)))
#     end
#     local _system_type = quote
#         Peridynamics.system_type(::$(esc(material))) = $(esc(system))
#     end
#     local _get_system = quote
#         function Peridynamics.get_system(body::Peridynamics.AbstractBody{$(esc(material))},
#                                          args...)
#             return $(esc(system))(body, args...)
#         end
#     end
#     return Expr(:block, _checks, _system_type, _get_system)
# end

# function macrocheck_input_system(system)
#     system isa Symbol && return nothing
#     (system isa Expr && system.head === :.) && return nothing
#     return throw(ArgumentError("argument `$system` is not a valid system input!\n"))
# end

# function typecheck_system(::Type{System}) where {System}
#     if !(System <: AbstractSystem)
#         msg = "$System is not a valid system type!\n"
#         throw(ArgumentError(msg))
#     end
#     return nothing
# end
