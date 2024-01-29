
@enum HaloExchangeType read_from_halo write_to_halo

struct ThreadsHaloExchange
    type::HaloExchangeType
    field::Symbol
    source::Int
    dest::Int
    send_from_to::Dict{Int,Int}
end

function get_halo_exchange_info(::T) where {T<:Material}
    @warn "material should specify the function `get_halo_exchange_info`!"
    return Dict{Symbol,HaloExchangeType}()
end
