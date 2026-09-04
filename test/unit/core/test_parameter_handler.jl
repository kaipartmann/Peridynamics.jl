# `ParameterHandler`: the per-point mapping to the parameter sets of a body with more than one
# `material!` call (`src/core/parameter_handler.jl`).

@testitem "get_param_spec: one or several parameter sets" setup=[BodyCase] begin
    body = BodyCase.tetra()
    material!(body; horizon=2, E=1, rho=1, Gc=1)
    @test Peridynamics.get_param_spec(body) isa Peridynamics.SingleParamChunk
    @test Peridynamics.parameter_setup_type(body, Peridynamics.SingleParamChunk()) ===
          Peridynamics.point_param_type(body.mat)
    point_set!(body, :a, 1:2)
    material!(body, :a; horizon=2, E=2, rho=2, Gc=2)
    @test Peridynamics.get_param_spec(body) isa Peridynamics.MultiParamChunk
    @test Peridynamics.parameter_setup_type(body, Peridynamics.MultiParamChunk()) ===
          Peridynamics.ParameterHandler{Peridynamics.point_param_type(body.mat),
                                        Vector{Peridynamics.point_param_type(body.mat)},
                                        Vector{Int}}
end

@testitem "ParameterHandler: maps the local points of a chunk to their parameters" setup=[Fixtures, BodyCase] begin
    body = BodyCase.tetra()
    point_set!(body, :a, [1, 4])
    material!(body; horizon=2, E=1, rho=1, Gc=1)
    material!(body, :a; horizon=2, E=2, rho=2, Gc=2)
    # a single chunk: the parameter setup is a handler with the map of the body
    chunk = Fixtures.chunk(body, VelocityVerlet(steps=1))
    (; paramsetup) = chunk
    @test paramsetup isa Peridynamics.ParameterHandler
    @test paramsetup.parameters === body.point_params
    @test paramsetup.point_mapping == body.params_map == [2, 1, 1, 2]
    @test Peridynamics.get_params(paramsetup, 1).E == 2.0
    @test Peridynamics.get_params(paramsetup, 2).E == 1.0
    @test Peridynamics.get_params(chunk, 4).E == 2.0
    # a single parameter set is its own setup
    p = body.point_params[1]
    @test Peridynamics.get_params(p, 3) === p
    # two chunks: the mapping is localized to the points of each chunk (incl. halo)
    dh = Fixtures.handler(body, VelocityVerlet(steps=1); n_chunks=2, init=false)
    for chunk in dh.chunks
        ids = chunk.system.chunk_handler.point_ids
        @test chunk.paramsetup.point_mapping == body.params_map[ids]
    end
end

@testitem "ParameterHandler: both containers are type parameters and move with Adapt" setup=[Fixtures, BodyCase] begin
    # a minimal stand-in for the array type of another backend, see `test_storage_fields.jl`
    struct WrappedArray{T,N} <: AbstractArray{T,N}
        a::Array{T,N}
    end
    Base.size(x::WrappedArray) = size(x.a)
    Base.getindex(x::WrappedArray, i...) = getindex(x.a, i...)
    Base.setindex!(x::WrappedArray, v, i...) = setindex!(x.a, v, i...)
    struct WrappedBackend end
    function Peridynamics.Adapt.adapt_storage(::WrappedBackend, a::Array{T,N}) where {T,N}
        return WrappedArray{T,N}(a)
    end

    body = BodyCase.tetra()
    point_set!(body, :a, [1, 4])
    material!(body; horizon=2, E=1, rho=1, Gc=1)
    material!(body, :a; horizon=2, E=2, rho=2, Gc=2)
    chunk = Fixtures.chunk(body, VelocityVerlet(steps=1))
    (; paramsetup) = chunk
    P = Peridynamics.point_param_type(body.mat)
    @test paramsetup isa Peridynamics.ParameterHandler{P,Vector{P},Vector{Int}}

    # both containers move to the other backend and the parameters themselves are unchanged
    adapted = Peridynamics.Adapt.adapt(WrappedBackend(), paramsetup)
    @test adapted isa Peridynamics.ParameterHandler{P}
    @test adapted.parameters isa WrappedArray{P,1}
    @test adapted.point_mapping isa WrappedArray{Int,1}
    @test adapted.point_mapping == paramsetup.point_mapping
    @test Peridynamics.get_params(adapted, 1).E == 2.0

    # a target that moves no array hands back the handler it was given
    @test Peridynamics.Adapt.adapt(Array, paramsetup) === paramsetup
end
