
@testitem "show DynamicRelaxation" begin
    io = IOBuffer()

    dr = DynamicRelaxation(steps=1)

    show(IOContext(io, :compact=>true), MIME("text/plain"), dr)
    msg = String(take!(io))
    @test contains(msg, "DynamicRelaxation(n_steps=1, Δt=1.0, Λ=1.0)")

    show(IOContext(io, :compact=>false), MIME("text/plain"), dr)
    msg = String(take!(io))
    @test contains(msg, "DynamicRelaxation:\n  n_steps  1\n  Δt       1\n  Λ        1")
end
