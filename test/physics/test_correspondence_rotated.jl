@testitem "CRMaterial" begin
    mat = CRMaterial()
    @test mat.constitutive_model isa LinearElastic

    @test_throws ArgumentError CRMaterial(model=NeoHookePenalty())

    io = IOBuffer()
    show(IOContext(io, :compact=>true), MIME("text/plain"), mat)
    msg = String(take!(io))
    @test startswith(msg, "CRMaterial")
end
