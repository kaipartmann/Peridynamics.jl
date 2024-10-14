@testitem "InterfaceError" begin
    ie1 = Peridynamics.InterfaceError(Float64, println)
    io = IOBuffer()
    Base.showerror(IOContext(io, :compact=>false), ie1)
    msg = String(take!(io))
    msg_correct = "interface method not correctly defined!"
    msg_correct *= "\n  type:    Float64\n  method:  println\n"
    @test msg == msg_correct
end
