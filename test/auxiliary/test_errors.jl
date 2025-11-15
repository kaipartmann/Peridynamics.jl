@testitem "InterfaceError" begin
    using Peridynamics.Printf
    ie1 = Peridynamics.InterfaceError(Float64, println)
    msg_correct = "interface method not correctly defined!"
    msg_correct *= "\n  type:    Float64\n  method:  println\n"
    @test sprint(showerror, ie1) == msg_correct
end

@testitem "NaNError" begin
    using Peridynamics.Printf
    ne1 = Peridynamics.NaNError(2.0, 2)
    msg_correct = "NaN values detected in simulation data!\n  time:    2.0\n  step:    2\n"
    @test sprint(showerror, ne1) == msg_correct
end
