@testitem "InterfaceError" begin # lint-ok: println is the example method of the error
    using Peridynamics.Printf
    ie1 = Peridynamics.InterfaceError(Float64, println)
    msg_correct = "interface method not correctly defined!"
    msg_correct *= "\n  type:    Float64\n  method:  println\n"
    @test sprint(showerror, ie1) == msg_correct
end

@testitem "NaNError" begin
    using Peridynamics.Printf
    ne1 = Peridynamics.NaNError(2.0, 2)
    msg_correct = "NaN values detected in force density field!\n"
    msg_correct *= "  time:    2.0\n  step:    2\n"
    @test sprint(showerror, ne1) == msg_correct
end

@testitem "InterfaceError: an optional hint is printed after the method" begin
    ie = Peridynamics.InterfaceError(Float64, "foo", "define `foo(::Float64)`")
    msg = sprint(showerror, ie)
    @test contains(msg, "type:    Float64")
    @test contains(msg, "method:  foo")
    @test endswith(strip(msg), "define `foo(::Float64)`")
    @test isempty(Peridynamics.InterfaceError(Float64, "foo").hint)
end

@testitem "StorageContractError: names every missing field and the reason" begin
    S = Peridynamics.storage_type(BBMaterial())
    missing_fields = [:residual => "the time solver `NewtonKrylov`",
                      :my_field => "the material `BBMaterial`"]
    err = Peridynamics.StorageContractError(S, typeof(BBMaterial()), NewtonKrylov, missing_fields)
    msg = sprint(showerror, err)
    @test contains(msg, "storage type does not contain all required fields!")
    @test contains(msg, "BBStorage")
    @test contains(msg, "residual   required by the time solver `NewtonKrylov`")
    @test contains(msg, "my_field   required by the material `BBMaterial`")
    @test contains(msg, "Peridynamics.@storage")
end

