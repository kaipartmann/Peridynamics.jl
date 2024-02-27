
@testitem "Aqua.jl" begin
    using Aqua
    Aqua.test_all(Peridynamics; stale_deps=(ignore=[:ThreadPinning],))
end
