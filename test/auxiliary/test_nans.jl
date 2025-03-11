@testitem "containsnan" begin
    A1 = [1.0, 2.0, NaN, 4.0]
    A2 = [0; 1; 2;; 3; 4; 5;; 6; 7; NaN]
    A3 = [A2;;; 2A2;;; 3A2]
    for A in [A1, A2, A3]
        @test Peridynamics.containsnan(A) == true
        a = replace(A, NaN=>0.0)
        @test Peridynamics.containsnan(a) == false
    end
end

@testitem "nancheck" begin
    ref_position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    volume = [1.0, 1.0]
    δ = 1.5
    E = 210e9
    body = Body(BBMaterial(), ref_position, volume)
    material!(body, horizon=δ, rho=1, E=E)
    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; b_int) = chunk.storage

    @test Peridynamics.nancheck(chunk, 0.0) === nothing

    b_int[3, end] = NaN
    @test_throws ErrorException Peridynamics.nancheck(chunk, 0.0)
end
