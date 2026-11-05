# Test const_one_kernel
@testitem "const_one_kernel returns 1 for any input" begin
    @test const_one_kernel(1.0, 2.0) == 1
    @test const_one_kernel(0.5, 0.1) == 1
    @test const_one_kernel(100.0, 0.001) == 1
end

# Test linear_kernel
@testitem "linear_kernel computes δ / L" begin
    @test linear_kernel(2.0, 1.0) ≈ 2.0
    @test linear_kernel(1.0, 2.0) ≈ 0.5
    @test linear_kernel(0.0, 1.0) ≈ 0.0
end

# Test cubic_b_spline_kernel
@testitem "cubic_b_spline_kernel basic values and support" begin
    δ = 1.0
    @test cubic_b_spline_kernel(δ, 0.0) ≈ 0 atol=eps()
    @test cubic_b_spline_kernel(δ, δ) ≈ 0 atol=eps()
    @test cubic_b_spline_kernel(δ, 0.5δ) ≈ (2/3 - 4*(0.5)^2 + 4*(0.5)^3)
    @test cubic_b_spline_kernel(δ, 0.25δ) ≈ (2/3 - 4*(0.25)^2 + 4*(0.25)^3)
    @test cubic_b_spline_kernel(δ, 0.75δ) ≈ (4/3 - 4*0.75 + 4*0.75^2 - 4/3*0.75^3)
end

# Test cubic_b_spline_kernel_norm
@testitem "cubic_b_spline_kernel_norm basic values and support" begin
    δ = 1.0
    C = 8/(π * δ^3)
    @test cubic_b_spline_kernel_norm(δ, 0.0) ≈ 0 atol=eps()
    @test cubic_b_spline_kernel_norm(δ, δ) ≈ 0 atol=eps()
    @test cubic_b_spline_kernel_norm(δ, 0.5δ) ≈ C * (1 - 6*(0.5)^2 + 6*(0.5)^3)
    @test cubic_b_spline_kernel_norm(δ, 0.25δ) ≈ C * (1 - 6*(0.25)^2 + 6*(0.25)^3)
    @test cubic_b_spline_kernel_norm(δ, 0.75δ) ≈ C * 2 * (1 - 0.75)^3
end
