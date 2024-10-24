function von_mises_stress(σ)
    σx, σy, σz = σ[1], σ[5], σ[9]
    τxy, τxz, τyz = σ[4], σ[7], σ[8]
    a = σx * σx + σy * σy + σz * σz
    b = - σx * σy - σx * σz - σy * σz
    c = 3 * (τxy * τxy + τxz * τxz + τyz * τyz)
    σvm = √(a + b + c)
    return σvm
end
