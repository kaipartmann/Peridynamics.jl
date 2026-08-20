# The exported strain energy density of every material under homogeneous deformations, checked
# against the closed forms of the continuum it is calibrated to. One item per material family;
# adding a (material, model) pair to `Fixtures.MODEL_MATERIALS` covers it here automatically.

@testmodule StrainEnergy begin
    using Peridynamics, Test
    using Peridynamics.StaticArrays

    mean(x) = sum(x) / length(x)

    # indices of the points of the unit cube further than one horizon from every surface;
    # points near a surface have incomplete families and a stiffness below the bulk value, so
    # only the interior error reflects the calibration of the bond constants
    function interior_ids(body, system)
        δ = body.point_params[1].δ
        X = system.position
        return [i for i in axes(X, 2) if all(abs(X[d, i]) < 0.5 - δ for d in 1:3)]
    end

    # the unit cube of `mat` with spacing `Δx`, `material!` called with `kwargs`
    function cube(mat, Δx; kwargs...)
        pos, vol = uniform_box(1, 1, 1, Δx)
        body = Body(mat, pos, vol)
        material!(body; kwargs...)
        return body
    end

    # check the exported strain energy density of `body` under `F_a` against the closed form
    # `Ψ_a`; `tols` is the two-sided band on the relative error over *all* points (dominated by
    # the surface effect), `interior_tols` the band on the interior points alone, where the
    # error is a deterministic constant and can be asserted tightly
    function check(body, F_a, Ψ_a, tols; interior_tols=nothing, atol0=eps())
        ts = VelocityVerlet(steps=1)
        dh = Peridynamics.threads_data_handler(body, ts, 1)
        Peridynamics.initialize!(dh, ts)
        (; mat, system, storage, paramsetup) = dh.chunks[1]
        field = Val{:strain_energy_density}()

        # no deformation -> no strain energy density
        storage.position .= system.position
        storage.strain_energy_density .= 0.0
        Ψ_0 = Peridynamics.export_field(field, mat, system, storage, paramsetup, 0.0)
        @test all(isapprox.(Ψ_0, 0.0; atol=atol0))

        # apply deformation gradient F_a
        for i in Peridynamics.each_point_idx(system)
            Xi = Peridynamics.get_vector(system.position, i)
            Peridynamics.update_vector!(storage.position, i, F_a * Xi)
        end
        Peridynamics.calc_force_density!(dh, 0.0, 0.0)

        Ψ_pd = Peridynamics.export_field(field, mat, system, storage, paramsetup, 0.0)
        ΔΨ = (Ψ_pd .- Ψ_a) ./ Ψ_a
        @debug "strain energy density" mat Ψ_a mean(Ψ_pd) extrema(ΔΨ)
        @test minimum(ΔΨ) > tols[1]
        @test maximum(ΔΨ) < tols[2]

        if !isnothing(interior_tols)
            ids = interior_ids(body, system)
            # a body without interior points would make the assertions below vacuous
            @test !isempty(ids)
            ΔΨ_int = ΔΨ[ids]
            @debug "strain energy density, interior" length(ids) extrema(ΔΨ_int)
            @test minimum(ΔΨ_int) > interior_tols[1]
            @test maximum(ΔΨ_int) < interior_tols[2]
        end
        return nothing
    end

    # the deformations: λ uniform extension (ε = λ - 1 the strain), β shear
    const λ = 1.1
    const ε = λ - 1
    const β = 0.1
    F_iso(λ=λ) = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    F_shear(β=β) = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    F_x(λ=λ) = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    F_y(λ=λ) = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    F_z(λ=λ) = @SMatrix [1 0 0; 0 1 0; 0 0 λ]

    # closed forms of the small strain continuum for the bond-based families, as functions of
    # the point parameters
    Ψ_small_iso(p) = 9/2 * p.K * ε^2
    Ψ_small_shear(p) = 1/2 * p.G * β^2
    Ψ_small_uniaxial_bb(p) = 3/5 * p.E * ε^2          # bond-based: nu = 1/4
    Ψ_small_uniaxial(p) = 1/2 * p.λ * ε^2 + p.μ * ε^2 # general isotropic

    # closed forms of the Saint-Venant-Kirchhoff continuum, which every constitutive model
    # reproduces at small strains up to its own nonlinearity
    Ψ_svk_iso(p, λ=λ) = 9/8 * p.K * (λ^2 - 1)^2
    Ψ_svk_shear(p, β=β) = (p.μ / 2) * β^2 + (p.λ / 8 + p.μ / 4) * β^4
    Ψ_svk_uniaxial(p, λ=λ) = (p.λ + 2 * p.μ) / 8 * (λ^2 - 1)^2
end

@testitem "strain energy density: bond-based materials" setup=[Fixtures, StrainEnergy] begin
    # Δx = 1/12 with a horizon of 3.015Δx leaves 216 of 1728 points more than one horizon from
    # every surface, so the bond constant can be checked away from the surface effect. Over all
    # points the error spans -76% to +10% and no tolerance in that range says anything. On the
    # interior it is a single number, identical for every interior point and unchanged between
    # Δx = 1/12 and Δx = 1/16: the lattice quadrature error of the bond constant at m = 3.015,
    # which `BBMaterial` carries because it does not normalize by the weighted volume of the
    # actual family. `DHBBMaterial` reproduces the same lattice constants exactly.
    S = StrainEnergy
    for (name, mat) in ("BB" => BBMaterial(), "DHBB" => DHBBMaterial())
        @testset "$name" begin
            Δx = 1 / 12
            body = S.cube(mat, Δx; horizon=3.015Δx, rho=8000, E=210e9)
            p = body.point_params[1]
            tols = (-0.9, 0.3)
            S.check(body, S.F_iso(), S.Ψ_small_iso(p), tols; interior_tols=(0.094, 0.104))
            S.check(body, S.F_shear(), S.Ψ_small_shear(p), tols; interior_tols=(0.191, 0.201))
            for F in (S.F_x(), S.F_y(), S.F_z())
                S.check(body, F, S.Ψ_small_uniaxial_bb(p), tols; interior_tols=(0.063, 0.073))
            end
        end
    end
end

@testitem "strain energy density: GBBMaterial and OSBMaterial" setup=[Fixtures, StrainEnergy] begin
    # both are exact for a homogeneous isotropic extension and carry larger errors in shear
    S = StrainEnergy
    Δx = 0.2
    @testset "GBB" begin
        body = S.cube(GBBMaterial(), Δx; horizon=3.01Δx, rho=8000, E=210e9)
        p = body.point_params[1]
        S.check(body, S.F_iso(), S.Ψ_small_iso(p), (-1e-10, 1e-10))
        S.check(body, S.F_shear(), S.Ψ_small_shear(p), (-0.4, 0.4))
        for F in (S.F_x(), S.F_y(), S.F_z())
            S.check(body, F, S.Ψ_small_uniaxial_bb(p), (-0.9, 0.3))
        end
    end
    @testset "OSB" begin
        body = S.cube(OSBMaterial(), Δx; horizon=3.01Δx, rho=8000, E=210e9, nu=0.25)
        p = body.point_params[1]
        S.check(body, S.F_iso(), S.Ψ_small_iso(p), (-1e-10, 1e-10))
        S.check(body, S.F_shear(), S.Ψ_small_shear(p), (-0.4, 0.4))
        for F in (S.F_x(), S.F_y(), S.F_z())
            S.check(body, F, S.Ψ_small_uniaxial(p), (-0.9, 0.3))
        end
    end
end

@testitem "strain energy density: CKIMaterial" setup=[Fixtures, StrainEnergy] begin
    # `CKIMaterial` has its own lattice constants, close to but not equal to the bond-based
    # ones; see the bond-based item for the resolution and the interior tolerances
    S = StrainEnergy
    Δx = 1 / 12
    body = S.cube(CKIMaterial(), Δx; horizon=3.015Δx, rho=8000, E=210e9, nu=0.25)
    p = body.point_params[1]
    S.check(body, S.F_iso(), S.Ψ_small_iso(p), (-0.8, 0.2); interior_tols=(0.103, 0.113))
    S.check(body, S.F_shear(), S.Ψ_small_shear(p), (-0.9, 0.3); interior_tols=(0.201, 0.211))
    for F in (S.F_x(), S.F_y(), S.F_z())
        S.check(body, F, S.Ψ_small_uniaxial(p), (-0.9, 0.3); interior_tols=(0.072, 0.082))
    end
end

@testitem "strain energy density: CKIMaterial with explicit C1, C2, C3" setup=[Fixtures, StrainEnergy] begin
    # With `nu = 0.25` the derived `C3` is zero and no three-neighbor interactions exist (see
    # `Fixtures.cki_kwargs`), so here the interaction parameters are given and the Lamé
    # parameters are derived from them; this is the only way to reach the two- and
    # three-neighbor force loops. The triplet search scales with the ninth power of the
    # horizon-to-spacing ratio, hence `2.01Δx`.
    S = StrainEnergy
    Δx = 0.2
    δ = 2.01Δx
    C1, C2, C3 = 210e9, 80e9, 50e9
    λlame = π * δ^4 / 30 * C1 + π^3 * δ^8 / 30 * C2 + π^4 * δ^12 / 32 * C3
    μlame = π * δ^4 / 30 * C1 + π^3 * δ^8 / 180 * C2
    pos, vol = uniform_box(1, 1, 1, Δx)
    body = Body(CKIMaterial(), pos, vol)
    @test_logs (:warn, r"specified manually") material!(body; horizon=δ, rho=8000, lambda=λlame,
                                                          mu=μlame, C1, C2, C3)
    p = body.point_params[1]
    # The 5³ cube has no point further than one horizon from every surface, so these are
    # bands on a surface-dominated error (measured: -0.82 to +0.21 at m = 2.01) that catch a
    # broken interaction loop, not the calibration; the calibration is checked on the interior
    # in the item above.
    S.check(body, S.F_iso(), S.Ψ_small_iso(p), (-0.9, 0.3))
    S.check(body, S.F_shear(), S.Ψ_small_shear(p), (-0.9, 0.3))
    for F in (S.F_x(), S.F_y(), S.F_z())
        S.check(body, F, S.Ψ_small_uniaxial(p), (-0.9, 0.3))
    end
end

@testitem "strain energy density: correspondence materials" tags=[:slow] setup=[Fixtures, StrainEnergy] begin
    # `:slow`: 12 material/model combinations, each compiled once; runs in the `extras` CI job
    # The correspondence formulations (classic, rotated, reproducing kernel, bond-associated)
    # reproduce the strain energy density of their constitutive model exactly for a homogeneous
    # deformation, because the deformation gradient is exact there. The Saint-Venant-Kirchhoff
    # and linear elastic models agree with the closed forms to rounding; the Neo-Hooke models
    # differ by their nonlinearity at these strains.
    S = StrainEnergy
    Δx = 0.2
    for (name, mat) in Fixtures.MODEL_MATERIALS
        # the bond-associated formulation does not export a strain energy density
        mat isa BACMaterial && continue
        @testset "$name" begin
            body = S.cube(mat, Δx; horizon=3.01Δx, rho=8000, E=210e9, nu=0.25)
            p = body.point_params[1]
            if mat.constitutive_model isa Union{NeoHooke,NeoHookePenalty}
                # small strains, where the Neo-Hooke models are close to Saint-Venant-Kirchhoff;
                # their energy at the reference configuration is zero only up to rounding of
                # the large moduli, hence `atol0`
                λ, β, tols, atol0 = 1.001, 0.05, (-0.02, 0.02), 1e-4
            else
                λ, β, tols, atol0 = S.λ, S.β, (-1e-12, 1e-12), eps()
            end
            S.check(body, S.F_iso(λ), S.Ψ_svk_iso(p, λ), tols; atol0)
            S.check(body, S.F_shear(β), S.Ψ_svk_shear(p, β), tols; atol0)
            for F in (S.F_x(λ), S.F_y(λ), S.F_z(λ))
                S.check(body, F, S.Ψ_svk_uniaxial(p, λ), tols; atol0)
            end
        end
    end
end
