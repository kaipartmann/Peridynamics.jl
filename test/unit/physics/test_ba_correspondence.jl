# `BACMaterial`: the bond-associated correspondence formulation of
# `src/physics/ba_correspondence.jl`. The bond-associated system and its invariants are tested
# in `test/unit/discretization/test_bond_associated_system.jl`, the force density invariants for all
# materials in `test/materials/`.

@testitem "BACMaterial: constructor and show" begin
    mat = BACMaterial()
    @test mat.kernel === linear_kernel
    @test mat.constitutive_model isa SaintVenantKirchhoff
    @test mat.dmgmodel isa CriticalStretch
    @test mat.maxdmg == 0.85
    mat2 = BACMaterial(kernel=cubic_b_spline_kernel, model=NeoHooke(), maxdmg=0.5)
    @test mat2.kernel === cubic_b_spline_kernel
    @test mat2.constitutive_model isa NeoHooke
    @test mat2.maxdmg == 0.5
    msg = sprint(show, mat2)
    @test contains(msg, "BACMaterial")
    @test contains(msg, "maxdmg=0.5")
    # the bond horizon is a material parameter of the bond-associated system
    @test :bond_horizon in Peridynamics.allowed_material_kwargs(mat)
end

@testitem "BACMaterial: bond-associated damage kills the family beyond maxdmg" setup=[Fixtures] begin
    # every bond of a point is deactivated once the point is damaged beyond `maxdmg`, so that
    # the deformation gradient of a family with too few bonds cannot blow up the forces
    pos, vol = uniform_box(1, 1, 1, 1 / 3)
    body = Body(BACMaterial(maxdmg=0.3), pos, vol)
    material!(body; horizon=3.015 / 3, rho=7850.0, E=210e9, nu=0.25, Gc=1e-6)
    solver = VelocityVerlet(steps=1)
    dh = Fixtures.handler(body, solver)
    chunk = dh.chunks[1]
    (; storage, system) = chunk
    # the center point keeps its family; a corner point is pulled so far that many of the
    # bonds of its neighbors break, which damages them beyond `maxdmg`
    storage.position[:, 1] .+= 2.0
    Fixtures.force_density!(chunk, solver.Δt, solver.Δt)
    Peridynamics.calc_damage!(chunk)
    @test storage.damage[1] ≈ 1.0
    @test storage.n_active_bonds[1] == 0
    @test all(iszero, storage.b_int[:, 1]) # no bonds, no force on the pulled point
end
