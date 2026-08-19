@testitem "required_fields_fracture" begin
    @test Peridynamics.required_fields_fracture(Peridynamics.AbstractMaterial) == ()

    rff_bond_system = (:damage, :n_active_bonds, :bond_active)
    @test Peridynamics.required_fields_fracture(BBMaterial) === rff_bond_system

    rff_interaction_system = (:damage, :n_active_one_nis, :one_ni_active)
    @test Peridynamics.required_fields_fracture(CKIMaterial) === rff_interaction_system
end

@testitem "get_frac_params: damage models without parameters" begin
    struct ParameterlessDamage <: Peridynamics.AbstractDamageModel end
    p = Dict{Symbol,Any}(:Gc => 1.0)
    @test Peridynamics.get_frac_params(ParameterlessDamage(), p, 1.0, 1.0) == (;)
    frac = Peridynamics.get_frac_params(CriticalStretch(), p, 1.0, 1.0)
    @test frac.Gc == 1.0 && frac.εc > 0
end
