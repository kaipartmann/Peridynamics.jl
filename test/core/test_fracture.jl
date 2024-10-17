@testitem "required_fields_fracture" begin
    @test Peridynamics.required_fields_fracture(Peridynamics.AbstractMaterial) == ()

    rff_bond_system = (:damage, :n_active_bonds, :bond_active)
    @test Peridynamics.required_fields_fracture(BBMaterial) === rff_bond_system

    rff_interaction_system = (:damage, :n_active_one_nis, :one_ni_active)
    @test Peridynamics.required_fields_fracture(CKIMaterial) === rff_interaction_system
end
