@testitem "required_fields" begin
    rf_abstractmat = (:position, :displacement, :velocity, :velocity_half, :acceleration,
                      :b_int, :b_ext, :velocity_half_old, :b_int_old, :density_matrix)
    @test Peridynamics.required_fields(Peridynamics.AbstractMaterial) === rf_abstractmat

    rf_bb = (:position, :displacement, :velocity, :velocity_half, :acceleration, :b_int,
             :b_ext, :velocity_half_old, :b_int_old, :density_matrix, :damage,
             :n_active_bonds, :bond_active)
    @test Peridynamics.required_fields(BBMaterial) === rf_bb

    rf_cki = (:position, :displacement, :velocity, :velocity_half, :acceleration, :b_int,
              :b_ext, :velocity_half_old, :b_int_old, :density_matrix, :damage,
              :n_active_one_nis, :one_ni_active)
    @test Peridynamics.required_fields(CKIMaterial) === rf_cki
end
