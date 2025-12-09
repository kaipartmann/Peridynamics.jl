@testsnippet StorageFieldSize begin
    function test_setup(mat, solver; n_chunks=2)
        position = zeros(3, 10)
        position[1, :] = 0.0:9.0
        volume = ones(10)
        body = Body(mat, position, volume)
        material!(body, horizon=1.5, rho=1.0, E=1.0, nu=0.25, Gc=1.0) # dummy parameters
        pd = Peridynamics.PointDecomposition(body, n_chunks)
        ps = Peridynamics.get_param_spec(body)
        chunk = Peridynamics.BodyChunk(body, solver, pd, 1, ps) # return only first chunk
        return chunk.storage, chunk.system
    end
    function get_numbers(system::Peridynamics.AbstractBondSystem)
        n_loc_points = Peridynamics.get_n_loc_points(system)
        n_points = Peridynamics.get_n_points(system)
        n_bonds = Peridynamics.get_n_bonds(system)
        n_loc_dof = Peridynamics.get_n_loc_dof(system)
        n_dof = Peridynamics.get_n_dof(system)
        return n_loc_points, n_points, n_bonds, n_loc_dof, n_dof
    end
    function get_numbers(system::Peridynamics.InteractionSystem)
        n_loc_points = Peridynamics.get_n_loc_points(system)
        n_points = Peridynamics.get_n_points(system)
        n_one_nis = Peridynamics.get_n_one_nis(system)
        n_loc_dof = Peridynamics.get_n_loc_dof(system)
        n_dof = Peridynamics.get_n_dof(system)
        return n_loc_points, n_points, n_one_nis, n_loc_dof, n_dof
    end
end

@testitem "BBMaterial storage field sizes" setup=[StorageFieldSize] begin
    mat = BBMaterial()

    # Test with VelocityVerlet
    solver = VelocityVerlet(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_loc_points)
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_length) == n_bonds
    @test length(s.bond_active) == n_bonds
    @test length(s.residual) == 0
    @test size(s.displacement_copy) == (0, 0)
    @test size(s.b_int_copy) == (0, 0)
    @test length(s.temp_force) == 0
    @test length(s.Δu) == 0
    @test length(s.v_temp) == 0
    @test length(s.Jv_temp) == 0

    # Test with DynamicRelaxation
    solver = DynamicRelaxation(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (3, n_loc_points)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_loc_points)
    @test size(s.b_int_old) == (3, n_loc_points)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (3, n_loc_points)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_length) == n_bonds
    @test length(s.bond_active) == n_bonds
    @test length(s.residual) == 0
    @test size(s.displacement_copy) == (0, 0)
    @test size(s.b_int_copy) == (0, 0)
    @test length(s.temp_force) == 0
    @test length(s.Δu) == 0
    @test length(s.v_temp) == 0
    @test length(s.Jv_temp) == 0

    # Test with Peridynamics.NewtonKrylov
    solver = Peridynamics.NewtonKrylov(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (0, 0)
    @test size(s.velocity_half) == (0, 0)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (0, 0)
    @test size(s.b_int) == (3, n_points)
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_length) == n_bonds
    @test length(s.bond_active) == n_bonds
    @test length(s.residual) == n_loc_dof
    @test size(s.displacement_copy) == (3, n_loc_points)
    @test size(s.b_int_copy) == (3, n_points)
    @test length(s.temp_force) == n_loc_dof
    @test length(s.Δu) == n_loc_dof
    @test length(s.v_temp) == n_loc_dof
    @test length(s.Jv_temp) == n_loc_dof
    @test length(s.precond_diag) == n_loc_dof
end

@testitem "DHBBMaterial storage field sizes" setup=[StorageFieldSize] begin
    mat = DHBBMaterial()

    # Test with VelocityVerlet
    solver = VelocityVerlet(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_points)  # @htlfield
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_length) == n_bonds
    @test length(s.bond_active) == n_bonds
    @test length(s.residual) == 0
    @test size(s.displacement_copy) == (0, 0)
    @test size(s.b_int_copy) == (0, 0)
    @test length(s.temp_force) == 0
    @test length(s.Δu) == 0
    @test length(s.v_temp) == 0
    @test length(s.Jv_temp) == 0

    # Test with DynamicRelaxation
    solver = DynamicRelaxation(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (3, n_loc_points)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_points)
    @test size(s.b_int_old) == (3, n_loc_points)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (3, n_loc_points)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_length) == n_bonds
    @test length(s.bond_active) == n_bonds
    @test length(s.residual) == 0
    @test size(s.displacement_copy) == (0, 0)
    @test size(s.b_int_copy) == (0, 0)
    @test length(s.temp_force) == 0
    @test length(s.Δu) == 0
    @test length(s.v_temp) == 0
    @test length(s.Jv_temp) == 0

    # Test with NewtonKrylov
    solver = Peridynamics.NewtonKrylov(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (0, 0)
    @test size(s.velocity_half) == (0, 0)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (0, 0)
    @test size(s.b_int) == (3, n_points)
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_length) == n_bonds
    @test length(s.bond_active) == n_bonds
    @test length(s.residual) == n_loc_dof
    @test size(s.displacement_copy) == (3, n_loc_points)
    @test size(s.b_int_copy) == (3, n_points)
    @test length(s.temp_force) == n_loc_dof
    @test length(s.Δu) == n_loc_dof
    @test length(s.v_temp) == n_loc_dof
    @test length(s.Jv_temp) == n_loc_dof
    @test length(s.precond_diag) == n_loc_dof
end

@testitem "GBBMaterial storage field sizes" setup=[StorageFieldSize] begin
    mat = GBBMaterial()

    # Test with VelocityVerlet
    solver = VelocityVerlet(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_loc_points)
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.weighted_volume) == n_loc_points
    @test length(s.bond_length) == n_bonds
    @test length(s.bond_active) == n_bonds
    @test length(s.residual) == 0
    @test size(s.displacement_copy) == (0, 0)
    @test size(s.b_int_copy) == (0, 0)
    @test length(s.temp_force) == 0
    @test length(s.Δu) == 0
    @test length(s.v_temp) == 0
    @test length(s.Jv_temp) == 0

    # Test with DynamicRelaxation
    solver = DynamicRelaxation(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (3, n_loc_points)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_loc_points)
    @test size(s.b_int_old) == (3, n_loc_points)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (3, n_loc_points)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.weighted_volume) == n_loc_points
    @test length(s.bond_length) == n_bonds
    @test length(s.bond_active) == n_bonds
    @test length(s.residual) == 0
    @test size(s.displacement_copy) == (0, 0)
    @test size(s.b_int_copy) == (0, 0)
    @test length(s.temp_force) == 0
    @test length(s.Δu) == 0
    @test length(s.v_temp) == 0
    @test length(s.Jv_temp) == 0

    # Test with Peridynamics.NewtonKrylov
    solver = Peridynamics.NewtonKrylov(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (0, 0)
    @test size(s.velocity_half) == (0, 0)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (0, 0)
    @test size(s.b_int) == (3, n_points)
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.weighted_volume) == n_loc_points
    @test length(s.bond_length) == n_bonds
    @test length(s.bond_active) == n_bonds
    @test length(s.residual) == n_loc_dof
    @test size(s.displacement_copy) == (3, n_loc_points)
    @test size(s.b_int_copy) == (3, n_points)
    @test length(s.temp_force) == n_loc_dof
    @test length(s.Δu) == n_loc_dof
    @test length(s.v_temp) == n_loc_dof
    @test length(s.Jv_temp) == n_loc_dof
    @test length(s.precond_diag) == n_loc_dof
end

@testitem "OSBMaterial storage field sizes" setup=[StorageFieldSize] begin
    mat = OSBMaterial()

    # Test with VelocityVerlet
    solver = VelocityVerlet(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_points)  # @htlfield
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_length) == n_bonds
    @test length(s.bond_active) == n_bonds
    @test length(s.residual) == 0
    @test size(s.displacement_copy) == (0, 0)
    @test size(s.b_int_copy) == (0, 0)
    @test length(s.temp_force) == 0
    @test length(s.Δu) == 0
    @test length(s.v_temp) == 0
    @test length(s.Jv_temp) == 0

    # Test with DynamicRelaxation
    solver = DynamicRelaxation(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (3, n_loc_points)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_points)
    @test size(s.b_int_old) == (3, n_loc_points)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (3, n_loc_points)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_length) == n_bonds
    @test length(s.bond_active) == n_bonds
    @test length(s.residual) == 0
    @test size(s.displacement_copy) == (0, 0)
    @test size(s.b_int_copy) == (0, 0)
    @test length(s.temp_force) == 0
    @test length(s.Δu) == 0
    @test length(s.v_temp) == 0
    @test length(s.Jv_temp) == 0

    # Test with NewtonKrylov
    solver = Peridynamics.NewtonKrylov(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (0, 0)
    @test size(s.velocity_half) == (0, 0)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (0, 0)
    @test size(s.b_int) == (3, n_points)
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_length) == n_bonds
    @test length(s.bond_active) == n_bonds
    @test length(s.residual) == n_loc_dof
    @test size(s.displacement_copy) == (3, n_loc_points)
    @test size(s.b_int_copy) == (3, n_points)
    @test length(s.temp_force) == n_loc_dof
    @test length(s.Δu) == n_loc_dof
    @test length(s.v_temp) == n_loc_dof
    @test length(s.Jv_temp) == n_loc_dof
    @test length(s.precond_diag) == n_loc_dof
end

@testitem "CMaterial storage field sizes" setup=[StorageFieldSize] begin
    mat = CMaterial()

    # Test with VelocityVerlet
    solver = VelocityVerlet(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_points)  # @htlfield
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test size(s.defgrad) == (9, n_loc_points)
    @test size(s.cauchy_stress) == (9, n_loc_points)
    @test length(s.von_mises_stress) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_active) == n_bonds
    @test length(s.residual) == 0
    @test size(s.displacement_copy) == (0, 0)
    @test size(s.b_int_copy) == (0, 0)
    @test length(s.temp_force) == 0
    @test length(s.Δu) == 0
    @test length(s.v_temp) == 0
    @test length(s.Jv_temp) == 0

    # Test with DynamicRelaxation
    solver = DynamicRelaxation(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (3, n_loc_points)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_points)
    @test size(s.b_int_old) == (3, n_loc_points)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (3, n_loc_points)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test size(s.defgrad) == (9, n_loc_points)
    @test size(s.cauchy_stress) == (9, n_loc_points)
    @test length(s.von_mises_stress) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_active) == n_bonds
    @test length(s.residual) == 0
    @test size(s.displacement_copy) == (0, 0)
    @test size(s.b_int_copy) == (0, 0)
    @test length(s.temp_force) == 0
    @test length(s.Δu) == 0
    @test length(s.v_temp) == 0
    @test length(s.Jv_temp) == 0

    # Test with NewtonKrylov
    solver = Peridynamics.NewtonKrylov(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (0, 0)
    @test size(s.velocity_half) == (0, 0)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (0, 0)
    @test size(s.b_int) == (3, n_points)
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test size(s.defgrad) == (9, n_loc_points)
    @test size(s.cauchy_stress) == (9, n_loc_points)
    @test length(s.von_mises_stress) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_active) == n_bonds
    @test length(s.residual) == n_loc_dof
    @test size(s.displacement_copy) == (3, n_loc_points)
    @test size(s.b_int_copy) == (3, n_points)
    @test length(s.temp_force) == n_loc_dof
    @test length(s.Δu) == n_loc_dof
    @test length(s.v_temp) == n_loc_dof
    @test length(s.Jv_temp) == n_loc_dof
    @test length(s.precond_diag) == n_loc_dof
end

@testitem "CRMaterial storage field sizes" setup=[StorageFieldSize] begin
    mat = CRMaterial()

    # Test with VelocityVerlet
    solver = VelocityVerlet(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.velocity_half) == (3, n_points)  # @lthfield
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_points)  # @htlfield
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test size(s.unrotated_stress) == (9, n_loc_points)
    @test size(s.defgrad) == (9, n_loc_points)
    @test size(s.cauchy_stress) == (9, n_loc_points)
    @test length(s.von_mises_stress) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test size(s.left_stretch) == (9, n_loc_points)
    @test size(s.rotation) == (9, n_loc_points)
    @test length(s.bond_active) == n_bonds

    # Test with DynamicRelaxation
    solver = DynamicRelaxation(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.velocity_half) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (3, n_loc_points)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_points)
    @test size(s.b_int_old) == (3, n_loc_points)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (3, n_loc_points)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test size(s.unrotated_stress) == (9, n_loc_points)
    @test size(s.defgrad) == (9, n_loc_points)
    @test size(s.cauchy_stress) == (9, n_loc_points)
    @test length(s.von_mises_stress) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test size(s.left_stretch) == (9, n_loc_points)
    @test size(s.rotation) == (9, n_loc_points)
    @test length(s.bond_active) == n_bonds
end

@testitem "RKCMaterial storage field sizes" setup=[StorageFieldSize] begin
    mat = RKCMaterial()

    # Test with VelocityVerlet
    solver = VelocityVerlet(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.defgrad) == (9, n_points)  # @lthfield
    @test length(s.weighted_volume) == n_points  # @lthfield
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_points)  # @htlfield
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.update_gradients) == n_loc_points
    @test size(s.cauchy_stress) == (9, n_loc_points)
    @test length(s.von_mises_stress) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_active) == n_bonds
    @test size(s.gradient_weight) == (3, n_bonds)
    @test size(s.bond_first_piola_kirchhoff) == (9, n_bonds)
    @test length(s.residual) == 0
    @test size(s.displacement_copy) == (0, 0)
    @test size(s.b_int_copy) == (0, 0)
    @test length(s.temp_force) == 0
    @test length(s.Δu) == 0
    @test length(s.v_temp) == 0
    @test length(s.Jv_temp) == 0

    # Test with DynamicRelaxation
    solver = DynamicRelaxation(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.defgrad) == (9, n_points)
    @test length(s.weighted_volume) == n_points
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (3, n_loc_points)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_points)
    @test size(s.b_int_old) == (3, n_loc_points)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (3, n_loc_points)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.update_gradients) == n_loc_points
    @test size(s.cauchy_stress) == (9, n_loc_points)
    @test length(s.von_mises_stress) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_active) == n_bonds
    @test size(s.gradient_weight) == (3, n_bonds)
    @test size(s.bond_first_piola_kirchhoff) == (9, n_bonds)
    @test length(s.residual) == 0
    @test size(s.displacement_copy) == (0, 0)
    @test size(s.b_int_copy) == (0, 0)
    @test length(s.temp_force) == 0
    @test length(s.Δu) == 0
    @test length(s.v_temp) == 0
    @test length(s.Jv_temp) == 0

    # Test with NewtonKrylov, but only with 1 chunk!
    solver = Peridynamics.NewtonKrylov(steps=1)
    s, system = test_setup(mat, solver; n_chunks=1)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.defgrad) == (9, n_points)
    @test length(s.weighted_volume) == n_points
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (0, 0)
    @test size(s.velocity_half) == (0, 0)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (0, 0)
    @test size(s.b_int) == (3, n_points)
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.update_gradients) == n_loc_points
    @test size(s.cauchy_stress) == (9, n_loc_points)
    @test length(s.von_mises_stress) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_active) == n_bonds
    @test size(s.gradient_weight) == (3, n_bonds)
    @test size(s.bond_first_piola_kirchhoff) == (9, n_bonds)
    @test length(s.residual) == n_loc_dof
    @test size(s.displacement_copy) == (3, n_loc_points)
    @test size(s.b_int_copy) == (3, n_points)
    @test length(s.temp_force) == n_loc_dof
    @test length(s.Δu) == n_loc_dof
    @test length(s.v_temp) == n_loc_dof
    @test length(s.Jv_temp) == n_loc_dof
    @test length(s.precond_diag) == n_loc_dof
end

@testitem "RKCRMaterial storage field sizes" setup=[StorageFieldSize] begin
    mat = RKCRMaterial()

    # Test with VelocityVerlet
    solver = VelocityVerlet(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.velocity_half) == (3, n_points)  # @lthfield
    @test size(s.defgrad) == (9, n_points)  # @lthfield
    @test size(s.defgrad_dot) == (9, n_points)  # @lthfield
    @test length(s.weighted_volume) == n_points  # @lthfield
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_points)  # @htlfield
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.update_gradients) == n_loc_points
    @test size(s.cauchy_stress) == (9, n_loc_points)
    @test length(s.von_mises_stress) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_active) == n_bonds
    @test size(s.gradient_weight) == (3, n_bonds)
    @test size(s.bond_first_piola_kirchhoff) == (9, n_bonds)
    @test size(s.left_stretch) == (9, n_bonds)
    @test size(s.rotation) == (9, n_bonds)
    @test size(s.bond_unrot_cauchy_stress) == (9, n_bonds)

    # Test with DynamicRelaxation
    solver = DynamicRelaxation(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.velocity_half) == (3, n_points)
    @test size(s.defgrad) == (9, n_points)
    @test size(s.defgrad_dot) == (9, n_points)
    @test length(s.weighted_volume) == n_points
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (3, n_loc_points)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_points)
    @test size(s.b_int_old) == (3, n_loc_points)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (3, n_loc_points)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test length(s.update_gradients) == n_loc_points
    @test size(s.cauchy_stress) == (9, n_loc_points)
    @test length(s.von_mises_stress) == n_loc_points
    @test length(s.strain_energy_density) == n_loc_points
    @test length(s.bond_active) == n_bonds
    @test size(s.gradient_weight) == (3, n_bonds)
    @test size(s.bond_first_piola_kirchhoff) == (9, n_bonds)
    @test size(s.left_stretch) == (9, n_bonds)
    @test size(s.rotation) == (9, n_bonds)
    @test size(s.bond_unrot_cauchy_stress) == (9, n_bonds)
end

@testitem "BACMaterial storage field sizes" setup=[StorageFieldSize] begin
    mat = BACMaterial()

    # Test with VelocityVerlet
    solver = VelocityVerlet(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_points)  # @htlfield
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test size(s.stress) == (9, n_loc_points)
    @test length(s.von_mises_stress) == n_loc_points
    @test length(s.bond_active) == n_bonds

    # Test with DynamicRelaxation
    solver = DynamicRelaxation(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (3, n_loc_points)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_points)
    @test size(s.b_int_old) == (3, n_loc_points)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (3, n_loc_points)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test size(s.stress) == (9, n_loc_points)
    @test length(s.von_mises_stress) == n_loc_points
    @test length(s.bond_active) == n_bonds

    # Test with NewtonKrylov
    solver = Peridynamics.NewtonKrylov(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (0, 0)
    @test size(s.velocity_half) == (0, 0)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (0, 0)
    @test size(s.b_int) == (3, n_points)
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_bonds) == n_loc_points
    @test size(s.stress) == (9, n_loc_points)
    @test length(s.von_mises_stress) == n_loc_points
    @test length(s.bond_active) == n_bonds
    @test length(s.residual) == n_loc_dof
    @test size(s.displacement_copy) == (3, n_loc_points)
    @test size(s.b_int_copy) == (3, n_points)
    @test length(s.temp_force) == n_loc_dof
    @test length(s.Δu) == n_loc_dof
    @test length(s.v_temp) == n_loc_dof
    @test length(s.Jv_temp) == n_loc_dof
    @test length(s.precond_diag) == n_loc_dof
end

@testitem "CKIMaterial storage field sizes" setup=[StorageFieldSize] begin
    mat = CKIMaterial()

    # Test with VelocityVerlet
    solver = VelocityVerlet(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_loc_points)
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_one_nis) == n_loc_points
    @test length(s.residual) == 0
    @test size(s.displacement_copy) == (0, 0)
    @test size(s.b_int_copy) == (0, 0)
    @test length(s.temp_force) == 0
    @test length(s.Δu) == 0
    @test length(s.v_temp) == 0
    @test length(s.Jv_temp) == 0

    # Test with DynamicRelaxation
    solver = DynamicRelaxation(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (3, n_loc_points)
    @test size(s.velocity_half) == (3, n_loc_points)
    @test size(s.velocity_half_old) == (3, n_loc_points)
    @test size(s.acceleration) == (3, n_loc_points)
    @test size(s.b_int) == (3, n_loc_points)
    @test size(s.b_int_old) == (3, n_loc_points)
    @test size(s.b_ext) == (3, n_loc_points)
    @test size(s.density_matrix) == (3, n_loc_points)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_one_nis) == n_loc_points
    @test length(s.residual) == 0
    @test size(s.displacement_copy) == (0, 0)
    @test size(s.b_int_copy) == (0, 0)
    @test length(s.temp_force) == 0
    @test length(s.Δu) == 0
    @test length(s.v_temp) == 0
    @test length(s.Jv_temp) == 0

    # Test with NewtonKrylov
    solver = Peridynamics.NewtonKrylov(steps=1)
    s, system = test_setup(mat, solver)
    n_loc_points, n_points, n_bonds, n_loc_dof, n_dof = get_numbers(system)
    @test size(s.position) == (3, n_points)
    @test size(s.displacement) == (3, n_loc_points)
    @test size(s.velocity) == (0, 0)
    @test size(s.velocity_half) == (0, 0)
    @test size(s.velocity_half_old) == (0, 0)
    @test size(s.acceleration) == (0, 0)
    @test size(s.b_int) == (3, n_points)
    @test size(s.b_int_old) == (0, 0)
    @test size(s.b_ext) == (3, n_points)
    @test size(s.density_matrix) == (0, 0)
    @test length(s.damage) == n_loc_points
    @test length(s.n_active_one_nis) == n_loc_points
    @test length(s.residual) == n_loc_dof
    @test size(s.displacement_copy) == (3, n_loc_points)
    @test size(s.b_int_copy) == (3, n_points)
    @test length(s.temp_force) == n_loc_dof
    @test length(s.Δu) == n_loc_dof
    @test length(s.v_temp) == n_loc_dof
    @test length(s.Jv_temp) == n_loc_dof
    @test length(s.precond_diag) == n_loc_dof
end
