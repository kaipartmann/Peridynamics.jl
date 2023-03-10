using Peridynamics, Test

c1 = Contact((1,2), 0.1)
@test c1.body_id_set == (1,2)
@test c1.search_radius == 0.1
@test c1.spring_constant == 1e12

c2 = Contact((3,4), 0.2, 2e13)
@test c2.body_id_set == (3,4)
@test c2.search_radius == 0.2
@test c2.spring_constant == 2e13
