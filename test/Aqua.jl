using Peridynamics, Aqua

Aqua.test_all(Peridynamics; stale_deps=(ignore=[:ThreadPinning],))
