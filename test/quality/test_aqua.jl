@testitem "Aqua.jl" tags=[:lint] begin
    # Aqua introspects methods and project files and contributes nothing to the line coverage of
    # `src`, and `persistent_tasks` precompiles the package in a fresh process. Both make it an
    # `:lint` item that runs in the `extras` CI job, off the critical path.
    using Aqua
    Aqua.test_all(Peridynamics; stale_deps=(ignore=[:ThreadPinning],))
end
