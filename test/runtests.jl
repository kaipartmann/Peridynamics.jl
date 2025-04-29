using TestItemRunner

# Only run tests that are not tagged with :mpi due to performance issues in CI
@run_package_tests verbose=true filter=ti->!(:mpi in ti.tags)
