using TestItemRunner

# Only run tests that are not tagged with :skipci due to performance issues in CI
@run_package_tests verbose=true filter=ti->!(:skipci in ti.tags)
