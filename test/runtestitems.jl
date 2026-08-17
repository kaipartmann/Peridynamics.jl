# julia -t 6 test/runtestitems.jl [name-fragment | file:<path> | tag:<tag> ...]
#
# Plain arguments match test item names, not filenames, and every selector has to match:
#
#     julia -t 6 test/runtestitems.jl
#     julia -t 6 test/runtestitems.jl BBMaterial dynamic
#     julia -t 6 test/runtestitems.jl file:test/core/test_halo_exchange.jl
#     julia -t 6 test/runtestitems.jl tag:mpi
#
# Needs TestEnv: julia -e 'import Pkg; Pkg.add("TestEnv")'

import Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

if Base.find_package("TestEnv") === nothing
	error("TestEnv.jl is not installed in Julia's global environment. Install it with:\n\n" *
		  "    julia -e 'import Pkg; Pkg.add(\"TestEnv\")'")
end

using TestEnv

TestEnv.activate()

using TestItemRunner

name_filters = String[]
file_filters = String[]
tag_filters = Symbol[]
for argument in ARGS
	if startswith(lowercase(argument), "tag:")
		tag = argument[5:end]
		isempty(tag) && error("A tag selector must use the form tag:<tag>.")
		push!(tag_filters, Symbol(lowercase(tag)))
	elseif startswith(lowercase(argument), "file:")
		file = argument[6:end]
		isempty(file) && error("A file selector must use the form file:<path>.")
		push!(file_filters, normpath(abspath(file)))
	else
		push!(name_filters, argument)
	end
end

@run_package_tests verbose=true filter=ti ->
	all(contains(lowercase(ti.name), lowercase(name_filter)) for name_filter in name_filters) &&
	all(normpath(abspath(ti.filename)) == file for file in file_filters) &&
	all(tag in ti.tags for tag in tag_filters)
