using Peridynamics, Test

es1 = ExportSettings()
@test es1.path == ""
@test es1.exportfreq == 0
@test es1.resultfile_prefix == ""
@test es1.logfile == ""
@test es1.exportflag == false

io = IOBuffer()
show(io, "text/plain", es1)
msg = String(take!(io))
@test msg == "ExportSettings" || msg == "Peridynamics.ExportSettings"

es2 = ExportSettings("test/path", 10)
@test es2.path == "test/path"
@test es2.exportfreq == 10
@test es2.resultfile_prefix == ""
@test es2.logfile == ""
@test es2.exportflag == true

io = IOBuffer()
show(io, "text/plain", es2)
msg = String(take!(io))
@test msg == "ExportSettings" || msg == "Peridynamics.ExportSettings"

