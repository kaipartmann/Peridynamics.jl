using Peridynamics, Test

es1 = ExportSettings()
@test es1.path == ""
@test es1.exportfreq == 0
@test es1.resultfile_prefix == ""
@test es1.logfile == ""
@test es1.exportflag == false

es2 = ExportSettings("test/path", 10)
@test es2.path == "test/path"
@test es2.exportfreq == 10
@test es2.resultfile_prefix == ""
@test es2.logfile == ""
@test es2.exportflag == true
