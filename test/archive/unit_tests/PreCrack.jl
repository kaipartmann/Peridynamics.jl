using Peridynamics, Test

precrack = PreCrack([1, 2, 3], [4, 5, 6])
io = IOBuffer()
show(io, "text/plain", precrack)
msg = String(take!(io))
@test msg == "PreCrack:\n 3 points in set a\n 3 points in set b" ||
      msg == "Peridynamics.PreCrack:\n 3 points in set a\n 3 points in set b"
