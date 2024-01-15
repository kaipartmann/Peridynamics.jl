is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")
