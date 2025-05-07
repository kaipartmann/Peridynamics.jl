function internal_api_warning()
    if VERSION < v"1.11"
        msg = """
        !!! warning "Internal use only"
            Please note that this is intended for internal use only. It is *not*
            part of the public API of Peridynamics.jl, and thus can be altered (or removed)
            at any time without it being considered a breaking change.
        """
    else
        msg = ""
    end
    return msg
end
