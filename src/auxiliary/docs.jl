# Peridynamics.jl has three API tiers, and every docstring says which one it belongs to:
#
#   1. the user API      -- exported, no marker needed
#   2. the extension API -- not exported, written as `Peridynamics.<name>`, marked with
#                           `extension_api_note()`
#   3. everything else   -- internal, marked with `internal_api_warning()`

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

function extension_api_note()
    msg = """
    !!! note "Extension API"
        This is part of the extension API of Peridynamics.jl. It is not exported, so use it
        as `Peridynamics.<name>` or import it explicitly. It is stable within a minor
        release, but while the package is at version 0.x it can still change in a minor
        version bump. Any such change is listed in `NEWS.md`.
    """
    return msg
end

function experimental_api_warning()
    msg = """
    !!! danger "Experimental feature"
        Please note that this is an experimental feature. It is *not*
        part of the public API of Peridynamics.jl, and thus can be altered (or removed)
        at any time without it being considered a breaking change. Also, the feature may be
        incomplete and/or contain bugs. Please use with caution.
    """
    return msg
end
