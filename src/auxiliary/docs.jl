# Peridynamics.jl has three API tiers, and every docstring says which one it belongs to.
# See `src/public_api.jl` for the declarations and `docs/src/api_stability.md` for what
# each tier promises.
#
#   1. the user API:      exported, no marker needed
#   2. the extension API: `public`, marked with `extension_api_note()`
#   3. everything else:   internal, marked with `internal_api_warning()`

function internal_api_warning()
    msg = """
    !!! warning "Internal use only"
        This is internal to Peridynamics.jl. It is not part of any API tier and can be
        changed or removed in any release without that being a breaking change. If you need
        it in an extension of your own, please open an issue.
    """
    return msg
end

function extension_api_note()
    msg = """
    !!! note "Extension API"
        This is part of the extension API of Peridynamics.jl. It is not exported, so write
        it as `Peridynamics.<name>` or import it explicitly. It is stable within a minor
        release. While the package is at version 0.x it can still change in a minor version
        bump, and every such change is listed in `NEWS.md`.
    """
    return msg
end

function experimental_api_warning()
    msg = """
    !!! danger "Experimental feature"
        This is an experimental feature. It is not part of the public API of
        Peridynamics.jl and can be changed or removed at any time without that being a
        breaking change. It may also be incomplete or contain bugs, so please use it with
        care.
    """
    return msg
end
