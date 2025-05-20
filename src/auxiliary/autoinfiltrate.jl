"""
    @autoinfiltrate
    @autoinfiltrate condition::Bool

$(internal_api_warning())

Invoke the `@infiltrate` macro of the package Infiltrator.jl to create a breakpoint for
ad-hoc interactive debugging in the REPL. If the optional argument `condition` is given, the
breakpoint is only enabled if `condition` evaluates to `true`.

As opposed to using `Infiltrator.@infiltrate` directly, this macro does not require
Infiltrator.jl to be added as a dependency to Peridynamics.jl. As a bonus, the macro will
also attempt to load the Infiltrator module if it has not yet been loaded manually.

Note: For this macro to work, the Infiltrator.jl package needs to be installed in your
current Julia environment stack.

See also: [Infiltrator.jl](https://github.com/JuliaDebug/Infiltrator.jl)
"""
macro autoinfiltrate(condition=true)
    pkgid = Base.PkgId(Base.UUID("5903a43b-9cc3-4c30-8d17-598619ec4e9b"), "Infiltrator")
    error_msg = "Cannot load Infiltrator.jl!"
    error_msg *= "Make sure it is included in your environment stack."
    if !haskey(Base.loaded_modules, pkgid)
        try
            Base.eval(Main, :(using Infiltrator))
        catch err
            @error error_msg
            Base.showerror(stderr, err)
        end
    end
    i = get(Base.loaded_modules, pkgid, nothing)
    lnn = LineNumberNode(__source__.line, __source__.file)
    isnothing(i) && return Expr(:macrocall, Symbol("@warn"), lnn, error_msg)
    return Expr(:macrocall, Expr(:., i, QuoteNode(Symbol("@infiltrate"))), lnn,
                esc(condition))
end
