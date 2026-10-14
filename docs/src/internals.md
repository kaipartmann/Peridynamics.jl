# Internals

Everything on this page is **internal**. It is documented because the documentation helps
when reading or contributing to the package, but it is not part of any API tier. It can be
renamed, resignatured or removed in any release without that being a breaking change. See
[API stability](@ref).

If you are extending Peridynamics.jl, work from the [Extension API](@ref) instead. If
something you need is only available here, please open an issue. Promoting a name to the
extension API is not a breaking change and can be done in a patch release.

This page is generated from every docstring in the package that is neither exported nor
declared `public`, so it cannot drift out of sync with the code.

```@meta
CollapsedDocStrings = true
```

```@autodocs
Modules = [Peridynamics]
Public = false
Private = true
```
