# Terminal formatting shared by `compare.jl` and `hpc/report.jl`.
#
# `printstyled` emits escape codes only when its output is a terminal, so piping a report into a
# file gives a clean file and nothing here needs a `--no-color` switch. The verdict is carried
# by a plain character (`▲ ▼ ·`) as well as by the color, so the tables stay readable in a log
# file and for a reader who cannot distinguish red from green.

using Printf

const WIDTH = 76

hline() = printstyled("─"^WIDTH, "\n"; color=:light_black)

function title(headline, subtitle)
    hline()
    printstyled(" ", headline, "\n"; bold=true)
    printstyled(" ", subtitle, "\n"; color=:light_black)
    hline()
    return nothing
end

section(name) = printstyled("\n ", name, "\n"; bold=true, color=:cyan)
colheads(cells...) = printstyled(cells..., "\n"; color=:light_black)

function verdict_style(verdict::Symbol)
    verdict === :regression && return (color=:red, bold=true, marker='▲')
    verdict === :improvement && return (color=:green, bold=true, marker='▼')
    return (color=:light_black, bold=false, marker='·')
end

# A relative change as a signed percentage, e.g. `+21.6 %`
pct(ratio::Real) = @sprintf("%+.1f %%", 100 * (ratio - 1))

# An integer with thousands separators, e.g. `2,700,000`
function prettycount(n::Integer)
    digits = reverse(string(n))
    return reverse(join((digits[i:min(i + 2, end)] for i in 1:3:length(digits)), ","))
end
