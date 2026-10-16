#=
Compare two time-to-first-step measurements and write a markdown table of the
change in each phase. Throws an error if the total time grows by more than
`TTFS_REGRESSION_THRESHOLD`, a fraction that defaults to 0.1.

```bash
julia .buildkite/time_to_first_step/compare.jl baseline.csv changed.csv comparison.md
```
=#

read_timings(path) = map(readlines(path)[2:end]) do line
    phase, seconds = split(line, ',')
    String(phase) => parse(Float64, seconds)
end

baseline_path, changed_path, output_path = ARGS
baseline = Dict(read_timings(baseline_path))
changed = read_timings(changed_path)
const DEFAULT_THRESHOLD = 0.1
threshold = parse(Float64, get(ENV, "TTFS_REGRESSION_THRESHOLD", string(DEFAULT_THRESHOLD)))

relative_change(phase, seconds) = seconds / baseline[phase] - 1

open(output_path, "w") do io
    upstream_packages = get(ENV, "UPSTREAM_PACKAGES", "")
    println(io, "### Time to first step: `$upstream_packages` vs. upstream main\n")
    println(io, "| Phase | Upstream main (s) | Changed (s) | Change |")
    println(io, "|---|---|---|---|")
    for (phase, seconds) in changed
        change = round(100 * relative_change(phase, seconds), digits = 1)
        sign = change >= 0 ? "+" : ""
        b, c = round(baseline[phase], digits = 1), round(seconds, digits = 1)
        println(io, "| $phase | $b | $c | $sign$change% |")
    end
end

total_change = relative_change("Total", Dict(changed)["Total"])
total_change > threshold && error(
    "Time to first step grew by $(round(100 * total_change, digits = 1))%, " *
    "more than the $(100 * threshold)% threshold",
)
