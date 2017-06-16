# RenewalProcesses.jl
Numeric computation of probabilities relevant in Renewal Processes using Julia.

Install using `Pkg.clone("url_of_this_repo")`.

## Functions

The exported functions that can be used are:
```julia
renewal_counting(n::Int, t, pdf, first_pdf = pdf) -> totalt, Qn, Pn
```
Recursively calculates the p.d.f. of the `n`-th event, `Qn`, happening at time
`totalt`, as
well as the probability to have `n` events, `Pn`, occuring up to time `totalt`. Use the "help" for more info.
