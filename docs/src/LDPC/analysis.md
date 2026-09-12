# [LDPC Ensemble Analysis](@id ldpc-analysis-api)

An LDPC ensemble is specified by its degree distributions rather than by a
particular parity-check matrix, and density evolution tracks the distribution
of decoder messages through the iterations to predict the ensemble's threshold:
the worst channel parameter for which the error probability still converges to
zero as the blocklength grows. The design of good irregular codes is the search
for degree distributions with a high threshold at a given rate, which is what
the optimization routines below do.

The ensemble types themselves, `LDPCEnsemble` and `METEnsemble`, are documented
with the [channels](@ref ldpc-channels-api). EXIT and protograph-EXIT charts visualize the same
convergence question one iteration at a time; the plotting functions require a
Makie backend to be loaded.

## Optimizing degree distributions

The optimization helpers support the following workflows:

* `optimal_lambda` and `optimal_rho`: given ``\rho`` (or ``\lambda``) and a
  threshold ``\epsilon^{BP}`` or target rate, find the other distribution with
  at least that threshold maximizing design rate, or with at least that target
  rate maximizing threshold.
* `optimal_lambda_and_rho`: given a target rate or threshold, find both
  ``\lambda`` and ``\rho`` maximizing the other quantity.
* `optimal_threshold`: given ``\lambda`` and ``\rho``, compute the threshold.

Example of using `optimal_lambda_and_rho` and `optimal_threshold`:

```julia-repl
julia> λ, ρ, r, ε = optimal_lambda_and_rho(8, 6, 0.4, :ε); 0.4 - optimal_threshold(λ, ρ)
2.849104958069226e-7

julia> λ, ρ, r, ε = optimal_lambda_and_rho(8, 6, 0.4, :ε, Δλ = 0.0001); 0.4 - optimal_threshold(λ, ρ)
1.0256726462598564e-7

julia> λ, ρ, r, ε = optimal_lambda_and_rho(8, 6, 0.4, :ε, Δλ = 0.0001); 0.4 - optimal_threshold(λ, ρ, Δ = BigFloat("1e-7"))
1.025672727266482436145720743991009459178786186514042575060182100591178904349331e-07

julia> λ, ρ, r, ε = optimal_lambda_and_rho(8, 6, 0.4, :ε, Δλ = 0.0001, Δρ = 0.001); 0.4 - optimal_threshold(λ, ρ, Δ = BigFloat("1e-7"))
1.025672727266482436145720743991009459178786186514042575060182100591178904349331e-07
```

This shows the accuracy of these functions and how to tune it. `optimal_lambda`
and `optimal_rho` also take a keyword `Δ` for the same purpose. Note that
`BigFloat` only behaves properly for `optimal_threshold`; any other `Δ` should
be a `Float64`, and even for `optimal_threshold` it is best to use `Float64`
unless you are specifically testing numerical stability.

```@autodocs
Modules = [CodingTheory]
Pages = ["LDPC/analysis.jl"]
Private = false
```
