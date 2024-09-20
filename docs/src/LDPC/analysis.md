# LDPC Ensemble Analysis

## Degree Distributions

```@docs
check_concentrated_degree_distribution
```

```@docs
optimal_lambda
```

```@docs
optimal_rho
```

```@docs
optimal_lambda_and_rho
```

```@docs
optimal_threshold
```

## Density Evolution
```@docs
density_evolution
```


# Misc

```@docs
multiplicative_gap_lower_bound
```

```@docs
density_lower_bound
```

See also: `plot_EXIT_chart`, `multiplicative_gap`

## Comments

(these are just temporary notes)

We have functions for the following:

* `optimal_lambda` (`optimal_rho`): Given $\lambda$ (or $\rho$) and a threshold $\epsilon^{BP}$ (or target rate), find the distribution $\rho$ (or $\lambda$) with at least that threshold maximizing design rate (with at least that target rate maximizing threshold).
* `optimal_lambda_and_rho`: Given a target rate (threshold), find distributions $\lambda$ and $\rho$ that maximize threshold (rate).
* `optimal_threshold`: Given distributions $\lambda$ and $\rho$, compute the threshold

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

This shows the accuracy of these functions and how to tune that accuracy. `optimal_lambda` and `optimal_rho` also have an optional keyword parameter `Δ` for tuning accuracy. Note that using `BigFloat` only behaves properly for `optimal_threshold`, any other `Δ` parameter should be `Float64`. Even for `optimal_threshold`, it's best to just use `Float64` unless you are specifically testing numerical stability.
