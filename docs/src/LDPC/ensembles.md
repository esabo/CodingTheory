# LDPC Ensemble Analysis

## Degree Distributions

```@docs
checkconcentrateddegreedistribution
```

```@docs
optimallambda
```

```@docs
optimalrho
```

```@docs
optimallambdaandrho
```

```@docs
optimalthreshold
```

## Density Evolution
```@docs
densityevolution
```

```@docs
plotEXITchart
```

# Misc
```@docs
multiplicativegap
```

```@docs
multiplicativegaplowerbound
```

```@docs
densitylowerbound
```

## Comments

(these are just temporary notes)

We have functions for the following:

* `optimallambda` (`optimalrho`): Given $\lambda$ (or $\rho$) and a threshold $\epsilon^{BP}$ (or target rate), find the distribution $\rho$ (or $\lambda$) with at least that threshold maximizing design rate (with at least that target rate maximizing threshold).
* `optimallambdaandrho`: Given a target rate (threshold), find distributions $\lambda$ and $\rho$ that maximize threshold (rate).
* `optimalthreshold`: Given distributions $\lambda$ and $\rho$, compute the threshold

Example of using `optimallambdaandrho` and `optimalthreshold`:

```julia-repl
julia> λ, ρ, r, ε = optimallambdaandrho(8, 6, 0.4, :ε); 0.4 - optimalthreshold(λ, ρ)
2.849104958069226e-7

julia> λ, ρ, r, ε = optimallambdaandrho(8, 6, 0.4, :ε, Δλ = 0.0001); 0.4 - optimalthreshold(λ, ρ)
1.0256726462598564e-7

julia> λ, ρ, r, ε = optimallambdaandrho(8, 6, 0.4, :ε, Δλ = 0.0001); 0.4 - optimalthreshold(λ, ρ, Δ = BigFloat("1e-7"))
1.025672727266482436145720743991009459178786186514042575060182100591178904349331e-07

julia> λ, ρ, r, ε = optimallambdaandrho(8, 6, 0.4, :ε, Δλ = 0.0001, Δρ = 0.001); 0.4 - optimalthreshold(λ, ρ, Δ = BigFloat("1e-7"))
1.025672727266482436145720743991009459178786186514042575060182100591178904349331e-07
```

This shows the accuracy of these functions and how to tune that accuracy. `optimallambda` and `optimalrho` also have an optional keyword parameter `Δ` for tuning accuracy. Note that using `BigFloat` only behaves properly for `optimalthreshold`, any other `Δ` parameter should be `Float64`. Even for `optimalthreshold`, it's best to just use `Float64` unless you are specifically testing numerical stability.
