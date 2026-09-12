# Quantum Weight Reduction

Weight reduction rewrites a quantum code so that its stabilizer generators and
qubit degrees meet a target, at the cost of extra qubits and checks. The
classical Hastings procedure is a composition of four steps, each exposed
separately here: copying splits a high-degree qubit into several, gauging
splits a high-weight check, thickening adds layers and chooses heights, and
coning removes the resulting high-weight checks along a cone. `copying` and
`gauging` also have formulations as instances of coning, which are provided for
comparison and testing.

`quantum_weight_reduction` runs the full pipeline. Note that weight reduction
preserves the number of logical qubits but generally reduces the relative
distance, so the reduced code is not a free improvement; see the
[weight reduction tutorial](@ref weight-reduction-tutorial) for the trade-offs.

```@autodocs
Modules = [CodingTheory]
Pages = ["Quantum/weight_reduction.jl"]
Private = false
```
