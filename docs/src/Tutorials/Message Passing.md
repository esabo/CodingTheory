# [Message-passing Decoding](@id message-passing-tutorial)

`CodingTheory.jl` provides reusable workspaces for binary soft-decision
decoding. A workspace stores the Tanner graph, messages, hard decisions, and
schedule, so simulations can decode many channel realizations without
rebuilding the graph.

Every example below continues the same Julia session and uses this
parity-check matrix of the ``[7,4,3]`` Hamming code:

```julia
using CodingTheory

H = UInt8[
    1 0 1 0 1 0 1
    0 1 1 0 0 1 1
    0 0 0 1 1 1 1
]
```

## Codeword decoding

The decoder consumes log-likelihood ratios (LLRs). The convention is

```math
L_i = \log\frac{\Pr(x_i=0\mid y_i)}{\Pr(x_i=1\mid y_i)},
```

so a positive LLR favors zero.

```julia
workspace = init_soft_workspace(H)
received = UInt8[1, 1, 0, 1, 0, 0, 1]
llrs = [bit == 1 ? -5.0 : 5.0 for bit in received]
decoded = zeros(UInt8, size(H, 2))

converged, iterations = decode!(
    workspace,
    llrs;
    algorithm=:sum_product,
    out=decoded,
)
```

`decoded` contains the hard decision. The same array is also available as
`workspace.current_bits` until the next call.

## Syndrome decoding

For error correction it is often more convenient to solve ``He=s`` directly.
Pass the measured syndrome while keeping the channel LLRs as priors on the
error. Loading a channel realization resets the messages, so the workspace
built above can be reused for an unrelated syndrome:

```julia
workspace = init_soft_workspace(H)

p = 0.05
error_llrs = fill(log((1 - p) / p), size(H, 2))
syndrome = UInt8[1, 0, 1]
correction = zeros(UInt8, size(H, 2))

converged, iterations = decode!(
    workspace,
    error_llrs;
    syndrome=syndrome,
    algorithm=:sum_product,
    max_iter=50,
    out=correction,
)
```

Success means that the returned correction has the requested syndrome. It does
not imply recovery of the exact sampled error: any representative of the same
syndrome class is valid for the classical decoding problem.

## Algorithms and schedules

The `algorithm` keyword supports:

- `:sum_product`
- `:min_sum`
- `:normalized_min_sum` with `attenuation`
- `:offset_min_sum` with `offset`, the default
- `:min_sum_correction`

The default flooding schedule updates all checks from the previous iteration.
Layered and serial schedules must be prepared when the workspace is created:

```julia
layered = init_soft_workspace(H; schedule=:layered)
decode!(
    layered,
    error_llrs;
    syndrome=syndrome,
    algorithm=:sum_product,
    schedule=:layered,
)
```

The examples here decode a tiny, dense, high-rate matrix from a constant
channel LLR vector. On such inputs the min-sum family under `:layered` or
`:serial` can stall at an exact stationary point and report no convergence,
which is why these examples request `:sum_product`. The approximations are
intended for large sparse codes with informative per-bit LLRs.

Use `layered_schedule(H)` to inspect the partition. For a large sparse matrix,
prefer the compressed-sparse-row constructor
`init_soft_workspace(row_ptr, col_ind, m, n; base=...)` to avoid a dense scan.
Its `base=0` form accepts SciPy CSR indices directly.

## Erasures and decimation

The `erasures` keyword sets selected channel LLRs to zero, which marks a
position as carrying no information rather than as likely correct. Manual
decimation pins selected variables to known values:

```julia
decode!(
    workspace,
    error_llrs;
    syndrome=syndrome,
    algorithm=:sum_product,
    erasures=[5],
    decimation=:manual,
    decimated_bits=[1],
    decimated_values=[0],
)
```

Erasures make some error patterns free, so they can remove the unique most
likely solution: with a constant prior on this matrix, erasing position 2
makes the weight-two pattern on positions 2 and 7 exactly as likely as the
weight-one pattern on position 5, and the decoder reports no convergence
rather than choosing arbitrarily.

Automatic and guided decimation are selected with `decimation=:auto` and
`:guided`. Set `oscillation=:active` to detect two-cycles in hard decisions; a
negative iteration count reports early termination caused by oscillation.

## Reusing workspaces

Allocate one workspace and output buffer per thread. Repeated `decode!` calls
reset messages and channel state in place:

```julia
workspace = init_soft_workspace(H; schedule=:layered)
output = zeros(UInt8, size(H, 2))

for syndrome in (UInt8[1, 0, 1], UInt8[0, 1, 1])
    converged, iterations = decode!(
        workspace,
        error_llrs;
        syndrome=syndrome,
        algorithm=:sum_product,
        schedule=:layered,
        out=output,
    )
    # consume `output` before the next iteration
end
```

Do not share a mutable workspace between threads.
