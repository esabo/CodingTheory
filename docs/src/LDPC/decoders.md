# [Decoding LDPC Codes](@id ldpc-decoders-api)

## Message passing

The soft-decision decoder keeps its Tanner graph, messages, hard decisions, and
schedule in a reusable workspace, so a simulation decodes many channel
realizations without rebuilding the graph. Sum-product and min-sum variants are
selected through the `algorithm` keyword to `decode!`, and the schedule must be
prepared when the workspace is created. Construct one workspace and one output
buffer per thread; see
[Message-passing Decoding](@ref message-passing-tutorial) for a complete
workflow.

```@autodocs
Modules = [CodingTheory]
Pages = ["LDPC/MP_decoders.jl"]
Private = false
```

## Linear programming

```@autodocs
Modules = [CodingTheory]
Pages = ["LDPC/LP_decoders.jl"]
Private = false
```

## Post-processing decoders

Message passing can terminate without satisfying the parity checks, most often
on a short cycle or a trapping set. The decoders below take the soft output of
that failed attempt, the total log-likelihood ratios, and search nearby for a
vector with the requested syndrome. They are used after `decode!` rather than
in place of it, and each one has its own workspace so the allocation happens
once.

*Ordered statistics decoding* (OSD) sorts the positions by reliability ``|L_v|``
and eliminates the parity-check matrix from the least reliable end, which leaves
an information set drawn from the most reliable positions. It then flips small
patterns inside that information set, re-encodes the remaining positions from
the syndrome, and keeps the candidate of least soft cost, the sum of ``|L_v|``
over the positions that disagree with the original hard decisions. The `order`
keyword bounds the weight of those patterns and so controls the cost. With the
default `method = :cs`, the weight-two sweep is restricted to the `cs_lambda`
least reliable positions of the information set instead of all of it.

*Guessing random additive noise decoding* (GRAND) enumerates low-weight error
patterns over the least reliable positions and returns the cheapest one whose
syndrome matches. Its search is bounded by the cost of the best match so far, so
the result is the minimum-cost pattern in the search space rather than the first
one found.

*Weighted bit flipping* (WBF) is the cheapest of the three. It repeatedly flips
the single position maximizing

```math
\mathrm{score}(v) = |\{c \sim v : c \text{ unsatisfied}\}| - \alpha\,|L_v|,
```

that is, the position explaining the most unsatisfied checks after discounting
the decoder's confidence ``|L_v|`` in that position. The weight ``\alpha`` is
the `alpha` keyword. Flipping stops when the residual syndrome vanishes or after
`max_iters` flips, so WBF is a local search and, unlike OSD and GRAND, offers no
optimality guarantee over its candidates.

```@autodocs
Modules = [CodingTheory]
Pages = ["LDPC/decoder_post.jl"]
Private = false
```

Generalized belief propagation, which passes messages between regions of the
Tanner graph instead of individual nodes, has its own page; see
[Generalized belief propagation](@ref ldpc-gbp-api).
