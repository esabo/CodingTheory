# Findings on `src/LDPC/MP_decoders.jl` and `src/LDPC/decoder_post.jl`

Branch: `weight_dist`. Baseline commit for everything below: `a4dfb8e`
("Fix OSD syndrome workspace and accept Flint matrices in the new decoders").

This note is written for someone maintaining CodingTheory who has no context on
the downstream project that surfaced these issues. It lives in a new top-level
`notes/` directory rather than under `docs/`, because `docs/` is a Documenter.jl
site driven by `docs/make.jl` and this is engineering notes, not user
documentation. Move it if the repo grows a better home for that.

Every claim is tagged:

* **[run]** — verified by executing code and reading the output.
* **[src]** — read directly from the source, not executed.
* **[hyp]** — hypothesis, stated as such.

Reproduction harness used for the **[run]** results in sections 2 and 4: a
throwaway module that stubs `fpMatrix`, `FqMatrix` and
`_Flint_matrix_to_Julia_support_matrix` and then `include`s
`src/LDPC/MP_decoders.jl` (and, where noted, `src/LDPC/decoder_post.jl`)
directly. That loads the two files with no Oscar dependency and makes the
internal functions (`_check_update!`, `_harden!`, `_syndrome_matches`) callable,
which is what the instrumented traces below need. Nothing about the results
depends on the harness; the section 1 and 5 results come from the real package.

---

## 1. Three bugs fixed in `a4dfb8e`

Read from `git show a4dfb8e`.

### 1.1 `OSDWorkspace` used an undeclared field `s_work`

`osd_decode!` and `_fast_osd!` in `src/LDPC/decoder_post.jl` both read and wrote
`W.s_work` — the syndrome that has to be row-reduced in lockstep with `H_work`
during the Gaussian elimination, and that `_evaluate_pattern!` reads at
`decoder_post.jl:86` to re-encode the least-reliable positions. The field was
never declared on the struct and never allocated.

Fix: added `s_work::Vector{UInt8}` to `OSDWorkspace` (now `decoder_post.jl:16`)
and initialised it as `zeros(UInt8, num_check)` in `init_osd_workspace`
(`decoder_post.jl:52`).

### 1.2 Stray `copyto!(W.s_work, syndrome)` inside `_fast_osd!`

`_fast_osd!` (`decoder_post.jl:171`) has no `syndrome` parameter, but its body
contained `copyto!(W.s_work, syndrome)`, so any call raised
`UndefVarError: syndrome`. Both `osd_decode!` methods already populate
`W.s_work` before delegating — `fill!(W.s_work, 0x00)` in the classical form
(`decoder_post.jl:249`) and `copyto!(W.s_work, syndrome)` in the syndrome form
(`decoder_post.jl:263`) — so the line was simply deleted and the comment above it
updated to record that the caller owns `s_work`.

### 1.3 Flint matrices could not reach the new decoders

`init_osd_workspace`, `init_grand_workspace`, `init_wbf_workspace`,
`init_soft_workspace`, `layered_schedule` and `csr_of` were typed to
`AbstractMatrix`. `Nemo.fpMatrix` and `Nemo.FqMatrix` are **not**
`AbstractMatrix` subtypes, so a CodingTheory caller holding a Flint matrix —
which is the normal situation inside this package — got a `MethodError`.

Fix: purely additive `Union{fpMatrix, FqMatrix}` convenience methods that convert
once and delegate, plus a new `_Flint_matrix_to_Julia_support_matrix` helper in
`src/utils.jl` that extracts the 0/1 support pattern through `iszero` (so it
handles `FqMatrix` as well as `fpMatrix`, unlike the existing `nmod_mat`-based
converters). The `AbstractMatrix` methods were deliberately left untouched.

> **Constraint to preserve.** The `AbstractMatrix` and compressed-sparse-row
> signatures in `MP_decoders.jl` exist so that numpy/scipy arrays can be passed
> in from Python. Do not narrow them to Flint or Oscar types. Any new type
> support must be additive convenience methods, as above.

Also in `a4dfb8e`: `csr_of` and `serial_schedule` are public-facing but were
missing from the `export` list; they were added at `src/CodingTheory.jl:489`.

---

## 2. Layered min-sum on small dense high-rate matrices

### 2.1 Scope of the observation **[run]**

A sweep over four parity-check matrices, syndrome decoding through
`init_soft_workspace` / `decode!` with a constant channel LLR vector
`fill(log((1-p)/p), n)`, `max_iter = 100`, 400 trials for the first matrix and
200 for the rest, counting only nonzero syndromes. `conv` is the fraction that
reported convergence; `exact` is the fraction where the returned vector equalled
the injected error.

```
classical Hamming(7,4)        3 x 7    p=0.050  rank=3
  sum_product           flooding conv=1.00 exact=0.80 | layered conv=1.00 exact=0.91
  min_sum               flooding conv=1.00 exact=0.80 | layered conv=0.36 exact=0.32
  normalized_min_sum    flooding conv=1.00 exact=0.80 | layered conv=0.71 exact=0.64
  offset_min_sum        flooding conv=1.00 exact=0.80 | layered conv=0.36 exact=0.32
  min_sum_correction    flooding conv=1.00 exact=0.80 | layered conv=0.59 exact=0.53

classical LDPC (3,6) n=120   58 x 120  p=0.030  rank=58
  sum_product           flooding conv=0.99 | layered conv=0.98
  min_sum               flooding conv=0.97 | layered conv=0.93
  normalized_min_sum    flooding conv=0.97 | layered conv=0.97
  offset_min_sum        flooding conv=0.98 | layered conv=0.92
  min_sum_correction    flooding conv=0.98 | layered conv=0.98

quantum BB [[30,6,5]] H_X    15 x 30   p=0.030  rank=12 (dependent rows)
  min_sum               flooding conv=0.80 | layered conv=0.80   (sum_product 0.84 / 0.93)
quantum BB [[42,12,5]] H_X   21 x 42   p=0.030  rank=15 (dependent rows)
  min_sum               flooding conv=0.89 | layered conv=0.91   (sum_product 0.96 / 0.97)
```

**Min-sum is healthy under `:flooding` on every matrix tested, and healthy under
`:layered` on the sparse classical LDPC matrix and on both quantum matrices.**
The only dramatic degradation is the 3x7 Hamming matrix under `:layered`.
If you have seen an earlier note claiming a general min-sum failure in this
decoder, that framing is wrong — the effect is specific to the regime described
below.

### 2.2 Verdict: not a bug **[run]**

`MP_decoders.jl` implements layered min-sum correctly. Four independent checks:

**(a) An independent reference implementation reproduces the failure exactly.**
A deliberately naive textbook serial min-sum was written from scratch — dense
arrays, dictionaries, an explicit product of signs times the minimum of the
*other* magnitudes, the syndrome applied as a sign flip — sharing no code with
`MP_decoders.jl`. Over all 7 nonzero syndromes of the Hamming matrix, for both
`:min_sum` and `:offset_min_sum` (14 cases), the reference agreed with the
library on the convergence flag, the iteration count and every output bit,
including all 5 failures. This is the decisive evidence.

**(b) The optimized check update is exact.** Over 20,000 random check updates
(degree 1 to 9, ~10% of inputs exactly `0.0`, both syndrome polarities, all
three min-sum-family algorithms), the `Val(:single_pass)` rule at
`MP_decoders.jl:785` agreed with both the quadratic `Val(:sequential)` fold and
the `Val(:forward_backward)` fold to **0.0 absolute difference** — bit-identical,
not merely close. The single-pass min1/min2/negative-count optimization is
sound, including the tie handling and the exact-zero handling.

**(c) The failure is a degenerate exact-tie effect.** Perturbing the channel LLR
vector by a *relative* `1e-12` of Gaussian noise (200 repetitions x 7 syndromes)
raises layered min-sum convergence on the Hamming matrix from 0.286 to 0.907,
and the number is flat from `1e-12` all the way to `1e-1`. A defect in the
scheduling or in the arithmetic would not be removed by noise 12 orders of
magnitude below the signal.

**(d) The failure mode is a stationary point, not oscillation or overflow.**
Instrumented trace, `:min_sum`, `:layered`, syndrome `(1,0,0)`, channel LLR
`2.9444` on all 7 bits. The greedy colouring produces three singleton layers,
`[[1], [2], [3]]`, so the sweep is fully serial.

```
it1 c1 vars=[1,3,5,7] V2C=[ 2.9444  2.9444  2.9444  2.9444] C2V=[-2.9444 -2.9444 -2.9444 -2.9444]
it1 c2 vars=[2,3,6,7] V2C=[ 2.9444  0.0000  2.9444  0.0000] C2V=[ 0.0000  0.0000  0.0000  0.0000]
it1 c3 vars=[4,5,6,7] V2C=[ 2.9444  0.0000  2.9444  0.0000] C2V=[ 0.0000  0.0000  0.0000  0.0000]
it1  total_llrs=[0.0000 2.9444 0.0000 2.9444 0.0000 2.9444 0.0000]  bits=0000000
it2  ... byte-for-byte identical to it1 ...
```

The mechanism, in three steps:

1. Every input to check 1 has the *same* magnitude, because the channel LLR is
   constant. Min-sum's outgoing magnitude is the minimum over the other edges,
   so it is exactly the incoming magnitude, `2.9444`. Check 1 is unsatisfied, so
   `flip && (agg = -agg)` at `MP_decoders.jl:822` negates it.
2. The layered engine folds that straight back into the posterior
   (`MP_decoders.jl:967-970`), giving `total_llrs == 0.0` *exactly* on all four
   of check 1's variables. This is min-sum's magnitude overshoot: sum-product
   emits only `1.8532` here, leaving the posterior at `1.0913`.
3. A min-sum check with an incoming `0.0` emits `0.0` on every other edge, which
   is correct (an uninformative input makes the parity uninformative). Checks 2
   and 3 each see two zeros, so `min1 == min2 == 0` and *all* their outgoing
   messages are zero. Nothing changes on the next sweep. It is an exact fixed
   point of the iteration.

`:flooding` escapes because all three checks read the same snapshot, so the
`-2.9444` from check 1 is added to `+2.9444` contributions from checks 2 and 3
in the same posterior update instead of being applied alone.

### 2.3 The `min_sum` == `offset_min_sum` coincidence is real **[run]**

The two algorithms produced identical `conv=0.36 exact=0.32`. The offset **is**
applied on the layered path — the traced `C2V` for `:offset_min_sum` is
`-2.4444`, i.e. `2.9444 - 0.5`, against `-2.9444` for `:min_sum`. What happens is
that the offset merely postpones the exact cancellation by one layer:

```
it1 c1 V2C=[2.9444 x4]           C2V=[-2.4444 x4]   -> posteriors become exactly 0.5
it1 c2 V2C=[2.9444 0.5 2.9444 0.5] C2V=[0 0 0 0]    -> min1 = 0.5, and 0.5 - 0.5 = 0
```

`max(0.0, abs(agg) - β)` at `MP_decoders.jl:514-515` maps the resulting `0.5`
magnitudes to exactly `0.0`, reaching the same stationary point one step later
and the same hard decisions. Not a missing offset, and not luck.

`:normalized_min_sum` (α = 0.75) and `:min_sum_correction` do better (5/7 and 4/7
of the Hamming syndromes rather than 2/7) because a multiplicative scaling
cannot produce an exact cancellation against the channel LLR. `:normalized_min_sum`
still fails on 2 of the 7 syndromes: the trace shows the relevant posterior
decaying `0.7361 -> 0.3220 -> 0.0115 -> ...`, converging on the same attractor
from above without ever crossing zero.

### 2.4 The convergence test is not at fault **[run]**

On the 5 failing Hamming syndromes the hard decision at `max_iter` is the
all-zero vector, whose syndrome is `000` and therefore genuinely does not match
the nonzero target. `_syndrome_matches` (`MP_decoders.jl:838`) is not rejecting a
valid solution.

`oscillation = :active` **does** catch it, because
`_check_oscillation(::Val{:active}, ...)` at `MP_decoders.jl:1052` compares
against the previous two hard-decision vectors and a fixed point trivially
matches. All 5 failing syndromes return `(false, -2)` instead of running the
full 100 iterations. If you want the cheap early exit, that keyword is the way
to get it; it does not improve the decode.

### 2.5 Guidance

No algorithm change was made. The only edit to `MP_decoders.jl` alongside this
note is documentation: the caveat below was added to the `:layered` engine
docstring and cross-referenced from `decode!`. The regime that triggers the
stall is narrow and identifiable:

* a **constant (or near-constant) channel LLR vector**, which is exactly the
  standard syndrome-decoding setup with a uniform prior — this is what creates
  the exact ties;
* a **sequential schedule** (`:layered` on a densely-overlapping graph, or
  `:serial`), so that one check's overshoot lands on the posterior alone;
* a **small, dense, high-rate** matrix, so that the resulting zeros reach every
  remaining check within one sweep.

Recommendations:

1. On small dense high-rate matrices, use `:flooding`, or use `:sum_product`.
   Both were 7/7 on every Hamming syndrome; `:min_sum` under `:layered` was 2/7.
2. If min-sum under a sequential schedule is required there, prefer
   `:normalized_min_sum` (or `:min_sum_correction`) over `:min_sum` and
   `:offset_min_sum`. A multiplicative correction cannot cancel exactly; an
   additive offset can, and does.
3. Breaking the tie works and is cheap: a relative perturbation of `1e-12` on
   the channel LLRs restored Hamming layered min-sum from 0.286 to 0.907. If you
   ever add a built-in dither option, that is the mechanism.
4. This is worth almost nothing on sparse graphs. The same `1e-6` perturbation
   moved a 58x120 column-weight-3 LDPC matrix from 0.951 to 0.977 **[run]** —
   real but marginal, which is consistent with the sweep table.
5. `oscillation = :active` turns a 100-iteration stall into a 2-iteration
   failure with no loss of accuracy on this matrix.

---

## 3. `MP_decoders.jl` and `decoder_post.jl` are not wired together

**Verified [run] + [src].** There is no BP+OSD path anywhere in `src/`:

* `decoder_post.jl` contains no reference to `SoftDecisionWorkspace`,
  `init_soft_workspace` or `decode!`.
* No other file in `src/` references `OSDWorkspace`, `osd_decode!`,
  `grand_decode!` or `wbf_decode!` except the export block at
  `src/CodingTheory.jl:458-460`.
* Grepping `src/` for `bp_osd`, `BP_OSD`, `bposd` returns nothing.

The two files are included independently at `src/CodingTheory.jl:457` and
`src/CodingTheory.jl:485`.

### What the hand-off needs

The good news is that the types already line up, so the gap is a missing
convenience function rather than a redesign. Verified by running it:

* `SoftDecisionWorkspace.total_llrs` is a `Vector{Float64}` of length `num_var`,
  holding the posterior LLRs (positive = bit 0). `osd_decode!` takes exactly
  `total_llrs::Vector{Float64}`.
* `SoftDecisionWorkspace.target_syndrome` is a `Vector{UInt8}` of length
  `num_check`, already normalised to 0/1 by `load_soft_channel!`.
  `osd_decode!`'s syndrome method takes exactly `syndrome::Vector{UInt8}`.
* So the entire hand-off, today, is:

  ```julia
  W_soft = init_soft_workspace(H; schedule = :layered)
  W_osd  = init_osd_workspace(H)
  conv, _ = decode!(W_soft, llr; algorithm = :min_sum, schedule = :layered,
                    syndrome = s, max_iter = 100)
  e = conv ? copy(W_soft.current_bits) :
      osd_decode!(W_osd, W_soft.total_llrs, W_soft.target_syndrome;
                  method = :cs, order = 2)
  ```

  Run on the Hamming matrix this takes layered min-sum from 2/7 syndromes to
  **7/7 valid solutions** **[run]**.

What a first-class integration would have to add:

1. **Two workspaces, one matrix.** `init_osd_workspace` re-densifies `H` into its
   own `H_dense`, independently of the CSR form the soft workspace holds. A
   combined constructor should build both from one description of `H`.
2. **A success flag from OSD.** `osd_decode!` returns a bare `Vector{UInt8}`
   (`decoder_post.jl:243`, `decoder_post.jl:257`), unlike `grand_decode!`
   (`Tuple{Bool, Vector{UInt8}}`, `decoder_post.jl:318`) and `wbf_decode!`
   (`Tuple{Bool, Vector{UInt8}, Int}`, `decoder_post.jl:454`). A BP+OSD wrapper
   that wants to report success has to re-verify `H*e == s` itself.
3. **Aliasing care.** `osd_decode!`'s first argument aliases `W_soft.total_llrs`
   in the sketch above. That is safe as written — `_fast_osd!` only reads it —
   but it is not documented, and a future in-place normalisation inside OSD would
   silently corrupt the soft workspace.
4. **A caveat that is not obvious.** Feeding *stalled* BP posteriors to OSD is
   worse than feeding the channel LLRs. On the Hamming syndrome `(1,0,0)`, OSD
   from the raw channel LLRs returns the weight-1 solution `1000000`, but OSD
   from the stalled layered-min-sum posteriors returns the weight-3 solution
   `0010101` **[run]** — still valid, but a worse estimate, because the exact
   zeros of section 2.2 make the stalled bits look maximally unreliable and drag
   them out of the most-reliable basis. Any BP+OSD wrapper should decide
   deliberately which LLRs it hands to OSD.

---

## 4. Other issues noticed

Only items actually checked are listed.

### 4.1 Four exported names do not exist **[run]**

`src/CodingTheory.jl:486-487` exports `DecoderWorkspace`, `HardDecisionWorkspace`,
`init_hard_workspace` and `load_hard_channel!`. None of them are defined in
`MP_decoders.jl`; they live in `src/LDPC/MP_decoders_old.jl`, which is **not**
included anywhere (`src/CodingTheory.jl:485` includes `MP_decoders.jl` only).
`using CodingTheory` succeeds — Julia allows exporting an undefined name — but
`CodingTheory.init_hard_workspace` raises
`UndefVarError: init_hard_workspace not defined in CodingTheory`.

### 4.2 `test/LDPC/MP_decoders_test.jl` is stale and all of it errors **[run]**

Running the two relevant test items on `weight_dist` HEAD:

```
CodingTheory/test/LDPC/decoder_post_test.jl   |  20 pass  1 fail
CodingTheory/test/LDPC/MP_decoders_test.jl    |               5 error
```

All five `MP_decoders_test.jl` testsets error because the file is written against
the `MP_decoders_old.jl` API: `layered_schedule(H, schedule = :layered)`,
`W_soft.layers`, a three-value `decode!` return, `decimated_bits_values`,
`init_hard_workspace`. The current API is `layered_schedule(H)`,
`W.layer_ptr`/`W.layer_checks`, a two-value `decode!` return, and
`decimated_bits`/`decimated_values`. **This file needs rewriting; there is
currently no working test coverage of `MP_decoders.jl` at all.**

### 4.3 `grand_decode!` ranked by Hamming weight, not soft cost — now FIXED **[run]**

`test/LDPC/decoder_post_test.jl:77` expects `grand_decode!` to return
`[1,1,1,0,0,0,0]` for `llrs = [0.5, 0.8, -5.0, 5.0, 5.0, 5.0, 5.0]`; it returns
the all-zero codeword.

The test's expectation is the maximum-likelihood answer: flipping bits 1 and 2
costs `0.5 + 0.8 = 1.3` in soft distance, flipping bit 3 costs `5.0`.
`grand_decode!` (`decoder_post.jl:318`) searches strictly by **Hamming weight**
— all weight-1 patterns over the least-reliable set, then all weight-2, then all
weight-3 — and returns the first pattern that matches the syndrome, without ever
comparing soft costs. The weight-1 flip of bit 3 is found before any weight-2
pattern is tried, so the higher-cost answer wins.

This was independent of `a4dfb8e`, which only added an `init_grand_workspace`
convenience method and did not touch `grand_decode!` **[src]**. The same weakness
made `grand_decode!` inconsistent with `osd_decode!`, which *does* rank
candidates by soft distance (`decoder_post.jl:98-104`).

**Resolution.** `grand_decode!` now returns the minimum-soft-cost pattern in its
search space, and the test expectation was left untouched. The sweeps still
enumerate by Hamming weight, but they track the cheapest match rather than
returning the first, and each loop carries a lower bound: because `sortperm!`
leaves the least-reliable array nondecreasing, the cheapest pattern still
reachable from any index takes the next consecutive bits, so a subtree can be
pruned once `rel[i] + rel[i+1] + rel[i+2] >= best_cost`. That yields the same
provable optimality as cost-ordered (ORBGRAND-style) enumeration without a heap
or a sort.

Cost-ordered enumeration was prototyped and rejected on measurement: it pays an
unconditional sort over every candidate pattern on each call, which costs more
than the syndrome checks it avoids (298 patterns generated, costed and sorted to
skip ~80 checks at `max_lrb = 12`, and it degrades as the space grows). On a
512-bit rate-1/2 code the shipped branch-and-bound version is *faster than the
original buggy code* — 4.10 vs 5.03 us/call on a match, 16.92 vs 21.36 on the
no-match worst case at `max_lrb = 20` — and the decode loop is now
allocation-free at **0 bytes/call**, down from 4320, after switching
`sortperm!` to the fully in-place `alg = QuickSort`. Note that omitting the
`syndrome` keyword still allocates 320 bytes from the `zeros(UInt8, W.num_check)`
default in the signature; that is pre-existing and callers crossing the Python
boundary pass a syndrome anyway.

Verified by 500 randomised trials against the cost-ordered prototype (0 soft-cost
disagreements) and 400 trials against exhaustive brute force on a 16-bit code.
`decoder_post_test.jl` went from 20 pass / 1 fail to **65 pass / 0 fail / 0
error**, all additions additive. The new coverage includes an exhaustive-brute-force
cross-check and a field-by-field comparison of the array and Flint
`GRANDWorkspace` constructors over `fieldnames`, which will catch any future
added-but-uninitialised field — the failure mode that caused bug 1 above.

Still open in the same file: `_fast_osd!` builds `mrb_indices = Int[]` by `push!`
on every call, so `osd_decode!` is not allocation-free.

### 4.4 The layer partition has no effect on the numerical result **[run]**

`_fast_decode!(..., ::Val{:layered}, ...)` at `MP_decoders.jl:939` iterates
layers, and within each layer iterates its checks, but the inner loop is
strictly sequential: each check subtracts, updates and re-adds before the next
one starts. Nothing in the engine exploits the fact that checks within a layer
are independent.

Consequence: only the **order** in which `layer_checks` lists the checks matters,
not how they are grouped. Verified on a 58x120 matrix over 200 syndromes:
putting all 58 checks in a single layer in natural order gives bit-identical
output to `:serial` (58 singleton layers, same order) in **200/200** trials,
while the real conflict-free colouring — which reorders the checks — differs in
2/200.

This is not a bug, and the docstring at `MP_decoders.jl:913-938` is accurate
about what layering *means*. But the layer partition is currently pure metadata:
it costs a graph colouring at construction time and buys nothing until someone
parallelizes the inner loop. That is now stated in the docstring, along with the
section 2 caveat, so nobody assumes `:layered` is already doing something
`:serial` is not.

### 4.5 `erasures` sets LLRs to exactly `0.0`, which min-sum propagates **[run]**

`load_soft_channel!` sets `W.channel_llrs[v] = 0.0` for each erasure
(`MP_decoders.jl:440-442`). Combined with the section 2 mechanism, one erased bit
zeroes **every** outgoing message of each check it touches, on every min-sum
iteration until the posterior moves off zero. Directly observed: with
`erasures = [3]` on the Hamming matrix, one min-sum update of check 1
(variables 1, 3, 5, 7) emits `[-0.0, -2.9444, -0.0, -0.0]` — the three
non-erased edges get nothing.

This is arithmetically correct box-plus behaviour and sum-product does the same,
but for min-sum the effect is total rather than partial, and `erasures` is the
one API in this file that manufactures exact zeros on purpose. Worth a docstring
warning at minimum. **[hyp]** A small epsilon instead of `0.0` would avoid it
without changing the semantics meaningfully, but I have not measured whether
that is a net improvement.

Related: `_harden!` (`MP_decoders.jl:850`) resolves a posterior of exactly `0.0`
to bit **0** via `total_llrs[v] < 0.0`. That is a defensible convention but it
means an undetermined bit always decodes as 0 rather than being flagged.

### 4.6 `layered_schedule(H)` output cannot be fed back to `init_soft_workspace` **[src]**

`layered_schedule(H::AbstractMatrix)` at `MP_decoders.jl:212` returns the nested
form `Vector{Vector{Int}}`, but `init_soft_workspace`'s `layer_ptr` /
`layer_checks` keywords (`MP_decoders.jl:281-282`) want the flat CSR-style form
produced by the four-argument `layered_schedule`. There is no converter between
the two, so the obvious "compute the schedule, inspect it, then build the
workspace with it" workflow does not typecheck. A `flatten`/`unflatten` pair, or
accepting the nested form in `_validated_layers` (`MP_decoders.jl:382`), would
close this.

### 4.7 Unsupported `algorithm`, `method` and `order` values surface as `MethodError` **[run]**

`decode!` dispatches on `Val(algorithm)` without validating it, so
`algorithm = :not_an_algorithm` produces a `MethodError` on an internal
`_apply_boxplus` rather than an `ArgumentError` naming the supported values.
Same for `osd_decode!`: `order = 3` and `method = :typo` both give `MethodError`
on `_run_osd_sweeps!`. Contrast `schedule`, which *is* validated —
`init_soft_workspace` throws `ArgumentError("Unknown schedule $schedule")`
(`MP_decoders.jl:336`), and running `:layered` against a flooding workspace
throws a clear `ArgumentError` (`MP_decoders.jl:942`). The dispatch-based
algorithm selection is a deliberate zero-cost design; a one-line `in` check at
the `decode!` entry point would keep that and still give a decent error.

### 4.8 `decoder_post.jl` public functions take concrete types **[src]**

`osd_decode!`, `grand_decode!` and `wbf_decode!` are all typed
`total_llrs::Vector{Float64}` and `syndrome::Vector{UInt8}`
(`decoder_post.jl:243`, `:257`, `:318`, `:454`). `MP_decoders.jl` uses
`AbstractVector{<:Real}` and `AbstractVector{<:Integer}` throughout. The concrete
types rule out views, `Vector{Float32}`, and `Vector{Int}` syndromes without a
copy. Worth aligning if `decoder_post.jl` is ever exposed across the same
Python boundary as `MP_decoders.jl`.

### 4.9 Minor, `decoder_post.jl`

* `@simd` at `decoder_post.jl:88` is applied to a loop whose body is guarded by
  `if W.is_mrb[m_col]`, which is not a form `@simd` is specified to accept.
  The reduction is `⊻=`, which is associative and commutative, so reordering is
  harmless; in practice the compiler will simply decline to vectorize. It is
  misleading rather than wrong. **[src]**
* `_evaluate_pattern!` re-encodes with an inner loop over all `num_var` columns
  for each parity column, so it is O(num_var^2) per candidate pattern, and
  `_run_osd_sweeps!(Val(:standard), Val(2), ...)` evaluates O(|MRB|^2) patterns.
  Order-2 OSD on a code of even moderate length will be very slow. Restricting
  the inner loop to the MRB columns (precomputed once per decode) is the obvious
  fix. **[src]**
* The Gaussian elimination at `decoder_post.jl:192` runs `for col in num_var:-1:1`
  over reliability-sorted columns, i.e. it pivots on the *least* reliable columns
  first and leaves the most reliable ones free. That is the correct
  most-reliable-basis construction, and OSD order 2 recovers the exact weight-1
  error for all 7 Hamming syndromes **[run]** — noted only because the descending
  loop reads like a mistake at first glance.

### 4.10 `Pkg.test()` cannot run on `weight_dist` at all **[run]**

Unrelated to the decoders, but it blocks whole-suite verification, so anyone
working on this branch will hit it immediately.

`ext/JLD2Ext/JLD2Ext.jl:4` does
`import CodingTheory: TriangularColorCode488, TriangularColorCode666, ...`, but
`include("Quantum/misc_known_codes.jl")` is **commented out** at
`src/CodingTheory.jl:559`, so those names are never defined in the module. The
extension therefore fails to precompile:

```
✗ CodingTheory → JLD2Ext
ERROR: LoadError: UndefVarError: `TriangularColorCode488` not defined in `CodingTheory`
  @ ext/JLD2Ext/Quantum/misc_known_codes.jl:11
```

Because `JLD2` is a direct dependency of `test/Project.toml`, the extension loads
during `Pkg.test()` and the failure aborts the entire suite before any test
runs. Individual test items still run fine through `TestItemRunner` against the
package environment, which is how the decoder work here was verified.

This is **not** caused by any change described in this document: it is already
present at `a4dfb8e^`, and the same `include` is commented out there. It is the
same *class* of defect as section 4.1 — a declared/imported name with no
definition behind it — and it appears to be fallout from the branch's
in-progress state, since `Quantum/weight_dist.jl` and the whole Quantum
known-codes export block are commented out alongside it.

Two candidate fixes, both a judgement call for whoever owns the refactor:
re-enable `include("Quantum/misc_known_codes.jl")`, or guard/trim `JLD2Ext` so
it does not import names the package does not currently define. Deliberately not
attempted here, since uncommenting a large disabled section on someone else's
in-flight branch is likely to cascade.
