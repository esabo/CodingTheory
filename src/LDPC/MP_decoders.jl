
# Copyright (c) 2023 - 2026 Eric Sabo, Benjamin Ide. BSD-style license, see the
# LICENSE file in the root of the upstream repository; a copy is vendored
# alongside this file as LICENSE.CodingTheory.

# Box-plus identity. Folding a message with this returns the message unchanged
# for every rule below, and unlike `Inf` it cannot produce `Inf - Inf == NaN`
# further down. Used as the seed of the forward-backward accumulation, so a
# degree-1 check emits a saturated (fully determined) message rather than the
# uninformative zero upstream produced.
const _LLR_ID = 1.0e12

# The magnitude of a fully determined belief: what decimation pins a bit to, and
# what a degree-1 check emits. Deliberately far below `_LLR_ID` so that a
# saturated message still adds normally to the floating-point sum of the rest.
# Messages are otherwise NOT clamped, so that every non-degenerate check produces
# bit-identical output to upstream.
const _LLR_MAX = 1.0e3

# How close to +/-1 the product form's `tanh` values are allowed to get (FIX-16). This
# caps any message it can emit at `2 * atanh(_TANH_CLAMP)` ~ 28.3, which is why that rule
# is opt-in rather than the default -- see `_check_update!(::Val{:product}, ...)`.
const _TANH_CLAMP = 0.999999999999

# Check degree at which forward-backward accumulation starts to pay for itself.
# Measured, not guessed: see FIX-12.
const _FB_MIN_DEGREE = 6

# Shared empty inputs, so the common "no erasures, no manual decimation" call
# allocates nothing for its defaults.
const _NO_INDICES = Int[]

"""
    SoftDecisionWorkspace{T}

Every buffer the soft-decision decoder needs, allocated once for one fixed
parity-check matrix and reused for every syndrome.

The Tanner graph is stored as a flat edge list in check-major order. Edge `e`
connects check `c` (the unique `c` with `chk_ptr[c] <= e < chk_ptr[c + 1]`) to
variable `edge_var[e]`. The edges incident on variable `v` are
`var_edges[var_ptr[v]:var_ptr[v + 1] - 1]`.
"""
struct SoftDecisionWorkspace{T <: AbstractFloat}
    num_var::Int
    num_check::Int
    num_edges::Int
    max_check_degree::Int

    # Messages, one entry per edge.
    V2C::Vector{T}
    C2V::Vector{T}

    # Flat Tanner graph (FIX-10).
    chk_ptr::Vector{Int}
    edge_var::Vector{Int}
    var_ptr::Vector{Int}
    var_edges::Vector{Int}

    # Per-variable state.
    channel_llrs::Vector{T}
    total_llrs::Vector{T}
    current_bits::Vector{UInt8}
    is_decimated::Vector{Bool}

    # Per-check state.
    target_syndrome::Vector{UInt8}

    # Layer partition, flat: layer `l` is `layer_checks[layer_ptr[l]:layer_ptr[l + 1] - 1]`.
    # Empty `layer_ptr` means no partition is loaded and only flooding can run.
    layer_ptr::Vector{Int}
    layer_checks::Vector{Int}

    # Forward-backward scratch, length `max_check_degree` (FIX-12).
    fwd::Vector{T}

    # Oscillation history.
    prev_bits_1::Vector{UInt8}
    prev_bits_2::Vector{UInt8}
end

"""
    csr_of(H) -> (row_ptr, col_ind)

Compressed sparse row description of `H`, 1-based, with column indices ascending
within each row. A convenience for Julia-side callers and tests; FlamingPy passes
scipy's CSR arrays straight to [`init_soft_workspace`](@ref) instead, and this
function is the only thing in the file that is quadratic in the matrix size.
"""
function csr_of(H::AbstractMatrix)
    num_check, num_var = size(H)
    row_ptr = Vector{Int}(undef, num_check + 1)
    col_ind = Int[]
    row_ptr[1] = 1
    for c in 1:num_check
        for v in 1:num_var
            iszero(H[c, v]) || push!(col_ind, v)
        end
        row_ptr[c + 1] = length(col_ind) + 1
    end
    return row_ptr, col_ind
end

"""
    csr_of(H::Union{fpMatrix, FqMatrix}) -> (row_ptr, col_ind)

Flint-native form, for CodingTheory callers holding an Oscar matrix. The plain
array and compressed-sparse-row entry points throughout this file are reserved
for the Python caller, which never sees a Flint matrix, so the conversion to a
dense Julia matrix is confined to this thin layer -- and, since every use of it
is a workspace constructor, is paid once per matrix rather than once per decode.
"""
csr_of(H::Union{fpMatrix, FqMatrix}) = csr_of(_Flint_matrix_to_Julia_support_matrix(H))

"""
    layered_schedule(row_ptr, col_ind, num_check, num_var; base = 1)

Partition the checks into layers such that no two checks in a layer share a
variable, so that the checks of one layer can be updated in any order -- or all
at once -- without changing the result.

Returns `(layer_ptr, layer_checks)`, the flat form described in
[`SoftDecisionWorkspace`](@ref), with 1-based check indices.

Greedy: each check takes the smallest admissible existing layer, and opens a new
one only if every existing layer already holds a neighbour. Admissibility is read
off a stamp array in one pass over the check's two-step neighbourhood, which is
linear in the graph rather than upstream's rescan of every layer per check.

Pass `base = 0` for 0-based input arrays, as produced by scipy.

Reference: Mansour and Shanbhag, "Turbo decoder architectures for low-density
parity-check codes" (2002).
"""
function layered_schedule(row_ptr::AbstractVector{<:Integer}, col_ind::AbstractVector{<:Integer},
                          num_check::Integer, num_var::Integer; base::Integer = 1)
    num_check > 0 && num_var > 0 || throw(ArgumentError("Input matrix of improper dimension"))
    shift = 1 - Int(base)

    # Variable-to-check adjacency, counted then filled: two linear passes, no push!.
    var_deg = zeros(Int, num_var + 1)
    for e in eachindex(col_ind)
        var_deg[Int(col_ind[e]) + shift + 1] += 1
    end
    var_ptr = Vector{Int}(undef, num_var + 1)
    var_ptr[1] = 1
    for v in 1:num_var
        var_ptr[v + 1] = var_ptr[v] + var_deg[v + 1]
    end
    fill_at = copy(var_ptr)
    var_chks = Vector{Int}(undef, length(col_ind))
    for c in 1:num_check
        for i in (Int(row_ptr[c]) + shift):(Int(row_ptr[c + 1]) + shift - 1)
            v = Int(col_ind[i]) + shift
            var_chks[fill_at[v]] = c
            fill_at[v] += 1
        end
    end

    layer_of = zeros(Int, num_check)     # 0 until assigned
    layer_size = Int[]
    blocked_by = zeros(Int, num_check)   # stamp: which check blocked this layer
    for c in 1:num_check
        # Stamp every layer that already holds a check sharing a variable with c.
        for i in (Int(row_ptr[c]) + shift):(Int(row_ptr[c + 1]) + shift - 1)
            v = Int(col_ind[i]) + shift
            for j in var_ptr[v]:(var_ptr[v + 1] - 1)
                other = var_chks[j]
                l = layer_of[other]
                l == 0 || (blocked_by[l] = c)
            end
        end

        # Smallest admissible layer, so the partition stays balanced.
        best, best_size = 0, typemax(Int)
        for l in eachindex(layer_size)
            if blocked_by[l] != c && layer_size[l] < best_size
                best, best_size = l, layer_size[l]
            end
        end
        if best == 0
            push!(layer_size, 1)
            layer_of[c] = length(layer_size)
        else
            layer_size[best] += 1
            layer_of[c] = best
        end
    end

    # Flatten, checks ascending within each layer.
    num_layers = length(layer_size)
    layer_ptr = Vector{Int}(undef, num_layers + 1)
    layer_ptr[1] = 1
    for l in 1:num_layers
        layer_ptr[l + 1] = layer_ptr[l] + layer_size[l]
    end
    layer_checks = Vector{Int}(undef, num_check)
    at = copy(layer_ptr)
    for c in 1:num_check
        l = layer_of[c]
        layer_checks[at[l]] = c
        at[l] += 1
    end
    return layer_ptr, layer_checks
end

"""
    layered_schedule(H) -> Vector{Vector{Int}}

Nested-vector form, for interactive use and for comparison against upstream.
"""
function layered_schedule(H::AbstractMatrix)
    row_ptr, col_ind = csr_of(H)
    layer_ptr, layer_checks = layered_schedule(row_ptr, col_ind, size(H, 1), size(H, 2))
    return [layer_checks[layer_ptr[l]:(layer_ptr[l + 1] - 1)] for l in 1:(length(layer_ptr) - 1)]
end

"""
    layered_schedule(H::Union{fpMatrix, FqMatrix}) -> Vector{Vector{Int}}

Flint-native form, as in [`csr_of`](@ref).
"""
layered_schedule(H::Union{fpMatrix, FqMatrix}) =
    layered_schedule(_Flint_matrix_to_Julia_support_matrix(H))

"""
    serial_schedule(num_check) -> (layer_ptr, layer_checks)

One check per layer: the fully serial schedule.
"""
function serial_schedule(num_check::Integer)
    return collect(1:(Int(num_check) + 1)), collect(1:Int(num_check))
end

"""
    balance_of_layered_schedule(layer_ptr)

Ratio of the largest layer to the smallest. 1 means every layer is the same size,
which is the best case for a parallel implementation of a layered schedule.

Reference: Layered decoding of quantum LDPC codes.
"""
function balance_of_layered_schedule(layer_ptr::AbstractVector{<:Integer})
    length(layer_ptr) >= 2 || throw(ArgumentError("Schedule cannot be empty"))
    smallest, largest = typemax(Int), 0
    for l in 1:(length(layer_ptr) - 1)
        len = Int(layer_ptr[l + 1]) - Int(layer_ptr[l])
        len > 0 || throw(ArgumentError("Schedule cannot contain an empty layer"))
        len < smallest && (smallest = len)
        len > largest && (largest = len)
    end
    return largest / smallest
end

function balance_of_layered_schedule(sch::AbstractVector{<:AbstractVector{<:Integer}})
    isempty(sch) && throw(ArgumentError("Schedule cannot be empty"))
    any(isempty, sch) && throw(ArgumentError("Schedule cannot contain an empty layer"))
    lengths = map(length, sch)
    return maximum(lengths) / minimum(lengths)
end

"""
    init_soft_workspace(row_ptr, col_ind, num_check, num_var; schedule, layer_ptr, layer_checks, base)

Allocate the decoder workspace for the parity-check matrix given in compressed
sparse row form. Linear in the number of edges (FIX-9).

`row_ptr` and `col_ind` are exactly scipy's `H.indptr` and `H.indices` for a
`csr_matrix`; pass `base = 0` for those, which is the default, since the Python
caller is the hot path.

`schedule` selects the layer partition to precompute: `:flooding` builds none,
`:layered` colours the graph, `:serial` puts one check per layer. Passing
`layer_ptr` and `layer_checks` supplies a partition directly -- computed once on
the Python side, cached, and reused across processes -- and overrides `schedule`.
"""
function init_soft_workspace(row_ptr::AbstractVector{<:Integer},
                             col_ind::AbstractVector{<:Integer},
                             num_check::Integer, num_var::Integer;
                             schedule::Symbol = :flooding,
                             layer_ptr::AbstractVector{<:Integer} = _NO_INDICES,
                             layer_checks::AbstractVector{<:Integer} = _NO_INDICES,
                             base::Integer = 0)
    num_check = Int(num_check)
    num_var = Int(num_var)
    num_check > 0 && num_var > 0 || throw(ArgumentError("Input matrix of improper dimension"))
    length(row_ptr) == num_check + 1 ||
        throw(ArgumentError("row_ptr must have num_check + 1 entries"))
    shift = 1 - Int(base)

    num_edges = length(col_ind)
    chk_ptr = Vector{Int}(undef, num_check + 1)
    edge_var = Vector{Int}(undef, num_edges)
    max_check_degree = 0
    for c in 1:(num_check + 1)
        chk_ptr[c] = Int(row_ptr[c]) + shift
    end
    for c in 1:num_check
        deg = chk_ptr[c + 1] - chk_ptr[c]
        deg > max_check_degree && (max_check_degree = deg)
    end
    for e in 1:num_edges
        v = Int(col_ind[e]) + shift
        1 <= v <= num_var || throw(ArgumentError("Column index $v out of range"))
        edge_var[e] = v
    end

    # Variable-side gather list, by counting sort on edge_var: two linear passes.
    var_ptr = Vector{Int}(undef, num_var + 1)
    counts = zeros(Int, num_var)
    for e in 1:num_edges
        counts[edge_var[e]] += 1
    end
    var_ptr[1] = 1
    for v in 1:num_var
        var_ptr[v + 1] = var_ptr[v] + counts[v]
    end
    var_edges = Vector{Int}(undef, num_edges)
    at = copy(var_ptr)
    for e in 1:num_edges
        v = edge_var[e]
        var_edges[at[v]] = e
        at[v] += 1
    end

    # Layer partition: supplied, precomputed, or absent.
    lay_ptr, lay_checks = if !isempty(layer_ptr)
        _validated_layers(layer_ptr, layer_checks, num_check)
    elseif schedule === :layered || schedule === :semiserial
        layered_schedule(chk_ptr, edge_var, num_check, num_var; base = 1)
    elseif schedule === :serial
        serial_schedule(num_check)
    elseif schedule === :flooding || schedule === :parallel
        Int[], Int[]
    else
        throw(ArgumentError("Unknown schedule $schedule"))
    end

    return SoftDecisionWorkspace{Float64}(
        num_var, num_check, num_edges, max_check_degree,
        zeros(Float64, num_edges), zeros(Float64, num_edges),
        chk_ptr, edge_var, var_ptr, var_edges,
        zeros(Float64, num_var), zeros(Float64, num_var), zeros(UInt8, num_var),
        zeros(Bool, num_var),
        zeros(UInt8, num_check),
        lay_ptr, lay_checks,
        Vector{Float64}(undef, max_check_degree),
        fill(0xFF, num_var), fill(0xFF, num_var),
    )
end

"""
    init_soft_workspace(H; schedule, layer_ptr, layer_checks)

Convenience constructor from any `AbstractMatrix`. Densely scans `H`, so prefer
the compressed-sparse-row form for anything large.
"""
function init_soft_workspace(H::AbstractMatrix; schedule::Symbol = :flooding,
                             layer_ptr::AbstractVector{<:Integer} = _NO_INDICES,
                             layer_checks::AbstractVector{<:Integer} = _NO_INDICES)
    row_ptr, col_ind = csr_of(H)
    return init_soft_workspace(row_ptr, col_ind, size(H, 1), size(H, 2);
                               schedule = schedule, layer_ptr = layer_ptr,
                               layer_checks = layer_checks, base = 1)
end

"""
    init_soft_workspace(H::Union{fpMatrix, FqMatrix}; schedule, layer_ptr, layer_checks)

Flint-native form, as in [`csr_of`](@ref).
"""
init_soft_workspace(H::Union{fpMatrix, FqMatrix}; schedule::Symbol = :flooding,
                    layer_ptr::AbstractVector{<:Integer} = _NO_INDICES,
                    layer_checks::AbstractVector{<:Integer} = _NO_INDICES) =
    init_soft_workspace(_Flint_matrix_to_Julia_support_matrix(H); schedule = schedule,
                        layer_ptr = layer_ptr, layer_checks = layer_checks)

"""
Validate and copy a caller-supplied layer partition: every check exactly once,
no empty layers.
"""
function _validated_layers(layer_ptr::AbstractVector{<:Integer},
                           layer_checks::AbstractVector{<:Integer}, num_check::Int)
    length(layer_checks) == num_check ||
        throw(ArgumentError("layer_checks must list all $num_check checks, got " *
                            "$(length(layer_checks))"))
    ptr = Vector{Int}(undef, length(layer_ptr))
    for i in eachindex(layer_ptr)
        ptr[i] = Int(layer_ptr[i])
    end
    first(ptr) == 1 || throw(ArgumentError("layer_ptr must start at 1"))
    last(ptr) == num_check + 1 || throw(ArgumentError("layer_ptr must end at num_check + 1"))
    for l in 1:(length(ptr) - 1)
        ptr[l + 1] > ptr[l] || throw(ArgumentError("Schedule cannot contain an empty layer"))
    end
    seen = falses(num_check)
    checks = Vector{Int}(undef, num_check)
    for i in 1:num_check
        c = Int(layer_checks[i])
        1 <= c <= num_check || throw(ArgumentError("Check index $c out of range"))
        seen[c] && throw(ArgumentError("Check $c appears in more than one layer"))
        seen[c] = true
        checks[i] = c
    end
    return ptr, checks
end

"""
    load_soft_channel!(W, LLR_in; syndrome, erasures, decimated_bits, decimated_values)

Load one channel realization into the workspace: channel LLRs, target syndrome,
erasures (neutral belief) and manually decimated bits (pinned belief). Resets the
messages, so a workspace can be reused for an unrelated syndrome.

Allocation-free: the defaults are shared empty constants and every argument is
copied into a pre-allocated buffer (FIX-14).
"""
function load_soft_channel!(W::SoftDecisionWorkspace{Float64}, LLR_in::AbstractVector{<:Real};
                            syndrome::AbstractVector{<:Integer} = _NO_INDICES,
                            erasures::AbstractVector{<:Integer} = _NO_INDICES,
                            decimated_bits::AbstractVector{<:Integer} = _NO_INDICES,
                            decimated_values::AbstractVector{<:Integer} = _NO_INDICES)
    length(LLR_in) == W.num_var ||
        throw(ArgumentError("Expected $(W.num_var) LLRs, got $(length(LLR_in))"))

    @inbounds begin
        copyto!(W.channel_llrs, LLR_in)

        if isempty(syndrome)
            fill!(W.target_syndrome, 0x00)
        else
            length(syndrome) == W.num_check ||
                throw(ArgumentError("Expected $(W.num_check) syndrome bits, " *
                                    "got $(length(syndrome))"))
            for c in 1:W.num_check
                W.target_syndrome[c] = iszero(syndrome[c]) ? 0x00 : 0x01
            end
        end

        for v in erasures
            W.channel_llrs[v] = 0.0
        end

        fill!(W.is_decimated, false)
        length(decimated_bits) == length(decimated_values) ||
            throw(ArgumentError("decimated_bits and decimated_values must have equal length"))
        for i in eachindex(decimated_bits)
            v = Int(decimated_bits[i])
            W.is_decimated[v] = true
            W.channel_llrs[v] = iszero(decimated_values[i]) ? _LLR_MAX : -_LLR_MAX
        end

        copyto!(W.total_llrs, W.channel_llrs)
        fill!(W.C2V, 0.0)
        fill!(W.V2C, 0.0)
        fill!(W.prev_bits_1, 0xFF)
        fill!(W.prev_bits_2, 0xFF)
    end
    return W
end

# ==============================================================================
# BOX-PLUS OPERATORS AND THEIR PER-ALGORITHM TRAITS
# ==============================================================================

"""
The exact box-plus operator (Jacobian logarithm) for sum-product. Mathematically
the tanh rule, but numerically bulletproof against NaNs.
"""
@inline function boxplus_exact(x::Float64, y::Float64)
    base = sign(x) * sign(y) * min(abs(x), abs(y))
    corr = log1p(exp(-abs(x + y))) - log1p(exp(-abs(x - y)))
    return base + corr
end

"""The min-sum box-plus operator."""
@inline boxplus_minsum(x::Float64, y::Float64) = sign(x) * sign(y) * min(abs(x), abs(y))

"""The min-sum box-plus operator with a low-complexity correction term."""
@inline function boxplus_minsum_correction(x::Float64, y::Float64)
    base = sign(x) * sign(y) * min(abs(x), abs(y))
    sum_abs = abs(x + y)
    diff_abs = abs(x - y)
    corr = 0.0
    if sum_abs < 2.0 && diff_abs > 2.0 * sum_abs
        corr = 0.5
    elseif diff_abs < 2.0 && sum_abs > 2.0 * diff_abs
        corr = -0.5
    end
    return base + corr
end

@inline _apply_boxplus(::Val{:sum_product}, x::Float64, y::Float64) = boxplus_exact(x, y)
@inline _apply_boxplus(::Val{:min_sum}, x::Float64, y::Float64) = boxplus_minsum(x, y)
@inline _apply_boxplus(::Val{:min_sum_correction}, x::Float64, y::Float64) =
    boxplus_minsum_correction(x, y)
# Normalized and offset min-sum fold exactly like min-sum and differ only afterwards.
@inline _apply_boxplus(::Val{:normalized_min_sum}, x::Float64, y::Float64) = boxplus_minsum(x, y)
@inline _apply_boxplus(::Val{:offset_min_sum}, x::Float64, y::Float64) = boxplus_minsum(x, y)

@inline _apply_boxplus(::Val{:sum_product_fast}, x::Float64, y::Float64) = boxplus_exact(x, y)
@inline _apply_boxplus(::Val{:min_sum_correction_fast}, x::Float64, y::Float64) =
    boxplus_minsum_correction(x, y)

@inline _apply_post_process(::Val{:sum_product}, agg::Float64, α::Float64, β::Float64) = agg
@inline _apply_post_process(::Val{:sum_product_fast}, agg::Float64, α::Float64,
                            β::Float64) = agg
@inline _apply_post_process(::Val{:min_sum}, agg::Float64, α::Float64, β::Float64) = agg
@inline _apply_post_process(::Val{:min_sum_correction}, agg::Float64, α::Float64, β::Float64) = agg
@inline _apply_post_process(::Val{:min_sum_correction_fast}, agg::Float64, α::Float64,
                            β::Float64) = agg
@inline _apply_post_process(::Val{:normalized_min_sum}, agg::Float64, α::Float64, β::Float64) =
    agg * α
@inline _apply_post_process(::Val{:offset_min_sum}, agg::Float64, α::Float64, β::Float64) =
    sign(agg) * max(0.0, abs(agg) - β)

# ==============================================================================
# CHECK-NODE UPDATE
# ==============================================================================

# Which fold a rule may use. Box-plus is associative for four of the five rules,
# so their all-but-one aggregates can be shared between edges. The low-complexity
# correction is NOT associative -- its +/-0.5 term is a threshold test on the pair
# being folded -- so re-associating it visibly changes the messages. It keeps the
# quadratic fold, and since it is only ever used on the low-degree checks this
# costs nothing in practice.
@inline _fold_trait(::Val) = Val(:forward_backward)
@inline _fold_trait(::Val{:min_sum_correction}) = Val(:sequential)
# The min-sum family needs no fold at all. One pass replaces the 2*deg folds.
# This is bit-identical, not merely equivalent:
# `min` and sign parity are exact in floating point, so the two routes agree on every input.
@inline _fold_trait(::Val{:sum_product_fast}) = Val(:product)
@inline _fold_trait(::Val{:min_sum_correction_fast}) = Val(:single_pass_corrected)
@inline _fold_trait(::Val{:min_sum}) = Val(:single_pass)
@inline _fold_trait(::Val{:normalized_min_sum}) = Val(:single_pass)
@inline _fold_trait(::Val{:offset_min_sum}) = Val(:single_pass)

"""
    _check_update!(W, algo, c, α, β)

Recompute every outgoing message of check `c` from the incoming ones, writing them
into `W.C2V`.
"""
@inline _check_update!(W::SoftDecisionWorkspace{Float64}, algo::Val, c::Int,
                       α::Float64, β::Float64) =
    _check_update!(_fold_trait(algo), W, algo, c, α, β)

"""
Forward-backward accumulation (FIX-12): `fwd[i]` folds the edges before `i` and a
running suffix folds those after it, so each message costs one further fold rather
than `deg - 1`. Exact by associativity, and bit-identical for the min-sum family,
where the fold is only a `min` and a product of signs.
"""
@inline function _check_update!(::Val{:forward_backward}, W::SoftDecisionWorkspace{Float64},
                                algo::Val, c::Int, α::Float64, β::Float64)
    @inbounds begin
        lo = W.chk_ptr[c]
        hi = W.chk_ptr[c + 1] - 1
        deg = hi - lo + 1
        deg > 0 || return nothing
        flip = W.target_syndrome[c] == 0x01

        # A degree-1 check determines its variable outright, so it emits a
        # saturated belief. Upstream emitted 0.0 here and threw that away.
        if deg == 1
            W.C2V[lo] = flip ? -_LLR_MAX : _LLR_MAX
            return nothing
        end

        # Below `_FB_MIN_DEGREE` the quadratic fold is measurably faster: it saves
        # no more than a couple of box-plus evaluations and the prefix array costs
        # a round trip to memory. It is also what upstream did, so taking it on the
        # low-degree checks that FlamingPy's FT codes are made of keeps this file
        # bit-identical to upstream there rather than merely equivalent.
        if deg < _FB_MIN_DEGREE
            return _check_update!(Val(:sequential), W, algo, c, α, β)
        end

        acc = _LLR_ID
        for i in 1:deg
            W.fwd[i] = acc
            acc = _apply_boxplus(algo, acc, W.V2C[lo + i - 1])
        end
        acc = _LLR_ID
        for i in deg:-1:1
            agg = _apply_boxplus(algo, W.fwd[i], acc)
            acc = _apply_boxplus(algo, acc, W.V2C[lo + i - 1])
            agg = _apply_post_process(algo, agg, α, β)
            # Inject the syndrome: an odd-parity check inverts its messages.
            flip && (agg = -agg)
            W.C2V[lo + i - 1] = agg
        end
    end
    return nothing
end

"""
Quadratic fold, kept bit-identical to upstream for the non-associative rule.
"""
@inline function _check_update!(::Val{:sequential}, W::SoftDecisionWorkspace{Float64},
                                algo::Val, c::Int, α::Float64, β::Float64)
    @inbounds begin
        lo = W.chk_ptr[c]
        hi = W.chk_ptr[c + 1] - 1
        deg = hi - lo + 1
        deg > 0 || return nothing
        flip = W.target_syndrome[c] == 0x01

        if deg == 1
            W.C2V[lo] = flip ? -_LLR_MAX : _LLR_MAX
            return nothing
        end

        for e_out in lo:hi
            agg = 0.0
            first_val = true
            for e_in in lo:hi
                e_in == e_out && continue
                if first_val
                    agg = W.V2C[e_in]
                    first_val = false
                else
                    agg = _apply_boxplus(algo, agg, W.V2C[e_in])
                end
            end
            agg = _apply_post_process(algo, agg, α, β)
            flip && (agg = -agg)
            W.C2V[e_out] = agg
        end
    end
    return nothing
end

"""
Single pass with an aggregate-level correction (FIX-17), for `:min_sum_correction_fast`.

`:min_sum_correction` is the one rule that cannot use a shared fold, because its +/-0.5 term
is a threshold test on the specific pair being folded, so re-associating changes the answer.
That forces `deg * (deg - 2)` box-plus evaluations per check, and it is the most expensive
rule in the file by a factor of four.

This rule moves the correction OUT of the fold: the magnitude comes from the same min1/min2
single pass as `:min_sum`, and the threshold test is applied once, comparing that magnitude
against the next-smallest one that the same edge sees. Three running minima instead of two.

This is a DIFFERENT DECODER, not an optimization of the existing one -- the messages it
produces differ, so it has its own symbol and needs its own logical-error-rate validation
before it is used in place of `:min_sum_correction`.
"""
@inline function _check_update!(::Val{:single_pass_corrected},
                               W::SoftDecisionWorkspace{Float64}, algo::Val, c::Int,
                               α::Float64, β::Float64)
    @inbounds begin
        lo = W.chk_ptr[c]
        hi = W.chk_ptr[c + 1] - 1
        deg = hi - lo + 1
        deg > 0 || return nothing
        flip = W.target_syndrome[c] == 0x01

        if deg == 1
            W.C2V[lo] = flip ? -_LLR_MAX : _LLR_MAX
            return nothing
        end

        min1 = Inf
        min2 = Inf
        min3 = Inf
        argmin1 = lo
        negatives = 0
        for e in lo:hi
            m = W.V2C[e]
            m < 0.0 && (negatives += 1)
            a = abs(m)
            if a < min1
                min3 = min2
                min2 = min1
                min1 = a
                argmin1 = e
            elseif a < min2
                min3 = min2
                min2 = a
            elseif a < min3
                min3 = a
            end
        end

        for e in lo:hi
            own = e == argmin1
            base = own ? min2 : min1
            other = own ? min3 : min2
            # The same threshold test as `boxplus_minsum_correction`, applied once to the
            # aggregate rather than inside each pairwise fold. `base` and `other` are both
            # magnitudes, so the sum and difference of the underlying signed pair reduce to
            # these two expressions.
            total = base + other
            spread = abs(other - base)
            corr = 0.0
            if total < 2.0 && spread > 2.0 * total
                corr = 0.5
            elseif spread < 2.0 && total > 2.0 * spread
                corr = -0.5
            end
            magnitude = max(0.0, base + corr)
            others_negative = negatives - (W.V2C[e] < 0.0 ? 1 : 0)
            agg = isodd(others_negative) ? -magnitude : magnitude
            agg = _apply_post_process(algo, agg, α, β)
            flip && (agg = -agg)
            W.C2V[e] = agg
        end
    end
    return nothing
end

"""
Product form of the exact box-plus, which is used for `:sum_product_fast`.

Algebraically the same rule as `:sum_product`: `tanh(boxplus/2)` is the product of the
incoming `tanh(m/2)`. Arithmetically much cheaper. The Jacobian-logarithm form spends
2 `exp` + 2 `log1p` on EVERY fold, and forward-backward performs about `2 * deg` folds per
check -- roughly `8 * deg` transcendentals. This form spends one `tanh` per edge going in and
one `atanh` per edge coming out, `2 * deg` total, with the forward-backward reduction running
over plain multiplies. Measured 4.3x faster per edge on a degree-6 check.

DYNAMIC RANGE -- the reason this is a separate rule and not a replacement. The `tanh` values
are clamped away from +/-1 so `atanh` stays finite, which caps any outgoing message at
`2 * atanh(_TANH_CLAMP)` ~ 28.3. The log-domain form has no such cap, and `_LLR_MAX` is 1e3:
a degree-1 check or a decimated bit injects a belief far above the cap and this rule will
flatten it. Harmless when beliefs stay small (measured bit-identical to `:sum_product` over
20 iterations on a degree-6 code at p = 0.15, where messages peaked at 0.90), wrong when they
do not. Use `:sum_product` when saturated beliefs are in play.
"""
@inline function _check_update!(::Val{:product}, W::SoftDecisionWorkspace{Float64},
                                algo::Val, c::Int, α::Float64, β::Float64)
    @inbounds begin
        lo = W.chk_ptr[c]
        hi = W.chk_ptr[c + 1] - 1
        deg = hi - lo + 1
        deg > 0 || return nothing
        flip = W.target_syndrome[c] == 0x01

        if deg == 1
            W.C2V[lo] = flip ? -_LLR_MAX : _LLR_MAX
            return nothing
        end

        for i in 1:deg
            W.fwd[i] = clamp(tanh(0.5 * W.V2C[lo + i - 1]), -_TANH_CLAMP, _TANH_CLAMP)
        end

        # Prefix products are stashed in the output slots, so no extra buffer is needed;
        # the backward sweep consumes each one before overwriting it.
        prefix = 1.0
        for i in 1:deg
            W.C2V[lo + i - 1] = prefix
            prefix *= W.fwd[i]
        end
        suffix = 1.0
        for i in deg:-1:1
            agg = 2.0 * atanh(W.C2V[lo + i - 1] * suffix)
            suffix *= W.fwd[i]
            agg = _apply_post_process(algo, agg, α, β)
            flip && (agg = -agg)
            W.C2V[lo + i - 1] = agg
        end
    end
    return nothing
end

"""
Single pass for the min-sum family.

The all-but-one min-sum aggregate on edge `e` is the product of every other edge's sign
times the smallest of every other edge's magnitude. Both are available from one sweep: the
two smallest magnitudes `min1 <= min2`, the index that attained `min1`, and the parity of
the negative count. The magnitude an edge sees is then `min2` if it owns `min1` and `min1`
otherwise -- so `deg` folds and `deg` further folds collapse into `deg` compares.

Bit-identical to the forward-backward fold, including ties (`min1` repeated makes
`min2 == min1`, which is what the fold returns for both copies) and exact zeros. Signs are
carried as a NEGATIVE COUNT rather than a product of `sign` calls, because `sign(0.0)` is
`0.0`: a product would zero out the aggregate that an incoming zero is excluded from, while
the fold -- which never sees that zero -- would not. Parity over the other edges is what the
fold actually computes.
"""
@inline function _check_update!(::Val{:single_pass}, W::SoftDecisionWorkspace{Float64},
                                algo::Val, c::Int, α::Float64, β::Float64)
    @inbounds begin
        lo = W.chk_ptr[c]
        hi = W.chk_ptr[c + 1] - 1
        deg = hi - lo + 1
        deg > 0 || return nothing
        flip = W.target_syndrome[c] == 0x01

        if deg == 1
            W.C2V[lo] = flip ? -_LLR_MAX : _LLR_MAX
            return nothing
        end

        min1 = Inf
        min2 = Inf
        argmin1 = lo
        negatives = 0
        for e in lo:hi
            m = W.V2C[e]
            m < 0.0 && (negatives += 1)
            a = abs(m)
            if a < min1
                min2 = min1
                min1 = a
                argmin1 = e
            elseif a < min2
                min2 = a
            end
        end

        for e in lo:hi
            # Exclude this edge from both reductions.
            magnitude = e == argmin1 ? min2 : min1
            others_negative = negatives - (W.V2C[e] < 0.0 ? 1 : 0)
            agg = isodd(others_negative) ? -magnitude : magnitude
            agg = _apply_post_process(algo, agg, α, β)
            flip && (agg = -agg)
            W.C2V[e] = agg
        end
    end
    return nothing
end

# ==============================================================================
# CONVERGENCE
# ==============================================================================

"""
    _syndrome_matches(W) -> Bool

Whether the current hard decisions reproduce the target syndrome.
"""
@inline function _syndrome_matches(W::SoftDecisionWorkspace{Float64})
    @inbounds for c in 1:W.num_check
        syn = W.target_syndrome[c]
        for e in W.chk_ptr[c]:(W.chk_ptr[c + 1] - 1)
            syn ⊻= W.current_bits[W.edge_var[e]]
        end
        syn == 0x00 || return false
    end
    return true
end

"""Refresh the hard decisions from the posteriors. Positive LLR means bit 0."""
@inline function _harden!(W::SoftDecisionWorkspace{Float64})
    @inbounds @simd for v in 1:W.num_var
        W.current_bits[v] = W.total_llrs[v] < 0.0 ? 0x01 : 0x00
    end
    return nothing
end

# ==============================================================================
# SOFT DECISION ENGINE: FLOODING
# ==============================================================================

function _fast_decode!(W::SoftDecisionWorkspace{Float64}, algo::Val, ::Val{:flooding},
                       decimation_type::Val, osc_type::Val, max_iter::Int,
                       α::Float64, β::Float64, dec_thresh::Float64, dec_rounds::Int)
    @inbounds for iter in 1:max_iter
        # 1. Variable nodes: every edge, from the same posterior snapshot.
        for v in 1:W.num_var
            tot = W.total_llrs[v]
            for i in W.var_ptr[v]:(W.var_ptr[v + 1] - 1)
                e = W.var_edges[i]
                W.V2C[e] = tot - W.C2V[e]
            end
        end

        # 2. Check nodes.
        for c in 1:W.num_check
            _check_update!(W, algo, c, α, β)
        end

        # 3. Posteriors and hard decisions.
        for v in 1:W.num_var
            if !W.is_decimated[v]
                tot = W.channel_llrs[v]
                for i in W.var_ptr[v]:(W.var_ptr[v + 1] - 1)
                    tot += W.C2V[W.var_edges[i]]
                end
                W.total_llrs[v] = tot
            end
        end
        _harden!(W)

        # 4. Hooks and convergence.
        _apply_decimation!(decimation_type, W, iter, dec_thresh, dec_rounds)
        _syndrome_matches(W) && return true, iter

        if _check_oscillation(osc_type, W)
            if decimation_type !== Val(:none) && decimation_type !== Val(:manual)
                _apply_decimation!(decimation_type, W, iter, dec_thresh, dec_rounds; force = true)
                _wipe_history!(osc_type, W)
            else
                return false, -iter
            end
        else
            _update_history!(osc_type, W)
        end
    end
    return false, max_iter
end

# ==============================================================================
# SOFT DECISION ENGINE: LAYERED
# ==============================================================================

"""
Layered schedule. Sweeps the layers of `W.layer_ptr`; within a layer the
checks share no variable, so their updates are independent, and each layer's new
messages are folded into the posteriors before the next layer reads them. That
immediate feedback is the point of layering: it typically halves the iteration
count against flooding for the same work per iteration.

The inner loop is sequential, so the partition itself does not change the
numbers -- only the ORDER `layer_checks` lists the checks in does. It is there
for a future parallel implementation of a layer, and to document independence.

MIN-SUM CAVEAT. The min-sum family can reach an exact stationary point here on a
SMALL DENSE HIGH-RATE matrix decoded from a CONSTANT channel LLR vector, which is
the usual syndrome-decoding setup. When every input to an unsatisfied check has
the same magnitude, min-sum's outgoing magnitude equals it exactly, this engine
folds it straight back, and the posterior lands on exactly 0.0; a min-sum check
with a 0.0 input emits 0.0 on every other edge, so the zeros spread and the state
repeats forever. `:offset_min_sum` reaches the same point one layer later, since
`max(0, |agg| - β)` maps the surviving `β` magnitudes to zero. On the 3x7
Hamming(7,4) matrix this costs 5 of the 7 nonzero syndromes.

Prefer `:flooding` or `:sum_product` in that regime; `:normalized_min_sum`, whose
correction is multiplicative and so cannot cancel exactly, degrades less. Sparse
graphs are essentially unaffected. `oscillation = :active` detects the stall on
the second iteration. See `notes/MP_and_OSD_decoder_findings.md`.
"""
function _fast_decode!(W::SoftDecisionWorkspace{Float64}, algo::Val, ::Val{:layered},
                       decimation_type::Val, osc_type::Val, max_iter::Int,
                       α::Float64, β::Float64, dec_thresh::Float64, dec_rounds::Int)
    isempty(W.layer_ptr) && throw(ArgumentError(
        "No layer partition in this workspace. Build it with " *
        "init_soft_workspace(...; schedule = :layered) or pass layer_ptr/layer_checks."))

    num_layers = length(W.layer_ptr) - 1
    @inbounds for iter in 1:max_iter
        for l in 1:num_layers
            for k in W.layer_ptr[l]:(W.layer_ptr[l + 1] - 1)
                c = W.layer_checks[k]
                lo = W.chk_ptr[c]
                hi = W.chk_ptr[c + 1] - 1

                # Withdraw this check's old contribution from the posteriors. What
                # is left is exactly the extrinsic belief, so it doubles as the
                # incoming message and no separate snapshot is needed. It already
                # includes every earlier layer of this iteration -- that feedback
                # is what layering buys.
                for e in lo:hi
                    v = W.edge_var[e]
                    W.is_decimated[v] || (W.total_llrs[v] -= W.C2V[e])
                    W.V2C[e] = W.total_llrs[v]
                end

                # Recompute the outgoing messages and fold them straight back in.
                _check_update!(W, algo, c, α, β)
                for e in lo:hi
                    v = W.edge_var[e]
                    W.is_decimated[v] || (W.total_llrs[v] += W.C2V[e])
                end
            end
        end
        _harden!(W)

        _apply_decimation!(decimation_type, W, iter, dec_thresh, dec_rounds)
        _syndrome_matches(W) && return true, iter

        if _check_oscillation(osc_type, W)
            if decimation_type !== Val(:none) && decimation_type !== Val(:manual)
                _apply_decimation!(decimation_type, W, iter, dec_thresh, dec_rounds; force = true)
                _wipe_history!(osc_type, W)
            else
                return false, -iter
            end
        else
            _update_history!(osc_type, W)
        end
    end
    return false, max_iter
end

# ==============================================================================
# DECIMATION, OSCILLATION
# ==============================================================================

# No decimation: the compiler erases the call.
@inline _apply_decimation!(::Val{:none}, W, iter, threshold, rounds; force::Bool = false) = nothing

"""Manual decimation is applied by the channel loader; nothing to do per iteration."""
@inline _apply_decimation!(::Val{:manual}, W, iter, threshold, rounds; force::Bool = false) =
    nothing

"""Auto decimation: pin any belief past the threshold."""
@inline function _apply_decimation!(::Val{:auto}, W::SoftDecisionWorkspace, iter::Int,
                                   threshold::Float64, rounds::Int; force::Bool = false)
    @inbounds for v in 1:W.num_var
        if !W.is_decimated[v]
            if W.total_llrs[v] > threshold
                W.is_decimated[v] = true
                W.total_llrs[v] = _LLR_MAX
                W.channel_llrs[v] = _LLR_MAX
            elseif W.total_llrs[v] < -threshold
                W.is_decimated[v] = true
                W.total_llrs[v] = -_LLR_MAX
                W.channel_llrs[v] = -_LLR_MAX
            end
        end
    end
    return nothing
end

"""Guided decimation: every `rounds` iterations, pin the most confident free bit."""
@inline function _apply_decimation!(::Val{:guided}, W::SoftDecisionWorkspace, iter::Int,
                                   threshold::Float64, rounds::Int; force::Bool = false)
    (force || iter % rounds == 0) || return nothing
    best_v = -1
    max_belief = -1.0
    @inbounds for v in 1:W.num_var
        if !W.is_decimated[v]
            belief = abs(W.total_llrs[v])
            if belief > max_belief
                max_belief = belief
                best_v = v
            end
        end
    end
    if best_v != -1
        @inbounds begin
            W.is_decimated[best_v] = true
            pinned = W.total_llrs[best_v] >= 0.0 ? _LLR_MAX : -_LLR_MAX
            W.total_llrs[best_v] = pinned
            W.channel_llrs[best_v] = pinned
        end
    end
    return nothing
end

@inline _check_oscillation(::Val{:none}, W) = false
@inline _update_history!(::Val{:none}, W) = nothing
@inline _wipe_history!(::Val{:none}, W) = nothing

@inline _check_oscillation(::Val{:active}, W) =
    (W.current_bits == W.prev_bits_1 || W.current_bits == W.prev_bits_2)

@inline function _update_history!(::Val{:active}, W)
    copyto!(W.prev_bits_2, W.prev_bits_1)
    copyto!(W.prev_bits_1, W.current_bits)
    return nothing
end

@inline function _wipe_history!(::Val{:active}, W)
    fill!(W.prev_bits_1, 0xFF)
    fill!(W.prev_bits_2, 0xFF)
    return nothing
end

# ==============================================================================
# ENTRY POINT
# ==============================================================================

"""
    decode!(W, LLR_in; kwargs...) -> (converged, iterations)

Decode one channel realization into `W`, and optionally copy the hard decisions
into `out`.

Only `(converged, iterations)` is returned, so that a decode moves no array
across a language boundary unless the caller asks for one (FIX-14). The decisions
also stay available as `W.current_bits` until the next decode overwrites them.
A negative `iterations` means the decoder stopped early on a detected
oscillation.

Keyword arguments:

  * `algorithm`: `:sum_product`, `:min_sum`, `:normalized_min_sum`,
    `:offset_min_sum` or `:min_sum_correction`.
  * `schedule`: `:flooding`, or `:layered`/`:serial` to sweep the partition the
    workspace was built with.
  * `max_iter`, `attenuation` (normalized min-sum), `offset` (offset min-sum).
  * `decimation`: `:none`, `:auto`, `:guided` or `:manual`.
  * `oscillation`: `:none`, or `:active` to detect a two-cycle in the hard
    decisions and either force a decimation step or give up.
  * `syndrome`, `erasures`, `decimated_bits`, `decimated_values`: forwarded to
    [`load_soft_channel!`](@ref).
  * `out`: a length-`num_var` integer buffer to receive the hard decisions.

On a small dense high-rate matrix decoded from a constant channel LLR vector, the
min-sum family under `:layered`/`:serial` can stall at an exact stationary point;
see the `:layered` engine above and `notes/MP_and_OSD_decoder_findings.md`.
"""
function decode!(W::SoftDecisionWorkspace{Float64}, LLR_in::AbstractVector{<:Real};
                 algorithm::Symbol = :offset_min_sum,
                 schedule::Symbol = :flooding,
                 decimation::Symbol = :none,
                 oscillation::Symbol = :none,
                 max_iter::Int = 100,
                 attenuation::Float64 = 0.75,
                 offset::Float64 = 0.5,
                 dec_thresh::Float64 = 10.0,
                 dec_rounds::Int = 10,
                 syndrome::AbstractVector{<:Integer} = _NO_INDICES,
                 erasures::AbstractVector{<:Integer} = _NO_INDICES,
                 decimated_bits::AbstractVector{<:Integer} = _NO_INDICES,
                 decimated_values::AbstractVector{<:Integer} = _NO_INDICES,
                 out::Union{Nothing, AbstractVector{<:Integer}} = nothing)
    load_soft_channel!(W, LLR_in; syndrome = syndrome, erasures = erasures,
                       decimated_bits = decimated_bits, decimated_values = decimated_values)

    target = schedule === :serial || schedule === :semiserial ? :layered : schedule
    converged, iters = _fast_decode!(W, Val(algorithm), Val(target), Val(decimation),
                                     Val(oscillation), max_iter, attenuation, offset,
                                     dec_thresh, dec_rounds)
    out === nothing || copyto!(out, W.current_bits)
    return converged, iters
end
