# Copyright (c) 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
# constructors
#############################

"""
    ConvolutionalCode(G::CTPolyMatrix, parity::Bool = false)

Return the convolutional code based on the matrix `G`.

# Notes
- Currently only accepting matrices with polynomial entries; any code with rational entries is equivalent to a code with polynomial entries by multiplying by the lcm of the denominators.
- Currently only accepting full-rank matrices.
"""
function ConvolutionalCode(G::CTPolyMatrix, parity::Bool = false)

    # TODO check base ring on G is appriorate
    D = gen(parent(G_new[1, 1]))

    G_new = deepcopy(G)
    G_new = _remove_empty(G_new, :rows)
    k = rank(G_new)
    k == nrows(G_new) || throw(ArgumentError("The input matrix must have `rank(G)` rows."))

    H = kernel(G_new, side = :right)
    rnk_H = rank(H)
    if ncols(H) == rnk_H
        H_tr = transpose(H)
    else
        # remove empty columns for flint objects https://github.com/oscar-system/Oscar.jl/issues/1062
        nr = nrows(H)
        H_tr = zero_matrix(base_ring(H), rnk_H, nr)
        for r = 1:nr
            for c = 1:rnk_H
                !iszero(H[r, c]) && (H_tr[c, r] = H[r, c];)
            end
        end
    end
    H = H_tr

    if parity
        temp = G_new
        G_new = H
        H = temp
    end

    n = ncols(G_new)
    vi = zeros(Int, k)
    row_degs = zeros(Int, k)
    for i = 1:k
        row_deg = [degree(G_new[i, c]) for c = 1:n]
        row_degs = sum(row_deg)
        vi[i] = maximum(row_deg)
    end
    m = maximum(vi)

    mnrs = minors(G_new, k)
    int_deg = maximum([degree(x) for x in mnrs])
    ext_deg = sum(row_degs)
    return ConvolutionalCode(
        base_ring(G_new),
        n,
        k,
        missing,
        D,
        m,
        vi,
        mnrs,
        int_deg,
        ext_deg,
        G,
        H,
        missing,
    )
end

#############################
# getter functions
#############################

"""
    field(C::AbstractConvolutionalCode)

Return the base ring of the generator matrix of `C`.
"""
field(C::AbstractConvolutionalCode) = C.F

"""
    length(C::AbstractConvolutionalCode)

Return the length of `C`.
"""
length(C::AbstractConvolutionalCode) = C.n

"""
    dimension(C::AbstractConvolutionalCode)

Return the dimension of the code.
"""
dimension(C::AbstractConvolutionalCode) = C.k

"""
    rate(C::AbstractConvolutionalCode)

Return the rate of `C`.
"""
rate(C::AbstractConvolutionalCode) = C.k / C.n

"""
    parity_check_matrix(C::AbstractConvolutionalCode)

Return the parity check matrix of `C`.
"""
parity_check_matrix(C::AbstractConvolutionalCode) = C.H

"""
    delay_operator(C::AbstractConvolutionalCode)

Return the delay operator of `C`.
"""
delay_operator(C::AbstractConvolutionalCode) = C.D

"""
    constraint_lengths(C::AbstractConvolutionalCode)

Return the constraint lengths of `C`.    
"""
constraint_lengths(C::AbstractConvolutionalCode) = C.vi

"""
    overall_constraint_lengths(C::AbstractConvolutionalCode)

Return the overall constraint length of `C`.    
"""
overall_constraint_length(C::AbstractConvolutionalCode) = sum(C.vi)

"""
    memory(C::AbstractConvolutionalCode)

Return the memory of `C`.
"""
memory(C::AbstractConvolutionalCode) = C.m

"""
    internal_degree(C::AbstractConvolutionalCode)

Return the internal degree of `C`.
"""
internal_degree(C::AbstractConvolutionalCode) = C.int_deg

"""
    external_degree(C::AbstractConvolutionalCode)

Return the external degree of `C`.
"""
external_degree(C::AbstractConvolutionalCode) = C.ext_deg

#############################
# setter functions
#############################

#############################
# general functions
#############################

"""
    is_realizable(C::AbstractConvolutionalCode)

Return `true` if the generator matrix of `C` is realizable; otherwise, `false`.

# Notes
- A generator matrix is realizable if it is a polynomial matrix or all the denominators of the rational functions are equal to one. Since only polynomial matrices are accepted at the moment, this function always returns true.
"""
is_realizable(C::AbstractConvolutionalCode) = true

"""
    is_polynomial(C::AbstractConvolutionalCode)

Return `true` if the generator matrix of `C` is polynomial; otherwise, `false`.

# Notes
- This would check all the denominators of the rational functions are equal to one; however, since only polynomial matrices are accepted at the moment, this function always returns true.
"""
is_polynomial(C::AbstractConvolutionalCode) = true

"""
    is_basic(C::AbstractConvolutionalCode)

Return `true` if the generator matrix of `C` is basic; otherwise, `false`.
"""
function is_basic(C::AbstractConvolutionalCode)
    SNF, _, _ = snf_with_transform(C.G)
    # diag elements of SNF are invariant factors
    # all invariant factors (SNF) are 1
    return all(is_one, [SNF[i, i] for i = 1:k])
end

# can always make minimal with row ops
"""
    is_minimal(C::AbstractConvolutionalCode)

Return `true` if the generator matrix of `C` is minimal; otherwise, `false`.
"""
function is_minimal(C::AbstractConvolutionalCode)
    is_basic(C) || return false
    G_h = zeros(UInt8, C.k, C.n)
    for r = 1:C.k
        for c in C.n
            degree(C.G[r, c]) == C.vi[r] && (G_h[r, c] = 1)
        end
    end
    return k == rank(G_h)
end

"""
    is_reduced(C::AbstractConvolutionalCode)

Return `true` if the generator matrix of `C` is reduced; otherwise, `false`.
"""
is_reduced(C::AbstractConvolutionalCode) = internal_degree(C) == external_degree(C)

"""
    is_canonical(C::AbstractConvolutionalCode)

Return `true` if the generator matrix of `C` is canonical; otherwise, `false`.
"""
is_canonical(C::AbstractConvolutionalCode) = is_basic(C) && is_reduced(C)

"""
    is_catastrophic(C::AbstractConvolutionalCode)

Return `true` if the generator matrix of `C` is catastrophic; otherwise, `false`.
"""
is_catastrophic(C::AbstractConvolutionalCode) = !is_one(sum(ZZ.(coefficients(gcd(C.mnrs)))))

"""
    generator_matrix(C::AbstractConvolutionalCode, systematic::Bool = false)

Return the generator matrix of `C` if `systematic` is `false` and a systematic generator matrix
otherwise.
"""
generator_matrix(C::AbstractConvolutionalCode, systematic::Bool = false) =
    systematic ? (return hnf(fraction_field(base_ring(C.G)).(C.G));) : (return C.G;)

function terminated_generator_matrix(C::AbstractConvolutionalCode, L::Int)
    L ≥ 1 || throw(DomainError(L, "The number of rows must be at least one."))

    G_mats = [zero_matrix(C.F, C.k, C.n) for _ = 1:C.m]
    for r = 1:C.k
        for c = 1:C.n
            for i = 1:m
                G_mats[m][r, c] = coeff(C.G[r, c], i - 1)
            end
        end
    end

    # this matrix has L shifts and truncates at row L
    G_L = zero_matrix(C.F, L * C.k, (m + L) * C.n)
    for r = 1:L
        for c = 1:m
            G_L[((r-1)*C.k+1):(r*C.k), ((c-1+r-1)*C.n+1):((c+r-1)*C.n)] .= G_mats[c]
        end
    end
    return G_L
end

function truncated_generator_matrix(C::AbstractConvolutionalCode, L::Int)
    L ≥ 1 || throw(DomainError(L, "The number of columns must be at least one."))

    G_mats = [zero_matrix(C.F, C.k, C.n) for _ = 1:C.m]
    for r = 1:C.k
        for c = 1:C.n
            for i = 1:m
                G_mats[m][r, c] = coeff(C.G[r, c], i - 1)
            end
        end
    end

    # this matrix truncates at column L
    G_L = zeros(C.F, L * C.k, L * C.n)
    for offset = 0:min(C.m-1, L)
        for i = 1:(L-offset)
            rows = (1:C.k) .+ ((i - 1) * C.k)
            cols = (1:C.n) .+ ((offset + i - 1) * C.n)
            G_L[rows, cols] .= G_mats[offset+1]
        end
    end
    return G_L
end

function tail_biting_generator_matrix(C::AbstractConvolutionalCode, L::Int)
    L ≥ 1 || throw(DomainError(L, "The number of columns must be at least one."))

    G_mats = [zero_matrix(C.F, C.k, C.n) for _ = 1:C.m]
    for r = 1:C.k
        for c = 1:C.n
            for i = 1:m
                G_mats[m][r, c] = coeff(C.G[r, c], i - 1)
            end
        end
    end

    # this matrix has L shifts and truncates at row L
    G_L = zero_matrix(C.F, L * C.k, L * C.n)
    for r = 1:L
        for c = 1:m
            if c - 1 + r - 1 ≤ L
                if c + r - 2 > L
                    cols = (((c-1+r-1)%L)*C.n+1):(((c+r-1)%L)*C.n)
                else
                    cols = ((c-1+r-1)*C.n+1):((c+r-1)*C.n)
                end
                G_L[((r-1)*C.k+1):(r*C.k), cols] .= G_mats[c]
            end
        end
    end
end

"""
    dual(C::AbstractConvolutionalCode)

Return the dual code of `C`.
"""
function dual(C::AbstractConvolutionalCode)
    # TODO do better in the future, avoid the recomputation
    return ConvolutionalCode(C.H)
end

function show(io::IO, C::AbstractConvolutionalCode)
    print(io, "($(C.n), $(C.k)")
    !ismissing(C.d) && print(io, ", $(C.d)")
    print(io, ")_$(order(C.F)) ")
    println(io, "convolutional code")

    if get(io, :compact, true)
        println(io, "Memory: $(C.m)")
        println(io, "Constraint length: $(overall_constraint_length(C))")
        if C.n ≤ 10
            println(io, "Generator matrix: $(C.k) × $(C.n)")
            println(
                io,
                "\t" * replace(
                    replace(
                        replace(repr(MIME("text/plain"), C.G), r"\n" => "\n\t"),
                        r"\[" => "",
                    ),
                    r"\]" => "",
                ),
            )
        end
        # if !ismissing(C.weight_enum)
        #     println(io, "\nComplete weight enumerator:")
        #     print(io, "\t", polynomial(C.weight_enum))
        # end
    end
end

"""
    encode(C::AbstractConvolutionalCode, v::CTPolyMatrix)

Return the encoding of `v` into `C`.
"""
function encode(C::AbstractConvolutionalCode, v::CTPolyMatrix)
    (size(v) != (1, C.k) && size(v) != (C.k, 1)) && throw(
        ArgumentError(
            "Vector has incorrect dimension; expected length $(C.k), received: $(size(v)).",
        ),
    )
    parent(v) == parent(C.G) ||
        throw(ArgumentError("Vector must have the same parent as the generator matrix."))
    nrows(v) ≠ 1 || return v * C.G
    return transpose(v) * C.G
end

# free distance
# canonical_generator_matrix
# minimum_distance
#     column_distance
#     distance_profile
#     row_distance
#     free_distance
#     also extended versions
# encode
#     interleave::Bool
#     scalar and polynomial versions
# syndrome
