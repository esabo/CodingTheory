# Copyright (c) 2021 - 2026 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
          # Misc
#############################

"""
$(TYPEDSIGNATURES)

Return the `[n, 0, 0]` zero code over `F`.
"""
function ZeroCode(F::CTFieldTypes, n::Integer)
    n > 0 || throw(ArgumentError("Code length must be positive (received n = $n)"))
    
    # A zero code has dimension 0, so G is 0 x n
    G = zero_matrix(F, 0, n)
    H = identity_matrix(F, n)
    
    # We can explicitly inject the weight distribution into the cache
    counts = Dict{Int, BigInt}(0 => BigInt(1))
    hwe = HammingWeightEnumerator(n, counts)
    
    cache = Dict{Symbol, Any}(:G => G, :H => H, :weight_enum => hwe)
    return LinearCode(F, n, 0, 0, 0, 0, cache)
end

"""
$(TYPEDSIGNATURES)

Return the `[n, 0, 0]` zero code over `GF(q)`.
"""
function ZeroCode(q::Integer, n::Integer)
    F = if is_prime(q) 
        Oscar.Nemo.Native.GF(q)
    else
        factors = Nemo.factor(q)
        length(factors) == 1 || throw(DomainError("There is no finite field of order $q"))
        p, t = first(factors)
        GF(p, t, :α)
    end
    return ZeroCode(F, n)
end

"""
$(TYPEDSIGNATURES)

Return the `[n, 0, 0]` binary zero code.
"""
ZeroCode(n::Integer) = ZeroCode(Oscar.Nemo.Native.GF(2), n)

"""
$(TYPEDSIGNATURES)

Return the `[n, n, 1]` identity code over `F`.
"""
IdentityCode(F::CTFieldTypes, n::Integer) = dual(ZeroCode(F, n))

"""
$(TYPEDSIGNATURES)

Return the `[n, n, 1]` identity code over `GF(q)`.
"""
IdentityCode(q::Integer, n::Integer) = dual(ZeroCode(q, n))

"""
$(TYPEDSIGNATURES)

Return the `[n, n, 1]` binary identity code.
"""
IdentityCode(n::Integer) = dual(ZeroCode(n))

"""
$(TYPEDSIGNATURES)

Return the `[n, 1, n]` repetition code over `GF(q)`.
"""
function RepetitionCode(q::Int, n::Int)
    F = if is_prime(q)
        Oscar.Nemo.Native.GF(q)
    else
        factors = Nemo.factor(q)
        length(factors) == 1 || throw(DomainError("There is no finite field of order $q"))
        (p, t) = first(factors) # BUG FIX: Safer destructuring
        GF(p, t, :α)
    end
    
    G = matrix(F, ones(Int, 1, n))
    H = hcat(matrix(F, ones(Int, n - 1, 1)), identity_matrix(F, n - 1))
    
    # Exact weight distribution is trivial for a repetition code
    counts = Dict{Int, BigInt}(0 => BigInt(1), n => BigInt(q - 1))
    hwe = HammingWeightEnumerator(n, counts)
    
    cache = Dict{Symbol, Any}(:G => G, :H => H, :weight_enum => hwe)
    return LinearCode(F, n, 1, n, n, n, cache)
end

"""
$(TYPEDSIGNATURES)

Return the `[n, n-1, 2]` single parity check code over `GF(q)`.
"""
function SingleParityCheckCode(q::Int, n::Int)
    iseven(q) && (return dual(RepetitionCode(q, n)))
    
    F = if is_prime(q)
        Oscar.Nemo.Native.GF(q)
    else
        factors = Nemo.factor(q)
        length(factors) == 1 || throw(DomainError("There is no finite field of order $q"))
        (p, t) = first(factors)
        GF(p, t, :α)
    end
    
    G = hcat(matrix(F, .-ones(Int, n - 1, 1)), identity_matrix(F, n - 1))
    
    cache = Dict{Symbol, Any}(:G => G)
    return LinearCode(F, n, n - 1, 2, 2, 2, cache)
end
SPCCode(q::Int, n::Int) = SingleParityCheckCode(q, n)

"""
$(TYPEDSIGNATURES)

Return the `[6, 3, 4]` hexacode over `GF(4)`.
"""
function Hexacode()
    # TODO: is Hexacode a Hamming code?
    F = GF(2, 2, :ω)
    ω = gen(F)
    G = matrix(F, [1 0 0 1 ω ω; 0 1 0 ω 1 ω; 0 0 1 ω ω 1])
    H = matrix(F, [1 ω ω 1 0 0; ω 1 ω 0 1 0; ω ω 1 0 0 1])
    
    # Omit the eager weight enumerator; let the lazy getter do the work later if asked
    cache = Dict{Symbol, Any}(:G => G, :H => H)
    return LinearCode(F, 6, 3, 4, 4, 4, cache)
end

#############################
         # Hamming
#############################
# unclear if this should be promoted to its own type so r can be extracted

# ==============================================================================
# HAMMING CODES
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the `[(q^r - 1)/(q - 1), (q^r - 1)/(q - 1) - r, 3]` Hamming code over `GF(q)`.
"""
function HammingCode(q::Int, r::Int)
    2 ≤ r || throw(DomainError(r, "Hamming codes require r ≥ 2; received r = $r."))
    r < 64 || throw(DomainError(r, "This Hamming code requires the implementation of BigInts. Change if necessary."))

    factors = Nemo.factor(q)
    length(factors) == 1 || throw(ArgumentError("There is no finite field of order $q."))

    F = if is_prime(q) 
        Oscar.Nemo.Native.GF(q)
    else
        p, t = first(factors)
        GF(p, t, :α)
    end

    n = div(q^r - 1, q - 1)
    k = n - r
    d = 3

    H = zero_matrix(F, r, n)
    
    # A valid parity-check matrix for a general Hamming code over GF(q) consists 
    # of exactly one representative from each 1-dimensional subspace.
    # We select all vectors where the first non-zero element is exactly `one(F)`.
    col_idx = 1
    for iter in Nemo.AbstractAlgebra.ProductIterator([F for _ in 1:r], inplace = true)
        for i in 1:r
            if !iszero(iter[i])
                if iter[i] == one(F)
                    for j in 1:r
                        H[j, col_idx] = iter[j]
                    end
                    col_idx += 1
                end
                break # Stop checking after finding the first non-zero element
            end
        end
    end

    cache = Dict{Symbol, Any}(:H => H)
    return HammingCode(F, n, k, d, d, d, r, cache)
end

# ==============================================================================
# EXTENDED HAMMING CODES
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the `[2^r, 2^r - 1 - r, 4]` extended binary Hamming code.

# Notes
* Extended Hamming codes are formed by adding an overall parity-check bit to a 
  binary Hamming code, creating a SEC-DED (Single Error Correction, Double Error Detection) code.
"""
function ExtendedHammingCode(r::Int)
    2 ≤ r || throw(DomainError(r, "Extended Hamming codes require r ≥ 2; received r = $r."))
    r < 64 || throw(DomainError(r, "This Extended Hamming code requires the implementation of BigInts. Change if necessary."))

    F = Oscar.Nemo.Native.GF(2)
    n = 2^r
    k = 2^r - 1 - r
    d = 4

    # The parity check matrix has dimension (r + 1) x 2^r
    H = zero_matrix(F, r + 1, n)
    
    # Fill the first r rows with the standard Hamming columns
    for i in 1:(n - 1)
        bin_vals = digits(i, base=2, pad=r)
        for j in 1:r
            H[j, i] = bin_vals[j]
        end
    end
    # The last column's first r rows remain 0
    
    # Fill the bottom row with 1s for the overall parity check
    for i in 1:n
        H[r + 1, i] = 1
    end

    cache = Dict{Symbol, Any}(:H => H)
    return LinearCode(F, n, k, d, d, d, cache)
end

# ==============================================================================
# TETRA CODE
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the `[4, 2, 3]` tetra code over `GF(3)`.

# Notes
* This is equivalent to the `Hamming(3, 2, 3)` code, but the construction here is
  based on the commonly presented generator and parity-check matrices.
"""
function TetraCode()
    F = Oscar.Nemo.Native.GF(3)
    G = matrix(F, [1 0 1 1; 0 1 1 -1])
    H = matrix(F, [-1 -1 1 0; -1 1 0 1])
    
    n, k, d = 4, 2, 3
    
    # The Tetra code has exactly 1 codeword of weight 0, and 8 codewords of weight 3
    counts = Dict{Int, BigInt}(0 => BigInt(1), 3 => BigInt(8))
    hwe = HammingWeightEnumerator(n, counts)
    
    cache = Dict{Symbol, Any}(:G => G, :H => H, :weight_enum => hwe)
    return LinearCode(F, n, k, d, d, d, cache)
end

# ==============================================================================
# SIMPLEX CODES
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the `[(q^r - 1)/(q - 1), r]` simplex code over `GF(q)`.

# Notes
* Generator matrices for the binary codes are constructed using the standard recursive definition.
* The higher fields return `dual(HammingCode(q, r))`.
* This is currently only implemented for binary codes.
"""
function SimplexCode(q::Int, r::Int)
    2 ≤ r || throw(DomainError(r, "Simplex codes require 2 ≤ r; received r = $r."))
    r < 64 || throw(DomainError(r, "The weight enumerator for the simplex codes for r > 64 require BigInts. Implement if necessary."))
    
    factors = Nemo.factor(q)
    length(factors) == 1 || throw(ArgumentError("There is no finite field of order $q."))

    q > 2 && return dual(HammingCode(q, r))

    # Binary simplex codes
    F = Oscar.Nemo.Native.GF(2)
    G2 = matrix(F, [0 1 1; 1 0 1])
    
    if r == 2
        Grm1 = G2
    else
        Grm1 = G2
        for i in 3:r
            zs1 = matrix(F, nrows(Grm1), 1, zeros(Int, nrows(Grm1), 1))
            bot = hcat(Grm1, zs1, Grm1)
            zs2 = matrix(F, 1, ncols(Grm1), zeros(Int, 1, ncols(Grm1)))
            os = matrix(F, 1, ncols(Grm1) + 1, ones(Int, 1, ncols(Grm1) + 1))
            top = hcat(zs2, os)
            Grm1 = vcat(top, bot)
        end
    end
    
    n = 2^r - 1
    k = r
    d = 2^(r - 1)
    
    # All nonzero codewords have weight 2^{r - 1}. There are 2^r - 1 of them.
    counts = Dict{Int, BigInt}(0 => BigInt(1), d => BigInt(2^r - 1))
    hwe = HammingWeightEnumerator(n, counts)
    
    cache = Dict{Symbol, Any}(:G => Grm1, :weight_enum => hwe)
    return SimplexCode(F, n, k, d, d, d, r, cache)
end

#############################
          # Golay
#############################

"""
$(TYPEDSIGNATURES)

Return the `[24, 12, 8]` extended binary Golay code if `p == 2` or the `[12, 6, 6]`
extended ternary Golay code if `p == 3`.
"""
function ExtendedGolayCode(p::Int)
    if p == 2
        F = Oscar.Nemo.Native.GF(2)
        A = matrix(F, [0 1 1 1 1 1 1 1 1 1 1 1;
                       1 1 1 0 1 1 1 0 0 0 1 0;
                       1 1 0 1 1 1 0 0 0 1 0 1;
                       1 0 1 1 1 0 0 0 1 0 1 1;
                       1 1 1 1 0 0 0 1 0 1 1 0;
                       1 1 1 0 0 0 1 0 1 1 0 1;
                       1 1 0 0 0 1 0 1 1 0 1 1;
                       1 0 0 0 1 0 1 1 0 1 1 1;
                       1 0 0 1 0 1 1 0 1 1 1 0;
                       1 0 1 0 1 1 0 1 1 1 0 0;
                       1 1 0 1 1 0 1 1 1 0 0 0;
                       1 0 1 1 0 1 1 1 0 0 0 1])
        G = hcat(identity_matrix(F, 12), A)
        H = hcat(-transpose(A), identity_matrix(F, 12))
        
        # Exact weight enumerator dictionary injection bypasses expensive math
        counts = Dict{Int, BigInt}(0 => 1, 8 => 759, 12 => 2576, 16 => 759, 24 => 1)
        hwe = HammingWeightEnumerator(24, counts)
        
        cache = Dict{Symbol, Any}(:G => G, :H => H, :weight_enum => hwe)
        return LinearCode(F, 24, 12, 8, 8, 8, cache)
        
    elseif p == 3
        F = Oscar.Nemo.Native.GF(3)
        A = matrix(F, [0  1  1  1  1  1;
                       1  0  1 -1 -1  1;
                       1  1  0  1 -1 -1;
                       1 -1  1  0  1 -1;
                       1 -1 -1  1  0  1;
                       1  1 -1 -1  1  0])
        G = hcat(identity_matrix(F, 6), A)
        H = hcat(-transpose(A), identity_matrix(F, 6))
        
        # Exact weight enumerator dictionary injection
        counts = Dict{Int, BigInt}(0 => 1, 6 => 264, 9 => 440, 12 => 24)
        hwe = HammingWeightEnumerator(12, counts)
        
        cache = Dict{Symbol, Any}(:G => G, :H => H, :weight_enum => hwe)
        return LinearCode(F, 12, 6, 6, 6, 6, cache)
        
    else
        throw(ArgumentError("Golay code not implemented for p = $p."))
    end
end

"""
$(TYPEDSIGNATURES)

Return the `[23, 12, 7]` binary Golay code if `p == 2` or the `[11, 6, 5]`
ternary Golay code if `p == 3`.
"""
function GolayCode(p::Int)
    C = puncture(ExtendedGolayCode(p), [1])
    d_val = p == 2 ? 7 : 5
    set_minimum_distance!(C, d_val) # Safely update d, l_bound, and u_bound
    return C
end

# ==============================================================================
# HADAMARD CODES
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the `[2^m, m, 2^{m - 1}]` binary Hadamard code.

# Notes
* Hadamard codes are generally constructed using Hadamard matrices `H_n`. When `n`
  is a power of two, the codes are linear. Constructing `H_{2^m}` then mapping `+/- 1`
  to `{0, 1}` gives a generator matrix which, up to permutation, is equivalent to
  using all `2^m` binary strings as column vectors. The construction here uses this
  latter definition.
* Note that some engineering fields define the Hadamard code to be the `[2^m, m + 1, 2^{m - 1}]` 
  augmented Hadamard code (equivalent to the first-order Reed-Muller code `RM(1, m)`).
"""
function HadamardCode(m::Int)
    m ≥ 1 || throw(DomainError(m, "Hadamard codes require m ≥ 1."))
    m < 64 || throw(DomainError(m, "This Hadamard code requires the implementation of BigInts. Change if necessary."))

    F = Oscar.Nemo.Native.GF(2)
    n = 2^m
    k = m
    d = 2^(m - 1)

    # Pre-allocate to avoid massive memory allocations on hcat
    G = zero_matrix(F, k, n)
    
    for i in 0:(n - 1)
        # Using digits natively populates the binary representations.
        # (Whether it's reversed or not doesn't matter, as it generates the complete set of columns).
        bin_vals = digits(i, base=2, pad=m)
        for j in 1:m
            G[j, i + 1] = bin_vals[j]
        end
    end

    # Each non-zero codeword has a Hamming weight of exactly 2^{m-1}
    counts = Dict{Int, BigInt}(0 => BigInt(1), d => BigInt(n - 1))
    hwe = HammingWeightEnumerator(n, counts)
    
    cache = Dict{Symbol, Any}(:G => G, :weight_enum => hwe)
    return LinearCode(F, n, k, d, d, d, cache)
end
WalshHadamardCode(m::Int) = HadamardCode(m)
WalshCode(m::Int) = HadamardCode(m)

# ==============================================================================
# BEST KNOWN LINEAR CODES (GAP Interface)
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the best known linear code of length `n` and dimension `k` over `GF(2)`
using the internal GAP database.
"""
function best_known_linear_code(n::Int, k::Int)
    C_GAP = GAP.Globals.BestKnownLinearCode(n, k, GAP.Globals.GF(2))
    
    G_GAP = GAP.Globals.GeneratorMat(C_GAP)
    dims_G = GAP.Globals.DimensionsMat(G_GAP)
    g = matrix(GF(2), [GAP.Globals.Int(G_GAP[i, j]) for i in 1:dims_G[1], j in 1:dims_G[2]])
    
    H_GAP = GAP.Globals.CheckMat(C_GAP)
    dims_H = GAP.Globals.DimensionsMat(H_GAP)
    h = matrix(GF(2), [GAP.Globals.Int(H_GAP[i, j]) for i in 1:dims_H[1], j in 1:dims_H[2]])
    
    # Uses our safe lazy constructor that supports G and H injection natively
    return LinearCode(g, h)
end

# ==============================================================================
# MACDONALD CODES
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the `[(q^k - q^u)/(q - 1), k, q^{k - 1} - q^{u - 1}]` MacDonald code over `GF(q)`.

# Notes
* MacDonald codes are the classic examples of two-weight codes.
* They are constructed by taking a Simplex code of dimension `k` and removing the columns 
  that correspond to a Simplex subcode of dimension `u`.
"""
function MacDonaldCode(q::Int, k::Int, u::Int)
    k > u >= 1 || throw(DomainError((k, u), "MacDonald codes require k > u ≥ 1."))
    
    factors = Nemo.factor(q)
    length(factors) == 1 || throw(ArgumentError("There is no finite field of order $q."))

    F = if is_prime(q) 
        Oscar.Nemo.Native.GF(q)
    else
        p, t = first(factors)
        GF(p, t, :α)
    end

    n = div(q^k - q^u, q - 1)
    d = q^(k - 1) - q^(u - 1)
    
    G = zero_matrix(F, k, n)
    
    col_idx = 1
    for iter in Nemo.AbstractAlgebra.ProductIterator([F for _ in 1:k], inplace = true)
        # Find the first non-zero element to ensure projective uniqueness (Simplex condition)
        first_nz_idx = 0
        for i in 1:k
            if !iszero(iter[i])
                first_nz_idx = i
                break
            end
        end
        
        if first_nz_idx > 0 && iter[first_nz_idx] == one(F)
            # Puncture condition: Keep the column ONLY IF at least one of 
            # the last (k - u) elements is non-zero.
            keep_column = false
            for i in (u + 1):k
                if !iszero(iter[i])
                    keep_column = true
                    break
                end
            end
            
            if keep_column
                for i in 1:k
                    G[i, col_idx] = iter[i]
                end
                col_idx += 1
            end
        end
    end

    cache = Dict{Symbol, Any}(:G => G)
    return MacDonaldCode(F, n, k, d, d, d, u, cache)
end

# ==============================================================================
# LEXICODES (Lexicographic Codes)
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the binary Lexicode of length `n` and minimum distance `d`.

# Notes
* Lexicodes are generated via a greedy algorithm. The algorithm iterates through all 
  `2^n` binary vectors in lexicographic order, adding a vector to the code if its 
  Hamming distance to all currently selected vectors is at least `d`.
* By the Conway-Sloane theorem, this greedy construction over GF(2) naturally 
  produces a linear code.
"""
function Lexicode(n::Int, d::Int)
    d <= n || throw(DomainError(d, "Minimum distance d must be ≤ n."))
    d > 0  || throw(DomainError(d, "Minimum distance d must be positive."))
    n <= 24 || @warn "Lexicode generation for n > 24 may take a significant amount of time."

    # Hardware-optimized greedy span generation using integers
    basis_int = Int[]
    codewords_int = [0]
    
    for i in 1:(2^n - 1)
        # Check distance to the span of the current basis
        min_wt = n + 1
        for cw in codewords_int
            # Hamming distance is the number of 1s in the XOR difference
            w = count_ones(i ⊻ cw)
            if w < d
                min_wt = w
                break
            end
        end
        
        if min_wt >= d
            push!(basis_int, i)
            # Add the new coset to the tracked codewords
            new_cws = [i ⊻ cw for cw in codewords_int]
            append!(codewords_int, new_cws)
        end
    end

    k = length(basis_int)
    F = Oscar.Nemo.Native.GF(2)
    
    if k == 0
        return ZeroCode(n)
    end

    # Translate the integer basis back to the formal generator matrix
    G = zero_matrix(F, k, n)
    for i in 1:k
        bin_vals = digits(basis_int[i], base=2, pad=n)
        # Lexicodes are traditionally read left-to-right (MSB to LSB), 
        # so we reverse the digits to match the standard lexicographic basis.
        reverse!(bin_vals) 
        for j in 1:n
            G[i, j] = bin_vals[j]
        end
    end

    cache = Dict{Symbol, Any}(:G => G)
    return LinearCode(F, n, k, d, d, d, cache)
end
