# Copyright (c) 2022 - 2026 Eric Sabo, Michael Vasmer
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
    concatenate(S_outer::AbstractStabilizerCode, S_inner::AbstractStabilizerCode)

Return a lazy `QuantumConcatenatedCode`. Computes parameters instantly 
without generating the massive concatenated stabilizer matrices.
"""
function concatenate(S_outer::AbstractStabilizerCode, S_inner::AbstractStabilizerCode)
    S_inner.k == 1 || throw(ArgumentError("The inner code must encode exactly 1 logical qubit (k=1)."))
    
    F = S_outer.F # Assuming F is accessible; if it's in a cache, use S_outer.cache[:F]
    
    # 1. Physical and Logical Qubits
    n_new = S_outer.n * S_inner.n
    k_new = S_outer.k 
    
    # 2. Distance Bounds (Multiplicative)
    l_bound = S_outer.l_bound * S_inner.l_bound
    u_bound = S_outer.u_bound * S_inner.u_bound
    
    d_exact = missing
    if !ismissing(S_outer.d) && !ismissing(S_inner.d)
        d_exact = S_outer.d * S_inner.d
    end
    
    cache = Dict{Symbol, Any}(:F => F)
    
    return QuantumConcatenatedCode(
        S_outer, S_inner, 
        n_new, k_new, 
        d_exact, l_bound, u_bound, 
        cache
    )
end

function stabilizers(S::QuantumConcatenatedCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    S_outer = S.outer_code
    S_inner = S.inner_code
    F = S.cache[:F]
    
    n_out = S_outer.n
    n_in = S_inner.n
    n_new = S.n
    
    stabs_out = stabilizers(S_outer)
    stabs_in = stabilizers(S_inner)
    
    L_in = logicals_matrix(S_inner)
    L_X_in = L_in[1:1, :]
    L_Z_in = L_in[2:2, :]
    
    num_stabs_in = nrows(stabs_in)
    num_stabs_out = nrows(stabs_out)
    
    is_sparse = stabs_out isa SparseMatrixCSC || stabs_in isa SparseMatrixCSC
    
    if is_sparse
        stabs_new_inner = spzeros(F, n_out * num_stabs_in, 2 * n_new)
        stabs_new_outer = spzeros(F, num_stabs_out, 2 * n_new)
    else
        stabs_new_inner = zero_matrix(F, n_out * num_stabs_in, 2 * n_new)
        stabs_new_outer = zero_matrix(F, num_stabs_out, 2 * n_new)
    end
    
    # 1. Inner Code Stabilizers (Applied to each block)
    curr_row = 1
    for i in 1:n_out
        offset_X = (i - 1) * n_in
        offset_Z = n_new + (i - 1) * n_in
        for r in 1:num_stabs_in
            for c in 1:n_in
                stabs_new_inner[curr_row, offset_X + c] = stabs_in[r, c]
                stabs_new_inner[curr_row, offset_Z + c] = stabs_in[r, n_in + c]
            end
            curr_row += 1
        end
    end
    
    # 2. Outer Code Stabilizers (Mapped through inner logicals)
    for r in 1:num_stabs_out
        for i in 1:n_out
            x_val = stabs_out[r, i]
            z_val = stabs_out[r, n_out + i]
            
            offset_X = (i - 1) * n_in
            offset_Z = n_new + (i - 1) * n_in
            for c in 1:n_in
                # Symplectic dot product applying the outer check to the inner logical block
                stabs_new_outer[r, offset_X + c] += x_val * L_X_in[1, c] + z_val * L_Z_in[1, c]
                stabs_new_outer[r, offset_Z + c] += x_val * L_X_in[1, n_in + c] + z_val * L_Z_in[1, n_in + c]
            end
        end
    end
    
    stabs_final = vcat(stabs_new_inner, stabs_new_outer)
    
    S.cache[:stabilizers] = stabs_final
    return stabs_final
end

"""
    HypergraphProductCode(C1::AbstractLinearCode, C2::AbstractLinearCode)

Return a lazy `HypergraphProductCode`. Computes parameters `n`, `k`, and bounds 
instantly without generating the quantum parity check matrices.
"""
function HypergraphProductCode(C1::AbstractLinearCode, C2::AbstractLinearCode)
    C1.F == C2.F || throw(ArgumentError("Codes must be over the same field."))
    
    # Parity check matrices for classical codes are generally already cached or cheap.
    # We only need them to grab the exact dimensions m1 and m2.
    H1 = parity_check_matrix(C1)
    H2 = parity_check_matrix(C2)
    m1, n1 = size(H1)
    m2, n2 = size(H2)
    
    # 1. Calculate Physical Qubits (n)
    n_new = n1 * n2 + m1 * m2
    
    # 2. Calculate Logical Qubits (k)
    k1_dual = m1 - (n1 - C1.k)
    k2_dual = m2 - (n2 - C2.k)
    k_new = C1.k * C2.k + k1_dual * k2_dual
    
    # 3. Handle Distance Bounds
    # The $X$-distance and $Z$-distance depend on both the primal distances ($d_1, d_2$) and the dual distances ($d_1^T, d_2^T$):$d_X = \min(d_1^T, d_2)$$d_Z = \min(d_1, d_2^T)$Because the overall distance is $d = \min(d_X, d_Z)$, the value $\min(d_1, d_2)$ is strictly an upper bound.
    d_u_bound = min(C1.u_bound, C2.u_bound)
    if !ismissing(C1.d) && !ismissing(C2.d)
        d_u_bound = min(C1.d, C2.d)
    end
    
    cache = Dict{Symbol, Any}(:F => C1.F)
    
    return HypergraphProductCode(
        C1, C2, 
        n_new, k_new, 
        missing, 1, d_u_bound, 
        cache
    )
end
HypergraphProductCode(C::AbstractLinearCode) = HypergraphProductCode(C, C)

function X_stabilizers(S::HypergraphProductCode)
    haskey(S.cache, :H_X) && return S.cache[:H_X]
    
    F = S.cache[:F]
    H1 = parity_check_matrix(S.C1)
    H2 = parity_check_matrix(S.C2)
    m1, n1 = size(H1)
    m2, n2 = size(H2)
    
    is_sparse = H1 isa SparseMatrixCSC || H2 isa SparseMatrixCSC
    F_one = F(1)
    
    if is_sparse
        H1_sp = H1 isa SparseMatrixCSC ? H1 : sparse(H1)
        H2_sp = H2 isa SparseMatrixCSC ? H2 : sparse(H2)
        I_n2 = spdiagm(0 => fill(F_one, n2))
        I_m1 = spdiagm(0 => fill(F_one, m1))
        
        H_X = hcat(kron(H1_sp, I_n2), kron(I_m1, transpose(H2_sp)))
    else
        I_n2 = identity_matrix(F, n2)
        I_m1 = identity_matrix(F, m1)
        
        H_X = hcat(kronecker_product(H1, I_n2), kronecker_product(I_m1, transpose(H2)))
    end
    
    S.cache[:H_X] = H_X
    return H_X
end

function Z_stabilizers(S::HypergraphProductCode)
    haskey(S.cache, :H_Z) && return S.cache[:H_Z]
    
    F = S.cache[:F]
    H1 = parity_check_matrix(S.C1)
    H2 = parity_check_matrix(S.C2)
    m1, n1 = size(H1)
    m2, n2 = size(H2)
    
    is_sparse = H1 isa SparseMatrixCSC || H2 isa SparseMatrixCSC
    F_one = F(1)
    
    if is_sparse
        H1_sp = H1 isa SparseMatrixCSC ? H1 : sparse(H1)
        H2_sp = H2 isa SparseMatrixCSC ? H2 : sparse(H2)
        I_n1 = spdiagm(0 => fill(F_one, n1))
        I_m2 = spdiagm(0 => fill(F_one, m2))
        
        if Int(characteristic(F)) == 2
            H_Z = hcat(kron(I_n1, H2_sp), kron(transpose(H1_sp), I_m2))
        else
            H_Z = hcat(kron(I_n1, H2_sp), kron(-transpose(H1_sp), I_m2))
        end
    else
        I_n1 = identity_matrix(F, n1)
        I_m2 = identity_matrix(F, m2)
        
        if Int(characteristic(F)) == 2
            H_Z = hcat(kronecker_product(I_n1, H2), kronecker_product(transpose(H1), I_m2))
        else
            H_Z = hcat(kronecker_product(I_n1, H2), kronecker_product(-transpose(H1), I_m2))
        end
    end
    
    S.cache[:H_Z] = H_Z
    return H_Z
end

function stabilizers(S::HypergraphProductCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    H_X = X_stabilizers(S)
    H_Z = Z_stabilizers(S)
    n = S.n
    F = S.cache[:F]
    
    is_sparse = H_X isa SparseMatrixCSC || H_Z isa SparseMatrixCSC
    
    # Assemble standard symplectic parity check matrix
    if is_sparse
        stabs = vcat(hcat(H_X, spzeros(F, size(H_X, 1), n)),
                     hcat(spzeros(F, size(H_Z, 1), n), H_Z))
    else
        stabs = vcat(hcat(H_X, zero_matrix(F, size(H_X, 1), n)),
                     hcat(zero_matrix(F, size(H_Z, 1), n), H_Z))
    end
    
    S.cache[:stabilizers] = stabs
    return stabs
end

function character_vector(S::Union{HypergraphProductCode, GeneralizedShorCode})
    haskey(S.cache, :char_vec) && return S.cache[:char_vec]
    
    # Default to an all-zero vector (+1 phases) lazily sized to the number of stabilizers
    stabs = stabilizers(S)
    F = S.cache[:F]
    char_vec = zero_matrix(F, 1, size(stabs, 1))
    
    S.cache[:char_vec] = char_vec
    return char_vec
end

# unable to yield quantum LDPC code families with non constant minimum distance
"""
    GeneralizedShorCode(C1::AbstractLinearCode, C2::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)
    BaconCasaccinoConstruction(C1::AbstractLinearCode, C2::AbstractLinearCode; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return the generalized Shor code of `C1` and `C2` with `C1⟂ ⊆ C2`.
"""
function GeneralizedShorCode(C1::AbstractLinearCode, C2::AbstractLinearCode)
    C1.F == C2.F || throw(ArgumentError("Codes must be over the same field."))
    
    # 1. Calculate Physical Qubits (n)
    n_new = C1.n * C2.n
    
    # 2. Calculate Logical Qubits (k)
    k_new = C1.k * C2.k
    
    # 3. Handle Exact Distance Bounds
    l_bound = min(C1.l_bound, C2.l_bound)
    u_bound = min(C1.u_bound, C2.u_bound)
    
    d_exact = missing
    if !ismissing(C1.d) && !ismissing(C2.d)
        d_exact = min(C1.d, C2.d)
        l_bound = d_exact
        u_bound = d_exact
    end
    
    cache = Dict{Symbol, Any}(:F => C1.F)
    
    return GeneralizedShorCode(
        C1, C2, 
        n_new, k_new, 
        missing, # r (gauge qubits, deferred until stabilizers are computed)
        d_exact, l_bound, u_bound, 
        cache
    )
end
BaconCasaccinoConstruction(C1::AbstractLinearCode, C2::AbstractLinearCode) =
    GeneralizedShorCode(C1, C2)

function gauge_operators(S::GeneralizedShorCode)
    haskey(S.cache, :gauge_operators) && return S.cache[:gauge_operators]
    
    F = S.cache[:F]
    H1 = parity_check_matrix(S.C1)
    H2 = parity_check_matrix(S.C2)
    
    m1, n1 = size(H1)
    m2, n2 = size(H2)
    n = S.n
    
    is_sparse = H1 isa SparseMatrixCSC || H2 isa SparseMatrixCSC
    F_one = F(1)
    
    if is_sparse
        H1_sp = H1 isa SparseMatrixCSC ? H1 : sparse(H1)
        H2_sp = H2 isa SparseMatrixCSC ? H2 : sparse(H2)
        
        I_n1 = spdiagm(0 => fill(F_one, n1))
        I_n2 = spdiagm(0 => fill(F_one, n2))
        
        G_X_mat = kron(H1_sp, I_n2)
        G_Z_mat = kron(I_n1, H2_sp)
        
        G_X_sym = hcat(G_X_mat, spzeros(F, size(G_X_mat, 1), n))
        G_Z_sym = hcat(spzeros(F, size(G_Z_mat, 1), n), G_Z_mat)
    else
        I_n1 = identity_matrix(F, n1)
        I_n2 = identity_matrix(F, n2)
        
        G_X_mat = kronecker_product(H1, I_n2)
        G_Z_mat = kronecker_product(I_n1, H2)
        
        G_X_sym = hcat(G_X_mat, zero_matrix(F, size(G_X_mat, 1), n))
        G_Z_sym = hcat(zero_matrix(F, size(G_Z_mat, 1), n), G_Z_mat)
    end
    
    gauge_ops = vcat(G_X_sym, G_Z_sym)
    
    S.cache[:gauge_operators] = gauge_ops
    return gauge_ops
end

function stabilizers(S::GeneralizedShorCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    # Subsystem codes extract stabilizers by finding the center of the gauge group.
    # We leverage your existing subsystem infrastructure (e.g., SubsystemCode constructor) 
    # to compute this dynamically from the gauge operators.
    gauge_ops = gauge_operators(S)
    
    # Assuming you have an internal function that computes the stabilizers from gauge operators.
    # If this logic natively lives inside your `SubsystemCode` constructor, we can invoke 
    # it temporarily to extract the required matrix.
    temp_code = SubsystemCode(gauge_ops)
    stabs = stabilizers(temp_code)
    
    # We can also lazily cache `r` now that we've done the math
    S.r = temp_code.r
    
    S.cache[:stabilizers] = stabs
    return stabs
end

function X_stabilizers(S::GeneralizedShorCode)
    haskey(S.cache, :H_X) && return S.cache[:H_X]
    stabs = stabilizers(S)
    H_X = stabs[:, 1:S.n]
    S.cache[:H_X] = H_X
    return H_X
end

function Z_stabilizers(S::GeneralizedShorCode)
    haskey(S.cache, :H_Z) && return S.cache[:H_Z]
    stabs = stabilizers(S)
    H_Z = stabs[:, S.n+1:end]
    S.cache[:H_Z] = H_Z
    return H_Z
end

"""
    HyperBicycleCodeCSS(a::Vector{CTMatrixTypes}, b::Vector{CTMatrixTypes}, χ::Int; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return the hyperbicycle CSS code of `a` and `b` given `χ`.

# Arguments
- a: A vector of length `c` of binary matrices of the same dimensions.
- b: A vector of length `c` of binary matrices of the same dimensions,
  potentially different from those of `a`.
- χ: A strictly positive integer coprime with `c`.

# Example

[[900, 50, 14]] CSS Hyperbicycle Code from Example 6 of [Kovalev_2013](@cite).

```jldoctest
julia> S, x = polynomial_ring(Oscar.Nemo.Native.GF(2), :x);

julia> l = 30; χ = 1;

julia> R, = residue_ring(S, x^l - 1);

julia> h = R(1 + x + x^3 + x^5);

julia> A = residue_polynomial_to_circulant_matrix(h);

julia> a1 = A[1:15, 1:15];

julia> a2 = A[1:15, 16:30];

julia> code = HyperBicycleCodeCSS([a1, a2], [a1, a2], χ);

julia> length(code), dimension(code)
(900, 50)
```
"""
function HyperBicycleCodeCSS(a::Vector{T}, b::Vector{T}, χ::Int) where T <: CTMatrixTypes
    χ > 0 || throw(ArgumentError("Required χ > 0."))
    
    c = length(a)
    gcd(c, χ) == 1 || throw(ArgumentError("The length of the input vectors must be coprime with χ."))
    c == length(b) || throw(ArgumentError("Input vectors must have same length."))
    
    k1, n1 = size(a[1])
    k2, n2 = size(b[1])
    F = base_ring(a[1])
    Int(order(F)) == 2 || throw(ArgumentError("Hyperbicycle codes require binary inputs."))
    
    for i in 1:c
        F == base_ring(a[i]) || throw(ArgumentError("Inputs must share the same base ring."))
        (k1, n1) == size(a[i]) || throw(ArgumentError("First set of matrices must all have the same dimensions."))
        F == base_ring(b[i]) || throw(ArgumentError("Inputs must share the same base ring."))
        (k2, n2) == size(b[i]) || throw(ArgumentError("Second set of matrices must all have the same dimensions."))
    end

    # Calculate physical qubits O(1)
    n_new = c * (k2 * n1 + n2 * k1)
    
    cache = Dict{Symbol, Any}(:F => F)    
    return HyperBicycleCodeCSS(
        a, b, χ, 
        n_new, missing, missing, 1, n_new, 
        cache
    )
end

"""
    HyperBicycleCode(a::Vector{CTMatrixTypes}, b::Vector{CTMatrixTypes}, χ::Int; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return the hyperbicycle non-CSS code of `a` and `b` given `χ`.

# Arguments
- a: A vector of length `c` of binary matrices of the same dimensions.
- b: A vector of length `c` of binary matrices of the same dimensions,
  potentially different from those of `a`.
- χ: A strictly positive integer coprime with `c`.

# Example

[[289, 81, 5]] non-CSS Hyperbicycle Code from Example 13 of [Kovalev_2013](@cite).

```jldoctest
julia> using CodingTheory, Oscar;

julia> S, x = polynomial_ring(Oscar.Nemo.Native.GF(2), :x);

julia> l = 17; χ = 1;

julia> R, = residue_ring(S, x^l - 1);

julia> h = R(x^4 * (1 + x + x^3 + x^6 + x^8 + x^9));

julia> A = residue_polynomial_to_circulant_matrix(h);

julia> code = HyperBicycleCode([A], [A], χ);

julia> length(code), dimension(code)
(289, 81)
```
"""
function HyperBicycleCode(a::Vector{T}, b::Vector{T}, χ::Int) where T <: CTMatrixTypes
    χ > 0 || throw(ArgumentError("Required χ > 0."))
    
    c = length(a)
    gcd(c, χ) == 1 || throw(ArgumentError("The length of the input vectors must be coprime with χ."))
    c == length(b) || throw(ArgumentError("Input vectors must have same length."))
    
    k1, n1 = size(a[1])
    k2, n2 = size(b[1])
    F = base_ring(a[1])
    Int(order(F)) == 2 || throw(ArgumentError("Hyperbicycle codes require binary inputs."))
    
    for i in 1:c
        F == base_ring(a[i]) || throw(ArgumentError("Inputs must share the same base ring."))
        (k1, n1) == size(a[i]) || throw(ArgumentError("First set of matrices must all have the same dimensions."))
        F == base_ring(b[i]) || throw(ArgumentError("Inputs must share the same base ring."))
        (k2, n2) == size(b[i]) || throw(ArgumentError("Second set of matrices must all have the same dimensions."))
    end

    # Calculate physical qubits O(1)
    n_new = c * k2 * n1 # For the non-CSS symplectic structure, column match dictates n
    
    cache = Dict{Symbol, Any}(:F => F)
    return HyperBicycleCode(
        a, b, χ, 
        n_new, missing, missing, 1, n_new, 
        cache
    )
end

# Helper function that executes the internal Hyperbicycle sum math
function _build_hyperbicycle_blocks(a, b, χ, F, c, k1, n1, k2, n2)
    H1 = zero_matrix(F, c * k1, c * n1)
    H2 = zero_matrix(F, c * k2, c * n2)
    HT1 = zero_matrix(F, c * n1, c * k1)
    HT2 = zero_matrix(F, c * n2, c * k2)
    
    Ic = identity_matrix(F, c)
    Sχ = Ic[mod1.(1:χ:c * χ, c), :]
    
    for i in 1:c
        Ii = vcat(Ic[i:end, :], Ic[1:i - 1, :])
        Iχi = Sχ * Ii
        ITχi = transpose(Sχ) * transpose(Ii)
        H1 += kronecker_product(Iχi, a[i])
        H2 += kronecker_product(b[i], Iχi)
        HT1 += kronecker_product(ITχi, transpose(a[i]))
        HT2 += kronecker_product(transpose(b[i]), ITχi)
    end
    return H1, H2, HT1, HT2
end

function X_stabilizers(S::HyperBicycleCodeCSS)
    haskey(S.cache, :H_X) && return S.cache[:H_X]
    
    a, b = S.a, S.b
    k1, n1 = size(a[1])
    k2, n2 = size(b[1])
    F = S.cache[:F]
    
    H1, H2, _, _ = _build_hyperbicycle_blocks(a, b, S.χ, F, length(a), k1, n1, k2, n2)
    
    Ek1 = identity_matrix(F, k1)
    Ek2 = identity_matrix(F, k2)
    
    H_X = hcat(kronecker_product(Ek2, H1), kronecker_product(H2, Ek1))
    S.cache[:H_X] = H_X
    return H_X
end

function Z_stabilizers(S::HyperBicycleCodeCSS)
    haskey(S.cache, :H_Z) && return S.cache[:H_Z]
    
    a, b = S.a, S.b
    k1, n1 = size(a[1])
    k2, n2 = size(b[1])
    F = S.cache[:F]
    
    _, _, HT1, HT2 = _build_hyperbicycle_blocks(a, b, S.χ, F, length(a), k1, n1, k2, n2)
    
    En1 = identity_matrix(F, n1)
    En2 = identity_matrix(F, n2)
    
    H_Z = hcat(kronecker_product(HT2, En1), kronecker_product(En2, HT1))
    S.cache[:H_Z] = H_Z
    return H_Z
end

function stabilizers(S::HyperBicycleCodeCSS)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    H_X = X_stabilizers(S)
    H_Z = Z_stabilizers(S)
    n = S.n
    F = S.cache[:F]
    
    stabs = vcat(hcat(H_X, zero_matrix(F, size(H_X, 1), n)),
                 hcat(zero_matrix(F, size(H_Z, 1), n), H_Z))
                 
    S.cache[:stabilizers] = stabs
    return stabs
end

function stabilizers(S::HyperBicycleCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    a, b = S.a, S.b
    k1, n1 = size(a[1])
    k2, n2 = size(b[1])
    F = S.cache[:F]
    
    H1, H2, HT1, HT2 = _build_hyperbicycle_blocks(a, b, S.χ, F, length(a), k1, n1, k2, n2)
    (H1 == HT1 && H2 == HT2) || throw(ArgumentError("H_i must equal H̃_i for i = 1, 2 for Non-CSS HyperBicycle."))
    
    Ek1 = identity_matrix(F, k1)
    Ek2 = identity_matrix(F, k2)

    stabs = hcat(kronecker_product(Ek2, H1), kronecker_product(H2, Ek1))
    S.cache[:stabilizers] = stabs
    return stabs
end

"""
    GeneralizedBicycleCode(A::CTMatrixTypes, B::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return the generealized bicycle code given by `A` and `B`.

# Example

[[254, 28, 14 ≤ d ≤ 20]] Generalized Bicycle Code from Appendix B, Example A1 of [panteleev2021degenerate](@cite).

```jldoctest
julia> using CodingTheory, Oscar;

julia> F = Oscar.Nemo.Native.GF(2);

julia> S, x = polynomial_ring(F, :x);

julia> l = 127;

julia> R, _ = residue_ring(S, x^l - 1);

julia> a = 1 + x^15 + x^20 + x^28 + x^66;

julia> b = 1 + x^58 + x^59 + x^100 + x^121;

julia> code = GeneralizedBicycleCode(R(a), R(b));

julia> length(code), dimension(code)
(254, 28)
```
"""
function GeneralizedBicycleCode(A::T, B::T) where T <: CTMatrixTypes
    F = base_ring(A)
    F == base_ring(B) || throw(ArgumentError("Arguments must be over the same base ring."))
    (iszero(A) || iszero(B)) && throw(ArgumentError("Arguments should not be zero."))
    
    m, n1 = size(A)
    m2, n2 = size(B)
    m == n1 == m2 == n2 || throw(ArgumentError("A and B must be square matrices of the same dimensions."))
    iszero(A * B - B * A) || throw(ArgumentError("Arguments must commute."))

    # Physical qubits n = 2m since H_X = [A | B]
    n_new = 2 * m
    
    cache = Dict{Symbol, Any}(:F => F)
    return GeneralizedBicycleCode(
        A, B, 
        n_new, missing, missing, 1, n_new, 
        cache
    )
end

"""
    GeneralizedBicycleCode(a::T, b::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: ResElem

Return the generealized bicycle code determined by `a` and `b`.

# Notes
- `l x l` circulant matrices are constructed using the coefficients of the polynomials
  `a` and `b` in `F_q[x]/(x^l - 1)` (`gcd(q, l) = 1`) as the first column
"""
function GeneralizedBicycleCode(a::T, b::T) where T <: ResElem
    parent(a) == parent(b) || throw(ArgumentError("Both objects must be defined over the same residue ring."))

    S = GeneralizedBicycleCode(residue_polynomial_to_circulant_matrix(a),
        residue_polynomial_to_circulant_matrix(b))
    S.cache[:polynomials] = (a, b)
    return S
end

"""
    GeneralizedBicycleCode(a::T, b::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: CTGroupAlgebra

Return the generealized bicycle code determined by `a` and `b`.

# Notes
- `|G| x |G|` circulant matrices are constructed using the coefficients of the elements in the group algebra `FG` as` the first column
"""
function GeneralizedBicycleCode(a::T, b::T) where T <: CTGroupAlgebra
    parent(a) == parent(b) || throw(ArgumentError("Both objects must be defined over the same residue ring."))

    S = GeneralizedBicycleCode(group_algebra_element_to_circulant_matrix(a),
        group_algebra_element_to_circulant_matrix(b))
    S.cache[:polynomials] = (a, b)
    return S
end

function X_stabilizers(S::GeneralizedBicycleCode)
    haskey(S.cache, :H_X) && return S.cache[:H_X]
    
    # H_X = [A | B]
    H_X = hcat(S.A, S.B)
    S.cache[:H_X] = H_X
    return H_X
end

function Z_stabilizers(S::GeneralizedBicycleCode)
    haskey(S.cache, :H_Z) && return S.cache[:H_Z]
    
    F = S.cache[:F]
    
    # H_Z = [B^T | ±A^T]
    if Int(order(F)) == 2
        H_Z = hcat(transpose(S.B), transpose(S.A))
    else
        H_Z = hcat(transpose(S.B), -transpose(S.A))
    end
    
    S.cache[:H_Z] = H_Z
    return H_Z
end

function stabilizers(S::GeneralizedBicycleCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    H_X = X_stabilizers(S)
    H_Z = Z_stabilizers(S)
    n = S.n
    F = S.cache[:F]
    
    is_sparse = H_X isa SparseMatrixCSC || H_Z isa SparseMatrixCSC
    
    if is_sparse
        stabs = vcat(hcat(H_X, spzeros(F, size(H_X, 1), n)),
                     hcat(spzeros(F, size(H_Z, 1), n), H_Z))
    else
        stabs = vcat(hcat(H_X, zero_matrix(F, size(H_X, 1), n)),
                     hcat(zero_matrix(F, size(H_Z, 1), n), H_Z))
    end
                 
    S.cache[:stabilizers] = stabs
    return stabs
end

"""
    BicycleCode(A::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return a lazy Bicycle code given by the square matrix `A`. 
This is equivalent to a Generalized Bicycle Code where `B = A^T`.
"""
function BicycleCode(A::CTMatrixTypes; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm)
    
    m, n = size(A)
    m == n || throw(ArgumentError("Input matrix must be square."))
    
    F = base_ring(A)
    Int(order(F)) == 2 || throw(ArgumentError("Bicycle codes are strictly defined over GF(2)."))
    
    # A Bicycle code is exactly a Generalized Bicycle Code with B = A^T
    return GeneralizedBicycleCode(A, transpose(A); char_vec = char_vec, logs_alg = logs_alg)
end

"""
    BicycleCode(a::ResElem; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return the lazy Bicycle code determined by the residue ring element `a`.
"""
function BicycleCode(a::ResElem; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm)
    
    A = residue_polynomial_to_circulant_matrix(a)
    
    S = BicycleCode(A; char_vec = char_vec, logs_alg = logs_alg)
    # Save the generating polynomial for researchers
    S.cache[:polynomials] = (a, )
    return S
end

"""
    BicycleCode(a::CTGroupAlgebra; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm)

Return the lazy Bicycle code determined by the group algebra element `a`.
"""
function BicycleCode(a::CTGroupAlgebra; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm)
    
    A = group_algebra_element_to_circulant_matrix(a)
    
    S = BicycleCode(A; char_vec = char_vec, logs_alg = logs_alg)
    # Save the generating group algebra element for researchers
    S.cache[:polynomials] = (a, )
    return S
end

# """
#     generalized_hypergraph_product_matrices(A::MatElem{T}, b::T) where T <: ResElem
#     GHGP_matrices(A::MatElem{T}, b::T) where T <: ResElem
#     lifted_product_matrices(A::MatElem{T}, b::T) where T <: ResElem

# Return the pre-lifted matrices `H_X` and `H_Z` of the generalized hypergraph product code of `A` and `b`.

# # Arguments
# - `A` - an `m x n` matrix with elements in `F_2[x]/(x^m - 1)`
# - `b` - a polynomial over the same residue ring

# # Notes
# - Use `LiftedProductCode` to return a quantum code over the base ring directly.
# """
# function generalized_hypergraph_product_matrices(A::MatElem{T}, b::T) where T <: ResElem

#     A * b == b * A || throw(ArgumentError("A and b must commute to form a valid lifted product."))
#     S = base_ring(b)
#     F = base_ring(S)
#     # Int(order(F)) == 2 || throw(ArgumentError("The generalized hypergraph product is only defined over GF(2)."))
#     R = parent(A[1, 1])
#     R == parent(b) || throw(ArgumentError("Both objects must be defined over the same residue ring."))
#     m, n = size(A)
#     (m != 1 && n != 1) || throw(ArgumentError("First input matrix must not be a vector."))
#     f = modulus(R)
#     l = degree(f)
#     f == gen(S)^l - 1 || throw(ArgumentError("Residue ring not of the form x^l - 1."))
#     # gcd(l, Int(characteristic(F))) == 1 || throw(ArgumentError("Residue ring over F_q[x] must be defined by x^l - 1 with gcd(l, q) = 1."))
    
#     A_tr = _CT_adjoint(A)
#     # b_coeffs = collect(coefficients(Nemo.lift(b)))
#     # for _ in 1:l - length(b_coeffs)
#     #     push!(b_coeffs, F(0))
#     # end
#     # b_coeffs[2:end] = reverse(b_coeffs[2:end])
#     # B_tr = R(S(b_coeffs))
#     B_tr = _CT_adjoint(matrix(R, 1, 1, [b]))[1, 1]
#     Mn = matrix_space(R, n, n)
#     H_Z = hcat(Mn(B_tr), A_tr)
#     Mm = matrix_space(R, m, m)
#     # TODO: check extending this function past F_2 makes sense
#     # branch for speed
#     if Int(order(F)) == 2
#         H_X = hcat(A, Mm(b))
#     else
#         H_X = hcat(A, -Mm(b))
#     end
#     return H_X, H_Z
# end
# GHGP_matrices(A::MatElem{T}, b::T) where T <: ResElem =
#     generalized_hypergraph_product_matrices(A, b)
# lifted_product_matrices(A::MatElem{T}, b::T) where T <: ResElem =
#     generalized_hypergraph_product_matrices(A, b)

# """
#     generalized_hypergraph_product_matrices(A::MatElem{T}, b::T) where T <: CTGroupAlgebra
#     GHGP_matrices(A::MatElem{T}, b::T) where T <: CTGroupAlgebra
#     lifted_product_matrices(A::MatElem{T}, b::T) where T <: CTGroupAlgebra

# Return the pre-lifted matrices `H_X` and `H_Z` of the generalized hypergraph product code of `A` and `b`.

# # Arguments
# - `A` - a matrix with elements in a group algebra
# - `b` - a polynomial over the same group algebra

# # Notes
# - Use `LiftedProductCode` to return a quantum code over the base ring directly.
# """
# function generalized_hypergraph_product_matrices(A::MatElem{T}, b::T) where T <: CTGroupAlgebra
#     FG = parent(A[1, 1])
#     parent(b) == FG || throw(ArgumentError("Inputs must be over the same group algebra"))

#     A_tr = _CT_adjoint(A)
#     B_tr = _CT_adjoint(matrix(FG, 1, 1, [b]))[1, 1]
#     m, n = size(A)
#     Mn = matrix_space(FG, n, n)
#     H_Z = hcat(Mn(B_tr), A_tr)
#     Mm = matrix_space(FG, m, m)
#     F = base_ring(FG)
#     # TODO: check extending this function past F_2 makes sense
#     # branch for speed
#     if Int(order(F)) == 2
#         H_X = hcat(A, Mm(b))
#     else
#         H_X = hcat(A, -Mm(b))
#     end
#     return H_X, H_Z
# end
# GHGP_matrices(A::MatElem{T}, b::T) where T <: CTGroupAlgebra =
#     generalized_hypergraph_product_matrices(A, b)
# lifted_product_matrices(A::MatElem{T}, b::T) where T <: CTGroupAlgebra =
#     generalized_hypergraph_product_matrices(A, b)

# """
#     GeneralizedHypergraphProductCode(A::MatElem{T}, b::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}
#     LiftedProductCode(A::MatElem{T}, b::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
#         logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

# Return the lifted (generalized hypergraph) product code of `A` and `b`.

# # Arguments
# - `A` - either an `m x n` matrix with elements in `F_2[x]/(x^m - 1)` or a group algebra
# - `b` - a polynomial over the same residue ring or group algebra
# """
# function GeneralizedHypergraphProductCode(A::MatElem{T}, b::T; char_vec::Union{Vector{zzModRingElem},
#     Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

#     H_X, H_Z = generalized_hypergraph_product_matrices(A, b)
#     return CSSCode(lift(H_X), lift(H_Z), char_vec = char_vec, logs_alg = logs_alg)
# end
# LiftedProductCode(A::MatElem{T}, b::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
#     logs_alg::Symbol = :stnd_frm) where T <: ResElem =
#     GeneralizedHypergraphProductCode(A, b, char_vec = char_vec, logs_alg = logs_alg)

# """
#     lifted_product_matrices(A::MatElem{T}, B::MatElem{T}) where T <: ResElem

# Return the pre-lifted matrices `H_X` and `H_Z` for the lifted quasi-cyclic lifted product code.

# # Arguments
# - `A` - an `m x n1` matrix with elements in `F_2[x]/(x^m - 1)`
# - `B` - an `m x n2` matrix with elements in the same residue ring

# # Notes
# - Use `LiftedProductCode` to return a quantum code over the base ring directly.
# """
# function lifted_product_matrices(A::MatElem{T}, B::MatElem{T}) where T <: ResElem
#     A * B == B * A || throw(ArgumentError("A and B must commute to form a valid lifted product."))
#     S = base_ring(A[1, 1])
#     F = base_ring(S)
#     Int(order(F)) == 2 || throw(ArgumentError("The quasi-cyclic lifted product is only defined over GF(2)."))
#     R = parent(A[1, 1])
#     R == parent(B[1, 1]) || throw(ArgumentError("Both objects must be defined over the same residue ring."))
#     f = modulus(R)
#     l = degree(f)
#     f == gen(S)^l - 1 || throw(ArgumentError("Residue ring not of the form x^l - 1."))
    
#     A_tr = _CT_adjoint(A)
#     B_tr = _CT_adjoint(B)

#     k1, n1 = size(A)
#     k2, n2 = size(B)
#     Ek1 = identity_matrix(R, k1)
#     Ek2 = identity_matrix(R, k2)
#     En1 = identity_matrix(R, n1)
#     En2 = identity_matrix(R, n2)

#     H_X = hcat(A ⊗ Ek2, Ek1 ⊗ B)
#     H_Z = hcat(En1 ⊗ B_tr, A_tr ⊗ En2)
#     return H_X, H_Z
# end

# """
#     lifted_product_matrices(A::MatElem{T}, B::MatElem{T}) where T <: ResElem

# Return the pre-lifted matrices `H_X` and `H_Z` for the lifted quasi-cyclic lifted product code.

# # Arguments
# - `A` - a matrix with elements in a group algebra
# - `B` - a matrix with coefficents in the same group algebra

# # Notes
# - Use `LiftedProductCode` to return a quantum code over the base ring directly.
# """
# function lifted_product_matrices(A::MatElem{T}, B::MatElem{T}) where T <: CTGroupAlgebra
#     FG = parent(A[1, 1])
#     parent(B[1, 1]) == FG || throw(ArgumentError("Inputs must be over the same group algebra"))

#     A_tr = _CT_adjoint(A)
#     B_tr = _CT_adjoint(B)

#     k1, n1 = size(A)
#     k2, n2 = size(B)
#     Ek1 = identity_matrix(FG, k1)
#     Ek2 = identity_matrix(FG, k2)
#     En1 = identity_matrix(FG, n1)
#     En2 = identity_matrix(FG, n2)

#     H_X = hcat(A ⊗ Ek2, Ek1 ⊗ B)
#     H_Z = hcat(En1 ⊗ B_tr, A_tr ⊗ En2)
#     return H_X, H_Z
# end

# """
#     LiftedProductCode(A::MatElem{T}, B::MatElem{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

# Return the lifted product code given by the matrices `A` and `B`.

# # Example

# [[882, 24, 18 ≤ d ≤ 24]] Lifted Product Code from Appendix B, Example B1 of [panteleev2021degenerate](@cite).

# ```jldoctest
# julia> using CodingTheory, Oscar;

# julia> F = Oscar.Nemo.Native.GF(2);

# julia> S, x = polynomial_ring(F, :x);

# julia> l = 63;

# julia> R, _ = residue_ring(S, x^l - 1);

# julia> A = matrix(R, 7, 7,
#            [x^27, 0   , 0   , 0   , 0   , 1   , x^54,
#             x^54, x^27, 0   , 0   , 0   , 0   , 1   ,
#             1   , x^54, x^27, 0   , 0   , 0   , 0   ,
#             0   , 1   , x^54, x^27, 0   , 0   , 0   ,
#             0   , 0   , 1   , x^54, x^27, 0   , 0   ,
#             0   , 0   , 0   , 1   , x^54, x^27, 0   ,
#             0   , 0   , 0   , 0   , 1   , x^54, x^27]);

# julia> b = R(1 + x + x^6);

# julia> code = LiftedProductCode(A, b);
# ┌ Warning: Commutativity of A and b required but not yet enforced.
# └ @ CodingTheory ~/Documents/GitHub/CodingTheory/src/Quantum/product_codes.jl:340

# julia> length(code), dimension(code)
# (882, 24)
# ```
# """
# function LiftedProductCode(A::MatElem{T}, B::MatElem{T}; char_vec::Union{Vector{zzModRingElem}, Missing} =
#     missing, logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

#     H_X, H_Z = lifted_product_matrices(A, B)
#     return CSSCode(lift(H_X), lift(H_Z), char_vec = char_vec, logs_alg = logs_alg)
# end

# """
#     bias_tailored_lifted_product_matrices(A::MatElem{T}, B::MatElem{T}) where T <: ResElem

# Return the pre-lifted stabilizer matrix for bias-tailored lifted product code of `A` and `B`.

# # Arguments
# - `A` - an `m x n1` matrix with elements in `F_2[x]/(x^m - 1)`
# - `B` - an `m x n2` matrix with elements in the same residue ring

# # Notes
# - Use `BiasTailoredLiftedProductCode` to return a quantum code over the base ring directly.
# """
# function bias_tailored_lifted_product_matrices(A::MatElem{T}, B::MatElem{T}) where T <: ResElem
#     A * B == B * A || throw(ArgumentError("A and B must commute to form a valid lifted product."))
#     S = base_ring(A[1, 1])
#     F = base_ring(S)
#     Int(order(F)) == 2 || throw(ArgumentError("The quasi-cyclic lifted product is only defined over GF(2)."))
#     R = parent(A[1, 1])
#     R == parent(B[1, 1]) || throw(ArgumentError("Both objects must be defined over the same residue ring."))
#     f = modulus(R)
#     l = degree(f)
#     f == gen(S)^l - 1 || throw(ArgumentError("Residue ring not of the form x^l - 1."))
    
#     A_tr = _CT_adjoint(A)
#     B_tr = _CT_adjoint(B)

#     k1, n1 = size(A)
#     k2, n2 = size(B)
#     Ek1 = identity_matrix(R, k1)
#     Ek2 = identity_matrix(R, k2)
#     En1 = identity_matrix(R, n1)
#     En2 = identity_matrix(R, n2)

#     A12 = A_tr ⊗ Ek2
#     A13 = En1 ⊗ B
#     A21 = A ⊗ En2
#     A24 = Ek1 ⊗ B_tr
#     return vcat(hcat(zero(A21), A12, A13, zero(A24)), hcat(A21, zero(A12), zero(A13), A24))
# end

# """
#     bias_tailored_lifted_product_matrices(A::MatElem{T}, B::MatElem{T}) where T <: CTGroupAlgebra

# Return the pre-lifted stabilizer matrix for bias-tailored lifted product code of `A` and `B`.

# # Arguments
# - `A` - a matrix with elements in a group algebra
# - `B` - a matrix with elements in the same group algebra

# # Notes
# - Use `BiasTailoredLiftedProductCode` to return a quantum code over the base ring directly.

# # Example

# [[882, 24, d ≤ 24]] BiasTailored Lifted Product Code from Appendix B of [roffe2023bias](@cite).

# ```jldoctest
# julia> using CodingTheory, Oscar;

# julia> F = Oscar.Nemo.Native.GF(2);

# julia> S, x = polynomial_ring(F, :x);

# julia> l = 63;

# julia> R, _ = residue_ring(S, x^l - 1);

# julia> A1 = matrix(R, 1, 1, [1 + x^1 + x^6]);

# julia> A2 = matrix(R, 7, 7,
#            [x^36, 0   , 0   , 0   , 0   , 1   , x^9 ,
#             x^9 , x^36, 0   , 0   , 0   , 0   , 1   ,
#             1   , x^9 , x^36, 0   , 0   , 0   , 0   ,
#             0   , 1   , x^9 , x^36, 0   , 0   , 0   ,
#             0   , 0   , 1   , x^9 , x^36, 0   , 0   ,
#             0   , 0   , 0   , 1   , x^9 , x^36, 0   ,
#             0   , 0   , 0   , 0   , 1   , x^9 , x^36]);

# julia> code = BiasTailoredLiftedProductCode(A1, A2);
# ┌ Warning: Commutativity of A and b required but not yet enforced.
# └ @ CodingTheory ~/Documents/GitHub/CodingTheory/src/Quantum/product_codes.jl:60

# julia> length(code), dimension(code)
# (882, 24)
# ```
# """
# function bias_tailored_lifted_product_matrices(A::MatElem{T}, B::MatElem{T}) where T <: CTGroupAlgebra

#     FG = parent(A[1, 1])
#     parent(B[1, 1]) == FG || throw(ArgumentError("Inputs must be over the same group algebra"))

#     A_tr = _CT_adjoint(A)
#     B_tr = _CT_adjoint(B)

#     k1, n1 = size(A)
#     k2, n2 = size(B)
#     Ek1 = identity_matrix(R, k1)
#     Ek2 = identity_matrix(R, k2)
#     En1 = identity_matrix(R, n1)
#     En2 = identity_matrix(R, n2)

#     A12 = A_tr ⊗ Ek2
#     A13 = En1 ⊗ B
#     A21 = A ⊗ En2
#     A24 = Ek1 ⊗ B_tr
#     return vcat(hcat(zero(A21), A12, A13, zero(A24)), hcat(A21, zero(A12), zero(A13), A24))
# end

# """
#     BiasTailoredLiftedProductCode(A::MatElem{T}, B::MatElem{T}; char_vec::Union{Vector{zzModRingElem}, Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

# Return the bias-tailored lifted product code of `A` and `B`.

# # Arguments
# - `A` - either an `m x n` matrix with elements in `F_2[x]/(x^m - 1)` or a group algebra
# - `B` - an `m x n2` matrix with elements in the same parent as `A`
# """
# function BiasTailoredLiftedProductCode(A::MatElem{T}, B::MatElem{T}; char_vec::Union{Vector{zzModRingElem},
#     Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

#     stabs = bias_tailored_lifted_product_matrices(A, B)
#     return StabilizerCode(lift(stabs), char_vec = char_vec, logs_alg = logs_alg)
# end

function LiftedProductCode(A::MatElem{T}, B::MatElem{T}; char_vec::Union{Vector{zzModRingElem}, Missing} =
    missing, logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

    A * B == B * A || throw(ArgumentError("A and B must commute to form a valid lifted product."))
    
    # Lazily determine the lifting size `l` by evaluating a single element
    R = parent(A[1, 1])
    test_lift = lift(matrix(R, 1, 1, [A[1, 1]]))
    l = size(test_lift, 1)
    F = base_ring(test_lift)
    
    k1, n1 = size(A)
    k2, n2 = size(B)
    
    # Calculate physical qubits O(1)
    n_pre = k2 * n1 + n2 * k1
    n_new = l * n_pre
    
    cache = Dict{Symbol, Any}(:F => F, :logs_alg => logs_alg)
    !ismissing(char_vec) && (cache[:char_vec] = char_vec)
    
    return LiftedProductCode(A, B, n_new, missing, missing, 1, n_new, cache)
end

function LiftedProductCode(A::MatElem{T}, b::T; char_vec::Union{Vector{zzModRingElem}, Missing} = missing,
    logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

    A * b == b * A || throw(ArgumentError("A and b must commute to form a valid lifted product."))
    
    # A Generalized Hypergraph Product is exactly a Lifted Product where B is a 1x1 matrix.
    R = parent(A[1, 1])
    B_mat = matrix(R, 1, 1, [b])
    
    return LiftedProductCode(A, B_mat; char_vec = char_vec, logs_alg = logs_alg)
end
GeneralizedHypergraphProductCode(A, b; kwargs...) = LiftedProductCode(A, b; kwargs...)

function BiasTailoredLiftedProductCode(A::MatElem{T}, B::MatElem{T}; char_vec::Union{Vector{zzModRingElem},
    Missing} = missing, logs_alg::Symbol = :stnd_frm) where T <: Union{ResElem, CTGroupAlgebra}

    R = parent(A[1, 1])
    test_lift = lift(matrix(R, 1, 1, [A[1, 1]]))
    l = size(test_lift, 1)
    F = base_ring(test_lift)
    
    k1, n1 = size(A)
    k2, n2 = size(B)
    
    # Calculate physical qubits O(1)
    n_pre = 2 * (n1 * n2 + k1 * k2)
    n_new = l * n_pre
    
    cache = Dict{Symbol, Any}(:F => F, :logs_alg => logs_alg)
    !ismissing(char_vec) && (cache[:char_vec] = char_vec)
    
    return BiasTailoredLiftedProductCode(A, B, n_new, missing, missing, 1, n_new, cache)
end

function X_stabilizers(S::LiftedProductCode)
    haskey(S.cache, :H_X) && return S.cache[:H_X]
    
    A, B = S.A, S.B
    R = parent(A[1, 1])
    k1, n1 = size(A)
    k2, n2 = size(B)
    F = S.cache[:F]
    
    Ek1 = identity_matrix(R, k1)
    Ek2 = identity_matrix(R, k2)
    
    if Int(characteristic(F)) == 2
        H_X_pre = hcat(kronecker_product(A, Ek2), kronecker_product(Ek1, B))
    else
        H_X_pre = hcat(kronecker_product(A, Ek2), kronecker_product(-Ek1, B))
    end
    
    H_X = lift(H_X_pre)
    S.cache[:H_X] = H_X
    return H_X
end

function Z_stabilizers(S::LiftedProductCode)
    haskey(S.cache, :H_Z) && return S.cache[:H_Z]
    
    A, B = S.A, S.B
    R = parent(A[1, 1])
    k1, n1 = size(A)
    k2, n2 = size(B)
    
    A_tr = _CT_adjoint(A)
    B_tr = _CT_adjoint(B)
    
    En1 = identity_matrix(R, n1)
    En2 = identity_matrix(R, n2)
    
    H_Z_pre = hcat(kronecker_product(En1, B_tr), kronecker_product(A_tr, En2))
    
    H_Z = lift(H_Z_pre)
    S.cache[:H_Z] = H_Z
    return H_Z
end

function stabilizers(S::LiftedProductCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    H_X = X_stabilizers(S)
    H_Z = Z_stabilizers(S)
    n = S.n
    F = S.cache[:F]
    
    stabs = vcat(hcat(H_X, zero_matrix(F, size(H_X, 1), n)),
                 hcat(zero_matrix(F, size(H_Z, 1), n), H_Z))
                 
    S.cache[:stabilizers] = stabs
    return stabs
end

function stabilizers(S::BiasTailoredLiftedProductCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    A, B = S.A, S.B
    R = parent(A[1, 1])
    k1, n1 = size(A)
    k2, n2 = size(B)
    
    A_tr = _CT_adjoint(A)
    B_tr = _CT_adjoint(B)
    
    Ek1 = identity_matrix(R, k1)
    Ek2 = identity_matrix(R, k2)
    En1 = identity_matrix(R, n1)
    En2 = identity_matrix(R, n2)
    
    A12 = kronecker_product(A_tr, Ek2) 
    A13 = kronecker_product(En1, B)    
    A21 = kronecker_product(A, En2)    
    A24 = kronecker_product(Ek1, B_tr) 
    
    # Safe zero blocks avoiding strict identity mirroring
    Z1 = zero_matrix(R, nrows(A12), ncols(A21))
    Z2 = zero_matrix(R, nrows(A13), ncols(A24))
    Z3 = zero_matrix(R, nrows(A21), ncols(A12))
    Z4 = zero_matrix(R, nrows(A24), ncols(A13))
    
    H_pre = vcat(hcat(Z1, A12, A13, Z2), 
                 hcat(A21, Z3, Z4, A24))
    
    stabs = lift(H_pre)
    S.cache[:stabilizers] = stabs
    return stabs
end

"""
    SPCDFoldProductCode(D::Int, s::Int = 1)
    SingleParityCheckDFoldProductCode(D::Int, s::Int = 1) = SPCDFoldProductCode(D, s)

Return the single-parity-check `D`-fold product code.

# Note
- This is defined in https://arxiv.org/abs/2209.13474

# Example

[512, 174, 8]] Symmetric 2-fold product CSS code from [ostrev2024classical](@cite)

```jldoctest
julia> using CodingTheory, Oscar;

julia> F = Oscar.Nemo.Native.GF(2);

julia> h = matrix(F, [1 1]);

julia> id = identity_matrix(F, 2);

julia> H_X = vcat(
             h ⊗ h ⊗ h ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id,
             id ⊗ id ⊗ id ⊗ h ⊗ h ⊗ h ⊗ id ⊗ id ⊗ id,
             id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ id ⊗ h ⊗ h ⊗ h);

julia> H_Z = vcat(
             h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id,
             id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id,
             id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h ⊗ id ⊗ id ⊗ h);

julia> code = SPCDFoldProductCode(3);

julia> length(code), dimension(code)
(512, 174)
```
"""
function SPCDFoldProductCode(D::Int, s::Int = 1)
    vec_S = Vector{AbstractCSSCode}() # Or whatever explicit CSS supertype you map it to
    for i in 1:D
        for l in 1:D
            if l == (i - 1) * D + i
                push!(vec_S, CSSCode(RepetitionCode(2, 2 * s)))
            else
                push!(vec_S, CSSCode(RepetitionCode(2, 2)))
            end
        end
    end
    
    S = symmetric_product(vec_S)
    
    # Inject exact distance boundaries manually 
    dist = 2^D
    S.d = dist
    S.l_bound = dist
    S.u_bound = dist
    S.cache[:pure] = true
    
    return S
end
SingleParityCheckDFoldProductCode(D::Int, s::Int = 1) = SPCDFoldProductCode(D, s)

#############################
      # getter functions
#############################

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

"""
    Quintavalle_basis(C::HypergraphProductCode)

Return a symplectic canonical basis for the logical operators of `C`.
    
# Note
- This implements https://doi.org/10.48550/arXiv.2204.10812.
"""
function Quintavalle_basis(C::HypergraphProductCode)
    H1 = parity_check_matrix(C.C1)
    H2 = parity_check_matrix(C.C2)

    # c - complement
    ker_H1, im_H1_tr_c = strongly_lower_triangular_reduction(H1)
    ker_H1_tr, im_H1_c = strongly_lower_triangular_reduction(transpose(H1))
    ker_H2, im_H2_tr_c = strongly_lower_triangular_reduction(H2)
    ker_H2_tr, im_H2_c = strongly_lower_triangular_reduction(transpose(H2))
    F = C.F
    lx = zero_matrix(F, C.k, C.n)
    lz = deepcopy(lx)

    l = 1
    temp = zero_matrix(F, 1, nrows(ker_H1_tr) * nrows(ker_H2_tr))
    tr_im_H1_tr_c = transpose(im_H1_tr_c)
    tr_ker_H2 = transpose(ker_H2)
    tr_ker_H1 = transpose(ker_H1)
    tr_im_H2_tr_c = transpose(im_H2_tr_c)
    for i in 1:nrows(tr_ker_H1)
        for h in 1:nrows(tr_ker_H2)
            lx[l:l, :] = hcat(tr_im_H1_tr_c[i:i, :] ⊗ tr_ker_H2[h:h, :], temp)
            lz[l:l, :] = hcat(tr_ker_H1[i:i, :] ⊗ tr_im_H2_tr_c[h:h, :], temp)
            l += 1
        end
    end

    temp = zero_matrix(F, 1, nrows(ker_H1) * nrows(ker_H2))
    tr_ker_H1_tr = transpose(ker_H1_tr)
    tr_im_H2_c = transpose(im_H2_c)
    tr_im_H1_c = transpose(im_H1_c)
    tr_ker_H2_tr = transpose(ker_H2_tr)
    for i in 1:nrows(tr_ker_H1_tr)
        for h in 1:nrows(tr_ker_H2_tr)
            lx[l:l, :] = hcat(temp, tr_ker_H1_tr[i:i, :] ⊗ tr_im_H2_c[h:h, :])
            lz[l:l, :] = hcat(temp, tr_im_H1_c[i:i, :] ⊗ tr_ker_H2_tr[h:h, :])
            l += 1
        end
    end
    return lx, lz
end

# TODO: present the stabilizers in the docs and mention how to switch X and Z by
# using the switch on the inputs beforehand
"""
    asymmetric_product(S1::T, S2::T) where {T <: AbstractSubsystemCode}

Return the asymmetric 2-fold product quantum CSS code of the CSS codes `S1` and `S2`.

# Note
- This is defined in https://arxiv.org/abs/2209.13474
"""
function asymmetric_product(::IsCSS, S1::AbstractSubsystemCode, S2::AbstractSubsystemCode)
    F = S1.F
    F == S2.F || throw(ArgumentError("Base rings must match."))
    
    # Calculate physical qubits O(1)
    n_new = S1.n * S2.n
    cache = Dict{Symbol, Any}(:F => F)
    
    return AsymmetricProductCode(S1, S2, n_new, missing, missing, 1, n_new, cache)
end
asymmetric_product(S1::T, S2::T) where {T <: AbstractSubsystemCode} = asymmetric_product(CSSTrait(T), S1, S2)
asymmetric_product(::IsNotCSS, S1::AbstractSubsystemCode, S2::AbstractSubsystemCode) = error("Only valid for CSS codes.")

"""
    symmetric_product(vec_S::Vector{T}) where {T <: AbstractSubsystemCode}

Return the symmetric `D`-fold product quantum CSS code, where `D` is
the square-root of the length of the vector of CSS codes `vec_S`.

# Note
- This is defined in https://arxiv.org/abs/2209.13474
"""
function symmetric_product(::IsCSS, vec_S::Vector{T}) where {T <: AbstractSubsystemCode}
    length(vec_S) >= 4 || throw(DomainError("The length of the input vector must be at least 4"))
    D_float = sqrt(length(vec_S))
    isinteger(D_float) ? (D = Int(D_float)) : throw(ArgumentError("The number of CSS codes must be D^2"))
    
    F = vec_S[1].F
    n_new = 1
    
    # Calculate physical qubits O(1) by multiplying all lengths
    for S in vec_S
        S.F == F || throw(ArgumentError("All codes must share the same base ring."))
        n_new *= S.n
    end
    
    cache = Dict{Symbol, Any}(:F => F)
    return SymmetricProductCode(vec_S, D, n_new, missing, missing, 1, n_new, cache)
end
function symmetric_product(vec_S::Vector{T}) where {T <: AbstractSubsystemCode}
    isempty(vec_S) && throw(ArgumentError("Input vector of CSS codes cannot be empty"))
    for S in vec_S
        if CSSTrait(typeof(S)) == IsNotCSS()
            return symmetric_product(IsNotCSS(), vec_S)
        end
    end
    return symmetric_product(IsCSS(), vec_S)
end
symmetric_product(::IsNotCSS, vec_S::Vector{T}) where {T <: AbstractSubsystemCode} = error("Only valid for CSS codes.")

function X_stabilizers(S::AsymmetricProductCode)
    haskey(S.cache, :H_X) && return S.cache[:H_X]
    F = S.cache[:F]
    H_X = vcat(kronecker_product(X_stabilizers(S.S1), identity_matrix(F, S.S2.n)), 
               kronecker_product(identity_matrix(F, S.S1.n), X_stabilizers(S.S2)))
    S.cache[:H_X] = H_X
    return H_X
end

function Z_stabilizers(S::AsymmetricProductCode)
    haskey(S.cache, :H_Z) && return S.cache[:H_Z]
    H_Z = kronecker_product(Z_stabilizers(S.S1), Z_stabilizers(S.S2))
    S.cache[:H_Z] = H_Z
    return H_Z
end

function stabilizers(S::AsymmetricProductCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    H_X = X_stabilizers(S)
    H_Z = Z_stabilizers(S)
    F = S.cache[:F]
    stabs = vcat(hcat(H_X, zero_matrix(F, nrows(H_X), S.n)),
                 hcat(zero_matrix(F, nrows(H_Z), S.n), H_Z))
    S.cache[:stabilizers] = stabs
    return stabs
end

# --- Symmetric Product Accessors ---
function X_stabilizers(S::SymmetricProductCode)
    haskey(S.cache, :H_X) && return S.cache[:H_X]
    D = S.D
    F = S.cache[:F]
    vec_S = S.vec_S
    
    H_X = nothing
    for j in 0:D - 1
        temp_row = nothing
        for l in 1:D^2
            mat = (j * D + 1 <= l <= (j + 1) * D) ? X_stabilizers(vec_S[l]) : identity_matrix(F, vec_S[l].n)
            temp_row = l == 1 ? mat : kronecker_product(temp_row, mat)
        end
        H_X = j == 0 ? temp_row : vcat(H_X, temp_row)
    end
    S.cache[:H_X] = H_X
    return H_X
end

function Z_stabilizers(S::SymmetricProductCode)
    haskey(S.cache, :H_Z) && return S.cache[:H_Z]
    D = S.D
    F = S.cache[:F]
    vec_S = S.vec_S
    
    H_Z = nothing
    for j in 0:D - 1
        temp_row = nothing
        for l in 1:D^2
            mat = ((l - 1) % D == j) ? Z_stabilizers(vec_S[l]) : identity_matrix(F, vec_S[l].n)
            temp_row = l == 1 ? mat : kronecker_product(temp_row, mat)
        end
        H_Z = j == 0 ? temp_row : vcat(H_Z, temp_row)
    end
    S.cache[:H_Z] = H_Z
    return H_Z
end

function stabilizers(S::SymmetricProductCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    H_X = X_stabilizers(S)
    H_Z = Z_stabilizers(S)
    F = S.cache[:F]
    stabs = vcat(hcat(H_X, zero_matrix(F, nrows(H_X), S.n)),
                 hcat(zero_matrix(F, nrows(H_Z), S.n), H_Z))
    S.cache[:stabilizers] = stabs
    return stabs
end

# has this been extended to subsystem codes?
"""
    homological_product(S1::AbstractStabilizerCode, S2::AbstractStabilizerCode, U::CTMatrixTypes = identity_matrix(S1.F, S1.n), V::CTMatrixTypes = identity_matrix(S2.F, S2.n))
    ⊠(S1::AbstractStabilizerCode, S2::AbstractStabilizerCode) = homological_product(S1, S2)

Return the single-sector homological product code of `S1` and `S2`.

# Note
- This is the single-sector homological product. Use ⊗ for the more general product.
"""
function homological_product(S1::AbstractStabilizerCode, S2::AbstractStabilizerCode,
    U::CTMatrixTypes = identity_matrix(S1.F, S1.n), V::CTMatrixTypes = identity_matrix(S2.F, S2.n))

    return homological_product(CSSTrait(typeof(S1)), CSSTrait(typeof(S2)), S1, S2, U, V)
end

function homological_product(::IsCSS, ::IsCSS, S1::AbstractStabilizerCode, S2::AbstractStabilizerCode, U::CTMatrixTypes, V::CTMatrixTypes)

    num_stabs1 = num_X_stabs(S1)
    num_stabs1 == num_Z_stabs(S1) || throw(ArgumentError("The first code didn't have the same number of X and Z stabilizers"))
    num_stabs2 = num_X_stabs(S2)
    num_stabs2 == num_Z_stabs(S2) || throw(ArgumentError("The second code didn't have the same number of X and Z stabilizers"))
    nrows(U) == ncols(U) == S1.n || throw(ArgumentError("U is the wrong size for the code S1"))
    nrows(V) == ncols(V) == S2.n || throw(ArgumentError("V is the wrong size for the code S2"))
    
    F = S1.F
    F == S2.F == base_ring(U) == base_ring(V) || throw(ArgumentError("S1, S2, U, and V should all have the same base ring"))

    # Physical qubits scale multiplicatively
    n_new = S1.n * S2.n
    
    cache = Dict{Symbol, Any}(:F => F)
    
    return HomologicalProductCode(
        S1, S2, U, V, 
        n_new, missing, missing, 1, n_new, 
        cache
    )
end

homological_product(::IsNotCSS, ::IsNotCSS, S1::AbstractStabilizerCode, S2::AbstractStabilizerCode, U, V) = throw(ArgumentError("This is only defined for CSS codes"))
homological_product(::IsNotCSS, ::IsCSS, S1::AbstractStabilizerCode, S2::AbstractStabilizerCode, U, V) = throw(ArgumentError("This is only defined for CSS codes"))
homological_product(::IsCSS, ::IsNotCSS, S1::AbstractStabilizerCode, S2::AbstractStabilizerCode, U, V) = throw(ArgumentError("This is only defined for CSS codes"))

@doc (@doc homological_product)
⊠(S1::AbstractStabilizerCode, S2::AbstractStabilizerCode) = homological_product(S1, S2)

function X_stabilizers(S::HomologicalProductCode)
    haskey(S.cache, :H_X) && return S.cache[:H_X]
    
    S1, S2 = S.S1, S.S2
    U, V = S.U, S.V
    F = S.cache[:F]
    
    num_stabs1 = num_X_stabs(S1)
    num_stabs2 = num_X_stabs(S2)
    
    X1, Z1 = X_stabilizers(S1), Z_stabilizers(S1)
    δ1 = zero_matrix(F, S1.n, S1.n)
    for i in 1:num_stabs1
        for j in 1:num_stabs1
            δ1 += U[i, j] * transpose(Z1[i:i, :]) * X1[j:j, :]
        end
    end

    X2, Z2 = X_stabilizers(S2), Z_stabilizers(S2)
    δ2 = zero_matrix(F, S2.n, S2.n)
    for i in 1:num_stabs2
        for j in 1:num_stabs2
            δ2 += V[i, j] * transpose(Z2[i:i, :]) * X2[j:j, :]
        end
    end

    i1 = identity_matrix(F, S1.n)
    i2 = identity_matrix(F, S2.n)
    
    # ∂ = δ1 ⊗ i2 - i1 ⊗ δ2
    ∂ = kronecker_product(δ1, i2) - kronecker_product(i1, δ2)
    
    S.cache[:H_X] = ∂
    return ∂
end

function Z_stabilizers(S::HomologicalProductCode)
    haskey(S.cache, :H_Z) && return S.cache[:H_Z]
    
    # For this construction, H_Z is just transpose(H_X)
    H_X = X_stabilizers(S)
    H_Z = transpose(H_X)
    
    S.cache[:H_Z] = H_Z
    return H_Z
end

function stabilizers(S::HomologicalProductCode)
    haskey(S.cache, :stabilizers) && return S.cache[:stabilizers]
    
    H_X = X_stabilizers(S)
    H_Z = Z_stabilizers(S)
    n = S.n
    F = S.cache[:F]
    
    stabs = vcat(hcat(H_X, zero_matrix(F, size(H_X, 1), n)),
                 hcat(zero_matrix(F, size(H_Z, 1), n), H_Z))
                 
    S.cache[:stabilizers] = stabs
    return stabs
end

function _rand_single_sector_boundary(n::Int, k::Int)
    num_stabs = divexact(n - k, 2)
    U = _rand_invertible_matrix(GF(2), n)
    d0 = zero_matrix(GF(2), n, n)
    d0[k + 1:k + num_stabs, k + num_stabs + 1:end] = identity_matrix(GF(2), num_stabs)
    return U * d0 * inv(U)
end

"""
   random_homological_product_code(n1::Int, k1::Int, n2::Int, k2::Int)

Return a random homological product code.

# Note
- This implements the construction in https://arxiv.org/abs/1311.0885.
"""
function random_homological_product_code(n1::Int, k1::Int, n2::Int, k2::Int)
    d1 = _rand_single_sector_boundary(n1, k1)
    d2 = _rand_single_sector_boundary(n2, k2)
    i1 = identity_matrix(GF(2), nrows(d1))
    i2 = identity_matrix(GF(2), nrows(d2))
    d = d1 ⊗ i2 - i1 ⊗ d2
    return CSSCode(d, transpose(d))
end
