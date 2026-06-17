# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

# ==============================================================================
# INVARIANT THEORY: DISTANCE BOUNDS
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the Mallows-Sloane upper bound on the minimum distance of a binary self-dual code.

# Notes
* `type = :TypeII` (Doubly-even): d ≤ 4 * ⌊n / 24⌋ + 4
* `type = :TypeI` (Singly-even): d ≤ 2 * ⌊n / 8⌋ + 2
* Codes meeting this bound exactly are called "Extremal".
"""
function Mallows_Sloane_bound(n::Int; type::Symbol=:TypeII)
    if type == :TypeII
        n % 8 == 0 || throw(ArgumentError("Type II doubly-even codes only exist for lengths divisible by 8."))
        return 4 * floor(Int, n / 24) + 4
    elseif type == :TypeI
        iseven(n) || throw(ArgumentError("Type I singly-even codes only exist for even lengths."))
        return 2 * floor(Int, n / 8) + 2
    else
        throw(ArgumentError("Code type must be :TypeI or :TypeII."))
    end
end

"""
$(TYPEDSIGNATURES)

Return the Bachoc-Gaborit upper bound on the joint geometric distance `2d + s` 
for a Type I code and its shadow.
"""
function Bachoc_Gaborit_bound(n::Int)
    iseven(n) || throw(ArgumentError("Type I singly-even codes only exist for even lengths."))
    
    if mod(n, 24) == 22
        return div(n, 2) + 6
    else
        return div(n, 2) + 4
    end
end

# ==============================================================================
# GLEASON INVARIANT GENERATORS
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the fundamental Gleason polynomial generators for binary self-dual codes.
Returns the tuple `(ϕ_2, ϕ_8, ϕ_24)` as polynomials in `R = QQ[x, y]`.
"""
function gleason_generators()
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    
    # Repetition code [2, 1, 2] generator
    ϕ_2 = x^2 + y^2
    
    # Extended Hamming code [8, 4, 4] generator
    ϕ_8 = x^8 + 14 * x^4 * y^2 + y^8
    
    # Extended Golay code [24, 12, 8] generator
    ϕ_24 = x^24 + 759 * x^16 * y^8 + 2576 * x^12 * y^12 + 759 * x^8 * y^16 + y^24
    
    return ϕ_2, ϕ_8, ϕ_24
end

# ==============================================================================
# SHADOW ENUMERATOR
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Return the shadow weight enumerator of a Type I binary self-dual code.

# Notes
* Transforms the primal weight enumerator `W(x, y)` via the MacWilliams-like 
  shadow transformation: `W_S(x, y) = 1/2^(n/2) * W(x+y, i*(x-y))`.
* If the returned polynomial contains negative or fractional coefficients, 
  the primal code is mathematically impossible (a "Ghost Code").
"""
function shadow_transform(W::MPolyRingElem, n::Int)
    iseven(n) || throw(ArgumentError("Shadow transform requires even code length n."))
    
    # 1. Establish the Cyclotomic Field Q(i) to handle the imaginary unit
    K, i_im = cyclotomic_field(4, "i")
    R_K, (xk, yk) = polynomial_ring(K, ["x", "y"])
    
    # 2. Lift the base weight enumerator into the complex polynomial ring
    W_K = map_coefficients(K, W)
    
    # 3. Apply the Conway-Sloane Shadow Substitution
    # x -> x + y
    # y -> i * (x - y)
    sub_x = xk + yk
    sub_y = i_im * (xk - yk)
    W_eval = evaluate(W_K, [sub_x, sub_y])
    
    # 4. Scale by the volume factor 1 / 2^(n/2)
    volume_factor = K(1) // K(2)^div(n, 2)
    W_shadow = volume_factor * W_eval
    
    # 5. Project back to QQ[x, y] (since physical shadow weights must be rational/integer)
    R_QQ, (x_q, y_q) = polynomial_ring(QQ, ["x", "y"])
    
    W_shadow_QQ = zero(R_QQ)
    for (c, exp_vec) in zip(coefficients(W_shadow), exponent_vectors(W_shadow))
        # In Oscar, cyclotomic field elements use coeff(c, index)
        # index 0 = real part (1), index 1 = imaginary part (i)
        real_coeff = coeff(c, 0)
        imag_coeff = coeff(c, 1)
        
        # Verify the imaginary part perfectly canceled out (which it must for a valid shadow)
        iszero(imag_coeff) || @warn "Imaginary artifacts detected in shadow transform. Input may not be self-dual."
        
        # Add the real part to the projected shadow enumerator
        W_shadow_QQ += real_coeff * (x_q^exp_vec[1]) * (y_q^exp_vec[2])
    end
    
    return W_shadow_QQ
end

# ==============================================================================
# SLOANE SANITY CHECKER
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Verify if a weight enumerator `W` is a mathematically valid Type I or Type II self-dual code.

# Notes
* Checks if `W` is in the Gleason invariant ring.
* Returns `(true, gleason_coeffs)` if valid, or `(false, [])` if it is a "Ghost Code".
"""
function is_valid_self_dual_enumerator(W::MPolyRingElem, n::Int; type::Symbol=:TypeII)
    R_QQ = parent(W)
    x, y = gens(R_QQ)
    
    # 1. Define the algebraically independent generators
    # Type I uses ϕ_2 and (x^2*y^2 * (x^2 - y^2)^2) as the basis to isolate degrees cleanly.
    # Type II uses ϕ_8 and (x^4*y^4 * (x^4 - y^4)^4).
    if type == :TypeI
        iseven(n) || return false, []
        max_i = div(n, 8)
        base_A = x^2 + y^2
        base_B = (x^2 * y^2) * (x^2 - y^2)^2
        deg_step = 8
        power_A = div(n, 2)
    elseif type == :TypeII
        n % 8 == 0 || return false, []
        max_i = div(n, 24)
        base_A = x^8 + 14 * x^4 * y^4 + y^8
        base_B = (x^4 * y^4) * (x^4 - y^4)^4
        deg_step = 24
        power_A = div(n, 8)
    else
        throw(ArgumentError("type must be :TypeI or :TypeII"))
    end
    
    # 2. Build the vector space basis for degree n
    basis = MPolyRingElem[]
    for i in 0:max_i
        # For Type I: base_A^(n/2 - 4i) * base_B^i
        # For Type II: base_A^(n/8 - 3i) * base_B^i
        power_diff = (type == :TypeI) ? 4*i : 3*i
        push!(basis, (base_A^(power_A - power_diff)) * (base_B^i))
    end
    
    # 3. Extract the coefficients to form a linear system A * c = b
    # We only need to check the first (max_i + 1) coefficients of y to fully determine the system
    target_coeffs = [QQ(coeff(W, [n - j, j])) for j in 0:(length(basis) - 1) * (type == :TypeI ? 2 : 4)]
    
    A_mat = zero_matrix(QQ, length(target_coeffs), length(basis))
    for (col, basis_poly) in enumerate(basis)
        for (row, j) in enumerate(0:(length(basis) - 1) * (type == :TypeI ? 2 : 4))
            A_mat[row, col] = QQ(coeff(basis_poly, [n - j, j]))
        end
    end
    
    b_mat = matrix(QQ, length(target_coeffs), 1, target_coeffs)
    
    # 4. Solve the system
    try
        c = solve(A_mat, b_mat, side = :right)
        
        # Verify the entire polynomial matches (handles cases where the degree was right but other terms were corrupted)
        W_test = zero(R_QQ)
        for col in 1:length(basis)
            W_test += c[col, 1] * basis[col]
        end
        
        if W_test != W
            return false, []
        end
        
        # If the coefficients of the basis are not integers, the code is a ghost.
        c_vals = [c[i, 1] for i in 1:nrows(c)]
        if !all(is_integer, c_vals)
            return false, []
        end
        
        return true, c_vals
    catch
        # If the system has no solution, it's not in the invariant ring
        return false, []
    end
end

# ==============================================================================
# EXTREMAL CODE GENERATOR
# ==============================================================================

"""
$(TYPEDSIGNATURES)

Generate the exact theoretical weight enumerator of an extremal self-dual code of length `n`.

# Notes
* Mathematically forces the low-weight coefficients to 0 to maximize minimum distance.
* If the resulting polynomial contains negative coefficients, no extremal code exists at this length!
"""
function extremal_weight_enumerator(n::Int; type::Symbol=:TypeII)
    R_QQ, (x, y) = polynomial_ring(QQ, ["x", "y"])
    
    if type == :TypeI
        max_i = div(n, 8)
        base_A = x^2 + y^2
        base_B = (x^2 * y^2) * (x^2 - y^2)^2
        power_A = div(n, 2)
        d_bound = Mallows_Sloane_bound(n, type=:TypeI)
        y_step = 2
    else
        max_i = div(n, 24)
        base_A = x^8 + 14 * x^4 * y^4 + y^8
        base_B = (x^4 * y^4) * (x^4 - y^4)^4
        power_A = div(n, 8)
        d_bound = Mallows_Sloane_bound(n, type=:TypeII)
        y_step = 4
    end
    
    basis = MPolyRingElem[]
    for i in 0:max_i
        power_diff = (type == :TypeI) ? 4*i : 3*i
        push!(basis, (base_A^(power_A - power_diff)) * (base_B^i))
    end
    
    # We want to force the coefficients of y^2, y^4, ... up to y^(d-2) to be 0.
    # The coefficient of y^0 must be 1.
    num_equations = max_i + 1
    A_mat = zero_matrix(QQ, num_equations, num_equations)
    b_mat = zero_matrix(QQ, num_equations, 1)
    
    b_mat[1, 1] = QQ(1) # The [n, 0, d] code must have exactly 1 word of weight 0 (the all-zeros word)
    
    for (col, basis_poly) in enumerate(basis)
        for (row, j) in enumerate(0:y_step:(num_equations - 1) * y_step)
            A_mat[row, col] = QQ(coeff(basis_poly, [n - j, j]))
        end
    end
    
    c = solve(A_mat, b_mat, side = :right)
    
    W_extremal = zero(R_QQ)
    for col in 1:length(basis)
        W_extremal += c[col, 1] * basis[col]
    end
    
    # Check for negative coefficients. If found, the extremal code does not physically exist.
    for coeff_val in coefficients(W_extremal)
        if coeff_val < 0
            @warn "Extremal weight enumerator contains negative coefficients. No such code exists for length $n."
            break
        end
    end
    
    return W_extremal
end
