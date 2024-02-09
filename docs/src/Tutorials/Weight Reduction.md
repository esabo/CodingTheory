# Weight Reduction
Weight reduction was first introduced for CSS codes in [hastings2016](@cite), [hastings2021quantum](@cite) and for classical codes in [hastings2021fiber](@cite). Here, we follow the finite-size analysis of [sabo2024weight](@cite). The arguments of the functions below are aligned with the terminology introduced in that paper.

## Classical Codes
Weight reduction applied to classical codes acts on parity-check matrices. To weight reduce a generator matrix instead, apply weight reduction to the dual code.
```
julia> C = ReedMullerCode(1, 3)
[8, 4, 4]_2 Reed-Muller code RM(1, 3)
Generator matrix: 4 × 8
        1 1 1 1 1 1 1 1
        0 1 0 1 0 1 0 1
        0 0 1 1 0 0 1 1
        0 0 0 0 1 1 1 1

julia> parity_check_matrix(C)
[1   1   1   1   1   1   1   1]
[0   1   0   1   0   1   0   1]
[0   0   1   1   0   0   1   1]
[0   0   0   0   1   1   1   1]

julia> C_wtred = weight_reduction(C)
[27, 4, 9]_2 linear code
Generator matrix: 4 × 27
        1 1 1 1 0 0 0 0 0 0 0 1 1 1 0 0 0 1 1 1 0 1 1 1 0 0 0
        1 1 0 0 1 1 0 0 0 0 0 1 1 0 1 0 0 0 0 0 1 0 0 0 1 0 0
        0 1 0 1 0 1 0 1 1 1 1 1 1 0 0 0 1 0 1 0 1 1 1 0 0 1 0
        1 1 0 0 0 0 1 1 1 1 1 1 0 0 1 1 0 0 0 1 0 0 1 0 0 0 1

julia> parity_check_matrix(C_wtred)
[0   1   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   1   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   1   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0]
[1   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   1   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0]
[0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0]
[0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0]
[0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0]
[0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0]
[0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0]
[0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0]
[0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0]
[0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0]
[0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0]
[0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0]
[0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0]
[0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1]
[0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1]
[0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
```
As described in [sabo2024weight](@cite), this function applies independent row and column permutations by default. These may be independently turned off using the optional arguments `permute_rows` and `permute_columns`, respectively.
```
julia> C_wtred = weight_reduction(C, permute_rows = false, permute_columns = false);

julia> parity_check_matrix(C_wtred)
[1   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   1   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   1   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   1   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   1   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0]
[0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0]
[0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0]
[0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0]
[0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0]
[0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0]
[0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0]
[0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0]
[0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0]
[0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0]
[0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1]
[0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1]
[0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
```
Use the optional arguments `rows = false` or `columns = false` to reduce only the columns or rows, respectively. Provide a vector of row or column indices to the optional arguments `row_indices` and `column_indices` to only reduce specific rows or columns, respectively. If the optional arguments `row_target` or `column_target` are set, then all rows and columns with weights greater than these values are weight reduced. Compressed weight reduction is available by setting `compressed = true`. Finally, the optional argument `seed` sets `Random.seed!(seed)`, which allows for reproducible permutations.

Weight reduction may also be applied to matrices directly without having to construct a code object. This may be used to reduce a generator matrix, if desired.
```
julia> H1 = matrix(GF(2), 2, 6, [1 1 1 1 0 0; 0 0 1 1 1 1]);

julia> H2 = H1[:, [1, 5, 3, 6, 4, 2 ]];

julia> weight_reduction(H1, permute_rows = false, permute_columns = false)
[1   0   0   0   0   0   1   0   0   0   0   0]
[0   1   0   0   0   0   1   1   0   0   0   0]
[0   0   1   0   0   0   0   1   1   0   0   0]
[0   0   0   1   0   0   0   0   1   0   0   0]
[0   0   1   0   0   0   0   0   0   1   0   0]
[0   0   0   1   0   0   0   0   0   1   1   0]
[0   0   0   0   1   0   0   0   0   0   1   1]
[0   0   0   0   0   1   0   0   0   0   0   1]

julia> weight_reduction(H2, permute_rows = false, permute_columns = false)
[1   0   0   0   0   0   1   0   0   0   0   0]
[0   0   1   0   0   0   1   1   0   0   0   0]
[0   0   0   0   1   0   0   1   1   0   0   0]
[0   0   0   0   0   1   0   0   1   0   0   0]
[0   1   0   0   0   0   0   0   0   1   0   0]
[0   0   1   0   0   0   0   0   0   1   1   0]
[0   0   0   1   0   0   0   0   0   0   1   1]
[0   0   0   0   1   0   0   0   0   0   0   1]
```

The easiest way to see the effect of the permutation `H2` of `H1` is to create code objects for the matrices. Since we have already applied the desired permutation, we will turn further permutations off. Since these codes are small, the `LinearCode` constructor will automatically compute their minimum distance. (This is Example 10 of [sabo2024weight](@cite).)
```
julia> C1 = LinearCode(H1, true);

julia> weight_reduction(C1, permute_rows = false, permute_columns = false)
[12, 4, 3]_2 linear code
Generator matrix: 4 × 12
        1 1 0 0 0 0 1 0 0 0 0 0
        0 0 1 1 0 0 0 0 1 1 0 0
        0 1 0 1 1 0 0 1 1 0 1 0
        0 0 0 0 1 1 0 0 0 0 0 1

julia> C2 = LinearCode(H2, true);

julia> C2_wtred = weight_reduction(C2, permute_rows = false, permute_columns = false)
[12, 4, 4]_2 linear code
Generator matrix: 4 × 12
        1 0 0 0 0 1 1 1 1 0 0 0
        1 1 1 0 0 0 1 0 0 1 0 0
        1 0 1 1 0 0 1 0 0 0 1 0
        1 0 0 1 1 0 1 1 0 0 0 1
```

We can check the parameters with a function like
```
function check_weights(C)
    w = maximum(count(!iszero, parity_check_matrix(C)[i, :]) for i in 1:nrows(parity_check_matrix(C)))
    q = maximum(count(!iszero, parity_check_matrix(C)[:, i]) for i in 1:ncols(parity_check_matrix(C)))
    @show (w, q)
    return nothing
end

julia> check_weights(C2_wtred)
(w, q) = (3, 2)
```

## Quantum Codes
Quantum weight reduction consists of four steps: copying, gauging, thickening and choosing heights, and coning. In addition to running the entire process on a pair of stabilizer matrices or code object, each step may be run individually.

### Coning
Example 1 of [sabo2024weight](@cite)
```
julia> F = GF(2);

julia> H_X = matrix(F, 4, 6, [
           1 1 1 0 0 0;
           1 1 0 0 1 1;
           1 0 1 1 1 0;
           1 0 0 0 0 1]);

julia> H_Z = matrix(F, 1, 6, [1 0 1 0 0 1]);

julia> tilde_H_X, tilde_H_Z = copying(H_X, H_Z);

julia> tilde_H_X
[1   0   0   0   1   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   1   0   0   0   1   0   0   0   0   0   0   0   0   0   0   1   0   0   0   1   0   0   0]
[0   0   1   0   0   0   0   0   0   1   0   0   1   0   0   0   0   1   0   0   0   0   0   0]
[0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0]
[1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0]
[0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0]
[0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0]
[0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1]

julia> tilde_H_Z
[1   1   1   1   0   0   0   0   1   1   1   1   0   0   0   0   0   0   0   0   1   1   1   1]
```

All of the examples in this section will also work using a code object.
```
julia> S = CSSCode(H_X, H_Z);

julia> S_copy = copying(S)
[[24, 1]]_2 CSS stabilizer code
X-stabilizer matrix: 22 × 24
         chi(0) 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
         chi(0) 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0
         chi(0) 0 0 1 0 0 0 0 0 0 1 0 0 1 0 0 0 0 1 0 0 0 0 0 0
         chi(0) 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0
         chi(0) 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
         chi(0) 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
         chi(0) 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
         chi(0) 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
         chi(0) 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
         chi(0) 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
         chi(0) 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
         chi(0) 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0
         chi(0) 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0
         chi(0) 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0
         chi(0) 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0
         chi(0) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
         chi(0) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0
         chi(0) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0
         chi(0) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0
         chi(0) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0
         chi(0) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0
         chi(0) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 
Z-stabilizer matrix: 1 × 24
         chi(0) 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1
```

We can check the parameters with a function like
```
function check_weights(S)
    w_X = maximum(count(!iszero, X_stabilizers(S)[i, :]) for i in 1:nrows(X_stabilizers(S)))
    q_X = maximum(count(!iszero, X_stabilizers(S)[:, i]) for i in 1:ncols(X_stabilizers(S)))
    w_Z = maximum(count(!iszero, Z_stabilizers(S)[i, :]) for i in 1:nrows(Z_stabilizers(S)))
    q_Z = maximum(count(!iszero, Z_stabilizers(S)[:, i]) for i in 1:ncols(Z_stabilizers(S)))
    @show (w_X, q_X, w_Z, q_Z)
    return nothing
end

julia> check_weights(S_copy)
(w_X, q_X, w_Z, q_Z) = (4, 3, 12, 1)
```

### Gauging
Example 2 of [sabo2024weight](@cite)
```
julia> S = Q15RM();

julia> H_X = X_stabilizers(S)[[2, 1], :];

julia> H_Z = Z_stabilizers(S)[[4, 3, 2, 1], :];

julia> tilde_H_X, tilde_H_Z = gauging(H_X, H_Z);

julia> tilde_H_X
[0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   1   1   0   0   0   0   0]
[0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   1   0   0   0   0   0]
[1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0]
[0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0]
[0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0]
[0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0]
[0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   1   1]
[0   0   0   0   0   0   0   0   0   0   0   0   1   0   1   0   0   0   0   0   0   0   0   0   1]

julia> tilde_H_Z
[0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   0   0   0   1   0   0   0   0   1   0]
[0   0   0   1   1   1   1   0   0   0   0   1   1   1   1   0   1   0   0   0   0   1   0   0   0]
[0   1   1   0   0   1   1   0   0   1   1   0   0   1   1   0   1   0   1   0   1   1   0   0   1]
[1   0   1   0   1   0   1   0   1   0   1   0   1   0   1   1   1   0   0   1   0   1   0   1   0]
```

### Thickening And Choosing Heights
For thickening and choosing heights, one must specify the thickening parameter `l` and `heights`. This is Example 3 of [sabo2024weight](@cite).
```
julia> F = GF(2);

julia> l = 3; heights = [1, 2];

julia> H_X = matrix(F, 1, 4, [1 1 1 1]);

julia> H_Z = matrix(F, 2, 4, [1 1 0 0; 1 0 1 0]);

julia> tilde_H_X, tilde_H_Z  = thickening_and_choose_heights(H_X, H_Z, l, heights);

julia> tilde_H_X
[1   0   0   1   0   0   1   0   0   1   0   0   1   0]
[0   1   0   0   1   0   0   1   0   0   1   0   1   1]
[0   0   1   0   0   1   0   0   1   0   0   1   0   1]

julia> tilde_H_Z
[1   0   0   1   0   0   0   0   0   0   0   0   0   0]
[0   1   0   0   0   0   0   1   0   0   0   0   0   0]
[1   1   0   0   0   0   0   0   0   0   0   0   1   0]
[0   1   1   0   0   0   0   0   0   0   0   0   0   1]
[0   0   0   1   1   0   0   0   0   0   0   0   1   0]
[0   0   0   0   1   1   0   0   0   0   0   0   0   1]
[0   0   0   0   0   0   1   1   0   0   0   0   1   0]
[0   0   0   0   0   0   0   1   1   0   0   0   0   1]
[0   0   0   0   0   0   0   0   0   1   1   0   1   0]
[0   0   0   0   0   0   0   0   0   0   1   1   0   1]
```

### Coning
This implementation uses the Decongestion Lemma [freedman2021building](@cite) to find a cycle basis (see [sabo2024weight](@cite)). This iteratively reduces the size of the graph, and any time the graph has no cycles of length one or two, a new edge is picked at random. Different cycle bases lead to different cellulations, which leads to different stabilizers. In this way, randomness is introduced into an any prodecure which uses coning as a subroutine. As with the classical case above, an optional `seed` argument is provided to control this.

At the moment, coning is only supported for reasonable codes [hastings2021quantum](@cite) in the case the `X` stabilizers have maximum weight two overlap with the support of the `Z` stabilizer being reduced. In other words, it is designed to be input the output of thickening and choosing heights and not a random code. We may extend this later, if desired.

```
julia> H_X = matrix(GF(2), 11, 10, [
           1 1 0 0 0 0 0 0 0 0;
           0 1 1 0 0 0 0 0 0 0;
           0 0 1 1 0 0 0 0 0 0;
           0 0 0 1 1 0 0 0 0 0;
           0 0 0 0 1 1 0 0 0 0;
           0 0 0 0 0 1 1 0 0 0;
           1 0 0 0 0 0 1 0 0 0;
           0 0 0 1 0 0 0 1 0 0;
           0 0 0 0 0 0 0 1 1 0;
           0 0 0 0 0 0 0 0 1 1;
           0 0 0 0 0 0 0 1 0 1
           ]);

julia> H_Z = matrix(GF(2), 1, 10, [1 1 1 1 1 1 1 1 1 1 ]);
```
To cone, we must specify which ``Z`` stabilizers to reduce.
```
julia> tilde_H_X, tilde_H_Z = coning(H_X, H_Z, [1]);

julia> tilde_H_X
[0   0   1   1   1   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0]
[0   1   0   0   0   1   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0]
[1   0   0   0   0   0   1   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0]
[0   0   0   0   0   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0]
[1   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0]
[0   1   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0]
[0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0]
[0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0]
[0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0]
[0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0]
[0   0   0   0   0   0   1   0   0   0   0   0   0   1   0   0   0   0   0   1   0   0   0]
[0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   1   0   0   0   1   0   0]
[0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   1   0]
[0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   1]
[0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   1   0   1]

julia> tilde_H_Z
[1   0   0   0   0   0   1   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0]
[1   1   0   0   0   0   0   0   0   0   0   0   1   0   1   0   0   0   0   0   0   0   0]
[0   1   1   0   0   0   0   0   0   0   0   1   0   0   0   1   0   0   0   0   0   0   0]
[0   0   1   1   0   0   0   1   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0]
[0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0]
[0   0   0   0   1   1   0   0   0   0   0   1   0   0   0   0   0   0   1   0   0   0   0]
[0   0   0   0   0   1   1   0   0   0   0   0   1   0   0   0   0   0   0   1   0   0   0]
[0   0   0   0   0   0   0   1   1   0   1   0   0   0   0   0   0   0   0   0   1   0   0]
[0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   1   0]
[0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   1]
```

### Improved Copying
The copying variants introduced in [sabo2024weight](@cite) are available via an optional argument to `copying`.
```
julia> S = Q15RM();

julia> copying(S)
[[60, 1]]_2 CSS stabilizer code

julia> copying(S, method = :reduced)
[[32, 1]]_2 CSS stabilizer code

julia> copying(S, method = :target, target_q_X = 3)
[[16, 1]]_2 CSS stabilizer code
X-stabilizer matrix: 5 × 16
         chi(0) 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0
         chi(0) 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0
         chi(0) 0 0 0 1 1 1 1 0 0 0 0 1 1 1 0 1
         chi(0) 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 1
         chi(0) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 
Z-stabilizer matrix: 10 × 16
         chi(0) 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1
         chi(0) 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 1
         chi(0) 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 1
         chi(0) 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1
         chi(0) 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 1
         chi(0) 0 0 0 0 1 0 1 0 0 0 0 0 1 0 1 1
         chi(0) 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 1
         chi(0) 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1 1
         chi(0) 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1
         chi(0) 0 0 0 0 0 0 0 0 1 0 1 0 1 0 1 1
```

### All Together
The functions `weight_reduction` and `quantum_weight_reduction` provide a wrapper for automatically running the entire quantum weight reduction process in order. The arguments provided to each step individually are be passed into these functions with the exception that copying's optional argument `method` is now `copying_type` and `target_q_X` is now `copying_target`. The optional parameter `target_q_X` now triggers a second round of thickening and choosing heights in coning. The first, manadatory round of thickening and choosing heights is controlled via `l1` and `heights`. The second, optional round is controlled via `l2` and `target_q_X`, where the second set of heights are determined by the target.

Note the randomness of the output induced by coning.
```
julia> S = Q15RM();

julia> l = 3; heights = [2, 1, 2, 1, 2, 3, 1, 3, 3, 1];

julia> quantum_weight_reduction(S, l, heights)
[[722, 1]]_2 CSS stabilizer code

julia> quantum_weight_reduction(S, l, heights)
[[721, 1]]_2 CSS stabilizer code

julia> quantum_weight_reduction(S, l, heights, copying_type = :reduced)
[[510, 1]]_2 CSS stabilizer code

julia> quantum_weight_reduction(S, l, heights, copying_type = :target, copying_target = 3)
[[315, 1]]_2 CSS stabilizer code
```

### Copying And Gauging As Coning
It was shown in [sabo2024weight](@cite) that copying and gauging can be thought of as mapping cones.
```
julia> F = GF(2);

julia> H_X = matrix(F, 4, 6, [
           1 1 1 0 0 0;
           1 1 0 0 1 1;
           1 0 1 1 1 0;
           1 0 0 0 0 1]);

julia> H_Z = matrix(F, 1, 6, [1 0 1 0 0 1]);

julia> tilde_H_X, tilde_H_Z = copying(H_X, H_Z);

julia> tilde_H_X_cone, tilde_H_Z_cone = copying_as_coning(H_X, H_Z);

julia> tilde_H_X == tilde_H_X_cone
true

julia> tilde_H_Z == tilde_H_Z_cone
true

julia> tilde_H_X, tilde_H_Z = copying(H_X, H_Z, method = :reduced);

julia> tilde_H_X_cone, tilde_H_Z_cone = copying_as_coning(H_X, H_Z, method = :reduced);

julia> tilde_H_X == tilde_H_X_cone
true

julia> tilde_H_Z == tilde_H_Z_cone
true

julia> tilde_H_X, tilde_H_Z = copying(H_X, H_Z, method = :target, target_q_X = 3);

julia> tilde_H_X_cone, tilde_H_Z_cone = copying_as_coning(H_X, H_Z, method = :target, target_q_X = 3);

julia> tilde_H_X == tilde_H_X_cone
true

julia> tilde_H_Z == tilde_H_Z_cone
true
```
While perhaps more elegant, solving a solution of equations is more time consuming.
```
julia> using BenchmarkTools

julia> @btime copying($H_X, $H_Z);
  6.675 μs (223 allocations: 17.71 KiB)

julia> @btime copying_as_coning($H_X, $H_Z);
  83.250 μs (1131 allocations: 132.41 KiB)
```

The results are similar for gauging, although now the mapping cone is slightly faster (on this example).
```
julia> S = Q15RM();

julia> H_X = X_stabilizers(S)[[2, 1], :];

julia> H_Z = Z_stabilizers(S)[[4, 3, 2, 1], :];

julia> tilde_H_X, tilde_H_Z = gauging(H_X, H_Z);

julia> tilde_H_X_cone, tilde_H_Z_cone = gauging_as_coning(H_X, H_Z);

julia> tilde_H_X == tilde_H_X_cone
true

julia> tilde_H_Z == tilde_H_Z_cone
true

julia> @btime gauging($H_X, $H_Z);
  45.167 μs (1582 allocations: 128.55 KiB)

julia> @btime gauging_as_coning($H_X, $H_Z);
  32.084 μs (722 allocations: 66.96 KiB)
```

## Classical Versus Quantum Weight Reduction
Consider the code from the first row of Table 1 in [sabo2024weight](@cite).
```
julia> C = best_known_linear_code(6, 3)
[6, 3, 3]_2 linear code
Generator matrix: 3 × 6
        1 0 1 0 1 0
        0 1 1 0 0 1
        0 0 1 1 1 1

julia> S = HypergraphProductCode(C)
[[45, 9, 3]]_2 subsystem code

julia> quantum_weight_reduction(S, num_Z_stabs(S), collect(1:l), seed = 5849772946347113199, copying_type = :target, copying_target = 3)
[[2892, 9]]_2 CSS stabilizer code
```

Weight reducing the classical codes before passing to the hypergraph product gives.
```
julia> C_wtred = weight_reduction(C)
[9, 3, 4]_2 linear code
Generator matrix: 3 × 9
        1 0 1 0 1 0 1 0 0
        0 1 1 0 0 1 0 1 0
        0 1 0 1 1 0 0 0 1

julia> HypergraphProductCode(C_wtred)
[[117, 9, 4]]_2 subsystem code

julia> C_wtred_com = weight_reduction(C, compressed = true)
[7, 3, 3]_2 linear code
Generator matrix: 3 × 7
        1 1 1 1 0 0 0
        0 1 1 0 0 1 0
        1 0 1 0 1 0 1

julia> HypergraphProductCode(C_wtred_com)
[[65, 9, 3]]_2 subsystem code
```

## Exploring The Cycle Structure
Classical weight reduction should not change the cycle structure of the code. We can test this. Recall that a parity-check matrix defines a Tanner graph, and the girth, ``g``, of the graph is defined to the length of the shortest cycle. Short cycles are defined to be cycles with length up to ``2g - 2``. The total number of short cycles are not preserved by weight reduction since the girth may not increase as much as the length of a cycle, pushing it beyond the ``2g - 2`` limit. Elementary cycles are cycles which do not pass through the same vertex twice. The total number of elementary cycles is invariant under classical weight reduction.

We will supress the plots output from the functions below.
```
julia> C = best_known_linear_code(6, 3)
[6, 3, 3]_2 linear code
Generator matrix: 3 × 6
        1 0 1 0 1 0
        0 1 1 0 0 1
        0 0 1 1 1 1

julia> L = LDPCCode(C)
[6, 3, 3]_2 irregular 4-limited LDPC code with density 0.5555555555555556.

Variable degree polynomial:
        3//10*x^2 + 2//5*x + 3//10
Check degree polynomial:
        2//5*x^3 + 3//5*x^2
Parity-check matrix: 3 × 6
        1 1 1 1 0 0
        0 1 1 0 1 0
        1 0 1 0 0 1

julia> girth(L)
4

julia> count_short_cycles(L)
(Plot{Plots.GRBackend() n=1}, Dict(4 => 4, 6 => 2))

julia> count_elementary_cycles(L)
(Plot{Plots.GRBackend() n=1}, Dict(4 => 4, 6 => 2))

julia> C_wtred = weight_reduction(C, permute_rows = false, permute_columns = false)
[9, 3, 4]_2 linear code
Generator matrix: 3 × 9
        1 1 0 0 1 1 1 0 0
        0 1 1 0 0 1 0 1 0
        0 0 1 1 1 1 0 0 1

julia> L_wtred = LDPCCode(C_wtred)
[9, 3, 4]_2 irregular 3-limited LDPC code with density 0.2962962962962963.

Variable degree polynomial:
        3//16*x^2 + 5//8*x + 3//16
Check degree polynomial:
        3//4*x^2 + 1//4*x
Parity-check matrix: 6 × 9
        1 0 0 0 0 0 1 0 0
        0 1 0 0 0 0 1 1 0
        0 0 1 0 0 0 0 1 1
        0 0 0 1 0 0 0 0 1
        0 1 1 0 1 0 0 0 0
        1 0 1 0 0 1 0 0 0

julia> girth(L_wtred)
6

julia> count_short_cycles(L_wtred)
(Plot{Plots.GRBackend() n=1}, Dict(6 => 2, 10 => 0, 8 => 4))

julia> count_elementary_cycles(L_wtred)
(Plot{Plots.GRBackend() n=1}, Dict(6 => 2, 8 => 4))
```
We see that the girth increased, as well as the cycle lengths, but the total number of elementary cycles is still six. The function `count_short_cycles` preallocates a dictionary with entries from ``g`` to ``2g - 2``, which in this case in ten. Since there are no length ten cycles, this entry still exists but with value zero.

The hypergraph product does not preserve cycle structure, and the maximum girth of the Tanner graph is now capped at eight. Ignoring the ``X``-``Z`` correlations, let's consider the ``X`` stabilizers of the following codes.
```
julia> S = HypergraphProductCode(C)
[[45, 9, 3]]_2 subsystem code


julia> S_wtred = HypergraphProductCode(C_wtred)
[[117, 9, 4]]_2 subsystem code

julia> L_X = LDPCCode(X_stabilizers(S))
[45, 27]_2 irregular 7-limited LDPC code with density 0.1111111111111111.

Variable degree polynomial:
        2//15*x^3 + 2//5*x^2 + 4//15*x + 1//5
Check degree polynomial:
        7//90*x^6 + 4//15*x^5 + 7//18*x^4 + 4//15*x^3


julia> girth(L_X)
4

julia> L_X_wtred = LDPCCode(X_stabilizers(S_wtred))
[117, 63]_2 irregular 6-limited LDPC code with density 0.03798670465337132.

Variable degree polynomial:
        33//80*x^2 + 19//40*x + 9//80
Check degree polynomial:
        1//10*x^5 + 11//24*x^4 + 11//30*x^3 + 3//40*x^2


julia> girth(L_X_wtred)
6

julia> _, D = count_elementary_cycles(L_X)
(Plot{Plots.GRBackend() n=1}, Dict(4 => 36, 6 => 92, 10 => 176, 12 => 104, 8 => 280, 14 => 68))

julia> _, D_wtred = count_elementary_cycles(L_X_wtred); print(ans[2])
Dict(16 => 356, 20 => 108, 12 => 544, 24 => 18, 28 => 2, 8 => 352, 22 => 196, 6 => 30, 14 => 1326, 10 => 960, 18 => 748, 26 => 22)

julia> sum(values(D))
756

julia> sum(values(D_wtred))
4662
```
Even though the weight-reduced code produced more cycles after the hypergraph product, the number of shorter cycles has decreased.
```
julia> count_short_cycles(L_X)
(Plot{Plots.GRBackend() n=1}, Dict(4 => 36, 6 => 92))

julia> count_short_cycles(L_X_wtred)
(Plot{Plots.GRBackend() n=1}, Dict(6 => 30, 10 => 960, 8 => 352))
```

The cycle structure is not preserved by quantum weight reduction [sabo2024weight](@cite).
```
julia> S_qwtred = weight_reduction(S, 4, rand(1:4, nrows(S.Z_stabs)))
[[2170, 9]]_2 CSS stabilizer code

julia> L_X_qwtred = LDPCCode(X_stabilizers(S_qwtred))
[2170, 1197]_2 irregular 9-limited LDPC code with density 0.002024658584234584.

Variable degree polynomial:
        7//1444*x^6 + 7//722*x^5 + 15//1444*x^4 + 34//1083*x^3 + 65//361*x^2 + 1451//2166*x + 203//2166
Check degree polynomial:
        9//722*x^8 + 2//57*x^7 + 259//4332*x^6 + 85//722*x^5 + 445//2166*x^4 + 484//1083*x^3 + 177//1444*x^2


julia> girth(L_X_qwtred)
4

julia> count_short_cycles(L_X_qwtred)
(Plot{Plots.GRBackend() n=1}, Dict(4 => 170, 6 => 300))

julia> _, D_qwtred = count_elementary_cycles(L_X_qwtred); print(D_qwtred)
Dict(78 => 28, 56 => 2894, 16 => 13428, 20 => 24654, 58 => 15616, 52 => 7728, 60 => 1062, 12 => 6768, 24 => 38458, 28 => 48938, 8 => 2990, 30 => 140804, 72 => 16, 22 => 86594, 32 => 52376, 6 => 300, 36 => 49460, 44 => 29098, 68 => 88, 14 => 27266, 74 => 260, 64 => 400, 46 => 94434, 66 => 3418, 76 => 4, 40 => 41326, 48 => 16998, 34 => 147230, 50 => 60238, 4 => 170, 54 => 32394, 70 => 1216, 10 => 8890, 18 => 53496, 26 => 118236, 38 => 142138, 42 => 125086, 62 => 7388)
```

## Lifted Products
Classical weight reduction also applies to other types of inputs, although with the current function, the row and column indices must be specified explicitly either as a vector or a range.
```
julia> F = GF(2);

julia> S, x = PolynomialRing(F, "x");

julia> l = 63;

julia> R = ResidueRing(S, x^l - 1);

julia> A = matrix(R, 7, 7,
               [x^27, 0, 0, 1, x^18, x^27, 1,
                1, x^27, 0, 0, 1, x^18, x^27,
                x^27, 1, x^27, 0, 0, 1, x^18,
                x^18, x^27, 1, x^27, 0, 0, 1,
                1, x^18, x^27, 1, x^27, 0, 0,
                0, 1, x^18, x^27, 1, x^27, 0,
                0, 0, 1, x^18, x^27, 1, x^27])
[x^27      0      0      1   x^18   x^27      1]
[   1   x^27      0      0      1   x^18   x^27]
[x^27      1   x^27      0      0      1   x^18]
[x^18   x^27      1   x^27      0      0      1]
[   1   x^18   x^27      1   x^27      0      0]
[   0      1   x^18   x^27      1   x^27      0]
[   0      0      1   x^18   x^27      1   x^27]

julia> b = R(1 + x + x^6)
x^6 + x + 1

julia> LiftedProductCode(A, b)
┌ Warning: Commutativity of A and b required but not yet enforced.
└ @ CodingTheory ~/Documents/GitHub/CodingTheory/src/Quantum/product_codes.jl:354
[[882, 48]]_2 CSS stabilizer code

julia> A_wtred = weight_reduction(A, row_indices = 1:4, column_indices = 1:4, permute_rows = false, permute_columns = false);

julia> LiftedProductCode(A_wtred, b)
┌ Warning: Commutativity of A and b required but not yet enforced.
└ @ CodingTheory ~/Documents/GitHub/CodingTheory/src/Quantum/product_codes.jl:354
[[4914, 48]]_2 CSS stabilizer code
```
