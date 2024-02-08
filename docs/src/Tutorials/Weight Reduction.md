# Weight Reduction
Weight reduction was first introduced for CSS codes in \cite{hastsings2016, hastings2023} and for classical codes in \cite{hastingsfiber}. Here, we follow the finite-size analysis of \cite{sabo2024}. The arguments of the functions below are aligned with the terminology introduced in that paper.

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
As described in \cite{sabo2024}, this function applies independent row and column permutations by default. These may be independently turned off using the optional arguments `permute_rows` and `permute_columns`, respectively.
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
Use the optional arguments `rows = false` or `columns = false` to reduce only the columns or rows, respectively. Provide a vector of row or column indices to the optional arguments `row_indices` and `column_indices` to only reduce specific rows or columns, respectively. If the optional arguments `row_target` or `column_target` are set, then all rows and columns with weights greater than these values are weight reduced. Compressed weight reduction is available by setting `compressed = true`. Finally, the optional argument `seed` sets `Random.seed!(seed)`, which allows for reproducable permutations.

Weight reduction may also be applied to matrices directly without having to construct a code object. This may be used to reduce a generator matrix, if desired.
```
julia> H1 = matrix(F, 2, 6, [1 1 1 1 0 0; 0 0 1 1 1 1]);

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

The easiest way to see the effect of the permutation `H2` of `H1` is to create code objects for the matrices. Since we have already applied the desired permutation, we will turn further permutations off. Since these codes are small, the `LinearCode` constructor will automatically compute their minimum distance. (This is Example 10 of \cite{sabo2024}.)
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

julia> weight_reduction(C2, permute_rows = false, permute_columns = false)
[12, 4, 4]_2 linear code
Generator matrix: 4 × 12
        1 0 0 0 0 1 1 1 1 0 0 0
        1 1 1 0 0 0 1 0 0 1 0 0
        1 0 1 1 0 0 1 0 0 0 1 0
        1 0 0 1 1 0 1 1 0 0 0 1
```

## Quantum Codes

# Example 1
julia> F = GF(2);

julia> H_X = matrix(F, 4, 6, [
           1 1 1 0 0 0;
           1 1 0 0 1 1;
           1 0 1 1 1 0;
           1 0 0 0 0 1
       ]);

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

# need a check parameters?

julia> S = CSSCode(H_X, H_Z);

julia> copying(S)
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


# Example 2
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

# Example 3
# this comes out wrong here
julia> l = 3; heights = [1, 2];

julia> H_X = matrix(F, 1, 4, [1 1 1 1]);

julia> H_Z = matrix(F, 2, 4, [1 1 0 0; 1 0 1 0]);

julia> tilde_H_X, tilde_H_Z  = thickening_and_choose_heights(H_X, H_Z, l, heights);

julia> tilde_H_X
[1   0   1   0   0   1   0   0   1   0   0   1   0   0]
[1   1   0   1   0   0   1   0   0   1   0   0   1   0]
[0   1   0   0   1   0   0   1   0   0   1   0   0   1]

julia> tilde_H_Z
[1   0   1   1   0   0   0   0   0   0   0   0   0   0]
[0   1   0   1   1   0   0   0   0   0   0   0   0   0]
[1   0   0   0   0   1   1   0   0   0   0   0   0   0]
[0   1   0   0   0   0   1   1   0   0   0   0   0   0]
[1   0   0   0   0   0   0   0   1   1   0   0   0   0]
[0   1   0   0   0   0   0   0   0   1   1   0   0   0]
[1   0   0   0   0   0   0   0   0   0   0   1   1   0]
[0   1   0   0   0   0   0   0   0   0   0   0   1   1]
[0   0   1   0   0   1   0   0   0   0   0   0   0   0]
[0   0   0   1   0   0   0   0   0   1   0   0   0   0]

# Improved Copying
julia> S = Q15RM();

julia> quantum_weight_reduction(S, 3, [2, 1, 2, 1, 2, 3, 1, 3, 3, 1])
[[722, 1]]_2 CSS stabilizer code

julia> quantum_weight_reduction(S, 3, [2, 1, 2, 1, 2, 3, 1, 3, 3, 1])
[[721, 1]]_2 CSS stabilizer code

# note randomness

julia> quantum_weight_reduction(S, 3, [2, 1, 2, 1, 2, 3, 1, 3, 3, 1], copying_type = :reduced)
[[514, 1]]_2 CSS stabilizer code

julia> quantum_weight_reduction(S, 3, [2, 1, 2, 1, 2, 3, 1, 3, 3, 1], copying_type = :reduced)
[[510, 1]]_2 CSS stabilizer code

julia> quantum_weight_reduction(S, 3, [2, 1, 2, 1, 2, 3, 1, 3, 3, 1], copying_type = :target, target_q_X = 3)
[[315, 1]]_2 CSS stabilizer code

# can we do these as coning maps?

# cycle counting examples

# Examples From Table 1
julia> C = best_known_linear_code(6, 3)
[6, 3, 3]_2 linear code
Generator matrix: 3 × 6
        1 0 1 0 1 0
        0 1 1 0 0 1
        0 0 1 1 1 1

julia> S = HypergraphProductCode(C)
[[45, 9, 3]]_2 subsystem code

julia> l = 6; 

julia> quantum_weight_reduction(S, l, rand(1:l, nrows(S.Z_stabs)), seed = 123)
[[3365, 9]]_2 CSS stabilizer code

julia> quantum_weight_reduction(S, l, rand(1:l, nrows(S.Z_stabs)), seed = 197)
[[3358, 9]]_2 CSS stabilizer code

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
