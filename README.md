# CodingTheory
A basic coding theory library for Julia.

The goal of this package (technically not yet a real package) is to develop a classical and quantum error-correcting codes package in as much pure Julia as possible, eliminating the need for external libraries such as GAP typically used in the backend of other libraries. The library is built around the AbstractAlgebra and Nemo packages, and many thanks to Tommy Hofmann of these packages for helping this repo get off the ground. Anyone is welcome to contribute, although the final form of any accepted code may be standardized to maintain intra-package consistency.

The quantum part of this library is designed to eventually replace the good package in MAGMA, which is not free. At the same time, the MAGMA library is limited in scope and the features of this library already surpass it. Significant work needs to be done to extend the quantum code to include more useful features to the community.

At the moment, all functions work as intended for test cases but have not been unit tested thoroughly enough to guarantee accuracy and error free usage. All results from this library should be mentally checked and any bugs reported (or fixed and pushed). This is particularly true for the quantum part where algorithms become increasingly more complicated.

As a temporary documentation, see the file testcases.jl for demonstrations on how to use some of the basic library functions.

## Current features
### Classical
- Basic linear codes
- Subcode checking
- Syndrome/encoding calculation
- Quotients of subcodes
- Duals, Hermitian duals
- Direct sums, tensor products, Schur product
- Extending, puncturing, expurgating, augmenting, shortening, lengthening (these definitions may vary from text to text)
- Plotkin-(u | u + v) construction, (u + w, v + w, u + v + w) construction, construction X, construction X3, juxtaposition
- Expanding codes over field extensions
- Cylotomic cosets, defining sets
- Cyclic, BCH, Reed Solomon codes (RS constructor temporarily broken)
- BCH Bound with Hartmann-Tzeng Bound refinement
- Complement, dual, intersection, and union of cyclic codes
- Reed Muller code (currently just binary)
- Extended binary and ternary Golay codes
- Weight enumerators, distribution, and minimum distance

### Quantum
- General and CSS codes with auto CSS detection
- Quadratic to symplectic, symplectic to quadratic functions
- Symplectic and trace inner product
- Pauli strings to numerical vector (qubit only)
- Characteristic vector/arbitrary phases for each bit (currently only +/- 1)
- CSS construction from linear codes or matrices
- General construction from Pauli strings or quadratic matrices
- Logical space construction (S⟂ \ S) (some special cases exist which make this fail due to additive only nature of QECC while this is technically a classical code)
- Logical operators - set manually and automated construction (untested - probably buggy)
- X-, Z-, and full syndrome computations
- Generate all stabilizers
- Weight enumerators (S⟂ temporarily disabled), distribution, and minimum distance
- The [[5, 1, 3]] perfect, Steane, Shor, and [[15, 1, 3]] Reed Muller codes
- Rotated surface/XZZX codes to any distance
- 4.8.8 and 6.6.6 triangular color codes for odd distances 3 to 21 (manual arrays - needs replacing)

## TODO:
- Hamming codes in miscknowncodes.jl or Hamming.jl depending on complexity
- Golay codes in miscknowncodes.jl are actually the extended Golay codes
- Special case of Schur product for Reed Muller and cyclic codes
- BCH codes are returning in the extension field and should recast the polynomial down to the subfield
- Remove dependence of Primes package in cycliccode.jl - try to just import single function
- Include more useful known quantum codes
- Gottesman algorithm for logicals as an alternative to the current solution
- Auto calculate the triangular color code lattices for any distance
- Disjointness
- Symplectic decomposition of gates
- Eventually replace the rref used here with an additive only version to replace the multiplication by omega hack used in computing the logical operators to avoid special case errors
- Whole library is based on row operations. Switch over to the transpose since Julia is column major for potential speedup.
- Use @views in more places
- Remove manual optimizations which might hinder the compiler, such as not doing something if the result is zero.
    Basic principals to use here:
    1. Branches are expensive, adds are not
    2. If assignments are unconditional, the compiler and processor microcode can optimize them
    3. If it thinks the old value may still be needed then it will have to keep it around in a register
- Make Julia package
- Finish documentation
- Do unit testing
