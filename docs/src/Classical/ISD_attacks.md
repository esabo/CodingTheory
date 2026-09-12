# Information-Set Decoding Attacks

Information-set decoding (ISD) is the fastest known generic attack on a random
linear code, and its cost is what sets parameters for code-based cryptography.
Each attack repeatedly guesses an information set, hoping the error is
distributed favorably with respect to it, and the variants differ in how much
work they do per iteration in exchange for a better success probability:
Prange's algorithm does the least, Lee-Brickell and Leon allow a few errors in
the information set, Stern and its DOOM variant collide partial sums across a
split, and MMT and BJMM apply representation techniques to the same collision
idea.

These functions are also useful outside cryptography. Solving the syndrome
decoding problem is exactly what a decoder for an arbitrary linear code must
do, and running an attack against a code whose distance is unknown yields an
upper bound on that distance; see
[Minimum distance](@ref classical-minimum-distance-api).

`required_ISD_iterations` estimates the number of iterations needed for a
target success probability, which is the usual way to convert an attack into a
security estimate.

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/ISD_attacks.jl"]
Private = false
```
