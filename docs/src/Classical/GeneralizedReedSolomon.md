# [Generalized Reed-Solomon Codes](@id generalized-reed-solomon-api)

This page covers the generalized Reed-Solomon (GRS) codes and the families
obtained from them by taking subfield subcodes: alternant, Goppa, Srivastava,
and generalized BCH codes. The twisted Reed-Solomon codes and the rank-metric
Gabidulin codes, both evaluation-code variants, are at the end.

The cyclic presentation of a Reed-Solomon code lives with the
[cyclic codes](@ref cyclic-codes-api); `GeneralizedReedSolomonCode` accepts a
`ReedSolomonCode` and returns the equivalent evaluation presentation.

## Generalized Reed-Solomon codes

Given ``n`` distinct evaluation points ``\gamma_1, \dots, \gamma_n`` in
``\mathbb{F}_q`` and nonzero scalars ``v_1, \dots, v_n``, the dimension-``k``
GRS code is

```math
\mathrm{GRS}_k(v, \gamma) = \{(v_1 f(\gamma_1), \dots, v_n f(\gamma_n))
    : f \in \mathbb{F}_q[x], \deg f < k\}.
```

These codes are MDS, so ``d = n - k + 1`` and the distance is set at
construction without a search. The dual is again a GRS code on the same
evaluation points; its scalars are computed by Lagrange interpolation and are
returned by `dual_scalars`. Because the length cannot exceed ``q``, a GRS code
over a small field is short, which is what motivates the subfield subcodes
below.

The generator and parity-check matrices are evaluated lazily and cached.

## Alternant codes

Fix a subfield ``\mathbb{F} \subseteq \mathbb{E}``. The alternant code
``A_r(v, \gamma)`` is the subfield subcode over ``\mathbb{F}`` of the GRS code
over ``\mathbb{E}`` with redundancy ``r``, that is, the codewords of the
``\mathbb{E}``-code whose coordinates all lie in ``\mathbb{F}``. In practice
the ``r \times n`` parity-check matrix of the parent is expanded over
``\mathbb{F}``, so the length is unchanged while the dimension satisfies
``k \geq n - rm`` for ``m = [\mathbb{E} : \mathbb{F}]``. The exact dimension
requires the rank of the expanded matrix and so is computed eagerly. The
designed distance is ``r + 1``, which is stored as a lower bound; the true
distance is generally unknown.

A generalized BCH code is the alternant code with ``r = \delta - 1`` and
scalars ``v_i = \gamma_i^b`` for offset ``b``.

Since an alternant code is defined through a parent GRS code, syndromes are
computed over the extension field by mapping back to that parent.

## Goppa codes

The Goppa code ``\Gamma(L, g)`` is the alternant code determined by a
polynomial ``g`` over ``\mathbb{E}`` and a support ``L \subseteq \mathbb{E}``
containing no root of ``g``, via the parity-check entries
``L_j^{i-1} / g(L_j)``. With ``t = \deg g`` the designed distance is
``t + 1``. Over ``\mathbb{F}_2`` the constructor sharpens this using the
factorization of ``g``, which recovers the familiar ``2t + 1`` for a separable
binary Goppa polynomial.

Goppa codes are the codes underlying the McEliece cryptosystem, which is why
`RandomGoppaCode` exists: it samples a random support and a random irreducible
Goppa polynomial, the usual private key.

## Srivastava codes

Generalized Srivastava codes are the alternant codes whose parity-check blocks
are ``z_j (a_j - w_l)^{-i}`` for ``i = 1, \dots, t`` and ``l = 1, \dots, s``,
giving designed distance ``st + 1``. A Srivastava code is the case ``t = 1``.

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/GRS_alternate.jl", "Classical/Goppa.jl"]
Private = false
```

## Twisted Reed-Solomon codes

Twisted Reed-Solomon codes perturb the monomial basis of a Reed-Solomon code by
adding higher-degree terms at selected hooks, which generally produces
non-MDS codes that are not equivalent to any Reed-Solomon code. A code is
specified by the evaluation points, the twist vector ``t``, the hook vector
``h``, and the coefficient vector ``\eta``. Taking the dual is an ``O(1)``
operation that tracks the twists, hooks, and coefficients exactly.

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/TwistedReedSolomon.jl"]
Private = false
```

## Gabidulin codes

Gabidulin codes are the rank-metric analogue of Reed-Solomon codes, evaluating
linearized rather than ordinary polynomials at points that are linearly
independent over the base subfield. Note that the distances reported by the
accessors are Hamming distances, not rank distances.

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/Gabidulin.jl"]
Private = false
```
