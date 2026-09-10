# Quantum-code implementation roadmap

Deferred quantum work should be split by mathematical responsibility rather
than added to `stabilizer_code.jl` or `subsystem_code.jl`.

## New codes from old

`src/Quantum/new_codes_from_old.jl` now provides additive direct sums,
puncturing, shortening, local Fourier/Hadamard transformations, and conversion
of stabilizer codes into subsystem codes by adding gauge generators.
Puncturing is allowed to return a subsystem code when projection makes the
generators noncommuting; shortening remains in the stabilizer category.

`augment` and `expurgate` now reconstruct codes through additive kernels and
the normal constructors. Augmentation implements a stabilizer measurement:
it retains the full additive subgroup centralizing the new row before imposing
that row, rather than assuming generator rows occur in adjacent pairs.

## Weight enumerators

Implemented in `Quantum/weight_enumerators.jl`: cache-backed Hamming
Shor--Laflamme `A` and `B` enumerators, the trace-symplectic MacWilliams
transform, stabilizer/normalizer/quotient distributions, additive
extension-field enumeration, and portable persistence. Enumeration streams
through one working word and stores coefficient dictionaries, never operator
lists. Complete, signed, and phase-sensitive enumerators remain deferred.

## Quantum bounds

Add a dedicated quantum-bounds file after deciding which bounds are both
mathematically applicable and computationally useful. Candidates to evaluate
include Singleton, Hamming (pure and impure qualifications), quantum
Gilbert--Varshamov, linear-programming/Rains bounds, and subsystem variants.
Keep theorem bounds separate from exact and heuristic distance routines.

## Interchange and persistence

`src/Quantum/io.jl` now defines a versioned, NumPy-friendly data schema,
portable TOML save/load, optional JLD2 save/load, and binary Pauli-string
import/export. The lossless formats preserve the base field, additive
generator coordinates, character vector, sparse preference, and certified
distance cache entries.

Future interchange work should add explicit adapters for Stim, QDistRnd, and
panqec only when their semantics can be represented without pretending that a
stabilizer presentation is a preparation circuit. Phase-aware Pauli text and
schema migration tests are also still needed.
