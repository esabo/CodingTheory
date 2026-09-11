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

`src/Quantum/bounds.jl` keeps parameter theorems separate from the witnessed
solver bounds in `min_dist_bounds.jl`. It implements the stabilizer and
subsystem Singleton bounds, quantum MDS detection, the pure quantum Hamming
bound, and additive, Fq²-linear, and pure Feng--Ma
Gilbert--Varshamov existence bounds.
Subsystem Singleton is seeded automatically only when prime-field linearity,
Fq-linearity, or purity makes its applicability rigorous.
The Hamming code-level API requires cached purity or an explicit
`assume_pure=true`; impure subsystem codes can violate sphere packing. The GV
bound is a parameter benchmark and never modifies a concrete code's certified
distance cache.

CSS constructors inherit proven lower bounds from their classical parent
codes. The restored product families use cache-backed distance metadata.
Hypergraph products retain both transposed classical codes and seed their
X/Z bounds from the corresponding primal/transposed distance pairs.

`QuantumBoundsExt` provides arbitrary-precision Shor--Laflamme/Rains,
low-generator-weight, general check-weight, and CSS split-enumerator LPs
through JuMP and Tulip. Exact `BigInt` construction and row scaling precede
the optimizer conversion. Numerical infeasibility remains explicitly
distinguished from exact certification. A `model_hook` supports experimental
family-specific lifts, including translation-orbit variables for BB-code
scripts; univariate weight enumerators alone cannot express BB invariance.
Selective exact rational/Farkas certification remains deferred.

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
