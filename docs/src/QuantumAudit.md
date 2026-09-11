# Quantum API audit

This report records the 2026 audit of the stabilizer/subsystem core and the
recommended implementation actions taken in response.

## Representation and correctness

Quantum generator matrices are treated as **additive spaces over the prime
field**, not automatically as linear spaces over `GF(q)`. If `GF(q)` has
degree `m` over its prime field and the additive stabilizer rank is `s`, the
reported logical dimension is

```math
k = n - \frac{s}{m}.
```

Consequently `k` may be rational. Subsystem dimensions include one prime-field
dimension for each gauge pair. Rank, containment, quotient, centralizer,
pairing, membership, random generation, and exhaustive group enumeration now
follow this additive convention. Commutation over extension fields uses the
trace-symplectic form.

The `GF(q^2)` constructor expands a Hermitian self-orthogonal linear code into
a complete prime-field additive generator set. Random constructors likewise
sample a prime-field symplectic basis and accept rational dimensions compatible
with the field degree.

## Sparse formats

Quantum constructors accept Oscar `SMat` and Julia `SparseMatrixCSC` inputs.
Binary Julia sparse matrices remain integer CSC matrices; extension-field
sparse matrices use Oscar `SMat`, because Julia CSC matrices cannot represent
implicit zeros for parent-dependent finite-field element types safely.
Algorithms densify only at explicit algebra boundaries that require Oscar
kernels or finite-field coordinate expansion, then restore sparse storage.

QLDPC degree and weight statistics have native paths for both sparse backends.
Round-trip tests cover dense binary, integer CSC, Oscar `SMat`, and additive
extension-field codes.

## Stabilizer and subsystem API

The audit repaired lazy cache use, quotient argument order, subsystem
dimensions, trait dispatch, and distance getters. Gauge fixing now materializes
a complete stabilizer code instead of returning an incomplete wrapper.

The public structural API now includes normalizer and gauge-centralizer
matrices, stabilizer/gauge/logical membership, symplectic weight, minimum group
weights, purity/degeneracy checks, CSS subsystem constructors, random
constructors, and QLDPC presentation statistics. Relic exports from the old
quadratic representation were removed rather than retained as undefined
promises.

## Phase 2: interchange and persistence

`Quantum/io.jl` provides:

- `quantum_code_data` / `quantum_code_from_data`, a versioned language-neutral
  dictionary schema;
- lossless TOML persistence suitable for Python's `tomllib` and NumPy array
  reconstruction;
- optional JLD2 persistence through `JLD2Ext`;
- binary Pauli-string import/export for phase-free interchange.

The lossless schema stores finite-field elements as prime-field coordinates,
plus field metadata, code kind, character vector, sparse preference, and
certified distance cache entries. Pauli text intentionally omits phase and
cache metadata and is not presented as a lossless format.

Direct Stim export was deferred. A stabilizer generator presentation is not a
Stim preparation circuit, and silently inventing a circuit would be a semantic
error. A future adapter should either export an explicit tableau or require a
chosen synthesis algorithm.

## Phase 2: new codes from old

`Quantum/new_codes_from_old.jl` provides:

- independent direct sums for stabilizer and subsystem codes;
- puncturing, which may correctly return a subsystem code if projected
  stabilizers become noncommuting;
- shortening by restricting to additive combinations trivial on removed
  coordinates;
- local Fourier transformations (Hadamard `X`/`Z` swaps in the binary case);
- promotion of additional Pauli generators into a gauge group;
- additive-kernel implementations of stabilizer measurement (`augment`) and
  stabilizer removal (`expurgate`) without row-pair assumptions.

Operations that cannot yet transport a nonempty character vector reject the
input instead of silently discarding phases. Local Clifford operations preserve
certified distance cache entries because they preserve Pauli weight.

## Phase 3: Hamming Shor--Laflamme enumerators

`Quantum/weight_enumerators.jl` implements the Shor--Laflamme Hamming pair:
`A` counts the additive stabilizer group by symplectic Hamming weight and `B`
counts its trace-symplectic normalizer. `B` is computed from `A` with the
quantum MacWilliams transform. Exact enumeration uses an additive prime-field
basis and a single mutable working word; no operator collection is stored.

The pair and its coefficient dictionaries are cached under
`:SL_weight_enum`, `:weight_enum_A`, `:weight_enum_B`, `:weight_dist_A`, and
`:weight_dist_B`. Portable quantum-code data stores cached coefficients as
decimal strings so arbitrarily large counts survive TOML and JLD2 round trips.
The five-qubit and Steane constructors seed their known `A` enumerators.

Portable payloads bind generators, phases, cached distance certificates, and
enumerator coefficients with a SHA-256 integrity fingerprint. Loading rejects
modified or incomplete payloads rather than attaching stale certified metadata
to a different code; callers may use `restore_cache=false` to reconstruct only
the code. Pauli files are explicitly lossy, require `allow_lossy=true` through
`save_quantum_code`, and reject signed input or nonempty character vectors.
`quantum_generator_array` provides a plain integer matrix suitable for direct
PythonCall/NumPy conversion, while `write_quantum_csv` writes the same
prime-field coordinate expansion as a header-free numeric CSV. CSV is a
matrix export, not a complete code serialization.

The shared `code_matrix_array` and `write_code_csv` APIs provide the same
NumPy-friendly export for classical linear and LDPC codes. Linear codes
default to their generator matrix; LDPC codes default to their parity-check
matrix. Either representation can be selected explicitly.

`save_code(path, code; type=:auto)` is the unified export entry point. It
normalizes the requested symbol and dispatches through `Val(type)` to CSV,
TOML, JLD2, NPZ, or Pauli backends. File extensions are used when
`type=:auto`, and `:nz` is accepted as an alias for `:npz`. NPZ remains a weak
dependency implemented by `NPZExt`; loading CodingTheory alone does not load
or precompile NPZ. The older format-specific functions remain available for
compatibility.

`misc_known_codes.jl` is restored. Its constructors now use cache-backed
distance properties. Stored triangular color-code matrices remain in the JLD2
extension; the incomplete procedural color-code generators are retained only
as internal experimental routines.

## Analytic bounds

`Quantum/bounds.jl` implements non-LP parameter bounds without conflating them
with solver certificates. `quantum_Singleton_bound` supports stabilizer and
subsystem parameters (including additive rational dimensions). Subsystem code
objects apply it automatically only for prime-field, Fq-linear, or certified
pure codes because the arbitrary additive impure case is not a general
theorem. Library `r` counts prime-field gauge pairs, so the subsystem formula
uses `k + r/degree(F)`. `is_quantum_MDS` compares it with a stored exact
dressed distance.
`quantum_Hamming_bound` implements sphere packing for pure stabilizer and
subsystem codes; its code method refuses to assume purity silently.
`quantum_Gilbert_Varshamov_bound` implements additive, Fq²-linear, and pure
Feng--Ma finite existence bounds and deliberately leaves concrete-code caches
unchanged.

The arbitrary-precision LP layer implements the binary Shor--Laflamme/Rains
system and the coarse and refined low-generator-weight constraints of
[Wei et al., *Theory of Low-Weight Quantum Codes*](https://arxiv.org/abs/2601.19848),
plus the CSS split-enumerator and general stabilizer constraints of
[Wang et al., *Check-Weight-Constrained Quantum Codes: Bounds and Examples*](https://arxiv.org/abs/2601.15446).
Combinatorial coefficients are built as `BigInt`, each row is scaled exactly,
and conversion to `BigFloat` occurs only inside a Tulip precision context.
Numerical infeasibility is reported explicitly as
`:infeasible_numerical`, never as an exact certificate. `model_hook` permits
scripts to lift these relaxations with family-specific variables and
constraints.

The implementation corrects three apparent errors in Wang et al. v1: the
Krawtchouk exponent is the polynomial degree minus the summation index, the
MacWilliams sum includes its weight-zero term, and the cumulative
check-weight inequality includes both weight-zero terms. Regression tests use
coefficients larger than `10^23`.

Raw stabilizer and subsystem constructors seed their Singleton upper bound.
CSS-from-classical constructors propagate certified parent lower bounds.
User-supplied exact distances that exceed Singleton are rejected.
Product-family hooks await restoration of the disabled product-code API;
exact rational LP certification remains a selective follow-up for boundary
cases.

## Deferred recommendations

1. Add complete or signed phase-sensitive enumerators only when a concrete
   research use requires them; do not cache operator collections.
2. Add exact rational/Farkas certification for LP boundary cases without
   replacing scalable arbitrary-precision solves.
3. Remove the quarantined legacy augmentation/expurgation implementations
   after downstream callers have migrated to the new semantics.
4. Add phase transport for puncturing, shortening, and local Clifford gates.
5. Add schema migration fixtures and external Python round-trip tests.
6. Add explicit tableau/circuit adapters for Stim, QDistRnd, and panqec rather
   than format-name-only exports.
