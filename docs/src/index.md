# ErrorCorrection.jl

Welcome to `ErrorCorrection.jl`, a coding theory library for Julia. The package is built on the [Oscar framework](https://docs.oscar-system.org/dev/), while using as much native Julia as possible. The library supports classical,  modern (LDPC), and quantum coding theory.

The main developers so far are
* [esabo](https://github.com/esabo)
* [benide](https://github.com/benide)
We are also grateful for contributions from
* [MikeVasmer](https://github.com/MikeVasmer)
* [kalmarek](https://github.com/kalmarek)
and various members of the Oscar project, especially [thofma](https://github.com/thofma) for helping get the initial package off the ground.

If you are interested in contributing to the library, see the [Developer Documentation](link) and feel free to contact us on the #codingtheory channel of the Julia Slack workspace.

## Structure

The following constants are refernced throughout this documentation
```
const CTFieldTypes = FinField
const CTFieldElem = FinFieldElem
const CTMatrixTypes = MatElem{<:CTFieldElem}
const CTPolyRing = PolyRing{<:CTFieldElem}
const CTPolyRingElem = PolyRingElem{<:CTFieldElem}
const CTGroupAlgebra = AlgGrpElem{fpFieldElem, AlgGrp{fpFieldElem, GrpAbFinGen, GrpAbFinGenElem}}
const CTChainComplex = Union{ComplexOfMorphisms{AbstractAlgebra.FPModule{fpFieldElem}}}
```
A code is defined by matrices of type `CTMatrixTypes`, which include `fpMatrix` and `fqPolyRepMatrix`. The former have base ring `GF(p)` and the latter `GF(p, l, :α)`. Due to the way finite fields are typically represented in a computer, matrices over `GF(p, l, :α)` are *considerably* larger and slower than those over `GF(p)`. This in turn *considerably* limits the length of the codes able to be handled by the library when using this type. Therefore, it is *strongly* encouraged to utilize the field constructor `GF(p)` instead of `GF(p, 1)` when building codes over $$\mathbb{F}_p$$.

It is recommended to avoid `deepcopy(C)` and instead use `copy(C)` to create a copy of the code `C`. The use of `deepcopy` on a code object will create a new Galois field object in the struct while the matrices in the struct will still be defined over the previous Galois field. Although these two Galois fields are mathematically identical, functions in the Oscar framework consider them to be different.

The various code families in the library are mathematically related in complex patterns which are unable to be faithfully represented in Julia's linear type hierarchy. As such, it is not recommended to rely on `typeof` to discern properties of codes. For quantum codes where this is more useful, the traits `HasLogicals/HasNoLogicals`, `HasGauges/HasNoGauges`, and `IsCSS/IsNotCSS` have been setup to detect graph states ($k = 0$), subsystem codes, and CSS codes, respectively.

It is often desirable to build a code with a specific matrix representation. While properties such as standard forms and correct parameters are computed and used throughout in the background, the original matrix (matrices) used to create the code is always kept and displayed. This is of particular importance in LDPC codes, where one wants a specific representation of the code, and quantum codes, where one often prefers an over complete set of stabilizers. The user should not use matrix dimensions to determine code parameters or code parameters to iterate over matrices. Use the function `is_overcomplete` to determine if any of the matrices representing the code are over complete, i.e., have more rows than its rank.

## Suppressing The Oscar Banner

The Oscar banner will display be default when calling `using ErrorCorrection`. This can be suppressed by running Julia with the `-q` flag: `julia -q`. Note that this will also suppress the Julia banner.
