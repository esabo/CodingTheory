# CodingTheory

It is often desirable in quantum error correction to work with a set of
overcomplete stabilizers. Therefore this constructor does not simplify any
provided set of stabilizers. The dimension of the code is computed based on the
rank, and the user should not use the matrix dimension of the stabilizers to
determine such quantities. Use `isovercomplete` to determine if an
`AbstractStabilizerCode` is overcomplete.

Credit to Tommy Hofmann of the AbstractAlgebra/Nemo/Hecke packages for help with
debugging and providing the most elegant implementation now used here.

