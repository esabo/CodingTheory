acknowledge help from Tommy Hofmann


****************************
linearcode.jl
****************************
include("linearcode.jl");

# To create a linear code, need to make generator matrix. We start by creating a field for this.
F, _ = FiniteField(2, 1, "α");
# one may either create a matrix space and then embedd the matrix into it as in
M = MatrixSpace(F, 4, 7);
G = M([1 0 0 0 0 1 1;
       0 1 0 0 1 0 1;
       0 0 1 0 1 1 0;
       0 0 0 1 1 1 1]);
# or use either of the two shorthand notations
G2 = matrix(F, [1 0 0 0 0 1 1;
       0 1 0 0 1 0 1;
       0 0 1 0 1 1 0;
       0 0 0 1 1 1 1]);
G3 = matrix(F, 4, 7, [1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1]);
G == G2
G == G3

# This is the generator matrix for the [7, 4, 3] binary Hamming code in standard form
C = LinearCode(G)
# We can get the basic information
length(C)
dimension(C)
cardinality(C)
rate(C)
# since we gave it a full rank matrix, the rank should equal the dimension
rank(G) == dimension(C)
# no minimum distance calculations were performed so the object so far does not know
# this value and the current state of the variable is set to 'missing'
# functions which rely on this value will report an error in this case
minimumdistance(C)
relativedistance(C)
genus(C)
# since we know the minimum distance of this code is 3, we can go ahead and set that
# no checking is done to verify the user's input
# since setting the minimum distance is a desirable thing to be able to do, all code
# objects are mutable (and hence live on the heap)
setminimumdistance!(C, 3)
minimumdistance(C)
relativedistance(C)
genus(C)
isMDS(C)
numbercorrectableerrors(C)
# can check the generator and parity check matrices and check their defining properties
GC = generatormatrix(C)
G == GC
H = paritycheckmatrix(C)
iszero(G * H')
iszero(H * G')
# if G is full rank, the matrix is stored directly as the generator matrix of C as written
# and H is computed as the right kernel of G
# the standard form of G and H may be retrieved by passing true into the optional parameter
Gstand = generatormatrix(C, true)
Hstand = paritycheckmatrix(C, true)
# G in the example was already in standard form, so we can check the computations
Gstand == GC
Hstand == H
# finally, we may compute the dual
# in order to save computational time, this command merely flipps the roles of G and H
# in the LinearCode constructor, and no new computation is done
dual(C)

# a vector v is in the code C if it has zero syndrome
r = generatormatrix(C)[1, :]
iszero(syndrome(r, C))
# a code C1 is a subset of C2 if every row of the generator matrix of C1 is in C2
C ⊆ C
C ⊆ dual(C)
issubcode(C, dual(C))
# two codes C1 and C2 are equivalent if C1 ⊆ C2 and C2 ⊆ C1
isequivalent(C, dual(C))
# a code is self dual if it is equivalent to its dual
isselfdual(C)
# code is self orthogonal if it is a subcode of its dual
isselforthogonal(C)

# finally, one may encode a length k vector by C using
encode(C.G[:, 1], C)
# here the vector to be encoded is an element of MatrixSpace(Nemo.GF(2), 4, 1)
# allowing for proper matrix multiplication between it and the generator matrix
# one may also pass in a Vector{Int64} of the proper length and the appropriate
# matrix space will automatically be created. this may be slow if this is done
# repeatedly as this object will be built from scratch each time. the same is
# true for syndrome(). either a (1, C.k) or a (C.k, 1) dimensional vector may
# be used as the code will automatically apply the transpose if necessary
# it will not however accept objects of type adjoint such as v'
generatormatrix(C)[:, 1]
v = [1, 0, 0, 0];
encode(v, C)
v2 = [1; 0; 0; 0];
encode(v2, C)
generatormatrix(C)[1, :]
v = [1, 0, 0, 0, 0, 1, 1];
syndrome(v, C)
v = [1; 0; 0; 0; 0; 1; 1];
syndrome(v, C)
# might want to remove these transpose examples and remove from code since these have same size

# need examples, not full rank
# passing in H
# passing in G and H directly

# we can preform several standard methods of building new codes from old ones
# direct sum
D = C ⊕ C;
generatormatrix(D)
D = directsum(C, C);
generatormatrix(D)
# direct product
D = C ⊗ C;
generatormatrix(D)
D = directproduct(C, C);
generatormatrix(D)
D = tensorproduct(C, C);
generatormatrix(D)
# the most common form of extending a code is to add parity check column to right of generator matrix
D = extend(C);
length(D)
generatormatrix(D)
F3, _ = FiniteField(3, 1, "α");
# M2 = MatrixSpace(Nemo.GF(3), 2, 4); #- ex page 15
G2 = matrix(F3, [1 0 1 1; 0 1 1 -1]);
C2 = LinearCode(G2);
generatormatrix(extend(C2))
paritycheckmatrix(extend(C2))
originalgeneratormatrix(C2) == G2
# in cases such as this it is often useful to keep track of the original generator matrix which
# was used to generate the code. can also do this for the parity check matrix. in cases where this
# does not make sense, the normal generator or parity check matrix is returned
originalgeneratormatrix(D)
originalgeneratormatrix(D) == generatormatrix(D)[:, 1:end-1]
paritycheckmatrix(D)
originalparitycheckmatrix(D)
# puncturing deletes columns from a generator matrix
C2 = puncture(D, [length(D)]);
isequivalent(C, C2)
generatormatrix(C) == generatormatrix(C2)
originalgeneratormatrix(C2) == generatormatrix(C)
# to expurgate a code is to delete rows from generator matrix and then remove any potentially resulting zero columns
D = expurgate(C, [dimension(C)]);
length(D)
dimension(D)
generatormatrix(D)
originalgeneratormatrix(D) == generatormatrix(C)
# to augment is to add rows to the generator matrix
# in this example we are going to first add the all 1's vector to the code, which is common
# we note that the all 1's vector is actually already in the code so the code
# shouldn't change
# M = MatrixSpace(Nemo.GF(2), 1, 7);
v = matrix(F, [1 for i in 1:7]')
v ∈ C
D = augment(C, v);
generatormatrix(D) == generatormatrix(C)
v = matrix(F, [1 0 1 0 1 1 1])
v ∈ C
D = augment(C, v);
generatormatrix(D)
originalgeneratormatrix(D) == generatormatrix(C)
# to shorten is to expurgate then puncture
D = shorten(C, [length(C)]);
generatormatrix(D)
originalgeneratormatrix(D) == generatormatrix(C)
# M2 = MatrixSpace(Nemo.GF(2), 3, 6);
G2 = matrix(F, [1 0 0 1 1 1; 0 1 0 1 1 1; 0 0 1 1 1 1]);
C2 = LinearCode(G2);
generatormatrix(shorten(C2, [5, 6]))
# the most common form of lengthening a code it to augment the all 1's vector and then extend
D = lengthen(C);
generatormatrix(D)
originalgeneratormatrix(D) == generatormatrix(C)
# the (u | u + v)- or Plotkin construction
# M2 = MatrixSpace(Nemo.GF(2), 3, 4);
# M3 = MatrixSpace(Nemo.GF(2), 1, 4);
G2 = matrix(F, [1 0 1 0; 0 1 0 1; 0 0 1 1]);
G3 = matrix(F, [1 1 1 1]);
C2 = LinearCode(G2);
C3 = LinearCode(G3);
generatormatrix(uuplusv(C2, C3)) #- ex p. 19

SingletonBound


****************************
cyclotomic.jl
****************************
include("cyclotomic.jl")

# ord_n(q)
q = 2;
n = 15;
b = 3;
δ = 4;
ord(n, q)
# at the moment this comes back sorted so one can easily record a coset representative
coset = cyclotomiccoset(3, q, n, false)
cyclotomiccoset(3, q, n, true)
allcyclotomiccosets(q, n, true)
dualqcosets(q, n, [coset])
complementqcosets(q, n, [coset])
# pairings
qcosetpairings(q, n)
# find all in a range
qcosettable(10, 13, q)


****************************
cycliccode.jl
****************************
include("cycliccode.jl")
q = 2;
n = 15;
b = 3;
δ = 4;
cosets = definingset([i for i = b:(b + δ - 2)], q, n, false)
C = CyclicCode(q, n, cosets)
basefield(C)
splittingfield(C)
polynomialring(C)
primitiveroot(C)
offset(C)
designdistance(C)
qcosets(C)
qcosetsreps(C)
definingset(C)
generatorpolynomial(C)
paritycheckpolynomial(C)
idempotent(C)
# CyclicCode is a subtype of LinearCode so the functions of the previous section also apply here
originalgeneratormatrix(C)
G = matrix(F, [1 0 0 0 0 1 1;
       0 1 0 0 1 0 1;
       0 0 1 0 1 1 0;
       0 0 0 1 1 1 1]);
D = LinearCode(G);
CD = directsum(C, D);
generatormatrix(CD)
# the generator matrix of a cyclic code is kept in cyclic form but the standard form is of course always available
generatormatrix(C)
generatormatrix(C, true)
# note how the previous code auto detected that it was a BCH code and returned an object of type BCHCode
complement(C)
# the above is equivalent to
B = BCHCode(q, n, δ, b)
C == B
# the optional offset parameter is set to 0 by default
BCHCode(q, n, δ)
# for cyclic codes, we can
C ∩ B == C
B = BCHCode(q, n, δ - 1, b + 4)
D = C ∩ B # check later that this is == repetition code
C + B

# Reed Solomon codes are BCH codes with n = q^m - 1
# this forces the cyclotomic cosets to each have size one
q = 16;
n = 15;
b = 3;
δ = 4;
allcyclotomiccosets(q, n, true)
# as with BCH codes, it will auto detect whether or not the code is Reed Solomon and call the appropriate constructor
cosets = definingset([i for i = b:(b + δ - 2)], q, n, false);
CyclicCode(q, n, cosets)
BCHCode(q, n, δ, b)

############## this doesn't work
ReedSolomonCode(q, δ, b)

come up with better examples from the book but basis functionality is working so far and interaction with LinearCode is good
fix C ∩ complement(C) example where defining set is an extreme (Vector{Any} so probably empty)
2nd form of CyclicCode constructor


****************************
ReedMuller.jl
****************************
include("ReedMuller.jl")

# so far only binary implemented
# fields are created via F, _ = FiniteField(2, 1, "α") by default
# codes are created using the recursive form of the generator matrix
q = 2;
ReedMullergeneratormatrix(q, 0, 3)
ReedMullergeneratormatrix(q, 2, 2)
ReedMullerCode(q, 1, 2)
C = ReedMullerCode(q, 1, 3)
ReedMullerCode(q, 2, 3)
isselfdual(C)
