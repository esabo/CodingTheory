# ****************************
linearcode.jl
# ****************************
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
# in cases such as this it is often useful to keep track of the original generator matrix which
# was used to generate the code. can also do this for the parity check matrix. in cases where this
# does not make sense, the normal generator or parity check matrix is returned
paritycheckmatrix(D)
# puncturing deletes columns from a generator matrix
C2 = puncture(D, [length(D)]);
isequivalent(C, C2)
generatormatrix(C) == generatormatrix(C2)
# to expurgate a code is to delete rows from generator matrix and then remove any potentially resulting zero columns
D = expurgate(C, [dimension(C)]);
length(D)
dimension(D)
generatormatrix(D)
# to augment is to add rows to the generator matrix
# in this example we are going to first add the all 1's vector to the code, which is common
# we note that the all 1's vector is actually already in the code so the code
# shouldn't change
# M = MatrixSpace(Nemo.GF(2), 1, 7);
v = matrix(F, [1 for i = 1:7]')
v ∈ C
D = augment(C, v);
generatormatrix(D) == generatormatrix(C)
v = matrix(F, [1 0 1 0 1 1 1])
v ∈ C
D = augment(C, v);
generatormatrix(D)
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
# the (u | u + v)- or Plotkin construction
# M2 = MatrixSpace(Nemo.GF(2), 3, 4);
# M3 = MatrixSpace(Nemo.GF(2), 1, 4);
G2 = matrix(F, [1 0 1 0; 0 1 0 1; 0 0 1 1]);
G3 = matrix(F, [1 1 1 1]);
C2 = LinearCode(G2);
C3 = LinearCode(G3);
generatormatrix(uuplusv(C2, C3)) #- ex p. 19



# ****************************
cyclotomic.jl
# ****************************
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


# ****************************
cycliccode.jl
# ****************************
G = matrix(
    F,
    [
        1 0 0 0 0 1 1;
        0 1 0 0 1 0 1;
        0 0 1 0 1 1 0;
        0 0 0 1 1 1 1
    ],
);
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

# Remove?
# We will only be concerned with cyclic Reed-Solomon codes, but the more general, and original, definition of Reed-Solomon codes will lead us into the final family of codes we will use in this work. Let $\mathcal{P}_k(x)$ denote the set of polynomials of degree less than $k$ in $\mathbb{F}_{p^m}[x]$. The Reed-Solomon code of length $n \leq p^m$ and dimension $k < n$ is given by
# 
# $$\mathrm{RS}_{p^m}(n, k) = \{ (f(\alpha_1), \dots, f(\alpha_n)) \mid f(x) \in \mathcal{P}_k(x)\},$$
# 
# where $\alpha_i \in \mathbb{F}_{p^m}$. The most common case $n = p^m$ is the extended code of the cyclic definition, but only the case $n = p^m -1$ is, in general, cyclic. The proof of this is direct application of the Chinese Remainder Theorem.
