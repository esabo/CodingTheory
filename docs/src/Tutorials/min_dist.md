# Minimum Distance Solvers

There are numerous 
most work by repeatedly generating codewords, some randomly, some deterministically and systematically

# Background
## Classical
Recall that for linear codes the minimum distance is equal to the minimum weight. It is therefore always valid to compute every element of the code one-by-one by constructing all possible combinations of the generators of the generator matrix, computing their weights, and then returning the minimum. This takes exponential time. Vardy showed that the computation of the minimum distance is NP-Hard and the corresponding decision problem is NP-Complete \cite{vardy1997intractability}, so we can't expect a polynomial-time algorithm but we can do better than brute force. There are two main minimum-weight algorithms in the literature: the Brouwer-Zimmermann (BZ) algorithm and Leon's algorithm \cite{leon1988probabilistic}. The latter is probabilistic and returns the minimum weight with high confidence for a binary linear code.

Consider an $[n, k, d]$ code. The BZ algorithm (roughly) works by enumerating every length-$k$ vector of weight $w$ for $1 \leq w \leq k$. After all vectors of a given weight are processed, a lower bound on the distance is increased. The algorithm also keeps an upper bound equal to the smallest weight codeword it has found so far. The algorithm terminates when the lower bound meets the upper bound. Usually this implies that a codeword of weight equal to the upper bound has been found. However, this library maintains an internal system of bounds on the distance for each code. Hence, it is possible that an upper bound has been previously computed by another method which had not yet been attained by a codeword during the BZ algorithm. Codewords acheiving the distance (or bound) are always returned with the weight whenever possible; otherwise they return the weight and `missing`.

In contrast to BZ's deterministic and systematic approach, some algorithms repeatedly generate random codewords. The lowest-weight codeword generated provides an upper bound on the distance. If the sampling is done well, a large number of iterations can provide a close bound on the distance with high probabiliity. This technique cannot provide a lower bound, so there is no way to terminate early or guarantee the answer is correct without obtaining a bound by other means. In particular, this approach can be combined with the BZ aglorithm to try to generate a codeword or lower an upper bound given BZ's lower bound. Bounds are automatically updated internally, although this may not be possible if a method is terminated prematurely by the user. Bounds from a previous run will be used in any subesquent runs.

## Quantum
The situation is much different in the quantum case. Recall that the minimum distance is given by the minimum weight of a non-trivial logical operator \eqref{dQECC}. This generally has nothing to do with the minimum distance of the corresponding stabilizer code considered as a classical, additive code. Note that $\mathcal{C}_{\mathcal{P}_n}(\mathcal{S}) \backslash \mathcal{S}$ is a set-difference of size $p^{n + k} - p^{n - k}$ and not a quotient module of size $p^{2k}$. Constructing a basis for the centralizer is an easy row reduction but enumerating its elements are not as easy as before. The brute-force method always works but the concept of an information set is now more complicated since one cannot row reduce down to the identity. White and Grassl tackled this in \cite{white2006enumeration, white2006new} where they map the additive code onto a linear code in a way that the minimum of distance of the additive code may be implied from that of the linear code. This mapping increases the parameters from $n$ to $3n$ and $d$ to $2d$, dramatically increasing the overall runtime of the BZ algorithm. Further complicating the quantum case, once a minimal weight vector is detected, one must check to see if it is an element of $S$.

To see how classical intuition can be harmful here, recall that the rotated surface code of distance $d$ has many weight-two elements. The Steane code also has minimum distance three despite having all elements of weight four. These apparent inconsistencies go back to the fact that stabilizer codes are specified by their parity-check matrices but the distances are determined by the dual. In general, low-weight elements are necessary for quantum codes to perform well against certain types of errors \cite{hu2021mitigating}.

often faster to bound dx/dz and then take min

# Methods
## Classical


Note that it is sometimes cheaper to compute the weight enumerator of the dual code then use the MacWilliams identities to compute the desired distance. The method `minimum_distance(C)` attempts to automatically determine which algorithm...


use subfield subcode to bound

## Quantum

exact algorithms
lower bounds w/ Gray code
upper bounds w/ random information sets
native QDistRnd but also interface to original GAP version
graph states
