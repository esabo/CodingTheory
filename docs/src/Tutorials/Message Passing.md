# Message Passing In CodingTheory
The following is intended to be a review and not an introductory lesson. Prior knowledge about all topics is assumed.

## Background

### The Decoding Problem
We will utilize common notation that capital letters denote random variables, lowercase letter represent specific instances of them, and will simplify expressions such that $\mathrm{Pr}_X(X = x)$ to $\mathrm{Pr}_X(x)$. Let $\mathcal{C}$ be a binary linear code. Let $X \in \mathcal{C}$ be a random variable chosen by the transmitter with probability $\mathrm{Pr}_X(x)$ and $Y$ be a random variable over the output space. The error channel is described by its *transition probability* $\mathrm{Pr}_{Y | X}(y \mid x)$. If we decode $y$ to $\hat{x}(y) \in \mathcal{C}$, then the probability of error is $1 - \mathrm{Pr}_{X \mid Y}(\hat{x}(y) \mid y)$. Thus, to minimize the probability of error, we should choose $\hat{x}(y)$ to maximize this probability. This gives the *maximum a posteriori* (MAP) decoding rule

$$\begin{align}
    \hat{x}^\text{MAP}(y) &= \argmax_{x \in \mathcal{C}} \mathrm{Pr}_{X \mid Y}(x \mid y)\\
    &= \argmax_{x \in \mathcal{C}} \mathrm{Pr}_{Y \mid X}(y \mid x) \frac{\mathrm{Pr}_X(x)}{\mathrm{Pr}_Y(y)}\\
    &= \argmax_{x \in \mathcal{C}} \mathrm{Pr}_{Y \mid X}(y \mid x) \mathrm{Pr}_X(x).
\end{align}$$

If all codewords are equally likely, then we get the *maximum likelihood* (ML) decoding rule

$$\hat{x}^\text{MAP}(y) = \argmax_{x \in \mathcal{C}} \mathrm{Pr}_{Y \mid X}(y \mid x) = \hat{x}^\text{ML}(y),$$

where the new terminology follows from the discussion below. The term $\mathrm{Pr}_{X \mid Y}(x \mid y)$ is the *a posteriori probability* (APP), so sometimes MAP decoders are referred to as APP decoders.

### Binary Symmetric Channel (BSC)
The probability of error, $p$, is called the *crossover probability.* Let $X_i$ be a random variable over the input $\{0, 1\}$ and $Y_i$ be a random variable over the output. The *channel transition probabilities* are

$$\begin{align*}
    \mathrm{Pr}_{Y \mid X}(1 \mid 0) &= \mathrm{Pr}_{X \mid Y}(0 \mid 1) = p,\\
    \mathrm{Pr}_{X \mid Y}(1 \mid 1) &= \mathrm{Pr}_{X \mid Y}(0 \mid 0) = 1 - p
\end{align*},$$

from which we find

$$L(x_i, y_i) = (-1)^{y_i} \log\left(\frac{1 - p}{p}\right).$$

This channel is memoryless, i.e., $\mathrm{Pr}_{X \mid Y}(Y | X) = \Pi_i \mathrm{Pr}_{X \mid Y}(y_i \mid x_i)$.

Consider decoding over the BSC using a $[n, k, d]$ linear code. The criterion for decoding the BSC is to minimize the probably of decoding error, which is equivalent to maximizing the *a posteriori* probability $\mathrm{Pr}_{X \mid Y}(x \mid y)$. Thus, the optimal solution is given by

$$\hat{v} = \argmax_v \mathrm{Pr}_{X \mid Y}(x \mid y) = \argmax_v \frac{\mathrm{Pr}_{Y \mid X}(y \mid x) \mathrm{Pr}_X(x)}{\mathrm{Pr}_Y(y)}.$$

The prior is uniform over the codewords, so $\mathrm{Pr}_X(x)$ is independent of $x$ and hence $v$, and $\mathrm{Pr}_Y(y)$ is also independent of $v$. Hence, the *maximum a posteriori* (MAP) criteria can be replaced by the *maximum-likelihood* (ML) criteria

$$\hat{v} = \argmax_v \mathrm{Pr}_{Y \mid X}(y \mid x).$$

Using the monotonicity of the $\log$ and letting $d_H(\cdot, \cdot)$ denote Hamming distance,

$$\begin{align*}
    \hat{v} &= \argmax_v \log \prod_i \mathrm{Pr}_{Y \mid X}(y_i \mid x_i)\\
    &= \argmax_v \sum_i \log \mathrm{Pr}_{Y \mid X}(y_i \mid x_i)\\
    &= \argmax \left[d_H(y, x) \log(p) + (n - d_H(y, x))\log(1 - p)\right]\\
    &= \argmax_v \left[d_H(y, x) \log \left( \frac{p}{1 - p}\right) + n \log(1 - p)\right]\\
    &= \argmin_v d_H(y, x).
\end{align*}$$

In addition to belief propagation, Gallager’s A, B, and bit-flipping algorithms can be used for the BSC. The capacity of the BSC is $C_\mathrm{Sha} = 1 - H(p)$, where $H(x) = - x \log_2 x - (1 - x) \log_2(1 - x)$ is the *binary entropy function*.

```
julia> ch = BinarySymmetricChannel(1/7)
BinarySymmetricChannel(0.14285714285714285, 0.40832722141767275)

julia> BSC(1/7)
BinarySymmetricChannel(0.14285714285714285, 0.40832722141767275)

julia> crossover_probability(ch)
0.14285714285714285

julia> capacity(ch)
0.40832722141767275
```

### Binary Additive White Gaussian Noise Channel (BAWGNC)
The binary code bits $v_i \in \{0, 1\}$ are mapped to channel inputs by $x_i = (-1)^{v_i}$. The channel transition *probability density function* (pdf) is

$$p(y_i \mid x_i) = \frac{1}{\sqrt{2 \pi} \sigma} \exp \left[-(y_i - x_i)^2/(2\sigma^2)\right],$$

where $\sigma^2$ is the variance of the zero-mean Gaussian noise sample $n_i$ such that $y_i = x_i + n_i$. This channel is memoryless. One can show that

$$\mathrm{Pr}_{X \mid Y}(x_i \mid y_i) = \left[ 1 + \exp\left(-2 y_i x_i / \sigma^2\right)\right]^{-1},$$

from which we find $L(v_i \mid y_i) = 2 y_i / \sigma^2$. In practice, one must first estimate the variance.

Similar to above, the channel decoding problem for the BAWGNC comes out to be

$$\begin{align*}
    \hat{v} &= \argmax_v \mathrm{Pr}_{Y \mid X}(y \mid x)\\
    &= \argmax_v \sum_i \log\left(\frac{1}{\sqrt{2\pi}\sigma} \exp \left[-(y_i - x_i)^2 / (2 \sigma^2)\right]\right)\\
    &= \argmin_v \sum_i (y_i - x_i)^2\\
    &= \argmin_v d_E(y, x),
\end{align*}$$

where $d_E(\cdot, \cdot)$ is the standard Euclidean distance.

The capacity of the BAWGNC is

$$C_\mathrm{Sha} = \frac{1}{2} \sum_{x_i = \pm 1} \int^\infty_{-\infty} p(y \mid x) \log_2 \left(\frac{p(y \mid x)}{p(y)}\right) dy,$$

where $p(y) = \frac{1}{2} \left[p(y \mid 1) + p(y \mid -1)\right]$.

```
julia> ch = BAWGNChannel(1/7)
BAWGNChannel(0.14285714285714285, missing)

julia> BAWGNC(1/7)
BAWGNChannel(0.14285714285714285, missing)

julia> standard_deviation(ch)
0.14285714285714285

julia> variance(ch)
0.02040816326530612

julia> capacity(ch)
ERROR: Not yet written
Stacktrace:
 [1] error(s::String)
   @ Base ./error.jl:35
 [2] capacity(Ch::BAWGNChannel)
   @ CodingTheory ~/Documents/GitHub/CodingTheory/src/LDPC/channels.jl:96
 [3] top-level scope
   @ REPL[71]:1
```

### Binary Erasure Channel (BEC)

WORDS HERE


```
julia> ch = BinaryErasureChannel(1/7)
BinaryErasureChannel(0.14285714285714285, 0.8571428571428572)

julia> BEC(1/7)
BinaryErasureChannel(0.14285714285714285, 0.8571428571428572)

julia> erasure_probability(ch)
0.14285714285714285

julia> capacity(ch)
0.8571428571428572
```

### Message Passing And Belief Propagation
The most common type of decoder for LDPC codes are called *message passing* decoders. These work directly on the Tanner graph of the code. Messages are passed from the variable nodes to the check nodes constituting a single round. This iterates until some stopping criteria is met, usually either a maximum number of iterations or some level of “certainty” is reached in the current answer. The most common and general message passing algorithm is a form of belief propagation, which is used throughout many fields. Lesser known, but also used, algorithms include Gallager A and Gallager B, which will be covered in the next section.

The messages passed from variable nodes to check nodes of the Tanner graph in belief propagation are probabilities, log-likelihoods, or otherwise some metric of belief that the bit has a certain value given the observed value. The messages passed from check to variable are beliefs that the variable node has a certain value given all the messages passed to the check node from all of the *other* variable nodes. The key idea behind this is that only *extrinsic* information is considered. In other words, the information a variable node holds about itself is never considered in the calculation at the check node; the check node makes an estimate of what it thinks the variable node should be based on all the other variable nodes connected to the check node. The only way this works is if the information about the given variable node in question does not propagate through the graph and make its way to another variable node via a second check node which then gets passed to the original check node making an estimate of the original variable node. When this happens, the probability of a bit being correct is now conditioned on other bits *and* itself, at which point the mathematical equations are no longer valid. Hence, belief propagation is only valid if the paths the messages travel in the Tanner graph unrolled in time form a tree. Once the path ceases to be a tree, i.e., a cycle was encountered, belief propagation is no longer valid. Thus, the success of the decoder is correlated to the length of the shortest cycle in the Tanner graph, which is called the *girth* of the graph. If $g$ is the girth, the first $g/2$ rounds of message passing are extrinsic (independent). In practice, one finds that message passing often works well in spite of cycles because the non-independent information is either too diluted at first to make a difference or it simply does not contradict the already formed belief. The latter of course depends on whether or not the variable node in question corresponds to an error.

One of the main ingredients in message passing is likelihoods. The difference between probabilities and likelihoods is deep and runs at the heart of interpretations of probability theory: frequentist, Bayesian, and likelihoodist. The key behind likelihoods is the *law of likelihood* which says that (the evidence) $E$ favors (hypothesis) $H_1$ over $H_2$ if and only if $L = \mathrm{Pr}(E \mid H_1) / \mathrm{Pr}(E \mid H_2) > 1$, which $L$ measuring the degree of favoring. It has been suggested to treat $L = 8$ as the threshold for declaring “fairly strong evidence” in favoring one hypothesis over another and $L = 32$ for “strong evidence”. This interpretation may be useful in assigning cutoffs in iterative decoding schemes, although we must apply this idea with caution without Bayesian priors. Belief propagation, however, is distinctively stemming from a Bayesian philosophy.

In Bayesian theory, if one learns (the evidence) $E$ with certainty and nothing else, then one should replace one’s prior degree of belief $\mathrm{Pr}(H)$ of any hypothesis with $\mathrm{Pr}(H \mid E)$. Applying Bayes’s,

$$\frac{\mathrm{Pr}(H_1 \mid E)}{\mathrm{Pr}(H_2 \mid E)} = \frac{\mathrm{Pr}(H_1)}{\mathrm{Pr}(H_2)}\frac{\mathrm{Pr}(E \mid H_1)}{\mathrm{Pr}(E \mid H_2)},$$

which says that the *a posterior* ratio for a pair of hypotheses following this belief update is their prior ratio times their likelihood ratio. This update holds for probability densities as well.

The problem with the MAP problem above is that it is global, considering all bits simultaneously. This is hard, so instead message passing chooses to solve the easier local problem of *bit-wise* MAP decoding. Let $\mathcal{C}$ be a binary linear code with random variable $X \in \mathcal{C}$. The $i$th element, $v_i$, is an element of $\mathbf{F}_2 = \{0, 1\}$, however, it is more convenient to work with the input $X_i \in \{\pm 1\}$ such that $X_i = (-1)^{v_i}$. We assume the channel is memoryless (i.i.d.) such that $\mathrm{Pr}_{Y \mid X}(y \mid x) = \prod_i \mathrm{Pr}_{Y_i \mid X_i}(y_i \mid x_i)$. Finally, let $\sim x_i$ denote every bit expect for the $i$th bit and $\mathbb{1}_{x \in \mathcal{C}}$ denote an indicator function on the code. The bit-wise MAP decoding rule is

$$\begin{align*}
\hat{x}^\text{MAP}_i(y) &= \argmax_{x_i \in \{\pm 1\}} \mathrm{Pr}_{X_i \mid Y}(x_i \mid y)\\
&= \argmax_{x_i \in \{\pm 1\}} \sum_{\sim x_i} \mathrm{Pr}_{X \mid Y}(x \mid y)\\
&= \argmax_{x_i \in \{\pm 1\}} \sum_{\sim x_i}\mathrm{Pr}_{Y \mid X}(y \mid x) \mathrm{Pr}_X(x)\\
&= \argmax_{x_i \in \{\pm 1\}} \sum_{\sim x_i} \left( \prod_j \mathrm{Pr}_{Y_j \mid X_j}(y_j \mid x_j)\right) \mathbb{1}_{x \in \mathcal{C}}.
\end{align*}$$

The last term here is a *marginal probability* and is able to be computed efficiently using factor graphs via the belief propagation algorithm. Fortunately, the Tanner graph happens to be a factor graph and so BP can be run on it directly. It is known that the bit-wise MAP problem is suboptimal, sometimes dramatically, compared to the full MAP problem.

!!! note "Example"
    For the Hamming code, $\mathbb{1}_{x \in \mathcal{C}} = \mathbb{1}_{x_1 + x_2 + x_4 + x_5 = 0} \, \mathbb{1}_{x_1 + x_3 + x_4 + x_6 = 0} \, \mathbb{1}_{x_2 + x_3 + x_4 + x_7 = 0}$.

We will need the following lemma to continue.

!!! note "Lemma"
    Consider a vector of $n$ independent binary random variables $a = [a_0, \dots, a_{n - 1}]$ such that $\mathrm{Pr}(a_i = 0) = p^{(i)}_0$ and $\mathrm{Pr}(a_i = 1) = p^{(i)}_1$. Then the probability that $a$ contains an even number of 1’s is $\displaystyle \frac{1}{2} + \frac{1}{2} \prod_{i = 0}^{n - 1} \left( 1 - 2 p^{(i)}_1\right)$ and the probability that $a$ contains an odd number of 1’s is $\displaystyle \frac{1}{2} - \frac{1}{2} \prod_{i = 0}^{n - 1} \left( 1 - 2 p^{(i)}_1\right)$.

The initial information in the Tanner graph at round 0 is at the variable nodes (the value of the received bit). We begin by computing the likelihood ratio of each value. To increase numerical stability (and speed), we shift to log-likelihood ratios $\ell = \ln L(x \mid y)$. This changes the belief update formula above from multiplication to addition. The check nodes now have log-likelihood ratios from all of its connected variable nodes. It now computes an estimate that the bit attached to the $j$th edge is correct using the information from all the edges but that one, i.e., it uses only extrinsic information. Each check node can be seen as a length $w$ random bit vector $a$, where $w$ is the weight of the parity check. In order for the check to be satisfied, this vector must have an even number of 1’s. Without loss of generality, consider the $0$th edge. The probability that this is a 0 while satisfying the check is given by the probability that the remaining edges have an even number of 1’s,

$$\mathrm{Pr}(a_0 = 0 \mid y) = \frac{1}{2} + \frac{1}{2}\prod_{i = 1}^{w - 1} (1 - 2 \mathrm{Pr}(a_i = 1 \mid y_i)).$$

Rearranging gives

$$2 \mathrm{Pr}(a_0 = 0 \mid y) - 1 = \prod_{i = 1}^{w - 1}(1 - 2 \mathrm{Pr}(a_i = 1 \mid y_i)) = \prod_{i = 1}^{w - 1}(2 \mathrm{Pr}(a_i = 0 \mid y_i) - 1).$$

A little algebra and the fact that the sum of all probabilities is 1 shows that $\mathrm{Pr}(a_i = 0 \mid y_i) = L(x_i \mid y_i) / (1 + L(x_i \mid y_i))$. Therefore, $2 \mathrm{Pr}(a_0 = 0 \mid y_i) - 1 = (L(x_i \mid y_i) - 1)/(L(x_i \mid y_i) + 1)$. Recalling that $\tanh x = \sinh x / \cosh x = (e^x - e^{-x})/(e^x + e^{-x}) = (e^{2x} - 1)/(e^{2x} + 1)$, we have that

$$2 \mathrm{Pr}(a_0 = 0 \mid y_i) - 1 = \tanh(\ell_i / 2),$$

where $\ell_i = \ln L(x_i \mid y_i)$. Therefore,

$$2 \mathrm{Pr}(a_0 = 0 \mid y) - 1 = \prod_{i = 1}^{w - 1} \tanh(\ell_i / 2)$$

Repeating for $a_0 = 1$, we get

$$L(a_0 \mid y) = \frac{\mathrm{Pr}(a_0 = 0 \mid y)}{\mathrm{Pr}(a_0 = 1 \mid y)} = \frac{1 + \prod_{i = 1}^{w - 1} \tanh(\ell_i / 2)}{1 - \prod_{i = 1}^{w - 1} \tanh(\ell_i / 2)}.$$

Sticking with the log domain, we pass the natural log of this value back to the variable node on edge 0, which then updates its belief in its own value. A check node repeats this for every edge, computing the likelihoods on the $i$th edge by using all the other edges. The message nodes update their beliefs using the messages passed from every edge (check node).

The more common way to see this equation written is by applying the same hyperbolic tangent trick to the full probability to find $\tanh(\ell / 2) = \prod^{w - 1}_{i = 1} \tanh(\ell_i / 2)$, where $\ell = \ln L(a_0 \mid y)$, in which case we find the update step

$$\ell = 2 \tanh^{-1} \left( \prod_{i = 1}^{w - 1} \tanh(\ell_i / 2)\right).$$

To summarize, let $m^{(i)}_{vc}$ denote the messages passed from the variable nodes to the check nodes at round $i$ and $m^{(i)}_{cv}$ denote similarly denote the messaged passed from check to variable. Then

$$m^{(i)}_{vc} = \begin{cases} \ln L(x_v \mid y) & i = 0,\\
\ln L(x_v \mid y) + \sum_{c^\prime \neq c} m^{(i - 1)}_{c^\prime v} & i \geq 1,
\end{cases}$$

and

$$m^{(i)}_{cv} = 2 \tanh^{-1} \left( \prod_{v^\prime \neq v} \tanh \left(m^{(i)}_{v^\prime c} / 2 \right)\right).$$

Expressed in this form, message passing is called the *sum-product* *algorithm*.

The standard stopping condition for this algorithm is obtaining a zero syndrome. This may never occur, so one can also set a max number of of iterations. Note that these criteria ignore the actual value of the likelihoods.

### Box-Plus And Min-Sum
The following is summarized from Chapter 5 of [ryan2009channel](@cite).


The check-node message as written above is numerically slow and unstable. Numerous papers have been written improving this in various directions. To improve this, first, we factor $m_{vc}$ into it's sign and magnitude, $m_{vc} = \alpha_{vc} \beta{vc}$, where $\alpha_{vc} = \mathrm{sign}(m_{vc})$ and $\beta_{vc} = |m_{vc}|$. Then,

$$\tanh(m_{vc} / 2) = \prod_{v^\prime \in \mathcal{N}(c) \setminus \{ v \}} \alpha_{v^\prime c} \prod_{v^\prime \in \mathcal{N}(c) \setminus \{ v \}} \tanh(\beta_{v^\prime c} / 2)$$

such that

$$\begin{align}
m_{cv} &= \prod_{v^\prime} \alpha_{v^\prime c} 2 \tanh^{-1} \left(\prod_{v^\prime} \tanh(\beta_{v^\prime c} / 2)\right)\\
&= \prod_{v^\prime} \alpha_{v^\prime c} 2 \tanh^{-1} \log^{-1} \log \left(\prod_{v^\prime} \tanh(\beta_{v^\prime c} / 2)\right)\\
&= \prod_{v^\prime} \alpha_{v^\prime c} 2 \tanh^{-1} \log^{-1} \sum_{v^\prime} \log \left(\tanh(\beta_{v^\prime c} / 2)\right)\\
&= \prod_{v^\prime} \alpha_{v^\prime c} \phi\left(\sum_{v^\prime} \phi(\beta_{v^\prime c})),$$

where

$$\phi(x) = - \log(\tanh(x / 2)) = \log \left( \frac{e^x + 1}{e^x - 1}\right)$$

and $\phi^{-1}(x) = \phi(x)$.

```
f = Figure()
ax = Axis(f[1, 1],
    aspect = 1,
    limits = ((0, 6), (0, 6)),
    xlabel = L"x",
    ylabel = L"\phi(x)"
)
x = range(0, 6, length = 10000);
# y = -log10.(tanh.(x / 2));
y = log.((exp.(x) .+ 1) ./ (exp.(x) .- 1));
lines!(ax, x, y, color = :blue)
lines!(ax, x, x, color = :black)
f

```
![phi](./../assets/images/phi.png)

Trying to plot this function with the $\tanh$ versus plotting this function with the exponentials makes the numerical instability of the $\tanh$ apparent. Fortunately, the exponential form makes it apparent that it is not only symmetric about $y = x$, but also that the function is dominated by smaller values of $x$:

$$\phi\left(\sum_{v^\prime} \phi(\beta_{v^\prime c})\right) \approx \phi\left( \phi \left(\min_{v^\prime} \beta_{v^\prime c}\right)\right) = \min_{v^\prime} \beta_{v^\prime c}.$$

This is called the *min-sum* approximation. It is faster, more hardware friendly, but often not a very good approximation. To compensate, a normalization constant called the *attenuation constant* is often included such that

$$m_{cv} = \prod_{v^\prime} \alpha_{v^\prime c} c_\mathrm{atten} \min_{v^\prime} \beta_{v^\prime c}.$$

In practice, the attenuation is often defaulted to $1/2$ and then changed until the min-sum decoding line comes within a desired distance of the full sum-product decoding line, if possible. Note that in practice it is generally innappropriate to use min-sum unless one has checked its accuracy against sum-product. In particular, min-sum is more sensitive to finite precision such as when run on real hardware like an FPGA.

Another way to derive the min-sum approximation is through the box-plus message formula. One can show using first prinicples that the log-likelihood ratio (LLR) of the sum of two binary random variables, $L_1$ and $L_2$, respectively, gives

$$ L_1 \boxplus L_2 = \log \left(\frac{1 + e^{L_1 + L_2}}{e^{L_1} + e^{L_2}}),$$

from which one can show

$$m_{cv} = \boxplus_{v^\prime} m_{v^\prime c}.$$

(There are a lot of details here; see the reference for a full explanation.) This is the *box-sum* update formula. It is equivalent to the full sum-product, but is more hardware friendly and numerically stable. Futhermore, one can show that

$$L_1 \boxplus L_2 = \mathrm{sign}(L_1) \mathrm{sign}(L_2) \left[\min\left(|L_1|, |L_2|\right) + s\left(|L_1|, |L_2|\right)\right],$$

where

$$s(x, y) = \log\left(1 + e^{-|x + y|}\right) - \log\left(1 + e^{-|x - y|}\right)$$

is the *correction term*. Ignoring the correction gives the min-sum approximation. The correction is expensive to compute; a lower complexity version is given by

$$\tilde{s} s(x, y) = \begin{cases} c && \text{if $|x + y| < 2$ and $|x - y| > 2|x + y|$}\\
-c && \text{if $|x - y| < 2$ and $|x + y| > 2|x - y|$}\\
0 && \text{otherwise} \end{cases}.$$

This is included in this library under the name `min_sum_with_correction` using $c = 0.5$.

### Syndrome-Based Belief Propagation
Let $\mathcal{C}$ be a binary linear code with parity-check matrix $H$. In the previous section on belief propagation, we viewed the problem as a codeword, $c \in \mathcal{C}$, being sent over a channel where an error, $e$, occurred, causing the vector $y = c + e$ to be received. Given $y$, the goal was to recover $c$. The error was never explicitly computed, although it could be easily inferred by subtraction. In syndrome-based belief propagation, we apply the parity-check matrix to get $s = Hy = Hc + He = He$, where we used the fact that $Hc = 0$ since $c \in \mathcal{C}$. Now, we are given the syndrome $s$ and the goal is to determine $e$. The exact codeword that was sent is irrelevant since it does not contribute to the syndrome.

The mathematical problem is to determine the most-likely error at each bit given the syndrome $s$,

$$\begin{align}
    \hat{e}^\text{MAP}_i(s) &= \argmax_{e_i \in \mathbf{F}_2} \mathrm{Pr}(e_i \mid s)\\
    &= \argmax_{e_i \in \mathbf{F}_2} \sum_{\sim e_i} \mathrm{Pr}(e \mid s)\\
    &= \argmax_{e_i \in \mathbf{F}_2} \sum_{\sim e_i} \mathbb{1}_{He = s} \mathrm{Pr}(e \mid s),
\end{align}$$

where we have used the fact that $\mathrm{Pr}(e \mid s) = 0$ if $He \neq s$. If $He = s$ and the errors are independent (memoryless), $\mathrm{Pr}(e \mid s) = \prod \mathrm{Pr}(e_i)$.

As before, belief propagation can be used to compute these marginals.
(give messages)



As written, the indicator function says we should fix an error on a variable node (edge), and if it combined with the current guess for the other variable nodes (edges) satisfy the check node, then the probability is non-zero and should be included in the message. This is repeated for all possible errors. A straightforward implementation would therefore include a `for` loop over all possible errors followed by an `if` statement matching the syndrome. However, for binary codes there is only one possible error which will work (either a 0 or a 1) and belief propagation automatically returns the desired probability. So the same implementation may be used for the regular and syndrome-based algorithms for binary codes. The difference between the two is the initialization and stopping criteria. For regular belief propagation, start with the received vector and check for the zero syndrome (a codeword). If $\mathrm{wt}(e)$ errors occurred, then it will need to flip $\mathrm{wt}(e)$ bits. For syndrome-based, start with the zero error (no error) and check for the syndrome $s$. If $\mathrm{wt}(e)$ errors occurred, then it will also need to flip $\mathrm{wt}(e)$ bits. The minus sign is necessary in the syndrome-based message formula to get the bit to flip from 0 to 1. Non-binary codes do not simplify as nicely, although quantum codes over $\mathbf{F}_{2^\ell}$ can use a similar trick since the trace inner product means the syndrome bits are elements of the prime subfield, $\mathbf{F}_2$.

### Scheduling Patterns
The order messages are passed around is referred to as a *schedule*. In the standard message passing algorithms described above, all check nodes are updated at the same time then all variable nodes are updated at the same time. This is called a *parallel* or *flooding* schedule. Here, a node uses the previous iteration's values to update the current values. If the value of a variable node is changed by a check node, the other check nodes will not know until the next iteration. In a *serial* schedule, we iterate over the check (variable) nodes and send all messages into the node and then all messages out of the node. If the value of a variable node is changed by a check node, all following check nodes use this new value immediately. In a *semi-serial* or *layered* schedule, the nodes are partitioned into disjoint sets. The first set of nodes performs one iteration in parallel, then the second set of nodes performs one iteration in parallel, and so on until all sets have performed one iteration. Each (parallel) set uses the new values of the previous (parallel) set as if it was run in serial. For this to work, the nodes in each set must not interact (be connected to) the nodes in any other set.

The advantage of a parallel schedule is that all messages are able to be computed and processed simultaneously. A serial schedule is slower because a node has to wait for all its predecessors to finish before proceeding. On the other hand, if a bit is going to flip, a serial schedule passes that information on faster than a parallel schedule, which is assumed to speed convergence. Fewer iterations is good for hardware but may still be slower than giving up parallelization. It can also be good for converging before the affect of short cycles causes the decoder to fail. A serial schedule can also be beneficial in situations where the value of a node is being pulled in opposite directions by different nodes causing no change to happen. In a serial schedule, the node will make a decision based on the first message, which could change the result of the second message and break the stalemate. A semi-serial schedule attempts to recover some of the speed of the parallel schedule while retaining some of the benefits of the serial schedule.

Parallel scheduling is default for all methods. Other scheduling may be set using the optional parameter `schedule`. No multi-threading is actually performed in these message passing implementations because it is assumed the goal is to run large scale simulations in which the threads would be better utilized in another manner.

## Belief Propagation In CodingTheory
To use the message-passing decoders coding theory, first we need a (code) parity-check matrix, a received vector, and a noise model. For simplicity, we will introduce a weight-one error.
```
using Oscar, CodingTheory
julia> H = matrix(GF(2), [1 1 0 1 1 0 0; 1 0 1 1 0 1 0; 0 1 1 1 0 0 1]);

julia> c = matrix(GF(2), 7, 1, [1, 1, 1, 0, 0, 0, 0]);

julia> e = matrix(GF(2), 7, 1, [0, 0, 1, 0, 0, 0, 0]);

julia> y = c + e;

julia> H * c
[0]
[0]
[0]

julia> syn = H * y
[0]
[1]
[1]
```

REDO NOISE MODELS TO MATCH CHNS ABOVE?

For the noise model, we will choose a BSC with the arbitrarily chosen crossover probability $p = 1/7$.
```
julia> nm = MPNoiseModel(:BSC, 1/7);
```
Finally, to decode
```
julia> flag, out, iter, = sum_product(H, y, nm)
(true, [1, 1, 1, 0, 0, 0, 0], 3)
```

If `flag` is true, then the decoder converged to the codeword `out` in `iter` iterations, at which point we declare it to be successful if `out` is equal to `c`.

This function uses the $\tanh$ formula. The box-plus update formula is available with `sum_product_box_plus`. The results should be equivalent, but one may be faster and more numerically stable. The optional parameters and their defaults are
- `max_iter::Int = 100` - sets the maximum number of iterations
- `chn_inits::Union{Missing, Vector{Float64}} = missing` - initial conditions on the bits in the form of log-likelihoods, used to pass soft information into the decoder
- `schedule::Symbol = :parallel` - choose either a `:parallel`, `:serial`, or `:layered` schedule
- `erasures::Vector{Int} = Int[]` - set the erasure locations with respect to the input vector `y`
These are also available with all of the belief-propagation-based algorithms described below.

An erasure is by definition a bit which we have no information about. For a binary code, the channel initialization on this bit is therefore $\mathrm{Pr}(0 \mid y) = \mathrm{Pr}(1 \mid y) = 1/2$ in the probability domain and 0 in the log domain. Since the symbol $e$ or $?$ are not valid finite-field elements, a value for the bit must be entered at erasures locations, although this will subsequently be ignored by the algorithm. This and the channel initialization are automatically handled by listing the erasure locations in the parameter `erasures::Vector{Int}`.

EXAMPLE HERE OF AN ERASURE, MAYBE DEMO BEYOND MIN D DECODING
```

```

For the syndrome-based version,
```
julia> flag, out, iter, = sum_product_syndrome(H, syn, nm)
(true, [0, 0, 1, 0, 0, 0, 0], 3)
```
Here, decoding is declared to be successful if `out` is equal to `e`.

The min-sum functions are similar; box-plus is the only formula implemented here.
```
julia> flag, out, iter = min_sum(H, y, nm)
(true, [1, 1, 1, 0, 0, 0, 0], 1)

julia> flag, out, iter = min_sum_syndrome(H, syn, nm)
(true, [0, 0, 1, 0, 0, 0, 0], 1)

julia> flag, out, iter = min_sum_with_correction(H, y, nm)
(true, [1, 1, 1, 0, 0, 0, 0], 13)

julia> flag, out, iter = min_sum_with_correction_syndrome(H, syn, nm)
(true, [0, 0, 1, 0, 0, 0, 0], 13)
```
The normalization constant is set to 0.5 by default and may be changed via the optional parameter `attenuation::Float64 = 0.5`.

### Gallager A & B And Bit-Flipping Decoding Algorithms
Gallager’s original work contained three message-passing algorithms. The messages are not beliefs, so they are distinct from the above work. In Gallager A, the check nodes perform a simple mod-2 sum of all the variable bits on all but a given edge and then pass back on that edge what that bit must be to satisfy the parity check. If *any* check node passes back the value received on the channel, it uses the channel value. If *all* check nodes disagree with the channel, then it uses the value that the check nodes agree on. Gallager B is almost identical but instead flips the channel bit if at least some fixed threshold $t$ of check nodes disagree with the channel bit. Gallager provided a formula for the optimal value of $t$ for regular LDPC codes. These algorithms actually perform quite well on some channels.

The bit-flipping algorithm first evaluates all parity-check equations (rows of $H$) and then flips any bits in the received word that are involved in more than some fixed number of failed parity checks. This is repeated until all the parity checks are satisfied or a fixed number of iterations has been reached.

The optional parameters and their defaults are
- `max_iter::Int = 100`
- `schedule::Symbol = :parallel`
and
- `threshold::Int = 2`
for Gallager B.

DOES NOT CONVERGE - HAD REWRITTEN FUNCTIONS, NEED TO FIX
```
julia> Gallager_A(H, y)
(true, [0, 0, 0, 0, 0, 0, 0], 2)

julia> Gallager_B(H, y)
(true, [0, 0, 0, 0, 0, 0, 0], 2)
```

## Decimation
Decimation was previously used in message-passing applications outside of error correction and was applied to stabilizer codes in [](@cite). The idea is to freeze the value of a variable node. We can either do this from the start to obtain a so-called *genie-aided* decoder, or we can periodically pause message passing to fix a bit. In *guided decimation*, we pause every fixed number of rounds and freeze the value of the variable node with the highest log-likelihood ratio. In *automated decimation*, we pause after every iteration and fix any bit whose absolute value has passed a certain threshold.

It is important to note that when a variable node is fixed to a specific value, the decoder is now sampling possible solutions with that fixed bit, which is different from the ML and MAP problems above. Furthermore, if there is a unique solution and a bit is fixed which does not match the solution, the decoder will fail instead of correcting that bit. For example, fixing the most reliable bit in guided decimation may mean fixing a bit which is still far from reliable and could go either way. On the other hand, fixing a bit could help the decoder converge faster and also break out of trapping sets. In this sense, decimation can be very helpful decoding degenerate stabilizer codes where there are many valid solutions and BP has a difficult time picking one to converge to.

## Trapping Sets

# Simulations
While the examples above are interesting learning tools, they are not particularly relevant. A more useful example is an FER simulation of a code. The way this library is written, the initialization and actual message passing is split into different functions. The goal of this is to only need to initialize a Tanner graph once per simulation. To decode multiple errors, we only need to reset the message values to either zero or the channel inputs btween runs. In this way, we can decode multiple times without allocating new memory. If the message passing function and simulation functions are written properly, the number of allocations will be independent of the number of samples taken.

The examples here will all use *direct sampling*. In general, this is highly inappropriate for research-level simulations; however, when importance sampling was included the number of lines of code tripled and became more difficult to follow. Hence for simplicity, since this is *not* a research paper, we will aim for easier to follow code. This has a number of downsides in terms of results. Most importantly, suppose we wish to demonstrate an FER of $10^{-6}$. If we sample 15k errors, the two best answers we can get are $0$ (no decoding failures) or $1/15000 = 6.67 * 10^{-5}$ (a single decoding failure). Hence, we may be dramatically over or underestimating the FER. In order to reach the desired accuracy, we'd need to direct sample at least 1 million times, which may or may not be feasible for certain codes without a large amount of computing resources. Having neither the time, desire, nor resources for this, these examples will stick to 15k samples. A black horizontal line will be drawn on plots at an FER of $1/15000$ to remind us of the artificial sampling-induced error floor none of our decoders can surpass. Lines which have an FER value of $0$ will abruptly terminate (going from right-to-left) since the $y$-axis is a plotted on a log scale. (This of course implies we need to sample more or switch to importance sampling.)

Another problem with direct sampling is the expected value of number of errors. Sampling errors on length 20 code under a BSC with crossover probability $p = 10^{-3}$ will produce on average no errors per sample. If 15k samples are taken but only a small number of them were not the all-zero error, then it's likely we have underestimated the FER. Again, since the sole goal here is to demonstrate the library functions, we will ignore these issues.

On the flip side, we don't want to waste time by oversampling. Unfortunately, the number of samples varies along the $x$-axis and is best determined after a preliminary plot has been created. We will come back to this later.

Finally, we wish to point out that the actual message-passing functions in this library are *not* parallelized. As a result, all schedules are going to take an equal amount of time. This is because it is a better use of limited resources to parallelize across a large number of errors in a simulation than to parallelize inside a function which may only take nanoseconds, which may actually slow it down.





EXAMPLE HERE COMPARING SYNDROME TO REGULAR BP
```

```

EXAMPLE HERE EXPLORING THE CONVERGENCE SPEED AND RESULTS OF THE DIFFERENT SCHEDULES
```

```


## Ensemble Decoding

CYCLE COUNTING?

# Belief Propagation Of Stabilizer Codes

# Simulations
## Example 1: Independent CSS Decoding With Bayes' Theorem
Although we haven't yet done decoding over $\mathbb{F}_4$, we are in a position to decode CSS codes using independent $X$ and $Z$ decoders. All the criticisms above about, for example, importance sampling, also hold for these simulations.

As before, we are going to start out by initializing the decoder parameters, noise model, and threading. Decoding failure here is measured by failure of any of the decoders to converge or if the error plus the correction implements a logical operator, and we want to keep track of these separately. We will also keep track of the number of iterations it took to converge, if it converges. To prevent slowdowns caused by appending iteration counts to a vector, we will allocate this ahead of time. As before, we will avoid race conditions by storing this separately for each thread and then combining all of the data at the end. 

```
function CSS_decoder_simple(S::AbstractStabilizerCode)
    # BP details
    n = length(S)
    num_runs = 15000
    max_iter = 20
    schedule = :layered
    X_layers::Vector{Vector{Int}} = layered_schedule(X_stabilizers(S))
    Z_layers::Vector{Vector{Int}} = layered_schedule(Z_stabilizers(S))
    erasures = Int[]

    # noise details
    noise = 10 .^(-2:0.1:-1);
    len_noise = length(noise)

    # multi-threading details
    num_threads = Threads.nthreads()
    runs_per_thread = cld(num_runs, num_threads)
    new_num_runs = runs_per_thread * num_threads
    println("Number of threads: $num_threads, runs per thread: $runs_per_thread, new number of runs: $new_num_runs")

    # preallocate places to store results
    FER = zeros(Float64, len_noise)
    local_log_failure_counts = zeros(Int, num_threads)
    local_conv_failure_counts = zeros(Int, num_threads)
    X_local_iters = [zeros(Int, max_iter) for _ in 1:num_threads]
    Z_local_iters = [zeros(Int, max_iter) for _ in 1:num_threads]
```

Next we are going to loop through the noise parameters and setup the noise model for the given probability of error, $p$. Since we will never write to these values, we can do it outside the threading. Here, we chose to simulate the depolarizing channel and apply Bayes' Theorem between the $X$ and $Z$ corrections. The choice to do $X$, Bayes', then $Z$ instead of $Z$, Bayes', then $X$ is arbitrary. Note that this necessarily slows down the decoder by having to wait for one to finish before starting the other instead of doing them both at the same time, but again, the way these functions are parallelized, they would not run at the same time anyway.

Each thread needs it's own local copy of the decoder and its initializations to avoid race conditions. Again note that this dramatically increases the memory footprint of the function since without the threading we would be able to do this once outside of the loop.

```
    for i in 1:len_noise
        p = noise[i]
        println("Starting p = $p")

        # setup noise model
        chn = Dict(0 => 1 - p, 1 => p /3, 2 => p /3, 3 => p / 3)
        # Pauli_types::Vector{Int} = collect(keys(chn))
        # Pauli_weights::Vector{Float64} = collect(values(chn))
        Pauli_types = Int[0, 1, 2, 3]
        Pauli_weights = Float64[1 - p, p / 3, p / 3, p / 3]
        PrII::Float64 = Pauli_weights[1] / (Pauli_weights[1] + Pauli_weights[4])
        PrIX::Float64 = Pauli_weights[2] / (Pauli_weights[1] + Pauli_weights[2])
        PrZI::Float64 = Pauli_weights[4] / (Pauli_weights[1] + Pauli_weights[4])
        PrZX::Float64 = Pauli_weights[3] / (Pauli_weights[2] + Pauli_weights[3])
        without_X_err::Float64 = log(PrII / PrZI)
        with_X_err::Float64 = log(PrIX / PrZX)
        Pauli_weights_W = Weights(Pauli_weights)
        Pauli_op::Vector{Int} = [0]

        Threads.@threads for th in 1:num_threads
            X_stabs, Z_stabs, logs, X_err, Z_err, X_var_adj_list,
                X_check_adj_list, Z_var_adj_list, Z_check_adj_list, X_chn_inits,
                Z_chn_inits, X_check_to_var_messages, X_var_to_check_messages,
                Z_check_to_var_messages, Z_var_to_check_messages, current_bits,
                totals, X_syn, Z_syn = CodingTheory._message_passing_init_CSS(S, chn, max_iter,
                missing, missing, schedule, erasures)

            nr_X = size(X_stabs, 1)
            nr_Z = size(Z_stabs, 1)
            nr_logs = size(logs, 1)
            X_meas_syn = zeros(Int, nr_X)
            Z_meas_syn = zeros(Int, nr_Z)
            skip = false
```

There are now four types of errors to sample: Pauli $I$, $X$, $Y$, and $Z$. Since we are direct sampling i.i.d. errors, we must do this for each qubit. Here we use the `StatsBase` function `sample!` to properly sample with respect to the given noise model. The resulting Pauli operator is stored in `Pauli_op[1]`, which is setup earlier to reduce allocations. Since we are using a CSS decoding scheme, the errors are separated in two length `n` vectors `X_err` and `Z_err` instead of a single length `2n` vector. Note that these vectors are allocated once before this loop and therefore require setting no error locations to zero in case there was an error there in the previous run.

Recall that $Z$ errors are detected by the $X$ stabilizers and vice versa. Earlier we allocated `X_meas_syn` of the proper size, and now we utilize `mul!` to store the result in this memory. Since BP checks if the syndrome of the current variable node estimates are equal to the passed in syndrome, we need to reduce the multiplication modulo two to match. Finally, we decode the $X$ syndrome. If this fails, we are going to skip to the end and declare that error as a failure since the result will never end up back in the codespace regardless of the outcome for $Z$. If the $X$ decoder converges, we are going to make a note of the number of iterations it took and apply the correction back to the current state. We proceed similarly for the $Z$ syndrome after improving the channel initialization with Bayes' Theorem. As before, we could add a check to decode if and only if the error vectors are non-zero to not only save time but also to not overestimate the results; however, in the range this becomes important, one should reconsider if direct sampling is appropriate in the first place.

```
            for _ in 1:runs_per_thread
                # sample
                @inbounds @simd for j in 1:n
                    sample!(Pauli_types, Pauli_weights, Pauli_op)
                    # println(Pauli_op)
                    if Pauli_op[1] == 0
                        # I
                        X_err[j] = 0
                        Z_err[j] = 0
                    elseif Pauli_op[1] == 1
                        # X
                        X_err[j] = 1
                        Z_err[j] = 0
                    elseif Pauli_op[1] == 2
                        # Y
                        X_err[j] = 1
                        Z_err[j] = 1
                    else
                        # Z
                        X_err[j] = 0
                        Z_err[j] = 1
                    end
                end

                # X syndrome
                LinearAlgebra.mul!(X_meas_syn, X_stabs, Z_err)
                @inbounds @simd for i in 1:nr_X
                    X_meas_syn[i] %= 2
                end

                # decode X
                X_flag, X_out, X_iter = CodingTheory._message_passing_layered(X_stabs, X_meas_syn,
                    X_chn_inits, CodingTheory._SP_check_node_message_box_plus, X_var_adj_list,
                    X_check_adj_list, max_iter, schedule, current_bits, totals,
                    X_syn, X_check_to_var_messages, X_var_to_check_messages, 0.0, X_layers)

                if !X_flag
                    # did not converge
                    skip = true
                else
                    # did converge
                    # shouldn't matter if I mod 2 this now or later
                    # a non-allocating version of Z_err += X_out
                    @inbounds @simd for i in 1:n
                        Z_err[i] = Z_err[i] + X_out[i]
                    end
                    X_local_iters[th][X_iter] += 1
                end

                # skip Z if X failed because the result will not end up in the codespace
                if !skip
                    # Z syndrome
                    LinearAlgebra.mul!(Z_meas_syn, Z_stabs, X_err)
                    @inbounds @simd for i in 1:nr_Z
                        Z_meas_syn[i] %= 2
                    end

                    # Bayes' Theorem: update priors on Z given the results from X
                    # always log(Pr(no error) / Pr(error)) so now
                    # log((Pr(I | I) + Pr(I | X)) / (Pr(Z | I) + Pr(Z | X))) etc
                    @inbounds @simd for i in 1:n
                        if X_out[i] == 0
                            Z_chn_inits[i] = without_X_err
                        else
                            Z_chn_inits[i] = with_X_err
                        end
                    end
                    
                    Z_flag, Z_out, Z_iter, = CodingTheory._message_passing_layered(Z_stabs, Z_meas_syn,
                        Z_chn_inits, CodingTheory._SP_check_node_message_box_plus, Z_var_adj_list,
                        Z_check_adj_list, max_iter, schedule, current_bits, totals,
                        Z_syn, Z_check_to_var_messages, Z_var_to_check_messages, 0.0, Z_layers)

                    if !Z_flag
                        # did not converge
                        skip = true
                    else
                        # did converge
                        # shouldn't matter if I mod 2 this now or later
                        # a non-allocating version of X_err += Z_out
                        @inbounds @simd for i in 1:n
                            X_err[i] = X_err[i] + Z_out[i]
                        end
                        Z_local_iters[th][Z_iter] += 1
                    end
                end
```

Finally, we determine if the decoder (error + correction) implemented a logical error and report the results. Here we return the frame (word) error rate and the convergence statistics of the two decoders. Different, more interesting, statistics may be collect, computed, or reported if desired.

```
                if !skip
                    # converged, check for logical errors
                    @inbounds for j in 1:nr_logs
                        iseven(dot(view(logs, j, 1:n), Z_err) - dot(view(logs, j, n +
                            1:2 * n), X_err)) || (local_log_failure_counts[th] += 1; break;)
                    end
                else
                    # did not converge
                    local_conv_failure_counts[th] += 1
                end

                # reset inputs for next run, but don't re-allocate new memory
                X_check_to_var_messages .= 0.0
                X_var_to_check_messages .= 0.0
                Z_check_to_var_messages .= 0.0
                Z_var_to_check_messages .= 0.0
                skip = false
            end
        end
        println("logical failures: $local_log_failure_counts")
        println("convergence failures: $local_conv_failure_counts")
        @inbounds FER[i] = (sum(local_conv_failure_counts) + sum(local_log_failure_counts)) / new_num_runs
        println("FER = $(FER[i])")
        println("Finished p = $p")

        local_conv_failure_counts[:] .= 0
        local_log_failure_counts[:] .= 0
    end

        # reduce iteration counts
    @inbounds for i in 2:num_threads
        @simd for j in 1:max_iter
            X_local_iters[1][j] += X_local_iters[i][j]
            Z_local_iters[1][j] += Z_local_iters[i][j]
        end
    end

    return FER, X_local_iters[1], Z_local_iters[1]
end
```


Let's test this with the following code [cite].

```
using LinearAlgebra, Oscar, Makie, CodingTheory
F = GF(2);
S, x = PolynomialRing(F, "x");
l = 127;
R = ResidueRing(S, x^l - 1);
a = 1 + x^15 + x^20 + x^28 + x^66;
b = 1 + x^58 + x^59 + x^100 + x^121;
aR = R(a);
bR = R(b);
A1 = GeneralizedBicycleCode(aR, bR);

julia> @time FER, X_iters, Z_iters = CSS_decoder_simple(A1);
Number of threads: 16, runs per thread: 938, new number of runs: 15008
Starting p = 0.010000000000000002
logical failures: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
convergence failures: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
FER = 0.0
Finished p = 0.010000000000000002
Starting p = 0.012589254117941675
logical failures: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
convergence failures: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
FER = 0.0
Finished p = 0.012589254117941675
Starting p = 0.015848931924611134
logical failures: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
convergence failures: [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
FER = 6.663113006396588e-5
Finished p = 0.015848931924611134
Starting p = 0.0199526231496888
logical failures: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
convergence failures: [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]
FER = 0.00013326226012793177
Finished p = 0.0199526231496888
Starting p = 0.025118864315095794
logical failures: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
convergence failures: [0, 0, 0, 1, 1, 0, 0, 2, 0, 1, 0, 0, 3, 0, 0, 0]
FER = 0.0005330490405117271
Finished p = 0.025118864315095794
Starting p = 0.03162277660168379
logical failures: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
convergence failures: [1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 0, 0, 2, 2, 2, 2]
FER = 0.0012659914712153518
Finished p = 0.03162277660168379
Starting p = 0.039810717055349734
logical failures: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
convergence failures: [5, 7, 2, 6, 6, 8, 7, 3, 6, 5, 3, 9, 4, 5, 8, 2]
FER = 0.0057302771855010665
Finished p = 0.039810717055349734
Starting p = 0.05011872336272722
logical failures: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
convergence failures: [22, 16, 27, 15, 20, 17, 27, 15, 33, 21, 16, 34, 21, 23, 20, 19]
FER = 0.023054371002132198
Finished p = 0.05011872336272722
Starting p = 0.06309573444801933
logical failures: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
convergence failures: [98, 95, 97, 74, 82, 95, 102, 98, 99, 78, 84, 108, 82, 103, 88, 98]
FER = 0.09868070362473348
Finished p = 0.06309573444801933
Starting p = 0.07943282347242814
logical failures: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
convergence failures: [303, 286, 287, 314, 288, 286, 295, 307, 299, 304, 284, 294, 298, 290, 308, 314]
FER = 0.3169642857142857
Finished p = 0.07943282347242814
Starting p = 0.1
logical failures: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
convergence failures: [633, 662, 644, 635, 601, 644, 652, 645, 627, 634, 645, 627, 642, 623, 642, 658]
FER = 0.6805703624733476
Finished p = 0.1
 105.068949 seconds (92.89 M allocations: 3.176 GiB, 0.40% gc time)
```


The first thing we immediately notice from the printout is that there were no logical errors. This is obviously not going to be true for other inputs such as the Steane code or surface codes; however, this is commonly claimed to be true throughout the quantum belief propagation literature. Second, the first two data points had no converge failures, implying that either they need more sample, or more likely, these error rates should be importance sampled. These will be automatically removed by the plotting software since we are using a log-log plot.

```
noise = 10 .^(-2:0.1:-1);
fig = Figure()
ax = Axis(fig[1, 1],
    limits = ((noise[1], noise[end]), (1e-6, 1)),
    xscale = log10,
    yscale = log10,
    yminorgridvisible = true,
    xminorgridvisible = true,
    title = "CSS Decoder - Simple",
    xlabel = L"Depolarizing Probability $p$",
    ylabel = "FER"
)
lines!(noise, FER, color = :red)
scatter!(noise, FER, color = :red)
current_figure()
```
![CSS_FER_test](./../assets/images/CSS_FER_test.png)


TODO REDO THIS FUNCTION

While the maximum number of iterations in this example was set to 20, the large majority of convergences happened quickly. To make the plots easier to read, let's filter out the first 15 iterations. 

```
X_iters_15 = filter(((k, v), ) -> k <= 15, X_iters)
Z_iters_15 = filter(((k, v), ) -> k <= 15, Z_iters)
X_values = collect(values(X_iters_15))
Z_values = collect(values(Z_iters_15))
m = max(maximum(X_values), maximum(Z_values))
fig2 = Figure()
ax1 = Axis(fig2[1, 1], xlabel = "Iteration Number", ylabel = "Number Of Times Converged",
        title = L"$X$ Convergence Distribution", limits = ((0, 16), (0, m + 500)))
ax2 = Axis(fig2[1, 2], xlabel = "Iteration Number", ylabel = "Number Of Times Converged",
    title = L"$Z$ Convergence Distribution", limits = ((0, 16), (0, m + 500)))
barplot!(ax1, collect(keys(X_iters_15)), X_values, bar_width = 1)
barplot!(ax2, collect(keys(Z_iters_15)), Z_values, bar_width = 1)
```
![CSS_test_iters](./../assets/images/CSS_test_iters.png)

The overwhelming majority of convergeneces, for both $X$ and $Z$, occurred within one or two iterations. This is plausible for a couple of reasons. First, note that all variable nodes have degree four but all check nodes have degree nine! This is large for an LDPC code. For low error rates when errors are sparsely distributed over 254 qubits, it may be common that a single check node does not connect to more than one inncorrect variable node and the degrees are high enough to immediately flip any bit.

```
julia> LDPCCode(X_stabilizers(A1))
[254, 141]_2 regular (5, 10)-LDPC code with density 0.03937007874015748.

Variable degree polynomial:
        x^4
Check degree polynomial:
        x^9


julia> LDPCCode(Z_stabilizers(A1))
[254, 141]_2 regular (5, 10)-LDPC code with density 0.03937007874015748.

Variable degree polynomial:
        x^4
Check degree polynomial:
        x^9
```

Second, note that the $Z$ stabilizers converged more often than the $X$ stabilizers. Unless one is familiar with the code family, it may be unclear that this code is not symmetric.

```
julia> CodingTheory._has_equivalent_row_spaces(X_stabilizers(A1), Z_stabilizers(A1))
false
```


Rerunning the simulation without using Bayes' Theorem returns almost symmetric $X$ and $Z$ iteration counts. For completeness, we include both runs on a single plot. The first blue point on the left is due to a single convergence error and the spike on the left further shows that direct sampling is either not appropriate for this error rate or it has not been sampled enough times for accuracy.

![CSS_FER_Bayes_test](./../assets/images/CSS_FER_Bayes_test.png)

## Example 2: Single-Shot Decoding With Metachecks
Next we're going to look at two single-shot decoding schemes. We will call the paper [quintavalle2021single](@cite) scheme one and [higgott2023improved](@cite) scheme two. We encourage the reader to check out both papers directly for details. Briefly, both schemes will consider data errors, as in the previous example, plus additional measurement errors (on the syndrome values). The code family we will look at has an extra matrix, $M$, with the property that $Ms = 0$ for any valid syndrome $s$ of the code. Then assuming that the measurement error $s_e$ didn't take us from a valid syndrome to another valid syndrome, $M(s + s_e) = Ms_e \neq 0$. Whether or not this happens depends on the properties of the classical code with $M$ as its parity-check matrix. To correct the syndrome, we decode using the Tanner graph based on $M$. Then we will use the corrected syndrome to decode the stabilizers.

The authors of scheme one noticed that sometimes correcting the syndrome produced an $s^\prime$ such that $Ms^\prime = 0$ but $s^\prime$ was not a valid syndrome. This occurs when $s^\prime = s + l$ where $s$ is a valid syndrome and $l$ is an element of the second homology. Using the fact that elements of homology and cohomology are paired [hatcher2005algebraic](@cite), we can generate a basis $L$ for this group and then use it to test if the syndrome contains an element of homology: $Ls^\prime \neq 0$. If this is true, then scheme one proceeds to redecode the measured syndrome with error using the Tanner graph corresponding to the matrix

$$\begin{pmatrix} M \\ L\end{pmatrix}.$$

The authors of scheme two noticed that the homology calculations of scheme one can be avoided by trying to solve for a correct and valid syndrome at the same time using the Tanner graph for

$$\begin{pmatrix} H & I \\ 0 & M\end{pmatrix}.$$

We will refer to this matrix as the single-shot matrix. Besides being simpler, one advatange of this scheme is that the matrix $L$ is often dense, making it a poor candidate for BP decoding. Along these lines, note that every stabilizer code has metachecks, but few are known which have sparse metachecks.

Here we compare scheme one and scheme two. Although it is a bit inconvenient, it is most accurate to compare on the same errors, including the same measurement errors. This means ininitalizing four sepearate BP decoders (scheme one: stabilizers. metachecks, $M + L$, scheme two: single-shot matrix) and using extra memory to not overwrite data in scheme one which needs to be used in scheme two. The stabilizer Tanner graph is initialized as before. The $M$ and $M + L$ Tanner graphs should use $\log((1 - p_\mathrm{meas}) / p_\mathrm{meas})$ to initialize the channel and the metacheck syndrome to decode. The single-shot matrix should use $\log((1 - p) / p)$ on the data qubits, $\log((1 - p_\mathrm{meas}) / p_\mathrm{meas})$ on the new columns corresponding to the metacheck, and the syndrome $(s \,\, Ms)^T$. Although CodingTheory makes it easy to compute homologies, we will eliminate a potential source of mismatch with the literature by using the matrices conveniently provided by the authors of scheme one [Vasmer_2020](@cite). A value of $p_\mathrm{meas} \in (10^{-3}, 10^{-2})$ is appropriate for many state-of-the-art quantum computing architectures. Note that both papers set $p_\mathrm{meas}$ equal to the depolarizing probability instead.

The failure condition of single-shot decoding is a bit different from the previous example. Due to measurement errors, there may be a small residue error left over after decoding, which will be handled in the next round. This proceeds for a fixed number of rounds, $L - 1$, and is followed by a round of perfect measurements which clears up any remaining residual errors. As a result, we don't care if the BP decoders converge until the last round, at which point we also check for logical errors. Note that this means that if $t$ samples are taken, $tL$ rounds are run, increasing both the time and memory requirements.

```

function _make_single_shot_Tanner_graph(H::T, M::T) where T <: CodingTheory.CTMatrixTypes
    # internal function, so skip checks on correctness
    F = base_ring(H)
    Fone = F(1)
    nr_H, nc_H = size(H)
    nc_M = ncols(M)
    single_shot_matrix = CodingTheory.direct_sum(H, M)
    row = 1
    @inbounds @simd for c in nc_H + 1:nc_H + nc_M
        single_shot_matrix[row, c] = Fone
        row += 1
    end
    # of form (H I; 0 M)
    return single_shot_matrix
end

function CSS_metachecks(X_stabs::T, X_logs::T, X_meta::T, X_meta_L::T) where T <: CodingTheory.CTMatrixTypes
    nr, n = size(X_stabs)
    n == ncols(X_logs) || throw(ArgumentError("Logs are the wrong size"))
    ncols(X_meta) == nr || throw(ArgumentError("Metacheck is the wrong size"))
    iszero(X_meta * X_stabs) || throw(ArgumentError("Either metachecks or stabilizers wrong"))
    X_ss = _make_single_shot_Tanner_graph(X_stabs, X_meta)

    # noise details
    p_meas = 1e-3
    chn_meas =  MPNoiseModel(:BSC, p_meas)
    chn_inits_meta_in = Float64[log((1 - p_meas) / p_meas) for _ in 1:nr]
    syn_err_weights = Weights(Float64[1 - p_meas, p_meas])
    noise = 10 .^(-3:.25:-1.5) ∪ 10 .^(-1.5:.1:-1)
    len_noise = length(noise)

    # BP details
    num_runs = 15000
    max_iter = 20
    # 8 w/ measurement error, 1 perfect
    single_shot_rounds = 9
    schedule = :layered
    X_layers_stabs = layered_schedule(X_stabs)
    X_layers_meta = layered_schedule(X_meta)
    X_layers_meta_L = layered_schedule(X_meta_L)
    X_layers_ss = layered_schedule(X_ss)
    erasures = Int[]
    nr_logs = size(X_logs, 1)
    Pauli_types = Int[0, 1]
    v = zeros(Int, n)
    v2 = zeros(Int, nr)
    v3 = zeros(Int, n + nr)
    X_logs_Int = CodingTheory._Flint_matrix_to_Julia_int_matrix(X_logs)

    # multi-threading details
    num_threads = Threads.nthreads()
    runs_per_thread = cld(num_runs, num_threads)
    new_num_runs = runs_per_thread * num_threads
    println("Number of threads: $num_threads, runs per thread: $runs_per_thread, new number of runs: $new_num_runs")

    # preallocate places to store results
    FER_scheme_1 = zeros(Float64, len_noise)
    FER_scheme_2 = zeros(Float64, len_noise)
    local_log_failure_counts_scheme_1 = [zeros(Int, nr_logs) for _ in 1:num_threads]
    local_log_failure_counts_scheme_2 = [zeros(Int, nr_logs) for _ in 1:num_threads]
    local_failure_counts_scheme_1 = zeros(Int, num_threads)
    local_failure_counts_scheme_2 = zeros(Int, num_threads)
    for i in 1:len_noise
        p = noise[i]
        println("Starting p = $p")

        # setup noise model
        chn = MPNoiseModel(:BSC, p)
        Pauli_weights = Weights(Float64[1 - p, p])
        chn_inits_stabs_in = Float64[log((1 - p) / p) for _ in 1:n]
        chn_inits_ss_in = Float64[chn_inits_stabs_in; chn_inits_meta_in]

        Threads.@threads for th in 1:num_threads
            # Tanner graph for the stabilizers
            X_stabs_Int, _, var_adj_list_stabs, check_adj_list_stabs, chn_inits_stabs,
                check_to_var_messages_stabs, var_to_check_messages_stabs, current_bits_stabs,
                totals_stabs, syn_stabs = CodingTheory._message_passing_init(X_stabs, v, chn, max_iter, :SP,
                chn_inits_stabs_in, schedule, erasures)
            # Tanner graph for the metachecks
            X_meta_Int, _, var_adj_list_meta, check_adj_list_meta, chn_inits_meta,
                check_to_var_messages_meta, var_to_check_messages_meta, current_bits_meta,
                totals_meta, syn_meta = CodingTheory._message_passing_init(X_meta, v2, chn_meas, max_iter, :SP,
                chn_inits_meta_in, schedule, erasures)
            # Tanner graph for the second metacheck matrix
            X_meta_L_Int, _, var_adj_list_meta_L, check_adj_list_meta_L, chn_inits_meta_L,
                check_to_var_messages_meta_L, var_to_check_messages_meta_L, current_bits_meta_L,
                totals_meta_L, syn_meta_L = CodingTheory._message_passing_init(X_meta_L, v2, chn_meas, max_iter, :SP,
                chn_inits_meta_in, schedule, erasures)
            # Tanner graph for scheme 2
            X_ss_Int, _, var_adj_list_ss, check_adj_list_ss, chn_inits_ss, check_to_var_messages_ss,
                var_to_check_messages_ss, current_bits_ss, totals_ss, syn_ss =
                CodingTheory._message_passing_init(X_ss, v3, chn, max_iter, :SP, chn_inits_ss_in, schedule,
                erasures)
            X_meas_syn_scheme_1 = zeros(Int, nr)
            X_meas_syn_scheme_2 = zeros(Int, nr)
            X_syn_err = zeros(Int, nr)
            X_syn_temp = zeros(Int, nr)
            Z_err = zeros(Int, n)
            state_scheme_1 = zeros(Int, n)
            state_scheme_2 = zeros(Int, n)
            m = zeros(Int, size(X_meta_Int, 1))
            m_2 = zeros(Int, size(X_meta_Int, 1))
            m2 = zeros(Int, size(X_meta_L_Int, 1))
            logs_check = zeros(Int, nr_logs, 1)

            for _ in 1:runs_per_thread
                X_flag_stabs = false
                X_flag_ss = false
                X_out_stabs = Int[0]
                X_out_ss = Int[0]
                for ss_iter in 1:single_shot_rounds
                    # sample
                    sample!(Pauli_types, Pauli_weights, Z_err)
                    @inbounds @simd for j in i:n
                        state_scheme_1[j] += Z_err[j]
                        state_scheme_2[j] += Z_err[j]
                    end

                    # syndrome
                    LinearAlgebra.mul!(X_meas_syn_scheme_1, X_stabs_Int, state_scheme_1)
                    LinearAlgebra.mul!(X_meas_syn_scheme_2, X_stabs_Int, state_scheme_2)
                    @inbounds @simd for j in 1:nr
                        X_meas_syn_scheme_1[j] %= 2
                        X_meas_syn_scheme_2[j] %= 2
                    end

                    # last measurement is considered perfect
                    if ss_iter != single_shot_rounds
                        # add measurement errors
                        sample!(Pauli_types, syn_err_weights, X_syn_err)
                        @inbounds @simd for j in 1:nr
                            X_meas_syn_scheme_1[j] = (X_meas_syn_scheme_1[j] + X_syn_err[j]) % 2
                            X_meas_syn_scheme_2[j] = (X_meas_syn_scheme_2[j] + X_syn_err[j]) % 2
                        end
                    end

                    ###########
                    # Scheme 1
                    ###########

                    LinearAlgebra.mul!(m, X_meta_Int, X_meas_syn_scheme_1)
                    LinearAlgebra.mul!(m_2, X_meta_Int, X_meas_syn_scheme_2)
                    @inbounds @simd for j in 1:length(m)
                        m[j] %= 2
                        m_2[j] %= 2
                    end

                    if !iszero(m)
                        # decode with metachecks
                        X_flag_meta, X_out_meta, _, = CodingTheory._message_passing_layered(X_meta_Int,
                            X_meas_syn_scheme_1, chn_inits_meta,
                            CodingTheory._SP_check_node_message_box_plus, var_adj_list_meta,
                            check_adj_list_meta, max_iter, schedule, current_bits_meta,
                            totals_meta, syn_meta, check_to_var_messages_meta,
                            var_to_check_messages_meta, 0.0, X_layers_meta)

                        # reset inputs for next run, but don't re-allocate new memory
                        check_to_var_messages_meta .= 0.0
                        var_to_check_messages_meta .= 0.0

                        if X_flag_meta
                            @inbounds @simd for j in 1:nr
                                X_syn_temp[j] = (X_meas_syn_scheme_1[j] + X_out_meta[j]) % 2
                            end

                            # check if has element of homology
                            LinearAlgebra.mul!(m2, X_meta_L_Int, X_syn_temp)
                            @inbounds @simd for j in 1:length(m2)
                                m2[j] %= 2
                            end

                            if !iszero(m2)
                                # println("decoding other metachecks")
                                X_flag_meta_L, X_out_meta_L, _, = 
                                CodingTheory._message_passing_layered(X_meta_L_Int, X_meas_syn_scheme_1,
                                    chn_inits_meta_L, CodingTheory._SP_check_node_message_box_plus,
                                    var_adj_list_meta_L, check_adj_list_meta_L, max_iter, schedule,
                                    current_bits_meta_L, totals_meta_L, syn_meta_L,
                                    check_to_var_messages_meta_L, var_to_check_messages_meta_L,
                                    0.0, X_layers_meta_L)
    
                                # reset inputs for next run, but don't re-allocate new memory
                                check_to_var_messages_meta_L .= 0.0
                                var_to_check_messages_meta_L .= 0.0
    
                                if X_flag_meta_L
                                    @inbounds @simd for j in i:nr
                                        X_meas_syn_scheme_1[j] += X_out_meta_L[j]
                                        X_meas_syn_scheme_1[j] %= 2
                                    end
                                else
                                    @inbounds @simd for j in i:nr
                                        X_meas_syn_scheme_1[j] = X_syn_temp[j]
                                    end
                                end    
                            else
                                @inbounds @simd for j in i:nr
                                    X_meas_syn_scheme_1[j] = X_syn_temp[j]
                                end
                            end
                        else
                            X_flag_meta_L, X_out_meta_L, _, = 
                            CodingTheory._message_passing_layered(X_meta_L_Int, X_meas_syn_scheme_1,
                                chn_inits_meta_L, CodingTheory._SP_check_node_message_box_plus,
                                var_adj_list_meta_L, check_adj_list_meta_L, max_iter, schedule,
                                current_bits_meta_L, totals_meta_L, syn_meta_L,
                                check_to_var_messages_meta_L, var_to_check_messages_meta_L,
                                0.0, X_layers_meta_L)

                            # reset inputs for next run, but don't re-allocate new memory
                            check_to_var_messages_meta_L .= 0.0
                            var_to_check_messages_meta_L .= 0.0

                            if X_flag_meta_L
                                @inbounds @simd for j in i:nr
                                    X_meas_syn_scheme_1[j] += X_out_meta_L[j]
                                    X_meas_syn_scheme_1[j] %= 2
                                end
                            end
                        end
                    end

                    # decode stabilizers
                    X_flag_stabs, X_out_stabs, _, = CodingTheory._message_passing_layered(X_stabs_Int, 
                        X_meas_syn_scheme_1, chn_inits_stabs, CodingTheory._SP_check_node_message_box_plus,
                        var_adj_list_stabs, check_adj_list_stabs, max_iter, schedule,
                        current_bits_stabs, totals_stabs, syn_stabs, check_to_var_messages_stabs,
                        var_to_check_messages_stabs, 0.0, X_layers_stabs)

                    # reset inputs for next run, but don't re-allocate new memory
                    check_to_var_messages_stabs .= 0.0
                    var_to_check_messages_stabs .= 0.0

                    if X_flag_stabs
                        @inbounds @simd for j in 1:n
                            state_scheme_1[j] += X_out_stabs[j]
                        end
                    end

                    ###########
                    # Scheme 2
                    ###########

                    # decode single-shot matrix
                    X_flag_ss, X_out_ss, _, = CodingTheory._message_passing_layered(X_ss_Int, Int[X_meas_syn_scheme_2; m_2],
                        chn_inits_ss, CodingTheory._SP_check_node_message_box_plus, var_adj_list_ss,
                        check_adj_list_ss, max_iter, schedule, current_bits_ss, totals_ss,
                        syn_ss, check_to_var_messages_ss, var_to_check_messages_ss, 0.0, X_layers_ss)

                    # reset inputs for next run, but don't re-allocate new memory
                    check_to_var_messages_ss .= 0.0
                    var_to_check_messages_ss .= 0.0

                    # only take the bits corresponding to the data error
                    @inbounds @simd for j in 1:n
                        state_scheme_2[j] += X_out_ss[j]
                    end
                end

                ###########
                # Scheme 1
                ###########

                # check if scheme 1 converged - returned to the code subspace
                if X_flag_stabs
                    # converged
                    # check for logical error
                    logs_flag = false
                    LinearAlgebra.mul!(logs_check, X_logs_Int, state_scheme_1)
                    @inbounds for j in 1:nr_logs
                        logs_check[j, 1] %= 2
                        if logs_check[j, 1] == 1
                            logs_flag = true
                            local_log_failure_counts_scheme_1[th][j] += 1
                        end
                    end
                    if logs_flag
                        # introduced a logical error
                        local_failure_counts_scheme_1[th] += 1
                    end
                else
                    # did not converge
                    local_failure_counts_scheme_1[th] += 1
                end

                ###########
                # Scheme 2
                ###########

                # check if scheme 2 converged - returned to the code subspace
                if X_flag_ss
                    # converged
                    # check for logical error
                    logs_flag = false
                    LinearAlgebra.mul!(logs_check, X_logs_Int, state_scheme_2)
                    @inbounds for j in 1:nr_logs
                        logs_check[j, 1] %= 2
                        if logs_check[j, 1] == 1
                            logs_flag = true
                            local_log_failure_counts_scheme_2[th][j] += 1
                        end
                    end
                    if logs_flag
                        # introduced a logical error
                        local_failure_counts_scheme_2[th] += 1
                    end
                else
                    # did not converge
                    local_failure_counts_scheme_2[th] += 1
                end

                # reset inputs for next run, but don't re-allocate new memory
                state_scheme_1 .= 0
                state_scheme_2 .= 0
            end
        end
        
        # reduce iteration counts
        @inbounds for i in 2:num_threads
            @simd for j in 1:nr_logs
                local_log_failure_counts_scheme_1[1][j] += local_log_failure_counts_scheme_1[i][j]
                local_log_failure_counts_scheme_2[1][j] += local_log_failure_counts_scheme_2[i][j]
            end
        end

        println("logical failures 1: $(local_log_failure_counts_scheme_1[1])")
        println("logical failures 2: $(local_log_failure_counts_scheme_2[1])")
        println("failures 1: $local_failure_counts_scheme_1")
        println("failures 2: $local_failure_counts_scheme_2")
        @inbounds FER_scheme_1[i] = sum(local_failure_counts_scheme_1) / new_num_runs
        @inbounds FER_scheme_2[i] = sum(local_failure_counts_scheme_2) / new_num_runs
        println("FER 1 = $(FER_scheme_1[i])")
        println("FER 2 = $(FER_scheme_2[i])")
        println("Finished p = $p")

        # reset inputs for next run, but don't re-allocate new memory
        @inbounds for i in 1:num_threads
            local_log_failure_counts_scheme_1[i] .= 0
            local_log_failure_counts_scheme_2[i] .= 0
        end
        local_failure_counts_scheme_1 .= 0
        local_failure_counts_scheme_2 .= 0
    end
    return FER_scheme_1, FER_scheme_2
end
```

Here we use eight rounds of noisy syndromes followed by a round of perfect measurement. As before, we will only do a rough sample of 15,000 points. The maximum number of iterations is set to 20, although the data suggests that raising this could improve convergence. Note that both references use BP+OSD, so our results are not going to be directly comparable. Also, the code family in [Vasmer_2020](@cite) is the 3D toric code. This is known for performing extremely poorly under BP decoding (due to trapping sets), and this is exactly what we expect to see that here. Later, we will apply other decoding to the same problem to show our progress.

The first thing we notice is that scheme one takes an incredibly long time due to the previously mentioned check-node degrees
```
function check_weights(H::CodingTheory.CTMatrixTypes)
    w = maximum(count(!iszero, H[i, :]) for i in 1:nrows(H))
    q = maximum(count(!iszero, H[:, i]) for i in 1:ncols(H))
    @show (w, q)
    return nothing
end

F = GF(2);
X_meta_L = matrix(F, load_alist("data/toric3D/toric3D_2_mxlmx.alist"));
julia> check_weights(X_meta_L)
(w, q) = (6, 3)

X_meta_L = matrix(F, load_alist("data/toric3D/toric3D_3_mxlmx.alist"));
julia> check_weights(X_meta_L)
(w, q) = (9, 3)

X_meta_L = matrix(F, load_alist("data/toric3D/toric3D_4_mxlmx.alist"));
julia> check_weights(X_meta_L)
(w, q) = (16, 3)

X_meta_L = matrix(F, load_alist("data/toric3D/toric3D_5_mxlmx.alist"));
julia> check_weights(X_meta_L)
(w, q) = (25, 3)

X_meta_L = matrix(F, load_alist("data/toric3D/toric3D_6_mxlmx.alist"));
julia> check_weights(X_meta_L)
(w, q) = (36, 3)

X_meta_L = matrix(F, load_alist("data/toric3D/toric3D_7_mxlmx.alist"));
julia> check_weights(X_meta_L)
(w, q) = (49, 3)

X_meta_L = matrix(F, load_alist("data/toric3D/toric3D_8_mxlmx.alist"));
julia> check_weights(X_meta_L)
(w, q) = (64, 3)

X_meta_L = matrix(F, load_alist("data/toric3D/toric3D_9_mxlmx.alist"));
julia> check_weights(X_meta_L)
(w, q) = (81, 3)
```

By distance eight, scheme one was difficult to run without a cluster. We did not attempt distance nine with this scheme. This is problematic for many reasons, the most important of which is that many code families do not "settle in" to their asymptotic behaviors until distances much higher than this (although the exact distance depends on the decoder being used). For example, for the surface codes under minimum-weight perfect-matching (MWPM), anything below distance 20 is considered the small-code regime (compare this to distance seven for the same code family using trellis decoding). 

![CSS_Single_Shot_test](./../assets/images/CSS_Single_Shot_test.png)

Logical errors were equally distributed among the cosets and only occurred on odd-distance codes.