# Message Passing

The following is intended to be a review and not an introductory lesson. Prior knowledge about the basic idea about the topic is assumed.

## Classical LDPC Codes
The most common type of decoder for LDPC codes are called *message passing* decoders. These work directly on the Tanner graph of the code. Messages are passed from the variable nodes to the check nodes constituting a single round. This iterates until some stopping criteria is met, usually either a maximum number of iterations or some level of “certainty” is reached in the current answer. The most common and general message passing algorithm is a form of belief propagation, which is used throughout many fields. Lesser known, but also used, algorithms include Gallager A and Gallager B, which will be covered in the next section.

The messages passed from variable nodes to check nodes of the Tanner graph in belief propagation are probabilities, log-likelihoods, or otherwise some metric of belief that the bit has a certain value given the observed value. The messages passed from check to variable are beliefs that the variable node has a certain value given all the messages passed to the check node from all of the *other* variable nodes. The key idea behind this is that only *extrinsic* information is considered. In other words, the information a variable node holds about itself is never considered in the calculation at the check node; the check node makes an estimate of what it thinks the variable node should be based on all the other variable nodes connected to the check node. The only way this works is if the information about the given variable node in question does not propagate through the graph and make its way to another variable node via a second check node which then gets passed to the original check node making an estimate of the original variable node. When this happens, the probability of bit $i$ being correct is now conditioned on other bits *and* itself, at which point the mathematical equations are no longer valid. Hence, belief propagation is only valid if the paths the messages travel in the Tanner graph unrolled in time form a tree. Once the path ceases to be a tree, i.e., a cycle was encountered, belief propagation is no longer valid. Thus, the success of the decoder is correlated to the length of the shortest cycle in the Tanner graph, which is called the *girth* of the graph. If $g$ is the girth, the first $g/2$ rounds of message passing are extrinsic (independent). In practice, one finds that message passing often works well in spite of cycles because the non-independent information is either too diluted at first to make a difference or it simply does not contradict the already formed belief. The latter of course depends on whether or not the variable node in question corresponds to an error.

One of the main ingredients in message passing is likelihoods. The difference between probabilities and likelihoods is deep and runs at the heart of interpretations of probability theory: frequentist, Bayesian, likelihoodist. for a good discussion. The key behind likelihoods is the *law of likelihood* which says that (the evidence) $E$ favors (hypothesis) $H_1$ over $H_2$ if and only if $L = \mathrm{Pr}(E \mid H_1) / \mathrm{Pr}(E \mid H_2) > 1$, which $L$ measuring the degree of favoring. According to the first above link, it has been suggested treating $L = 8$ as the threshold for declaring “fairly strong evidence” in favoring one hypothesis over another and $L = 32$ for “strong evidence”. This interpretation may be useful in assigning cutoffs in iterative decoding schemes, although the second and third links above shows we must apply this idea with caution without Bayesian priors. Belief propagation, however, is distinctively stemming from a Bayesian philosophy.

In Bayesian theory, if one learns (the evidence) $E$ with certainty and nothing else, then one should replace one’s prior degree of belief $\mathrm{Pr}(H)$ of any hypothesis with $\mathrm{Pr}(H \mid E)$. Applying Bayes’s,

$$
\frac{\mathrm{Pr}(H_1 \mid E)}{\mathrm{Pr}(H_2 \mid E)} = \frac{\mathrm{Pr}(H_1)}{\mathrm{Pr}(H_2)}\frac{\mathrm{Pr}(E \mid H_1)}{\mathrm{Pr}(E \mid H_2)},
$$

which says that the *a posterior* ratio for a pair of hypotheses following this belief update is their prior ratio times their likelihood ratio. This update holds for probability densities as well.

Before returning to message-passing, it is necessary to comment on exactly what problem it is trying to solve. We will utilize common notation that capital letters denote random variables, lowercase letter represent specific instances of them, and will simplify expressions such that $\mathrm{Pr}_X(X = x)$ to $\mathrm{Pr}_X(x)$. Let $\mathcal{C}$ be a binary linear code. Let $X \in \mathcal{C}$ be a random variable chosen by the transmitter with probability $\mathrm{Pr}_X(x)$ and $Y$ be a random variable over the output space. The error channel is described by its *transition probability* $\mathrm{Pr}_{Y | X}(y \mid x)$. If we decode $y$ to $\hat{x}(y) \in \mathcal{C}$, then the probability of error is $1 - \mathrm{Pr}_{X \mid Y}(\hat{x}(y) \mid y)$. Thus, to minimize the probability of error, we should choose $\hat{x}(y)$ to maximize this probability. This gives the *maximum a posteriori* (MAP) decoding rule

$$
\begin{align*}
\hat{x}^\text{MAP}(y) &= \argmax_{x \in \mathcal{C}} \mathrm{Pr}_{X \mid Y}(x \mid y)\\
&= \argmax_{x \in \mathcal{C}} \mathrm{Pr}_{Y \mid X}(y \mid x) \frac{\mathrm{Pr}_X(x)}{\mathrm{Pr}_Y(y)}\\
&= \argmax_{x \in \mathcal{C}} \mathrm{Pr}_{Y \mid X}(y \mid x) \mathrm{Pr}_X(x).
\end{align*}
$$

If all codewords are equally likely, then we get the *maximum likelihood* (ML) decoding rule

$$
\hat{x}^\text{MAP}(y) = \argmax_{x \in \mathcal{C}} \mathrm{Pr}_{Y \mid X}(y \mid x) = \hat{x}^\text{ML}(y),
$$

where the new terminology follows from the discussion above. The term $\mathrm{Pr}_{X \mid Y}(x \mid y)$ is the *a posteriori probability* (APP), so sometimes MAP decoders are referred to as APP decoders.

The problem with the MAP problem above is that it is global, considering all bits simultaneously. This is hard, so instead message passing chooses to solve the easier local problem of *bit-wise* MAP decoding. Let $\mathcal{C}$ be a binary linear code with random variable $X \in \mathcal{C}$. The $i$th element, $v_i$, is an element of $\{0, 1\}$, however, it is more convenient to work with the input $X_i \in \{\pm 1\}$ such that $X_i = (-1)^{v_i}$. We assume the channel is memoryless (i.i.d.) such that $\mathrm{Pr}_{Y \mid X}(y \mid x) = \prod_i \mathrm{Pr}_{Y_i \mid X_i}(y_i \mid x_i)$. Finally, let $\sim x_i$ denote every bit expect for the $i$th bit and $\mathbb{1}_{x \in \mathcal{C}}$ denote an indicator function on the code. The bit-wise MAP decoding rule is

$$
\begin{align*}
\hat{x}^\text{MAP}_i(y) &= \argmax_{x_i \in \{\pm 1\}} \mathrm{Pr}_{X_i \mid Y}(x_i \mid y)\\
&= \argmax_{x_i \in \{\pm 1\}} \sum_{\sim x_i} \mathrm{Pr}_{X \mid Y}(x \mid y)\\
&= \argmax_{x_i \in \{\pm 1\}} \sum_{\sim x_i}\mathrm{Pr}_{Y \mid X}(y \mid x) \mathrm{Pr}_X(x)\\
&= \argmax_{x_i \in \{\pm 1\}} \sum_{\sim x_i} \left( \prod_j \mathrm{Pr}_{Y_j \mid X_j}(y_j \mid x_j)\right) \mathbb{1}_{x \in \mathcal{C}}.
\end{align*}
$$

The last term here is a *marginal probability* and is able to be computed efficiently using factor graphs via the belief propagation algorithm. Fortunately, the Tanner graph happens to be a factor graph and so BP can be run on it directly. It is known that the bit-wise MAP problem is suboptimal, sometimes dramatically, compared to the full MAP problem.

!!! note "Example"
💡 For the Hamming code above,
    $\mathbb{1}_{x \in \mathcal{C}} = \mathbb{1}_{x_1 + x_2 + x_4 + x_5 = 0} \, \mathbb{1}_{x_1 + x_3 + x_4 + x_6 = 0} \, \mathbb{1}_{x_2 + x_3 + x_4 + x_7 = 0}$.

We will need the following lemma to continue.

!!! note "Lemma"
    Consider a vector of $n$ independent binary random variables $a = [a_0, \dots, a_{n - 1}]$ such that $\mathrm{Pr}(a_i = 0) = p^{(i)}_0$ and $\mathrm{Pr}(a_i = 1) = p^{(i)}_1$. Then the probability that $a$ contains an even number of 1’s is $\displaystyle \frac{1}{2} + \frac{1}{2} \prod_{i = 0}^{n - 1} \left( 1 - 2 p^{(i)}_1\right)$ and the probability that $a$ contains an odd number of 1’s is $\displaystyle \frac{1}{2} - \frac{1}{2} \prod_{i = 0}^{n - 1} \left( 1 - 2 p^{(i)}_1\right)$.

The initial information in the Tanner graph at round 0 is at the variable nodes (the value of the received bit). We begin by computing the likelihood ratio of each value. To increase numerical stability (and speed), we shift to log-likelihood ratios $\ell = \ln L(x \mid y)$. This changes the belief update formula above from multiplication to addition. The check nodes now have log-likelihood ratios from all of its connected variable nodes. It now computes an estimate that the bit attached to the $j$th edge is correct using the information from all the edges but that one, i.e., it uses only extrinsic information. Each check node can be seen as a length $w$ random bit vector $a$, where $w$ is the weight of the parity check. In order for the check to be satisfied, this vector must have an even number of 1’s. Without loss of generality, consider the $0$th edge. The probability that this is a 0 while satisfying the check is given by the probability that the remaining edges have an even number of 1’s,

$$
\mathrm{Pr}(a_0 = 0 \mid y) = \frac{1}{2} + \frac{1}{2}\prod_{i = 1}^{w - 1} (1 - 2 \mathrm{Pr}(a_i = 1 \mid y_i)).
$$

Rearranging gives

$$
2 \mathrm{Pr}(a_0 = 0 \mid y) - 1 = \prod_{i = 1}^{w - 1}(1 -2 \mathrm{Pr}(a_i = 1 \mid y_i)) = \prod_{i = 1}^{w - 1}(2 \mathrm{Pr}(a_i = 0 \mid y_i) - 1).
$$

A little algebra and the fact that the sum of all probabilities is 1 shows that $\mathrm{Pr}(a_i = 0 \mid y_i) = L(x_i \mid y_i) / (1 + L(x_i \mid y_i))$. Therefore, $2 \mathrm{Pr}(a_0 = 0 \mid y_i) - 1 = (L(x_i \mid y_i) - 1)/(L(x_i \mid y_i) + 1)$. Recalling that $\tanh x = \sinh x / \cosh x = (e^x - e^{-x})/(e^x + e^{-x}) = (e^{2x} - 1)/(e^{2x} + 1)$, we have that

$$
2 \mathrm{Pr}(a_0 = 0 \mid y_i) - 1 = \tanh(\ell_i / 2),
$$

where $\ell_i = \ln L(x_i \mid y_i)$. Therefore,

$$
2 \mathrm{Pr}(a_0 = 0 \mid y) - 1 = \prod_{i = 1}^{w - 1} \tanh(\ell_i / 2)
$$

Repeating for $a_0 = 1$, we get

$$
L(a_0 \mid y) = \frac{\mathrm{Pr}(a_0 = 0 \mid y)}{\mathrm{Pr}(a_0 = 1 \mid y)} = \frac{1 + \prod_{i = 1}^{w - 1} \tanh(\ell_i / 2)}{1 - \prod_{i = 1}^{w - 1} \tanh(\ell_i / 2)}.
$$

Sticking with the log domain, we pass the natural log of this value back to the variable node on edge 0, which then updates its belief in its own value. A check node repeats this for every edge, computing the likelihoods on the $i$th edge by using all the other edges. The message nodes update their beliefs using the messages passed from every edge (check node).

The more common way to see this equation written is by applying the same hyperbolic tangent trick to the full probability to find $\tanh(\ell / 2) = \prod^{w - 1}_{i = 1} \tanh(\ell_i / 2)$, where $\ell = \ln L(a_0 \mid y)$, in which case we find the update step

$$
\ell = 2 \tanh^{-1} \left( \prod_{i = 1}^{w - 1} \tanh(\ell_i / 2)\right).
$$

To summarize, let $m^{(i)}_{vc}$ denote the messages passed from the variable nodes to the check nodes at round $i$ and $m^{(i)}_{cv}$ denote similarly denote the messaged passed from check to variable. Then

$$
m^{(i)}_{vc} = \begin{cases} \ln L(x_v \mid y) & i = 0,\\
\ln L(x_v \mid y) + \sum_{c^\prime \neq c} m^{(i - 1)}_{c^\prime v} & i \geq 1,
\end{cases}
$$

and

$$
m^{(i)}_{cv} = 2 \tanh^{-1} \left( \prod_{v^\prime \neq v} \tanh \left(m^{(i)}_{v^\prime c} / 2 \right)\right).
$$

Expressed in this form, message passing is called the *sum-product* *algorithm*.

The standard stopping condition for this algorithm is obtaining a zero syndrome. This may never occur, so one can also set a max number of of iterations. Note that these criteria ignore the actual value of the likelihoods.

The check-node message as written above is numerically slow and unstable. Numerous papers have been written improving this in various directions. Chapter 5 of [ryan2009channel](@cite) contains a good summary of the following modifications:

- Box-Plus Decoder
- Min-Sum Decoder
- Attenuated & Offset Min-Sum Decoder
- Min-Sum With Correction Decoder
- Approximate $\min^*$ (a-$\min^*$) Decoder
- Richardson/Novichkov Decoder
- Reduced-Complexity Box-Plus Decoder

As described, all variable nodes compute their messages in parallel and pass them in parallel. Likewise, all check nodes compute their messages in parallel and pass them in parallel. This is called a *parallel (flooding) update schedule*. Sometimes other schedules can help break the decoder out of local minima it is stuck in (trapping sets). In a *serial schedule*, all of the check nodes connected to the first variable node sends their messages to the variable node, which then sends its updated value back to all check nodes. The same procedure is repeated for all variable nodes. A *layered schedule* where some set of variable nodes is updated in parallel followed by the next set of variable nodes is also popular.

### Gallager A & B And Bit-Flipping Decoding Algorithms

Gallager’s original work contained three message-passing algorithms. The messages are not beliefs, so they are distinct from the above work. In Gallager A, the check nodes perform a simple mod-2 sum of all the variable bits on all but a given edge and then pass back on that edge what that bit must be to satisfy the parity check. If *any* check node passes back the value received on the channel, it uses the channel value. If *all* check nodes disagree with the channel, then it uses the value that the check nodes agree on. Gallager B is almost identical but instead flips the channel bit if at least some fixed threshold $t$ of check nodes disagree with the channel bit. Gallager provided a formula for the optimal value of $t$ for regular LDPC codes. These algorithms actually perform quite well on some channels.

The bit-flipping algorithm first evaluates all parity-check equations (rows of $H$) and then flips any bits in the received word that are involved in more than some fixed number of failed parity checks. This is repeated until all the parity checks are satisfied or a fixed number of iterations has been reached.

### Binary Symmetric Channel (BSC)

The probability of error, $p$, is called the *crossover probability.* Let $X_i$ be a random variable over the input $\{0, 1\}$ and $Y_i$ be a random variable over the output. The *channel transition probabilities* are

$$
\begin{align*}
\mathrm{Pr}_{Y \mid X}(1 \mid 0) &= \mathrm{Pr}_{X \mid Y}(0 \mid 1) = p,\\
\mathrm{Pr}_{X \mid Y}(1 \mid 1) &= \mathrm{Pr}_{X \mid Y}(0 \mid 0) = 1 - p
\end{align*},
$$

from which we find

$$
L(x_i, y_i) = (-1)^{y_i} \log\left(\frac{1 - p}{p}\right).
$$

This channel is memoryless, i.e., $\mathrm{Pr}_{X \mid Y}(Y | X) = \Pi_i \mathrm{Pr}_{X \mid Y}(y_i \mid x_i)$.

Consider decoding over the BSC using a $[n, k, d]$ linear code. The criterion for decoding the BSC is to minimize the probably of decoding error, which is equivalent to maximizing the *a posteriori* probability $\mathrm{Pr}_{X \mid Y}(x \mid y)$. Thus, the optimal solution is given by

$$
\hat{v} = \argmax_v \mathrm{Pr}_{X \mid Y}(x \mid y) = \argmax_v \frac{\mathrm{Pr}_{Y \mid X}(y \mid x) \mathrm{Pr}_X(x)}{\mathrm{Pr}_Y(y)}.
$$

The prior is uniform over the codewords, so $\mathrm{Pr}_X(x)$ is independent of $x$ and hence $v$, and $\mathrm{Pr}_Y(y)$ is also independent of $v$. Hence, the *maximum a posteriori* (MAP) criteria can be replaced by the *maximum-likelihood* (ML) criteria

$$
\hat{v} = \argmax_v \mathrm{Pr}_{Y \mid X}(y \mid x).
$$

Using the monotonicity of the $\log$ and letting $d_H(\cdot, \cdot)$ denote Hamming distance,

$$
\begin{align*}
\hat{v} &= \argmax_v \log \prod_i \mathrm{Pr}_{Y \mid X}(y_i \mid x_i)\\
&= \argmax_v \sum_i \log \mathrm{Pr}_{Y \mid X}(y_i \mid x_i)\\
&= \argmax \left[d_H(y, x) \log(p) + (n - d_H(y, x))\log(1 - p)\right]\\
&= \argmax_v \left[d_H(y, x) \log \left( \frac{p}{1 - p}\right) + n \log(1 - p)\right]\\
&= \argmin_v d_H(y, x).
\end{align*}
$$

In addition to belief propagation, Gallager’s A, B, and bit-flipping algorithms can be used for the BSC. The capacity of the BSC is $C_\mathrm{Sha} = 1 - H(p)$, where $H(x) = - x \log_2 x - (1 - x) \log_2(1 - x)$ is the *binary entropy function*.

### Binary Additive White Gaussian Noise Channel (BAWGNC)

The binary code bits $v_i \in \{0, 1\}$ are mapped to channel inputs by $x_i = (-1)^{v_i}$. The channel transition *probability density function* (pdf) is

$$
p(y_i \mid x_i) = \frac{1}{\sqrt{2 \pi} \sigma} \exp \left[-(y_i - x_i)^2/(2\sigma^2)\right],
$$

where $\sigma^2$ is the variance of the zero-mean Gaussian noise sample $n_i$ such that $y_i = x_i + n_i$. This channel is memoryless. One can show that

$$
\mathrm{Pr}_{X \mid Y}(x_i \mid y_i) = \left[ 1 + \exp\left(-2 y_i x_i / \sigma^2\right)\right]^{-1},
$$

from which we find $L(v_i \mid y_i) = 2 y_i / \sigma^2$. In practice, one must first estimate the variance.

Similar to above, the channel decoding problem for the BAWGNC comes out to be

$$
\begin{align*}
\hat{v} &= \argmax_v \mathrm{Pr}_{Y \mid X}(y \mid x)\\
&= \argmax_v \sum_i \log\left(\frac{1}{\sqrt{2\pi}\sigma} \exp \left[-(y_i - x_i)^2 / (2 \sigma^2)\right]\right)\\
&= \argmin_v \sum_i (y_i - x_i)^2\\
&= \argmin_v d_E(y, x),
\end{align*}
$$

where $d_E(\cdot, \cdot)$ is the standard Euclidean distance.

The capacity of the BAWGNC is

$$
C_\mathrm{Sha} = \frac{1}{2} \sum_{x_i = \pm 1} \int^\infty_{-\infty} p(y \mid x) \log_2 \left(\frac{p(y \mid x)}{p(y)}\right) dy,
$$

where $p(y) = \frac{1}{2} \left[p(y \mid 1) + p(y \mid -1)\right]$.