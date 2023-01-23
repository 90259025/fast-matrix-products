Given an unknown polynomial $p \in \mathbb{Z}[x]$ of degree $d$ and the values of $p(0), p(1), ..., p(d)$, we find $p(a), p(a + 1), ..., p(a + d)$. There are two cases:

Case 1. $a \neq d + 1$. Define

$$
    \delta(i, d) = \prod_{j = 0, j \neq i}^d (i - j)
$$

and

$$
    \Delta(a, i, d) = \prod_{j = 0}^d (a + i - j).
$$

Lagrange Interpolation states that

$$
    P(a + k) = \sum_{i = 0}^d \tilde{P_i} \prod_{j = 0, j \neq i}^d (a + k - j),
$$

where $\tilde{P_i} = \frac{P(i)}{\delta(i, d)}$.

Define

$$
    Q_k = \sum_{i = 0}^d \tilde{P_i} \frac{1}{a + k - i}
$$

and

$$
    \tilde{P}(x) = \sum_{i = 0}^d \tilde{P_i} x^i
$$

and

$$
    S(x) = \sum_{i = 0}^d \frac{1}{a + i - d} x^i.
$$

Then for all $k$ between $0$ and $d$, $Q_k$ is the coefficient of $x^{k + d}$ in $\tilde{P}S$, and $P(a + k) = Q_k \prod_{j = 0}^d (a + k - j)$.

So the algorithm first computes all $d + 1$ different $\delta(i, d)$, then all the different $\tilde{P_i}$. Then it multiplies $\tilde{P}$ and $S$ and extracts the $d + 1$ coefficients that correspond to the values of $P(a + k)$.

Case 2. $a == d+1$.