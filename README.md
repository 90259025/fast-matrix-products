# fast-matrix-products

This is a library to compute matrix products of the form

$$
  \prod_{i=0}^{k-1} M(i)
$$

where $M$ is an $n\times n$ matrix over $\mathbb Z[x]$ (or $\mathbb Z[x]/n\mathbb Z$) of bounded degree. 

The algorithm is inspired by [BGS07], which covers the case when $M$ has degree at most 1. However, much of the heavy lifting is actually done with Theorem 3.1 of [Sh91].

The case of interest is in computing zeta functions of nondegenerate hypersurfaces in toric varieties via p-adic cohomology. The linear case is used for Corollary 3.2 of [CHK18] and the higher degree case for Proposition 1.18 of [Cos15].


- [BGS07] https://specfun.inria.fr/bostan/publications/BoGaSc07.pdf
- [CHK18] https://arxiv.org/pdf/1806.00368.pdf
- [Cos15] https://math.mit.edu/~edgarc/files/EdgarCosta-PhDthesis.pdf
- [Sh91] https://dl.acm.org/doi/pdf/10.1145/120694.120697
