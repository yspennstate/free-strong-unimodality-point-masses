# Only point masses are freely strongly unimodal

Yitzchak Shmalo

Paper: [paper/main.pdf](paper/main.pdf) (source [paper/main.tex](paper/main.tex)).

A probability measure $\mu$ on $\mathbb R$ is freely strongly unimodal if $\mu\boxplus\nu$ is unimodal for every
unimodal $\nu$. Hasebe and Ueda (ALEA 15 (2018) 353-374, arXiv:1710.08240) showed that a measure with finite
variance has this property only if it is a point mass, and asked in their Problem 3.11 whether any other measure
has it. The paper shows that none does: for every $\mu$ that is not a point mass,
$\mu\boxplus\bigl(\tfrac1{100}\delta_0+\tfrac{99}{100}U[0,L]\bigr)$ is not unimodal for all large $L$. No moment
condition is needed.

## Code

Python 3.14 with NumPy 2.4 and mpmath 1.3. Each script writes its output next to itself.

- `code/check_semicircle.py`: the density of $\sigma\boxplus\kappa$, $\kappa=\tfrac1{100}\delta_0+\tfrac{99}{100}U[0,1]$,
  for $\sigma$ a semicircle law of variance $10^{-6}$ or $10^{-4}$ and for the symmetric Bernoulli law on
  $\{\pm10^{-3}\}$, by fixed-point iteration of the subordination equations; prints $f'(0.01)$, $f'(0.25)$ and the
  interior local minima (Section 4 of the paper, `code/check_semicircle.out`).
- `code/constants.py`: exact rational and 50-digit checks of the constants in Section 3 ($g'$, $g''$ at the two
  points, the bounds on $\Delta$, $\rho$ and $2r/\Delta$, and the two slopes $-200000/99$ and $704/375$)
  (`code/constants.out`).

## Funding

The research presented in this paper was supported by the European Research Council (ERC) under the European
Union's Horizon Europe research and innovation programme (grant agreement No. 101041711), by the Simons Foundation
as part of the Collaboration on the Mathematical and Scientific Foundations of Deep Learning, by Heights Labs, by
Convex Nexus Capital, by the Israel Science Foundation (grant number 2258/19), and by the Israel Science Foundation
(ISF Grant 4101/25). Computational resources were provided in part by Google Cloud through its education credits
program (credit 530550103).
