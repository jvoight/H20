# Checks reconstructed from the manuscript

These scripts were reconstructed on 15 November 2026 from the statements
and numerical tables in `rewrite-v11.tex` and the preserved surface models.
They replace the missing scripts cited there; they are not recovered copies
of the original implementations. Both pass with Magma 2.30-1.

## Run

From this directory, with Magma on your path:

```sh
magma -b herfurtner_checks.m
magma -b modular_dimensions.m
```

On the local Apple Silicon installation, use `MAGMA_ARCH=arm64 magma`.
The scripts need no network access, external packages, or downloaded
eigenvalues. The Herfurtner check loads `herfurtner_models.m` and the shared
utilities at `../202-candidates/scripts/euler_tools.m`.

Saved output is in `results/`. Require the final `PASS` line **and no errors
anywhere in the log**: Magma may exit with status zero after an error. Each
main check runs in a procedure, so a failed assertion aborts before its
`PASS` message. Neither script overwrites saved results automatically.

## Herfurtner checks

`herfurtner_models.m` transcribes the five relevant pairs of polynomials
from [`../legacy-surfaces/maximal.m`](../legacy-surfaces/maximal.m).
For each pair `(g2,g3)`, the untwisted equation is
`y^2 = x^3 - g2*x + 2*g3`. The supplied quadratic twist is the product of
the monic finite bad-place polynomials, exactly as in that source. The
verification recomputes this product and checks the explicit twist; its
constant squareclass is not adjusted to obtain a match.

For H7, H13, H15, H25, and H28, the script checks the characteristic-zero
geometric fiber types, `b1=0`, `b2=34`, and trivial-lattice rank 30. Thus
`pg=2`, `h11=30`, and the geometric Mordell–Weil rank is zero. The root
lattice discriminants give orientation squareclass 3 for H7 and 1 for the
other four models. The possible torsion index changes a discriminant by
a square and does not affect this check.

The finite-field computations use Magma's
[`LFunction`](https://magma.maths.usyd.edu.au/magma/handbook/text/1581).
Because the characteristic-zero Néron–Severi space is the trivial lattice
here, the resulting degree-four polynomial retains the required
transcendental contribution, including any extra algebraic classes arising
after reduction. The script also checks that the geometric fiber types and
second Betti number agree at each tested prime.

- **H13, H15, H25, H28:** all seven primes `5,7,13,17,19,31,37` in v11's
  “Full-Euler diagnostics” table. All five polynomial coefficients and the
  complete sets of compatible QM component orders are checked: 28 entries.
- **H7:** the two displayed polynomials at 13 and 37. Both are recomputed,
  checked coefficient by coefficient, and shown irreducible with Galois
  group of order 8, hence `D4`. These are the two primes used for H7;
  in particular, 5 is bad for the displayed model, where the finite bad
  places `t=0` and `t=-10` collide.

The saved log records all 30 full polynomials, with coefficients in ascending
degree. This reproduces the manuscript's diagnostics; it does not identify
an Artin representation, prove a noncyclic QM component group, or exclude
every CM possibility for H7.

## Modular dimensions

`modular_dimensions.m` computes `dim S_3(Gamma_1(N))` for every level from
1 through 40 in two ways: Magma's dimension formula and the dimension of
full cuspidal modular symbols divided by two. It also compares the index
and regular-cusp formula for levels at least 5. The use of sign zero and the
factor of two follow the
[modular-symbols conventions](https://magma.maths.usyd.edu.au/magma/handbook/text/1721).

For levels 9, 10, and 12, it independently enumerates projective cosets as
primitive bottom rows modulo the level, up to sign, and computes cusp widths
as translation orbits. In particular, Gamma_1(12) has SL2 index 96, ten
regular cusps, and cusp-form dimension 3; levels 9 and 10 have dimension 2.

Finally it checks all six odd quadratic refinements in the manuscript:
`(12,-4), (14,-7), (15,-3), (15,-15), (16,-8), (18,-3)`.
Each has cusp-form dimension 2. The level-14 and level-16 spaces have new
dimension zero; the other four have new dimension 2.

The final call `VerifyModularDimensions(40)` sets the finite level bound.
The run checks the stated small-level dimensions; it does not turn a
finite enumeration into a proof of the assertion for all larger levels.
