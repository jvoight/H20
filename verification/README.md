# Verification checks

The Herfurtner and modular-dimension scripts were reconstructed on 15 November 2026 from the statements
and numerical tables in `rewrite-v11.tex` and the preserved surface models.
They replace the missing scripts cited there; they are not recovered copies
of the original implementations. Both pass with Magma 2.30-1.

The independent verification of the rational rank-four QM counterexample
289.3.e.o was added on 20 November 2026; its reproduction instructions,
mathematical scope, and references appear below.

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
`PASS` message. No script overwrites saved results automatically.

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


## Rational rank-four QM counterexample: 289.3.e.o

Added on 20 November 2026. Run from this directory:

```sh
MAGMA_ARCH=arm64 magma -b qm_289_counterexample.m
```

Saved output: [`results/qm_289_counterexample.log`](results/qm_289_counterexample.log),
produced with Magma 2.30-1. Require the final
`PASS: qm_289_counterexample` message and no errors anywhere in the log.
The whole verification runs inside one procedure, so a failed assertion
prevents that final message. The script does not overwrite the saved log.
No network access, external packages, or downloaded Fourier coefficients
are required.

The input is the unique degree-16 newform orbit in weight 3, level 289,
with primitive character of conductor 17 and order 16, labeled
[289.3.e.o in the LMFDB](https://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/289/3/e/o/).
The script chooses the character with chi(3)=zeta_16 and selects the unique
relative-degree-2 constituent over its character field. It constructs the
presentation E=Q(zeta_16,b), b^2=22*zeta_16, a_2=zeta_16^7, a_3=b.
It does not rely on the ordering of Magma's newform constituents or on the
name of a number-field generator.

The verification has five parts:

1. **Exact inner twists.** All 16 automorphisms
   zeta_16 -> zeta_16^r, b -> s*zeta_16^((r-1)/2)*b are tested against the
   character chi^((r-1)/2+8*(1-s)/2), through q^89. The common-level Sturm
   bound is 76. Together with the twist-level bound this certifies the
   identities globally and gives rational twist-invariant field. The
   script does not call Magma's numerically implemented `InnerTwists`.
2. **Quaternion and sign cocycles.** The scalar alpha has square
   chi(g)*22^e(g), where (-1)^e=chi_17. Its positive cocycle part contributes
   (17,22), ramified at 11 and 17. The sign class gamma_chi is ramified at
   17 and infinity, so the descent algebra is B=(-22,-1), ramified at 11
   and infinity. The script checks these symbols and the identity
   gamma_chi=(17,-1)+(-3,-1)+(17,-3), as well as the multiplication
   cocycles of alpha and I^e J^h, with I^2=-22, J^2=-1, IJ=-JI and
   (-1)^h=chi_-3. The description of gamma_chi in terms of local parity is
   a cited theorem, not a claim inferred from a prime sample.
3. **A concrete Artin factor.** The ray class group of Q(sqrt(-3)) modulo
   17*p_3 is C_96. The script chooses an order-32 character theta for which
   det(Ind theta)=chi^-1. It checks that identity on all units modulo 51,
   and checks theta^conj/theta=chi_17 restricted to the quadratic field on
   a ray class generator. The induced representation has projective
   field Q(sqrt(17),sqrt(-3)) and projective group D_2=C_2 x C_2. The ray
   groups for the proper divisor moduli have orders 48 and 1, certifying
   conductor 17*p_3 for theta and conductor 2601 for its induction.
   The powers of a chosen ray class generator printed in the log can
   depend on Magma's generator convention; the defining character
   conditions are checked independently of that convention.
4. **Full Euler polynomials.** For all 44 primes p<200 other than 3 and
   17, the script tensors companion matrices for rho_f(Frob_p) and
   tau(Frob_p). It checks all five Euler coefficients, rational
   integrality, the orthogonal reciprocal relations, and the component
   criterion r=4 exactly when both 17 and -3 are squares modulo p.
   These are consistency checks for the constructed factors, with no
   surface point counts as input.
5. **Local twist obstruction.** On mu_16 in Z_17^times, multiplying chi
   by an arbitrary character square preserves an odd exponent and hence
   order 16. This finite character-group calculation is checked exactly.
   Since nebentypus values lie in the coefficient field, every Dirichlet
   twist has coefficient field containing Q(zeta_16), hence degree at
   least 8. This includes characters with additional ramification.

The global descent argument is mathematical, not a conclusion drawn from
part 4. Ribet--Momose descent puts rho_f(g)/alpha(g) in B tensor Q_l.
The sign-cocycle identity in part 2 supplies a sign cochain u such that
u(g)*alpha(g)^(-1)*I^e(g)*J^h(g) is a finite-image representation tau.
The tensor then acts on B tensor Q_l by left and right multiplication;
B_l tensor B_l is M_4(Q_l), including at l=11. Thus the tensor has an
actual four-dimensional Q_l-model for every l, with reduced-norm pairing,
Hodge--Tate weights 0,0,2,2, and QM connected monodromy. The concrete
induction in part 3 differs from this tau at most by a quadratic twist,
which preserves rational descent.

This is a counterexample to a coefficient-degree-at-most-four assertion
based solely on a rational rank-four QM compatible system. It does not
identify a maximal surface whose entire transcendental system is this
one. A detailed proof is included in the Codex comment in the separately
maintained manuscript draft dated 20 November 2026.

References:

- A. J. Best et al., *Computing classical modular forms*, in *Arithmetic
  Geometry, Number Theory, and Computation*, Simons Symposia (2022),
  [Lemma 11.2.1 and Theorem 11.2.4](https://jvoight.github.io/articles/cmfs-mf-101222-scrubbed.pdf):
  twist levels and exact inner-twist certification.
- E. Ghate, E. Gonzalez-Jimenez, J. Quer, *On the Brauer class of modular
  endomorphism algebras*, IMRN 2005, no. 12, 701-723,
  [Theorem 4.1 and proof of Theorem 4.3](https://arxiv.org/pdf/2602.13401)
  ([published DOI](https://doi.org/10.1155/IMRN.2005.701)):
  Brauer formula in all weights and local parity of the sign class.
- D. Banerjee, E. Ghate, *Adjoint lifts and modular endomorphism algebras*,
  Israel J. Math. 195 (2013), 507-543,
  [Section 3, especially Proposition 8](https://maths-people.anu.edu.au/~banerjeed/adjointdebarghaghate.pdf):
  the scalar alpha, its trace relation, and its cocycle.
- F. Momose, *On the l-adic representations attached to modular forms*,
  J. Fac. Sci. Univ. Tokyo Sect. IA Math. 28 (1981), 89-109;
  see also J. Lang, *On the image of the Galois representation associated
  to a non-CM Hida family*, Algebra Number Theory 10 (2016), 155-194,
  [Theorem 6.1](https://msp.org/ant/2016/10-1/ant-v10-n1-p06-p.pdf),
  for the quaternion-valued representation on the inner-twist kernel.
