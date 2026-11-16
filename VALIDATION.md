# Import validation — 15 November 2026

The import was checked in a disposable copy using Magma 2.30-1, a C++17
compiler, and Python 3. All twelve computation logs contained the expected
PASS messages and no Magma runtime, assertion, syntax, or user errors.

| Check | Result |
| --- | --- |
| Congruence enumeration | 64 classes and 558 actual subgroups; independent cusp and level checks passed. |
| Standard congruence examples | All six quadratic refinements and Γ₁(9), Γ₁(10) matched. |
| Lifts of projective congruence groups | 1,745 actual subgroups, including 1,187 noncongruence subgroups. |
| Rational-lift enumeration and identification | All 64 geometric classes checked at exact and proof moduli; every computed lift matched a unique LMFDB fine label. |
| Adelic obstructions and higher-level checks | The three exceptional classes passed the obstruction tests and the doubled/quadrupled-level checks. |
| Congruence comparison actions | Projective classes and all 64 signed actions regenerated. |
| C++ exhaustive enumeration | Compiled with `-O3 -std=c++17 -Wall -Wextra -pedantic`; regenerated `maps.txt` and `summary.csv` byte for byte. |
| Independent noncongruence audit | 47,244 total classes, 47,180 noncongruence classes; Hamilton–Loeffler and Burnside checks passed. |
| Explicit surfaces | Γ(4) twist and Calegari example checks passed. |

The three Python exporters, and reconstruction of the LMFDB cache from its
saved API pages, reproduced all checked-in generated data byte for byte.
Both copies of the upstream Cummins–Pauli tables matched their recorded
SHA-256 checksums. Fresh computation logs are saved in the corresponding
directories.

The numerical continuation used for the saved Calegari row identification
was not rerun during this import; its existing action data and provenance
are retained, and the exact comparison with the complete census was rerun.
The historical surface summary files were preserved without mathematical
reverification.

## Additional surface and candidate code

The subsequently supplied `H20/Code/` and `newvers/202-candidates/` folders
contributed 26 and 7 original files respectively. All 33 imported files
matched their source SHA-256 checksums. Editor backups and macOS metadata
were excluded. Each destination directory contains its import manifest.

The two candidate-package checks were rerun with Magma 2.30-1:

| Check | Result |
| --- | --- |
| `202-candidates/scripts/qm_candidates.m` | All nine rational rotation groups, element-order counts, and symbolic tensor, Asai, and trace identities passed. |
| `202-candidates/scripts/example_S1_9.m` | All 25 complete Euler polynomials matched at primes at most 101 other than 3, including 2. |

Both complete logs contained the expected `PASS` line and no runtime,
assertion, syntax, or user errors. Their output matched the supplied logs
exactly after normalizing the Magma version string. Fresh logs are saved
under `202-candidates/results/validation-magma-2.30-1/`; the original logs
remain unchanged.

The legacy surface calculations were preserved without mathematical
reverification. Their README records known missing dependencies and
unfinished scripts; they are not claimed to form a passing test suite.

## Reconstructed manuscript checks

On 15 November 2026, the missing `herfurtner_checks.m` and
`modular_dimensions.m` were reconstructed in `verification/` from the v11
claims and the preserved Herfurtner models. Both passed with Magma 2.30-1,
with one final PASS marker per log and no reported errors.

- Herfurtner: five characteristic-zero models have the asserted fibers,
  second Betti number 34, and trivial-lattice rank 30. All 28 diagnostic
  entries for H13/H15/H25/H28 match at the seven stated primes. H7's full
  quartics at 13 and 37 match and are irreducible with Galois group D4.
  All 30 full polynomials are recorded, including their orientation signs.
- Modular dimensions: independent dimension-formula and full modular-symbol
  computations agree at levels 1 through 40. Explicit projective coset
  orbits reproduce the cusp widths for Gamma_1(9), Gamma_1(10), and
  Gamma_1(12), including SL2 index 96, ten cusps, and dimension 3 at level 12.
  All six quadratic refinements have the asserted cusp and new dimensions.

The scripts, explicit models, scope, and reproduction instructions are in
[`verification/README.md`](verification/README.md); fresh logs are under
`verification/results/`. This reconstruction verifies the stated finite
checks. It does not identify the noncyclic Artin factors or prove a
modularity claim from the Herfurtner diagnostics.
