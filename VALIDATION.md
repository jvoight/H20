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
