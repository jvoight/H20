# Original surface calculations

These 26 files were imported unchanged on 15 November 2026 from the
supplied `H20/Code/` directory. `import-sha256.json` records their original
checksums. The directory structure, saved output, and mathematical comments
are preserved; macOS metadata and editor backup files were omitted.

These are working calculations and research notes, not a single executable
verification suite. They have not been mathematically reverified during
this import. Some expect variables from an interactive Magma session, such
as `Q := Rationals()`, and some contain unfinished code. The directory's
original filenames and example numbers need not match the current paper.

| Files | Contents and status |
| --- | --- |
| `maximal.m`, `maximal-summary.txt` | Herfurtner models, quadratic twists, and preliminary modular-form matches. `maximal.m` contains a missing list separator after the fourth model. |
| `lfs.m`, `lfs.out` | Exploratory analysis of Euler factors and saved output. The source expects earlier session data and contains an unfinished conditional and uncommented working notes. |
| `singular-k3/` | Fibrations on singular K3 surfaces, arranged by discriminant; see the original `readme` there. `rank20-k3.mag` is an unfinished list. |
| `good_twists.mag` | Helpers for choosing quadratic twists using fiber data and discriminants. |
| `x0_15.mag`, `count_nr.mag`, `pols_15.mag` | Calculations over X₀(15), point-counting helpers, and saved Euler polynomials. The relative load of `count_nr.mag` assumes this directory is the working directory. |
| `steiner.m` | Steiner-quartic and other double-cover calculations. Later calculations load `k3pic.mag`, which was not supplied. |
| `20.3.f.a.fail.mag` | Exploratory obstruction and twist calculations. Requires `test_j.mag` and a local `ModFrmGL2` package, neither supplied here. |
| `h02.m`, `hgm-212.m` | Hodge-number and exploratory hypergeometric calculations; these are working fragments rather than standalone checked examples. |

The missing dependencies have not been replaced with guessed substitutes.
For the self-contained quaternionic and S₁(9) checks, see
[`../202-candidates/`](../202-candidates/README.md). For previously collected
surface summary tables, see [`../historical-data/`](../historical-data/README.md).
