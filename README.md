# H20

Computational code and data accompanying *Maximal surfaces of geometric
genus 2 and modular forms* by Asher Auel, Adam Logan, and John Voight.

## Contents

| Directory | Computation and saved results |
| --- | --- |
| [`congruence-enumeration/`](congruence-enumeration/README.txt) | 64 congruence classes in SL₂(ℤ), comprising 558 actual subgroups; signed cusp data and finite matrix generators. |
| [`rational-lifts/`](rational-lifts/README.txt) | 61 geometric classes admitting determinant-surjective adelic lifts; 457 finite lifts at the exact geometric levels and 636 at the proof moduli. |
| [`noncongruence-enumeration/`](noncongruence-enumeration/README.txt) | All 47,244 classes, including 47,180 noncongruence classes; permutation actions, an independent Magma audit, and explicit examples. |
| [`historical-data/`](historical-data/README.md) | Original numerical summary tables for quadratic twists and double and triple covers. |

The enumeration directories include source data, machine-readable results,
and logs from the original computations. The main manuscript and its CM
companion are maintained separately.

## Requirements and reproducibility

- Magma (the saved computations used version 2.30-1).
- A C++17 compiler for the exhaustive map enumeration.
- Python 3.9 or later, using only the standard library, for data exports.

Run Magma and C++ commands from the indicated subdirectory. On the local
Apple Silicon installation, use `MAGMA_ARCH=arm64 magma` in place of
`magma`. Require the expected `PASS` messages and check that no errors appear:
Magma can return exit status zero and continue after reporting an error.

Each directory's README lists the commands and their outputs. For a full
regeneration, work in this order:

1. Regenerate the congruence census and its JSON/CSV exports.
2. Regenerate and identify the rational lifts, then export their tables.
3. Export the congruence comparison actions, enumerate all oriented maps,
   audit the signed lifts in Magma, and export the complete census.

The later computations read the single canonical generator list at
`congruence-enumeration/representatives.m`. Regenerating it may change
representative choices or row ordering, so regenerate dependent data
together. No network access is needed with the included caches.

The raw mathematical computations are distinct from `record_results.py`:
those Python scripts export and check saved output; they do not replace
the Magma or C++ calculations. The scripts write their generated data in
their own directories. Run checks in a disposable copy when you want to
preserve the saved results unchanged.

## Sources and scope

The Cummins–Pauli genus-zero and genus-one congruence tables were downloaded
on 15 November 2026 from
[their congruence subgroup tables](https://mat112.uncg.edu/tables/congruence/).
The unchanged extracts and their SHA-256 checksums accompany both
enumerations. The complete noncongruence test uses Hamilton–Loeffler
relations independently of those tables.

The rational-lift directory includes public LMFDB fine modular-curve records
retrieved on 15 November 2026. `fetch_lmfdb.py` rebuilds the consolidated
cache from the saved API pages, retrieving missing pages when necessary.
Use a fresh cache directory to obtain a new snapshot.

Calegari's example is attributed in its code and data. The supplied
unpublished note itself is not included. The historical surface tables
contain preliminary matches; their inclusion does not assert a proved
identification of compatible systems.

The repository's [MIT license](LICENSE) applies to the project code.
Upstream mathematical data retain their original attribution.

## Citation

John Voight, *H20: Code for maximal surfaces of genus 2 and their modular
forms*, software and data repository, 2026,
<https://github.com/jvoight/H20>.

See [`CITATION.cff`](CITATION.cff) for machine-readable citation metadata.
For reproducibility, also record the commit used.
