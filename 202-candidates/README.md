# Candidate groups and Euler-polynomial checks

This directory contains the working package supplied as
`newvers/202-candidates/`, dated 23 October 2026. All seven supplied files
are preserved unchanged; `import-sha256.json` records their checksums at
import on 15 November 2026.

## Reproduction

Run from this directory with Magma on your path:

```sh
magma -b scripts/qm_candidates.m
magma -b scripts/example_S1_9.m
```

On the local Apple Silicon installation, use `MAGMA_ARCH=arm64 magma`.
Both commands load `scripts/euler_tools.m` using a relative path. Neither
requires the other repository directories, network access, or external
data. Check the entire output for errors as well as the final `PASS` line.

- `qm_candidates.m` verifies nine rational quaternionic rotation groups,
  element-order counts, and symbolic tensor, Asai, and trace identities.
- `example_S1_9.m` compares all five coefficients of the geometric and
  modular Euler polynomials at the 25 primes at most 101 other than 3,
  including 2. The bound is the argument of `RunS19(101)` at the end.
- `euler_tools.m` supplies the exact polynomial utilities used by both.

The original logs in `results/` were produced with Magma 2.29-1. Fresh
Magma 2.30-1 checks are stored in `results/validation-magma-2.30-1/`.

## Notes and scope

`classification-notes.tex` and `README.txt` are the original working
documents. They describe the package at its creation; remarks there about
an unavailable code archive or temporary local setup are historical.
The classification note is not the current manuscript or CM companion.

The rotation groups are algebraic candidates, not asserted Galois or surface
realizations. The S₁(9) computation is finite Euler-factor evidence and does
not by itself prove an identification of compatible systems. The script
uses the known characteristic-zero Néron–Severi contribution; it does not
discard extra factors caused by Picard-rank jumps after reduction.
