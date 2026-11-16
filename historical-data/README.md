# Historical surface tables and construction notes

These are the original numerical summaries accompanying the surface
constructions in the paper:

- `summary-quadratic-twists.txt`: twists of Herfurtner's rational elliptic surfaces.
- `summary-twisted-double-covers.txt`: twisted double covers.
- `summary-triple-covers.txt`: triple covers.
- `summary`: preliminary search by degree of the modular j-map, including
  unsuccessful cases and unexplained point counts.
- `k3-double-covers.txt`: constructions and numerical matches for twists of
  singular K3 fibrations.
- `double-covers.txt`: degree-10 double covers of P(1,1,2), including the
  cuspidal sextic and five-conic constructions.

The files preserve the historical labels and base-surface numbering. That
numbering is not consistent between the three cover/twist tables or with the paper's
Herfurtner table. Some form labels and twists remain to be audited; these
tables are records of computations, not proofs of the proposed arithmetic
identifications. Consult the paper's census and realization-status discussion
for the current qualifications.

## Reading the records

For the double-cover table, identify a row by filename and `(i,n)`.
For the triple-cover table, also specify its section: total, mixed, or
partial ramification. A repeated form label or a repeated row does not
establish an isomorphism of surfaces or compatible systems. In particular,
these tables alone do not certify counts of distinct systems or twist classes.

The three construction notes added on 15 November 2026 are byte-for-byte
copies from the manuscript directory; `notes-import-sha256.json` records
their source paths and checksums. The notes are historical records,
not a newly verified computational supplement. They contain instructions
and suggestions from the original research process; these are retained as
source text, not current reproduction instructions.

The main manuscript (`H20-modfrms.tex`, formerly v11) incorporates selected
constructions and qualifications:

- The K3 note's level-512 associations at discriminant -8 and its level-1800
  association at discriminant -12 are numerical records, not proved
  compatible-system identifications. The level-128/256 discrepancy at
  discriminant -4 remains unresolved.
- The Legendre cubic twist `t(t-1)y^2 = x(x-1)(x-t^3)` is an explicit
  maximal example with no form identified in the note.
- In the five-conic note, the listed singularities of the `t=6` example,
  `2D6 + D4 + 3A3 + 2A1`, have total root rank 27, not the asserted 28.
  That list and the canonical class alone do not prove maximality.
  The recorded list itself still requires verification.
- The cuspidal-sextic/quartic construction and the level-200 five-conic
  association remain candidates pending geometric and arithmetic checks.

The full exploratory tables remain here to keep the manuscript compact.
