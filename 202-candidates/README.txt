Working package: (2,0,2) motives, 23 October 2026

CONTENTS
classification-notes.tex: endomorphism definitions and comparison hypotheses;
  non-CM compact candidate classification; rational quaternion embedding
  conditions; initial realization-status table; CM boundary; evidence policy.
scripts/euler_tools.m: exact rank-four orthogonal, tensor and Asai helpers.
scripts/qm_candidates.m: nine rational quaternionic rotation groups, component
  element counts, and a symbolic tensor-polynomial identity.
scripts/example_S1_9.m: independent full Euler comparison for the S_1(9)
  elliptic surface and the weight-three modular form of level 9.
results/: logs from Magma 2.29-1.

RUN
From this directory, with Magma on PATH:
  magma -b scripts/qm_candidates.m
  magma -b scripts/example_S1_9.m

The example uses all primes <= 101 except 3, including 2. Change RunS19(101)
at the end of the script to increase the bound. It recomputes modular symbols
and finite-field elliptic-surface L-functions. It does not need LMFDB access,
an external data file, or the authors' existing code archive.

VERIFIED
All 25 complete Euler polynomials agree, including their constant, linear,
quadratic, cubic and quartic coefficients. The tensor description uses
tau = 1 direct-sum chi^(-1), with chi a generator of the order-six character
group modulo 9. The corresponding projective Artin image is C6.
The nine rational rotation-group orders and element-order counts pass, as
does the symbolic 2-by-2 tensor characteristic-polynomial identity.

LIMITS
This is one pilot, not a conversion or audit of all existing examples.
The supplied Code.zip could not be downloaded in this session. No existing
manuscript or code was overwritten. Other candidate-table rows remain
unassigned in this audit; this does not assert nonexistence of examples.
The rational quaternion algebra and explicit geometric correspondences for
the pilot remain to be determined. The CM arithmetic list is not completed.
The scripts provide numerical evidence, not a proof of a compatible-system
identification or of algebraicity of Hodge endomorphisms.

LOCAL TEST NOTE
The installed Magma needed to compile a missing modular-symbols package data
file. Tests used an unchanged writable copy of its standard package directory
in temporary storage. No installed Magma source files were changed. A normal
Magma installation with its package data available runs the commands above.

EVIDENCE POLICY
Do not remove every root-of-unity-times-p factor from a reduction: Picard
rank can jump. Remove only the known specialization of the characteristic-
zero Neron--Severi part. LFunction(E_p) removes the trivial lattice; known
characteristic-zero Mordell--Weil contributions need separate removal.
The S_1(9) pilot has characteristic-zero geometric MW rank zero.
Keep geometric construction, group/endomorphism identification, and modular
Euler-factor matching as distinct status fields.
