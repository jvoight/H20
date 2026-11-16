RATIONAL LIFTS: VERIFIED RESULTS, 15 NOVEMBER 2026

61 of 64 geometric congruence classes admit determinant-surjective lifts.
The three failures are IDs 2, 6, 7 (CP 4D0, 5F0, 6A0), even at adelic level.
There are 457 GL2-conjugacy classes at the exact geometric levels, and
636 at the corresponding proof moduli 2Nprime. The latter include the former.

labels-by-geometric-group.csv gives ALL labels, grouped by geometric class.
all-lifts.csv has one row for each of the 636 distinct fine LMFDB labels.
all-lifts.json includes computed generators, LMFDB generators, and conjugators.
geometric-groups.json includes the geometric generators and cusp signatures.

Run locally from this directory (Magma 2.30-1 was used):
  MAGMA_ARCH=arm64 magma -b enumerate_lifts.m
  MAGMA_ARCH=arm64 magma -b identify_lifts.m
  MAGMA_ARCH=arm64 magma -b obstructions.m
  MAGMA_ARCH=arm64 magma -b higher_lifts.m
  python3 record_results.py
The cached LMFDB Magma data suffice; no network access is needed to verify.
fetch_lmfdb.py refreshes public API records when desired.

The genus-zero and genus-one counts of successful geometric classes are
29 and 32, respectively. These are geometric SL2(Z) conjugacy classes;
fine labels count finite GL2 conjugacy classes. Counts at larger unrestricted
adelic levels are not claimed to be finite.

The Magma scripts read ../congruence-enumeration/representatives.m.
The exporter also reads that directory's results.json; regenerate the
congruence census first if its representatives have changed.
