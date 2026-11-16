COMPUTATIONAL SUPPLEMENT TO SECTION 3.3

Congruence census: 64 SL_2(Z)-conjugacy classes, comprising 558 actual
subgroups; 32 classes in each base genus 0 and 1.
Partial noncongruence count: 1187 actual noncongruence subgroups with
congruence projective image, not a complete noncongruence census.

From this directory run (Magma 2.30-1 used for the saved logs):
  magma -b enumerate_congruence.m
  magma -b verify_examples.m
  magma -b noncongruence_lifts.m
  python3 record_results.py
On the local Apple Silicon installation, prefix with MAGMA_ARCH=arm64.
No additional package or network access is required. Check the final PASS
line in each run; Magma can return a zero exit status after an error.

The output representatives.m gives finite matrix generators modulo twice
the projective level. Their full inverse images define the subgroups.
The projective index mu in the computational files equals n/2 in the
manuscript. The full list is in results.csv/json and representatives.m.

The first script regenerates representatives.m, results.txt, and audit.txt;
the third script reads audit.txt. The Python exporter regenerates JSON/CSV from results.txt; generator choices and row ordering can change on rerunning Magma.

Original Cummins--Pauli data, unchanged extracts from:
https://mat112.uncg.edu/tables/congruence/csg24.dat.tar.gz
https://mat112.uncg.edu/tables/congruence/
Downloaded 15 November 2026. SHA-256 checksums: source-data-hashes.json.
Lift-level theorem: Kiming--Schuett--Verrill, Proposition 3,
https://doi.org/10.1112/jlms/jdq062 (arXiv:0905.4798).
