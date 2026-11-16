COMPLETE NONCONGRUENCE CENSUS -- 15 NOVEMBER 2026

47244 SL2(Z) conjugacy classes in total, of which 64 are congruence and
47180 noncongruence. Actual subgroup counts: 1466010, 558, 1465452.
Noncongruence by base genus: 44071 genus zero, 3109 genus one.
All 64 congruence classes were matched individually to the earlier census.

Build/run from this directory:
  MAGMA_ARCH=arm64 magma -b export_congruence.m
  c++ -O3 -std=c++17 enumerate_maps.cpp -o enumerate_maps
  ./enumerate_maps
  MAGMA_ARCH=arm64 magma -b verify_census.m
  python3 record_results.py
The cached congruence_actions.txt is only used as a preliminary sanity check;
the actual congruence classification uses Hamilton--Loeffler relations on
every signed action and is independent of the Cummins--Pauli tables.

maps.txt and maps.json encode the projective pairs and lift bits.
subgroups.csv/json identify all 47244 lifted classes and their signed cusps.
R and S in maps.json are zero-based images. S is the image of [0,1;-1,0],
R that of S*[1,1;0,1]. For each free edge i, a lift bit changes both S-signs;
the baseline S-sign is 0 at i<S(i) and 1 at i>S(i), and all R-signs are 0.
Thus a signed action can be recovered without any group-presentation search.

Calegari's explicit surface is rigorously within the census by its geometry;
its individual row identification uses numerically continued monodromy,
repeated with two path resolutions, followed by exact permutation comparison.
The supplied private note is cited but not redistributed in this package.

Additional example checks (run from this directory):
  MAGMA_ARCH=arm64 magma -b gamma4_example.m
  MAGMA_ARCH=arm64 magma -b calegari_example.m
The Gamma(4) example agrees, up to constant twist by -1, with row (1,15)
of the earlier twisted-double-cover computation.
The CP source tables are included. export_congruence.m reads the original
64 representatives from ../congruence-enumeration/representatives.m and
regenerates the comparison data.
