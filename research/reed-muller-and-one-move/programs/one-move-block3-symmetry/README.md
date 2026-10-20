# One move at block 3 under symmetry (tier 2)

What these programs show: RM14 and the refuted reach of section 4, item 4. `lib3.py` is the shared library; `s1_levi.py`
to `s14_fixedmove_check.py` are the screens (Levi subgroups, character splits, lifts, the fixed-move rule, the
21-fold group, the permutation symmetry of the three address bits and its six 5-in-14 shapes). The skeptic's programs
are in `independent-check/` (`vlib.py`, `v1_basic.py` to `v7_traces_gauge.py`); its gauge-symmetry table, which shows
the surviving 20 in 58, 20 in 59 and 15 in 42, is `../../data/one-move-block3-symmetry/independent-check/out_v4_gauge_table.txt`.
The proofs are in `../../proofs/one-move-block3-symmetry.md`; their exclusions hold under exact symmetry only.

Each screen runs in seconds to a minute: `python3 s5_table.py > out_s5_table.txt`, and so on; all 13 were re-run for this revision
and print the same results (s7 and s14 also print numerical warnings from numpy, which the recorded files hold
mixed into the output).

Run everything from a working copy made by `sh ../workdir.sh` (from this folder), inside `work/one-move-block3-symmetry/`, where the data files of `../../data/one-move-block3-symmetry/` sit next to the programs. Python programs need Python 3.9 or later and numpy.

## Files

Each line is the file's own description (from its header), lightly shortened; the files themselves keep the working labels of the rounds, which the section README's Sources maps to these folders.

| File | What it does |
|---|---|
| `independent-check/v1_basic.py` | Levi lifts: P_A * diag phase * Pauli, free, with K^-1 g K free Singer: A of order 7 |
| `independent-check/v2_f21.py` | all elements g = P_A D_S X^a Z^b (B mod scalars) that are in the Levi and satisfy g c g^-1 = c^2 exactly (up to scalar 1) n-fixed cells at distance 1 |
| `independent-check/v3_cells.py` | exact commuters in the second cell, chunked over B |
| `independent-check/v4_gauge_table.py` | Z7 per-character screens (the tier-2 block-3 track's own: hook 80 m_k >= 224 a + 32 b; Theorem D trivial m_0 >= a+b+2; non-trivial m_k >= a+2 (a>0); depth m_k >= 3 when any piece is present) in ... |
| `independent-check/v5_s3.py` | helper (no description in the file) |
| `independent-check/v7_traces_gauge.py` | gauge-symmetric single move, f = 1 and f = 2: Z7 regular orbit of arrays, c Singer, c' = K^-1 c^-1 K |
| `independent-check/vlib.py` | Block-3 engine letters as 8x8 complex matrices; symplectic images on labels (x/z); Lagrangians of F2^6; the f = 5 hook module V (x) Lambda^4 V with Y = N = span{e_x (x) e |
| `lib3.py` | block-3 engine letters as 8x8 complex matrices, symplectic images, membership in the free group B (monomial Clifford letters), the Levi group B cap K B K^-1, and the f = 5 hook module M = V (x) ... |
| `s10_f21.py` | the F21 refinement of the trivial-character problem of the Z7 survivor 15 in 42. F21 = <c, n>, c the Singer lift (order 7), n an order-3 Levi element with n c n^-1 = c^2 (s7_equiv.py). Array side: ... |
| `s11_f21_table.py` | F21-symmetric block-3 shapes (PROOF.md S9). Array orbits: alpha of size 21, beta of size 7 (stabiliser Z3), gamma of size 3 (stabiliser Z7), delta fixed. Move orbits: r of size 21, q of size 7 ... |
| `s12_s3.py` | S3 = the coordinate permutations (an exact symmetry of K and of the letters) acting on W in M at block 3, per irreducible representation of S3 (PROOF.md S10). Array orbit with stabiliser S ... |
| `s13_fixedmoves.py` | robustness of the fixed-move cap splits used by s5/s11/s12. For T = <transposition (01)>, <3-cycle>, S3 (coordinate permutations) and <n> (the F21 element), collect every move letter m = b dir z ... |
| `s14_fixedmove_check.py` | calibration of PROOF.md S4 for fixed moves. c = 3-cycle (commutes with K, so c' = c^-1); a fixed array y restricted to N_k = { y in N : c' y = zeta^k y }; the raw input register is tau-fixed and ... |
| `s1_levi.py` | the Levi group B cap K B K^-1 at block 3, its Singer element, and the per-character hook demands of the input-side Singer action on N/E_N (f = 5). Also calibrates the hook numbers (cap 80, demand ... |
| `s2_split.py` | eigen-splits that decide the per-character hook test of symmetric block-3 shapes. For a symmetry element c of exact order n acting through the Levi (c and c' = K^-1 c^-1 K in B): - demand split: ... |
| `s3_decomp_check.py` | numerical calibration of the per-character decomposition (PROOF.md Lemma S1). A random Z7-symmetric program (Singer lift c acting through the Levi, c'_h = K^-1 c_h^-1 K) is run on random inputs in ... |
| `s4_lifts.py` | demand and fixed-move cap splits for EVERY Pauli-twisted lift c P of the 3-cycle and of the transposition (P runs over the 64 Paulis; c P is rescaled to exact order n when (c P)^n is a scalar), ... |
| `s5_table.py` | which symmetric block-3 shapes W in M survive the per-character screens (PROOF.md S1-S5). Group H of order n acting through letters (Levi lifts for fixed arrays). Shape data: a regular array ... |
| `s6_rf.py` | the filtered rung (ladder Theorem RF) in the f = 5 hook module, for full arrays and for Singer-projected pieces, to screen the per-character programs of the Z7 survivors. V_0 = N, V_(j+1) = ... |
| `s7_equiv.py` | checks behind PROOF.md Lemma S6 (only two distinct per-character problems for the Singer). (a) an order-3 Levi element n with n c n^-1 = c^2 exactly (so n' = K^-1 n K is free and n' c' n'^-1 = ... |
| `s8_hyperplanes.py` | PROOF.md S3' says a symmetry element can fix a move only if its linear part on L_Z fixes a hyperplane. Count, for every element of GL(m,2) (m = 3, 4), by order and conjugacy-invariant data (order, ... |
| `s9_levigroup.py` | do the 168 canonical Levi lifts c_A = P_A diag(i^q_A) form a group (exactly, or up to scalars / Paulis)? Prints the set of correction factors c_A c_B c_AB^-1. |
