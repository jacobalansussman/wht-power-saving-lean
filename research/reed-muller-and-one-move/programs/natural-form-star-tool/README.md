# The natural form at block 5: the star tool

What these programs show: RM7. `xbase.py` (field, forms, canonical Singer targets in a self-dual basis) and
`gen_tables.py` (writes `tables3.bin`, `tables4.bin`, `tables5.bin`) set up the system; `scanrc.c` runs deletion,
the ten-piece rank cover and the one-level star test; `starrec.c` and `starrec2.c` run the star recursion with pooled
zero blocks and the clique and matching tests; `rcdp.c`, `rc_run.py` and `rc_ilp.py` are the general rank cover;
`classes.py` enumerates the classes of target sets (Burnside counts 50,895, 130,065, 286,155); `cal2.sh` runs
`starrec2` on lines `m k1,k2,... label`; `merge8.sh` assembles the 8-in-31 verdict file `all8.txt` from the job
outputs. Second programs written by the same track, sharing no code with the scanner, are in `second-programs/`
(`vscan.py`, `vstar2.py`, `vcert.py`, `vcert2.py`). The skeptic's own field, forms, flats and recursion are in
`independent-check/` (`vf.py`, `vrec.py`, `vgen.py`, `vlines.py`, `coverage.py`, `vcert8a.py`, `vcert10a.py`,
`vcert2me.py`, controls `ctrl.py`, `ctrl2.py`, `u3true.py`).

Commands re-run for this revision:

```
cc -O2 -o starrec2 starrec2.c
python3 gen_tables.py 3; python3 gen_tables.py 4; python3 gen_tables.py 5     # seconds; byte-identical to the originals
printf '5 0,2,4,8,18,24,25,27 first-set\n5 0,1,2,3,4,5,6,7 a\n5 0,1,2,3,4,5,7,8 b\n5 0,1,3,4,5,6,7,8 c\n4 0,1,2,4 b4\n' | ./cal2.sh
      # about 6 minutes in all: the first set is excluded by the clique test (11 > 8), the three sampled classes
      # get exactly the verdicts recorded in all8.txt, and the block-4 class is excluded
(cd independent-check && python3 coverage.py)                                  # 10 s, reads jobs/: identical to out_coverage.txt
```

A class that passes (a plain control, or one of the 1,266 survivors) takes several minutes, because every form is
tried. Output lines: `S ks P` passes, `D` deletion, `Z` star recursion, `C` clique, `M` matching, `X`/`Y` the rank
cover and one-level star test of `scanrc`.

Not re-run: the full scans (`scanrc` and `starrec2` over every class, many hours on several 8-processor Linux
machines). Their outputs are in `../../data/natural-form-star-tool/job-outputs.tar.gz` (the work copy unpacks them
into `jobs/`), with the merged first-pass verdicts in `all8.txt` (1,669 classes pass there; the deep pass, jobs deep to deep3,
excluded 403 of them) and the 1,266 survivors in `surv8_final.txt`.

Run everything from a working copy made by `sh ../workdir.sh` (from this folder), inside `work/natural-form-star-tool/`, where the data files of `../../data/natural-form-star-tool/` sit next to the programs. Python programs need Python 3.9 or later, numpy, and scipy for the rank-cover ILP.

## Files

Each line is the file's own description (from its header), lightly shortened; the files themselves keep the working labels of the rounds, which the section README's Sources maps to these folders.

| File | What it does |
|---|---|
| `cal2.sh` | run starrec2 on "m ks [label]" lines from stdin (test forms = every class with r = 0 for these targets) |
| `classes.py` | orbit representatives of W-subsets of Z_n (canonical Singer copy indices k, target mult(alpha^2k)) under k -> 2^j k + c (Singer cycle and Frobenius; they act on the canonical targets and fix the ... |
| `dump_cert.py` | <case> <F> [dmax]: run rcdp on one form and write its certificate (raw rcdp text) to stdout. |
| `gen_tables.py` | tables for scanrc.c (block 5, one move per direction, canonical Singer targets T_k = mult(alpha^2k)). Writes tables5.bin (little-endian int32): nf, then mask[nf] (bit z-1 <-> direction z), wt[nf], ... |
| `independent-check/cmp.py` | <expect.txt> <out files...>: compare the star-tool track verdicts with vrec verdicts (letter and, for Z, round). |
| `independent-check/coverage.py` | own orbit canonicalisation (min 31-bit mask over all 155 maps k -> 2^j k + c), Burnside counts, and coverage / verdict bookkeeping of the the star-tool track output files. |
| `independent-check/ctrl.py` | genuine plain programs (canonical frames, disjoint) at block 5 with one move per direction; build the TRUE E = I + Gamma - K of the program (Lemma 1 with identity mixings), check T10 on all 32,767 ... |
| `independent-check/ctrl2.py` | all frames (5 directions z with sum z z^T = T_k) by brute force; disjoint frame packings for the plain controls; TRUE E of each plain program; T10 on all forms; vrec.run with trueE (u(A) >= true ... |
| `independent-check/u3true.py` | a genuine slice unit U_3 2 in 7 (all 7 directions once, targets I, I) found numerically by an earlier families round's solver (refl2.System2, as in t_theorem10.py); its E = I + H - K (refl ... |
| `independent-check/vcert10a.py` | helper (no description in the file) |
| `independent-check/vcert2me.py` | own line-by-line check of starrec2 CERT=1 output lines (C, M, Z; round 1). Each cited zero cell x,y:F/S/round/d/i: r_F = 0; S <= D_F; if d >= 0: S = D_F n A for a d-flat A on which OWN round-1 ... |
| `independent-check/vcert8a.py` | line-by-line check of a starrec2 CERT clique/matching certificate with own code (vf + vrec flats). For each cited cell x,y:F/S/round/d/i : r_F = nu (<= 0 here) for the target set; S = D_F n A for ... |
| `independent-check/vf.py` | skeptic's own GF(32) / forms / canonical-target parities, written from scratch. Field: GF(2)[x]/(x^5+x^2+1) by carry-less multiplication (no log tables). Directions z in 1..31, bit j = coordinate ... |
| `independent-check/vgen.py` | block-m version (m = 3, 4, 5) of vf.py + vrec.py, own code: field by carry-less multiplication, own self-dual basis search, canonical targets mult(alpha^2k), forms, flats, star recursion with ... |
| `independent-check/vlines.py` | verify scanrc lines (D, X ten-piece, Y one-level star) from scratch with vf.py tables. Every bound is recomputed: min over ALL nonzero forms G with D_G >= piece of s_G = /D_G/ - r_G (and the piece ... |
| `independent-check/vrec.py` | skeptic's own re-implementation of the star recursion + pooled zero blocks + clique + matching + /Y/-W tests for the per-array-target system at block 5, one move per direction, canonical targets. ... |
| `merge8.sh` | assemble the W = 8 starrec2 scan from its index slices (each slice file has one line per class, in |
| `rc_ilp.py` | piece + term-rank bound on rank E[D_F, D_F], optimised exactly by an ILP (HiGHS). Bound: E[D_F,D_F] = sum over chosen disjoint pieces R x C (each inside D_G x D_G for a form G, so rank <= s_G) + ... |
| `rc_run.py` | build a multi-target system, run rcdp (guillotine rank-cover bound) on chosen forms, and verify every reported violation INDEPENDENTLY (partition of D_F x D_F, each piece's bound re-derived with a ... |
| `rcdp.c` | guillotine rank-cover bound for the multi-target system, fixed move counts. For a form F with support D_F (a set of directions; moves = directions with multiplicity counts[z]): rank E[D_F, D_F] = ... |
| `scanrc.c` | fast ten-piece rank-cover filter for the block-5 multi-target system with one move per nonzero direction and canonical Singer targets T_k = mult(alpha^(2k)) (tables from gen_tables.py). For a test ... |
| `second-programs/vcert.py` | verify a rank-cover certificate for the multi-target system with one move per nonzero direction and canonical Singer targets T_k = mult(alpha^(2k)), alpha a root of x^5 + x^2 + 1. Everything is ... |
| `second-programs/vcert2.py` | (construction part copied from vstar2.py) : an independent re-derivation of the exclusions made by starrec2 at block 5 (one move per nonzero direction; array l owed T_l = mult(alpha^(2 k_l)), al |
| `second-programs/vscan.py` | verify every line of scanrc output. Lines: "S k1,..,kW X F phi psi bound sF G_H G_H' G_psi0 G_x0 G_x1 G_psi1" (ten-piece exclusion) "S k1,..,kW Y F phi psi c1 c2 2/Psi/ sum G1 G2 G3" (star-lemma |
| `second-programs/vstar2.py` | an independent re-derivation of the exclusions made by starrec2 at block 5 (one move per nonzero direction; array l owed T_l = mult(alpha^(2 k_l)), alpha a root of x^5 + x^2 + 1). It |
| `star_explore.py` | <k list>: for one canonical target set, the star-lemma recursion on every even form F: u(S) for S = D_F n (affine flat), dims 4 -> 0, starting from u(S) = min(/S/, s_G over forms G != F with S in ... |
| `starrec.c` | the star-lemma recursion on the block-5 multi-target system, one move per nonzero direction, canonical Singer targets T_k = mult(alpha^(2k)) (tables5.bin). For an even form F (r_F = 0), E[D_F, D_F] is |
| `starrec2.c` | star-lemma recursion + pooled zero blocks + clique and matching tests, for the multi-target system at block m (3, 4, 5) with ONE move per nonzero direction and canonical Singer targets T_k = ... |
| `vmask.py` | helper (no description in the file) |
| `x_probe.py` | helper (no description in the file) |
| `xbase.py` | own field / forms / targets code for the multi-target system (Theorem B form). Coordinates: a self-dual basis f_1..f_m of GF(2^m) (poly x^5+x^2+1 at m=5, x^4+x+1 at m=4, x^3+x+1 at m=3), so the ... |
| `zpool.py` | <zb dump> <k list>: pool the zero blocks E[Y,Y] = 0 derived by the star recursion (global facts about E), then (a) clique test: a set Y with E[Y,Y] = 0 and /Y/ > W is impossible (K[Y,Y] = I + ... |
