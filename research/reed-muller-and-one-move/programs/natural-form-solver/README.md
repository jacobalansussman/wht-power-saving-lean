# The natural form at block 5: the calibrated solver, and sharing next to plain

What these programs show: section 5 of the section README (the search found nothing) and section 4, item 2 (the
corrections to this track). `inst.py` builds instances of the multi-target system; `glm.py` is the Gauss-Newton /
Levenberg-Marquardt solver judged on every form; `run.py`, `warm.py`, `warm2.py`, `homotopy.py`, `homchain.py`,
`tri_direct.py` drive it; `packing.py` and `packing2.py` count frames and plain packings (every canonical target at
block 5 has exactly 6 frames); `scan_frames.py`, `scan_plain.py`, `scan_sets.py` run the exact tests on the shapes
next to plain. The skeptic's programs are in `independent-check/` (`vtests.py`, `vcert.py`, `vtrcert.py`,
`vscan3.py`, `vw4b5.py` and others); its scan found the 9 four-array sets at block 5 that pass every test
(`survivors_b5w4.txt` in the data folder).

Re-run for this revision: `python3 packing2.py` (about a minute; it prints `frames per target: {0: 6, 1: 6, 2: 6}`). The solver runs (remote
8-processor machines, hours) were not re-run.

Run everything from a working copy made by `sh ../workdir.sh` (from this folder), inside `work/natural-form-solver/`, where the data files of `../../data/natural-form-solver/` sit next to the programs. Python programs need Python 3.9 or later and numpy.

## Files

Each line is the file's own description (from its header), lightly shortened; the files themselves keep the working labels of the rounds, which the section README's Sources maps to these folders.

| File | What it does |
|---|---|
| `gf.py` | small GF(2^m) toolkit for the general-m track (field model of the Singer no-phase game). Elements are ints 0..q-1 in the polynomial basis of a primitive polynomial. Everything is exact (bit ... |
| `glm.py` | Gauss-Newton / Levenberg-Marquardt for the multi-target slice, with J^H J built from per-form Gram matrices (no Jacobian is ever stored). F(Q) = A T (I + H T)^-1 B^T = sum over time-increasing ... |
| `homchain.py` | deletion CHAIN from a solved plain program, one move at a time. usage: homchain.py SPEC ORDER NSUB MAXDEL CANDS SEED TLIMIT [X0.npy] [EXTRA] SPEC 'k:0,2,4,8,...' canonical Singer targets T_k ... |
| `homotopy.py` | deletion homotopy from a solved plain W-in-5W program to W in M_target. usage: homotopy.py KS ORDERSEED NDEL DELSEED nsub mu0 mu1 steps per [x40.npy] 1. plain instance (frames of the targets T_k, ... |
| `independent-check/vcases4.py` | helper (no description in the file) |
| `independent-check/vcert.py` | write a line-by-line zero-block certificate for named instances , then CHECK it from scratch in pure Python (no numpy): every claimed tight form is recomputed, the clique pairs are checked, and ... |
| `independent-check/vchain_tests.py` | exact tests (mine) on the 8-in-39 deletion candidates of the the solver track chains, in the solver track's coordinates. |
| `independent-check/veval.py` | independent residual of a saved point. Data (targets T_l, directions z_i in time order) are taken from the solver track's inst.canonical; the forms, parities and F = A_D (I + H_D)^-1 B_D^T (direct ... |
| `independent-check/vfield.py` | independent GF(2^m) model. Targets as symmetric F2 matrices (bitmask of upper triangle), directions as ints 1..2^m-1 in a self-dual basis found by exhaustive search (not greedy-random). |
| `independent-check/vframes.py` | check symmetric and self-dual also: any multiset decomposition with fewer than m distinct directions? rank argument says no; check sets of size < m |
| `independent-check/vpack.py` | pairwise: do any two frames of I intersect? additive set: check that alpha^(2k) for k in (1,9,16,18,21,25) and 0 ... form subspace |
| `independent-check/vplainctl.py` | control: every plain multiset (W frames, no deletion) must pass S, D, C, Z, TR (they are realised programs). |
| `independent-check/vscan3.py` | helper (no description in the file) |
| `independent-check/vscan3f.py` | block-m W-array 'plain minus one move' scan; S+D screened by the fast form, then every S+D pass re-checked directly and run through C, Z, TR (vtests.py). |
| `independent-check/vtests.py` | independent exact tests for a multi-target commutative program shape (multiset of directions, list of targets): S (span), D (deletion), C (tight-pair clique <= W), Z (zero-block), TR (term rank of ... |
| `independent-check/vtight31.py` | helper (no description in the file) |
| `independent-check/vtrcert.py` | term-rank certificate for the the general-m track set's 8 in 39 , in the solver track's coordinates. Finds the form F with term rank < /D_F/ - r_F, and a Konig cover (rows Rc, columns Cc covering eve |
| `independent-check/vw4b5.py` | exact tests (mine, incl. Z and TR) on a random sample of the solver track's block-5 W=4 S+D passes (their coordinates). |
| `inst.py` | instances of the multi-target commutative slice (Theorem B form of the no-phase Singer game). A program at block m: W arrays, array l owed [T_l] (T_l symmetric over F_2), moves [z_i z_i^T] in time ... |
| `packing.py` | for a canonical target set at block 5, the largest number of arrays whose frames are pairwise disjoint (the plain packing number in the one-move-per-direction game), and the disjointness graph. |
| `packing2.py` | ALL frames of each canonical target (5-sets of directions z with sum z z^T = T_k over F2) and the true plain packing number of a target multiset in the one-move-per-direction game (pairwise ... |
| `parse_runs.py` | LOG...: for each RESTART line of run.py logs, which arrays are served exactly (per-array < 1e-10), whether their frames are pairwise disjoint (a plain packing), and the residual. |
| `probe_border.py` | m TSPEC ORDER seed iters: long LM run; print cost and parameter norm along the way. |
| `run.py` | LM runs on one instance, several random restarts. usage: run.py m TSPEC DIRS ORDER restarts iters nsub seed tlimit [sa sh] [tag] TSPEC: 'k:0,2,4,...' canonical Singer targets T_k; 'I:8' slice (all ... |
| `scan_frames.py` | W [maxclasses]: every array uses ANY of its target's 6 frames (plain programs W in 5W). For each target class (up to k -> 2^s k + c) and each frame choice: which single deletions pass the exact ... |
| `scan_plain.py` | W: canonical target sets {0 = k_1 <= k_2 <= ...} of size W at block 5, up to the affine maps k -> 2^s k + c; for each, the plain frames' exact tests and how many single deletions pass span + deletion. |
| `scan_sets.py` | exact tests (span, deletion, tight forms) for structured target sets at block 5, one move per direction (U_5 with per-array canonical Singer targets T_k). |
| `served.py` | SPEC ORDER XFILE: restrict a solved run to its exactly served arrays; for each move, the residual of that sub-program when the move is dropped (no re-fit); the essential moves of each served array. |
| `t_check.py` | check glm.py's F, gradient and J^H J against direct formulas (block 3, small sizes). |
| `t_check2.py` | helper (no description in the file) |
| `table.py` | summary of all 8/9/10-in-31 random-start runs (run.py logs): per target set, runs, best residual on all forms, most arrays served, plain packing number. |
| `tri_direct.py` | CASES_FILE NRESTART ITERS NSUB: random-start LM on the multisets 'plain frames minus one move' that pass the exact tests (lines 'PASS ks=(..) frames=(..) delete dir z ...' of scan_frames.py). Each ... |
| `warm.py` | KS ORDER SEED ITERS NSUB TLIMIT [eps]: 8 in 31 (one move per direction) started from the LARGEST plain packing: the arrays of a maximum family of pairwise disjoint frames (packing2.py) get the ... |
| `warm2.py` | KS ORDER SEED ITERS NSUB TLIMIT [eps]: 8 in 31 (one move per direction) started with EVERY array on an XOR chain over one of its frames, frames chosen to minimise the number of moves used by two ... |
