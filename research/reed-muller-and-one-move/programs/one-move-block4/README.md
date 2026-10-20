# One move at block 4: pairs, triples and more arrays (tier 2)

What these programs show: the pair half of RM5, RM13, and section 5. `certify_t.py` certifies the multi-target
linear-program bound for each kind of target pair (`certpairs.sh` runs all 20; output `out_certpairs4.txt`, every pair
at least 56/15 moves per array); `deltri.py` and `tclasses.py` count the kinds of target triple and their deletion
minima; `lp_triples.py` and `mtlp_t.py` are the linear programs for triples; `scan.py`, `scan_s.py` scan shapes;
`lms.py`, `lmh.py`, `long.py` are the numerical solvers (evidence only). The skeptic's programs are in
`independent-check/`: `vpairs.py` (the 420 targets, the 20 kinds of pair, an ILP deletion minimum for every target,
and the exhaustive enumeration of all 116,280 multisets of 7 directions), `vlp.py`, `vlp3.py`, `vglclass.py` (the 758
kinds of triple), `vtri.py`, `vlm.py`.

Command re-run for this revision (same counts):

```
cd independent-check
python3 vpairs.py     # 20 s: 20 pair classes; deletion minimum 8 for all 420 targets; 0 of the 116,280 multisets pass
```

Not re-run: `certpairs.sh` (3 minutes; it also writes three 5 MB row dumps that are not included) and the solver runs.

Run everything from a working copy made by `sh ../workdir.sh` (from this folder), inside `work/one-move-block4/`, where the data files of `../../data/one-move-block4/` sit next to the programs. Python programs need Python 3.9 or later, numpy and scipy.

## Files

Each line is the file's own description (from its header), lightly shortened; the files themselves keep the working labels of the rounds, which the section README's Sources maps to these folders.

| File | What it does |
|---|---|
| `certify_t.py` | certify.py : exact, independent check of a cell-partition multi-target LP bound (mtlp.py --dump). Rebuilds |
| `deltri.py` | the exact DELETION minimum for multi-target slice programs at block 4. For a target multiset {T_1..T_W}, deletion says: for every form Q, (moves with Q(z) = 1) >= (arrays with <Q,T> = 1). Minimise ... |
| `forms.py` | block-m binary quadratic forms, GL(m,2) action, the virtual-target distribution of Lemma C (targets h h^T, h in H), and per-orbit support size s(F) and odd fraction a(F). Conventions: direction x ... |
| `independent-check/lmsumm.py` | summary of my LM trial files (count, lowest costs, histogram of floors). |
| `independent-check/vcore.py` | my own block-4 machinery, written without reading the tier-2 block-4 track's helpers' internals. Conventions (mine): a direction is a 0/1 numpy vector of length m (index 1..15 <-> integer with bit ... |
| `independent-check/vglclass.py` | merge the O_4-classes of {I, T2, T3} into classes of target multisets under the full relabelling {T_l} -> {g T_l g^T} (g in GL(4,2)), which leaves the multi-target game (and deletion) unchanged. ... |
| `independent-check/vlm.py` | my own solver for the multi-target commutative slice at block 4 (exact formulation). Derivation (mine): a program's registers are, as formal elements, the inputs plus scalar combinations of the ... |
| `independent-check/vlp.py` | my own build of the multi-target rank-cover LP , in my own conventions (vcore.py), plus my own exact dual certificate. Per array normalisation: c_x = (#moves in direction x)/N, p_xy (x <= y) = m |
| `independent-check/vlp3.py` | my own build of the multi-target rank-cover LP , in my own conventions (vcore.py), plus my own exact dual certificate. Per array normalisation: c_x = (#moves in directio |
| `independent-check/vpairs.py` | my own check of the pair claims at block 4. (1) number of targets h h^T, /O_4/, number of O_4-classes of T (pairs {I, T}); (2) for EVERY one of the 420 targets T: integer deletion minimum of {I, ... |
| `independent-check/vpairs_sanity.py` | (a) the float32 matmul used in vpairs.py equals the exact integer count on a sample of multisets; (b) positive control for the deletion test: plain 2 in 8 (directions e1..e4 and the columns of h, ... |
| `independent-check/vscan.py` | independent LM sampling (my solver vlm.py) for 3 in 11 at block 4, multi-target slice. Mode 'perm': random orders of a fixed direction multiset; mode 'rand': random multisets of 11 directions that ... |
| `independent-check/vtri.py` | my own classification of target triples {I, T2, T3} (T2, T3 any targets h h^T, repetition allowed) up to O_4 acting by T -> g T g^T, and the exact integer DELETION minimum of each class (ILP, ... |
| `lmh.py` | complex Levenberg-Marquardt with an analytic Jacobian for the NO-PHASE engine game with a finite address group H, any H (not only cyclic). Exact formulation by an earlier construction track ... |
| `lms.py` | the MULTI-TARGET commutative slice at block m, exact formulation, complex LM (analytic Jacobian). Algebra C[Sym_m(F_2)] : moves multiply a register by [z z^T]; free = scalar mixing; array l must ... |
| `long.py` | one long LM run of the multi-target slice (lms.py) on a FIXED order of directions; prints cost and /x/ every 50 iterations (border behaviour: cost -> 0 with /x/ -> infinity would be an approximate ... |
| `lp_triples.py` | the multi-target rank-cover LP (mtlp_t.py's rows, built once) for random target sets {I, T2, ..., TW} at block 4; only the right-hand side a_F changes. Usage: lp_triples.py W count seed |
| `mtlp_t.py` | \nmtlp.py : the multi-target rank-cover LP. By Lemma C of an earlier construction track, a unit of C[Sym_m x/ H] (moves z_t) gives, in the channe |
| `scan.py` | many LM trials of lmh.py's exact system for one address group H and one shape W in M. Each trial: a random sequence of direction ORBIT representatives (directions matter only up to H, an earlier ... |
| `scan_s.py` | LM trials of the multi-target slice (lms.py) at block m for W arrays in M moves. Targets: --targets fixed list, or --randT: T_1 = I and T_2..T_W random among the h h^T (all 420 at m = 4). ... |
| `summ.py` | helper (no description in the file) |
| `tclasses.py` | classes of target multisets {I, T2, ..., TW} for the multi-target slice at block m, up to the orthogonal group (g g^T = I), which fixes I and acts by T -> g T g^T. T ranges over h h^T (h in ... |
