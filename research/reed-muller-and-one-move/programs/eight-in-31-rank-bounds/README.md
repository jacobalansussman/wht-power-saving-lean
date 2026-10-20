# 8 in 31 at block 5: rank bounds from below (round 3)

What these programs show: RM8. `starrec4.c` (stage 1: lower bounds that use the derived zeros), `starrec6.c` (stages
2 and 3: pooling, `NU=2` for forms with up to two odd arrays, `DUMPB=1` prints every bound), `starrec3.c` (lower
bounds without zeros), `starrec5.c` (block-plus-matching covers inside the recursion), `linkrec.c` (subsets of the
arrays). `bpm.c` with `vbpm.py` and `vrec_z.py` is the independent stage-1 check written for this round; `chkL.py`
is the hand check of the worked example in plain Python; `u3bounds.py` checks every derived bound on the true E of
genuine U_3 units (it uses `../shared/`). `cover.py`, `blockdef.py`, `defilp.py` are tests that did not bite.
The skeptic's own program is `independent-check/vr3.c` (with `mytab.py` for its tables, `u3check.py` and `u4mt.py`
for bounds on genuine units, `vmerge.py` and `vsummary.py` for the comparison). The verdicts are
`../../data/eight-in-31-rank-bounds/verdicts.txt` (one line per class), `undecided.txt` (the 127),
`hyperplane_classes.txt` (the 24), and the skeptic's `vverdicts.txt` and `extra_exclusions.txt` (the 31).

Command re-run for this revision (same verdicts as recorded):

```
cc -O2 -o starrec4 starrec4.c
python3 gen_tables.py 5                           # writes tables5.bin (the same file as the star-tool folder's)
printf '0,1,10,12,13,15,18,21\n0,1,10,12,14,15,17,21\n0,3,6\n' | MAXD=16 ROUNDS=3 LB=2 VMASK=1 ./starrec4 tables5.bin list
      # 15 s: the two survivors of round 2 are excluded ('L', lower bound above upper bound), the plain control passes
```

`chkL.py` takes the zero dump of a `starrec4` run as arguments (form, block and zero cells); the worked example's
numbers are in the section README (RM8). Not re-run: the stages over all classes and the skeptic's cascade (remote
machines; outputs in `job-outputs.tar.gz`). The skeptic's 49 MB dumps of bounds on genuine block-4 units are not
included; `independent-check/u4mt.py` regenerates them.

Run everything from a working copy made by `sh ../workdir.sh` (from this folder), inside `work/eight-in-31-rank-bounds/`, where the data files of `../../data/eight-in-31-rank-bounds/` sit next to the programs. Python programs need Python 3.9 or later, numpy, and scipy for two of the tests that did not bite.

## Files

Each line is the file's own description (from its header), lightly shortened; the files themselves keep the working labels of the rounds, which the section README's Sources maps to these folders.

| File | What it does |
|---|---|
| `blockdef.py` | block-diagonal deficiency test. If moves X split into parts P_1..P_t with every cross cell E[P_i,P_j] (i != j) a derived zero, then rank E[X,X] = sum rank E[P_i,P_i] <= sum u(P_i), while L[X,X] = ... |
| `bpm.c` | block-plus-matching test. Input (stdin, text): nf, then nf lines "mask weight slack" (mask: bit z-1 <-> direction z; slack = /D_Q/ - r_Q), then 31 hex masks Z[x] (bit y-1 set <-> |
| `chkL.py` | sanity check of one LB2 contradiction with plain Python (tables5.bin data, final zero dump). |
| `cover.py` | forced-full-rank test from matching equality . For a form Q with nu(D_Q) = s_Q (term rank of E[D_Q,D_Q] avoiding derived zeros equals the required rank), every minimum vertex cover (R', C') (rows ... |
| `defilp.py` | deficiency ILP, a strengthening of the clique test (CLQ). Fact: L = E + K is unit lower triangular in time order and rank K <= W. So for EVERY set X of moves /X/ = rank L[X,X] <= rank E[X,X] + ... |
| `fill.py` | fill the NOTE placeholders from verdicts.txt (after merge.py) and stage counts. |
| `gen_tables.py` | tables for scanrc.c (block 5, one move per direction, canonical Singer targets T_k = mult(alpha^2k)). Writes tables5.bin (little-endian int32): nf, then mask[nf] (bit z-1 <-> direction z), wt[nf], ... |
| `independent-check/mytab.py` | my own GF(2^m) field, quadratic forms and parities, written from the definitions only. Moves = nonzero field elements z (1..2^m-1, polynomial basis bits); D_Q = {z : Q(z) = 1}. |
| `independent-check/u3check.py` | my vr3 bounds (DUMP=1 POOL=1, block 3, targets 0,0) against the TRUE E of genuine numerical U_3 = 2 in 7 units . Direction z (refl2 coords) <-> my field element |
| `independent-check/u4check.py` | my vr3 bounds (DUMP=1) against the TRUE E of genuine numerical programs found by an earlier families round's solver refl2.System2 (targets all I = "0,...,0"), at block m with S = all 2^m - 1 ... |
| `independent-check/u4mt.py` | genuine numerical programs with a TARGET PER ARRAY (T_k = mult by alpha^(2k) in my self-dual basis f; frame alpha^k f_j) at block m, all 2^m-1 directions once, found by an earlier families round's ... |
| `independent-check/vmerge.py` | merge my vr3 cascade outputs (jobs/vr3a/remote-out/*/out_[ABC]_*.txt) with the claim's verdicts.txt. Writes ../vverdicts.txt: class, claimed verdict+stage, my verdict (stage A/B/C or PASS-all / ... |
| `independent-check/vr3.c` | a skeptic's from-scratch re-implementation of the round-3 class track's stages 1-3. Shares no code with starrec*.c / bpm.c / vrec*.py. Tables from mytab.py (own field, forms, parities). Facts used ... |
| `independent-check/vr3_v1.c` | vr3.c : a skeptic's from-scratch re-implementation of the round-3 class track's stages 1-3. Shares no code with starrec*.c / bpm.c / vrec*.py. Tables from mytab.py (own field, forms, parities). ... |
| `independent-check/vsummary.py` | table of claimed stage x my cascade stage from ../vverdicts.txt; writes ../extra_exclusions.txt |
| `linkrec.c` | starrec2.c : star-lemma recursion + pooled zero blocks + clique and matching tests, for the multi-target system at block m (3, 4, 5) with ONE move p |
| `merge.py` | one verdict per survivor from the stage outputs (first verdict found wins in stage order). usage: merge.py surv8_final.txt stage1_files... -- stage2_files... -- ... prints per stage counts, writes ... |
| `starrec2.c` | star-lemma recursion + pooled zero blocks + clique and matching tests, for the multi-target system at block m (3, 4, 5) with ONE move per nonzero direction and canonical Singer targets T_k = ... |
| `starrec3.c` | starrec2 + LOWER-BOUND propagation (env LB=1, default on): ell(A) <= rank E[S_A,S_A] from (i) the top, rank E[D_F,D_F] = /D_F/ - nu, and deletion: rank E[S_A] >= rank E[S_B] - 2/S_B \ S_A/; (ii) ... |
| `starrec4.c` | starrec3 with LB=2: the top-down lower bounds use the zero cells: rank E[X,X] <= rank E[S,S] + (largest matching of X x X minus S x S avoiding derived zeros), for S = S_A inside X = D_F (4- and ... |
| `starrec4d.c` | starrec4.c : starrec3 with LB=2: the top-down lower bounds use the zero cells: rank E[X,X] <= rank E[S,S] + (largest matching of X x X minus S x S avoiding derived zeros), for S = S_A inside X = ... |
| `starrec5.c` | starrec4 + BLOCK-PLUS-MATCHING covers after each round (env BPM=1/2): F with r_F <= numax (rank E[D_F] = /D_F/ - nu), blocks S_i = D_F n D_G_i (rank E[S_i] <= s_G_i), disjoint for two blocks: ... |
| `starrec6.c` | starrec4 + MULTI-FORM POOLING of rank bounds (env POOL=1): every upper bound u(A) a form's recursion derives below its start value, and every lower bound ell(A) > 0, is a fact about the one matrix ... |
| `u3bounds.py` | soundness check of starrec6's bounds on GENUINE programs. Numerical U_3 2-in-7 units (block 3, targets I, I, all 7 directions once; an earlier families round solver) give a true E. Every bound ... |
| `u3trace.py` | in genuine numerical U_3 2-in-7 units , print diag(E) and per-array traces a_l.b_l. |
| `vbpm.py` | INDEPENDENT check of the new exclusions. Zero cells of E are re-derived by the skeptic's own code (vrec_z.py = the skeptic vrec.py unchanged except that it keeps its zero set; own field, self-dual ... |
| `vf.py` | skeptic's own GF(32) / forms / canonical-target parities, written from scratch. Field: GF(2)[x]/(x^5+x^2+1) by carry-less multiplication (no log tables). Directions z in 1..31, bit j = coordinate ... |
| `vmw.py` | tables ks t: hex mask of parity classes v with w(v) = #{l : <class v, T_k_l> = 1} <= t. |
| `vrec_z.py` | vrec.py : skeptic's own re-implementation of the star recursion + pooled zero blocks + clique + matching + /Y/-W tests for the per-array-target |
| `xbase.py` | own field / forms / targets code for the multi-target system (Theorem B form). Coordinates: a self-dual basis f_1..f_m of GF(2^m) (poly x^5+x^2+1 at m=5, x^4+x+1 at m=4, x^3+x+1 at m=3), so the ... |
| `zstruct.py` | direction z (coordinate int) <-> field element: coord(x) = sum_j Tr(x f_j) 2^j classify by the field ratio / sum of the two directions |
| `zsum.py` | summarise linkrec output: per subset, number of zero cells and zero-diagonal moves. |
