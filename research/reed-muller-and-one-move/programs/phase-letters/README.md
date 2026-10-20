# Phase letters: the symplectic game, the angle barrier, and phase letters at the ends

What these programs show: RM10 and the phase-letter half of RM5. `spg.py` builds the game on Sp(2m, 2) with exact
bookkeeping and `t_setup.py` checks it at blocks 2 to 5; `thresholds.py` computes, in exact rationals, what a
certificate must reach to refuse each family shape; `t_avg.py`, `t_lifts.py`, `t_lift4.py` test frame families
(lifts) and their refusal calibration (`t_lifts.py` reads `notes/network-coding/cert_E162_bitmodel.json` of the
repository root by default); `t_alt.py` checks the phase-extended targets at even blocks. In `ends-only-block4/`:
`t_small_alt.py` (the exhaustive small-lift tables with alternating targets), `fam_enum.py`, `tile_alt.py`,
`congr_classes.py` (congruence classes of symmetric matrices), with the ILP and rank-cover helpers copied from the
star-tool track. The skeptic's programs are in `independent-check/`. The proofs are in
`../../proofs/phase-letters-and-angle-barrier.md`.

Commands re-run for this revision (identical):

```
python3 thresholds.py > out_thresholds.txt    # under 1 s
python3 t_setup.py    > out_t_setup.txt       # about 4 minutes
```

Not re-run: the 4-in-15 enumeration with alternating targets (`ends-only-block4/`, remote machines).

Run everything from a working copy made by `sh ../workdir.sh` (from this folder), inside `work/phase-letters/`, where the data files of `../../data/phase-letters/` sit next to the programs. Python programs need Python 3.9 or later, numpy and scipy.

## Files

Each line is the file's own description (from its header), lightly shortened; the files themselves keep the working labels of the rounds, which the section README's Sources maps to these folders.

| File | What it does |
|---|---|
| `ends-only-block4/congr_classes.py` | check by orbit enumeration that the nonsingular symmetric m x m matrices over F_2 fall into congruence classes S -> A S A^T as follows: one class (that of I) for odd m; two classes (I and the ... |
| `ends-only-block4/fam_enum.py` | exhaustive S+D+C+Z pruning for ONE direction multiset with a large stabiliser, using the stabiliser to fix one target (t_small_alt.py enumerates every target multiset; this does the same job for ... |
| `ends-only-block4/glm.py` | Gauss-Newton / Levenberg-Marquardt for the multi-target slice, with J^H J built from per-form Gram matrices (no Jacobian is ever stored). F(Q) = A T (I + H T)^-1 B^T = sum over time-increasing ... |
| `ends-only-block4/ilp_alt.py` | the star-tool track's ILP piece + term-rank bound (rc_ilp.py, copied unchanged) on classes with targets of any kind; the verifier's second parity is replaced by the decomposition-based one of ... |
| `ends-only-block4/lm_alt.py` | the solver track's calibrated Levenberg-Marquardt (glm.py, copied unchanged) on per-array targets that may be ALTERNATING (Lemma A). The model only needs t_i = Q(z_i) and p_l = <Q, T_l> over all ... |
| `ends-only-block4/lm_batch_alt.py` | Levenberg-Marquardt on every class that the tiling test lets through, read from tile_alt.py's output lines "class i: dirs (...) targets (...) ...: passes ...". Per class: several random time or |
| `ends-only-block4/orbit_count.py` | sum of orbit sizes (under the direction multiset's stabiliser) of the survivor classes in a fam_enum dump, split by the number of alternating targets; cross-check against t_small_alt.py tuple counts. |
| `ends-only-block4/rc_ilp.py` | piece + term-rank bound on rank E[D_F, D_F], optimised exactly by an ILP (HiGHS). Bound: E[D_F,D_F] = sum over chosen disjoint pieces R x C (each inside D_G x D_G for a form G, so rank <= s_G) + ... |
| `ends-only-block4/rc_run.py` | build a multi-target system, run rcdp (guillotine rank-cover bound) on chosen forms, and verify every reported violation INDEPENDENTLY (partition of D_F x D_F, each piece's bound re-derived with a ... |
| `ends-only-block4/rcdp.c` | guillotine rank-cover bound for the multi-target system, fixed move counts. For a form F with support D_F (a set of directions; moves = directions with multiplicity counts[z]): rank E[D_F, D_F] = ... |
| `ends-only-block4/split_tight.py` | split a survivor dump into classes with and without a tight form (/D_Q/ = r_Q >= 1). |
| `ends-only-block4/t_free.py` | checks for the free-direction lemma. Lemma (free directions). In the multi-target slice, if the DISTINCT move directions y have F_2-independent rank-ones y y^T, then deletion alone forces M >= ... |
| `ends-only-block4/t_small_alt.py` | exhaustive exact pruning of SMALL multi-target slice programs (the small lift of Theorem B). For W arrays at block m with M moves: every M-multiset of directions up to G |
| `ends-only-block4/tile_alt.py` | the star-tool track's rank-cover tiling test (rcdp.c, Theorem 10 with a target per array, the general-m track 8.5) run on per-array targets that may be ALTERNATING (Lemma A of ../PROOF.md: phase ... |
| `ends-only-block4/xbase.py` | own field / forms / targets code for the multi-target system (Theorem B form). Coordinates: a self-dual basis f_1..f_m of GF(2^m) (poly x^5+x^2+1 at m=5, x^4+x+1 at m=4, x^3+x+1 at m=3), so the ... |
| `independent-check/v_a2.py` | helper (no description in the file) |
| `independent-check/v_arith.py` | helper (no description in the file) |
| `independent-check/v_cls.py` | the skeptic: own class-based D + C + Z enumeration at block 4, any W. usage: v_cls.py W M mode [shard nshards] mode 'alt': target W-multisets from all 448 with at least one alternating; 'noalt': ... |
| `independent-check/v_congr.py` | m = 5: count nonsingular symmetric and the orbit of I by random sampling m = 6: nonsingular alternating |
| `independent-check/v_congr5.py` | helper (no description in the file) |
| `independent-check/v_recheck.py` | Re-test the phase-letter track's survivor lists with MY D + C + Z (v_cls leaf logic) and count tight forms. Targets are read through the phase-letter track's own index -> matrix map ... |
| `independent-check/v_small.py` | the skeptic: own exhaustive D + C + Z enumeration at block 4 . usage: v_small.py W M first pool [chunk_index nchunks] first: 'alt' (T_1 = a fixed alternating target A0) or 'I' (T_1 = I) pool: ... |
| `independent-check/v_termrank.py` | Own exact test (valid, time-order free): rank E[D_F, D_F] <= term rank of the positions not forced to zero by tight supports (diagonal at live moves; pairs sharing a tight support). Exclusion if ... |
| `independent-check/vb.py` | the skeptic: own base code (written from the definitions in the general-m track 8.5-8.7, not from the phase-letter track's files). Block m. Directions z in F_2^m minus 0 (index = integer z-1). ... |
| `spg.py` | Lemma N's block as a game on Sp(2m,2), with exact F_2 bookkeeping. Vectors u in F_2^(2m) are ints: bits 0..m-1 = x, bits m..2m-1 = b (Pauli label X^x Z^b). omega(u,v) = x.b' + b.x'. A group ... |
| `t_alt.py` | the phase-extended per-array targets at even blocks, exact checks at m = 2, 4 (and 6 for (b)-(c)). With phases free, an array may be owed x(S) for ANY nonsingular symmetric S (two-sided free ... |
| `t_avg.py` | the Averaging Lemma written directly in Lemma N's block C[Sp(4,2)] (block 2), as a code check. Enumerate Sp(4,2) (720 elements) and P = Stab(L0) (48). Take a random family Y_U (U = the 15 frames) ... |
| `t_avg_lib.py` | exact integer rank (fraction-free Bareiss), shared by the phase-letter track scripts. |
| `t_free.py` | checks for the free-direction lemma. Lemma (free directions). In the multi-target slice, if the DISTINCT move directions y have F_2-independent rank-ones y y^T, then deletion alone forces M >= ... |
| `t_lift4.py` | an explicit frame family at block 4, built and checked here, as a test in C[Sp(8,2)]. A lift gives each of the 2295 frames U a 4-dim subspace Y_U of Q^8 spanned by integer rows reducing mod 2 to a ... |
| `t_lifts.py` | refusal calibration of the frame-family test in Lemma N's block, exact, from the rows only. block 2: search a lift (Y_U spanned by integer rows reducing mod 2 to a basis of U) of the 15 frames ... |
| `t_setup.py` | exact checks of Lemma N's Sp(2m,2) game at m = 2..5. (1) every generator is symplectic; (2) K = prod_j t_(e_j) = x(I) (the plain schedule, one array); (3) y(I) x(I) y(I) = w; (4) lev(A) and y(C) ... |
| `t_small_alt.py` | exhaustive exact pruning of SMALL multi-target slice programs (the small lift of Theorem B). For W arrays at block m with M moves: every M-multiset of directions up to G |
| `thresholds.py` | what a Theorem-R certificate (= a frame family, = a test in Lemma N's block) must reach to refuse the family's shapes W in 2^m - 1, with phases. Exact rationals. need(W) = (2^m - 1)/W (the ... |
