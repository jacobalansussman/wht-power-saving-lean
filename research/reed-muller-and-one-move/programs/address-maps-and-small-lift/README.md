# Address maps as per-array targets, block 3, and the small lift

What these programs show (section 3 of the section README): RM1 (address maps without phase letters are per-array
targets; `t_lift.py` is an exact group-algebra check of Theorem B at blocks 3 and 4), RM3 (Theorem C3: `t_face3.py`
recomputes the block-3 face certificate in the form the proof needs; `independent-check/vk_face3.py` is the skeptic's
exact recomputation over Q(i)), RM5 (the exhaustive small-lift tables: `t_small.py`, with `t_free.py`, `t_orbits.py`,
the cross-checks `t_crosscheck*.py` and the term-rank re-test `tr_check.py`), RM11 (the field model and the placement
theorem: `gf.py`, `t_field.py`, `t_place.py`), the 2.706 transfer at block 4 (`t_cert4.py`), and block-5 scans with
the zero-block test (`mt_scan5z.py`, `t_tightcode.py`). The `mt_*.py` and `lm_batch.py` files are numerical solvers
(evidence only).

Commands (each was re-run for this revision on an Apple M4 Mac mini; outputs identical apart from timings):

```
python3 t_face3.py                 > out_t_face3.txt                 # 6 s, 25 MB: caps 8 x 7, demand 24 in all 28 bases
python3 t_small.py 3 3 8           > out_t_small_b3_W3_M8.txt        # under 1 s: no 3 in 8 at block 3
python3 t_small.py 4 3 10          > out_t_small_b4_W3_M10_nofilter.txt   # 15 s: no 3 in 10 at block 4
(cd independent-check && python3 vk_face3.py > out_vk_face3.txt)     # 14 s, exact over Q(i)
```

Heavy, not re-run for this revision: `python3 t_small.py 4 3 11` and `python3 t_small.py 4 4 13` (hours on an
8-processor Linux machine; outputs `out_W3_M11_nofilter.txt`, `out_W4_M13.txt`); `t_small.py 4 4 14` was stopped
unfinished (`out_W4_M14_partial.txt`, at least 10 surviving classes, below the bar); `t_cert4.py` (about an hour,
output `out_cert4.txt`). Note: a docstring says 420 bases at block 4; there are 840 unordered bases of the 4-bit address
space (20,160 / 24), as the program itself counts, which a skeptic pointed out.

Run everything from a working copy made by `sh ../workdir.sh` (from this folder), inside `work/address-maps-and-small-lift/`, where the data files of `../../data/address-maps-and-small-lift/` sit next to the programs. Python programs need Python 3.9 or later and numpy.

## Files

Each line is the file's own description (from its header), lightly shortened; the files themselves keep the working labels of the rounds, which the section README's Sources maps to these folders.

| File | What it does |
|---|---|
| `gf.py` | small GF(2^m) toolkit for the general-m track (field model of the Singer no-phase game). Elements are ints 0..q-1 in the polynomial basis of a primitive polynomial. Everything is exact (bit ... |
| `independent-check/vk_b5.py` | the skeptic: own check of the general-m track 5b 'block 5' (mt_scan5z / t_tightcode claims) on MY OWN random draws. Block 5, one move per direction (31 moves), W distinct canonical Singer targets ... |
| `independent-check/vk_cert4.py` | the skeptic: own recomputation of the general-m track 5c (block-4 transfer, row 3 of wedge4's cert_4col.txt). N = span of b(P) = sum over the 4 cosets t+P of e_{t+p1} ^ e_{t+p2} ^ e_{t+p3} ^ ... |
| `independent-check/vk_cert4_float.py` | float (complex SVD) check that row-r caps of wedge4's cert equal the mod-p values (mod-p rank <= true rank, so equality in float means the file's exact caps are reproduced as true ranks). usage: ... |
| `independent-check/vk_face3.py` | the skeptic: independent recomputation of the block-3 face certificate (Theorem C3 inputs). Own code. Exact ranks over Q(i) via realification (integer matrices, exact Fraction elimination), plus ... |
| `independent-check/vk_field.py` | the skeptic: own check of the general-m track's canonical-placement claims (section 2, 3, 8.2-8.4, 8.8) at m = 3..6. Own GF(2^m) code (primitive polys chosen here), own self-dual basis search, own ... |
| `independent-check/vk_field_m3.py` | the skeptic: own check of the general-m track's canonical-placement claims (section 2, 3, 8.2-8.4, 8.8) at m = 3..6. Own GF(2^m) code (primitive polys chosen here), own self-dual basis search, own ... |
| `independent-check/vk_small.py` | the skeptic: independent small-lift enumerator , own code and own enumeration order. Instead of direction classes up to GL with all target tuples, it FIXES T_1 = I (GL is transitive on targets ... |
| `independent-check/vk_small2.py` | the skeptic: independent small-lift enumerator , own code and own enumeration order. Instead of direction classes up to GL with all target tuples, it FIXES T_1 = I (GL is transitive on targets ... |
| `independent-check/vk_survtest.py` | the skeptic: independent small-lift enumerator , own code and own enumeration order. Instead of direction classes up to GL with all target tuples, it FIXES T_1 = I (GL is transitive on targets ... |
| `lm_batch.py` | Levenberg-Marquardt (mt_lm.MTSystem) on a file of small-lift survivors 'dirs/targets' (t_small.py DUMP format; targets are 'gl:' indices). Each line: orders lex and rand<k>, `restarts` restarts ... |
| `mt_b3.py` | block 3 probe of target diversity. Multi-target slice = no-phase Singer game (Theorem B). Shapes: 5 in 14 (2.8 = the hook bound, not excluded by any test in the canonical placement), controls 5 in ... |
| `mt_b3neg.py` | find block-3 multi-target instances that are impossible (below the engine hook 2.8 moves per array: 3 in 8, 4 in 11) yet pass the exact tests: these are the NEGATIVE controls for the LM probe. |
| `mt_b3sweep.py` | block 3, the whole no-phase game (Theorem B: targets = any of the 28 matrices congruent to I). Random 5-in-14 instances (random targets, random direction multisets using every direction, random ... |
| `mt_cover.py` | the slice's covering test for the MULTI-TARGET slice with one move per direction. For a fully even form F (no array odd) the block E[D_F, D_F] is nonsingular; if every set of <= s+ |
| `mt_floor.py` | where does the LM residual sit at its floor? Per form: residual norm^2, the form's weight on the moves, n_odd, slack, and whether the form is linear / quadratic (b and c parts in the canonical ... |
| `mt_lm.py` | the MULTI-TARGET commutative slice, which by Theorem B is exactly the no-phase Singer game. A multi-target slice program at block m: W arrays, array l owed [T_l] (T_l symmetric over F_2), moves [z ... |
| `mt_scan5.py` | block 5, the free-target U_5 (8, 9, 10 arrays in 31 moves, one per direction), canonical Singer targets T_k = sigma^(k) sigma^(k)T (k in Z_31). Exact tests only (deletion, tight forms, tight ... |
| `mt_scan5z.py` | mt_scan5.py's random canonical-Singer target sets for the uniform family at block 5 (W in 31, one move per direction), now with the zero-block test of NOTE 8.7 added (t_small.leaf_tests). |
| `t_cert4.py` | the block-4 certificate of an earlier round (wedge4/cert_4col.txt: four kept subspaces of Lambda^4(C^16), multipliers 438/6095, 392/6095, 45/1219, 22/1219, bound 3728/1219 for the STANDARD ... |
| `t_crosscheck.py` | compare t_small.py's tests (span, deletion, cliques) with mt_lm.py's independent exact_tests on random small-lift configurations (block m, W arrays, M moves; targets drawn from those in the span ... |
| `t_crosscheck2.py` | for one direction multiset, every target W-multiset that passes deletion (t_small's code) is re-tested with mt_lm.Setup.exact_tests (independent deletion/tight/clique code); verdicts compared. |
| `t_face3.py` | an independent recomputation of the block-3 FACE certificate (an earlier round), in the form needed for Theorem C3 (NOTE 8.9): f = 3 digit columns, restricted to the antisymmetric part ... |
| `t_field.py` | exact checks of the field model of the Singer no-phase game . Coordinates: a self-dual basis f_1..f_m of F_q (q = 2^m) over F_2, so the dot product is Tr(xy) and the target I is the trace form. ... |
| `t_free.py` | checks for the free-direction lemma. Lemma (free directions). In the multi-target slice, if the DISTINCT move directions y have F_2-independent rank-ones y y^T, then deletion alone forces M >= ... |
| `t_lift.py` | exact check of the "boundary lift" : a multi-target commutative-slice program whose array l has target T_l = h_l h_l^T (h_l in H) becomes a unit of C[Sym_m x/ H] with the same moves, by ... |
| `t_nested.py` | the NESTED-SUPPORT test for programs given by a direction list and per-array targets: for forms F, G with D_F a subset of D_G (as multisets of moves), /D_F/ + r_F <= /D_G/ + r_G. Proof: rank ... |
| `t_orbits.py` | completeness check of t_small.py's enumeration: the GL(m,2)-orbits of M-multisets of nonzero directions must have sizes summing to C(2^m - 1 + M - 1, M). Same orbit code as t_small.py. |
| `t_place.py` | every placement of the Singer group relative to the target, at block m . Fix the Singer cycle sigma: x -> alpha x in self-dual coordinates. A placement is a target T (symmetric, invertible, ... |
| `t_quot.py` | single-array cost in the Gold quotient games of the trace placement. Quotient of C[Sym_m x/ Z_n] by the Singer-invariant part of Alt_m outside one cyclotomic class 1 + 2^i: the group F_q (+) F_q ... |
| `t_small.py` | exhaustive exact pruning of SMALL multi-target slice programs (the small lift of Theorem B). For W arrays at block m with M moves: every M-multiset of directions up to GL(m,2); every W-multiset of ... |
| `t_tightcode.py` | block 5, canonical Singer targets, one move per direction. Check that 'some tight form' <=> 'the trace code {(Tr(b gamma_l))_l : b} on the targets has a word of weight 2^(m-2)', and that every ... |
| `tr_check.py` | re-test a survivors file 'dirs/targets' with the leaf tests of t_small.py (TERMRANK=1 adds the term-rank test). Prints how many survive and why the others fail. |
