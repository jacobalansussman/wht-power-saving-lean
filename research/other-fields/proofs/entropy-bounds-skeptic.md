> **What this is.** The record of the independent AI skeptic that checked the entropy-bounds round on 15 October 2026, with its own code (private paths removed). "I" in the file is the skeptic. Its programs are in [../programs/entropy-bounds/skeptic/](../programs/entropy-bounds/skeptic/), the outputs of its Colab runs in [../data/entropy-bounds/skeptic-colab-runs/](../data/entropy-bounds/skeptic-colab-runs/) (file names start with the run name, `job_a` to `job_w`).

# v-nc-gain: skeptic's check of nc-gain (entropy and cut bounds for coherent units), session 2

15 October 2026, 17:56 to 19:12 EDT. I tried to refute nc-gain's claims with my own code: my own frame
graphs, my own cut LPs and Shannon LP, and direct checks of the proofs. Everything is in `calc/` and the `job_*`
folders here. Predictions are in STATUS.txt (written before any computation), the log in RUNLOG.txt.

## Verdict

**Upheld, with wording corrections (not refuted).** Every number I recomputed agrees with nc-gain, both new
lemmas check out (by reading and by direct exact checks on all three graphs), and the LP and inequality encodings
pass known-answer tests. The corrections are about scope sentences, listed at the end.

In one paragraph: my own Γ_2, Γ_3, Γ_4 coincide with nc-gain's frame for frame. Lemma W and Lemma W+ are correct,
and I checked their certificates directly, exactly, on random and hill-climbed layouts at blocks 2 to 4. The near-polygon
claim holds by brute force. The exact check in `w3screen.exact_value` is sound, and the rigorous PDLP bound is valid.
My own LPs reproduce every quoted value I tried: block 2 (9), block 3 (9, 9, 9, 9, 11, 13, 14, image 8.5), block 4
images (17 layouts), and a full Γ_4 solve (12). An undirected solve of the same layout also gives 12, so that result
does not depend on orientation. My own Shannon LP finds the butterfly's gain (7 against 8) and a second known gain
(combination network, 9 against 10), and it reproduces session 1's three-cell contraction of Γ_2 (80). nc-gain's
Ingleton and DFZ rows pass known-answer tests. Going further than nc-gain:
- **Unsent hard layouts.** I solved all of the 16 hard block-4 layouts that nc-gain did not send to the full LP.
  Every one has cut value 12 (two by the directed LP, fourteen by the undirected one, which covers every orientation; every rigorous bound is at least 11.999997).
- **Independent samples.** At block 3 I drew my own sample (15 random non-compatible layouts and 5 hill-climbed ones
  with e = 2): all 20 have cut value exactly 9 (my LP, rigorous bound 9.000000). At block 4 I ran my own image screen on 80 random non-compatible layouts, drawn from 326 (75.5 % were compatible). Small images certified 78 of them (11.25 to 12). The other two stopped at exactly 11.0, and both have cut value 12 on all of Γ_4: [(1477,944),(1123,424),(1700,2056)] (e = 2) and [(2220,467),(859,547),(468,334)] (e = 1, undirected, bound 11.999999).
- **A cut gap with four arrays at block 4.** nc-gain had a W = 4 cut gap only at block 3. I found two block-4
  triangle-plus-pendant layouts with cut values 14 and 15, against 16.

## What I checked, and how

1. **My own Γ_m** (`calc/vgamma.py`). I built the Lagrangian subspaces of F_2^(2m) by extending isotropic
   subspaces, with adjacency "meet in dimension m − 1" and distances by breadth-first search. Γ_2, Γ_3, Γ_4 have 15,
   135, 2,295 frames, degrees 6, 14, 30, diameter m, and 120, 8,640, 2,350,080 ordered opposite pairs; the distance
   is m − dim(E ∩ F) everywhere. Mapped onto nc-gain's numbering (by equal subspaces), the edge sets are identical
   at m = 2, 3, 4. So every frame number in nc-gain's files means the same frame in mine.

2. **The symmetry behind orbit variables and averaging** (`calc/vsym.py`). The transvections act as graph
   automorphisms; the group they generate is transitive on ordered opposite pairs (120 of 120; 8,640 of 8,640) and
   on arcs (90 of 90; 1,890 of 1,890). This is what Lemma 1 (the uniform layout is the worst) and Lemma 2 (averaging
   over the group) need. Enumerating all 2^15 cuts of Γ_2 gives the uniform cut value 3/2 per session, as quoted.

3. **Lemma W** (two-class stretch). The proof is correct: max(|x|, |y|) = (|x + y| + |x − y|)/2, so the max of two
   1-Lipschitz line metrics is a sum of two line metrics, each a cut packing with every edge load at most 1, and
   each pair is stretched to m by its class's function. The converse (one function exists only for a stretchable
   set) is also right: for two pairs i, j the two cross inequalities add to d(low_i, high_j) + d(low_j, high_i) ≥ 2m.
   Direct check (`calc/vlemmas.py`): for every two-class stretchable layout among 300 random ones for each
   W = 3, 4, 5, 6 at each block 2, 3, 4, I built the packing from φ_A, φ_B and checked in exact arithmetic that every
   edge load is ≤ 1 and every pair gets exactly m. No failure. At block 2 I also solved the cut LP with every one of
   the 2^14 cuts listed (exact rational dual) on 160 layouts with W = 3..6 (shared frames allowed): Lemma W never
   failed, and both cut gaps I found (W = 6) are on layouts it does not cover.

4. **Lemma W+** and its near-polygon claim (b).
   - Claim (b) by brute force: every triple (u opposite y, d(w, y) < m) has a neighbour y' of y opposite u with
     d(w, y') = d(w, y) + 1. All 840 triples at block 2, all 613,440 at block 3, and 381,300 at block 4 (300 random
     opposite pairs (u, y), every w; by transitivity on opposite pairs this is every case up to symmetry). No failure.
   - I read the proof of (b) and of (a), (c): correct.
   - Direct check of the construction: for each random layout with e(D) ≥ 1 (482, 691, 824 at blocks 2, 3, 4) and 40
     hill-climbed ones per block (up to e = 3 at block 4), I moved the ends by claim (b), confirmed the moves took
     exactly a and b steps and gave a compatible layout, and evaluated Lemma W's packing for the new layout on the
     old pairs: always ≥ 3m − e(D), with every edge load ≤ 1. No failure.
   - The e(D) distribution in my own sample (3,000 per block): block 2: 0 84 %, 1 16 %; block 3: 0 77 %, 1 22.8 %,
     2 0.23 %; block 4: 0 72.5 %, 1 26.8 %, 2 0.63 %. nc-gain: 85/15; 78/21/0.3; 74.5/24.9/0.55. Consistent.
   - Two-class stretchable shares: block 4, W = 4: 35 % (nc-gain 38 %); block 3, W = 5: 12.4 % of 2,000 (12 %).

5. **The pull-back screen and its exact check.** I read `w3screen.exact_value` and `w3image.cut_value`. The dual
   multipliers are rationalised, the edge loads are summed in exact fractions over the generated cut rows, the packing
   is scaled down if any load exceeds 1, and the value is Σ λ·sep. That is a valid cut packing of the image graph;
   pulled back along the distance map (adjacent frames map to equal or adjacent images; pairs keep their sep counts)
   it is a valid packing of Γ_m. When constraint generation stops at its round cap, the last dual still covers exactly
   the rows of the last LP, so the check stays valid. With M an integer, a value above 3m − 1 gives M ≥ 3m. Sound.
   - My own image graphs and my own LP (compact flow form, HiGHS, with my own Lagrangian lower bound) on 14 of
     nc-gain's certified block-4 layouts (the 6 lowest and 8 random, with nc-gain's orientation) give the same
     values (11.125 to 12); on one, [(1451,1722),(338,243),(579,532)], my LP gives 11.3333 where nc-gain's exact
     value is 11.3043 (theirs is a slightly weaker but valid packing). On three uncertified hard layouts I get the
     same best image values: 10.5, 11.0, 10.667.
   - Count check: the seven JSON files hold 701 layouts (524 compatible, 177 certified, 0 open); replaying seed 18
     with nc-gain's sampler on my graph gives seven compatible layouts and then [(2051,594),(2031,712),(777,952)] as
     the eighth, as claimed: 709, 531, 177 + 1. The hill-climb files hold 46 layouts, 22 certified, 24 not. The nine
     full-LP outputs are as quoted.

6. **The full directed LP and its rigorous bound.** I read `layout_cut_pdlp.py`: one flow per terminal pattern under
   shared arc capacities is exactly the directed cut condition (max-flow min-cut), and the bound
   b·y + Σ min(0, reduced cost)·U with clipped signs is valid for any multipliers, given that some optimal solution
   lies in the box [0, Wm] (it does: capacities can be capped at W, flows taken acyclic). The directed rows are the
   ones Lemma 4.2's proof gives (entering moves ≥ inputs whose end is inside and start outside).

7. **Γ_2 and Γ_3 values** (my LP, job_a and job_g). Block 2 gap layout [(0,7),(3,1),(4,10),(0,10),(10,7)]: 9
   (directed, undirected and all-cuts exact). Block 3: [(4,123),(25,33),(78,48)] 9 (undirected and directed),
   [(34,15),(65,28),(126,113)] 9, [(114,97),(17,41),(16,29)] 9, triangle plus pendant [(0,95),(0,80),(80,95),(84,60)]
   11. Block 3, W = 5 (directed): 13 and 14, as quoted. Three-end image of [(4,123),(25,33),(78,48)]: 8.5; of
   [(114,97),(17,41),(16,29)]: 9. My own block-3 sample (job_g): all 20 have cut value exactly 9 (my LP, rigorous bound 9.000000).

7b. **Block 4 on all of Γ_4** (my graph, my flow model, HiGHS PDLP, my own bound; job_c .. job_u):
   - nc-gain's [(331,751),(1319,596),(837,1071)] (e = 3): directed 11.999999 (bound 11.999999), as quoted (12);
     undirected 12.000000 (bound 11.999999). So this one holds for every orientation of the pairs.
   - nc-gain's 16 unsent hard layouts (images 10.5 to 11.0): every one 12, with bounds 11.999997 to 12.000000. Directed: [(787,1231),(576,1200),(1677,534)] and [(954,702),(1849,1516),(1697,1353)] (both e = 3). Undirected: [(422,72),(2002,1363),(1426,1720)], [(626,1671),(1513,184),(510,884)] (both e = 3), [(270,639),(1060,611),(2028,1218)], [(1991,864),(1730,2166),(1997,1197)], [(2272,2070),(1433,1339),(891,934)], [(220,1768),(2005,1343),(1231,1113)], [(974,904),(259,1045),(971,1410)], [(492,513),(1511,762),(1883,1768)], [(1538,1305),(2077,511),(736,1392)], [(1659,1830),(390,143),(1585,487)], [(1907,1125),(1625,789),(878,530)], [(2093,1504),(1307,1407),(1020,960)], [(2188,1587),(1972,1363),(1657,1877)], [(842,730),(1084,234),(1769,1372)]
   - W = 4 triangle plus pendant (three mutually opposite frames a, b, c with pairs a–b, a–c, b–c; pendant (x, y)
     with x next to b, y next to c, x opposite y; not two-class stretchable), undirected:
     [(0,2039),(0,2084),(2039,2084),(1204,750)] 14.000000 (bound 13.999999) and
     [(0,2197),(0,2224),(2197,2224),(1245,1526)] 14.999999 (bound 14.999999). So with four arrays cuts cannot
     exclude a coherent 4 in 14 or 4 in 15 at block 4 on these layouts.
   - My own random image screen at block 4 (job_h): 80 random non-compatible layouts, drawn from 326 (75.5 % were compatible). Small images certified 78 of them (11.25 to 12). The other two stopped at exactly 11.0, and both have cut value 12 on all of Γ_4: [(1477,944),(1123,424),(1700,2056)] (e = 2) and [(2220,467),(859,547),(468,334)] (e = 1, undirected, bound 11.999999).

8. **The entropy LP encodings.**
   - My own Shannon LP (`calc/vshannon.py`, no shared code): path 3; directed butterfly multicast 7 against
     routing 8 (unchanged by Ingleton rows); butterfly as two unicasts 6 (no false gain); combination network
     C(3,2) (two symbols to three sinks) 9 against routing 10, a second known coding gain (10/9); exchange on one
     edge 2.
   - nc-gain's `nclp.py` on the combination network: 9, correct.
   - On Γ_2 itself, I rebuilt session 1's three-cell contraction from my own graph (`calc/vshannon2.py`, job_k2): a
     line through frame 0 (3 frames), the frames opposite frame 0 (8), and the rest (4). Grouped arcs and sessions
     (rates 16, 8, 8, 16, 16, 16), full Shannon rows, and closure rows for every union of cells give 80.000000,
     which is 1.07143 per session after averaging, as nc-gain found.
   - nc-gain's Ingleton row: never positive on 150 random GF(2)/GF(3) rank functions, and +1 on the Vámos matroid
     (the classical violation). Its DFZ (65) row: never positive on 120 random GF(3) rank functions, 0 at most on
     the seven vectors of the Fano configuration over GF(3), and +1 on the Fano matroid over GF(2). So it is a
     characteristic-not-2 inequality encoded correctly (`calc/vineq.py`).

9. **Translation to savings.** s = ln(Wm/M)/ln m reproduces every quoted figure: 0.415, 0.491, 0.547 (cuts on the
   uniform layout), 0, 0.0628, 0.1714 (certificates), and the table's thresholds m^0.1, m^0.01. No known unit is
   below a stated ceiling: the known savers (E8 unit, U8X0HQ, T74, the block-110 unit) are at blocks 45 to 110 and
   pack frames (fp-measure: 5 to 55 extra labels per paid move); the small exact units of fp-measure all have Δ ≤ 0;
   the block-2 stand-ins are not engine units (right in one representation only, v-slack), so Lemma W does not
   apply to them.

## Corrections to wording (none changes a number)

1. **"Cut gaps exist with W = 4, so cuts cannot settle 4 in 13..15"** (ceilings; NOTE item 5: "cuts alone will not
   settle the four-array targets"). When nc-gain wrote this, it had a W = 4 gap only at block 3 (11 against 12, below
   the floor there). A block-3 gap does not show a block-4 one. My two block-4 layouts with cut values 14 and 15
   now support the sentence for 4 in 14 and 4 in 15. For 4 in 13, a layout with cut value at most 13 is still
   needed, and none is known. Suggested wording: "cuts cannot exclude 4 in 14 or 4 in 15 (block-4 layouts with cut
   values 14 and 15, v-nc-gain); 4 in 13 is not settled either way".
2. **"No layout reached leaves room for a coherent 3 in 10 or 3 in 11"** and NOTE item 6 ("as every one I checked
   does"). As nc-gain filed it, 16 hill-climbed layouts stopped at image values 10.5 to 11.0. Those values exclude
   3 in 10 (M ≥ 11) but not 3 in 11. "Reached" and "checked" should read "settled". After my runs this point is
   closed: I solved all 16 on Γ_4, and every one has cut value 12, so the sentence is now true as written for every layout either of us reached. With nc-gain's 19 solves, my 16 and my 2 random ones, 37 three-pair layouts at block 4 have been solved on all of Γ_4, and every one has cut value 12.
3. **Orientation.** nc-gain's nine full-Γ_4 solves use the directed LP, which holds for the sampled orientation
   (each input goes from the first frame of its pair to the second; Lemma 4.2's proof gives exactly these rows). The
   image certificates are undirected and cover every orientation. The text should say which one it uses. My
   undirected solves cover every orientation: [(331,751),...] and the unsent layouts I ran undirected.
4. **Minor.** On [(1451,1722),(338,243),(579,532)], nc-gain's stored exact value (11.3043) is a little below the image
   LP's optimum (11.3333). It is still valid and still above 11, so nothing changes.

## Files

- `calc/vgamma.py` (my Γ_m), `calc/vlemmas.py` (index map, claim (b), Lemma W and W+ checks, e(D)),
  `calc/vsym.py` (symmetry, uniform cut at block 2), `calc/vcut.py` (my cut LPs and bound), `calc/vcut_tests.py`
  (block 2 all-cuts tests), `calc/vrun.py` (Colab runner), `calc/vshannon.py`, `calc/nclp_comb.py`, `calc/vineq.py`,
  `calc/vmisc.py`, `calc/vseed18.py`, `calc/vgen.py`, `calc/vtri.py`; outputs `calc/out_*.txt`.
- Colab jobs (each folder: `tasks.json`, the programs and graph files sent, output in `out.txt`):
  job_a (Γ_2/Γ_3 values), job_b (17 block-4 image values), job_c / job_d (nc-gain's [(331,751),...] directed and
  undirected), job_e, job_f, job_l .. job_u (the 16 unsent hard layouts), job_g (my block-3 sample), job_h (my block-4
  image screen), job_v, job_w (its two layouts left at 11.0), job_i, job_j (W = 4 triangle plus pendant at block 4),
  job_k2 (Shannon LP on the three-cell contraction of Γ_2; job_k, the same LP by simplex, was stopped as too slow).
  All jobs are finished; nothing of mine is running on the Colabs. Nothing was posted.

## What this leaves open

- At block 4 with three arrays, 37 layouts have now been solved on all of Γ_4 (nc-gain's 19, my 18). Most were
  hill-climbed toward the worst e(D), and every one has cut value exactly 12. A proof for the non-compatible layouts
  (nc-gain's next step 2) is still missing, so "every coherent three-array unit needs 3m" is evidence, not a theorem.
- With four arrays at block 4, cuts are not enough: two layouts have cut values 14 and 15. Excluding a coherent
  4 in 14 or 4 in 15 needs more than cut rows, or a proof that such layouts cannot come from a unit. 4 in 13 is open.
- Everything here is about coherent units. It reaches every unit only through Proposition 6.
