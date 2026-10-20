> **Status now (19 October 2026).** This is the round's working report as the AI agent wrote it (private paths removed). Its first lines say no skeptic had read it; afterwards one skeptic upheld it with wording corrections ([entropy-bounds-skeptic.md](entropy-bounds-skeptic.md)): the four-array cut gaps at block 4 that item 5 speaks of were found by the skeptic, not here; "reached" and "checked" in item 6 and section 8 should read "settled"; and the full Γ_4 solves here use the directed LP, which holds for the sampled orientation. The skeptic then solved the 16 hard layouts this report did not send to the full LP (all have cut value 12). Where this report and the skeptic disagree, the skeptic's later verdict wins. Files it names that are not in this repository (run logs, status files, briefs, notes to other rounds) are the rounds' private working files.

# Can entropy bounds on network coding beat our certificates for coherent units? (nc-gain)

Packing round, added track, 15 October 2026, in two sessions: 11:31 to about 14:00, and (after a restart) 15:20 to
about 18:00 EDT (Mac clock). The second session finished the first one's pending runs and added section 8 (few
arrays) and items 5 and 6 below. Proofs are in `PROOF.md`, the run log is in `RUNLOG.txt`, the dated predictions are
in `STATUS.txt`. No skeptic has read this. Nothing was posted.
Heavy LPs ran on the Colabs (at most two of my processes per Colab); the Mac ran only short scripts. Labels:
[proved] means a proof is in PROOF.md; [computed] means a program and its output; [search] means a heuristic search
(what it finds exists, but it may miss better); [cited] means a published theorem cited from memory; [read] means
taken from one of our notes.

## The short answer

1. **No. On the frame graphs Γ_2, Γ_3 and Γ_4, every entropy LP small enough to solve proves no more than the cut
   bound, and our certificates are far stronger.** Coherent units (more generally, units with no frame packing on
   live moves in some gauge) are linear network codes on Γ_m, with one session per array from its source frame to
   the opposite frame.

   | block | cut bound (moves per array) | best entropy-LP bound within reach | our certificates (all units) |
   |---|---|---|---|
   | 2 | 1.5 | 1.5 | 2 (no coding gain) |
   | 3 | 1.75 | 1.75 | 2.8 |
   | 4 | 1.875 | 1.875 | 3.154 (41/13) |

   As ceilings on the saving of coherent units at that block, the cut and entropy columns allow 0.415, 0.491 and
   0.547. Our certificates allow 0, 0.063 and 0.171. The entropy route improves no ceiling.

2. **Why.** The full LP on Γ_m has one variable per set of sessions and arcs: 2^210 sets at block 2, before
   symmetry. Any LP small enough to solve must lump frames together, and I prove that lumping caps the result
   [proved, PROOF.md Lemma 3]. Whatever inequalities the LP uses (Shannon, Ingleton, DFZ, or any other valid for
   linear codes), it can never certify more than plain routing would cost on the lumped network. On these graphs
   that routing is cheap:
   - lumping into four groups or fewer never beats the cut bound [proved, Lemma 4];
   - block 2: the smallest lumping that beats cuts has five groups and reaches 1.6. Its LP has 26 variables (2^26
     sets), already beyond reach. Reaching the certificate's 2 needs all 15 frames kept apart [search];
   - block 3: no lumping into fewer than 14 groups beat the cut bound 1.75 in my searches. Reaching 2.8 needs
     about 124 of the 135 frames kept apart [search];
   - block 4: lumpings into up to 300 groups stay below 1.875. Reaching 41/13 needs about 2,027 of the 2,295 frames
     kept apart [search].

3. **The method works on the smallest networks, and already stops short on slightly larger ones** [computed]:
   - Okamura–Seymour network (K_{2,3} with four one-way sessions): cuts allow 6, routing needs 8, and the entropy LP
     proves 8 (rigorous bound 7.9999986; costs are whole numbers). This matches Harvey, Kleinberg and Lehman [cited].
     Over all rates of its four session types it proves "no coding gain" (ratio ≥ 0.99999988).
   - Directed butterfly with two multicast sources: the LP gives 7 against 8 for the best routing, so it detects the
     coding gain 8/7. Paths and two-way exchanges come out exact.
   - But the same K_{2,3} with sessions in both directions (the core of the five-group lumping of Γ_2: 20
     variables, 3.3 million rows) gives only 70.8 against routing 80 (cuts: 64). So even where a lumping could beat
     cuts, the Shannon LP does not reach the lumping's own routing value.
   - Ingleton and DFZ rows changed nothing where I could add them (Okamura–Seymour, a 3-group lumping of Γ_2).

4. **Useful by-products.**
   - **The worst layout for coherent units is the uniform one** [proved, Lemma 1]. In it every ordered pair of
     opposite frames appears equally often. So one number per block caps every coherent unit at that block, and
     fp-capacity's coding-gain term γ_m(h) is at most W(m − L_m), with one L_m for every gauge.
   - **A screen for small fixed layouts** (section 5). The targets have few arrays, where the uniform layout is far
     too pessimistic. For one fixed layout the cut LP is exact and takes about a minute.
     - At block 3, every random layout of five pairs has cut value above 14, so a coherent 5 in 14 cannot use it.
     - But layouts built on three mutually opposite frames reach 14, and one reaches 13. So cuts cannot exclude a
       coherent 5 in 14, and our certificate (14 moves for five arrays) is the stronger bound.
     - At block 4 with three arrays, the layouts solved all have cut value 12 = routing, so none of them allows a
       coherent 3 in 10 or 3 in 11. That is only a few layouts, at 10 to 20 minutes each. (Session 2 replaces this
       sampling by a proof for most layouts: item 5.)

5. **New in session 2: with few arrays, cuts alone exclude a coherent saver on most layouts, and I can prove it.**
   - **Lemma W** [proved, section 8; any block, any number of arrays]. Call a set of input pairs *stretchable* if
     each pair can be oriented so that every start frame in the set is opposite every end frame in the set. If the
     pairs of a layout split into two stretchable sets, every coherent unit with that layout needs W·m moves: no
     saving. The proof is two lines: one distance function per set, and the fact that the plane with the max-norm
     is the plane with the sum-norm.
   - **Three arrays.** If two of the three pairs share a frame, or two pairs are *compatible* (both cross pairs
     opposite, or both parallel pairs opposite), Lemma W applies. That is 74 to 85 % of random layouts at blocks 2
     to 4. In general a layout that falls short of compatibility by e steps needs at least 3m − e moves
     (**Lemma W+**, proved by moving one end of a pair one step at a time). At block 4, e ≤ 1 for about 99 % of
     random layouts, so a coherent 3 in 10 is excluded by proof there.
   - **An exact screen for the rest** [computed]: map every frame to its distances from one end of each pair. This
     shrinks Γ_m to at most (m+1)^W points, and a cut bound on the small picture is a cut bound on Γ_m. A rational
     certificate is checked exactly, and since moves come in whole numbers a value above W·m − 1 means W·m moves.
     - At block 4 I screened 709 random layouts of three pairs with six distinct frames. 531 are compatible
       (proved). Every one of the other 178 was certified: 177 by a small image (exact values 11.125 to 12), and one
       by the full LP of Γ_4 (value 12).
     - Layouts hill-climbed toward the worst case (e = 2 or 3) are harder: the small images certify 22 of 46. I solved
       eight of the others on all of Γ_4 with session 1's exact LP, the two worst included: every one has cut value
       12 (two came out at 11.99999, the solver's accuracy, still above 11). The other sixteen were not sent to the full LP.
   - **What this means for the targets.** A coherent 3 in 10 or 3 in 11 at block 4 (the fp-build target) must use a
     layout with six distinct frames and no compatible pair; none of the layouts I settled allows it (all 19 solved
     on Γ_4 have cut value 12). A builder who wants a three-array saver must pack frames on live moves in every
     gauge, or find one of the rare clustered layouts that nobody has solved.
   - **Four and five arrays.** Lemma W covers 38 % of random four-pair layouts at block 4 (so a coherent 4 in 13 to
     4 in 15 cannot use them) and 12 % of five-pair layouts at block 3 (no coherent 5 in 14 on them). But cut gaps
     already occur with four arrays: at block 3 a triangle of mutually opposite frames plus one more pair has cut
     value 11 against 12. So cuts alone will not settle the four-array targets.
   - **Calibration.** No layout covered by Lemma W ever showed a cut gap (374 layouts at block 2 with three to six
     arrays), and the known gap layouts of section 5 are not covered by it. In random graphs with three pairs at
     the largest distance, all 52 gaps I found were on non-compatible layouts. So the lemma uses only "opposite =
     farthest apart"; whatever closes the remaining layouts on Γ_m must use more of Γ_m's structure.

6. **What these bounds cover.** Every bound here is about *coherent* units: units with a gauge in which no paid
   move acts on an array whose terms sit on two frames (common-frame units). For them they are refund-type
   ceilings in fp-capacity's sense. They bound every unit only if Proposition 6 (`../../prop6/`) holds: "every unit
   costs at least the best common-frame network with the same number of arrays". If every three-pair layout has a cut
   value above 3m − 1 (as every one I checked does), the best common-frame three-array network costs exactly 3m at every block, and
   Proposition 6 for three arrays becomes "every three-array unit needs 3m", which would settle 3 in 10 and 3 in 11
   at block 4 the way the two-array case settles 2 in 7.

## 1. Coherent units as network codes

- Take a unit with an Sp-gauge in which no paid move acts on terms sitting on two frames. Then it is a linear network
  code on Γ_m (sl-combinatorics Lemma 4.1, Remark 4.7) [read]:
  - every frame is a node that keeps everything it receives and combines it for free;
  - every paid move sends one symbol along one edge;
  - input a enters at its source frame and must be read out at the opposite frame.
- The symbols may be matrix rows (Pauli twins are allowed). That makes it a vector-linear code over C, and every
  bound below holds for it unchanged (PROOF.md §0).
- Routing never adds two inputs, and it costs exactly m moves per input. A coherent saver is a code that beats
  routing. The *coding gain* is routing cost ÷ code cost.
- **The worst layout** [proved, Lemma 1]:
  - Sp(2m, 2) moves any pair of opposite frames to any other.
  - Copying a code once for every group element gives a code for the uniform layout, at the same cost per input.
  - So the uniform layout has the largest coding gain, and any lower bound for it holds per input for every layout.

## 2. The entropy LP, and its calibration

The LP has one variable for each set of sessions and arcs: the dimension of the span of their symbols. Its rows:
- Shannon's inequalities (the polymatroid bound);
- the Ingleton inequality and DFZ's inequality (65) (arXiv:1311.4601, Theorem 8.6), added as cutting planes. Both
  hold for linear codes over C; the second needs characteristic ≠ 2;
- the network rows:
  - the inputs are independent;
  - everything sent is made from the inputs;
  - for every set Q of nodes, what enters Q and what starts in Q determine what leaves Q, and what Q must read out.

Each arc direction is a separate variable, so a symbol sent both ways is paid twice. The code (`calc/nclp.py`) uses
symmetry classes of sets as variables, and merges every set with what it determines (the closure reduction, which
halves the Okamura–Seymour LP). It solves with HiGHS. For the large instances it uses HiGHS's first-order method,
together with a rigorous lower bound computed from the approximate dual (PROOF.md §4).

| network | routing | cut LP | entropy LP |
|---|---|---|---|
| path of length 3, one session | 3 | 3 | 3 |
| directed butterfly, 2 sources multicast to 2 sinks | 8 | — | 7 (coding gain 8/7 detected) |
| one edge, sessions both ways | 2 | 2 | 2 |
| Okamura–Seymour K_{2,3}, 4 sessions | 8 | 6 | ≥ 7.9999986, so 8 |
| same, worst case over the rates of its 4 session types | ratio 1 | — | ≥ 0.99999988 |
| 3-group lumping of Γ_2 | 80 | 80 | 80 |
| distance-2 skeleton of the 5-group lumping of Γ_2 (K_{2,3}, sessions both ways, 20 variables) | 80 | 64 | 70.8: **short of routing** |

## 3. Why no solvable LP beats the cut bound on Γ_2, Γ_3, Γ_4

**Lumping caps every LP** [proved, PROOF.md Lemma 3].
- Suppose an LP writes its network rows only for some family of node sets. Its frames fall into the cells of that
  family's Venn diagram.
- Route every session along a shortest path between cells, sending nothing inside a cell. This routing satisfies
  all the LP's rows. It also satisfies every inequality that real linear codes satisfy, because it *is* a linear
  code (on the lumped network; it skips the travel inside cells, which no row of such an LP can see). So the LP's
  value is at most that routing's cost.
- After averaging over the symmetry group (PROOF.md Lemma 2), the bound per input is at most
  ρ = (average distance between the cells of opposite frames) ÷ (share of edges between different cells).
- An LP small enough to solve must lump:
  - each pair of cells joined by sessions needs its own source variable;
  - each pair of adjacent cells needs two arc variables;
  - full Shannon LPs stop at about 16 to 20 variables (2^16 to 2^20 sets).

**How much lumping costs** [proved for four cells or fewer; search beyond]:
- With four cells or fewer, ρ never exceeds the cut bound. Every metric on four points is a sum of cut metrics
  (PROOF.md Lemma 4).
- Best ρ found, by number of cells (programs listed under Files):

| cells | block 2 (cut 1.5, certificate 2) | block 3 (cut 1.75, certificate 2.8) |
|---|---|---|
| ≤ 4 | 1.5 | 1.75 |
| 5 | **1.6** (one frame, the 3 lines through it, the 8 frames opposite it: a K_{2,3}) | below 1.75 |
| 6 | 1.676 | below 1.75 |
| 9 | 1.846 | below 1.75 |
| 13 | 1.849 | below 1.75 |
| 14 | 1.909 (all 105 one-pair merges tried) | 1.784 |
| 15 | 2 (no lumping) | 1.828 |
| 20 | — | 1.933 |
| 30 | — | 2.256 |
| 60 | — | 2.426 |
| 100 | — | 2.729 |
| 124 | — | 2.811 (first ≥ 2.8) |
| 134 | — | 2.966 (one pair merged) |

- Block 4 (2,295 frames, cut 1.875, certificate 41/13 = 3.154):
  - lumpings by distance to a few reference frames, up to 300 cells: at most 1.803;
  - merging pairs one at a time from the full graph: 3.993 after one merge, 3.531 after 100, 3.267 after 200;
  - the first fall below 41/13 comes at 2,026 cells.
- At block 2, no lumping that beats cuts has fewer than 26 LP variables [search]. So even the best-placed lumping
  is beyond full LPs, and if it were solved it could not pass 1.6, against the certificate's 2.

## 4. Ceilings on the saving of coherent units

Saving s ≤ ln(γ)/ln m for a coding gain γ. A 0.1 unit needs γ ≥ m^0.1, and a 0.01 unit needs γ ≥ m^0.01.

| block | γ allowed by cuts | by entropy LPs here | by our certificates | saving ceilings (cut / entropy / certificates) | 0.1 needs γ ≥ | 0.01 needs γ ≥ |
|---|---|---|---|---|---|---|
| 2 | 4/3 | 4/3 | 1 | 0.415 / 0.415 / 0 | 1.072 (excluded) | 1.007 (excluded) |
| 3 | 12/7 | 12/7 | 15/14 | 0.491 / 0.491 / 0.0628 | 1.116 (excluded) | 1.011 |
| 4 | 32/15 | 32/15 | 52/41 | 0.547 / 0.547 / 0.1714 | 1.149 | 1.014 |

A coherent 0.1 unit is not excluded at block 4. It would need a coding gain of about 15 % on its layout, and nothing
computable here rules that out. Other proved facts still apply to it: it needs three or more arrays (Theorem 4.3),
and a layout that is neither single-source nor single-sink (Corollary 4.4) [read].

## 5. Small fixed layouts (the targets have W = 3, 4, 5)

A unit with W arrays that is coherent in some gauge is a code for some layout of W opposite pairs. For one fixed
layout, the cut rows alone already give an exact, cheap LP: `calc/layout_cut2.py`, by constraint generation with
max-flow separation (PROOF.md §7). If every layout of W pairs had a cut value above M, no coherent W in M would
exist. A search can find layouts where the cut value is low; it cannot prove that none exists.

Which targets are still open after the certificates [read]:
- block 3: W ≤ 4 cannot save at all; W = 5 allows 5 in 14; W = 6 allows 6 in 17;
- block 4: W = 3 allows 3 in 10 and 3 in 11; W = 4 allows 4 in 13 to 4 in 15.

**Calibration at block 2** [computed]. There Theorem 2W proves that routing is optimal for every layout.
- Every sampled layout with W = 3 or W = 4 has cut value = routing (6, 8), about 90 layouts each with local search.
- From W = 5 cuts fall short:
  - layout (0→7), (3→1), (4→10), (0→10), (10→7) has cut value 9 against routing 10. It is a triangle of mutually
    opposite frames with two pendant sessions;
  - a W = 6 layout has 11 against 12, matching sl-combinatorics' 1.83 per pair.
- So beyond four arrays a low cut value says nothing about the truth.
- Lumping that W = 5 layout so that the entropy LP could try to close the gap keeps 33 variables or more: out of
  reach.

**Block 3, five arrays (the 5 in 14 shape)** [computed, search]:
- Random layouts of five opposite pairs: all eight tried have cut value above 14 (certified 14.005 to 14.036). So
  no coherent 5 in 14 uses those layouts.
- Layouts shaped like the block-2 failure (three mutually opposite frames a, b, c with sessions a→b, a→c, c→b,
  plus two pendant sessions) can go lower:
  - (0→128), (0→124), (124→128), (8→124), (76→44) has cut value exactly 14;
  - (0→95), (0→80), (80→95), (12→80), (84→60) has cut value 13.
  - Of 30 layouts of this shape, 22 are above 14, 7 are exactly 14 and 1 is 13.
  - The frames of the first two are written out in `calc/layouts_block3.json`.
- So cut arguments cannot exclude a coherent 5 in 14 at block 3, and the certificate (14 moves for five arrays) is
  stronger than cuts there.
- What this does give fp-build: a coherent 5 in 14 must use a layout whose cut value is at most 14.
  `calc/layout_cut2.py` screens a layout in about a minute on a Colab. The screen excludes random
  layouts; it does not exclude triangle-with-pendants layouts.

**Block 4, three arrays (the 3 in 10 shape)** [computed, few layouts]:
- For three sessions, cuts can fail only when the three pairs are pairwise disjoint. Every other shape of three
  demand edges lies inside K4 or inside two stars, where the cut condition already forces routing
  (Papernov–Lomonosov, sl-combinatorics Corollary 4.4(b)) [cited, from memory].
- Ten such layouts have cut value exactly 12 = routing (rigorous bounds, 12.00000 by PDLP): one random, and three
  each of "parallel" (each pair's ends adjacent to the previous pair's ends), "crossed" (the next pair's start next
  to this pair's end) and "mixed". So no coherent 3 in 10 or 3 in 11 uses those layouts.
  - Each solve takes 10 to 60 minutes on a Colab (2.6 million variables).
- Section 8 replaces this sampling by a proof for most layouts and a fast exact screen for the others.

## 6. For fp-capacity

Their Theorem A(d) carries the coding-gain term γ_m(h) ≤ W(m − L_m).
- By Lemma 1, one L_m serves every gauge and every layout: the per-input floor of codes for the uniform layout.
- The best L_m for codes is still the certificate floor: 2, 2.8 and 41/13 at blocks 2, 3 and 4. Cuts give only 1.5,
  1.75 and 1.875, and no entropy LP within reach does better than cuts.
- So network-coding entropy bounds cannot sharpen a packing ceiling at blocks 2 to 4.

The same is in `FOR-fp-capacity.txt`.

## 7. What would change this

- **A non-lumped LP.** Keep every frame separate, and restrict instead which *sets* get variables: a partial
  polymatroid over the full ground set (1,890 arcs and 8,640 sessions at block 3). Lemma 3 does not cap it. It is
  the only entropy route left, and nobody knows how strong such LPs are. I did not try it.
- **A proof that coding on K_{2,3} never beats routing, for every choice of rates.** That would make the 5-group
  lumping of Γ_2 certify 1.6 without an LP. It is still below the certificate's 2, so it would only illustrate
  the method.
- **A smarter argument specific to these graphs.** Theorem R already is one, built on the representation theory of
  Sp(2m, 2). Theorem 2W even proves that the uniform layout on Γ_2 has no coding gain at all.
- **For few arrays (session 2): a proof for the non-compatible three-pair layouts.** Lemma W cannot reach them (one
  distance-like function can stretch two pairs only if they are compatible), and in general graphs such layouts do
  have cut gaps. Every one I solved on Γ_3 has cut value exactly 9, and the ones on Γ_4 that the screen certifies
  have values between 11.1 and 12. A proof would have to use the near-polygon structure, as Lemma W+ does.

## 8. Few arrays: when cuts alone exclude a coherent saver (session 2)

The targets have three to five arrays, and for a fixed layout the cut rows of sl-combinatorics Lemma 4.2 already
give an exact LP (section 5). Session 2 asked when that LP is forced to equal Wm, so that no coherent unit with that
layout saves anything. Proofs are in PROOF.md §8.

**Lemma W (two-class stretch)** [proved].
- Orient each pair (low, high). A set of pairs is *stretchable* if every low in it is opposite every high in it.
- If the pairs split into two stretchable sets A and B, take φ_A(v) = the distance from v to the nearest low of A,
  and φ_B the same for B. Both change by at most 1 along an edge, and φ_A rises from 0 to m along every pair of A.
- The map v ↦ (φ_A(v), φ_B(v)) into the plane with the max-norm moves adjacent frames at most 1 apart and keeps every
  pair m apart. The max-norm plane is the sum-norm plane turned by 45 degrees, so this is a weighted family of cuts,
  and the cut condition then gives M ≥ Wm.
- Only "opposite = farthest apart" is used. Conversely, one distance-like function can stretch a set of pairs only
  if the set is stretchable (a negative-cycle argument), so layouts outside the lemma need genuinely different cuts.

**Three arrays** [proved].
- Two pairs that share a frame are compatible (orient both with the shared frame as low). So every three-pair
  shape except three disjoint pairs is covered, without the Papernov–Lomonosov theorem.
- For three disjoint pairs, Lemma W applies when two pairs are compatible: both cross pairs (start of one, end of
  the other) opposite, or both "parallel" pairs (start-start, end-end) opposite.
- **Lemma W+.** If the best choice falls short by e steps in total (the two cross distances add up to 2m − e), then
  M ≥ 3m − e. The proof moves one end of a pair one step at a time, away from the other pair's start while staying
  opposite its own partner. That such a step always exists is the near-polygon property of Γ_m: on every line the
  point nearest to a given frame is unique.

| block | random layouts (six distinct frames) | e = 0: M ≥ 3m proved | e = 1: M ≥ 3m − 1 proved | e = 2 | worst e found by hill-climbing |
|---|---|---|---|---|---|
| 2 | 3,000 | 85 % | 15 % | — | 1 |
| 3 | 3,000 | 78 % | 21 % | 0.3 % | 2 |
| 4 | 20,000 | 74.5 % | 24.9 % | 0.55 % | 4 |

So at block 4 a coherent 3 in 10 (3m − 2 moves) is excluded by proof on every layout with e ≤ 1, which is 99.45 %
of random layouts; only clustered layouts with e ≥ 2 need the screen.

**The exact screen for the remaining layouts** [computed].
- Map each frame v to ψ(v) = (distances from v to one end of each pair). Adjacent frames go to points at most one
  step apart, and each pair stays m apart. So any cut bound on the image (at most (m+1)^W points) is a cut bound on
  Γ_m.
- I solve the cut LP on the image by constraint generation, take its dual, round it to fractions and check it
  exactly. A value above Wm − 1 proves M ≥ Wm, because moves come in whole numbers.
- If the image built from one end of each pair falls short, the program tries the other ends, then all 2W ends,
  then Γ_m itself.

Results:

| block, arrays | layouts | covered by Lemma W (proof) | others certified by the screen | others solved on all of Γ_m | still open |
|---|---|---|---|---|---|
| 3, three (random) | 100 | 78 | 22, all with cut value exactly 9 | — | 0 |
| 4, three (random) | 709 | 531 (75 %) | 177 (values 11.125 to 12) | 1, cut value 12 | 0 |
| 4, three (session 1's structured layouts) | 10 | 7 | — | 3, cut value exactly 12 | 0 |
| 4, three (hill-climbed to e = 2 or 3) | 46 | 0 | 22 (values 11.125 to 12) | 8 of the other 24, all 12 | 16 not sent to the full LP |
| 4, four (random) | 3,000 | 1,147 (38 %) | not finished | — | the rest |
| 3, five (random) | 2,000 | 239 (12 %) | — | — | the rest |

- "Solved on all of Γ_m" uses session 1's directed cut LP (PDLP with a rigorous bound, 20 to 60 minutes a layout).
  It settled every hard layout I sent (small-image value in brackets): e = 3: [(331,751), (1319,596),
  (837,1071)] (10.75) and [(2115,615), (1327,2101), (641,486)] (10.5); e = 2: [(2030,1063), (1397,1085), (1249,145)]
  (11.0; LP 11.99999), [(487,1544), (722,1494), (2264,721)] (10.67), [(1291,349), (1486,674), (108,1841)] (10.67),
  [(1139,722), (987,1071), (2122,28)] (11.0), [(649,2184), (60,1126), (1376,646)] (10.8; LP 11.99999),
  [(1080,93), (3,906), (632,1269)] (10.67), and the random [(2051,594), (2031,712), (777,952)]. All 12.
- So every three-pair layout at block 4 that I solved on all of Γ_4 (19 of them, counting session 1's ten) has cut
  value 12, and every other one I screened has a certificate above 11. Where the small images fall short, the
  layout is clustered (most cross distances 2 or 3) and the shortfall is the images' coarseness: Γ_4 itself gives 12.
- With four arrays cut gaps do occur: at block 3 the layout made of a triangle of mutually opposite frames and one
  pendant pair has cut value 11 against 12. (At block 3 the certificate floor is already 12, so nothing opens.)

**Calibration.**
- Lemma W: at block 2, 374 random stretchable layouts with three to six arrays all have cut value Wm. The cut gaps
  (five arrays: 5 of 64 other layouts; six arrays: 11 of 114) are all on layouts outside the lemma, as are the
  known gap layouts of section 5.
- Random graphs on 9 to 12 vertices with three pairs at the largest distance: 52 of 1,070 have a cut gap (7.67 to
  8.5 against 9), every one on a non-compatible layout. So Lemma W is sharp in general graphs, and closing the
  non-compatible layouts on Γ_m needs more of Γ_m's structure.
- The screen gives lower bounds, never more than the full cut LP. On the two non-compatible block-3 layouts where
  I ran both, the full LP is exactly 9; the three-end image gives 9 on one and 8.5 on the other (enough, by whole
  numbers), and the six-end image gives 9 on both.

**What is still open.** "Every three-pair layout on Γ_4 has cut value above 11" is proved for 74.5 % of random
layouts (e = 0), checked layout by layout on the 765 block-4 layouts above (19 of them on all of Γ_4, every one
12), and not proved for all. A proof for the
non-compatible layouts would make "every coherent unit with W ≤ 3 needs Wm moves" a theorem at every block,
extending Theorem 4.3 from two arrays to three.

## Files

All in `nc-gain/`.
- `PROOF.md`: Lemmas 1 to 4; calibration, ceiling curves and the translation (§4 to §6).
- `calc/dpg.py`: the frame graphs Γ_1 to Γ_4 (exact).
- `calc/nclp.py`: the entropy LP (orbit variables, closure reduction, Shannon rows, Ingleton and DFZ cutting
  planes, lazy rows, PDLP with a rigorous dual bound).
- `calc/calib_small.py`, `calc/test_cp.py`, `calc/test_closed.py`, `calc/test_pdlp.py`: calibration on small
  networks.
- `calc/quot.py`, `calc/quot2.py`, `calc/quot3.py`, `calc/quot4.py`, `calc/quot5.py`, `calc/small5.py`: lumping
  searches. Outputs: `calc/quot3_m2_prof+agg.json`, `calc/quot3_m3_agg.json`, `job_quot*/`,
  `job_m4/quot5_m4.json`.
- `calc/star.py`: the 5-group lumping of Γ_2 (routing 128, cut LP 112, 26 variables).
- `calc/k23gap.py`: the smallest demand sets on K_{2,3} with a gap between cuts and routing (four sessions).
- `calc/cells_lp.py`: the entropy LP on actual lumpings of Γ_2 and Γ_3.
- `calc/layout_cut.py` (one big LP, too slow), `calc/layout_cut2.py` (exact, constraint generation),
  `calc/layout_cut_pdlp.py` (fast, rigorous dual bound), `calc/layout_search.py`: the cut LP for a fixed small layout
  and searches over layouts; `calc/lay5_g2.py` (the failing block-2 layout), `calc/dump_layouts.py` →
  `calc/layouts_block3.json` (the block-3 layouts with cut value 13 and 14, frames written out).
- Session 2 (section 8): `calc/w3cut.py` (undirected cut LP with its dual), `calc/w3probe.py` (which cuts the duals
  use), `calc/w3ident.py`, `calc/w3compat.py` (compatibility), `calc/stretch2.py` (Lemma W for any W),
  `calc/lemma_calib.py` (calibration of Lemma W at block 2), `calc/w3rand.py`, `calc/w3rand2.py` (random graphs),
  `calc/edef.py` (cross deficiency e), `calc/w3image.py` (image graphs and the king graph), `calc/w3screen.py`
  (three-array screen with exact dual check), `calc/wscreen.py` (any W). Colab jobs: `job_w3s3/`, `job_w3s4a/` to
  `job_w3s4d/` (three arrays), `job_w4s4a/`, `job_w4s4b/` (four arrays), `job_w3full4/` (full Γ_4 LP),
  `job_edef4/` (e at block 4 and the hill-climbed hard layouts), `job_hardpd/` to `job_hardpd5/` (session 1's directed
  LP on all of Γ_4 for the hard layouts), `job_hard1/`, `job_hard2/` (six-end images, stopped), `job_lpd4/` (session
  1's structured block-4 layouts, with `compat_check.txt`); `calc/kingcfg.py` (the king-graph relaxation, too coarse).
- `job_*/`: the Colab jobs and their outputs (``).
- `STATUS.txt`, `RUNLOG.txt`, `FOR-fp-capacity.txt`, `FOR-prop6-and-fp-build.txt`.

## Sources

- Read:
  - `../BRIEF.md`; `../../exist3/BRIEF.md` (standing rules); `../../slack/BRIEF.md`;
  - `../../slack/sl-combinatorics/PROOF.md` and `NOTE.md` (Lemma 4.1, Theorem 4.3, Corollary 4.4, Remark 4.7, the
    cut constants); `../../slack/sl-algebra/NOTE.md`;
  - `../../slack/v-slack/NOTE.md` and `STATUS.txt` (the cut constants 3/2, 7/4 and 15/8, recomputed there);
  - `../../anatomy/NOTE.md`; `../../meta/MAP.md` §0–1.4;
  - `../../proofs/duality/PROOF.md` §0–1 and `dfz_check.py` (the DFZ encoding);
  - `../fp-capacity/PROOF.md` §3 (Theorem A(d)).
  - Session 2: the round brief again; `../../prop6/BRIEF.md`;
    `../../prop6/p6-refute/NOTE.md` (top); `../../slack/sl-combinatorics/PROOF.md` §4 again (Lemma 4.2 in its
    undirected form, Remark 4.3b); `../fp-build/NOTE.md` (top); `../../iso/NOTE.md` §4 and
    `../../../wht33/hitmap/l-bounds.md` §2 (the Li–Li conjecture and the network-coding cases on file).
- Cited from memory:
  - T. C. Hu, Operations Research 11 (1963);
  - N. Harvey, R. Kleinberg, A. Lehman, IEEE Trans. Inf. Theory 52 (2006) (the Okamura–Seymour network);
  - H. Okamura, P. Seymour, J. Combin. Theory B 31 (1981);
  - M. Deza, M. Laurent, *Geometry of Cuts and Metrics* (1997) (on four points the metric cone equals the cut cone);
  - R. Dougherty, C. Freiling, K. Zeger, arXiv:1311.4601 (inequality (65)).
