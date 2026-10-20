<!-- header: frame-packing-two-arrays -->
# The skeptic's step-by-step check of the single-term game

**About this file.** This is the working proof file as an AI agent wrote it during the round of 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent. Its opening lines say that no skeptic had read it; that was true when it was written. What happened afterwards is below, and where a later skeptic's verdict differs from this file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean. File names inside it refer to the programs and outputs now in `../programs/` and `../data/`; [README.md](README.md) maps the working names.

**Status now.** This is the skeptic's own record. Its verdict: not refuted; Proposition 4.4 (no single-term 2 in 7 in the block-4 cube) re-run independently on all 107.0 million six-move states.

---

# v-p6-refute: step-by-step check of p6-refute's second session (proofs, computations, corrections)

15 October 2026, 17:39 to about 18:50 EDT (Mac clock). One skeptic. Labels: [checked] (I re-derived the step or
re-computed it with my own code), [corrected] (the step is wrong or overstated as written; the fix is given),
[not checked]. Everything I ran is in this folder; nothing was posted and nothing is in Lean.

Object under review: `../p6-refute/PROOF.md` section 4, `../p6-refute/NOTE.md` (second-session section),
`../p6-refute/STATUS.txt` (session 2), `../p6-refute/stg.py`, and the Colab outputs `job12` to `job19`.

## 1. The model: is the single-term class really the position game that stg.py plays?

The class: two inputs, any number of scratch registers, every register a paid move acts on holds at most one term of
each input. p6-prove's Lemma M turns a unit into the module game V <- V + C[B] p v (v single-term), and Lemma G
describes the single-term vectors of V through the cover graph.

- **Lemma G, position level** [checked, by hand]. Cover edges over one copy form a matching on the lifts (at a cover
  vertex (2, h) the only edge of a copy with twist t comes from (1, h t^-1)). So a closed walk in a TREE component
  of the position graph lifts to a backtracking walk: every cover component over a locked tree contains exactly one
  lift of each position, and nothing in it is available. A cycle in the position graph may or may not lift to a
  non-degenerate cycle. Treating every cycle as an unlock ("permissive") therefore over-approximates the module and
  never under-approximates it. One lineage's copies form a matching, so it cannot make a cycle by itself; two
  lineages with the same twist at the same position make a 2-cycle (a real unlock when their ratios differ).
- **Lemma G, completeness** [checked, `vmod.py`]. My own explicit modules in C[G]^2 (annihilator basis, every
  single-term vector of V found by zero or parallel columns), against my own position-level prediction (singles at
  available positions; pairs with both ends available; pairs inside a locked tree with the path product as twist).
  Random single-term plays, biased to pairs, in D_5, D_7, D_9, S_4/Z_3, B_3, J(4,2), H(2,3), the block-4 cube B_4 and
  H(3,3): 2,732 checks, **0 missed, 0 over-predicted**. Controls: without the cycle rule 8 to 16 misses per 300
  checks; without composites 5 to 11. So the check discriminates, and no single-term vector outside the four
  move types (single, copy, composite of a path of three or more edges, created pair) exists in V.
- **stg.py implements that game** [checked, `vcross.py`]. On one game object shared by both programs (built from my
  own groups), for every state of my exhaustive search: stg's closure = my closure on every successor; stg's
  successor set = my unpruned successor set, up to the parent state itself; stg's `goal_in_one` = brute force over
  my unpruned moves, both with the real target and with an artificial lower target (so that goals are frequent);
  every brute-force goal state passes bfs3's level filter. Games and depths: D_7 6, D_5 4, J(4,2) 4, H(2,3) 4,
  S_4/Z_3 4, D_9 7, B_3 5, B_4 3. **0 disagreements** (B_3 alone: 2,977,092 closure comparisons, 94,290 goal
  tests).
- **Lemma 4.1 (symmetries)** [checked]. Right multiplication of one input's entries by a free element keeps
  V_0 = C[B]^2 and commutes with left multiplication, so it maps plays to plays; levels and the target orbit are
  invariant. stg's key is the transformed state itself, so equal keys mean equivalent states. `vsym.py`: 1,800
  random symmetries on states of B_3, S_4/Z_3, H(2,3), D_7, B_4, J(4,2): canon invariant, successor classes equal,
  L_1 + L_2 and the goal preserved, 0 failures.
- **Lemma 4.2 (pruned moves)** [checked]: the move-set comparison above is exactly this statement.

## 2. Lemma 4.3 (cube and Hamming games inside the block-m frame graph)

- [checked, `vlem43.py`, actual symplectic matrices over F_2]. Cube m = 2, 3, 4 and Hamming m = 2, 3, 4: every
  group element is a symplectic matrix and the point action is a homomorphism; the sub-game's free group equals
  sub-group ∩ Stab(L0); P g -> g^-1 L0 is well defined and injective; coset-graph distance = m − dim(A ∩ B) for every
  pair of positions; every paid letter has d(L0, h L0) = 1; the twist action P g -> P g t equals frame -> t^-1 frame;
  the level-m positions lie at frame distance m (the cube's top is span(e_i), not L_T; the Hamming top orbit
  contains L_T). 0 failures.
- So at the **position level** (which is all the permissive single-term game uses) both are sub-games of the block-m
  single-term game, with the same moves and the same permissive closure. A violation of (I) found there whose
  releases are all by a loop or a catalyst (no cycle) would be a real engine single-term saver; one that needs a
  cycle would still need an exact check, because the exact rule depends on the cycle's holonomy in the real group.
  This is how PROOF 4.3 should read; as written ("the same moves, availability and levels") it is right for the
  permissive rule only [corrected, wording].
- **[corrected] The multi-term remark is wrong as stated.** NOTE ("a scratch 2 in 7 in the cube game with general
  contents ... would be an engine (frame-game) 2 in 7"), PROOF 4.3 and the summary's item (5) and next step (c)
  treat a multi-term unit of the cube game C[B_4] as an engine unit. It is not automatically one:
  - B_m does not lift to the engine's group G. G lies in GL(2^m, Q(i)) (all free letters and C have entries in
    Q(i)), so the Hadamard is not in G; the elements of G with the Sp-image of an e_i <-> f_i swap are n·h with
    h = e^(iπ/4)H, h^2 = i·I, n in the Paulis times {1, i, −1, −i}, and (n h)^2 is never 1 (it is ±i, ±Y-type or
    ±i·Y-type; computed by hand for n = I, X, Y, Z times scalars). The same holds for the Z <-> Y and X <-> Y
    apartments (C^2 = X, S^2 = Z). So a relation such as swap^2 = 1, which a C[B_4] identity may use, is false in
    C[G].
  - A C[B_4] unit IS a unit of the Sp(8,2) frame game C[Sp(8,2)] (B_4 is a genuine subgroup). But C[Sp(2m,2)] is
    C[G/N] (N = Paulis times scalars), the summand of C[G] where every Pauli acts trivially; it sits inside the
    f ≡ 0 (mod 4) block, not "the f = 1 mod 4 part" as NOTE says. An engine unit must hold in every summand. So a
    multi-term cube 2 in 7 would be a frame-game 2 in 7, i.e. an image-level object, not an engine 2 in 7.
  - The single-term results are not affected: they use positions only, and positions are the same in C[G],
    C[Sp(2m,2)] and C[B_m].

## 3. The exhaustive results

- **Outputs** [checked]: every output named in the summary exists and says what is claimed. job16 parts:
  1,395 + 1,395 + 1,394 × 4 = 8,366 layer-4 states (the full layer); depth-6 successors sum to 93,790,253; filtered
  states sum to 1,161,509; goal-in-one "no" in all six parts. job14: 4 × ~22.8 M = 91.2 M depth-6 successors, no
  goal at 6. job13: 55.57 M + 55.66 M (H(4,3) depth 5). job12: H(3,3) depth 6, 79,829,304 successors, (I) holds,
  goal first at 6. job19: 19 m = 3 cases, each "goal-in-one no" (out_f is empty because those groups have no m = 3
  class, as RUNLOG says). B_3 to 6 and D_9 to 8 in job10. Sp(4,2) to 3 in job10/outB_final.
- **bfs3's filter and goal test are admissible** [checked]. A move creates frames only next to existing frames of
  the same input (single: next to an available frame; copy: x' next to x and t(x') next to t(x); created pair and
  composite likewise), and a release only makes existing frames available; so a goal at depth D needs a frame of
  level >= m − 1 of each input at depth D − 1. The goal test's skip rule is exact because a new generator merges at
  most the components of its own ends. Both were also confirmed by brute force (section 1).
- **Independent re-computation** [checked, my own solver `vstg.py`/`vb4.py`: unpruned moves generated from group
  elements, global dedup (each state at its least depth), my own canonical form and my own goal-in-one test by
  generator ends, itself checked against brute force on every tested state with real and lowered targets]:
  - (I) at every state, goal first at 2m: D_5 (4), D_7 (6), D_9 (8), J(4,2) (4), H(2,3) (4); S_4/Z_3 to 5 and
    B_3 to 4 without goal (exhaustive, no symmetry).
  - Proposition 4.5 (no single-term 2 in 5 at m = 3): reproduced in all 19 cases (D_7 three classes, S_4/Z_2,
    S_4/Z_3, B_3, A_5/Z_2 class 2, A_5/Z_3 class 0, A_5/Z_5 two classes, A_5/V_4 three classes, S_5/D_5, S_5/S_3,
    A_6/3^2:2 two classes (Colab, `vp6r_a6/`), J(6,3) two classes; my class enumeration gives the same classes and
    the same m). Calibration: the same decider finds the plain 2 in 6 in A_5/V_4, A_5/Z_5, S_4/Z_2, A_5/Z_2, D_7.
  - Proposition 4.4 (block-4 cube, no single-term 2 in 7): see section 5 (Colab run vp6r_b4).
- **Conjecture S, my own test** [checked as data, `vinv.py`]: every reachable state (least depth) of D_7 (6),
  D_9 (8), D_5 (4), J(4,2) (4), H(2,3) (4), S_4/Z_3 (5), B_3 (4), S_4/Z_2 (5), A_5/V_4 class 1 (5; 917,930 states):
  maximum excess 0. Pairs witnessed by one lineage only reach excess 3 (D_7, S_4, B_3, A_5) and 4 (D_9 at 8 moves),
  more than the "+2" of `../p6-refute/PROOF.md` §4.6 (which came from shallower runs); this is consistent with
  Lemma 4.7, not a counterexample.

## 4. The small lemmas

- **Lemma 4.7** [checked; the proof is right]. But it is proved **given Conjecture S (and T) at earlier times**:
  the root bound lv(s_1) + lv(s_2) <= t_R − 1 is Conjecture T for a created pair and S's "two ends of a locked chain
  of three or more copies" clause for a composite. PROOF.md says so; the summary lists it under "Proved" without the
  hypothesis [corrected, wording]. As an induction step inside a proof of S it is fine; on its own it bounds
  nothing.
- **§4.5 (the natural crossing twist ties)** [argued in the source, not a proof]: the creation bound uses (I) at the
  creation time, so it is circular as a proof of (I); it is labelled [argued], correctly.
- **Correction to session 1 (parity)** [checked]: Lemma D's determinant argument needs every array, scratch included,
  to end with the transform (square words); with free scratch it gives nothing. The correction is right.
- **Equivalence (I) <=> "single-term 2 in 2m − e"** and the sizes: 2 in 7 at block 4 saves ln(8/7)/ln 4 = 0.096;
  2 in 8 at block 5 saves 0.139 [checked].

## 5. The independent block-4 cube run [checked, `vb4.py`, Colab job `vp6r_b4/`]

My decider: unpruned moves and closure from `vstg.py`; my canonical form; global dedup (each canonical state at its
least depth) to depth 5, split six ways at layer 4; every depth-6 successor closed and checked for (I) and the goal;
then, for every depth-6 state whose inputs both own a frame of level >= 3, an exact goal-in-one test by generator
ends. Coverage: every state of least depth 6 has a predecessor of least depth exactly 5, so the depth-6 phase covers
all of them; a goal at least depth 7 has its six-move prefix there.
- Layer 4: 8,074 canonical states (1,346 × 4 + 1,345 × 2 over the parts); layer 5 per part 125,981 to 128,063.
- Depth 6: (I) holds on 17,949,286 + 17,949,852 + 17,852,110 + 17,651,571 + 17,858,536 + 17,742,341 = 107,003,696
  successors; no goal at 6.
- Depth 7: 262,959 + 245,345 + 273,404 + 223,219 + 193,943 + 240,735 = 1,439,605 states tested; goal-in-one 0.
- **Proposition 4.4 is reproduced.**

Also: `vbeam.py` in H(4,3) to depth 7 (width 300, two scorings): no (I) violation (a search, not a proof).
