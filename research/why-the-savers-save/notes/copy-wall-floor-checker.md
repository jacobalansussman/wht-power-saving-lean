> Working note written by an AI agent during the round of 13 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# u-floor: the copy wall made exact (13 October 2026)

Standard library only (Python 3.9 on this Mac). Run everything with `python3 -I`.

## The statement

A unit has v ports. Port T has a frame label a_T in K^h (K = Q on the bit side, F_2 on the complex side) and the
form N on K^h. Source T starts at frame <a_T>, and target S lives at frames inside a_S^perp. Target S *knows*
source T when N(a_S, a_T) = 0, because it can then cancel x_T at a shared frame with no copy. A *decoder* over the
payload field F_p is a set of terms (r_j, g_j). Term j is one scratch *total* r_j.x, and it is added with weight
g_j[S] into target S. The decoder must *fit*: with B = sum_j g_j r_j^T, B[S,S] != 0, and B[S,T] = 0 whenever S does
not know T.

A total carries sources whose labels span P_j = span{a_T : r_j[T] != 0}. By the frame sandwich (hub C.1, #446
Theorem A) it sits at a frame containing P_j. It can reach target S without a further drop only at a frame inside
a_S^perp. So it has to be copied and dropped.

- A **copy child** is the paid block that performs that drop. A child of rank r costs r copy moves.
- A **star** decoder drops each total in one child, from P_j all the way to frame 0. The child's rank is
  beta_j = dim P_j = k * w_in, where w_in is the **input width** (the dimension spanned by the labels read) and k is
  the label rank (k = 1 here).
- **cst** is the sum of beta_j per invocation and stage. **x = cst / v** is the copy charge per port. It enters
  the deficit D = T[(s-1) v k - s cst].
- The **star floor** is the least cst over fitting decoders with star copies. The weaker **pairing floor** (#446
  Theorem C) charges each term only rank N[P_j, span a_{served}] = dim P_j - dim(P_j cap common cut), which is a
  drop to the common cut instead of to 0. Term by term, pairing <= star <= w_in.

The hub's wall, and what this folder recomputes from scratch (`python3 -I copyfloor.py e8`, 4 s):

| unit | decoder | star | pairing | hub |
|---|---|---|---|---|
| bit P120 (E8, Q frames) | 9 totals = columns of U, U U^T = B0 = (1+G) mod 2 | 72 | 72 | 72 |
| bit P120 | every one of the 511 functionals q reads a set spanning Q^8 | min 8 | | "all 511 span" |
| bit P111 (9 lines with independent u_y punctured; our line numbering) | greedy min-weight basis of functionals | 63 = 9 x 7 | 63 | 63 |
| complex P112 (frame decoder, the hub's cxp.py formula) | 8 totals A_a = lines not orthogonal to 2e_a | 64 = 8 x 8 | 64 | 64 |
| E7 (63 lines of E8 orthogonal to a line), bit-type | 7 label totals | 49 = 7 x 7 | | x = 0.78 (hub screen) |

## Hypotheses the floor uses (the to-do list for breaking it)

Each one is a door: dropping it either moves the floor or would need a new proof.

1. **Scalar decoder.** Each port's symbol is one register, copied whole (block length 1). Vector-linear decoders
   (repeated labels, t columns per port) are outside the hub's Lemma Q. In index coding, vector-linear codes can
   beat scalar ones. Test: `--dup t`.
2. **Payload field F_2** on the bit side. The 9 is the F_2 minrank of the E8 pattern (8 is UNSAT). Over another
   field, the minrank of the same pattern is bounded below only by the largest set of pairwise non-orthogonal
   lines. Test: `--p 3` (Dijkstra, small v).
3. **Rank-one labels** (k = 1). X depends on (s, h/k, x/k) only.
4. **Fitting = the endpoint (X, Y, H) -> (X, Y + X, H)** with targets inside a_S^perp ("y avoids u_t", hub C.2 and
   C.3). Complement sinks are legal, but they are forced on today's words.
5. **Side information only through frames.** S cancels x_T only when N(a_S, a_T) = 0. Helpers, timing and B_K
   windows decide whether that cancellation is actually realised. The floor treats it as free, so this hypothesis
   only makes the floor generous.
6. **Star shape.** One child per total, dropped to frame 0, so beta_j = w_in. The pairing floor relaxes this to
   the common cut. On every E8 unit checked (hub s2-decoder) the two are equal, because the targets served by any
   total span the whole space.
7. **Totals are formed at P_j** from sources at their labels (the sandwich). The floor assumes no total is
   assembled from content that has already been dropped.
8. **Copy trees, not DAGs.** Theorem C is stated for copy trees. On a modular lattice a single-content tree costs
   exactly a drop to its meet. Merging two dropped packets and dropping the sum again is dominated in the cases
   worked by hand (see the report), but this is not proved in general.
9. **Nesting.** Frames climb along a register's life, and drops happen only as copies. A non-nested step f -> g
   costs d(f, g) = dim f + dim g - 2 dim(f cap g).
10. **Common frame at every gate** (spec rule 1).
11. **Read-only sources.** Sources are restored and only copied from. No source is temporarily turned into a total.
12. **The decoder pattern.** B0 is unique at window <= 9 (hub C.6). A decoder that differs off the pattern needs
    window >= 10, which the floor does not cover.
13. **One decoder per invocation, paid at every stage.** That gives s * cst in D. Copies shared between stages,
    replicas or calls are not considered.
14. **The unit.** The E8 line system, punctured. A different unit, or one with smaller x, moves the floor (hub open
    problem 2).
15. **Pricing.** One copy child per total, priced as a child of rank beta_j in the pressure identity (P) with
    batching (Bat).
16. **Whole-array frames.** One frame per register. Coordinate copies, where part of a register drops, are
    excluded.

## What the exact search shows (decoder level, all exact)

`python3 -I copyfloor.py floors <instance> [--p P] [--dup t] [--pair]`

- Star floor over F_2: an exact branch and bound. The star cost depends only on the read vectors, and a family
  of read vectors admits fitting served vectors iff, for every target S, their span contains b with b_S = 1 and
  b_T = 0 for every T that S does not know. So the star floor is the least weight of an independent read family
  that passes this test.
- Pairing floor, and any p: Dijkstra over the constrained entries of B. Practical for v <= 5, and slow above that.

Results are in `out/floors_small.txt` and `out/floors_more.txt`. For every tiny unit tried (A2, B2, A3, the odd
vectors of F_2^3 and F_2^4) both floors are v, which is one drop per port. Decoders gain nothing there, so x = 1
and the copy question is trivial. The smallest root-system unit with x < 1 is E7 (v = 63), and a plain exhaustive
search does not reach it. The model therefore splits into two parts:

- tiny units, where the star and pairing floors are exact and equal (and equal v);
- E7 and E8, where the checker evaluates any candidate decoder exactly but the minimum is known only within the
  label-decoder family (hub SAT).

## Testing a candidate decoder

Write a JSON file `{"p": 2, "dup": 1, "terms": [{"r": [...v entries...], "g": [...]}, ...]}` and run
`python3 -I copyfloor.py eval <instance> decoder.json`. It prints whether the decoder fits, its star cost and its
pairing cost. `example_decoder.json` is a worked example on `f2h3-odd` (fits, star 4, pairing 4). For E8-sized units,
import the module and use `Unit(...).delivered / fits / cost` directly, as `e8_calibration()` does.

## Runs recorded on 13 October 2026 (`out/`)

- `e8_calibration.txt`: 72 / 72 (bit P120), 63 / 63 (bit P111), 64 / 64 (complex P112, the hub's formula), and
  E7 at 7 totals x 7 = 49.
- `floors_small.txt`: A2 3, B2 4, f2h3-odd 4 (also 4 with payload F_3), A3 6, f2h4-odd 8. In every case the star
  floor equals the pairing floor, which equals v.
- `floors_dup.txt` (vector-linear, t copies of each port): A2 x2 6, f2h3-odd x2 8, so no gain per port. B2 x2
  and A2 x3 did not finish within the 60 s limit.
- D4 (v 12) and cube13 (v 13) did not finish within the 120 s limit: the plain branch and bound is too slow above
  v of about 10. The next step is a better bound, or an ILP once a solver is installed in a venv here.
