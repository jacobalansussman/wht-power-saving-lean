<!-- header: frame-packing-two-arrays -->
# Where a counterexample to Proposition 6 can live, and the single-term game (the refute track of the prop6 round)

**About this file.** This is the working proof file as an AI agent wrote it during the round of 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent. Its opening lines say that no skeptic had read it; that was true when it was written. What happened afterwards is below, and where a later skeptic's verdict differs from this file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean. File names inside it refer to the programs and outputs now in `../programs/` and `../data/`; [README.md](README.md) maps the working names.

**Status now.** Section 4 (the second session: the single-term game, Propositions 4.4 and 4.5) was upheld by one AI skeptic, who rebuilt the model, the solver and the large runs with its own code ([single-term-game-skeptic.md](single-term-game-skeptic.md)). Corrections: a multi-term 2 in 7 in the block-4 cube would be a frame-game program, not automatically an engine unit (the cube's symmetries do not lift to the engine's letters); Lemma 4.7 assumes Conjecture S at earlier moments; "same availability" in Lemma 4.3 holds for the solver's permissive cycle rule; one-lineage pairs reach an excess of 4, not 2. Sections 1 to 3 (the first session: Theorem H and numerical frame-game searches) were not checked by any skeptic, apart from one parity correction.

---

# p6-refute: where a counterexample to Proposition 6 can live, and abstract frame games (proofs and computations)

15 October 2026, from 11:59 EDT (Mac clock). One agent; no skeptic has read this. Labels: [proved] (proof here),
[computed] (program in this folder and its output), [read] (from a named file), [argued]. Moves are counted one per
direction. Delta := W*m - M.

Proposition 6 (sl-algebra PROOF.md section 4.3) says: every engine unit with W arrays costs at least as much as the
cheapest common-frame network with W arrays at the same block. A counterexample is a unit U with M(U) < CF(m, W), where
CF(m, W) is the least cost of a common-frame network with W arrays at block m.

## 1. The cells where a counterexample can exist

**Lemma A (no counterexample at blocks 2 and 3 with few arrays)** [proved from Theorem R, read].
- Every engine unit at block 2 has M >= 2W, and at block 3 M >= 14W/5 (Theorem R, `meta/MAP.md` section 1.2 and the
  ceilings table in `slack/NOTE.md`). The plain schedule (cost W*m) is common-frame (sl-algebra Theorem 1, twists e).
- So CF(2, W) = 2W and every block-2 unit costs >= CF(2, W): no counterexample at block 2, for any W.
- At block 3, ceil(14W/5) = 3W exactly when W <= 4 (14W/5 > 3W - 1 iff W < 5). So CF(3, W) = 3W and no counterexample
  exists at block 3 with W <= 4.
- At block 4, ceil(41W/13) < 4W for every W >= 2 (82/13 = 6.31 <= 7), so the floor leaves room from W = 2 on.

*Consequence.* An exhaustive search at blocks 2-3 with W <= 4 can never refute the cost form of Proposition 6; it can
only test its structural form ("every unit is common-frame after a change of basis or gauge"), which p6-prove's item
(D) covers. The cells where a cost counterexample can live are:
1. **W = 2, block m >= 4**, M <= 2m - 1, because CF(m, 2) = 2m at every block (sl-algebra Corollary 2; Theorem 4.3 of
   sl-combinatorics even with scratch). These are exactly the open two-array targets: 2 in 7 at block 4; 2 in 7 and
   2 in 8 at block 5; 2 in 8 to 2 in 10 at block 6 (2 in 6 at block 5 and 2 in 7 at block 6 are closed by mw-big56's
   two-layer bound M >= W + m, unchecked). So for W = 2, refuting Proposition 6 IS finding a two-array saver: the
   proposition opens no new cell there. Of these, 2 in 7 at block 5 and 2 in 8 at block 6 have M = m + 2 = W + m,
   the first size the two-layer bound leaves open; section 3 tests that size in frame games (none found at m = 3
   and m = 4).
2. **W >= 5 at block 3, or W >= 3 at block >= 4**, but only together with a proof that CF(m, W) is larger than the
   unit's cost. No such proof exists for any W >= 3 (cuts give at most 1.5, 1.75, 1.875 moves per array at blocks
   2, 3, 4: sl-algebra Proposition 4, v-slack). For example, a 5 in 14 at block 3 (exactly on the hook floor) would
   refute Proposition 6 only if no common-frame 5 in 14 exists, which nobody can currently decide (nc-gain's
   entropy LP is the only tool aimed at it).

**Lemma B (mixed targets give no smaller cell)** [proved].
Let a unit at block m deliver to array r the element K_(S_r) x_pi(r) up to a free factor, where K_S is the product of
the kernels dir e_j over j in S and |S_r| = k_r (a rank-k_r target). Then there is a full unit (every array K x) with
M + sum_r (m - k_r) moves. So its deficit sum_r k_r - M is unchanged.
*Proof.* Drop the final free factor u_r of array r, apply dir e_j for j not in S_r (the dir letters on different
coordinates commute, so K_([m] \ S_r) K_(S_r) = K), and apply u_r. Each array costs m - k_r more moves. In frame-game
terms (section 3): if g_r has paid length k_r and the target T has length m with g_r on a geodesic to T (l(T g_r^-1)
= m - k_r), append a geodesic word for T g_r^-1; the initial free factors absorb any right P-factor. QED.
- Consequence: a mixed-target two-array saver at block m gives a full two-array saver at block m. At blocks 2 and 3
  both are excluded (Lemma A); for W = 2 mixed targets open no new cell.

## 2. Two facts about the first and last moves

Notation for this section and the next: a unit is A(end) = F_M Y_M F_(M-1) ... Y_1 F_0 = T * I_W in M_W(C[G]) with F_k
in M_W(C[P]) (free mixing; any final free monomial is absorbed into F_M) and Y_k = diag(h_k, 1, ..., 1) (wlog on slot 0,
since swaps are free). pi_T denotes the projection of C[G] onto the span of the right coset P T (the frame L_T); it is
left C[P]-linear, since left multiplication by P preserves right cosets.

**Lemma C (one move never delivers both inputs)** [proved]. In every frame game, with any number of scratch arrays:
let R(t) be the left C[P]-module spanned by all rows (data and scratch) after step t, and Lambda(t) = pi_T(R(t)).
Free steps, mixes, copies and overwrites never enlarge R, so never enlarge Lambda. A paid move on row v gives
Lambda(t+1) contained in Lambda(t) + C[P] pi_T(h v): it adds at most one cyclic module. At the end Lambda contains
(C[P] T)^2, which is not cyclic (a cyclic submodule of C[P]^2 containing (1, 0) and (0, 1) would be all of C[P]^2,
of composition length twice that of C[P]). So at least two different paid moves create the delivered top-frame
content. In particular **2 in 1 (m = 1) is impossible in every frame game, with any scratch**: the single move would
have to create both.
- Direct form, no scratch: T I_2 = F_1 diag(h, 1) F_0 forces row 2 of F_0 into C[P]T intersected with C[P] = 0.
- By time reversal (the inverse word takes T I to I), at least two different moves take content out of the raw
  frame L_0.
- This is the precise obstruction to the "XOR in the air" pattern (two-way relay with broadcast, which saves one of
  four transmissions in wireless networks): a move on an array whose terms sit on several labels at one frame does
  send pieces to several neighbouring frames at once (broadcast), but whatever arrives at one frame from one move is a
  single element, so the two deliveries still need two arrivals, and the two arrivals need their own last moves.
  Separating two deliveries that sit in one array is impossible: if a row holds (alpha T, beta T) at L_T with alpha,
  beta in C[P], then a C[P]-combination giving (T, 0) needs c alpha = 1, hence c is a unit and c beta = 0 forces
  beta = 0.

**Corollary C' (one move per array, every frame game)** [proved]. The argument of Lemma C with W inputs: the
top-frame module must reach (C[P]T)^W, of composition length W len(C[P]), and each move adds at most len(C[P]).
So M >= W in every frame game, with any scratch. For m = 1 this is M >= Wm: **frame games with m = 1 never save**,
for any W. (In the engine the stronger total-rank bound M >= W + m - 1 of ladder Theorem D is on file.)

**Lemma Z (what a 2 in 3 at m = 2 must look like)** [proved]. Every frame game, m = 2, two arrays, no scratch.
Absorb the first free matrix into the input basis, so the first move acts on the row (1, 0) and the target becomes
T G0 with G0 in GL_2(C[P]) (rows of the final state generate (C[P]T)^2 G0). Write the state after the first move as
rows (h, 0), (0, 1); the second move acts on v = (alpha h, beta), the other row is z = (gamma h, delta), with
A = [[alpha, beta], [gamma, delta]] in GL_2(C[P]); the third move acts on y = eps h v + zeta z and leaves
done = eta h v + theta z, with [[eps, zeta], [eta, theta]] in GL_2(C[P]). Then:
- the second coordinate of done is eta h beta + theta delta, a level-1 term plus a level-0 term, and done lies in
  the level-2 span, so **eta h beta = 0 and theta delta = 0**;
- done is a free generator of the final module, so done != 0, hence eta != 0, hence **beta is a nonzero non-unit
  of C[P]** (if beta were a unit, eta h = 0; if beta = 0, then delta is a unit, theta = 0, eta and zeta are units,
  and the second coordinate h zeta delta of the third move's output is a nonzero level-1 term: impossible);
- theta = 0 is impossible by the same chain, so **theta and delta are nonzero zero divisors with theta delta = 0**.
So the row moved by the middle move must carry the second input with a zero-divisor coefficient (an idempotent-type
element such as (1 + p)/2), and the bookkeeping row must cancel the other row's second input through a product of
zero divisors. In particular, **if P is trivial (a Cayley-graph game) two arrays need 4 moves at m = 2**, since C has no
zero divisors. This is the frame-game form of "a two-array saver must split an input" (Corollary 2'): the split is
by idempotents of the free algebra.

**Theorem H (two inputs at distance 2 need 4 moves, in every frame game, with any scratch)** [proved; NOT NEW: it is
the case W = m = 2 of mw-big56's two-layer bound, Theorem T, M >= W + m (`../../moonshot/mw-big56/PROOF.md`, proved
there, unchecked by a skeptic), whose proof is the same top-layer / next-layer argument and whose scope remark says it
needs only a layer function kept by free letters on both sides and raised by at most 1 per paid letter, i.e. it holds
in every frame game. I derived the argument below independently at 13:00 and found Theorem T on file at 13:10.]
Let (G, P, H, T) be a frame game with l(T) = 2 (any diameter, any number of paid double cosets). Every unit for two
inputs, with any number of scratch arrays and any contents, has M >= 4. The same proof works verbatim in the engine at
block 2 with C[B] in place of C[P] (consistent with Theorem R there), and in any game whose paid length is bi-invariant
under the free group and whose free group is finite.

*Setting.* Use the module game (p6-prove Lemma M; equivalently R(t) of Lemma C): V_0 = C[P]^2 (row vectors over
C[G]), and a paid move adds C[P] w with w = h v, v in the current V, h a paid letter (free factors of the letter are
absorbed, since C[P] q h q' v = C[P] h (q'v)). A unit with M moves exists iff T e_1, T e_2 lie in V_M. Grade C[G] by
level (paid length): pi_j is the projection onto the span of the elements of level j. Levels are invariant under left
and right multiplication by P, so pi_j is a C[P]-bimodule map; paid letters have level 1; h times level 0 is level 1.
L denotes the composition length of C[P] as a left module.

*Proof.* Suppose M = 3 (M <= 2 is excluded by Lemma C). Write the moves as
  w1 = h1 x1,   w2 = h2 (x2 + y w1),   w3 = h3 (x3 + y' w1 + z' w2),   x_k in C[P]^2, y, y', z' in C[P],
and the final memberships as T e_l = a_l + b_l w1 + c_l w2 + d_l w3 (a_l in C[P]^2, b_l, c_l, d_l in C[P]).
1. *Level 2.* w1 has level 1 and a_l level 0, so T e_l = c_l pi_2(w2) + d_l pi_2(w3). Hence (C[P]T)^2, of length 2L,
   lies in C[P] pi_2(w2) + C[P] pi_2(w3), of length <= 2L; so the two are equal, pi_2(w2) = sigma T and
   pi_2(w3) = tau T with sigma, tau in C[P]^2, and with Sigma := [sigma; tau] and C' := [[c_1, d_1], [c_2, d_2]] we get
   C' Sigma = I: **C' is invertible.**
2. *The level-0 vectors x1, x2 form a basis.* pi_2(w2) = pi_2(h2 y h1) x1 (h2 x2 has level 1) and
   pi_2(w3) = Psi' x1 + pi_2(h3 z' h2) x2 with Psi' := pi_2(h3 y' h1) + pi_2(h3 z' h2 y h1) (h3 x3 has level 1). So
   [[pi_2(h2 y h1), 0], [Psi', pi_2(h3 z' h2)]] N = Sigma T with N := [x1; x2] in M_2(C[P]). The right side is
   invertible in M_2(C[G]), a finite-dimensional algebra, so **N is invertible.**
3. *Level 1.* The level-1 part of T e_l is zero: b_l w1 + c_l pi_1(w2) + d_l pi_1(w3) = 0. Since C' is invertible,
   pi_1(w2) and pi_1(w3) lie in C[P] w1 = C[P] h1 x1. In particular
     h2 x2 + pi_1(h2 y h1) x1 = kappa h1 x1   for some kappa in C[P].
4. *Contradiction.* Multiply on the right by N^-1 (allowed: right multiplication by a matrix over C[P]); x1 N^-1 =
   (1, 0) and x2 N^-1 = (0, 1). The left side becomes (pi_1(h2 y h1), h2) and the right side (kappa h1, 0). The second
   coordinates give h2 = 0, impossible. QED.
- The proof is a junk argument: input 2 must enter the moves through a level-0 vector (step 2), which leaves level-1
  junk h2 x2 that only another move could cancel, and with three moves there is none to spare (step 3).
- What it does not give: for m >= 3 junk can be carried down to level 0 by a later move and absorbed there, and that
  move can be packed with useful work; this is exactly the packing-minus-branching balance P - B of sl-combinatorics.
  So the argument does not extend directly. It does say where a two-array saver must spend its effort: on moves
  that carry junk down while carrying data up.
- Consistency: every frame-game search at m = 2 (section 4) found no 2 in 3, with or without scratch; Theorem H says
  none exists.
- Scope. The proof uses only: P finite (so C[P] has finite composition length; step 1), paid length invariant under
  P on both sides (so pi_j is a bimodule map), paid letters at level 1, and group elements forming a basis of C[G]
  (levels are supports). G may be infinite: M_2(C[G]) is stably finite (Kaplansky), and an element of M_2(C[P])
  invertible in M_2(C[G]) is invertible in M_2(C[P]) (apply the projection onto the span of P, which is C[P]-linear).
  So it covers the engine at block 2 (already known from Theorem R) and the widened game at block 2. It does not
  cover games whose free group is infinite or whose letters mix digit columns (XW), and it says nothing about a
  single fixed representation (the block-2 stand-in's setting), where group elements are not independent.
- Relation to Theorem T. Theorem T gives M >= m + 2 for two arrays in every frame game; that equals 2m only at m = 2.
  So in the frame-game laboratory the open window for two arrays is m + 2 <= M <= 2m - 1: from 2 in 5 at m = 3 on.
  The m = 3 searches of section 4 test exactly the smallest point of that window (2 in 5) and found nothing.
- Where m = 2 is special: at level m the producing block is square exactly when W = m = 2 (the two producers w2, w3
  against the two level-0 vectors x1, x2), which is what makes C' and N invertible. For m >= 3 with 2m - 1 moves
  the corresponding block is 2 x m against m x 2, and the level-1 junk of input 2 can be cancelled by later moves'
  level-1 parts (back steps), so the contradiction in step 4 disappears. A proof for m >= 3 needs a new idea.

## 3. Abstract frame games: a laboratory for Proposition 6 with W = 2

**Definition.** A frame game is (G, P, H, T): a group G (finite in all computations below; D_infinity in
Proposition E), a finite "free" subgroup P, paid letters H = a union of
double cosets P h P closed under inverses, a target T. Frames are the right cosets P g; a paid move takes P g to P h g;
the paid length l(g) is the distance from P to P g in this graph, and m := l(T). Contents are elements of C[G]; free
mixing is by M_W(C[P]); a unit is an identity in M_W(C[G]) as in section 2. The engine at block m maps onto the case
G = Sp(2m, 2), P = Siegel parabolic, H = P T_z P (sl-combinatorics Lemma 1.1(d); fp-build's frame game, which records
that every engine unit gives a frame-game program with the same W and M) [read].

The frame-level lower bounds on file use only frame-game facts: free steps fix frames, a paid move moves a term by
one edge, right multiplication is an isometry, and the target is at distance m (sl-algebra section 2a lists exactly
these for the widened game). So Corollary 2 (common-frame networks: M >= 2m for W = 2), Corollary 2' (single-term
contents), Theorem 4.3 (no frame packing in some gauge, labels modulo the core of P), Theorem 5, the ladder's
Theorem D (M >= W + m - 1) and mw-big56's Theorem T (M >= W + m) hold in every frame game (the last two by their own
scope remarks; I have not re-checked them line by line). The question this section tests: **is "two arrays need 2m" true in every
frame game?** If some small frame game had a two-array saver, the W = 2 case of Proposition 6 would not follow from
the frame-game axioms, and any proof would have to use Clifford-specific structure (as Theorem R does). If none does,
the abstract statement is a candidate theorem that p6-prove could attack without Clifford specifics.

**Lemma D (characters trivial on P force a congruence)** [proved]. If psi is a one-dimensional character of G with
psi(P) = 1, then every unit with one paid double coset has M = W m modulo the order of psi(h).
*Proof.* Apply psi and the trivial character epsilon entrywise. They agree on M_W(C[P]); so with f_k := psi(F_k) =
epsilon(F_k) in GL_W(C), the epsilon-image gives prod det f_k = 1 and the psi-image gives prod det f_k * psi(h)^M =
psi(T)^W = psi(h)^(mW) (T lies in (P h P)^m P). QED.
- Examples: dihedral groups D_n with n even (psi(r) = -1, psi(s) = 1) and the hyperoctahedral groups B_m with
  P = S_m (psi = product of the signs, psi(flip) = -1): there M must have the parity of W m, so 2 in 2m - 1 is
  impossible and 2m - 2 is the first open case.
- In the engine no such psi exists (the frame graph has triangles: three paid letters with product in P force
  psi(h)^3 = 1, and psi(h)^2 = 1 since (dir z)^2 = X^z is free), consistent with the odd deficits on file.

**Proposition E (the infinite dihedral game: no saving at all)** [proved]. G = D_infinity = <r, s>, P = <s>, H = PrP,
T = r^m. Every unit, any W, has M >= W m.
*Proof.* In the two-dimensional representations rho_z (z on the unit circle), in the eigenbasis of r, rho(r) =
diag(z, 1/z) and rho(s) swaps the basis; rho(C[P]) consists of constant matrices. A unit's rho_z-image is a Laurent
polynomial matrix in z of size 2W; free factors are constant in z, and each Y = diag(z, 1/z, 1, ..., 1) has every
k x k minor of z-degree at most 1. By Cauchy-Binet, the largest z-degree among the W x W minors of a product is at most
the sum over the factors, so it is at most M after M paid moves. The target diag(z^m, z^-m) (x) I_W has the W x W
minor z^(mW) (rows and columns of the z^m entries). So M >= mW. QED. (The frame graph is a path; the argument is the degree valuation, which a finite group cannot have.)

**Lemma F (D_4, by hand)** [proved]. In D_4 (frame graph a 4-cycle, T = r^2, m = 2), two arrays need 4 moves: M = 3
is excluded by Lemma D, and M = 2 because both moves must then act on one array (the first-move/last-move argument
forces F_1 to be monomial), so the only top-frame term of the result is a T (central) times a scalar a, and
pi_T(Pi_11) = a F2_11 F0_11 T must be a unit times T while pi_T(Pi_12) = a F2_11 F0_12 T = 0 forces F0_12 = 0 and then
pi_T(Pi_22) = a F2_21 F0_12 T = 0, a contradiction.

**Computations** [computed: `toylab.py`, `sweep.py`, `sweep_scr.py`, `t_*.py`; outputs `out_*.txt`, `sweep*/`].
- Groups from permutation or F_2-matrix generators; all irreducible representations numerically, from the eigenspaces
  of a random Hermitian element of the right-regular action (sum of squared dimensions checked to equal |G|, and the
  homomorphism property checked). A unit is checked irrep by irrep, i.e. as a formal identity in C[G].
- Solver: unknowns F_0..F_M in M_W(C[P]) (W^2 |P| complex numbers each); the identity is multilinear, each F_k enters
  linearly; Levenberg-Marquardt on all F_k jointly from random complex starts (with a stagnation stop), after an
  alternating-least-squares variant proved too slow. Moves wlog on slot 0 (free swaps). With scratch (`SolverS`):
  W_t arrays, 2 data arrays, the scratch rows unconstrained at the end, copies and overwrites allowed.
- What a negative line means: the best relative residual over the starts stayed far from 0 (0.5 to 0.7). It is not
  a proof of impossibility; savers are expected in special position (MAP section 1.3), which random starts can miss.
  The calibrations show the solver finds the plain 2 in 2m from random starts in the same games (most often in the
  first start), so the negatives are not a failure to converge at that size.
- Results: the table in NOTE.md ("Results of the frame-game searches"), with every line traceable to RUNLOG.txt and
  to the output files: `sweepjob/remote-out/outA.txt` and `outB2.txt` (no scratch, random starts),
  `scrjob/remote-out/outS12.txt` (one scratch array), `job3/out3ab.txt` (structured starts), `job4/remote-out/
  out4ab.txt` and `job5/` (module game), `out_sp4.txt`, `out_bld.txt`, `out_obstr_J42.txt`, and the local `t_*.py`
  runs quoted in RUNLOG.txt.
- Layer relaxations (`modgame.py`: the module game in the group basis, scipy least squares, only chosen levels of the
  final identity enforced; calibrated: D_5 2 in 4 and D_7 2 in 6 found at 1e-15, D_7 2 in 5 with all levels 0.50,
  equal to the representation-based value). D_7, 2 in 5: levels {3, 2} only: 1.3e-6 (decreasing; exact or border
  solution of the relaxation); all levels: 0.50. So at m = 3 the top two layers alone admit 2 in m + 2 and the third
  layer obstructs it. S_4 with P of order 3 (class 0), 2 in 5, levels {3, 2} only: exact (2e-15; `t_layer.py`: the
  solution's level-2 and level-3 errors are 2e-15, its level-1 error has norm 15, level 0 is absorbed), all levels:
  0.41. (S_4 with P of order 2 is not comparable: `modgame.py` uses one letter of its inverse pair only.) At m = 2
  the same holds one level up: 2 in 3 with only level 2 enforced is solvable (J(4,2) exact, 5e-16; D_5 1.4e-6),
  and Theorem H's level-1 junk step is exactly what kills it. Pattern so far: with 2m - 1 moves the upper levels
  alone can be met, and the bottom level (leftover one-step terms) is the obstruction. Further games in RUNLOG.txt.
- Single-representation relaxations (`t_obstr.py`, J(4,2), 2 in 3): every irreducible representation alone admits
  the identity or comes arbitrarily close to it (error creeping to 1e-3 or below: border behaviour), every pair of
  irreducibles already obstructs it (0.34 to 0.41), all five together give 0.61. The obstruction is a coupling of
  representations through the shared free coefficients, which is how Theorem H's argument (supports of group
  elements) sees it.

## 4. Session 2 (15:21 to 17:40 EDT): the single-term class as an exact game

The 13:20 lead of the brief: two inputs, any scratch, every MOVED register holding at most one term per input
(p6-prove's §4; Conjecture T open from block 4). Labels as above.

### 4.1 The game

By p6-prove's Lemma M and Lemma G (`../p6-prove/PROOF.md` §1, §4) [read], the class is the following game on the
frames (positions) of any frame game (G, P, paid class), which `stg.py` implements.
- State: the available positions A_1, A_2 of each input, and the locked pair copies. A pair copy is a vector
  (c_1 g_1, c_2 g_2); it sits at input-1 position x = P g_1 and input-2 position y = P g_2 = x·t, with twist t = g_1^-1 g_2
  acting on positions by right multiplication. Copies of one lineage share t and the ratio c_2/c_1.
- Moves (cost 1): a single from an available position to a neighbour; a copy of a lineage to a neighbour x' of x
  (its input-2 end goes to x'·t); a CREATED pair (x_0 in A_1, x_0·t in A_2, any t; new lineage); a COMPOSITE (the
  vector of V supported on (1, u), (2, w) given by a path of at least three copies in a locked component; twist the
  alternating product along the path; new lineage).
- Closure: in the bipartite position graph (one edge per copy), a component is unlocked if it contains an available
  vertex or a cycle; then all its vertices are available. Lemma G's exact rule asks for a cycle whose lineage ratios
  do not cancel; "any cycle" over-approximates it, which is the safe direction for lower bounds and for (I) below.
- Checked at every state: (I) max level(A_1) + max level(A_2) <= M (Conjecture T), and the goal (both inputs available
  at a frame of the target's free orbit).

**Lemma 4.1 (symmetries)** [proved]. For free elements a, b in P, the map (A_1, A_2, copies (t, x)) -> (A_1 a, A_2 b,
(a^-1 t b, x a)) and the input swap (A_1 <-> A_2, (t, x) -> (t^-1, x t)) map plays to plays and preserve levels and the
target orbit. *Proof.* Right multiplication by a in P fixes the raw frame, is an automorphism of the frame graph and
commutes with the left (paid and free) actions; it is the effect of applying the free letter a to input 1 before
the first move (entries of contents are right-multiplied by a). Likewise b for input 2; the swap relabels inputs. ∎
States are deduplicated by the lexicographically least image under these maps (sound for any subset of them).

**Lemma 4.2 (pruned moves are dominated)** [proved]. A single moved onto an available position adds nothing; a
created pair moved so that one of its ends lands on an available position gives exactly the state of a single move
of the other input (catalyst release), and with both ends available it adds nothing. So omitting these moves loses
no state up to a wasted move, which cannot help (I) or the goal.

### 4.2 Calibration against the algebra

`t_stgcal.py` [computed]: random single-term plays in D_7, S_4/Z_3, J(4,2), D_5 (40 plays of 8 to 10 moves each),
every generator an explicit vector of C[G]^2 with random complex coefficients, V the span of all P-translates. At
every step the true availability (is e_(l,g) in V?) was compared with stg's closure: 1,440 checks, 0 mismatches. The
plays contained 20 releases by a cycle and about 300 releases by a catalyst; every composite was found by linear
algebra (a vector of V supported on the two named vertices) and had the twist stg predicts (11 of 11). Control:
the closure without the cycle rule misses exactly the 20 cycle releases, so the check discriminates. A second run biased
towards composites (B_3, S_4/Z_2, D_9, S_4/Z_3; 1,200 checks): the module never had an available single that the
solver missed. In D_9 the solver four times declared a single available that the module lacks: a cycle that the
permissive rule unlocks but whose ratios cancel (degenerate in Lemma G's sense). So the solver can over-predict,
never under-predict; over-prediction is the safe direction for every negative result here.

### 4.3 Sub-games of the engine

**Lemma 4.3 (cube and Hamming games are sub-games of the block-m frame game)** [proved].
- (Hamming) Fix a splitting of F_2^(2m) into m hyperbolic planes, adapted to L0 = {(0|b)} and L_T = {(z|z)}. The
  frames that are direct sums of lines of the planes form Γ_1^m = H(m, 3) (3^m frames), isometrically: the
  intersection dimension adds over the planes, so d is the Hamming distance. Sp(2,2)^m = S_3^m lies in Sp(2m,2) and
  acts coordinatewise; its elements at distance 1 from P (one coordinate moved) lie in the paid double coset (the
  double cosets of P in Sp(2m,2) are indexed by d(L0, g L0)); its intersection with P is Z_2^m.
- (Cube) Inside it, the apartment of the frames spanned by e_i or f_i (2^m frames, the m-cube, isometric) with the
  Weyl group B_m of monomial symplectic matrices; B_m ∩ P = S_m; a swap e_i <-> f_i is paid. Its level-m frame
  is the span of the e_i, not L_T; both lie in the distance-m double coset of P, so a program for it becomes one for
  the transform by free letters at the two ends.
- So every play of the single-term game in H(m,3) with twists S_3^m, or in the m-cube with twists B_m, is a play of
  the single-term frame game at block m with the same moves, availability and levels (the closure rule is the same
  position-level rule). Hence a violation of (I) there at depth t with excess e = L_1 + L_2 - t >= 1, finished
  plainly (each coordinate still at 0 changed once; P transitive on the level-m frames of the sub-game), is a
  single-term two-array frame-game program with 2m - e moves at block m. Non-existence in a sub-game proves nothing
  for the engine.
- Correction (also stated in NOTE.md): the character argument (Lemma D) is a
  statement about square matrix words without scratch; with scratch no parity is proved in cube games, so the cube
  checks are genuine tests of single-term 2 in 2m - 1, not calibrations. In the Hamming games no character of S_3^m
  is trivial on Z_2^m except the trivial one, so Lemma D gives nothing there either.

### 4.4 Results [computed; outputs in job10/ to job19/ (remote-out) and RUNLOG.txt]

(I) held at every reachable state in every case below; wherever 2m was within the completed depth, the goal was
first reached at exactly 2m. "Complete to M"
means every state up to M moves (canonical deduplication up to M - 1, the last layer checked successor by successor).

| game | frames | m | complete to M | largest layer |
|---|---|---|---|---|
| D_5 (2 classes), J(4,2), J(5,2) (2 classes), S_4/Z_4, GL(3,2) Borel (2 classes), A_5/S_3 (2 classes) | 5 to 21 | 2 | 4 = 2m | 82,000 states |
| D_7 (3 classes) | 7 | 3 | 6 = 2m | 451 states |
| S_4/Z_3, S_4/Z_2 | 8, 12 | 3 | 6 | 187,203 states |
| A_5/Z_2 (class 2), A_5/V_4 (3 classes) | 30, 15 | 3 | 6 | 475,636 states |
| A_5/Z_3 (class 0), A_5/Z_5 (class 0) | 20, 12 | 3 | 5 | state cap at 6 |
| B_3 (cube: block-3 apartment sub-game) | 8 | 3 | 6 | 63,080 states |
| H(3,3) (block-3 Hamming sub-game) | 27 | 3 | 6 = 2m | 79,829,304 successors |
| D_9 (3 classes) | 9 | 4 | 8 = 2m | 2,536 states |
| B_4 (cube: block-4 apartment sub-game) | 16 | 4 | 6, and the goal at 7 | 91,211,451 successors (job14); 93,790,253 (job16) |
| H(4,3) (block-4 Hamming sub-game) | 81 | 4 | 5 | 111,232,729 successors (2 parts) |
| A_5/Z_2 (classes 1, 3), A_5/Z_3 (classes 1, 2) | 30, 20 | 5 | 6, 7 | state cap |

Reading: no toy frame game, and neither engine sub-game, has a single-term state with L_1 + L_2 > M at these depths.
In particular no single-term 2 in 2m - 1 exists in the toys at m = 2, 3, 4 (D_9). In H(4,3) a single-term 2 in 7
would need its first violation of (I) at depth 6 or 7.

**Proposition 4.4 (the block-4 cube has no single-term 2 in 7)** [computed, exhaustive, `job16/`]. In the cube
sub-game B_4 (the 16 frames of one apartment of the block-4 frame graph, twists the Weyl group B_4, free group S_4),
no single-term two-input play delivers both inputs in 7 moves. With Theorem R (no 2 in 6 at block 4), the cheapest
single-term two-array program there costs exactly 8 = 2m.
- Method (`stg.py` bfs3): canonical deduplication to depth 5 (six parts, split at layer 4), every successor at
  depth 6 checked for (I) (93,790,253 in all: it holds), and for every depth-6 state whose highest input-1 and
  input-2 frames (available or locked) both reach level 3, an exact test whether one more move delivers both inputs
  (1,161,509 states tested: none). The level filter is admissible because a move creates frames only next to
  existing ones; the test skips only moves whose new generator neither ends on a level-4 frame nor touches a
  component holding one, which cannot release a level-4 frame.
- Check of the test itself: with and without the move pruning it gives the same goal counts at 2m and none at
  2m - 1 on D_7, S_4/Z_3, B_3, J(4,2), S_4/Z_2 (the first version, which identified the new generator by comparing
  sets of (twist, frame) pairs, missed goals when two lineages share a twist; fixed before job16 ran).
- Since B_4 is a sub-game (Lemma 4.3), this excludes nothing in the full engine: it says the apartment with its Weyl
  twists is too small to carry a single-term saver at block 4.

**Proposition 4.5 (no single-term 2 in 5 in the m = 3 games)** [computed, exhaustive, `job19/`, `t_stg_2m1.py`].
For every paid class with m = 3 of D_7 (3 classes), S_4/Z_2, S_4/Z_3, B_3, A_5 with free subgroups Z_2, Z_3, Z_5
(2 classes), V_4 (3 classes), S_5 with D_5 and with S_3, A_6 with 3^2:2 (2 classes), and J(6,3) (S_6 with S_3 x S_3,
classes 0 and 1; 19 cases in all), there is no single-term two-input 2 in 5: dedup to 3 moves, every 4-move state
checked for (I) (it holds), and the exact one-move goal test on every 4-move state that passes the level filter
(from 89 to 121,051 states per case). 2 in 5 = 2 in m + 2 is the first size the two-layer bound allows at m = 3,
so in the single-term class these games add nothing beyond the bound's sharpness question: they all need 2m.

### 4.5 Why the natural crossing twist ties [argued]

Let t swap L0 and L_T (for example the label of the transform itself). A t-pair at (z, z·t) carries input 2 towards
L0 when input 1 goes towards L_T, and the other way round. Branching it delivers both inputs for m moves in all:
moved to (L_T, L0) it releases input 1 at L_T against the free single of input 2 at L0, and moved to (L0, L_T) it
releases input 2 at L_T. But a t-pair can only be created from available singles at (x_0, x_0·t), whose levels add
up to lv(x_0) + d(x_0, L_T) >= m, or composed from earlier t-copies; so (I) at the creation time makes the total
at least 2m. Twist-1 pairs do not help to create it: a twist-1 copy at (z, z) is locked until a second object
reaches z. This is the pattern seen in every search: potential is cheap to build (a pair moved k times holds levels
(k, k) after k moves), and unlocking it costs exactly what it would have saved.

### 4.6 Candidate invariants for a proof of Conjecture T (data, and one small lemma)

`t_inv.py` [computed] evaluates candidate invariants on EVERY reachable state of the exhaustive runs (local: D_7 to
6 moves, S_4/Z_3, B_3, S_4/Z_2 to 5, J(4,2), H(2,3) to 4, D_9 to 7; Colab, `job15/`, `job17/`: S_4/Z_2, B_3, S_4/Z_3
to 6, D_9 to 8, D_7 (class 1) to 6, A_5/Z_2 (class 2) to 5, H(3,3), J(5,2), GL(3,2) Borel to 4). Write W(u) for the set of objects that witness an
input-1 vertex u (the mark A if u is available, and every lineage with a locked copy whose input-1 end is u), and
likewise W(w) for input-2 vertices. Maximum excess (left side minus M) over all states:
- C0: L_1 + L_2 <= M: excess 0 (Conjecture T).
- C1 and its component form: L_1 + (highest input-2 level of a locked component) <= M, and symmetrically: 0.
- C2: lv(x_E) + lv(y_F) <= M for locked copies E, F of different lineages: 0.
- C2c: two vertices of one locked component at path distance >= 3 (whatever lineages the end copies have): 0.
- S: every pair (u, w) witnessed by two different objects (|W(u) ∪ W(w)| >= 2, or both available): 0.
- C2s: the same for two copies of ONE lineage in different places: excess up to +2 (expected: branches of one
  lineage share their prefix).

So the data support **Conjecture S**: at every moment, lv(u) + lv(w) <= M for every input-1 vertex u and input-2
vertex w, except when both are witnessed only by one and the same lineage and are not joined inside one locked
component by a path of three or more copies. Conjecture S contains Conjecture T (u, w available) and would be the
natural induction hypothesis: a release of a component C by a single at one of its vertices is covered by S
applied to the catalyst's starting point and the other end (checked by hand for the cases "C is one copy", "the
two ends share a vertex with another lineage" and "path distance >= 3"). The one case I could not close by hand is
a release that makes a same-lineage pair (x_E, y_F) with E inside C and F outside C a two-object pair; the next
lemma bounds its excess but does not show it is paid for.

**Lemma 4.7 (the excess of one lineage is its shared prefix)** [proved, given S at earlier times]. Let the root copy
of a lineage be produced at time t_R from a created pair at available (s_1, s_2) or a composite at (u_0, w_0) (path
>= 3 in a locked component), so lv(s_1) + lv(s_2) <= t_R - 1 by S at time t_R - 1. For copies E, F of the lineage
with depths δ(E), δ(F) in its move tree (the root has depth 1): lv(x_E) + lv(y_F) <= M + δ(LCA(E, F)).
*Proof.* Each move of the lineage changes each end's level by at most 1, so lv(x_E) <= lv(s_1) + δ(E) and
lv(y_F) <= lv(s_2) + δ(F); the moves on the two root paths are δ(E) + δ(F) - δ(LCA) distinct moves after time
t_R - 1. ∎
- For E = F this says a copy's excess is at most its depth: the potential of a pair moved k times is k, which is
  exactly the cost of a second object reaching it in the plain way. A proof of Conjecture T along these lines
  must show that every event turning a same-lineage pair into a two-object pair costs at least the shared prefix.

### 4.7 A dimension potential (Conjecture Φ)

For a module V (Lemma M) and k >= 1 let G_k(V) = {v in V : every term of v has level >= k} (a C[P]-submodule) and
Φ(V) = Σ_{k=1..m} dim ε(G_k(V)), ε the augmentation of each entry (ε: C[G]^W -> C^W).
- Φ(V_0) = 0 (V_0 = C[P]^W lies at level 0); at the end K e_l is in G_k for every k <= m, so Φ >= W m. An available
  single of input l at level L_l puts e_l into ε(G_k) for k <= L_l, so Φ >= L_1 + L_2 when W = 2. Hence
  "Φ(V_M) <= M" implies Conjecture T, and "every move raises Φ by at most 1" implies M >= W m for every unit.
- So Φ cannot be a potential for all units: frame-game savers exist at large blocks (Lemma N images of the engine's
  savers), and for them some move raises Φ by 2 or more. Theorem C (p6-prove) is the homogeneous case, where the
  chain E_k ⊆ E_(k-1) shows each move raises Φ by at most 1.
- **Dual description** [proved, linear algebra]: dim ε(G_k(V)) = W - dim Λ_k with Λ_k = {λ in C^W : some functional
  vanishing on V takes the value λ_l on every coordinate (l, g) of level >= k}. For single-term V this is a statement
  about the position graph: a functional vanishing on V is zero on unlocked components and, on a locked tree
  component, is one free multiple of a fixed potential ψ (ψ(2, y) = -ψ(1, x)/r along a copy with ratio r). So level k
  carries the constraint λ_l = 0 for every available position of input l at level >= k, and for every locked
  component: λ_1 = 0 if it has two input-1 vertices of level >= k with different ψ (and likewise λ_2); λ in the line
  λ_1/λ_2 = ψ(u)/ψ(w) if it has high vertices u, w of both inputs. Φ = Σ_k rank of these constraints. Λ_k grows
  with k, so Φ = κ_1 + κ_2 with κ_j the highest level carrying constraints of rank j.
- **Data** [computed]: `t_phi.py` (module Φ, explicit vectors): random single-term plays in D_7, S_4/Z_3, J(4,2),
  B_3, Sp(4,2) (the block-2 engine frame game): 1,140 moves; random multi-term moves (D_7, S_4/Z_3, J(4,2)): 691;
  "packed" multi-term moves (random elements of G_k(V), moved; D_7, S_4/Z_3, J(4,2), B_3): 941. The largest
  increase of Φ in one move was 1 in every case. `t_inv.py` (position-level Φ, an upper bound for the module Φ
  because it treats every lineage ratio as independent): Φ <= M at every reachable state of the exhaustive runs
  (D_7 to 6 moves, S_4/Z_3, B_3, S_4/Z_2 to 6, D_9 to 8, A_5/Z_2 to 5, H(3,3), J(4,2), J(5,2), GL(3,2) to 4), with
  equality attained. `t_phi.py` per-move data from `job15/`: 1,680 more single-term moves (H(3,3), S_4/Z_2, D_9, B_3,
  A_5/Z_2) and 1,813 more packed multi-term moves (the same games and Sp(4,2)), largest increase 1.
  `job18/`: the block-4 sub-games H(4,3) and B_4 with explicit vectors (2,592 and 768 coordinates): 700 single-term,
  1,100 multi-term moves, largest increase 1.
- **Conjecture Φ.** For single-term plays, Φ(V_M) <= M. It implies Conjecture T. A one-move proof cannot work as
  stated: joining two locked components that each hold a high vertex of only one input creates a line constraint
  at all levels up to the lower of the two at once. In the plays I could construct, the crossing twists needed for
  such components were paid for earlier by singles that already fill those levels (creating a pair with twist t
  needs available ends whose levels add up to at least d(L0, L0 t^-1)), so the proof has to carry that cost along,
  for example through Conjecture S and Lemma 4.7.
