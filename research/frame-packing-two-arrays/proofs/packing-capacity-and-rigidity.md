<!-- header: frame-packing-two-arrays -->
# How much one paid move can carry (the capacity track of packing round 1)

**About this file.** This is the working proof file as an AI agent wrote it during the round of 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent. Its opening lines say that no skeptic had read it; that was true when it was written. What happened afterwards is below, and where a later skeptic's verdict differs from this file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean. File names inside it refer to the programs and outputs now in `../programs/` and `../data/`; [README.md](README.md) maps the working names.

**Status now.** Upheld by one AI skeptic with its own code, with four corrections, none load-bearing: (1) counting packing by image frames is not always the sharpest count (an exact block-2 unit is counted more sharply by labels); (2) one parenthetical of Corollary A'3 holds only in raw coordinates; (3) the reason given for step 1 of Proposition D is wrong inside a stretch where both arrays hold the same input, but the step holds (every step is an invertible 2 x 2 matrix); (4) the skeptic said the two-array table left out 2 in 7 at block 6 (worth 0.301). That fourth correction is itself mistaken: 2 in 7 at block 6 is already excluded by the two-layer bound (every unit with two or more arrays needs at least W + m moves, upheld on 14 October), which the skeptic did not take into account, so the original table was right on that point. The step "every ceiling is a refund ceiling" is argued, not proved.

---

# PROOF.md: packing capacity on the frame graph (fp-capacity)

Packing round, 15 October 2026. One agent; no skeptic has read this. Labels: [proved] (proof written here),
[computed] (a program in this folder and its output), [cited] (a published theorem, cited from memory), [read]
(taken from one of our notes, not re-derived). Programs: `vflow.py` (Theorem A on actual units, exact), `geom.py`
(Lemma S, joint ascent, rigid pairs), `vflowB.py` (Theorem A′), `tables.py`; outputs `vflow_s2.out`, `vflow_s1.json`
(first run, older code), `vflow_s2.json`, `vflowB_s3.out`, `vflowB_s4.out` (and .json), `geom.json`, `geom.out`,
`tables.out`. The engine (exact contents, frame labels) is sl-combinatorics' `ledger.py`,
imported read-only.

Everything is in the engine game at block m, with rank-1 moves (a grouped rank-r move counts as r consecutive rank-1
moves on one array, as in sl-combinatorics Remark 4.3b). The setting, letters, contents, gauges, labels, positions
and levels are those of `../../slack/sl-combinatorics/PROOF.md` sections 0 and 1, which I use without change:
- 𝒩 = Paulis·scalars (normal in G, inside B); frame labels are elements of G/𝒩 = Sp(2m, 2);
- the *position* of a label λ is the coset Bλ, i.e. the Lagrangian X = λ̄⁻¹L0, a vertex of Γ_m = DSp(2m, 2);
- the *image frame* of λ is F = λ̄L0; level δ(λ) = d(F, L0) = d(X, L0);
- a gauge h = (h_1, …, h_W) ∈ G^W; the label of a term (g, l) is g h_l mod 𝒩; S_a = position of h_a, T_a = position
  of K h_a, d(S_a, T_a) = m. The *raw gauge* is h = 1: then S_a = L0 and T_a = L1 := K̄L0 for every input.
- Δ = Wm − M is the saving in moves.

## 1. Frame packing, exactly

**Definition 1.1 (frame classes, lineage, liveness).** Fix a gauge h. At time t the *frame classes* of register r
are the distinct labels Λ_r(t) of its terms; the part of c_r with class κ is c_{r,κ}. The *class graph* has a node
(t, r, κ) for every class present, and an edge to time t+1 when the step carries it:
- a move or free letter x on r: (t, r, κ) → (t+1, r, x̄κ) (left multiplication; a bijection on classes);
- a mix c_r := Σ_s μ_s b_s c_s (μ_s ≠ 0): (t, s, κ) → (t+1, r, b̄_s κ) whenever that class is present at t+1;
- a gate (scalar matrix A on registers R): (t, s, κ) → (t+1, r, κ) for A_rs ≠ 0, when present;
- untouched registers: (t, r, κ) → (t+1, r, κ).
A node is *live* if a path leads from it to a final class of a delivering register. (`vflow.py` computes exactly
this, by a backward pass; it can only over-count liveness, which weakens nothing below.)

**Definition 1.2 (virtual moves and packing measures).** For the j-th paid move (on register r_j at time t_j) let
k_j(h) be the number of live classes of r_j just before it. The pairs (j, κ) with κ live are the *virtual moves*;
each moves one class from its position X to the position X' of x̄κ, and d(X, X') = 1 (sl-combinatorics Lemma 1.3).
- gross frame packing P_fr(h) = Σ_j (|Λ_{r_j}| − 1)^+ (sl-combinatorics' P with frame labels);
- live packing P_live(h) = Σ_j (k_j(h) − 1)^+;
- **irreducible packing** P* = min over all Sp-gauges h of P_live(h). (For W = 2, P_live depends only on the
  relative gauge h_2 h_1⁻¹, and `vflow.py`'s `pstar_W2` computes P* exactly: every relative gauge outside the finite
  set {ḡ_2⁻¹ ḡ_1 : terms (g_1, 1), (g_2, 2) of one register at one time} merges nothing, like a generic one.)
- raw-gauge weights: f(X) = d(L0, X) − d(L1, X) ∈ [−m, m]; a virtual move X → X' has weight w = max(0, f(X') −
  f(X))/2 ∈ {0, ½, 1}. Weight 1 (a *progress* step: away from the raw frame and toward the transformed frame) is
  exactly an upward step of the interval between L0 and L1 or of a parallel layer; weight ½ is a step that changes only
  one of the two distances. P_w := Σ_j (Σ_{live κ} w(j, κ) − 1). With the level potential instead,
  P_asc := Σ_j (number of live classes the move ascends − 1).

P_live ≤ P_fr always. A *coherent* unit (sl-combinatorics Definition 3.2) has P_fr(h) = 0 in its coherent gauge, so
P* = 0. The order-3 control words have P* = m (computed below).

## 2. The virtual network and its cut lemma

**Lemma 2.1 (virtual contents)** [proved]. Choose a representative λ_κ ∈ G for every class. Then
c_{r,κ} = λ_κ Σ_l A_{κ,l} h_l⁻¹ y_l with A_{κ,l} in the algebra 𝒜 = C[𝒩] spanned by Paulis and scalar elements
(𝒜 is semisimple, of finite composition length n as a left module over itself). Write v_{r,κ} = (A_{κ,l})_l ∈ 𝒜^W.
Changing λ_κ multiplies v by a unit of 𝒜 on the left, so the cyclic module 𝒜v_{r,κ} is well defined.
- (i) A move or free letter x on r maps class κ to x̄κ, with representative xλ_κ, and leaves v unchanged.
  Every class of the register moves, each along its own edge of Γ_m (free letters: no move, B x λ = Bλ).
- (ii) In a mix, the class κ' of the result has v_{r,κ'} = Σ μ_s n_s v_{s,κ} over the summand classes with
  b̄_s κ = κ' (n_s ∈ 𝒩, from b_s λ_κ = λ_{κ'} n_s). All those summand classes have the result's position, since
  B b_s λ_κ = B λ_κ. Gates: the same with b_s = 1.
- (iii) Input a starts as one class at S_a with v = e_a. A delivering register r ends as one class at T_{π(r)} with
  v = u e_{π(r)}, u a unit scalar.
*Proof.* g h_l ∈ λ_κ 𝒩 gives g = λ_κ n h_l⁻¹ with n ∈ 𝒩, so grouping the terms of class κ gives the form. (i) x(λ_κ A
h⁻¹ y) = (xλ_κ)(A h⁻¹ y). (ii) b_s λ_κ and λ_{κ'} lie in one 𝒩-coset, and 𝒩 is normal. (iii) is the definition of a
unit. ∎

So the live virtual sub-registers form an 𝒜-linear network on Γ_m: data enter at S_a, are combined only at a common
position, and travel only by virtual moves. A paid move pays for one virtual move of every live class of its array.

**Lemma 2.2 (directed cut lemma with packing)** [proved]. For every set Q of positions,

  #{virtual live moves entering Q} ≥ #{a : S_a ∉ Q, T_a ∈ Q},  and  #{leaving Q} ≥ #{a : S_a ∈ Q, T_a ∉ Q}.

*Proof.* Let V_Q(t) ⊆ 𝒜^W be the left submodule generated by v_{r,κ} over live classes positioned in Q at time t.
By Lemma 2.1, V_Q(t+1) ⊆ V_Q(t) + Σ (𝒜 v over the live virtual moves entering Q at step t): a live class in Q at
t+1 comes from a live class in Q (free letters, untouched registers, and mixes, whose summands sit at the result's
position and are live because they have an edge to a live node), or from a live class outside Q by an entering move.
Let π project 𝒜^W onto the coordinates a with S_a ∉ Q. Then π V_Q(0) = 0, π V_Q(end) ⊇ ⊕_{S_a ∉ Q, T_a ∈ Q} 𝒜 e_a
(Lemma 2.1(iii); final delivering classes are live), of length n·#{…}, and each entering virtual move adds at most n
to the length. Leaving: the same for the complement. ∎

This is sl-combinatorics' Lemma 4.2 with one change: a paid move on an array holding k live frame classes is
counted k times, once per class, and the hypothesis "no frame packing" is dropped. The directed form is new.

## 3. Theorem A: each saved move needs one extra live class carried by a paid move

**Theorem A (packing capacity)** [proved; computed on 152 (unit, gauge) rows, `vflow_s2.out`].

(a) *Raw gauge, every unit, any W, any scratch.* For every function φ on positions with φ(L1) − φ(L0) = m,

  Σ over virtual live moves of max(0, φ(X') − φ(X)) ≥ W m.

With φ = level: Δ ≤ P_asc = Σ_j (live classes the move ascends − 1). With φ = (f + m)/2: Δ ≤ P_w ≤ P_live(raw).

(b) *Two arrays, every Sp-gauge h* (any block, any scratch): M + P_live(h) ≥ 2m. Hence **Δ ≤ P*.**

(c) *Every gauge whose demand layout has one source position, or one sink position, for all inputs* (any W):
Δ ≤ P_live(h). The same holds for every layout for which the cut condition implies a fractional multicommodity flow
(Papernov–Lomonosov layouts [cited, not re-read]).

(d) *Every gauge, any W:* Δ ≤ P_live(h) + γ_m(h), where γ_m(h) = Wm − (least number of edge uses of an 𝒜-linear
network code on Γ_m for the unicast demands S_a → T_a) is the coding gain of that layout. γ_m(h) = 0 in the cases of
(a)–(c), and γ_m(h) ≤ W(m − L_m) for the proved floor L_m (Theorem R and descendants).

*Proof.* (a) In the raw gauge every S_a = L0 and T_a = L1. Give each arc X → X' of Γ_m the capacity c(X → X') =
the number of virtual live moves along it. For every Q with L1 ∈ Q ∌ L0, Lemma 2.2 gives at least W units of capacity
entering Q. By directed max-flow min-cut there is a flow F ≤ c of value W from L0 to L1. Decompose F into paths and
cycles; along any path from L0 to L1 the increments of φ sum to m, so the positive parts sum to at least m, and cycles
contribute ≥ 0. Hence Σ max(0, Δφ)·c ≥ Σ max(0, Δφ)·F ≥ Wm. Each paid move contributes Σ_{live κ} max(0, Δφ) ≤ k_j
(|Δφ| ≤ 1 on an edge for the two potentials used, since each distance changes by at most 1), so
Wm ≤ Σ_j Σ_κ max(0, Δφ) = M + Σ_j (Σ_κ max(0, Δφ) − 1).
(b) Undirected capacities c(e) = virtual live moves along e in either direction. Lemma 2.2 is the cut condition for
the two demands (S_1, T_1), (S_2, T_2). Hu's two-commodity theorem (T. C. Hu, Oper. Res. 11, 1963) [cited] gives a
fractional flow f_1 + f_2 ≤ c of value 1 each; each unit travels at least d(S_a, T_a) = m. So Σ_e c(e) ≥ 2m, and
Σ_e c(e) = Σ_j k_j ≤ M + P_live(h). (W = 1: max-flow min-cut.)
(c) One source (or one sink): a super-sink (super-source) reduces it to directed single-commodity max-flow, and
Lemma 2.2 gives the cut condition; then as in (b).
(d) The live virtual network is an 𝒜-linear network code on Γ_m with Σ_j k_j ≤ M + P_live(h) edge uses
(Lemma 2.1). Conversely every 𝒜-linear code is a frame-coherent unit (an 𝒜-combination is a mix with Pauli letters;
an edge from a class is reachable by a free letter followed by one paid letter), so its edge uses obey every proved
floor: γ_m(h) ≤ W(m − L_m). ∎

**Remarks.**
- (a) is "every saving is packing" as an *inequality on live packing*, valid in raw coordinates for every unit. It
  never mentions branching. sl-combinatorics Corollary 2.3 (Δ = P − B − 2D − Z − I0) is the exact ledger with gross
  packing and branching; neither statement implies the other.
- (b) is the quantitative form of sl-combinatorics Theorem 4.3 (which is the case P_live(h) = 0) and of
  sl-algebra Corollary 2 (common-frame networks). **A two-array unit saving k moves needs at least k units of live
  frame packing in every gauge.**
- The bound in (a) can be loose by capacity that no flow uses: surplus classes (the order-3 controls carry a whole
  unused path), and [argued] bundles of three or more classes that step back and forth, which can raise P_w by one
  per round trip without moving anything. The sharp form keeps only the capacity used by an acyclic flow:
  Δ ≤ Σ_j (Σ_κ w(j,κ) φ_{j,κ} − 1) with φ_{j,κ} ∈ [0, 1] the used fraction. For a pair this cannot happen
  (Corollary B1 telescopes).
- Per move: **a paid move contributes to the saving at most (its live classes that progress) + ½(those that
  half-progress) − 1.** A move on an array with one live class contributes at most 0, and −1 when that class does
  not progress.

**Calibration** [computed, `vflow_s2.out`; 152 rows at blocks 2–4: plain W = 1, 2, 3; pair trick; copy-move-copy;
gauge exchange; packing gadget; ±1 halving controls; order-3 controls C6O3, C8O3; 20 random order-3 packers with
W = 2, 3; each W = 2 case in the raw gauge, its coherent gauge if any, two random gauges, and the exact P* search].
- Every virtual move is an edge of Γ_m; raw max-flow ≥ W; Σ w·c ≥ Wm; Δ ≤ P_w ≤ P_live; Δ ≤ P_asc; for W = 2 the
  two-commodity LP with per-class capacities is feasible in every gauge. 0 violations.
- Negative control: with one unit of capacity per paid move (packing ignored), the two-commodity LP is infeasible in
  28 rows (random gauges of the pair trick, the halving controls, the order-3 controls and packers). Counting classes
  is necessary there. In the raw gauge the one-per-move flow still reaches W for these words: their raw packing is
  surplus.
- Exact P*: 0 for every coherent case; **m for C6O3, C8O3 and every random order-3 packer** (P* = 3 at block 3, 4 at
  block 4), with Δ = 0. These are the only known words with irreducible packing, and they refund all of it.
- No known engine unit at blocks 2–4 saves, so the inequality Δ ≤ P* is not tested at equality. What is tested is the
  lemma behind it.

## 3a. Coarser classes: image frames (the definition I recommend)

The proofs of Lemmas 2.1 and 2.2 need only that 𝒩 is a subgroup of B (in a mix, write b_s λ_κ = λ_{κ'} u with u in
the subgroup; normality is not used) and that the quotient keeps every source m steps from its sink. So the classes
can be coarsened.

**Theorem A′ (classes modulo any subgroup of the free group)** [proved; computed, `vflowB_s3.out`]. Let H be any
subgroup of B and h a gauge. Call the class of a term (g, l) the coset g h_l H, its *place* the double coset
B g h_l H, and let Γ_m/H̄ be the graph of places (H̄ acting on Γ_m through X ↦ ū⁻¹X). Suppose every place of a
source S_a is at distance m from every place of its sink T_a; this holds in the raw gauge for every H ⊆ B, because H̄
fixes L0. Then Lemma 2.2 holds on Γ_m/H̄ with 𝒜 replaced by the finite group algebra C[H], and Theorem A(a)–(c) hold
with live classes counted modulo H.
*Proof.* As for Lemmas 2.1 and 2.2: a term g y_l with g ∈ λ_κ H is λ_κ u y_l with u ∈ H, so c_{r,κ} = λ_κ Σ_l z_l y_l
with z_l ∈ C[H]; a letter acts on λ_κ from the left and leaves z alone; a mix with b_s ∈ B keeps the place, since
B b_s λ H = B λ H, and combines z's by left multiplication with elements of H; C[H] is semisimple of finite length
because H is finite. Places along a virtual move change by one edge of the quotient graph, and a path of places from
S_a to T_a lifts to a path of Γ_m of the same length, so it has length ≥ m under the hypothesis. ∎

Three choices matter, all in the raw gauge:
- H = 𝒩: Sp labels on Γ_m (sections 1–3).
- H = the free letters fixing both L0 and L1 (the Levi of the pair; it contains every coordinate permutation, since
  permutations commute with K): the class of g is the opposite pair (ḡL0, ḡL1), and the weight f is defined on
  places (it is invariant), so the progress weights of Theorem A(a) still apply.
- **H = B: the class of g is its image frame F = ḡL0** (G/B = Sp(2m,2)/P̄ is the set of Lagrangians). The places are
  the levels 0, 1, …, m, and the quotient graph is a path. Theorem A′ then reads:

  **for every level step i − 1 → i, at least W pairs (paid move, live image frame of the moved array at level
  i − 1 that the move ascends) occur. Hence Δ ≤ P_img := Σ_j (live image frames that move j ascends − 1).**

This is the strongest of the three in the raw gauge, since two Sp labels with the same image frame are one image-frame
class. I recommend it as the exact definition: **an array packs frames when its live terms have two or more distinct
image frames ḡL0, and a paid move carries packing when it ascends two or more of them.** Two terms that differ by a free
letter on the right (Pauli twins, (1 + Z^v)g, and the order-3 controls' (1 + s)g) have one image frame and are not
packing.

**Calibration** [computed, `vflowB_s3.out`, 264 rows at blocks 2–4: H = B in the raw gauge (image frames), H = B in
the gauge K⁻¹ (co-image frames, Corollary A′2) and the pair Levi]. The named cases,
12 random order-3 packers per block, a "crossing pair" (r0 = D1 x0 + x1, r1 = D1 x0 − x1, both moved along D2…Dm,
then (r0 ± r1)/2 and a last D1), and 36 random *crossing units* (input a moves along its own prefix set, a random
invertible Gaussian gate mixes all arrays, every array moves along a common set, the inverse gate unmixes, input a
moves along its own suffix set). 26 of the 88 words carry live packing (78 rows over the three class choices, up to
10 extra classes). Every level step is crossed at least W times by live classes, in both directions; Δ ≤ P_img and
Δ ≤ P_w hold; 0 violations. A second run (`vflowB_s4.out`: 75 random crossing units and 50 random order-3 packers,
459 rows, 50 words with live packing) also has 0 violations.
- The order-3 controls have **no image-frame packing at all** (their two Sp labels per array share the image frame),
  although their Sp-label packing is P* = m. Their "pack, move, unpack" is a free relabelling, as o3-build's Lemma 1
  says, and Theorem A′ sees it that way.
- The crossing units pack image frames live (their arrays hold inputs at different levels) and save nothing: the
  packing is surplus capacity, as in the hub's design without its centre gain.

**Corollary A′1 (exact accounting by chosen ascents)** [proved]. Choose, for each level step, exactly W of its live
image-frame ascents (possible by Theorem A′), and let u_j be the number chosen on move j. Then Σ_j u_j = Wm, so

  Δ = Σ_j (u_j − 1)  exactly:

the saving is the number of extra chosen ascents that shared moves carry, minus the number of moves that carry none.

**Corollary A′2 (co-image frames; the two end steps are never doubled)** [proved].
- The gauge h_a = K⁻¹ for every input, with H = B, also satisfies the hypothesis of Theorem A′: the class of a term g
  is its *co-image frame* F′ = ḡL1, the source place is level m and the sink place level 0. So, for every unit,
  every co-level step i → i − 1 (of d(ḡL1, L0)) is crossed at least W times by live co-image frames, and
  Δ ≤ Σ_j (live co-image frames that move j brings one step closer to L0 − 1).
- L0 is the only frame at level 0. So an array holds at most one live image frame at level 0 and at most one live
  co-image frame at co-level 0. Hence **no paid move lifts two live image frames off the raw frame, and no paid move
  brings two live co-image frames onto it** (a co-image frame is L0 exactly when the term is finished, g ∈ B·K up to
Paulis and scalars): at least W different moves do each of the two end steps. (A move can still
  combine an end-step crossing with crossings of other steps on other classes. At block 1 this is M ≥ W.)

**Corollary A′3 (units that never pack image frames never save, at any W)** [proved]. If every paid move acts on an
array with at most one live image frame (or at most one live co-image frame, or one live class modulo any subgroup
H ⊆ B with the hypothesis of Theorem A′), then M ≥ Wm. This is a class of units, of any width, for which "every unit
costs at least the best common-frame network" (sl-algebra Proposition 6) holds, trivially: they cost at least the
plain schedule. Common-frame networks with non-trivial twists are generally not in it, since in raw coordinates their
arrays hold inputs on different image frames.

**Lemma B′ (rigid image frames)** [proved]. The image frames of the classes of one register are all moved by the same
isometry: p̄ for a paid letter, b̄ ∈ P̄ for a free one. So their pairwise distances are fixed while they stay together,
and a pair at distance δ contributes at most (2m − δ)/2 to P_img over any chain of moves during which they are the
moved array's only live image frames (per move a_1 + a_2 − 1 ≤
(Δℓ_1 + Δℓ_2)/2, which telescopes, and δ ≤ ℓ_1 + ℓ_2 ≤ 2m).

## 3b. Relation to sl-algebra's Theorem 5 and Proposition 6

- **Theorem 5** (every unit: Δ ≤ Σ_j (μ_j − 1), μ_j the number of inputs whose chosen contributing history passes
  move j) counts *inputs* per move; Theorem A′ counts *live image frames* per move. Neither implies the other: an array
  holding x_1 + x_2 has one image frame and may serve two histories; an array holding (1 + D_1)x_1 has two image frames
  and serves one input. Both are valid at once, so a saver must satisfy both.
- **Proposition 6** ("every unit costs at least the cheapest common-frame network with the same W"). Theorem A(d)
  proves the relaxed form: for every gauge h, M ≥ NC(h) − P_live(h), where NC(h) is the least cost of an 𝒜-linear
  (frame-coherent) network for h's layout. For two arrays NC(h) = 2m (Hu), so **M ≥ 2m − P*: the whole gap between
  Proposition 6 and the truth for two arrays is irreducible live packing.** A unit refuting Proposition 6 with two
  arrays and saving k must carry at least k units of live packing in every Sp gauge (Theorem A(b)) *and* at least k
  extra live image-frame ascents in raw coordinates (Theorem A′), and its packed frames must not be refunded. The order-3
  controls meet the first condition (P* = m) and fail the second (no image-frame packing).
- Proving Proposition 6 for two arrays therefore amounts to proving that irreducible live packing is always paid back.
  I have no proof. The rigidity results (Lemma B, Corollary B1, Lemma B′) bound what one packed pair can carry, not how
  many pairs a unit can form. And payback is not forced class by class: a packed class can be split off by a canceller
  that is itself live through another copy (the crossing pair does exactly this: (r0 + r1)/2 and (r0 − r1)/2 cancel
  each other's packed classes and every class stays live). That mutual cancellation is coding, so a refund theorem
  would have to be a coding-gain bound for bundled networks [argued].

**Proposition D (two arrays, no scratch: the split must be a frame split)** [proved; this is sl-algebra's episode
argument (its Corollary 2′, upheld by v-slack) run with Sp labels in place of group elements]. Let a unit have exactly
two registers (no scratch), and suppose that no register ever holds one input on two different Sp labels (Paulis and
scalars may differ, and contents may have any 𝒜-coefficients). Then M ≥ 2m.
*Proof.*
1. *Contents that hold one input alone are full.* If at some time an array's whole content is λ·A·y_a (one label,
   A ∈ 𝒜), then A is invertible in 𝒜: every later content of input a lies in C[G]·A y_a, the output K y_a must be one of
   them, so A has a left inverse in C[G]; C[G] is a free C[𝒩]-module, so A is already invertible in 𝒜 = C[𝒩].
2. *Separate periods and episodes.* With two registers and both inputs present, at each time either each register
   holds exactly one input (a separate period) or some register holds both (an episode).
3. *One twist per episode.* In an episode let a register holding both inputs have labels λ¹, λ²; its twist is
   τ = (λ¹)⁻¹λ² ∈ Sp(2m, 2). Moves and free letters leave τ unchanged (left multiplication). In a mix or gate whose result
   receives input a from two summands, the two contributions sit on labels λ and b̄λ′, and the hypothesis forces them to
   be equal (two different labels cannot cancel each other's terms, since Sp labels of terms are a property of the
   group elements). So every register holding both inputs during one episode has the same twist, from the merge that
   opens the episode to its end.
4. *Inside an episode the unit is frame-coherent.* Take the gauge (h_1, h_2) = (1, τ⁻¹). Every register then holds at
   most one label (registers holding one input by hypothesis, registers holding both by step 3). Start the episode just
   before its opening merge and end it just after it closes; at both ends each input sits alone in one register on one
   label with full content (step 1). The cut lemma (Lemma 2.2, with these start and end states in place of sources and
   sinks) and Hu's theorem give: moves in the episode ≥ d_1 + d_2, where d_a is the distance between input a's positions
   at the two ends. Right multiplication by the gauge is an isometry, so d_a = ℓ(λ_end^a (λ_start^a)⁻¹), the distance
   input a's label travels, measured without any gauge.
5. *Separate periods:* each input moves alone, so the moves are at least the distance its label travels.
6. *Sum.* For each input the travelled distances add up, by the triangle inequality, to at least ℓ(K) = m. So M ≥ 2m. ∎
- This strengthens sl-algebra Corollary 2′ (one group element per input per array) to one Sp label per input per
  array: Pauli twins and 𝒜-coefficients are now allowed and still cannot save. It is incomparable with
  sl-combinatorics Theorem 4.3 (scratch allowed there, one global gauge needed) and with Theorem A(b).
- So every two-array, no-scratch saver must, at some moment, hold one input on two different Sp labels in one array
  (a *same-input frame split*). The order-3 controls do (two labels, one image frame) and save nothing; Theorem A′ then
  adds that the unit must also carry live *image-frame* packing at least equal to its saving.
- The same argument with image frames in place of Sp labels does not go through: a twist is then defined only up to
  free letters on both sides, which leaves no common gauge for an episode.

## 4. What one paid move can serve

**Lemma S (serving shares)** [proved; computed exactly at m = 2, 3, 4, `geom.json`]. Let a label have image frame F
at level k (U = F ∩ L0, dim U = m − k), and let a paid move have neighbour N of L0 (on the line through L0 whose
points contain the hyperplane H = N ∩ L0). The level after the move is d(F, N), and:
- it ascends iff U ⊄ H; this happens for 2(2^m − 2^k) of the 2(2^m − 1) neighbours;
- otherwise one of the two neighbours on the line descends and the other is lateral: 2^k − 1 each.
*Proof.* The Lagrangians containing H are L0 and two others, L = H + ⟨v⟩ with v̄ one of the three isotropic lines of
H^⊥/H. Since F is Lagrangian, F ∩ H^⊥ is the kernel of F → H*, of dimension 1 + dim(F ∩ H); its image in the
2-dimensional H^⊥/H is isotropic, hence exactly one line. So exactly one point L* of the line has
dim(F ∩ L*) = dim(F ∩ H) + 1, and the others have dim(F ∩ H) (a near-polygon). If U ⊄ H, dim(F ∩ H) = dim U − 1 and
L* = L0, so both other points are one step farther: ascent. If U ⊆ H, dim(F ∩ H) = dim U, L* ≠ L0 is a descent and the
third point is lateral. Hyperplanes containing U: 2^k − 1. ∎
In the position picture, a position on the interval between L0 and L1 at level k has 2^(m−k) − 1 progress
neighbours (Δf = 2) and 2^(m−k) − 1 half-progress neighbours (Δf = 1), out of 2(2^m − 1) [computed, m = 2, 3, 4];
the full profiles by (d(L0, X), d(L1, X)) are in `geom.json`.

**Proposition G (geometry does not cap packing)** [proved]. For any finite list of edges (X_i, X'_i) of Γ_m there
is a content whose classes are moved along exactly these edges by one paid letter. *Proof.* Fix a paid letter p, with
N = p̄⁻¹L0. Sp(2m, 2) is transitive on ordered edges, so there is σ_i with σ_i L0 = X_i and σ_i N = X'_i; take a label
λ_i ∈ G with λ̄_i = σ_i⁻¹ and the content Σ_i λ_i y. ∎ In particular a move can ascend every label at level 0 at
once, and by Lemma S a bundle of image frames at levels k_i < m is ascended by one common move whenever
Σ_i (2^{k_i} − 1) < 2^m − 1 (union bound over hyperplanes). The bound is sharp in kind: three frames at level m − 1
block every move (computed: the smallest blocking set is 3 at m = 2 and 3).

**Lemma B (bundles are rigid)** [proved]. If κ_1, κ_2 are classes of one register, τ = κ_1⁻¹κ_2 is unchanged by moves
and free letters on that register, and by mixes and gates that carry both from one summand. While it lasts, the
positions satisfy X_2 = τ̄⁻¹ X_1 and the image frames F_2 = κ̄_1 τ̄ L0. *Proof.* (xκ_1)⁻¹(xκ_2) = κ_1⁻¹κ_2; and
X_2 = (κ̄_1τ̄)⁻¹L0 = τ̄⁻¹X_1. ∎ So a register with several classes moves a rigid configuration: one choice of
neighbour for the first class fixes the step of every other one.

**Corollary B1 (what a pair can carry)** [proved; computed]. Along any chain of consecutive moves applied to one
lineage of a pair of live classes with relative label τ, while the moved register holds exactly these two live
classes, those moves contribute at most

  ρ(τ) = (max g − min g)/2,  g(X) = ½[f(X) + f_A(X)],  f_A = d(A, ·) − d(A', ·),  (A, A') = (τ̄L0, τ̄L1),

to the bound P_w of Theorem A(a). ρ(τ) = m exactly when (A, A') = (L0, L1), i.e. when τ fixes both the raw frame and
the transformed frame (τ in the Levi of the pair: free letters that commute with K up to free letters, such as the
coordinate permutations, like the order-3 controls' s).
*Proof.* For one move, with class increments x, y ∈ {−2, …, 2} of f and f(X_2) = f_A(X_1): the contribution is
x⁺/2 + y⁺/2 − 1 ≤ (x + y)/4 = Δg/2 (check the three sign cases: both ≥ 0, since x + y ≤ 4; one negative, since
x − 4 ≤ y; both negative, since x + y ≥ −4). Telescope. g ≤ m with equality only at X = L1 = A'; g ≥ −m with equality
only at X = L0 = A. ∎
Computed over all ordered opposite pairs (A, A'), which is the uniform distribution over τ ∈ Sp(2m, 2) (`geom.json`):

| block | pairs | ρ = m | mean ρ | longest path progressing for both frames: distribution |
|---|---|---|---|---|
| 2 | 120 (all) | 1 | 1.18 | 0: 30, 1: 89, 2: 1 |
| 3 | 8,640 (all) | 1 | 1.77 | 0: 127, 1: 7,329, 2: 1,183, 3: 1 |
| 4 | 300 (sampled) | 0 | 2.31 | 0: 1, 1: 10, 2: 287, 3: 2 |

A random pair bundle can pay for at most about two moves over its whole joint life at blocks 3 and 4. Only one
relative position in 120 (block 2) or 8,640 (block 3), the Levi type, can carry a full m; the order-3 controls use
exactly that one and refund it all at their readout (and their two labels share an image frame, so Theorem A′ does not
count them as packing at all).

## 5. Ceilings: what packing capacity gives, and what the proved floors say in packing terms

**Proposition C (no ceiling from capacity alone)** [proved for the first two sentences; argued for the last].
Geometry puts no bound on how many live classes one move can carry or ascend (Proposition G), and Theorem A bounds
the saving by the packing, not the packing by anything. If every moved array is limited to at most k live classes,
Theorem A(a) gives Wm ≤ Σ_j k_j ≤ kM, i.e. a saving of at most ln k / ln m (0.5 at block 4 for pairs), far above
every proved ceiling. So every per-block ceiling must come from the cost of *unpacking* (the refund), and the only
proved measure of that cost is Theorem R's potential, whose cap holds for packed and unpacked contents alike.

**Corollary C1 (the proved floors as packing ceilings)** [proved from the floors; `tables.out`]. For every unit at
block m, Δ/M ≤ κ_m := m/L_m − 1:

| block | floor L_m | largest saving | κ_m: net saved moves per paid move | 0.1 needs (per move) | 0.01 needs (per move) |
|---|---|---|---|---|---|
| 2 | 2 | 0 | 0 | 0.072, excluded | 0.0070, excluded |
| 3 | 14/5 | 0.0628 | 1/14 = 0.0714 | 0.116, excluded | 0.0110 |
| 4 | 41/13 | 0.1714 | 11/41 = 0.268 | 0.149 | 0.0140 |
| 5 | 2.99172 | 0.3191 | 0.671 | 0.175 | 0.0162 |
| 6 | 2.97290 | 0.3919 | 1.018 | 0.196 | 0.0181 |
| 8 | 2.9447 | 0.4806 | 1.717 | 0.231 | 0.0210 |
| 16 | 2.9447 | 0.6105 | 4.43 | 0.320 | 0.0281 |

"Needs" is m^s − 1: by Theorem A(a) a unit with saving s must carry, per paid move, at least that many extra live
ascending classes in raw coordinates, and all of it must survive the readout.

## 6. Sources

- `../../slack/sl-combinatorics/PROOF.md` (sections 0–4: setting, Lemma 1.3, Lemma 4.1, Lemma 4.2, Theorem 4.3,
  Remark 4.3b, Corollary 4.4); `../../slack/sl-algebra/PROOF.md` (Corollary 2, Theorem 5); `../../proofs/NOTE.md`,
  `../../proofs/ladder/PROOF.md` rung 1 (Theorem R), `../../meta/MAP.md` §1.2 (floors at blocks 5, 6 and ≥ 4).
- T. C. Hu, Multi-commodity network flows, Operations Research 11 (1963) [cited]. Papernov–Lomonosov classification
  (Schrijver, Combinatorial Optimization, multicommodity flows) [cited, not re-read].
