> **Status now (19 October 2026).** This is the working proof file of the entropy-bounds round of 15 October, as the AI agent wrote it, with private paths removed. Its own first lines say no skeptic had read it; afterwards one independent AI skeptic read it and checked Lemma W, Lemma W+ (with its near-polygon claim, by brute force), the symmetry facts that Lemma 1 needs, the exact screen and every quoted LP value it recomputed, and upheld them with wording corrections ([entropy-bounds-skeptic.md](entropy-bounds-skeptic.md)). It wrote no separate check of the proofs of Lemmas 2 to 4 (lumping). So: Lemma W and W+ argued by hand, not checked, one skeptic re-derived them; Lemmas 1 to 4 argued by hand, not checked; the numbers checked by computer. No person has read this file and nothing in it is checked in Lean. "I" in the file is the agent that wrote it.

# PROOF.md: entropy-LP bounds on the coding gain of coherent units (nc-gain)

nc-gain, packing round, 15 October 2026. No skeptic has read this. Labels: [proved] (proof written here),
[computed] (a program in `calc/` or a `job_*` folder and its output), [cited] (a published theorem, cited from
memory), [read] (from one of our notes), [search] (a heuristic search; what it finds is achievable, but it may miss
better objects).

## 0. Setting

- Γ_m is the frame graph (dual polar graph of Sp(2m, 2)): N_V frames (15, 135, 2295 at m = 2, 3, 4), degree
  2(2^m − 1), arc set A (two arcs per edge). Opposite = distance m. U is the *uniform layout*: every ordered opposite
  pair (S, T) once; |U| = N_U = 120, 8640, 2,350,080 at m = 2, 3, 4 [computed, `calc/dpg.py`].
- **Codes.** A layout D is a list of sessions (S_a, T_a) with S_a, T_a opposite. A *code* for D is what
  sl-combinatorics Lemma 4.1 calls the virtual network of a unit with no frame packing on live moves: every node
  (frame) keeps everything it has received and combines it linearly for free; each move sends one element of 𝒜^W
  (𝒜 = M_n(C), n = 2^m) along one arc; input a starts at S_a and must be recovered at T_a. By Morita equivalence a
  left 𝒜-submodule of 𝒜^W is the same as a subspace of C^(nW) (its rows), with composition length = dimension, so
  this is a vector-linear network code over C with symbols of dimension n and sources of dimension n. Scalar codes
  (strictly coherent units, sl-combinatorics Remark 4.7) are the case n = 1. Normalising every dimension by n,
  each move carries capacity 1 and each source has rate 1.
- cost*(D) = the least number of moves of a code for D. Routing (no mixing of different inputs) costs exactly
  m·|D| (each input travels distance m; sl-combinatorics Remark 4.7). The *coding gain* is γ(D) = m|D| / cost*(D),
  and γ_m = sup over layouts D of γ(D).
- Every unit of the engine game that has an Sp-gauge without frame packing on live moves is such a code with
  cost = M (sl-combinatorics Lemma 4.1 and Remark 4.7) [read], so **M ≥ W·m/γ_m for all of them**, and their saving
  is at most ln(γ_m)/ln m.

## 1. The worst layout is the uniform one

**Lemma 1** [proved]. Let LB be any lower bound on cost* that is homogeneous (LB(c·D) = c·LB(D) for c copies of each
session) and valid for every code. Then for every layout D,

  cost*(D) ≥ |D| · LB(U) / N_U.

In particular γ_m ≤ m·N_U / LB(U), and γ_m = m·N_U / cost*-rate of U (the best per-session cost of codes for many
copies of U).

*Proof.* Sp(2m, 2) acts on Γ_m by graph automorphisms, transitively on ordered opposite pairs (the stabiliser of
L0 acts simply transitively, through its unipotent radical, on the frames opposite L0). For g in the group G =
Sp(2m, 2), moving every position of a code for D by g gives a code for gD with the same number of moves. Running the
|G| codes side by side on disjoint inputs gives a code for ⊎_g gD, which is c copies of U with c = |D|·|G|/N_U, at
cost |G|·cost*(D). So |G|·cost*(D) ≥ cost*(c·U) ≥ LB(c·U) = c·LB(U). ∎

So one LP per block, on U, bounds every layout. Every LP below is homogeneous in (entropies, rates), and its
optimum may be taken G-invariant (average an optimal point over G), which is what licenses orbit variables.

## 2. The entropy LP and why each row is valid

Ground set: one element Y_a per session (or group of sessions with the same source node and sink node) and one
element X_e per arc (or group of arcs). For a code, h(S) = dimension of the span of the symbols in S (normalised by
n). This h is the rank function of a subspace arrangement over C, so it satisfies:
- the polymatroid (Shannon) inequalities: h(∅) = 0, monotone, submodular;
- the Ingleton inequality (all linear rank functions over any field);
- DFZ inequality (65) of arXiv:1311.4601, Theorem 8.6 (linear rank functions over fields of characteristic ≠ 2;
  our field is C). The encoding is the one checked in `../../proofs/duality/dfz_check.py` [read].

Network rows, valid for every code:
- (src) the sources are independent: h(Y_S) = Σ_{a∈S} r_a;
- (fun) every symbol sent is a function of all sources: h(everything) = h(all Y);
- (closure) for every node set Q: h(In(Q)) = h(In(Q) ∪ Out(Q) ∪ Dec(Q)), where In(Q) = sources at Q and arcs
  entering Q, Out(Q) = arcs leaving Q, Dec(Q) = sessions with a sink in Q. Valid because everything a node of Q
  ever holds was a source in Q or came in through an arc entering Q (induction over time);
- (cost) the moves along arc e number at least h(X_e), so cost ≥ Σ_e h(X_e).

LP(N, D) = min Σ_e h(X_e) over all h satisfying the rows; it is a lower bound on cost*(D). With an undirected
graph each direction of an edge is a separate element: a symbol sent both ways is paid twice (the exchange
u ↔ v with sessions both ways costs 2, not 1: checked, `calc/calib_small.py`, case 4).

## 3. Contractions, and the ceiling of every LP that uses only some node sets

**Lemma 2 (contraction)** [proved]. Let π be a partition of the frames into parts, Γ/π the network whose nodes are
the parts and whose arcs are the arcs of Γ_m between different parts, and U/π the sessions of U between different
parts. A code for U on Γ_m induces a code for U/π on Γ/π whose cost is the number of moves along crossing arcs.
Averaging over G (Sp(2m, 2) is arc-transitive on Γ_m, being distance-transitive):

  cost*(U) ≥ cost*(U/π on Γ/π) · |A| / |A_π|,   A_π = arcs crossing π.

*Proof.* Pool the nodes of each part (their stores merge, internal arcs become free). For the average: for each g,
the moves along arcs crossing gπ number at least cost*(U/π) (U is G-invariant); summing over g, each arc is counted
|G|·|A_π|/|A| times. ∎

**Definition.** R(π) = Σ_{(s,t) ∈ U} d_π(π(s), π(t)), d_π the distance in the quotient graph (parts adjacent when
an edge of Γ_m joins them), and ρ(π) = (R(π)/N_U) · (|A|/|A_π|): the routing value of the contraction, per session,
after averaging.

**Lemma 3 (ceiling of closure-family LPs)** [proved]. Let 𝒬 be any family of node sets of Γ_m, and π_𝒬 the partition
of the frames into the cells of its Venn diagram. Consider any LP whose network rows are (src), (fun), (cost) and
the (closure) rows for Q ∈ 𝒬 only, on any grouping of the ground set for which those rows can be written, together
with ANY further inequalities valid for all linear rank functions over C (Shannon, Ingleton, DFZ, all linear rank
inequalities, or even only entropic ones). Its value is at most R(π_𝒬). If it is used through Lemma 2 (averaged
over G), the bound it certifies is at most ρ(π_𝒬) moves per session.

*Proof.* Route every session of U/π_𝒬 along a shortest path of the quotient graph, each quotient step realised on
one Γ_m-arc between those cells, and send nothing on arcs inside cells. This is a code on the quotient network (a
routing; it is not a code on Γ_m, because it skips the travel inside cells). Put on the original ground set (arcs
inside cells get the zero subspace, each session its own coordinate), it is the rank function h of a subspace
arrangement, so it satisfies every inequality valid for linear rank functions. It satisfies (src) and (fun), and
(closure) for every Q that is a union of cells, in particular for Q ∈ 𝒬: a routed symbol leaving Q either started
in Q or entered it earlier, and a session with its sink in Q ends there. (A closure row for a set that splits a cell
could fail, which is why the cells are those of the Venn diagram of all sets used.) Its cost Σ_e h(X_e) is the
number of quotient steps, R(π_𝒬). Grouping ground elements only removes rows. With averaging, the certified bound is
LP·|A|/(|A_π| N_U) ≤ ρ(π_𝒬). ∎

So an entropy LP can certify β moves per session on Γ_m only if its closure family separates the frames finely
enough that the routing value of the Venn contraction is at least β. The size of the LP grows with the number of
cells: every pair of cells joined by sessions is a separate source element (merging sources with different source
or sink cells breaks the decoding rows), and every pair of adjacent cells gives two arc elements.

**Lemma 4 (at most four cells: never better than cuts)** [proved, citing a classical fact]. If π has at most four
parts, ρ(π) ≤ c_m, the cut value (best single cut: 3/2, 7/4, 15/8 at m = 2, 3, 4 [read: sl-combinatorics `cutconst`,
v-slack]).

*Proof.* Every metric on at most four points is a non-negative combination of cut metrics (the metric cone and the cut
cone coincide for n ≤ 4; Deza–Laurent, *Geometry of Cuts and Metrics*) [cited]. Write d_π = Σ_S λ_S δ_S over cuts S
of the parts. Pulling back, each δ_S is the cut metric of the union of the parts in S. Crossing arcs have d_π = 1 and
internal arcs 0, so |A_π|/|A| = Σ_S λ_S·(share of arcs cut by S), and R(π)/N_U = Σ_S λ_S·(share of opposite pairs
split by S). The ratio of the two sums is at most the largest ratio of the terms, which is at most c_m. ∎

The first contractions that can beat the cut bound have five cells, with a non-l1 quotient metric (a K_{2,3}
pattern: the pentagonal inequality fails).

## 4. Calibration of the LP [computed]

Programs: `calc/nclp.py` (orbit variables; closure-reduced form, in which every set is identified with its
functional closure, valid because h(S) = h(cl S) follows from the closure rows by submodularity; full Shannon rows;
Ingleton and DFZ (65) cutting planes; HiGHS simplex, or HiGHS PDLP with a rigorous dual bound: for any multipliers
y, with the ≤-rows clipped to y ≤ 0, LB = b_eq·y_eq + Σ_j min(0, (c − Aᵀy)_j)·R_tot is a valid lower bound because
every variable is some h(S) ≤ h(all) = R_tot).

| network | routing | cut LP | Shannon LP | + Ingleton, DFZ | where |
|---|---|---|---|---|---|
| path of length 3, one session | 3 | 3 | 3 (simplex; PDLP bound 3) | — | `calib_small.py`, `test_pdlp.py` |
| directed butterfly, 2 sources multicast to 2 sinks | 8 (Steiner) | — | 7 (simplex; PDLP bound 7) | — | same |
| one edge, sessions both ways | 2 | 2 | 2 | — | `calib_small.py` (case 4, by hand: 2 arcs, each ≥ 1) |
| Okamura–Seymour K_{2,3}, sessions x→y, y→z, z→x, a→b | 8 | 6 | ≥ 7.9999986 (PDLP rigorous bound; full LP 7.9999998), so every code costs ≥ 8 | ≥ 7.9999992 (8 Ingleton rows violated at the PDLP point by 1.4e-2, none after; no DFZ row violated) | `job_os_closed/`, `job_os_pdlp/` |
| same four session types, rates free | cost/routing 1 | — | ≥ 0.99999988: no coding gain for any rates | — | `job_os_closed/` |
| 3-cell contraction of Γ_2 (a line, the far set, the rest) | 80 | 80 | ≥ 79.99999 | Ingleton rows violated by ≤ 1.4e-7 (round-off), value unchanged | `job_cells/` |
| K_{2,3} with A→B 8, B→A 8, X_i→X_j 4 each (the distance-2 skeleton of the 5-cell lumping of Γ_2; 20 elements, symmetry of order 12; 64,065 variables, 3,305,746 rows) | 80 | 64 | 70.79995 (rigorous ≥ 70.79996), **short of routing** | not tried (one solve took 46 min) | `job_k23two/` |

The last row shows the limit of the method even on K_{2,3}: with sessions in both directions, Shannon's inequalities with
the closure rows stop at 70.8 against 80. So even the 5-cell lumping of Γ_2, if it could be solved, would probably
not certify its routing value 1.6. The OS row is the calibration that matters: the smallest undirected network on which cuts are not enough (four
sessions are needed, `calc/k23gap.py`), where the entropy LP certifies the routing value, as Harvey, Kleinberg and
Lehman proved [cited].

## 5. The ceiling curves [computed, search]

ρ(π) for the best contractions found (lower estimates of the best possible ρ with that many cells):

- Γ_2 (cut 3/2, certificate 2): cells 2–4: 1.5 (Lemma 4: no more possible); 5: 1.6 (the K_{2,3} star: A = one frame
  v, X_1..X_3 = the other two frames of each line through v in GQ(2,2), B = the 8 frames opposite v; quotient
  routing 128, quotient cut LP 112, LP elements 12 arcs + 14 session groups = 26, `calc/star.py`); 6: 1.676;
  8: 1.737; 9: 1.846; 13: 1.849; 14: 1.909 (every one of the 105 contractions merging one pair; exhaustive); 15: 2.
  No contraction with ρ > 1.5 and at most 20 LP elements was found (`calc/small5.py`: the smallest has 26).
- Γ_3 (cut 7/4, certificate 14/5): cells ≤ 13: nothing found above 1.75; 14: 1.784; 15: 1.828; 20: 1.933;
  30: 2.256; 60: 2.426; 100: 2.729; 120: 2.777; 123: 2.795; 124: 2.811 (first ≥ 2.8); 134 (one pair merged): 2.956
  to 2.966. The four-cell contractions that reach 1.75 have quotient cut LP = quotient routing (as Lemma 4 says).
- Γ_4 (cut 15/8, certificate 41/13 = 3.1538): distance-profile contractions with up to 300 cells: at most 1.803
  (at 198 cells), below the cut value; greedy merging from the full graph (30 sampled candidates per merge): 3.993
  after 1 merge, 3.531 after 100, 3.267 after 200, 3.1548 at 2,027 cells, 3.1534 at 2,026 cells (first below the
  certificate), 3.024 at 1,895 cells.

## 6. Ceilings on the saving of coherent units

For a coherent unit (or any unit without frame packing on live moves in some gauge), M ≥ W·m/γ_m, so the saving is
at most ln(γ_m)/ln m. Bounds on γ_m: from cuts m/c_m = 4/3, 12/7, 32/15 (saving ≤ 0.415, 0.491, 0.547); from every
entropy LP of solvable size, the same (§3, §5); from Theorem R (all units) m/L_m = 1, 15/14, 52/41 (saving ≤ 0,
0.0628, 0.1714) [read: the certificate floors 2W, 14W/5, 41W/13].

## 7. Fixed small layouts: the cut LP per layout

The targets are units with few arrays, and for a fixed layout D the uniform-layout bound is far too pessimistic.
For a fixed D, the *directed cut LP*

  cut(D) = min Σ_arcs c  subject to  c(arcs entering Q) ≥ #{a : T_a ∈ Q, S_a ∉ Q}  for every node set Q

is a lower bound on the moves of every code for D (the cut-set rows of §2; equivalently sl-combinatorics Lemma 4.2
in its directed form) [proved: each row holds for every code]. It is computed exactly by constraint generation
(`calc/layout_cut2.py`): an LP over c alone, and for each pattern P of endpoints inside Q a max-flow from the
endpoints outside Q to P; a min cut below the requirement is added as a row. The max-flow runs on ⌊10^6 c⌋, so the
final "no violated row" is rigorous up to 10^-6 per row.

**Consequence.** If cut(D) > M for every layout D of W opposite pairs on Γ_m, then no unit with W arrays and M moves
at block m that is coherent (no frame packing on live moves) in some gauge exists: in that gauge it is a code for
some such D. Searches over layouts (`calc/layout_search.py`) can only find layouts where cuts fail; they cannot
prove the "every layout" statement.

Which (m, W, M) matter, given the certificates (M ≥ 14W/5 at block 3, M ≥ 41W/13 at block 4) [read]:
- block 3: W ≤ 4 cannot save at all (⌈14W/5⌉ = 3W for W ≤ 4); W = 5 allows 5 in 14; W = 6 allows 6 in 17;
- block 4: W = 3 allows 3 in 10 and 3 in 11; W = 4 allows 4 in 13 to 4 in 15.

Results [computed; exact constraint generation `layout_cut2.py`, or PDLP with a rigorous dual bound
`layout_cut_pdlp.py`: every variable boxed in [0, W·m], which some optimal solution respects]:
- block 2: W = 3, 4: cut value = routing on every sampled layout (about 90 each with local search); W = 5:
  (0→7), (3→1), (4→10), (0→10), (10→7) has 9 < 10; W = 6: 11 < 12. Theorem 2W makes routing optimal at block 2, so
  from five arrays on the cut rows are not enough there.
- block 3, W = 5: random layouts: all 8 samples certified > 14 (14.005 to 14.036). Layouts with three mutually opposite frames
  a, b, c, sessions a→b, a→c, c→b and two pendant sessions: of 30, 22 > 14, 7 exactly 14, 1 exactly 13 (also 13.00000
  by PDLP). So the cut rows cannot exclude a coherent 5 in 14; the certificate floor 14 is the stronger bound.
- block 3, W = 6: one random layout certified > 17.
- block 4, W = 3: a random layout and a "parallel" one (each pair's ends adjacent to the previous pair's ends): cut
  value 12.00000 = routing (rigorous), 10–20 minutes each on 2.6 million variables.
- For W = 3, only three pairwise disjoint pairs can have a cut gap: every other demand graph with three edges is a
  subgraph of K4 or of a union of two stars, where the cut condition implies a fractional multiflow
  (Papernov–Lomonosov) [cited from memory], hence cost ≥ 3m.

## 8. Session 2: few arrays. Lemma W (two-class stretch) and the pull-back screen

Session 2 of nc-gain, 15 October 2026, from 15:20 EDT. No skeptic has read this section.

Setting as in §7, with the UNDIRECTED cut rows of sl-combinatorics Lemma 4.2: for a unit with no frame packing on
live moves in some Sp-gauge, with layout D = ((s_1, t_1), ..., (s_W, t_W)) (opposite pairs, so d(s_j, t_j) = m =
diameter of Γ_m), let c(e) be the number of paid moves along the edge e (either direction; a grouped move counts its
rank, Remark 4.3b). Lemma 4.2: c(δ(Q)) ≥ sep(Q) := #{j : exactly one of s_j, t_j in Q} for every set Q of frames.

**Fact (cut packings)** [proved, LP duality made explicit]. Let δ = Σ_Q λ_Q δ_Q (λ ≥ 0, δ_Q the cut metric of Q) with
δ(e) ≤ 1 on every edge. Then M ≥ Σ_j δ(s_j, t_j).
*Proof.* M = Σ_e c(e) ≥ Σ_e c(e) δ(e) = Σ_Q λ_Q c(δ(Q)) ≥ Σ_Q λ_Q sep(Q) = Σ_j δ(s_j, t_j). ∎
Every line metric |φ(u) − φ(v)| of a real function φ is such a combination (λ = length of the level interval), and so
is any sum of them.

**Definition.** A set A of pairs of D is *stretchable* if each pair j ∈ A can be oriented (low_j, high_j) so that
d(low_i, high_j) = m for all i, j ∈ A (for i = j this is the pair itself). Two pairs are *compatible* when {i, j} is
stretchable: d(s_i, t_j) = d(s_j, t_i) = m, or d(s_i, s_j) = d(t_i, t_j) = m.

**Lemma W (two-class stretch)** [proved]. If the pairs of D split into two stretchable classes A and B, then every
unit with layout D and no frame packing on live moves in some Sp-gauge has M ≥ W·m. (Any graph of diameter m with
pairs at distance m; nothing about Γ_m is used.)

*Proof.* φ_A(v) = min_{i ∈ A} d(low_i, v) is 1-Lipschitz on Γ_m (a minimum of 1-Lipschitz functions), φ_A(low_j) = 0,
and φ_A(high_j) = min_i d(low_i, high_j) = m for j ∈ A. The same for φ_B. Put
  δ(u, v) = max(|φ_A(u) − φ_A(v)|, |φ_B(u) − φ_B(v)|) = ½|Δ(φ_A + φ_B)| + ½|Δ(φ_A − φ_B)|
(the l∞ plane is isometric to the l1 plane by (x, y) ↦ ((x + y)/2, (x − y)/2)). It is a sum of two line metrics,
δ(e) ≤ 1 on every edge, and δ(s_j, t_j) ≥ m for every j (by φ_A if j ∈ A, by φ_B if j ∈ B). Apply the Fact. ∎

**Converse for one class** [proved]. A single 1-Lipschitz φ with |φ(s_j) − φ(t_j)| = m for all j ∈ A exists iff A is
stretchable: φ(high_j) − φ(low_j) = m with φ(v) − φ(u) ≤ 1 on edges is a system of difference constraints; it is
feasible iff the constraint graph has no negative cycle, and a cycle through the pairs i_1, ..., i_k weighs
Σ_r (d(low_{i_r}, high_{i_{r+1}}) − m) ≤ 0, with equality only if every cross distance is m. So any construction of
the form "max of two 1-Lipschitz functions" needs Lemma W's hypothesis; non-compatible layouts need other metrics.

**Corollary W3** [proved]. W = 3:
- if two of the three pairs share a frame, they are compatible (orient both with that frame as low): M ≥ 3m. This
  covers every three-pair demand graph except three disjoint pairs, without Papernov–Lomonosov (cf. Cor. 4.4(b));
- if the six frames are distinct and some two pairs are compatible: M ≥ 3m.
- With Theorem 4.3 (W ≤ 2): every unit with W ≤ 3 arrays, no live frame packing in some gauge, and a layout in which
  two pairs are compatible, has M ≥ W·m, at every block m.

**The pull-back screen** [proved method; computed results]. For frames u_1, ..., u_k let ψ(v) = (d(u_1, v), ...,
d(u_k, v)) and H the image graph (vertices ψ(Γ_m), an edge for each pair of distinct images of adjacent frames).
ψ is non-expanding and, if u_j is an end of pair j, keeps d_H(ψ s_j, ψ t_j) = m. Any cut packing of H pulls back to
a cut packing of Γ_m with the same values on the pairs, so the cut LP of H (at most (m+1)^k vertices) is a lower
bound for Γ_m. The dual of the final LP is rationalised and scaled so that every edge load is ≤ 1 in exact
arithmetic (`calc/w3screen.py: exact_value`); since M is an integer, an exact value > W·m − 1 proves M ≥ W·m.
Levels tried: ψ to one end of each pair (2^W choices), then to all 2W ends, then Γ_m itself.

Calibration [computed]:
- Lemma W never fails where it applies: block 2, 374 stretchable layouts with W = 3..6 (shared frames allowed):
  undirected cut LP = W·m on every one (`calc/lemma_calib.py`); the cut gaps (W = 5: 5 of 64 non-stretchable
  layouts; W = 6: 11 of 114) are all on non-stretchable layouts. The known gap layouts of §7 (block 2 cut 9; block 3
  cut 13 and 14) are not two-class stretchable.
- Random graphs, three disjoint pairs at the diameter (`calc/w3rand2.py`, 1,070 instances on 9 to 12 vertices): 52
  have a cut gap (values 7.67 to 8.5 against 9); every one is non-compatible. So Lemma W uses only the diameter,
  and non-compatible layouts can have gaps in general graphs: whatever excludes them on Γ_m uses Γ_m.
- Block 3, the non-compatible layout [(4,123), (25,33), (78,48)]: best pull-back to three frames 8.5 (> 8, so
  M ≥ 9 by integrality), to all six frames 9, Γ_3 itself 9. Block 3 [(34,15), (65,28), (126,113)]: 9 at every level.

**Lemma W+ (moving an end)** [proved]. For a layout D of three opposite pairs let the *cross deficiency* be
  e(D) = min over two of the pairs {i, j}, and ends u of i and w of j (partners u', w'), of (m − d(u, w')) + (m − d(w, u')).
Then every unit with layout D and no live frame packing in some Sp-gauge has M ≥ 3m − e(D). (e(D) = 0 is Lemma W.)

*Proof.* (a) If D' differs from D by replacing one end x of a pair by a frame x' still opposite the other end, then
for every cut packing δ' (δ' ≤ 1 on edges, hence δ' ≤ d): Σ_k δ'(D_k) ≥ Σ_k δ'(D'_k) − δ'(x, x') ≥ Σ_k δ'(D'_k) − d(x, x').
(b) Claim: if y is opposite u and d(w, y) = k < m, then y has a neighbour y' opposite u with d(w, y') = k + 1.
  - Lines through y are ℓ_K = {the three Lagrangians containing K} for hyperplanes K of y. For any frame z,
    dim(z ∩ K^⊥) = 2m − dim(z + K) = dim(z ∩ K) + 1, so (z ∩ K^⊥ + K)/K is one of the three Lagrangian lines of the
    plane K^⊥/K. Hence ℓ_K has exactly one point nearest to z, at distance m − 1 − dim(z ∩ K); the other two are
    one step farther (the near-polygon property).
  - w ∩ y has dimension m − k ≥ 1; take a hyperplane K of y not containing it (there are 2^m − 2^k of them). Then y
    is the point of ℓ_K nearest to w, and the other two are at distance k + 1 from w.
  - u ∩ K = 0 (u is opposite y), so the point of ℓ_K nearest to u is at distance m − 1, and it is not y. Of the
    other two points, one is at distance m − 1 from u and one, y', is opposite u. ∎ (claim)
(c) Take the pairs and ends realising e(D) = a + b, a = m − d(w, u'), b = m − d(u, w'). Move u' a steps by (b) (each
    step keeps it opposite u and increases its distance from w), then w' b steps (opposite w, away from u; this does
    not change d(w, u')). The new layout D' is compatible, so its cut value is 3m (Lemma W), and by (a) the cut value
    of D is at least 3m − a − b. Apply the Fact. ∎

Distribution of e(D) for random layouts with six distinct frames [computed, `calc/edef.py`]: block 2 (3,000):
e = 0 85 %, 1 15 %; block 3 (3,000): 0 78 %, 1 21 %, 2 0.3 %; block 4 (20,000): 0 74.5 %, 1 24.9 %, 2 0.55 %.
Adversarial hill-climbing reaches e = 1 at block 2, 2 at block 3, 4 at block 4. So Lemma W+ alone gives M ≥ 3m − 1
for about 99 % of random layouts but only M ≥ 3m − 4 in the worst case found at block 4; the pull-back screen is
what decides the rest.

**Screen results** [computed; every certificate is an exact rational cut packing; files under `job_w3s*/`,
`job_edef4/`, `job_lpd4/`, `job_hard*/`].
- Block 3, three arrays, 100 random layouts with six distinct frames (`job_w3s3`, fractional mode): 78 compatible;
  the 22 others all have cut value exactly 9 (19 by a three-end image, 3 by the six-end image). The hill-climbed
  layout [(114,97), (17,41), (16,29)], all cross distances 2 (e = 2): three-end image 9, Γ_3 itself 9.
- Block 4, three arrays, 709 random layouts with six distinct frames (`job_w3s4a`–`d`, seeds 11–18): 531 compatible;
  177 of the other 178 certified M ≥ 12 by a three-end image (exact values 11.125 to 12, 25 of them exactly 12); the
  last one, [(2051,594), (2031,712), (777,952)] (e = 2), defeated the three-end images and was solved on all of Γ_4:
  directed cut value 12.00000 (`job_hardpd3/`). None open.
- Block 4, session 1's ten layouts solved on all of Γ_4 (directed cut LP, PDLP with a rigorous bound; the directed
  rows of Lemma 4.2 are valid for coherent units and at least as strong as the undirected ones): all 12.00000.
  Seven are compatible; three are not (e = 1: [(0,1847),(1273,27),(87,1938)], [(0,1754),(26,457),(188,47)],
  [(0,1357),(11,1348),(1263,66)]), so non-compatible layouts with cut value exactly 12 exist (`job_lpd4/compat_check.txt`).
- Block 4, hill-climbed layouts with e = 2 or 3 (`job_edef4/h4b.out`, `h4c.out`): of 46 (40 with e = 2,
  6 with e = 3), 22 certified by three-end images (11.125 to 12); 24 not (best three-end values 10.5 to 11.0). Eight of
  those 24 were solved on
  all of Γ_4 with the directed LP (`job_hardpd*/hardpd*.out`): [(331,751),(1319,596),(837,1071)] 12.00000,
  [(2115,615),(1327,2101),(641,486)] 12.00000 (both e = 3), [(2030,1063),(1397,1085),(1249,145)] 11.99999,
  [(487,1544),(722,1494),(2264,721)] 12.00000, [(1291,349),(1486,674),(108,1841)] 12.00000,
  [(1139,722),(987,1071),(2122,28)] 12.00000, [(649,2184),(60,1126),(1376,646)] 11.99999, [(1080,93),(3,906),(632,1269)] 12.00000 (rigorous PDLP bounds).
  Sixteen were not sent.
- Four arrays, block 3: [(0,95), (0,80), (80,95), (84,60)] (a triangle of mutually opposite frames and one pendant;
  not two-class stretchable) has undirected cut value 11 < 12: cut gaps begin at W = 4 (below the certificate floor
  ⌈56/5⌉ = 12 there, so this opens nothing). Block 4: 1,147 of 3,000 random four-pair layouts are two-class
  stretchable (M ≥ 16 proved); the pull-back screen did not finish any other four-pair layout in the time I had.
- Five arrays, block 3: 239 of 2,000 random layouts are two-class stretchable (no coherent 5 in 14 on them).
- The king graph {0..m}^3 (which needs only nine distances) is too coarse: 2 of 15 random admissible configurations
  at m = 3 fall to 7.67 and 8 (`calc/kingcfg.py`), so the image graphs themselves are needed.
