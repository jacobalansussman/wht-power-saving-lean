> Working note written by an AI agent during the round of 15 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# an-hub: why the hub's complex units save, and what caps them

an-hub, 15 October 2026, 08:11 to 08:30 EDT (Mac clock). Source: the community hub
CrocSwap/integer-mult-bounds. I read its local mirror at commit 0b1cccf8 (14 October, 02:49 EDT), the #1839
package (T74 supplier, the record chain), the T74 words, and OpenAI's paper "Integer multiplication below n log n",
section 3. No skeptic has read this note. Nothing was posted.

Labels: [read] = read first-hand in the hub or OpenAI sources; [computed] = my own small script, output in `out/`;
[argued] = my reasoning, checked by nobody; [derived here] = arithmetic done for this note.

## In short

1. **The object.** The hub's complex unit is one fixed local scalar program, the "word": about 12,500 additions and
   dyadic scalings with R = 480 helper roles, at block 73. Its geometry is E8: 112 ports, which are the 112 lines through E8
   roots of the shape (±1,±1,±1,±1,0,0,0,0) on the 14 tetrads of the Hamming code, plus 8 centres. The template
   runs the word in 9 stages ("windows") over 2 × 4 × 112 = 896 data families per replica. T74 uses 74 replicas,
   666 calls in all. Every helper is "hosted": it runs on the stream of an idle data family, so helpers add no
   arrays. [read, counted]
2. **Where the saving comes from.** Exactly one place: **each data family gets one of its 73 coordinates for
   free.** That coordinate comes from OpenAI's endpoint-label trick. Each call then pays 64 back for its centre path
   (8 centres × rank 8). Per call that leaves 99.56 − 64 = **35.56 unpaid coordinates out of 7,267.6**, so the
   unit is "99.56 arrays in 7,232 moves at block 73" and saves 0.36 of a move per array (0.49 %). [read; computed]
3. **The other third is grouping.** Children (recursive calls) of rank r cost (r/73)^(1−a), so large children are
   cheaper per coordinate. In T74 the mass-weighted geometric mean child rank is **3.94**: the local word's
   children have rank 1 to 8, and the idle remainder sits in rank-48 children. With rank-1 children only, the same
   deficit would give 1.143e-3. Grouping lifts it to **1.6797e-3**, so **32 % of the saving is grouping**. The hub's
   own diagnosis says the same: a ≈ D/E (autoresearch.md, "revision 3"). [read; computed]
4. **Since #588 the deficit has not moved at all.** In the hub's own checkpoint files, 1 − bill(0) = 4.8924e-3
   exactly for #588, #617, #741, #799, #890, #1042 and T74, and the implied W is exactly 896 × replicas each time.
   The coverage law fixes it: paid mass = 72 × families + 64 × calls. The whole climb of the complex saving, from
   1.6049e-3 (#588) to 1.6798e-3 (the T74 line), +4.7 %, is grouping: the effective child rank rose from 3.44 to
   3.95. The outer reduction added another +0.17 % to κ. The pivot-566
   helper-donor kernel, for example, changes the child counts by +1 at rank 1, −2 at rank 2 and +1 at rank 3 per
   call. That leaves the mass unchanged and raises Σ r ln r. My first-order prediction for its gain is 4.16e-8;
   the hub reports 4.14e-8. [computed]
5. **How it escapes our failure modes.** It never shares a paid move between contents that must later be pulled
   apart. Paid moves act on one stream between gates. Gates are free and act only on wires that sit at one common
   frame, so every gate is "clean". Leftover terms are removed by a second route with the opposite coefficient
   (centres versus side helpers, with halving), never by a zero-divisor factor. Junk (dirty helper content) cancels
   for every value of the junk through the schedule's undo rows. When a hosted helper does share a paid move with
   its host, both need exactly that coordinate factor, so the tie that the pattern identity describes is harmless.
   [argued from the read construction]
6. **What caps it.** The saving is about (u − c/w) / (m · ln(m/r_eff)), with u unpaid coordinates per family
   (u = 1 in every hub word), centre charge c per call over w families per call (c/w = 0.643), block m = 73,
   effective child rank r_eff ≤ R < m/√2 (R = 51, from the outer halving guard). For T74's local word, packing
   alone caps a at 1.6876e-3, +0.47 % (the hub's method; I reproduce its 1.6901e-3 for #1094's word). Local frames
   live in one 9-dimensional window, so local and hosted children have rank ≤ 9. With every local child at rank 9
   and every idle child at 51, the cap is **2.88e-3**. Even the most favourable relaxation of the template (centre
   free, every child at rank 51) stops at **3.98e-2**. [read; computed]
7. **How it was found.** First theory: OpenAI's construction (subspace labels, Alon's design, restoration of
   arbitrary scratch). Then the compiler work that shrank the side roles, to 2^-15. Then a jump from a search:
   Sussman's E8 unit of 120 labels of width 9 (#352 integrates it). Then structure: 7 and 9 windows, hosting,
   removal of every dedicated bank. Since #588 it has been optimisation: owner maps, slots, exact fill, local
   kernels found by search, finer certificates. One theory step in the outer reduction (PR1673) turned κ = a/(1+a)
   into κ ≈ a, worth +0.168 %. [read]
8. **Knobs.** Three knobs exist: more families per centre, children closer to rank 51, and more unpaid coordinates
   per family. The first two are bounded by the caps table in section 3: at block 73 nothing exceeds 4e-2. Only the
   third knob scales, and it is the "fat label" knob. My E8 note records that our own search for fat labels found
   no family that saves. A smaller block helps only if the exchange network stays cheap there. That needs many
   lines per centre, but small label spaces have few lines, and hosting needs idle windows, which needs many
   stages. Those are the walls. [argued]

## Words used here

- **Family / array / role.** A family is one data array that must receive the transform (X or Y bank of a pair). A
  role is a fixed wire of the scalar program; a helper is an auxiliary role whose initial content may be
  arbitrary ("dirty") and must be restored.
- **Block m, coordinates.** The transform acts on m address coordinates of each array (m = 73 here). The plain way
  pays m per array.
- **Frame, label, child.** Each wire carries a frame: a subspace label inside a fixed m-dimensional space. Labels
  only grow along a wire, apart from the centre return. A step of the label by r dimensions is a child: a recursive
  call of rank r. Paid mass = Σ ranks.
- **Bill.** At saving a the bill is Σ_r n_r (r/m)^(1−a) / W, with W the number of families. A unit certifies saving
  a when the bill is ≤ 1 (the hub's model; the glossary and ceiling.py state it). In our words a "W in M" unit
  with rank-1 moves has bill 1 at a = ln(mW/M)/ln m.
- **Deficit D, grouping E.** D = mW − Σ n_r r (unpaid coordinates); E = Σ n_r r ln(m/r). To first order a ≈ D/E.
  The effective child rank r_eff is the mass-weighted geometric mean of r.
- **Coverage law (PR1012).** Paid mass = (m − 1) × families + 64 × calls: every live family keeps exactly one
  unpaid coordinate, and every call pays its centre path.
- **Hosting.** A helper's frame climb is carried out on the stream of a data family that is idle at that time and
  still has those coordinates pending. The helper's work then shares the host's own paid children.
- **Idle remainder.** A family's coordinates that no gate needs at a given time are climbed in large children
  (rank 48 in T74).
- **Halving guard.** The outer recipe requires m² > 2R² for every child rank R (R ≤ 51 at m = 73).

## 1. The object [read]

**Game.** It is OpenAI's linear-network game. A fixed scalar network exchanges two banks, (X_a, Y_a) → (−Y_a, X_a),
with prescribed signs and restoration of every auxiliary value. The scalar program is lifted to arrays: each
scalar gate acts pointwise on whole arrays, and each wire carries a frame. At a gate all wires share one frame, so
the gate commutes with it. The implemented operator is D_out · S · D_in^{-1}: once the scalar routing is right,
only the endpoint frames matter (OpenAI, section 3, "From scalar values to arrays"). Costs are dimension changes
along edges; the paper's budget is Σ|Δdim| = Wm − 2N + 2L.

**The hub's unit (T74, #1839 line).**
- Two local words, `gauge44` and `extra`, dispatched to 249 and 417 calls. Each has 112 ports, 8 centre registers
  with centre frames of rank 8 (cst = 64), R = 480 helper roles (the hub names words by this count, e.g. "R480") and h = 9 stages. The local profiles per call are
  {1:2348, 2:500, 3:416, 4:130, 5:82, 8:8} (mass 5,590) and {1:2340, 2:508, 3:420, 4:126, 5:82, 8:8} (mass
  5,594). [read in the words' `report`]
- The 112 root lines are vectors with four entries ±1 on 14 supports, the weight-4 words of the extended Hamming
  code. These are the 224 E8 roots of that shape, up to sign. The scatter has 112 rows with entries ±1/4.
  [counted in the words]
- Template: K = 4 pairs, 9 stages, m = 73. Live families = 74 × 896 = 66,304; calls = 666. The global profile is
  {1:1,559,100; 2:339,000; 3:277,392; 4:84,912; 5:54,612; 8:5,328; 48:22,748}. [read, SPEC/PROFILE]
- Free steps used: all scalar gates (additions and dyadic scalings over Z[i, 1/2]); binary phase frames with
  Clifford/Z4 phase corrections; address isometries between windows ("bridges"); halving (the 1/4 scatter, the
  1/2 of H0). Helpers are dirty and restored. The hub lists all-size compiler, realisation, hosting,
  precision/tape and analytic transfer as conditional interfaces; they are not proved there. [read]
- The block: by analogy with #386's bit side (nine windows in Q^8 ⊕ (Q^7)^8, so m = 64), the complex block is
  73 = 9 + 8 × 8. [argued: inferred, not read directly]

## 2. Where the saving comes from

**Per call, averaged over T74** [computed, `out/anatomy.txt`]:

| item | coordinates |
|---|---:|
| plain cost, 73 × 99.556 families | 7,267.56 |
| free coordinate, one per family | −99.56 |
| centre path, 8 centres × rank 8 | +64 |
| paid mass (= 72 × 99.556 + 64, checked exactly: 4,816,512 in all) | 7,232.00 |
| of which local-word children, ranks 1–8 (rank 1: 2,341) | 5,592.5 |
| of which idle remainder, 34.16 children of rank 48 | 1,639.5 |
| helpers (480 helper roles per word, about 4.8 per family): extra mass | 0 (hosted) |
| unpaid in the end, D | 35.56 (0.489 %) |

**Move by move, what does the work** [argued, from OpenAI's lemma and the words]:

1. *The free coordinate.* An X family starts at a 1-dimensional label, its own line U_a, at no cost: the input is
   defined as D_in^{-1} f, so no physical step is made. Its Y partner ends at U_a^⊥, of codimension 1. The value
   that enters on wire X_a leaves on wire Y_a. The scalar exchange carries it across wires, and only at gates,
   where frames agree. Every wire's label only climbs, yet the value's route goes "down" from U_a to U_a^⊥ for
   free. For a self-inverse transform, H(U_a^⊥) · H(U_a) is the full transform. This is the −2N of OpenAI's
   budget, which the hub states as "one unpaid coordinate per family".
2. *Doing the exchange exactly costs the centre path.* The exchange runs through gather → centres → scatter.
   Copies of neighbouring sources reach each target through side helpers and cancel the unwanted cross terms. In
   our five-stage E8 unit (120 labels of width 9) a target gets, after the scatter, its own source plus ±1/2 of each of its 56 neighbours, and the
   helpers take those away (my E8 note, section 2). The centres must return from the larger label to the smaller
   one; that return is the only paid decrease. It is OpenAI's 2L, and here it is 8 × 8 = 64 per call.
3. *Helpers cost nothing extra.* The schedule restores arbitrary helper content: subtract the old effect, add,
   then undo. This is the linear form of Buhrman–Cleve–Koucký–Loff–Speelman's transparent computation, cited by
   OpenAI. So a helper may borrow any array, including an idle data family of another call. The host family must
   climb the same coordinates later anyway, so the helper's residual rides on the host's own pending children. The
   hub's idle-axis-hosting package (κ ≈ 1.2376e-3, +9.2 %) and the T21 controller ("removes every dedicated bank")
   did this. The coverage law then holds with families as the only arrays. [read]
4. *Grouping turns the small deficit into a 1.47-times larger saving.* With all-rank-1 children, D gives 1.143e-3
   (−ln bill(0) / ln 73). With r_eff = 3.94 the denominator becomes ln(73/3.94) = 2.92 instead of 4.29, which gives
   1.680e-3. The idle remainder at rank 48 alone carries 1,091,904 of the 4,816,512 mass.

**Why it does not hit our near-miss obstructions** [argued]:

- *Zero-divisor dilemma (Lemma Z).* Our near misses share a paid move between two arrays and must then kill a
  leftover with a factor like 1 + u (u² = 1), which makes the readout non-invertible. The hub has no shared paid
  move that needs separating. Contents are mixed only by free gates at a common frame, and they are separated
  again by exact scalar subtraction. The readout is the signed bank exchange, invertible by construction. Leftovers
  are cancelled additively, by halved coefficients from a second route, never multiplicatively.
- *One route per target.* Each target receives every neighbour twice: through the centres, with coefficient
  (|S∩T| − 1)/2 in OpenAI's version, and through a side helper with the opposite coefficient. The two routes pass
  through different frames: the centre wire's A⊗F returns to B⊗F, while side wires go through P⊗⟨t_X⟩ ⊂ P⊗t_Y^⊥.
  Our moonshot wave-2 note found that the 2-in-7 near miss lacked exactly this second route into each target.
- *Pattern identity (a2-nag).* That identity says junk rows are touched only by moves that also reach the target,
  with tied coefficients. In the hub, junk is dirty helper content, and it never needs a paid move of its own: it
  cancels for every value of the junk. Where a hosted residual and its host do share a paid child, both need that
  same coordinate factor, so the tie is the wanted coefficient. The hub states the window condition: the host is
  untouched between the hosting call and its scheduled advance.
- *Cleanliness at depth 3 (o3).* Every gate's inputs sit at one common label, so there are no mixed-depth
  contents to clean. OpenAI's lemma, "every edge joins comparable labels", is the hub's "literal monotone local
  frames" check.
- *MAP §1.3 conditions.* The hub's unit meets every one of them: one finite word; exact; coefficients in
  Z[i, 1/2]; special position (E8); divides by 2 (scatter 1/4); phases (Clifford/Z4 corrections); multi-term
  contents (centres gather many taps; 28 per centre in #626's words, read in complex-center-quarter); rank
  ≥ W + m − 1, trivially. Its port geometry is the Hamming code built from the Fano plane, with halving over a
  shared base (the centres), which is the pf-duality blueprint. [argued link]

**The price of all this is that the unit is almost plain.** It saves 0.36 coordinate per array out of 73. Our
small-block targets need roughly 15 % to 40 % fewer moves per array (3 in 10 at block 4; 2 in 6 at block 5). The hub's mechanism buys exactness and cleanliness
by spending a whole exchange network to gain one coordinate per family.

## 3. What caps it

No tight certificate in our sense (Theorem R) is known for the hub's units. Every cap below is a relaxation
inside the hub's own model and template: fixed coverage law, halving guard, window. None is a lower bound for
multiplication or for the transform. [read; argued]

1. **The coverage law fixes D.** D/(mw) = (1 − c/w)/m with u = 1. For T74, c/w = 64/99.556 = 0.643: the centre
   path takes 64 % of the raw deficit. The hub says this too. [read; computed]
2. **The local word fixes most of E.** Local children carry 77.3 % of T74's mass, in ranks 1 to 8, and 42 % of the
   local mass is rank 1. By cost share over all children: rank 1 is 32 %, ranks 2 and 3 are 31 %, rank 48 is
   23 %. Child ranks are set by the dimension steps between consecutive gates on a wire, and in this word those
   steps are small. The hub's fixed-word ceiling for #1094's word is a < 1.690054e-3; I reproduce 1.6900541e-3.
   The same ceiling for T74's own local word (idle mass at rank 51) is 1.687577e-3, +0.47 % over T74.
   [read; computed]
3. **The window caps the local ranks.** Every frame in the T74 words is a subspace of a 9-bit space: basis
   vectors up to 256; cleanup at [256, 128, …, 1]. So no local child exceeds rank 9. Hosted helper residuals have
   width ≤ 8 (hub: "eight rank-8 residuals" per bank; "max child 48 → 8" after hosting). Within this window
   architecture, with D fixed: all local mass at rank 8 and idle at 48 gives 2.716e-3; all local at rank 9 and
   idle at 51 gives **2.884e-3**. The hub's "rank-16 local children" lever (×2.31) therefore needs a wider window,
   which means a different port geometry. [computed; argued]
4. **The halving guard caps every child at R = 51.** So ln(m/r_eff) ≥ ln √2 = 0.347, and
   **a ≤ (u − c/w) / (m · 0.347)** (first order). [read guard; derived here]
5. **The pricing itself is concave and mass-neutral.** Any move that keeps the coverage law can only change E. All
   kernels, owner maps and exact fills since #588 are of this kind. [computed]
6. **Search size.** The hub's census (#1264) found no covered open claim that could pass κ ≈ 0.0016873 by packing.
   Its leading claims have lower caps. The hub's ideal lever multipliers are 2.80× (centre free), 1.90× (centre
   shared by two calls), 2.31× (rank-16 local children), 3.77× (two unpaid coordinates per family) and 6.5× (centre
   free with rank-16 children, κ ≈ 1.09e-2). On T74 my multipliers are 2.80, 1.90, 2.30, 3.85 and, with every
   child at rank 51, 8.14 (1.37e-2). [read; computed]

**Caps table for the one-unpaid-coordinate architecture** [derived here, `out/caps.txt`]. Entries are
a = −ln(1 − (1 − c/w)/m) / ln(m/r_eff):

| block m | rank-1 only, c = 0 | r_eff = 4, c/w = 0.643 | r_eff = 4, c = 0 | r_eff = m/√2, c/w = 0.643 | r_eff = m/√2, c = 0 |
|---:|---:|---:|---:|---:|---:|
| 4 | 0.2075 | 0.270 | 0.830 | 0.270 | 0.830 |
| 6 | 0.1018 | 0.151 | 0.450 | 0.177 | 0.526 |
| 8 | 0.0642 | 0.066 | 0.193 | 0.132 | 0.385 |
| 16 | 0.0233 | 0.016 | 0.047 | 0.065 | 0.186 |
| 28 | 0.0109 | 0.007 | 0.019 | 0.037 | 0.105 |
| 45 | 0.0059 | 0.003 | 0.009 | 0.023 | 0.065 |
| 73 | 0.0032 | **0.0017 (T74)** | 0.0047 | 0.014 | 0.040 |

(At m ≤ 5, r_eff = 4 exceeds the guard, so R is used.) Reading: at m = 73 this architecture cannot pass 4e-2
under any packing. The 1e-1 bar needs **m ≤ 28 even with the centre free and every child at the maximal rank**,
and m ≤ 6 or so with today's grouping. These columns are relaxations: they omit the cost of making the exchange
network exist at that block.

## 4. How it was found, and what each step changed

| step (PR, author) | κ or complex a | what changed in D/E terms |
|---|---|---|
| OpenAI #109 framework | κ = 2^-182 | theory: triple-subset labels, h = 100, m = 10^6; D tiny next to the huge side-wire count |
| hub to PR #39 / #62 (compact control, rectangles, shared DAGs, cross-point sharing, compilers) | 5.1e-5 | fewer side roles per unit of deficit (D/(mW) up), by structure and compilers |
| #144 paired cube (icekylinx) | 4.61e-4 | ×9; shared cores (an664) |
| #186, #210, #279 (Dugongue; eumemic; DreamingOfClouds, Rohan Arun), on Grant and Sussman's five-stage construction | 6.62e-4 → 7.24e-4 | banks, reused sources, retiming |
| #352 (DaysSky): "Sussman's E8 complex unit" | (our E8: 8.76e-4 at block 45) | found by **search**: a fat-label search turned up thin families; E8 was recognised; an integer programme chose the shared sums. Copy ratio 0.6, 6.53 helpers per pair |
| #363, #375, #377 (DaysSky; Zhihao Chen): 7 and 9 windows, bridges, halving degree two | ~9.99e-4 (#386 κ) | more stages: m = 73, more families per centre path |
| #442/#434 (DaysSky and lineage): whole-residual and idle-axis hosting (hub's hosting lemma: a dirty invocation acts on any stream as an exact coordinate factor, and idle children are products of commuting coordinate factors) | complex 1.184e-3 → 1.267e-3 | helpers moved onto idle data families. If the bill's W is the "live stock": D/(mW) 2.69e-3 → 4.07e-3 (W 12,720 → 8,400 streams, deficit 2,496), r_eff about 7.5 → 2.9 (max child 48 → 8) [derived here] |
| #531, #535 (icekylinx, DreamingOfClouds) | complex a 1.520e-3, 1.590e-3 | direct delivery, post-centre kernels |
| #588 (Kothe, mixed complex providers ck/jp8/atm3), 8 replicas | complex a 1.604850e-3 | **D/(mW) = 4.8924e-3 (W = 7,168 = 8 × 896), r_eff 3.442**: banks gone, D at its final value [computed] |
| #617 cold helper endpoints (Kothe) | 1.637862e-3 | r_eff 3.661, D unchanged |
| #741, #799 X4 (Boukhalfa; Kothe) | 1.663967e-3, 1.664400e-3 | r_eff 3.837, 3.840 |
| #890 U8/T21 (NewBArch #869/#870, Hao Yan #858 controller), 21 replicas | 1.664724e-3 | r_eff 3.843 (W = 18,816 = 21 × 896) |
| #1012 (thakrarsagar) | — | coverage identity stated (it was already true from #588 on) |
| #1042 (Kothe), #1094 (Rohan Arun), 16 replicas | a 1.674989e-3 → 1.675463e-3 | owner packing only; r_eff 3.913 → about 3.92 |
| PR1297/1319 T74 exact fill (Rohan Arun), PR1625/1672 kernels (mdumitrean, LeverRank) | a 1.679672e-3 → 1.679755e-3 | all idle at rank 48; mass-neutral local reshaping; r_eff 3.945 |
| PR1673 outer (Vojdanovski), #1839 (Kothe) | κ = b − 3.4e-103 = 1.6797550e-3 | **outer theory**: κ ≈ a instead of a/(1+a), +0.168 % |

κ figures and the complex a of #531 to #890 are the hub's [read: checkpoint ASSEMBLY.json files]. The D/(mW) and
r_eff columns are mine [computed, `tools/lineage.py`, `out/lineage.txt`]. W is implied by solving bill = 1 at the
hub's own a, and it came out as an exact multiple of 896 every time, which checks the model. The #1042 and T74 a
values are also mine: 1.674989e-3 (its a/(1+a) = 1.672188e-3 equals the hub's #1042 κ) and 1.6796717955e-3
(equal to the T74 PRICE.json root). The #434 split is [derived here], assuming its "live stock" is the W of the
bill.

**The lineage in one sentence.** Up to #588 the structural work raised D: smaller labels (E8), more windows,
helpers hosted, banks removed. That ended at D/(mW) = (1 − 64/w)/73. From #588 on, everything has
been grouping inside a fixed deficit (+4.7 % in a), plus the outer reduction (+0.17 % in κ).

## 5. Knobs, and what stops each

1. **More families per call at the same centre path** (more ports V or pairs K). This gives ×1.90 at 2× and ×2.35
   at 4×, with a limit of ×2.80 (centre free). It is stopped by the label geometry. E8 gives 112 tetrad lines in
   dimension 8, and the decoder needs about 8 centres of rank 8 (k² for a k-dimensional window). The hub's scoped
   negative results exclude better E8 star copy ratios through 56 deleted lines (e8-first-rank-drop: no ratio
   below 63/111). [read; computed]
2. **Larger children.** This gives ×1.62 at rank 8 and ×2.30 at rank 16, with a limit of ×8.1 (all at rank 51). It
   is stopped by the gates: a family's frame must sit at each gate's common label, so between gates it can only
   step by the gap, which is small. It is stopped as well by the 9-dimensional window: local and hosted children
   have rank ≤ 9, which caps this knob at 2.88e-3. Only idle stretches give rank 48. [computed; argued]
3. **More unpaid coordinates per family** (labels of dimension u ≥ 2, "fat labels"). This gives ×3.8 at u = 2,
   and it is the only knob that scales without bound in u/m. It is stopped because fat labels need neighbour
   structures with good copy ratios. Our E8 note (section 1) records that the overnight search for fat labels "found
   no fat family that saves". The hub lists "changed coverage law" as an open lever. [read]
4. **A smaller block.** By the caps table this is the only way to the bar. It is stopped by three things at once,
   pulling in opposite directions [argued]:
   - small label spaces hold few lines, so w per centre falls and c/w passes 1. In F_2^3, for example, there are
     7 lines (Fano) against about 9 centre coordinates;
   - hosting needs idle windows, so many stages; but more stages make m larger;
   - side corrections grow with the number of neighbours (56 per pair in E8) and must fit in the idle capacity.
5. **Nesting or mixing with other units.** Nesting is already used: children are recursive calls. Mixing is
   already used as well: mixed providers, two local words. ex-laser's Φ barrier says products and powers do not
   beat the best single unit. [read MAP]

## Predictions (written 08:12) against outcomes

- P1, many small credits accumulated rather than one cancellation: **half right.** The deficit is one structural
  credit (one coordinate per family). The accumulation is all in grouping E.
- P2, sharing rather than cancellation: **right in substance.** Hosting is sharing of paid children between a
  host and a dirty helper. The scalar level does use designed cancellation through a second route. The "exact
  copies" detail is wrong: dirt is restored, not copied.
- P3, a model-internal cap and lineage steps ≤ ~1 %: **mostly right.** The hub's own packing ceiling is +0.87 %.
  The steps since #588 are 0.03 % to 2 % each (#617, +2.1 %, is the largest), 4.7 % in all. The earlier
  structural jumps are bigger: ×9 at #144, +9 % from hosting.
- P4, structure-guided search: **right for the recent half; wrong for the start.** The early jumps were
  structural or came from our search (E8).
- P5, no knob to 1e-1 at this block: **right**, with the cap made exact: 4e-2 at m = 73 in the most favourable
  relaxation.

## For an-synth: three scaled constructions this suggests, each with its deciding test

1. **Fat-label hub unit (u = 2) at a nine-window template.** Predicted ×3.8, giving a ≈ 6.5e-3 (still far below
   1e-1). Test: find a family of 2-dimensional nondegenerate labels in an 8- or 9-dimensional binary space whose
   neighbour graph has a scatter with copy ratio below 1, then price it with `ceiling.ceiling(local, 73, 896, 9,
   64, unpaid=2)` in the hub's model.
2. **Few-window, many-line unit at m ≈ 16 to 28.** Predicted up to 0.05 to 0.1, only with c/w ≤ 0.1 and r_eff near
   m/√2. Test: count lines and copy ratio for a candidate geometry (E8 ⊕ E8, or Barnes–Wall in dimension 16) and
   check hosting capacity per call (idle pending mass ≥ helper residual mass) before any word is built.
3. **The hub's mechanism transplanted to our small blocks** (one free coordinate per array, so W in (m − 1)W:
   0.2075 at m = 4, 0.1386 at m = 5, 0.1018 at m = 6, all with rank-1 moves). At m = 4 our certificates already
   exclude it: MAP §0 item 5 lists the block-4 grouped 2,2,2,1 at 0.1675 as tight. At m = 5 and 6 it is not
   excluded by any bound I know of (MAP: at least 2.99 and 2.97 moves per array). Test: an exact bank-exchange
   network at block 5 whose centre path plus non-hosted helper cost stays below 1 move per family, with free
   endpoint labels (source line, sink complement). This is the same question as our open small-block targets,
   now with a named mechanism (endpoint labels plus restored dirty helpers) to search in. [argued] A caution
   [guess]: with one window at m = 5, the 16 odd-weight lines of F_2^5 and E8-style centres (about k centres of
   rank k, so 25), c/w ≈ 25/32 = 0.78 before any helper, and a ≈ 0.028. The mechanism needs a much cheaper
   exchange at small blocks than E8's ratio.

Rank: 3 first. It is the only one that reaches the bar on paper. 1 is the cheapest decisive test. 2 is the
largest unknown.

## Sources

- Hub mirror: `checks/wht43/bound-check/bh-state/hub.git`, commit 0b1cccf8 (2026-10-14
  02:49 −0400). I read: README.md; docs/research/bottleneck-brief.md, research-glossary.md, current-status.md,
  result-history.md, community-reconciliation-third100-review.md (#352, #363, #375, #377 entries);
  autoresearch.md (revision-3 diagnosis and lever table); research/complex-local-word-ceiling/{README.md,
  ceiling.py}; research/idle-axis-hosting/README.md; research/frontier1094-idle-star-ceiling/README.md;
  research/{e8-r480-nine-stage, complex-center-quarter, complex-t21-owners, complex-t16-owner-exchanges,
  complex-t37-owners, e8-first-rank-drop}/README.md.
- #1839: PR body in `bound-check/bh-state/pr1839.json` (created 2026-10-14T11:50:52Z, jeremykothe2); package README
  and PROFILE.json and PRICE.json in `bound-check/rl-price/hub1839/`; the words `gauge44-word-566.json.gz` and
  `extra-word-566.json.gz` and SPEC.json in `bound-check/bh-gap/dl1839/`.
- #386 body: `bound-check/a-e8-seven-stage-word/bodies/pr386.md`.
- OpenAI, "Integer multiplication below n log n" (23 Sep 2026), sections 0 and 3:
  `upstream/openai-math/preprints/Integer-multiplication-below-n-log-n-September-23-2026/build/sections/`.
- Ours: `notes/how-the-e8-unit-was-found.md`; meta/MAP.md §1.1 to
  §1.4; moonshot/NOTE-wave2.md (Lemma Z, one route); meta2/NOTE.md (pattern identity); moonshot/o3/NOTE.md;
  exist2/NOTE.md.
- Hub checkpoint records `certificates/checkpoints/{family588, cold617, frontier741, frontier799, frontier890,
  frontier1042, direct531, post-center535}.json` and the ASSEMBLY.json files they name (complex_profile,
  complex_decimal), extracted to `data/checkpoint-profiles.jsonl`. Also read: research/idle-axis-hosting/PROOF.md
  §1–2 (the hosting lemma) and EXPECTED.json.
- My scripts: `tools/anatomy.py`, `tools/caps.py`, `tools/lineage.py`, `tools/window.py`; outputs
  `out/anatomy.txt`, `out/anatomy.json`, `out/caps.txt`, `out/lineage.txt`, `out/window.txt`.
