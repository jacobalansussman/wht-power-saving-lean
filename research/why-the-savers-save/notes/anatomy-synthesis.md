> Working note written by an AI agent during the round of 15 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# an-synth: one mechanism behind every saver we know, what caps it, and three scaled constructions

Anatomy round, synthesis track, 15 October 2026, 08:35 to 09:03 EDT and, after the orchestrator restart, 09:13 to 09:45 (Mac clock). Input: the three
dissections (`../an-hub/NOTE.md`, `../an-openai/NOTE.md`, `../an-ours/NOTE.md`), read in full, plus the sources
named at the end. I added four small programs of my own (`calc/`, stdlib, `python3 -I`, each under 3 s on the
Mac). No skeptic has read this note. Nothing was posted. No Colab was used. No third-party code was run.

Labels: [computed] = a program here and its output; [read] = taken from a named report or source;
[argued] = my reasoning, checked by nobody; [guess]. Savings are quoted as 1.68e-3 and so on.

## In short

1. **There is one mechanism, not three.** OpenAI's unit, the hub's T74 unit and our Lean-checked words all save in
   exactly one way: *crossing slack*.
   - A data array starts on a small frame (its source label, a line) and its data ends on another array whose
     frame is the complement of that line.
   - The route then receives the full transform although no array's frame ever descends.
   - Exact identity (an-openai, checked there by brute force): paid moves = W·m − 2·(sum of source-label
     dimensions) + 2·(dimension that bridges drop).
   - Our words are literally the hub's object, and OpenAI's is the same mechanism with different parts.
     [read; argued]
2. **One first-order law prices all three.** a ≈ (D/(W·m)) / l̄, where D is the free coordinates net of bridges and
   l̄ is the mean of ln(m/r) over the paid children. With each source's own inputs the law gives:
   - OpenAI 2.654e-13 (reported 2.65e-13);
   - T74 1.677e-3 (1.680e-3);
   - U8X0HQ 1.623e-3 (1.626e-3).

   [computed, `calc/cap.out` (a)]
3. **The bridge (the "centre") is an index code.** The low-rank centre broadcasts t sums. Each target decodes its
   own source using "side information": copies of the junk terms, delivered by carriers or helpers. Those copies
   are cheap only between ports whose binary labels are orthogonal. Two lower bounds follow:
   - the number of centres t is at least the largest set of mutually non-orthogonal ports;
   - the binary window is at least the largest set of mutually junk ports.

   With junk coefficients ±1/2, ports per centre n/t ≤ 1 + (junk per port)/4. E8 meets that bound exactly.
   [computed, `calc/ade.out`, `calc/clique.py`]
4. **Among all two-angle designs, E8 is the best centre there is** [computed for every simply-laced root system;
   argued for all designs, through the Cameron–Goethals–Seidel–Shult classification of lines at 60° and 90°].
   - Every simply-laced root system is a ready-made design: binary labels from B = J + Gram mod 2.
   - Copy ratios (centre cost per port): E8 0.600, E7 0.778 (both exactly the ratios our 9–10 Oct search
     flagged), E6 1.167, every D_n ≥ 1, every A_n ≥ 1.5.
   - In the hub's template only E7 and E8 save at all, and the first saving block is 25 (E8) or 31 (E7).
   - The complex (Gaussian) form of E8 halves the centres, but 12 of its lines are mutually junk. So its window
     is at least 12 and its copy ratio at least 0.8: worse.
   - OpenAI's Alon triples tie E8's best template figure only at ground set 20 to 25, with 13 to 29 times E8's
     junk per port.
5. **How real savers escape the failure modes of our near misses.** Every target gets two singular routes whose
   sum is the identity: the low-rank centre, plus junk carriers that run only between binary-orthogonal ports.
   - Separation is additive, at one common frame, by helpers that are dirty and restored exactly.
   - No paid move ever carries two contents that must be pulled apart later.
   - The price is break-even arrays (OpenAI: about 220,000 per saved move) or narrow mass plus a large block
     (hub: helpers carry 77 % of the narrow mass, and idle cover forces block 73).
   - Our block-3 to block-5 units cannot borrow this. With 2 or 3 arrays the bridge alone costs at least 2m, an
     index code needs about 100 ports before the centre is cheaper than one move per port, and there is no idle
     cover. So crossing slack is not what our near misses lack.

   [argued from the three reports]
6. **What caps the saving in general: three walls.**
   - (i) The index-code floor on the centre (item 3).
   - (ii) Junk cleaning: every junk entry needs a route, paid either as dilution or as narrow mass with idle
     cover.
   - (iii) The stage–window wall: data arrays never descend and gates need one common frame, so every shear
     stage costs each data array k − 1 fresh coordinates of a k-bit window. Hence m ≥ 1 + (k − 1)(2P + 1).

   Over every design I computed, one-sided lines give D/(W·m) ≤ 7.6e-3 even with helpers free, so a ≲ 3e-3.
   With one free coordinate per array, the 1e-1 bar needs m ≤ 6 with rank-1 children, or about 11 with
   r_eff = 4. At m = 25 to 73 it needs 4 to 18 free
   coordinates per array. **No construction in this class reaches 1e-1.** The bar needs a different mechanism.
   [computed tables, argued walls]
7. **Three scaled constructions, ranked** (section 4):
   1. **Two-sided slack**, a second unpaid coordinate per family: the Y bank also starts on a line. Predicted
      ×2.0 over the record (a ≈ 3.4e-3). The range is ×0.2 to ×3.9, depending on how many extra centre returns
      the template needs. A check I ran on our U8 word: no helper or total ever goes from a target back to a
      source inside a call (0 of 616 slots), so the local word adds no cost. The decisive exact test is a
      template-level frame solve.
   2. **A port design beyond E8** (multi-angle tight designs in windows ≤ 12): at most ×1.3 net, because halving
      the centre roughly doubles the junk. A cheap screen decides it; I expect it to fail.
   3. **A crossing network at block ≤ 16 whose stages share window coordinates.** This is the only route in the
      class to the bar on paper (0.096 at m = 16 with free bridges). I expect the stage–window wall to exclude
      it (85 %). An exact label ILP at m = 8 to 12, or a proof of the wall, decides it.

## Words used here

- **Array, block m, paid move, child, W in M, saving**: as in `../../exist2/NOTE.md` and `../an-hub/NOTE.md`.
  - A paid child of rank r is one recursive call on r of the m directions.
  - The saving a solves Σ_children (r/m)^(1−a) = W. First-order law: a ≈ D/E, where D = mW − Σ r and
    E = Σ r ln(m/r).
- **Frame, label, source label, sink label.** Each array carries a subspace of F_2^m (its frame). The frame only
  grows, except on bridge descents. A data route receives the full transform when its sink label is the
  complement of its source label (an-openai 2.2).
- **Crossing slack.** The free coordinates a route gets because its data leaves its starting array. In every
  known saver this is one line per family: u = 1.
- **Bridge / centre / total.** Arrays that carry data from frames containing the source lines to frames avoiding
  them. They must descend. That descent is the centre charge c per call, over w families per call.
- **Junk, carrier, helper.** The low-rank centre delivers unwanted ±1/2 (or ±1/4) multiples of neighbouring
  sources. Carriers (OpenAI, one per junk entry) or helpers (hub, shared partial sums) cancel them. They are
  cheap only between ports whose binary labels are orthogonal.
- **Two-field design.** n ports, each with a real or complex vector f_S and an odd-weight binary label t_S, such
  that A_ST = ⟨f_S, f_T⟩ (normalised, A_SS = 1) vanishes whenever t_S · t_T = 1 over F_2.
  - Centre count: t = rank A.
  - Window: k = the length of the labels.
  - Copy ratio: about t·k/n, the centre charge per port.
- **Idle cover.** The coordinates that a family must climb in stages where its pair is inactive. Hosted helpers
  park their leftover climbs there (an-ours).
- **ρ (dilution).** Data families divided by all arrays: 4.5e-6 for OpenAI, 1 for the hub.

## 1. The common mechanism

**The identity** [read: an-openai 2.1, brute-force checked at h = 7, 8, 9; the hub's coverage law is its form
with u = 1]:

  paid = Σ_arrays (dim sink − dim source) + 2·(descents) = W·m − 2·Σ dim(source labels) + 2·L.

- Gates cost nothing, because every array at a gate shares its frame.
- An array whose data stays home gains nothing: its sink would be the complement of its own source, and the two
  are not nested.
- So **saving needs data to change arrays**, and it is worth exactly the source dimension that is carried
  across, minus what the bridges drop.

**The parameters of the three savers, side by side** [read from the three reports; the last column is mine]:

| | OpenAI (explicit DFT) | hub T74 | our U8X0HQ |
|---|---|---|---|
| block m | 10^6 (= h³, h = 100) | 73 (= 1 + 8·9) | 73 |
| free coordinates per family u | 1 (X bank only) | 1 | 1 |
| ports per invocation, design | 161,700 triples (Alon), t = 101 | 112 E8 root lines, t = 8 | same as hub |
| centre charge c/w | 0.187 | 0.643 | 0.643 (+ 323 banks) |
| junk per port | 147,731 | 52 [computed: E8 minus 8 mutually orthogonal lines; tight, 112/8 = 14 = 1 + 52/4] | same |
| junk cleaning | one private carrier per junk entry | about 4.5 shared-sum helpers per pair, hosted on idle families | same, same-run hosting |
| dilution ρ | 4.5e-6 | 1 | 0.957 |
| pricing | per dimension (l̄ = ln m = 13.8) | per edge, r_eff 3.94 (l̄ 2.92) | per edge (l̄ 2.885) |
| saving | 2.65e-13 | 1.6797550e-3 | 1.6256629e-3 (Lean) |
| first-order law [computed] | 2.654e-13 | 1.677e-3 | 1.623e-3 |

**Where they differ** [argued]:
- *Design.* OpenAI buys a cheap centre (1,601 ports per centre) with enormous junk. E8 has 15 ports per centre
  and 56 junk terms per port.
- *Junk cleaning.* OpenAI uses one array per junk entry, which makes 99.9995 % of its arrays break-even
  carriers. The hub uses shared sums chosen by search and integer programmes, at 0.087 helpers per junk entry
  (504 helpers per call for 112 × 52 = 5,824 junk entries).
- *Where the helpers' cost lands.* For OpenAI it lands in W, as dilution. For the hub it lands in narrow mass,
  with no extra rank (an-ours' identity: local rank − leftover widths = 224·8 + 64 on every word).
- *Pricing.* Per dimension versus per edge: ×80 to ×360 on OpenAI's own network (an-openai).

So the community's whole climb, from OpenAI's 2^-182 to 1.68e-3, came from three changes:
- a better design (E8);
- cheaper junk cleaning: shared sums and hosting, which turned dilution into narrow mass;
- a better price (per edge, then grouping).

The deficit itself has not moved since #588 (an-hub's lineage computation). [read]

## 2. How real savers escape our failure modes, and why our small units cannot borrow the escape

| our failure mode | what all three savers do | why it is out of reach for W = 2 or 3 at blocks 3 to 5 |
|---|---|---|
| zero-divisor dilemma (absorb4: one route must be invertible and also kill leftovers) | the target is the sum of two singular routes: the low-rank centre RG plus the junk routes JV, with RG + JV = I. Neither route has to be invertible | two routes into every target need a bridge; with 2 arrays any bridge already costs ≥ 2m (an-openai 2.3) |
| one route per target | a target gets the centre and many carriers or helpers (36 to 97 terms in U8; 101 + 147,731 in OpenAI) | there are not enough arrays to host the second route |
| pattern identity (a2-nag: junk rows are touched only by moves that also reach the target) | junk does ride with the target through the centre, but every junk term also has a route of its own that never touches the target's data | the same |
| cleanliness at depth 3 (o3) | no mixed-depth contents: every gate joins arrays at one common frame, and dirty scratch is restored exactly (Bennett; Buhrman–Cleve–Koucký–Loff–Speelman; the hub's hosting lemma). In our U8 word every one of the 616 helpers and totals visits sources before targets and never goes back [computed, RUNLOG 08:58] | the restoration needs idle arrays or idle coordinates |
| MAP §1.3 (halving, phases, address maps, multi-term contents, special position) | all three use all of them: E8 and Alon's two-field identities are the special position | not a separating condition |

**The price of the escape is the same everywhere: many arrays per saved coordinate.** In OpenAI they are
break-even carriers; in the hub they are a large block and hosted helpers. A saving of this kind needs, at the
very least:
- tens of ports per centre (an index code with c/w < 1 needs n/t > k·(stage factor); the smallest that works is
  E7 with 56 ports);
- an idle stock to host the junk cleaning.

This is the plain reason why no crossing-slack unit exists at our small blocks [argued]. What our near misses
lack is not this mechanism. If a 1e-1 unit exists, it uses something else.

## 3. What caps saving in general

**Wall 1: the index-code floor on the centre** [computed for root systems; argued in general].
- *Floors.*
  - t ≥ (largest set of ports with mutually non-orthogonal labels), because A is diagonal on such a set.
  - The binary window k ≥ (largest set of mutually junk ports), because their labels must be pairwise
    orthogonal with odd weight.
  - For a tight design, n/t = 1 + Σ_T |A_ST|², which is at most 1 + s/4 when the junk entries have size ≤ 1/2.
- *E8 meets the centre floor exactly* (t = 8, and E8 has 8 mutually orthogonal root lines). Its window 9 is one
  above the junk-clique floor 8. It is tight, and so is its 112-port form (112/8 = 14 = 1 + 52/4).
- The table (`calc/ade.out`, copy ratio = t·k/n):

| design | ports n | centres t | junk per port | n/t | window k | copy ratio |
|---|---:|---:|---:|---:|---:|---:|
| A_n (n ≤ 16) | 6 to 136 | n | 4 to 30 | 2 to 8.5 | n or n+1 | 1.5 to 2.0 |
| D_n (n ≤ 17) | 12 to 272 | n | 8 to 60 | 3 to 16 | n or n−1 | 1.00 to 1.25 |
| E6 | 36 | 6 | 20 | 6 | 7 | 1.167 |
| E7 | 63 | 7 | 32 | 9 | 7 | **0.778** |
| E8 | 120 | 8 | 56 | 15 | 9 | **0.600** |
| E8 as 60 Gaussian lines in C^4 | 60 | 4 | 44 | 15 | ≥ 12 (junk clique) | ≥ 0.8 |

- By the classification of line systems at 60° and 90° (Cameron, Goethals, Seidel and Shult, 1976; cited from
  memory, not re-read), every real design with junk ±1/2 is a subset of one of these. The hub's port-deletion
  study (no ratio below 63/111) covers E8's subsets [read]. So **E8 is the best centre among two-angle designs**,
  and the community's and our searches had already found it.

**Wall 2: junk.** The trade-off is n/t ≤ 1 + s·max|A|², so a cheaper centre means more junk per port.
- Junk costs either dilution (one break-even array per junk entry: OpenAI) or narrow mass with idle cover
  (shared-sum helpers: the hub, 77 % of the narrow mass).
- OpenAI's Alon triples placed in the hub's template show the trade exactly [computed, `calc/cap.out` (b)]:

| design | copy ratio | junk per port | best P | m | D/(W·m), helpers ignored |
|---|---:|---:|---:|---:|---:|
| E8 (P112) | 0.571 | 56 | 2 (4 in T74) | 41 (73) | 6.97e-3 (4.89e-3) |
| E7 | 0.750 | 32 | 3 | 43 | 2.91e-3 |
| Alon triples h = 16 | 0.486 | 325 | 1 | 49 | 5.54e-3 |
| Alon triples h = 20 | 0.368 | 731 | 1 | 61 | 7.33e-3 |
| Alon triples h = 25 | 0.283 | 1,606 | 1 | 76 | 7.58e-3 |

  No design beats E8 by more than 9 % in D/(W·m) even with helpers free, and the ones that tie need 13 to 29
  times the junk.

**Wall 3: stage and window** [argued; the template geometry is read from an-ours and pr717].
- Data arrays never descend, and a gate needs one common frame. So in every shear stage each participating data
  array climbs k − 1 coordinates of that stage's k-bit window, and the next stage needs fresh coordinates.
- A signed swap needs at least three shears in all (two per pair, interleaved, in the hub's template), so
  m ≥ 1 + (k − 1)(2P + 1).
- Any saving two-angle design needs k ≥ 7 (E7) or 9 (E8). That gives m ≥ 25.

**Together.** Write a ≈ ρ·(2·avg source dim − c/w)/(m·l̄).
- With one-sided lines (u = 1) the computed ceiling over these designs is D/(W·m) ≤ 7.6e-3. With l̄ ≥ ln(m/9),
  that gives a ≲ 3e-3.
- With two-sided lines (u = 2), about twice that.
- The 1e-1 bar needs these free coordinates per array [computed, `calc/cap.out` (d)]:

| m | r_eff = 1 | r_eff = 4 | r_eff = m/√2 |
|---:|---:|---:|---:|
| 8 | 1.50 | 0.54 | 0.27 |
| 16 | 3.87 | 2.07 | 0.55 |
| 25 | 6.88 | 4.19 | 0.85 |
| 73 | 25.5 | 18.4 | 2.49 |

T74 has 0.357 at m = 73. **The crossing class cannot reach the bar at any block its walls allow.** The bar needs
a mechanism that frees a fifth to a quarter of each array's coordinates, which no known object does [argued].

## 4. Three scaled constructions, ranked

### 4.1 First: two-sided slack, a second unpaid coordinate per family

**Construction.** Keep the hub's template and the U8 (or T74) local word.
- Give every Y family a source line v_a as well. The v-lines are a mirrored copy of the port lines (same binary
  Gram, coordinates disjoint from the u-window: one extra shared coordinate, so m = 74).
- Keep every X family inside v_a^⊥. Then X: U_a → V_a^⊥ and Y: V_a → U_a^⊥, so each family pays m − 2 instead of
  m − 1.
- Both routes still receive the full transform, because Φ(V^⊥)Φ(V)^-1 = C^⊗m·R_v is a free shift, as for u.
- This is not "fat labels": both banks keep line labels. It is the hub's open card "Coverage law: explain the
  source of the single unpaid coordinate per family; stop if it is forced by the data endpoint contract"
  (research-priority-log, read). It is not forced by the endpoint contract [argued].

**Predicted saving** [computed from T74's per-call numbers, `calc/cap.out` (c)]:

| what the template needs | D per call | a | × record |
|---|---:|---:|---:|
| no extra centre return (the hub's ideal lever, ×3.77 in its own table) | 135.1 | 6.5e-3 | 3.88 |
| one extra rank-8 return per centre per call (the v-span) | 71.1 | **3.4e-3** | **2.03** |
| two extra returns (a dirty target–source–target–source pattern, as in OpenAI's restoration) | 7.1 | 3.4e-4 | 0.20 |

- **What I already checked** [computed, `calc/twosided.py`, RUNLOG 08:58]: inside the U8 local word the second
  line costs nothing.
  - All 616 helpers and totals meet sources before targets and never go back.
  - No gate joins a source and a target whose ports are not orthogonal.
  - So with mirrored v-lines no slot needs a v-drop inside a call.
- Any extra cost therefore sits at the template level: how dirty helpers are restored across calls, and how the
  stage in which Y's data flows to X drops the v-span.
- My central prediction is ×2.0 (one extra return per call), with ×3.9 possible if hosting restores the helpers
  without a target-to-source visit [argued].

**Decisive exact test (for an-build).**
1. *Calibrate.*
   - Rebuild the template frame schedule of the U8 line with our own code (`an-ours/identity.py` already closes
     the per-run identity: 72 paid per data array, 64 per call, helpers neutral).
   - Reproduce D = 2560 and the exact-fill price 1.6647327e-3 (`an-ours/knobs.out`).
2. *Scaled version.* Add the constraints:
   - Y_a ⊇ ⟨v_a⟩ from the start;
   - X_a ⊆ v_a^⊥ throughout;
   - in the stages where Y is the source, the v-lines are its ports.

   Solve the frame assignment exactly for the fixed gate order (an ILP of the #698 / pr717 P5 kind: frames are
   optimal for a fixed gate order). Read off the exact extra paid rank per call.
3. *Price* with our fold root.

The outcome decides the hub's card either way:
- **Kill:** the extra charge is ≥ 2 × 64 per call.
- **Pass:** it is ≤ 64 per call, which gives D ≥ 2(w − 64) and a ≥ 3.3e-3.

Only after a pass would the scalar network be built and checked exactly in Lean (heavy Lean on Colab, as
always).

**Why first.** It is the largest realistic gain in the class (×2 to ×3.9 over the record), it has a cheap exact
test on our own word, and it answers an open question of the community. It cannot reach 1e-1: at best 6.5e-3.

### 4.1a Addendum from the resumed run (09:15 to 09:40): where the second line comes from, and a sharper test

I re-checked proposal 1 against the hub's card and our template notes. The ranking and the numbers above stand.
Three refinements:

- **The second line need not be new** [read: an-ours; argued]. In the template, every call already has its sources
  on their port lines and its targets ending on their port hyperplanes, inside that call's window (8 levels per
  active stage; stage order PAIR = [1,1,1,2,2,3,3,4,4] in `../an-ours/identity.py`).
  - So each call has 224 in-window slack coordinates, and 2,016 per replica. The coverage identity keeps 896 of
    them, one per data array.
  - Two-sided slack keeps one more per family: the Y family's port line *in its own source stage* goes into its
    start frame, and the X family's end frame then avoids it.
  - No extra coordinate is needed, so m stays at 73, not the 74 I assumed above. The gain from that is only
    about 2 % (74/73 × ln(74/3.94)/ln(73/3.94) = 1.018).
  - One condition has to be checked family by family: the two lines must be orthogonal over F_2. Two lines in
    different windows can fail this only through a coordinate the windows share.
- **The decisive test, made exact** [argued]. Instead of an ILP, use the frame-legality propagation that
  pr717's README states (line 361; read as text, not run):
  - G(r, t) is the join of start frames along time-respecting gate journeys.
  - U(r, t) is the meet of the locked end frames reachable forward.
  - A frame assignment is legal if and only if G ⊆ U.

  Run it on one whole replica with the two-sided endpoint frames, allowing the centres their present per-call
  descents. Every place where G ⊄ U needs a descent, at 2 moves per dimension, so the sum is the exact extra
  charge per call at a fixed gate order.
  - Our gcert/1 tools hold frames for one call only. Hosting across calls exists only in the pricing tools
    (`x-hosting/tools/host_price.py`). So an-build must first write the replica frame model and check that it
    reproduces D = 320 per replica (the calibration step).
- **Where a conflict can arise**, which is what the test must look at. Inside one call there is none (computed
  above). Across calls there are three places:
  1. a hosted helper that has touched y_S carries v_S onto its host, whose frame must stay orthogonal to the
     host's own lines;
  2. a centre carries the v-span into the next call;
  3. pair 1 is active in three consecutive stages, so its X sources carry their free line through a stage that
     is not the line's own.

  Point 3 already arises in the one-sided template. Yet the identity says every call pays one return (64), not
  two. Whatever lets pair 1's third stage avoid a second return is exactly what would make a second free line
  cheap. Finding it in the replica frames is the first thing to read off [argued].
- **An upper variant** [guess]. Keep every active stage's line, so 2,016 free coordinates per replica.
  - If each carried line cost one extra 64 in each call whose sources carry it, my crude count gives
    D ≈ 1,312 per replica against 320, about ×4.
  - But that same crude model also predicts the double return of point 3, which the identity rules out. So I
    do not trust this number. It only marks how much the coverage-law question is worth.

### 4.2 Second: a port design beyond E8

**Construction.** Replace the 112 E8 ports by a tight design with inner products in {0, ±1/4, ±1/2} and a window
≤ 12. Candidates [guess]:
- the minimal vectors of the Barnes–Wall lattice BW16 (2,160 lines in R^16);
- the Coxeter–Todd lattice;
- sub-designs of the Leech lattice's minimal vectors (98,280 lines, n/t = 4,095, but 51,704 junk terms per port).

**Predicted saving.** Wall 2 says that halving the copy ratio roughly doubles the junk per port.
- With cst halved, D per call rises ×1.9 (to 67.6). If the helpers' narrow mass grows ×1.5 to ×2, E rises ×1.4
  to ×1.8.
- Net ×1.05 to ×1.35: a ≈ 1.8e-3 to 2.3e-3. The free-centre ceiling is ×2.8 (4.7e-3).

I expect every candidate to fail the screen (70 %) [guess].

**Decisive exact test.** A screen in seconds per design, by generalising `calc/ade.py` and `calc/clique.py`:
1. Compute n/t and tightness, the junk clique ω, the window rank_F2(B), and the copy ratio t·k/n against the
   bound ω/(1 + Σ|A|²).
2. Compute the junk per port.

A design passes if:
- its copy ratio is ≤ 0.35;
- its window is ≤ 12;
- its junk per port is ≤ 120 (about twice E8's).

Only a passing design goes on to a shared-sum helper ILP (`tools/e8/xmilp.py`) and to pricing.

**Why second.** The screen is cheap and exact. The upside is small, and the floor argument suggests that E8 is
already optimal.

### 4.3 Third: a small-block crossing network whose stages share window coordinates

**Construction.** At block m ≤ 16, use two-sided lines on all data arrays, a centre of t shared arrays and junk
carriers. The shears re-use window coordinates instead of taking fresh ones.

**Predicted saving if it existed.** With u = 2 and free bridges at m = 16, r_eff = 4: a = 0.096. This is the only
route to the bar inside the crossing class [computed, `calc/cap.out` (d)].

**Decisive exact test.** An exact label ILP:
- Variables: subspaces of F_2^m for every array at every gate, at m = 8, 10, 12.
- Fixed network: a three-shear signed swap of n ≤ 16 pairs through t shared centres plus one carrier per junk
  entry.
- Objective: minimise Σ|Δdim|, with monotone frames except on centres and one common frame per gate.
- Calibration: the same ILP must give exactly 2m per pair with no bridge (an-openai 2.3) and reproduce OpenAI's
  label table at h = 7 (`an-openai/calc/labels_check.py`, whose sum is W·m − 2N + 2L).

Alternatively, prove the stage–window lemma ("each shear stage costs every participating data array k − 1 fresh
coordinates").

I predict no network with D > 0 at m ≤ 16 (85 %) [guess].

**Why third.** It is the only proposal that touches the bar, but its predicted outcome is negative. Its value is
to close the crossing class, which would make "find a different mechanism" the only open direction for 1e-1.

## 5. Predictions written before reading (STATUS.txt), and outcomes

| prediction | outcome |
|---|---|
| P1: one mechanism (endpoint labels minus bridge); hub and ours the same object; OpenAI the same mechanism at a worse ratio | **held.** Its form is crossing slack with an index-code centre; OpenAI's worse ratio is dilution, not the centre (its c/w of 0.187 is better than E8's 0.643) |
| P2: a cap of the form u_eff/(m·ln(m/r_eff)); with u_eff ≤ 1 the bar needs m ≤ 6 to 8 | **held.** The first-order law reproduces all three savers; the bar needs m ≤ 6 (rank-1 children) to about 11 (r_eff = 4) |
| P3: rank first a small-block exchange; refuted at m ≤ 6; best realistic ≤ 1e-2 | **half held.** After the analysis I rank it third (stage–window wall); the best realistic proposal is 3.4e-3, under 1e-2 |
| P4: fat labels second; the bridge charge scales with u, so u_eff stays ≤ 1 | **held in substance, wrong in form.** The second coordinate comes from two-sided slack, not fat labels. Its charge does scale (one more return), giving ×2 rather than the hub's ideal ×3.8; at the hub's c/w, u_eff = 0.71 |
| P5: mixing will not help | **not tested.** Nothing in the three reports contradicts it (Φ barrier; τ-transfer) |

## 6. Files and checks

All in `checks/wht43/anatomy/an-synth/`:
- `calc/ade.py` → `calc/ade.json`, `calc/ade.out`. Two-field designs from A_n, D_n, E6, E7, E8 and the 60
  Gaussian lines of E8: n, t, junk, exact tightness (A² = (n/t)A), window rank_F2(J + G), copy ratio.
- `calc/clique.py` → `calc/clique.out`. Exact maximum cliques of the junk graphs: E8 8, Gaussian E8 12. Also, by an inline check in
  RUNLOG, 8 mutually orthogonal E8 root lines.
- `calc/cap.py` → `calc/cap.out`, `calc/cap.json`. (a) the law calibrated on the three savers; (b) design ×
  template table; (c) two-sided slack on T74; (d) the free coordinates per array that the bar needs.
- `calc/twosided.py` → `calc/twosided.out`, `calc/twosided.json`. Our U8 lowering under mirrored v-lines: 0
  conflicts. The inline follow-up in RUNLOG: 0 of 616 slots meet a source after a target.
- `STATUS.txt` (dated predictions, progress) and `RUNLOG.txt`.

## Sources

- The three dissections: `../an-hub/NOTE.md` (and its `out/`), `../an-openai/NOTE.md` (and `calc/`),
  `../an-ours/NOTE.md` (and `identity.out`, `knobs.out`, `stages.out`, `u8word.out`).
- `../../meta/MAP.md` sections 0 to 1.4; `../../moonshot/NOTE-wave2.md` (opening: two-layer bound, Lemma Z, the
  second route); `../../meta2/NOTE.md` (v-nag, nag-rerun: pattern identity); `../../moonshot/o3/NOTE.md`
  (opening); `../../exist2/NOTE.md` (opening).
- `notes/how-the-e8-unit-was-found.md` sections 1 to 3 (the
  flagged copy ratios 0.600, 0.778, 0.781), and `notes/open-directions-3.md` section 7 (the copy floor as a
  minrank problem).
- Hub mirror `bound-check/bh-state/hub.git` at 0b1cccf8 (14 Oct, 02:49 EDT), read as text: `autoresearch.md`
  (lever table: two unpaid coordinates ×3.77); `docs/research/research-priority-log.md` (the "Coverage law"
  card); `research/bit-revision3-levers-seed03/README.md` (opening).
- `bound-check/g-p112-decoder/pr717/README.md` P4 (seven stages +0.2 % only) and the claims table.
- Our U8 lowering `bound-check/k-lower/out/U8X-parsed.json.gz`; its format from our own
  `tools/gx/gxcore.py` of this repository.
- Cameron, Goethals, Seidel and Shult, "Line graphs, root systems, and elliptic geometry", J. Algebra 43 (1976):
  cited from memory for the classification of real line systems at 60° and 90°; not re-read.
