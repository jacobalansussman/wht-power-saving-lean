> Working note written by an AI agent during the round of 15 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# cx-mirror: the mirrored second line passes its counts and fails on cleanliness

Anatomy round 3, track cx-mirror, 15 October 2026, 12:24 to 13:31 EDT (Mac clock). Task (BRIEF3.md): give each
family its own second line v_b = s + g_b, count each line's crossings per call with an2-word's formulas, and only if a
gain is left, design the gather and a mirrored helper circuit, verify it exactly and price it.

Labels: **[computed]** a program in `calc/` (or `job-twistall/`) and its output; **[proved]** a complete argument
written here, checked by me only; **[argued]** reasoning with a named gap; **[read]** from a named file. All programs
are our own (stdlib Python). Mac runs went through a private memory-capping wrapper; one exhaustive job ran on a Colab. Nothing was posted.
No third-party code was run.

## In short

1. **The lines exist, and the natural choice is very good.** Take g_b = j + u_b, where j is the all-ones vector of
   the 9-coordinate window, so g_b is the complement of the port. Then [computed]:
   - every g_b has even weight, and every family's two lines are orthogonal (u_b · v_b = 0);
   - the second lines have exactly the same binary Gram matrix as the ports (v_a · v_b = u_a · u_b), so the second
     centre can be a second copy of the hub's 112-port code, at 64 per call;
   - the round-2 orthogonality test (g_a · g_b = 1 on every junk pair) passes automatically.

   The same holds for every "twisted mirror" v_b = s + j + A·u_b with A an orthogonal map of the window. Every such A
   turns out to be a symmetric involution [proved, with `qcheck.py`]; there are 69,616 of them [computed].
2. **The counts the brief asked for leave a gain.** With our template's per-family gluing, every call starts with
   x_b on its port line and y_b on its own second line. So each line needs exactly 1 crossing per call, in both
   directions. That is a floor of 9 per replica against a gain of 896 [computed]. Read with s as one fixed global
   coordinate instead, the vector s itself needs 112 crossings in 5 of the 9 calls, a floor of 564, which still does
   not reach 896 on its own [computed].
3. **Two of the round-2 requirements are wrong as worded.**
   - R2, the "common scatter frame", cannot exist [proved, computed]. A y array never leaves its own u^⊥. So a
     frame shared by several y's lies inside the meet of their u^⊥. For the 56 readers of one U8 total that meet is
     2-dimensional, and with the junk condition it forces all 56 lines to equal s: the shared line, which round 2
     closed.
   - The forward-call formula counts only climbs. A forward call also needs at least one drop across each family's
     second line [proved]. This changes nothing in our template, because one frame word serves both directions.
4. **What the second line may cost.** If it cost only a second P112 centre (64 per call), the price would be
   3.2824847e-3, ×1.954 the hub [computed, calibrated on 1.6647327e-3]. To beat the hub at all, everything beyond the
   two centres must stay under **34 levels per call** [computed].
5. **What kills it: the helpers can no longer be cleaned.**
   - In U8, every x ends on the same frame, the whole window. One end gate can therefore clean every helper.
   - With family-specific lines, each x_b ends on its own hyperplane v_b^⊥. A register that has met x_a can meet
     x_b later, without a paid drop in between, only if g_b · u_a = 0 [proved].
   - For the natural mirror that fails exactly on the junk pairs. So the pairs whose side information makes the
     index code cheap are exactly the pairs that no register can visit one after the other.
   - Inside each of U8's 8 totals (56 x's, 784 junk pairs), at most **4** x's can be visited in a row [computed,
     exact]. For a twisted mirror the proved bound is 16, and the exhaustive run over all 69,616 twists finds at most 8
     [computed].
6. **Priced at first order, the cleaning alone is many times the room.** Each total's content has to be put
   together from its 56 x's and later taken apart against them, a few x's at a time, with a paid drop between runs.
   - Natural mirror: 462 levels per call. Best twisted mirror: 272 per call. The room is 34 [computed, under a
     deliberately optimistic model].
   - Counting only the number of runs, at one level per drop, the natural mirror still needs at least 208 per call
     (52 if one register served all four totals of each x). A twist needs at least 96 (24 with that sharing). So for
     the twists the verdict rests on the drop sizes, about 2.8 levels each, not on the run count alone [argued].
   - Even if one register could serve all four totals of each x, the cost stays above the room [argued].
7. **U8's own gates cannot carry it either.** At U8's gate order, 2,453 of its 7,948 gates are frame-illegal under the
   natural mirror, and 2,360 under a twisted one. Our one-sided calibration of the same test gives 0 [computed, exact].
8. **Verdict: excluded at first order for every second line that keeps the second centre a copy of P112.** That
   covers the natural mirror and all 69,616 twisted mirrors. For the natural mirror the run count alone excludes it.
   For the twists the drop sizes are needed as well. No word was built, so there was nothing to verify on raw arrays
   and no real price to state. What is not covered: second lines outside the twisted family (my search found none, but
   it was too short to settle it), and centres that are not organised as one register per total.

## 1. The setting, made exact

**The call contract** [read: gxcore.py chk_replay; an-build calib.py; our HDef.lean gluing]. In the one-sided
template, each call starts with x_b on its port line ⟨u_b⟩ and y_b at frame 0. It ends with x_b on the whole window
and y_b on u_b^⊥. The gluing x01..y78 of our Lean layout routes each family into a group of the next stage with its
own frame map. So the identification between calls is per family.

The two-sided version adds a local coordinate s (block 74):
- x_b runs from ⟨u_b⟩ to v_b^⊥;
- y_b runs from ⟨v_b⟩ to u_b^⊥;
- at the end of a call, x_b sits on C_b + ⟨u_b⟩ and y_b on C_b + ⟨v_b⟩, with C_b = u_b^⊥ ∩ v_b^⊥.

The next call re-identifies C_b with its background and the two lines with its own local lines. This is the same
form of equation as the one-sided gluing, with two re-identified lines instead of one [argued: I did not write the
Lean]. Each data array pays 72 of 74, so D per replica = 1792 − 9C, where C is the whole centre charge per call.

**The ports** [computed, `lines.py`]. U8's 112 ports are odd vectors of weights 3 and 7. Its centre consists of 8
totals, each the sum of 56 x's. Each y reads 4 totals, and the delivered map has diagonal 1. There are 5,824 junk
entries, 52 per target. **The junk pairs are exactly the binary-orthogonal pairs**, and all 59 other partners of a
port are non-orthogonal.

## 2. The second lines

**Natural mirror** [computed]. g_b = j + u_b gives:
- v_b = s + j + u_b, with v_b · v_b = 1;
- u_b · v_b = 0;
- v_a · v_b = u_a · u_b for all pairs;
- the cross terms u_a · v_b = 1 + u_a · u_b.

The map ψ(w) = w + wt(w)(s + j) on the window, with ψ(s) = j, is an isometric involution of F₂¹⁰. It swaps u_b and v_b
for every family [proved: a three-line check; also checked by `psi.py`]. So the two line systems are exact mirror images inside one extra
coordinate.

**All mirrors that keep the second centre a P112 copy** [proved, computed: `twist.py`, `perms.py`, `qcheck.py`].
v_b = s + j + A·u_b keeps the v-Gram equal to the u-Gram for every orthogonal A. The remaining condition is
u_b · A·u_b = 1 on all 112 ports.
- An orthogonal A fixes j.
- So the condition reads ⟨e, A·e⟩ = 0 on the 112 vectors e = u_b + j of the even space E.
- The quadratic forms that vanish on the ports form a 9-dimensional space: the 8-dimensional space of forms that
  vanish on every odd vector, plus one form that is 1 on odd weights 1, 5 and 9.
- That last form is 1 at j, and ⟨0, A·0⟩ = 0, so it cannot occur. Hence ⟨e, A·e⟩ = 0 on all of E.
- So A + Aᵀ has rank at most 1 and is alternating, hence zero: **A is a symmetric involution.**
- These are A = I + M, with M symmetric, of zero diagonal, M² = 0 and M·j = 0. The image of M is a totally isotropic
  subspace of E. Counting: 1 + 5,355 + 2,295 × 28 = **69,616 twisted mirrors**.
- Of the 9! coordinate permutations, only the identity passes. 200,000 random products of Eichler maps and
  permutations gave no asymmetric A, as the proof says.

## 3. The counts the brief asked for

an2-word's formulas, for every vector of F₂¹⁰ and every call of the schedule (GF at stages 0, 2, 4, 6, 8; GB at 1, 3,
5, 7) [computed, `counts.py`]:

| reading of the chain between calls | largest count per line, per call | floor per replica | gain |
|---|---|---:|---:|
| per-family gluing (our Lean form) | 1 in every call | 9 | 896 |
| s one fixed global coordinate | 1 in each pair's first call; 112 on s in the other 5 calls | 564 | 896 |

Neither reading kills the lead on its own, so by the brief's rule I went on to the design.

**A correction to the formulas** [proved]. Take, for one family t, the registers whose frame is not inside v_t^⊥.
- y_t is one of them from the start.
- At the end, their contents span {x_c for c ≠ t, y_a + x_a for all a}, which does not contain y_t.
- Gates and climbs never shrink the span of what has ever been on that side. Only a drop out of it can lower it, and
  by at most one.
- So **a forward call needs at least one drop across each v_t**, not only climbs. By the mirror ψ, a backward call
  likewise needs drops across the u-lines.

So every call needs drops of both line systems. In our template one frame word serves both readings, so its drops
count once.

## 4. R2 and the scatter

[proved, computed: `lines.py`, `design.py`]
- A y array's frame always lies inside its own u^⊥: it ends there and never drops.
- So a frame shared by a set of y's lies inside the meet of their u^⊥.
- The 56 readers of each total have ports spanning 8 dimensions, so the meet is ⟨s, w_k⟩, of dimension 2. Over all
  112 y's it is ⟨s⟩.
- A shared frame that also holds the readers' lines forces g_t ∈ {0, w_k}. The junk pairs among the readers then need
  w_k · w_k = 1, which no even g allows. So **every reader's line is s: the shared-s case that round 2 closed.**
- For the natural mirror, at most 8 y's can share a frame that holds all their lines (a non-orthogonal clique).

R2's sufficient condition, "contains none of the readers' port lines", is necessary but far from enough. A working
scatter needs copies that climb from reader to reader. Its cost is the second centre: the merged copy holds all 56
readers' lines and drops 8 levels per total, 64 per call [computed]. That is the "gather" of R3 in the backward
reading.

## 5. Room in the price

[computed, `price.py`: the Lean fold hypothesis, the U8 profile ×72, exact fill at rank 48, calibrated on
1.6647327e-3]

| case | C per call | D per replica | price |
|---|---:|---:|---:|
| one-sided calibration, block 73 | 64 | 320 | 1.6647327e-3 |
| second P112 centre and nothing else, block 74 | 128 | 640 | 3.2824847e-3 (×1.954 the hub) |
| largest C that still beats our exact fill | 163 | 325 | |
| largest C that still beats the hub (1.6797550e-3) | 162 | 334 | |

So everything beyond the two centres (extra drops, helpers that cannot be hosted) must stay under 34 levels per call.
That is about 4 levels per total. The same holds whether the extra pieces are priced as centre pieces or absorbed into
the fill.

## 6. Why the cleaning breaks

**The visiting lemma** [proved]. Between two of its drops, a register's frame only grows, and each gate puts it on
the frame of the arrays it joins. So:
- if it joins x_a and later x_b, then u_a lies in X_b's frame, which lies inside v_b^⊥. Hence N[b][a] := g_b · u_a = 0;
- if it joins y_t and later y_t′, then N[t][t′] = 0;
- if it joins x_b and later y_t, then u_b · u_t = 0;
- if it joins y_t and later x_b, then g_t · g_b = 1.

For the natural mirror, N[a][b] = 1 exactly on junk pairs. For every twisted mirror, N is symmetric. So a drop-free
"run" of x-visits is a set with no N-edges.

**Run lengths** [computed, `design.py`, `twist3.py`, `job-twistall`; proved bounds]:
- Natural mirror: inside each total's support the longest run is 4 (exact). Covering a support needs at least 14
  runs.
- Twisted mirrors:
  - a run of any twist spans at most 5 dimensions, so it has at most 16 members [proved: the differences of a run are
    orthogonal to A applied to its span];
  - exhaustively, over all 69,616, the best is 8 (table below).
- Any second line at all: a run has at most P7 + 1 + 8 = 45 < 56 members, where P7 = 36 is the most ports in a
  7-dimensional subspace [proved, computed]. So **for no choice of second lines can one register clean a U8 total in
  a single run.** That is exactly what U8's end gate does.

Exhaustive run over all 69,616 twisted mirrors [computed, `job-twistall/twistall.py` on a Colab, exact maximum
runs by branch and bound; the natural mirror's 4 checked independently by `calc/vclique.py`]:

| longest run in the worst support | number of twisted mirrors |
|---:|---:|
| 4 | 3,333 (the natural mirror among them) |
| 6 | 56,539 |
| 8 | 9,744 |
| more than 8, in any support | 0 |

So **no twisted mirror lets a register visit more than 8 x's of a U8 total in a row**, and every total needs at least 7
runs, under every twist.

**First-order cost** [computed, `chaincost.py`; the model is argued and deliberately optimistic].
- After its U-drop, a total's content σ_k sits on one register at frame 0. The register removes σ_k run by run.
- Before each new run it must drop out of whatever its frame holds that the next run's lines forbid. That costs the
  codimension.
- The x's are assumed to stay on their own lines, which is optimistic. A greedy choice of runs (exact maximum runs,
  cheapest drop first) gives:

| lines | runs per total | cleaning drops per total | gathering + cleaning per call | room |
|---|---:|---:|---:|---:|
| natural mirror | 14 | 27 to 31 | 462 | 34 |
| best twist of the exhaustive run | 10 | 17 | 272 | 34 |

Gathering is counted the same as cleaning, by the time-reversal symmetry of the frame rules [argued]. Sharing could
lower these figures, if one register served the four totals that contain each x: at most by ×4, which still leaves
115 and 68 per call [argued].

How robust is this, if only the number of runs is counted, at one level per drop (the absolute minimum)?
- Natural mirror: at least 13 drops per total for cleaning and 13 for gathering, so at least 208 per call. That is
  52 per call with maximal sharing, still above the room of 34.
- Every twist (runs of at most 8, exhaustively): at least 6 + 6 drops per total, so 96 per call. With maximal sharing that falls to
  24, inside the room.
- So for the twists the verdict rests on the drop sizes. The greedy finds about 2.8 levels per drop (17 over 6
  drops), not on the run count alone. The model also leaves out three costs that would only add more: junk helpers obey the
same visiting lemma, the x's frames grow as they meet helpers, and helpers that cannot be hosted pay for their climbs.

**U8 at its own gate order** [computed, exact, `screen.py`]. This is pr717's legality propagation in subspace form:
the join of the start frames that reach each gate must lie inside the meet of the ceilings reachable from it. The
dirt gates g1 and g8 are left out (two-sided frames forbid them outright). Slots start at 0. The scatter is modelled
as fresh copies after the paid drop.

| case | illegal gates of 7,948 | on them |
|---|---:|---|
| one-sided calibration | 0 | |
| natural mirror, one copy per total | 2,453 (all 448 scatter reads, 2,005 in phase B) | all 112 y's, 416 scratch |
| natural mirror, one copy per (total, reader) | 2,429 | all 112 y's, 856 scratch |
| a twisted mirror (best of `twist.py`'s 40), one copy per (total, reader) | 2,360 | all 112 y's, 856 scratch |

## 7. What is left open

- **Second lines outside the twisted family.** These are g with the junk condition whose second Gram differs from
  the ports'. The second centre would then not be a P112 copy, and its price is unknown.
  - The general run bound (45 < 56) still holds.
  - My local search found none in 90 seconds; it is inconclusive.
  - An exhaustive search would be a constraint problem: odd h_b with u_b · h_b = 1 and h_a · h_b = 0 on every
    orthogonal pair. Colab-sized.
- **Centres not organised as one register per total.** The visiting lemma and the run bounds hold for every word. The
  cost table assumes U8's code (8 totals of 56) and cleaning against the x's. A code with smaller supports pays more
  centre charge; a code built only on non-orthogonal sets has no side information, so it saves nothing. That
  trade-off is [argued], not proved.
- **The Lean side** of the two-sided gluing was not written.

## 8. Predictions (STATUS.txt, written 12:40) against outcomes

| prediction | outcome |
|---|---|
| M1 (85 %): natural mirror passes; v-Gram = u-Gram | **held**, exactly; and every twisted mirror is a symmetric involution |
| M2 (80 %): R2 wrong as worded; shared frame forces shared s | **held**: meet 2-dim per total, forced to s by the junk condition |
| M3 (70 %): formulas give 1 per line (gauged), and 112 on s (fixed), which kills | **half held**: the counts are as predicted, but the fixed reading's floor, 564 per replica, stays below the gain of 896, so it does not kill on its own |
| M4 (60 %): forward calls need v-drops; floor small, gain left | **held** |
| M5 (60 %): "counts leave a gain; no word; obstruction named" | **held**, sharper than expected: runs of 4 inside a support, cleaning ≥ 272 per call against a room of 34 |
| M6 (97 %): nothing near 1e-1, nothing above 3.4e-3 | **held**: ceiling 3.28e-3, unreachable |

## Files

All in `checks/wht43/anatomy/cx-mirror/`:
- `STATUS.txt` (dated predictions, progress, outcomes), `RUNLOG.txt`, this `NOTE.md`.
- `calc/lines.py` → `lines.out/json`: ports, junk, the natural mirror, R2 meets.
- `calc/counts.py` → `counts.out/json`: the formulas, call by call, in both readings.
- `calc/design.py` → `design.out/json`: junk inside supports, run covers, P7, the V-drop.
- `calc/price.py` → `price.out/json`: ceiling and break-even.
- `calc/twist.py`, `twist_lib.py`, `twist2.py`, `twist3.py`, `perms.py`, `qcheck.py` (with `.out/.json`): the twisted
  mirrors.
- `calc/screen.py` → `screen.out/json`: U8 at its own gate order; `calc/psi.py` → `psi.out`: the mirror map.
- `calc/chaincost.py` → `chaincost.out/json`: the first-order cost of cleaning.
- `calc/gsearch.py` → `gsearch.out/json`: the inconclusive search outside the twisted family; `calc/vclique.py` → `vclique.out`: an independent check of the run bound 4.
- `job-twistall/twistall.py`, `data.json`, results in `job-twistall/twistall.json`: all 69,616
  twisted mirrors; `calc/chaincost2.py` → `chaincost2.out/json`: the cleaning cost of the exhaustive best.

## Sources

- Briefs: `../BRIEF3.md` (cx-mirror), `../BRIEF.md`, `../BRIEF2.md`; `../NOTE.md` rounds 1 and 2.
- `../an2-word/NOTE.md` section 4 (R1 to R5, the formulas), `../v-an2-word/NOTE.md` (the whole-run bound, prices),
  `../an2-smallnet/NOTE.md` sections 2 to 4 (side information, the second-line orthogonality), `../v-an2-smallnet/`.
- `../an-build/` (calib.py, twosided_local.py, dirtjoin.py, bwdpaths.py; the dirt gates), `../an-synth/NOTE.md` 4.1
  and 4.1a, `../v-an-build/calc/vports.out`.
- Read as text or data: `tools/gx/gxcore.py` (gcert/1 frame rules);
  `community/u8xmhh/overlay/Work/XHost/HDef.lean` (gluing x01 to y78);
  `bound-check/k-lower/out/U8X-parsed.json.gz`; `bound-check/fsp-price/data/words.json` (the U8 profile).
