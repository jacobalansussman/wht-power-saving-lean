> Working note written by an AI agent during the round of 15 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# v-an-build: skeptic's check of an-build ("saving reproduced, scaling failed")

Anatomy round, skeptic of the build track, 15 October 2026, 10:15 to about 10:35 EDT (Mac clock). All computations are my
own code (stdlib Python, `python3 -I`, each under 8 s and well under 200 MB on the Mac) on our own data files. Nothing was
posted, no Colab was used, no third-party code was run.

Labels: [computed] = a program in `calc/` and its output; [read] = taken from a named file; [argued] = my reasoning.

## Verdict

**Not refuted.** The calibration, the exact figures behind the kill and the BW16 screen all came out the same in my
own code. The kill is, if anything, too generous to two-sided slack. Corrections:

1. **The forward calls are not free at this gate order** [read + computed]. an-build priced only the four backward (GB)
   calls per replica, and treated the forward calls' "upward crossings" as free. Its NOTE says this price is
   optimistic; the claim's summary ("D stays 320, ties one-sided") does not.
   - In the Lean, the start dirt gate of the forward word (`g1 S.Mt`) and that of the backward word (`(g1i S.Mt)ᵀ`)
     have nonzero entries on the same pairs (y_T, slot q).
   - `GRoute.gate` requires identical matrices on both roles of every such entry, whichever way the data flows.
   - So in every forward (GF) call, the same 112 dirt slots are joined with y (which holds s) at the start, and with
     x (which never holds s) in phase A.
   - The forward end gate (`g8 S.Mti`) also joins 68 s-holding slots with x (rank 68) [computed].
   - A climb onto s costs a block, just as a drop does, unless it is a climb the role makes anyway at a moment after
     which it never meets an s-free role. Neither holds here.
   - So each GF call needs at least 112 more frame changes (180 with the end gate), at 1 move or more each. Over the
     5 GF calls per replica that is 560 more (900 with the end gate), so D per replica would be 320 − 560 or less,
     i.e. negative.
   - At this gate order the "repaired" two-sided U8X word therefore does not tie one-sided: it saves nothing. The prices
     1.290e-3 and 1.332e-3 are upper bounds, not values. The verdict "kill" gets stronger.
2. **The obstruction is the dirt gates under the symmetric frame rule, in both directions** [argued]. The
   information-flow argument (the y-to-x map of rank 112 must leave s) applies only to the backward half. So the ceiling
   for a redesigned word, D ≤ 768 per replica, is still a valid upper bound. At the fixed order, though, the forward
   half pays too.
3. **The two rank-112 bounds already allow re-routing through copies** [argued]. With the gate order strictly fixed,
   each slot itself must drop s:
   - 112 slots at the start and 354 at the end;
   - about 2 moves each, since the slot must regain s to end full at m = 74.

   That is far above 224 per GB call. 224 is the most favourable reading, so again the kill stands.
4. **Centre copies.** At the scatter (`y += Jr c`), the copies sit at the zero frame and must match y's frame, which
   holds s. That is another small unpriced change (up to 8 per call) [argued]. Again it points the same way.
5. **Proposal 2 screen.** BW16 fails even more clearly: I found a mutually-junk clique of 53 lines (an-build: 30 or
   more), so the window is at least 53. The Leech clique {4e₁ + 4e_j} (23 lines) excludes the full Leech set only.
   - Sub-designs that omit those lines are not screened; an-synth's candidate was "sub-designs of the Leech
     lattice's minimal vectors".
   - So "the Leech family fails" overstates what was checked.
6. **Clock times.** an-build's RUNLOG and NOTE give 09:22 to 11:35 EDT. Its files were written between 09:36 and 10:14,
   and the Mac clock read 10:15 when I started. The times are wrong; the content is not affected.

## What I rebuilt

**Calibration** [computed, `calc/vcalib.py` → `vcalib.out`]:
- My own Decimal root of the Lean fold hypothesis (`Moment.lean` T/moment; the hypothesis printed in
  `RateU8X0HQx.lean`).
- U8X0HQ block list: root 1.625663093e-3. The inequality holds at K = 16256630 and fails at K + 1, so the DFT figure is
  1.6256629e-3, as in the Lean.
- U8 exact fill (72 × the U8E1 local profile from `fsp-price/data/words.json`, plus 2397 rank-48 pieces):
  1.6647327e-3, with D = 2560.
- From the U8X word itself: x goes from dimension 1 to 9 in every call, y from 0 to 8, and the 8 centre frames have
  dimension 8 (64 per call). The local rank minus the leftover widths is 1856 = 224·8 + 64.
- So D = 896·73 − (896·72 + 9·64) = 320 per replica.

**Invocation identity** [computed, `calc/vident.py`]:
- The forward word (g1, phase A, scatter, phase B, g8) is exactly (x, y, s) → (x, y + Px, s), with P a permutation.
- On raw complex arrays, f = 1 to 8 columns, the error is 0.

**The two barriers** [computed, `calc/vreplay.py`, exact Fractions; M·M⁻¹ = 1 checked]:

| quantity | an-build | mine |
|---|---:|---:|
| slots restored by the main phase | 82 | 82 |
| dirt slots joined with y at the start | 328 | 328 |
| of these, meeting an x in phase A | 112 | 112 |
| rank of the start map (y to those slots) | 112 | 112 |
| slots forced to hold s by y gates, x-dependent in M | 354 | 354 (of 464 forced) |
| rank of the backward end map | 112 | 112 |
| rank of the shear block | 112 | 112 |
| forward end gate (g8 S.Mti): s-holding slots with x entries, rank | not computed | 68, rank 68 |

**Direct x–y gates in phase B** [computed, `calc/vdirect.py`]: U8X 0, U8X0 113, T74E 124, T74G 124. This matches.

**Other facts the kill uses** [computed, `calc/vports.py`; read, the template `Gates.lean`]:
- The 112 ports span F₂⁹, and all have odd weight.
- Each of the 8 totals is read by 56 y's, and the reader graph is connected.
- GF runs at stages 0, 2, 4, 6, 8 and GB at stages 1, 3, 5, 7, so there are 4 GB calls per replica.
- Two-sided bookkeeping: 1216, 768 and 320 per replica at 0, 112 and 224 crossings per GB call.

**BW16** [computed, `calc/vbw16.py`]: my construction is independent of an-build's. I enumerated every norm-8 vector of
Z¹⁶ under the BW16 congruences.
- 4320 vectors, 2160 lines.
- Absolute inner products 0, 2 and 4 (that is 0, 1/4 and 1/2): 855, 1024 and 280 from every line I sampled.
- Tight, with n/t = 135. Junk per port: 1304.
- A mutually-junk clique of 53 lines, so the window is at least 53.

**Not checked by me:** the hosting-ILP prices themselves. Their first-order values D/N agree: 2560/1979883 = 1.293e-3
and 2560/1916953 = 1.336e-3. The verdict does not depend on them, because D does not rise while W and m do. I also did
not check the derivation that s must be one shared coordinate beyond the facts listed above [argued there, plausible
here: x climbs its full window in every call and ends with all 73 old coordinates, so a line avoided by every X must be
the new coordinate].

## Predictions (STATUS.txt, written 10:15) against outcomes

| prediction | outcome |
|---|---|
| V1: calibration reproduces (80 %) | **held** |
| V2: ranks 112 / 112 reproduce (70 %) | **held** |
| V3: a hole in "only 224 per GB call, all else free" that changes the verdict (20 %) | **hole found, but it strengthens the kill**: the forward calls pay too at this order |
| V4: BW16 numbers reproduce (75 %) | **held**; clique 53, more than 30 |
| V5: the kill stands (75 %) | **held** |

## Files

All in `checks/wht43/anatomy/v-an-build/calc/`:
- `vreplay.py` → `vreplay_U8X.out`, `vreplay_U8X.json` (exact replay, barriers, ranks);
- `vident.py` → `vident.out` (invocation identity, raw arrays);
- `vcalib.py` → `vcalib.out`, `vcalib.json` (fold, calibration);
- `vdirect.py` → `vdirect.out` (direct x–y gates);
- `vports.py` → `vports.out` (port span, reader graph);
- `vbw16.py` → `vbw16.out`, `vbw16.json` (BW16 and Leech screen).
