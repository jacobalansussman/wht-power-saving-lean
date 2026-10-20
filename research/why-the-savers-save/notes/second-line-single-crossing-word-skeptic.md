> Working note written by an AI agent during the round of 15 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# v-an2-word: skeptic's check of an2-word ("the single-crossing local word does not exist")

Anatomy round 2, skeptic of the an2-word track, 15 October 2026, 11:19 to 11:36 EDT (Mac clock). All
computations are my own code (stdlib Python, `python3 -I`, under a private memory-capping wrapper, each under 70 s and well under 1 GB). No
Colab was needed. Nothing was posted. No third-party code was run.

Labels: **[computed]** = a program in `calc/` and its output; **[read]** = taken from a named file; **[argued]** = my
reasoning, which no program has checked.

## Verdict

**Not refuted in its outcome.** The lead is excluded: a two-sided word with one shared extra line s can never beat
its own one-sided projection, so the single-crossing local word cannot reach the record. an-build's ceiling (D ≤ 768,
about 3.2e-3) is rightly withdrawn. The Coxeter–Todd screen reproduces with an independent construction and fails.

**One quantitative statement is wrong as worded, and the verdict does not depend on it.**

1. **"The best any two-sided word with a shared s can do is D ≤ 208 per replica" is false without an extra
   hypothesis.** The per-call count (at least 112 paid crossings in each of the 9 calls, 1,008 per replica) holds only
   if no X array ever holds s, that is, if every data array keeps its s-status at every call boundary. The claim's
   list of hypotheses does not name this one.
   - Without it, only a whole-run information bound holds: per pair of banks, 112 x channels must enter the s side
     and 112 y channels must leave it, so 896 per replica (plus one climb per bank, which the banks' own extra
     direction cancels). That gives **D(two-sided) ≤ D(projection)**: a tie at best, not 112 worse.
   - The tie is reached by a trivial word: run the one-sided word with s added to every frame. Every X array climbs s
     at the start and drops it at the end (2 moves each); Y arrays start on ⟨s⟩ free. That is D = 320 per replica at
     m = 74 (the "trivial lift").
   - The difference comes from the three calls on the first pair (stages 0, 1, 2): 336 crossings in the per-call
     count against 224 in the whole-run count. The other pairs have two calls each, and there the two counts agree.
2. **The tie still loses.** Priced with the same model (U8 profile, exact fill, the Lean fold root), the trivial lift
   gives **1.6246e-3** with its crossings as rank-1 pieces, and **1.6538e-3** if they were absorbed into the idle fill
   (the most optimistic case). Both are below the one-sided exact fill, 1.6647327e-3, and below the hub, 1.6797550e-3
   **[computed]**. To first order a two-sided word is dominated child by child by its projection [argued; below].
3. **Smaller points.**
   - "Priced at most 1.0545261e-3" is the figure with the crossings as rank-1 pieces. Absorbing them into the fill
     gives 1.0758e-3. So 1.0545e-3 is not an upper bound in the model, though nothing turns on it **[computed]**.
   - "D ≤ 208" as an absolute figure also assumes the projected one-sided word has D ≤ 320, which rests on the centre
     floor of 64 per call. That floor has been checked for one class of decoders only. The robust statement is the
     relative one: no better than the word's own projection.

## 1. The crossing count, rebuilt

**The status argument** [argued, then computed on toys]. Shared s gives each register one status bit: its frame
holds s or not. A free gate joins equal frames, so it never mixes statuses. Gates within one status act on the
register tuple, so they cannot enlarge the span of that status's contents; throwing a scratch copy away can only
shrink it. Only a paid frame change of a register holding nonzero content moves a vector across. So:

- climbs across s ≥ dim(span(S_start ∪ needed S contents)) − dim(S_start);
- drops ≥ the same with the two sides swapped.

This is register by register. It does not care what a register is used for (data, bank, copy, hosted helper), so the
hosting question in the claim does not open a hole: a host carrying its helper's dirt is one register with one
vector, and the final s side needs the host's own content and the channel separately [argued]. `vbank.py` confirms it
on a toy: a bank's single climb never carries a channel.

**Two counts, depending on what is fixed at the call boundaries** [computed, `calc/vmulticall.py`]. This is an exact
0-1 search on the toy (1 x, 1 y, copies, GF(3)). Each call must leave the lifting-chain contents on the data registers,
so it is the identity on arbitrary inputs.

| run of calls | statuses fixed at every call boundary (claim's setting), copies discarded | statuses free between calls |
|---|---:|---:|
| GF | 1 | 1 |
| GB | 1 | 1 |
| GF, GB (or GB, GF) | 2 | 2 |
| GF, GB, GF (pair 1 in the template) | **3** | **2** |
| GF, GB, GF, GB | 4 | 2 |

- With statuses fixed, the minimum is one crossing per call (the claim's lemma).
- With statuses free, it is 2 whatever the number of calls: one channel in and one out.
- If copies must be returned to zero instead of discarded, the fixed-status counts double (2 per call). That only
  strengthens the lemma in its own setting.

With a bank register as well (`calc/vbank.py`, GF(2)): fixed = calls + 1, free = 3, as the whole-run bound predicts.

**The template's schedule** [read: `bound-check/h-L3-design/gen/template/Gates.lean`, `chainQ`]:
- GF at stages 0, 2, 4, 6, 8 and GB at 1, 3, 5, 7.
- The twin gate layers join X with X and Y with Y only, so they are status-clean.
- The run's result is X_p = −y_p, Y_p = x_p, with the banks unchanged.

So per replica the whole-run bound is 4 × 224 = 896 crossings, plus one per bank. The per-call bound with fixed
boundaries is 9 × 112 = 1,008, plus one per bank. In the fixed-boundary setting, each call's 112 x contents are
independent modulo the s side, because the X and Y contents together span all inputs. So the claim's 112 per call is
right there [argued].

**The projection identity** [computed, `calc/vproj.py`]. I tested 120,000 random elementary climbs and drops in F₂⁹,
including frames that hold s + a vectors and not s. Cost minus projected cost equals [the step crosses s] with 0
failures. The two-sided endpoint frames project onto the one-sided ones. So D(two) = D(projection) + W − (crossings)
is exact, and with crossings ≥ 896 + banks it gives D(two) ≤ D(projection).

**Why a tie cannot price higher** [argued, first order; computed for the U8 profile]. A child of rank r that does not
cross s costs r·ln(74/r) > r·ln(73/r). One that crosses costs (r+1)·ln(74/(r+1)) ≥ r·ln(73/r): the difference is
decreasing on [0, 73] and equals 0 at r = 73. So the narrow mass can only grow, and D can only fall.

## 2. Prices [computed, `calc/vprice.py`]

Model: the Lean fold hypothesis, my own Decimal bisection, U8 local profile × 72, the rest of the fill at rank 48.
Calibration: 1.6647327e-3.

| case | m | D per replica | price |
|---|---:|---:|---:|
| one-sided U8, exact fill (calibration) | 73 | 320 | 1.6647327e-3 |
| claim's floor, 1,008 rank-1 crossings | 74 | 208 | 1.0545261e-3 (reproduced) |
| claim's floor, crossings absorbed into the fill | 74 | 208 | 1.0758216e-3 |
| trivial lift, 2 rank-1 crossings per X array | 74 | 320 | 1.6246236e-3 |
| trivial lift, crossings absorbed (most optimistic) | 74 | 320 | 1.6537530e-3 |
| the hub's record (T74) | | | 1.6797550e-3 |

Every two-sided case is below the one-sided exact fill and below the hub.

## 3. Coxeter–Todd, rebuilt independently [computed, `calc/vct12.py`, `calc/vct12_phase.py`]

My construction is not the claim's. The claim used x mod 2 in the hexacode. I used E⁶ with all coordinates congruent
mod θ = √−3 and coordinate sum ≡ 0 mod 3. Its minimum norm is 6 (norms 1 to 5 are checked empty by brute force).

| quantity | claim | mine |
|---|---:|---:|
| minimal vectors | 756 | 756 (270 of shape (θu, −θu′, 0⁴), 486 of shape (u₁ … u₆)) |
| real lines; products 0 / ¼ / ½ from each line | 378; 135 / 160 / 82 | 378; 135 / 160 / 82 |
| real junk per port; tight | 242; yes | 242; yes |
| real junk clique (window ≥) | 45 | 45 (2,000 greedy restarts) |
| complex lines; products 0 / ¼ from each line | 126; 45 / 80 | 126; 45 / 80 |
| complex junk per port; clique; tight | 80; 15; yes | 80; 15; yes |
| copy ratios | 1.43 / 0.71 | 1.43 / 0.71 |

The phases claim also holds, and is a little stronger than stated. Over all mutually non-orthogonal triples of complex
lines, the triple products have phases 0°, ±60° and 180°, and 40,320 triples are at ±60°. Triple products do not change
when line representatives are rephased, so no choice of phases makes every product a power of i. Both forms fail the
screen (real form on all three marks; complex form on window and copy ratio).

## 4. What I did not check

- I did not rebuild a local word or check a forward identity on raw arrays. an2-word built no word (nothing passed),
  and the trivial lift is the one-sided word with one extra coordinate, so there was nothing new to verify.
- I did not check the claim's mirrored-variant estimate (D ≈ 640, about ×1.9) or its per-line count formulas. Both
  are labelled [argued] and unbuilt. The whole-run argument above does not cover family-specific lines, because one
  paid step can cross several lines at once. Section 4 of the claim says this too.
- I did not check that s must be a single shared coordinate (an-build's derivation). I took it as given.
- I did not check an-build's LP figures beyond confirming that an2-word's rerun output is byte-identical to an-build's
  (`dirtjoin`, `mintrans2`, `rank0`).

## Predictions (STATUS.txt, written 11:23) against outcomes

| prediction | outcome |
|---|---|
| P1: core exclusion right; projection identity holds for all subspace steps (85 %) | **held** |
| P2: the per-call 1,008 needs "X never holds s"; with free boundary statuses the bound is 896 per replica, a tie at best; pair 1's three calls make the difference (65 %) | **held**: toy 3 against 2; the trivial lift reaches D = 320 at m = 74 |
| P3: verdict still "excluded"; the tie prices below one-sided and the hub (85 %) | **held**: 1.6246e-3 to 1.6538e-3 |
| P4: schedule GF 0,2,4,6,8 / GB 1,3,5,7 (85 %) | **held** |
| P5: Coxeter–Todd figures reproduce independently; fails (75 %) | **held**, every figure |
| P6: 1.0545e-3 reproduces, but "at most" is loosely stated (50 %) | **held**: absorbed crossings give 1.0758e-3 |
| P7: not refuted in outcome; D ≤ 208 overstated (65 %) | **held** |

## Files

All in `checks/wht43/anatomy/v-an2-word/`:
- `STATUS.txt` (dated predictions, outcomes), `RUNLOG.txt`;
- `calc/vmulticall.py` → `vmulticall.out`, `vmulticall.json` (per-call against whole-run crossing counts);
- `calc/vbank.py` → `vbank.out`, `vbank.json` (the bank case);
- `calc/vproj.py` → `vproj.out` (the projection identity and endpoints);
- `calc/vprice.py` → `vprice.out`, `vprice.json` (prices);
- `calc/vct12.py` → `vct12.out`, `vct12.json`; `calc/vct12_phase.py` → `vct12_phase.out` (Coxeter–Todd).
