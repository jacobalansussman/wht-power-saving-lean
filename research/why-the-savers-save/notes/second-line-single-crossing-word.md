> Working note written by an AI agent during the round of 15 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# an2-word: the single-crossing local word does not exist

Anatomy round 2, track an2-word, 15 October 2026, 10:39 to 11:20 EDT (Mac clock). Task (BRIEF2.md): design a
local word for the hub's template whose backward (GB) stage crosses the extra line s only once, test it with an-build's
cheap test, price any pass; and screen the Coxeter–Todd lattice as a port design.

All computations are my own code, stdlib Python except the rerun of an-build's own scripts, each on the Mac under
a private memory-capping wrapper, under 35 s and well under 1 GB. No Colab was needed. Nothing was posted. No third-party code was run.

Labels: **[computed]** = a program in `calc/` and its output; **[read]** = taken from a named file; **[argued]** = my
reasoning, which no program has checked.

## In short

1. **No local word can pass the brief's test, and the lead is closed.** The test asks for at most 112 crossings per
   backward call and none in the forward calls. The forward calls can never be free:
   - A free gate joins only arrays at one frame. So an array's content can move from "frame without s" to "frame
     with s" only when some array changes its frame across s. That change is a paid move of at least 1.
   - A forward call must deliver 112 independent channels of x (which never holds s) into y (which always holds s).
     So it needs at least 112 crossings.
   - None of them can ride on a climb that is paid anyway. No X array ever climbs s and every Y array holds it from
     the start. A bank does climb s once, but its own content needs that climb.
   - A backward call needs at least 112 crossings the other way (drops), which are never free.

   So every call pays at least 112, whatever the gate order, the dirt gates, the copies or the hosting. **[argued;
   checked exhaustively on small words, below]**
2. **The deficit can only fall.** Per replica the extra charge is at least 9 × 112 = 1,008. The second line gains 896.
   - So D ≤ 1216 − 1008 = **208 per replica**, against 320 one-sided, at the larger block m = 74.
   - Equivalently: removing s from every frame turns any two-sided word into a one-sided word of the template. That
     one-sided word has at least 112 more free moves per replica, at m = 73.
   - Priced as optimistically as possible (U8 profile, exact fill), the best two-sided word is worth at most
     **1.05e-3**. That is ×0.63 of the one-sided exact fill (1.66e-3) and below the hub's 1.68e-3. **[computed]**
3. **an-build's ceiling is withdrawn.** It allowed D ≤ 768 per replica (about 3.2e-3), because it counted only the
   backward calls and took the forward calls' crossings as free. The skeptic of an-build had already found that the
   forward calls pay at U8X's gate order. The lemma shows they pay at every gate order. **[argued, computed]**
4. **The Coxeter–Todd lattice fails the screen** in both forms **[computed]**:
   - real form: 378 lines, 242 junk terms per port (mark ≤ 120), a mutually-junk clique of 45 (so window ≥ 45, mark
     ≤ 12), copy ratio ≥ 1.43 (mark ≤ 0.35);
   - complex form: 126 lines in C⁶, 80 junk terms per port, clique 15 (window ≥ 15), copy ratio ≥ 0.71, worse than
     E8's 0.60. Its phases are cube roots of unity, which are not among the game's free phases (powers of i)
     **[argued]**.
5. **What is left of two-sided slack** is a second line that differs from family to family (v_b for family b), not
   one shared s. This is a mirrored second index code. The lemma does not cover it. In the U8 word the start dirt gate
   and the scatter force one shared s (an-build), so it needs a new local word built around it. Section 4 lists what
   such a word must satisfy and the count it must beat. Its first-order value is about ×1.9 (an-synth's ×2.03 row),
   but nothing is built **[argued]**.

## 1. The crossing lemma

**Setting** [read: an-build, v-an-build; the frame rule of our engine]:
- In the two-sided template the second line s is one coordinate shared by every y and held by no x (an-build,
  `twosided_local.py`: the ports span F₂⁹ and the scatter's reader graph is connected).
- A free gate needs one common frame, so all its arrays agree on "holds s" or "does not".
- A frame change that adds or removes s is a climb or a drop of at least one move.
- Helpers live on host arrays and their frame is the host's frame. A host is an idle X array (never holds s), an idle
  Y array (always holds s) or a bank (frame 0 at the start and full at the end, so it climbs s exactly once in its
  life). Copies are fresh arrays; every frame change they make while holding content is paid.

**Forward call (y += P x, P a permutation, rank 112).** Track the span of the contents of the arrays that hold s:
- A gate never enlarges that span, because it only recombines arrays on one side.
- Only a climb across s can bring in something new, and each climb brings in at most one new direction.
- At the end, the arrays holding s must contain all 112 values y_b + x_b and the restored contents of every bank that
  climbed s in the call. That needs 112 + (number of such banks) new directions.
- Each bank climbs once, so at least 112 climbs are not a bank's own climb. They are paid.

**Backward call (x −= Pᵀ y).** The same argument with drops. A drop is never part of any array's planned climb, so at
least 112 drops are paid.

**Per replica:** 5 forward and 4 backward calls give at least 1,008 paid crossings. Each costs at least 1 move. The
second line gains one free move per data array, 896 in all. So D falls from 1216 to at most 208.

**Projection form.** Delete s from every frame (map each frame f to (f + ⟨s⟩)/⟨s⟩):
- Equal frames stay equal.
- X: ⟨u⟩ → s^⊥ becomes ⟨u⟩ → everything, and Y: ⟨s⟩ → u^⊥ ⊕ ⟨s⟩ becomes 0 → u^⊥. That is exactly the one-sided
  contract at m = 73.
- Every climb or drop costs exactly 1 less in the projection when it crosses s, and the same otherwise
  **[computed, `quotient.py` (a): 40,000 random climbs and drops in F₂⁷, 0 failures]**.
- Hence D(two-sided, m = 74) = D(projection, m = 73) + W − (crossings) ≤ D(projection) − 112 per replica.

In this template the second line is therefore never worth having, however the local word is built.

**Checked exhaustively on small words** **[computed, `toybfs.py`]**. The model has two statuses (holds s or not),
free gates between equal statuses, and paid status changes, searched exactly (0-1 BFS over all contents in GF(3) or
GF(2)). In all 13 cases the minimum number of paid changes is at least rank P + (slots whose status differs from start
to end); in one of them (no copy and no slot) no word exists at all, since nothing can cross. Equality holds when
copies are available:

| case | minimum | bound |
|---|---:|---:|
| 1 x, 1 y, one copy: forward / backward | 1 / 1 | 1 / 1 |
| a slot that must climb s anyway (a bank), plus one copy | 2 | 2 (its own climb + 1 paid) |
| the same slot, no copy | 3 | 2 |
| 2 x, 2 y, P = I, one copy: forward / backward | 3 / 3 | 2 / 2 (one copy must go back) |
| 2 x, 2 y, P of rank 1, one copy | 1 | 1 (the count is the rank, not the number of targets) |

The slot case is the one that matters: a bank's free climb never pays for an x-to-y channel.

## 2. The cheap test, rerun, and why "no cost in the forward calls" cannot be met

an-build's scripts reproduce exactly on U8X **[computed, `calc/rerun/`; the output JSON equals an-build's]**:
- 328 dirt slots, 112 of them meeting an x in phase A;
- an integral LP with 112 crossing copies per backward call and 112 switching slots;
- ranks 112 at the start and 112 at the end, so at least 224 per backward call at this gate order.

The LP's 112 "switching" slots climb s inside the forward word, and an-build counted that climb as part of the slot's
hosted width. No host can absorb it:
- X hosts never hold s, and Y hosts hold s from the start.
- an-build's `hosty.py` therefore sent these slots to banks.
- But a bank's single climb of s is used by the bank's own content (section 1, and the toy case above).

So those 112 climbs are paid in every call. The fixed-order word costs more than the 224 per backward call that the
kill used, as the skeptic of an-build found.

There is no other gate order, dirt-gate merge or copy placement to search. The lemma uses none of them: 112 per call in
each direction is a floor for every word.

**Prices** (`quotient.py` (c); U8 local profile ×72, idle remainder at rank 48, crossings as rank-1 pieces; the Lean
fold root; calibrated on the one-sided exact fill) **[computed]**:

| case | m | D per replica | price |
|---|---:|---:|---:|
| one-sided U8, exact fill (calibration) | 73 | 320 | 1.6647327e-3 |
| two-sided, the lemma's floor | 74 | 208 | **1.0545261e-3** |
| two-sided, an-build's fixed-order tie (backward crossings only) | 74 | 320 | 1.6246236e-3 (not reachable) |
| two-sided, an-build's ceiling (forward calls free) | 74 | 768 | 3.92e-3 (impossible) |

The 1.05e-3 is an upper bound. The real hosting would lose more.

## 3. Proposal 2: the Coxeter–Todd lattice

**Construction** **[computed, `ct12.py`]**. Eisenstein integers E = Z[ω], with K12 = {x ∈ E⁶ : x mod 2 lies in the
hexacode} (E/2E = F₄). The checks:
- the hexacode has 64 words, of weights 1 × 0, 45 × 4 and 18 × 6;
- the minimum norm is 4: no lattice vector of norm 1, 2 or 3 exists (brute force);
- there are exactly 756 minimal vectors, enumerated rather than counted by formula.

The screen is the one used for BW16 and Leech (pass marks: copy ratio ≤ 0.35, window ≤ 12, junk per port ≤ 120).

| form | lines | centres t | n/t | tight | junk per port | junk clique (window ≥) | copy ratio ≥ | verdict |
|---|---:|---:|---:|---|---:|---:|---:|---|
| real, in R¹² | 378 | 12 | 31.5 | yes | 242 (160 at 1/4, 82 at 1/2) | 45 | 1.43 | fails all three |
| complex, in C⁶ | 126 | 6 | 21 | yes (1 + 80/4 = 21) | 80 (all \|A\|² = 1/4) | 15 | 0.71 | fails window and copy ratio |

- The cliques are lower bounds found by 400 randomised greedy runs. The plain label choice B = [A = 0] gives windows
  of 287 and 35. Either way both forms fail.
- The complex form is the closest to E8 of any design screened so far: 80 junk terms against E8's 52 to 56. But its
  copy ratio, 0.71, is worse than E8's 0.60.
- Leech sub-designs that leave out the 23-line clique remain unscreened.

## 4. What is left of two-sided slack: a second line per family

The lemma rests on one fact: s is a single status bit shared by every y. If each Y family b instead starts on its own
line v_b (X family b avoiding v_b), the count changes. Its exact form [argued]:
- In a forward call, the paid climbs of a line v number at least the number of families b such that y_b holds v at
  the start of the call and x_b never holds v during the call.
- In a backward call, the paid drops of v number at least the number of families b such that y_b holds v at the start
  of the call and x_b lacks v at its end.

For shared s both counts are 112. For family-specific lines both can be as small as the number of families sharing a
line. And one paid move that changes a frame by several lines crosses all of them at once. That is the index-code
effect that makes today's centre cheap (the one-sided backward call passes this count with its 64 per call).

So the mirrored proposal of an-synth (section 4.1, "mirrored v-lines") is not excluded. It needs a new local word,
because U8 forces one shared line:
- **R1.** The lines must differ by family. Shared s is closed by section 1.
- **R2.** No gate may join y's of different lines while each holds only its own line.
  - U8's start dirt gate does exactly that at time 0, and so does its scatter (each copy is read by 56 y's at one
    frame; the reader graph is connected) [read: an-build].
  - The scatter can be kept if the y's first climb a common frame F that contains all their lines and none of the
    readers' port lines. Such an F exists, because every port has odd weight [read: v-an-build `vports.out`]. F = ⟨s⟩
    plus any even-weight subspace avoids all 112 ports, so v_b = s + g_b with g_b of even weight would do [argued].
- **R3.** In U8's backward reading, the y information reaches the x's only through the start dirt gate (an-build
  `bwdpaths.py`: 0 of 112 without it).
  - A mirrored word must let that information drop F first and then reach the x's: a "gather" after the scatter.
  - My reachability check shows that in U8X no source's information ever reaches another source (`reach.py`: each x
    is reached only by itself) [computed]. So the sources themselves impose no constraint between lines. The
    constraints are R2 and R3.
- **R4.** If a second centre costs what the first does (64 per call) and the helpers stay hosted, then
  D ≈ 1792 − 2 × 576 = 640 per replica at m = 73 or 74. That is about ×1.9, near an-synth's ×2.03 row (3.4e-3)
  [argued, first order].
- **R5.** If the second lines share the port window (m = 73), the skeptic's orthogonality check applies: with the
  natural matching, at least 52 of 112 families fail it. If they take a new shared coordinate, m = 74.

## 5. Predictions (STATUS.txt, written 10:52) against outcomes

| prediction | outcome |
|---|---|
| W1: excluded for every local word; GF ≥ 112 per call; D ≤ 208; an-build's 768 wrong (85 %) | **held**: the lemma; toy search; projection identity; D ≤ 208, price ≤ 1.05e-3 |
| W2: at U8X's order an LP charging climbs gives ≥ 112 per GF call (80 %) | **held, without a new LP**: an-build's own solution climbs s in 112 slots per call, which no host can absorb (section 2) |
| W3: no word passes (90 %) | **held** (a proof, not a search) |
| W4: Coxeter–Todd fails; real-form junk 250 to 350 (90 %) | **held**: fails in both forms; real junk 242, just under my range |
| W5: what survives is a family-specific second line, worth at most ×1.2 (55 %) | **half held**: that is what survives, but its first-order value is about ×1.9, not ≤ ×1.2; unbuilt |
| W6: nothing near 1e-1 (99 %) | **held** |

## 6. Next steps

1. **Stop the single-crossing word.** It cannot beat one-sided slack in this template, whatever its gate order.
2. **If two-sided slack is pursued, the test to run next** is cheap and exact:
   - choose second lines v_b = s + g_b (g_b of even weight);
   - check R2 on U8's scatter with a common frame F;
   - check that a gather after the scatter exists (R3).

   Then count each line's crossings per call with the formulas of section 4 before any pricing. Only after that would
   a mirrored helper circuit be worth searching for.
3. **Proposal 2:** BW16, full Leech and now Coxeter–Todd all fail. Only Leech sub-designs without the 23-line clique
   are left to screen. I expect them to fail too.

## Files

All in `checks/wht43/anatomy/an2-word/`:
- `STATUS.txt` (dated predictions, outcomes) and `RUNLOG.txt`.
- `calc/toybfs.py` → `toybfs.out`, `toybfs.json`: the exhaustive check of the crossing lemma (13 small words).
- `calc/quotient.py` → `quotient.out`, `quotient.json`: (a) the projection identity, (b) the replica bound, (c) prices.
- `calc/ct12.py` → `ct12.out`, `ct12.json`: the Coxeter–Todd screen, both forms.
- `calc/reach.py` → `reach.out`, `reach.json`: reachability between the sources of U8X.
- `calc/rerun/`: an-build's `dirtjoin.py`, `mintrans2.py` and `rank0.py` rerun on U8X (calibration; identical
  output).

Sources read:
- `../BRIEF.md`, `../BRIEF2.md`, `../NOTE.md`;
- `../an-synth/NOTE.md` (sections 4.1, 4.1a, 4.2);
- `../v-an-synth/NOTE.md`;
- `../an-build/NOTE.md`, `STATUS.txt`, `RUNLOG.txt` and `calc/` (`mintrans2.py`, `rank0.py`, `dirtjoin.py`,
  `endstep.py`, `bwdpaths.py`, `hosty.py`, `bw16.py`, `twosided_local.py` and their outputs);
- `../v-an-build/NOTE.md` and `calc/` (`vreplay.py`, `vident.py`, `vcalib.py`, `vbw16.py`, `vports.out`);
- `../an-ours/NOTE.md`;
- the U8X lowering `bound-check/k-lower/out/U8X-parsed.json.gz` and `bound-check/fsp-price/data/words.json`, read as
  data.
