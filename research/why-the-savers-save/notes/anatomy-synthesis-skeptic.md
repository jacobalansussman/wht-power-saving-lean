> Working note written by an AI agent during the round of 15 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# v-an-synth: a skeptic's check of the anatomy synthesis

Anatomy round, skeptic of an-synth, 15 October 2026, 10:15 to 10:35 EDT (Mac clock). I read the synthesis, the three
dissections and the sources they rest on first-hand, and I re-derived the key figures with three small programs of my
own (`calc/`, stdlib, `python3 -I`, each under 2 s). Nothing was posted. No Colab was used. No third-party code was run.

Labels: [computed] = a program here and its output; [read] = from a named file; [argued] = my reasoning.

**Verdict: refuted.** The design numbers and the account of how the savers work are right. But the general cap does
not follow from its own stated assumptions, and it rests on walls nobody has proved. The top proposal's decisive test
has its kill mark set in the wrong place, and its central prediction is an assumption that the exact attempt did not
bear out.

## In short

What I confirmed:

1. **The design figures are exact** [computed, `calc/roots.py`, written from scratch].
   - E8: 120 lines, 8 centres, 56 junk terms per line, tight, window rank_F2(J + Gram) = 9, junk clique 8, copy ratio 0.600.
   - E7: 63 lines, 7 centres, 32 junk terms, tight, window 7, junk clique 7, copy ratio 0.778.
   - P112 (E8 minus an orthogonal frame): 112 ports, 8 centres, 52 junk terms, tight, 112/8 = 14 = 1 + 52/4.
2. **The mechanism is described as the dissections describe it** [read].
   - The identity paid = W·m − 2·Σ(source dims) + 2·(drops) is an-openai's.
   - One free coordinate per array, 72 of 73, is in an-ours.
   - The deficit has been 4.8924e-3 at every hub checkpoint since #588, per an-hub's lineage table.
3. **The proposal tables are arithmetic-correct within their own model** [computed, `calc/capcheck.py`]: ×3.88, ×2.03
   and ×0.20; 6.5e-3; ×1.05 to ×1.35; 0.096; and table (d).
4. **The proposals are consistent with the proved lower bounds, but those bounds test nothing here.**
   - The bounds: Theorem 2W (M ≥ max(W + m, 2W)), the Hecke bound (M ≥ max(2W, W + m − 1)), MAP §1.3 (total rank
     ≥ W + m − 1), and 14W/5 and 41W/13 at blocks 3 and 4.
   - Every proposal sits at block 8 or more with many arrays, where none of them binds. For example, proposal 3 pays
     14 of 16 per array: M = 14W, far above max(W + 16, 2W).

What I could not confirm, or found wrong:

5. **The cap figure "a ≲ 3e-3" does not follow from the derivation the note gives** [computed].
   - The note takes D/(W·m) from its design × template table and bounds l̄ ≥ ln(m/9).
   - With those two inputs, the best case is P112 at P = 1, m = 25: D/(W·m) = 5.71e-3, which gives a ≤ **5.59e-3**.
     At P = 2, m = 41, it gives 4.60e-3.
   - "≲ 3e-3" appears only if l̄ = ln(m/3.94), that is, if T74's own child grouping (r_eff = 3.94) is assumed.
   - That is not a bound. an-hub's ceiling row already places the idle mass at rank 51.
6. **The "general cap" is not a cap** [argued].
   - Wall 3 (stage and window, m ≥ 25) is argued only. It is exactly what proposal 3 would test.
   - "Over all designs" means a finite list (A_n up to n = 16, D_n up to n = 17, E6 to E8, Alon triples up to h = 30),
     all placed in one template (the hub's 2P + 1 stages).
   - The 1e-1 conclusion depends on Wall 3. The note's own table (d), column r_eff = m/√2, shows that one free
     coordinate per array would meet the bar up to m = 29 (0.85 needed at m = 25).
   - So the arithmetic alone does not exclude the class. Only the unproved wall, together with the listed designs'
     D/W, does: P112 has D/W = 0.143 at m = 25.
   - The synthesis does say these are "relaxation ceilings and arguments", and the summary itself is consistent with
     that. But "No construction in this class reaches 1e-1 at any block its walls allow" is stated as a result.
7. **Proposal 1's decisive test has a miscalibrated kill mark** [computed].
   - In the note's own model (`cap.py` (c)), two-sided slack beats the record only while the extra charge stays below
     **99.9 per call**, that is, about w.
   - An extra of 110 gives ×0.71, and 128 gives ×0.20, yet the note kills only at 128 or more.
   - an-build's measured average, 99.6 per call (224 per backward call), fell exactly in that gap, with zero net gain.
8. **The ×2.03 central prediction is an assumption, and it did not survive the exact attempt** [read; computed in
   part].
   - It assumes one extra rank-8 return per centre per call.
   - The note's in-call check (`calc/twosided.py`) models only the forward local word of the copy lowering U8X. It does
     not model the backward stages (x −= y, on the same frames) or the end dirt gate.
   - an-build found that exactly these force 224 crossings per backward call, for a net gain of 0 (×0.79 to ×0.82 at
     m = 74). an-build is not yet checked by a skeptic, and its RUNLOG times run ahead of the Mac clock.
   - Either way, "inside one call it costs nothing" is not established for the network as it runs.
9. **The note is inconsistent about where the second line lives** [read; computed in part].
   - Section 4.1 and `twosided.py` use mirrored v-lines in coordinates disjoint from the u-window, plus one extra
     coordinate (m = 74). The summary and 4.1a say no new coordinate is needed (m = 73).
   - The 4.1a condition u_a · v_a = 0 was never checked.
   - In U8X, the 112 port labels span F₂⁹, and every bit, including the coordinate all windows share, is set in 52 to
     56 of them.
   - So if X and Y families keep the same port index in two windows that share that coordinate, at least 52 of 112
     families fail the condition [counts computed; the matching argued].
10. **"One first-order law prices all three" overstates what was shown.**
    - The inputs D, W and r_eff (or l̄) are each unit's own. For OpenAI, whose children are all rank 1, the law is
      simply the definition of the saving.
    - The close agreement is a consistency check of the first-order expansion, not independent evidence for a common
      mechanism.

Smaller points:
- Proposal 3 is called "the only proposal that touches 1e-1", but its own figure is 0.096.
- "With 2 or 3 arrays any bridge costs at least 2m": an-openai argues this for W = 2 only.
- The "77 % of the narrow mass" figure is from our U8 word (an-ours), not from the hub's T74.

## 1. Design figures (confirmed)

`calc/roots.py` builds the 240 roots of E8 and takes lines. E7 is the set of lines orthogonal to one root; P112 is E8
minus the frame e1±e2, e3±e4, e5±e6, e7±e8. For each, it computes:
- the exact rational rank t of the Gram matrix;
- the junk per line;
- tightness, A² = (n/t)·A;
- the window rank_F2(J + Gram);
- the maximum cliques of the junk and orthogonality graphs (Bron–Kerbosch).

| design | n | t | junk | n/t = 1 + s/4 | window | junk clique | copy ratio |
|---|---:|---:|---:|---:|---:|---:|---:|
| E8 | 120 | 8 | 56 | 15 | 9 | 8 | 0.600 |
| E7 | 63 | 7 | 32 | 9 | 7 | 7 | 0.778 |
| P112 | 112 | 8 | 52 | 14 | 9 | 8 | 0.643 |

All of an-synth's figures match.

- The window 9 is a property of the label choice B = J + Gram. The floor set by the junk clique is 8, so the true
  minimum copy ratio of E8 lies between 0.533 and 0.600. The note states this correctly.
- "E8 is the best two-angle centre" rests on two things: the classification of line systems at 60° and 90° (cited
  from memory; the theorem is real) and the hub's study of deleting E8 ports.
- Subsets of large D_n, and label choices other than J + Gram, are not covered. So that claim stays argued.

## 2. The cap (not confirmed; the stated arithmetic is wrong)

`calc/capcheck.py` uses the note's own template model: m = 1 + (k − 1)(2P + 1), and D/W = 1 − (2P + 1)·cst/(2P·v).

| P112 | m | D/W | D/(W·m) | a with l̄ = ln(m/9) | a with l̄ = ln(m/3.94) |
|---|---:|---:|---:|---:|---:|
| P = 1 | 25 | 0.143 | 5.71e-3 | **5.59e-3** | 3.10e-3 |
| P = 2 | 41 | 0.286 | 6.97e-3 | 4.60e-3 | 2.98e-3 |
| P = 4 (T74) | 73 | 0.357 | 4.89e-3 | 2.34e-3 | 1.68e-3 |

- The note states l̄ ≥ ln(m/9) and concludes a ≲ 3e-3. That inference is wrong by almost a factor of 2.
- At P = 1 there is no idle cover to host the helpers (`an-ours/stages.out`), so the real figure there is lower. But
  the note's ceiling explicitly sets helpers free.
- With the note's third pricing column (r_eff = m/√2), u_eff = 1 meets 1e-1 for m up to 29. The conclusion that the
  bar is unreachable therefore comes from Wall 3 (m ≥ 25) together with the listed designs' D/W, not from the pricing.

## 3. Proposal 1 (not confirmed)

| extra charge per call | a (note's model) | × record |
|---:|---:|---:|
| 0 | 6.51e-3 | 3.88 |
| 64 (pass mark) | 3.41e-3 | 2.03 |
| 99.6 (an-build's average) | 1.69e-3 | 1.01 |
| 110 | 1.20e-3 | 0.71 |
| 128 (kill mark) | 3.38e-4 | 0.20 |

The break-even extra is 99.9 per call. That is where the kill mark should have been.

The 113 direct x–y gates (U8X0) and 124 (T74E, T74G) that an-build reports in phase B do exist
[computed, `calc/xyprobe.py`]. But every one of them joins orthogonal ports.
- Under an-synth's mirrored-line model they are legal.
- an-build's "illegal outright" holds only under its own single-shared-coordinate model.
- So that point does not refute the synthesis by itself.
- The backward-stage crossings do bear on it, if they are confirmed. The synthesis never modelled them.

## Files

All in `checks/wht43/anatomy/v-an-synth/`:
- `STATUS.txt`: predictions written before reading, and their outcomes.
- `RUNLOG.txt`.
- `calc/roots.py` → `calc/roots.json` (the design figures).
- `calc/capcheck.py` → `calc/capcheck.out` and `calc/capcheck.json` (the cap arithmetic, the proposal 1 break-even, and
  proposal 3).
- `calc/xyprobe.py` → `calc/xyprobe.out` and `calc/xyprobe.json` (direct x–y gates in our lowerings).

## Sources

- `../an-synth/NOTE.md`, `STATUS.txt`, `calc/cap.py`, `cap.out`, `ade.out`, `clique.out`, `twosided.py` and
  `twosided.out`.
- `../an-ours/identity.py`, `identity.out`, `knobs.out`, `stages.out` and `u8word.out`.
- `../an-openai/NOTE.md` sections 2 and 3.
- `../an-hub/out/anatomy.txt`, `caps.txt` and `lineage.txt`.
- `../an-build/NOTE.md`, `STATUS.txt`, `RUNLOG.txt` and `calc/viol_probe.py`, read as text.
- `../../meta/MAP.md` sections 1.2 and 1.3.
- `../../../wht33/HANDOFF.md`, the lines on Theorem 2W and the Hecke bound.
- Our lowerings in `../../bound-check/`: `k-lower/out/U8X-parsed.json.gz`, `fsp-cc/out/U8X0.gcert.json.gz`,
  `lg-a1/out/T74E.gcert.json.gz` and `T74G.gcert.json.gz`, read as data.
