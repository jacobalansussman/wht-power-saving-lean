> Working note written by an AI agent during the round of 15 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# v-cx-mirror: skeptic's check of "the mirrored second line is excluded"

Anatomy round 3, 15 October 2026. Two sessions: 13:30 to 13:51, then (after the first session died) 15:21 to 15:38,
Mac clock. All programs are my own (stdlib Python, run with `python3 -I` under a private memory-capping wrapper). They read the U8X
lowering and cx-mirror's exhaustive twist file as data. Nothing was posted and no Colab was used.

## Verdict

**Refuted as stated, but the lead stays closed in practice.**

- Every fact the exclusion is built on reproduces with my own code: the ports, the junk pairs, the mirror, the twist
  count, the run bounds, the price room, the visiting lemma, the R2 refutation and the forward-call v-drops.
- The step that does the excluding does not hold up. The claim prices the cleaning at "462 / 272 levels per call,
  under a deliberately optimistic model", and its STATUS says "cleaning + gathering >= 272 per call". Those numbers
  are the cost of one greedy scheme, not a lower bound.
  - Within the claim's own pricing, an explicit scheme with more registers costs **120 levels per call for every
    second line**, natural or twisted. That is 112 cleaning drops plus 8 extra centre levels. I checked it gate by
    gate, on frames and exact values.
  - The claim's own single-register model also has cheaper orders than its greedy found: 364 and 248 per call.
- For the twisted mirrors, the only floor the claim proves is its run-count floor with sharing, 24 per call. That is
  inside the room of 34. So the twist exclusion is not established even in the claim's own model.
- No scheme I found fits the room. The cheapest is 120 per call against a room of 34, which leaves D below today's.
  So the mirrored second line looks dead at first order, but by a margin of about ×3.5, not ×8 to ×13. And no proof
  covers every scheme.

## 1. What reproduces

From `calc/v1.py`, `v2price.py` and `v4.py` (first session) [computed]:

- U8X's phase A gives 8 totals of 56 x's, with coefficients ±1/2 and ±1/4. The delivered diagonal is 1.
- There are 5,824 junk entries, 52 per target, and they are exactly the binary-orthogonal port pairs.
- The natural mirror v_b = s + j + u_b behaves as the claim says: g_b is even, u_b·v_b = 0, and the v-lines have the
  same Gram matrix as the ports. The map ψ is an involutive isometry that swaps u_b and v_b.
- The longest drop-free run inside each total's support is 4. The longest in the whole port set is 8. P7 = 36.
- Quadratic forms that vanish on the 112 vectors e_b vanish on all of E (rank 36 = 36). That is the step the
  symmetric-involution proof needs.
- There are 5,355 isotropic planes and 2,295 Lagrangians, giving 69,616 twisted mirrors.
- The price model reproduces: the calibration gives 1.6647327e-3 and the ceiling 3.2824847e-3. The break-even is
  C ≤ 162 per call, which leaves a room of 34.
- No set of the readers' frames can serve as R2's single common frame.
- I sampled 251 random twisted mirrors. All pass, every one gives a symmetric N, and none has a run longer than 8.
  This agrees with the claim's exhaustive Colab result: 3,333 / 56,539 / 9,744 mirrors with a worst-support run of
  4 / 6 / 8.

**Checked by hand** [proved, line by line]:
- **The visiting lemma.** Between drops a register's frame only grows. A gate puts it on the frame of the x it
  joins. An x's frame stays inside v_b^⊥, because data never drops. So meeting x_a and then x_b forces
  u_a · v_b = g_b · u_a = 0. The three mixed cases (x then y, y then x, y then y) follow the same way.
- **The forward-call v-drops.** Take the registers whose frame is outside v_t^⊥. The span of their contents can
  fall only through a drop. At the start y_t is in that span; at the end it is not, because only x_t sits inside.
- **The R2 forcing.** All the readers' lines would have to lie in ⟨s, w_k⟩, which forces g_t ∈ {0, w_k}. Then the
  junk condition fails.

I did not rerun the claim's crossing counts (9 and 564 per replica) or its gate screen (2,453 illegal gates). Neither
is needed for the verdict.

## 2. Where the exclusion breaks

**(a) The claim's cost figures are not lower bounds, even in its own model.** `calc/v3order.py` uses the claim's
single-register model with lazy drops: a blocked x costs exactly 1 level, since a hyperplane section always has
codimension at most 1.
- Natural mirror: 21 to 24 drops per total per phase, against the claim's 27 to 31.
- Best twist: 15 to 16 per total, against 17.
- Per call (gather and clean): 364 and 248, against 462 and 272.

The gather order also ends at a frame of dimension 3 or 4. So the U-drop it feeds costs 3 or 4, not the 8 that the
claim charges on top.

**(b) A multi-register scheme is cheaper still** [computed, `calc/v5perx.py`]. It uses the claim's pricing: climbs
are free, a drop costs its codimension, extra scratch registers are free, gates join only arrays at one frame, data
never drops, and each x frame stays between ⟨u_b⟩ and v_b^⊥.
1. Copy each x_b once into its own register, at x_b's start frame ⟨u_b⟩. This needs no drop.
2. Climb all the copies to the window W. Build each total R_k there from the copies.
3. Make the U-drop: a copy of R_k goes from W to 0, which costs 9 instead of 8. Afterwards that copy is cancelled
   against R_k, and R_k is cancelled against the copies. Both steps are free.
4. Clean each copy with one 1-level drop, from W to W ∩ v_b^⊥, where x_b meets it.

The simulation checks every gate for equal frames and every climb for nesting. It checks every x frame against
⟨u_b⟩ ≤ F ≤ v_b^⊥, that each R_k holds exactly σ_k at its U-drop, and that every scratch register is exactly zero at
the end.

Result, the same for the natural mirror and for the claim's best twist: 112 cleaning drops plus 72 U-drop levels
(64 + 8) per call. So the extra cost beyond the two centres is **120 per call** for any second lines with u_b·v_b = 0.
That is ×3.9 below the claim's 462 and ×2.3 below its 272.

Grouping x's into run registers does not help with U8's coefficients. No two of the 112 x's have proportional
coefficient vectors across the 8 totals, so the best grouping is the singletons [computed]. Within this family of
schemes the cost is at least 0.5 per x, that is at least 56 per call [argued]: a run of L ports needs v-rank at least
1 + log2 L.

**(c) The twist exclusion rests on a figure that is not a bound.** The claim's one proved floor for the twists, the
run count with ×4 sharing, is 24 per call, inside the room of 34. It then says the verdict "rests on the drop sizes,
about 2.8 levels each". That figure comes from its greedy. The greedy drops eagerly, but lazy drops always cost 1
level each. So the claim's model does not establish the twist exclusion. What does keep the twists out is only that
no scheme anyone has built fits the room. My best is 120 per call.

## 3. What this means for the lead

- With the best scheme I found, C = 64 + 64 + 120 = 248 per call. That gives D = 1792 − 9C < 0.
- Even if the V-centre cost nothing in the forward reading (C = 184), D would be 136 per replica, below the one-sided
  320. So the scheme gains nothing.
- The mirrored second line therefore stays closed at first order. Its correct status is narrower than the claim's.
  The best cleaning anyone has costs 120 per call against a room of 34, and no lower bound above 34 is proved for
  every scheme.
  - Not covered by anyone: hybrids of chains and copies, sharing across totals, cleaning by way of the y arrays, and
    dirt carried into the next call.
  - Hosting the 112 extra registers is not counted either. It could only add cost.
- The fix for cx-mirror's NOTE is small:
  - replace "≥ 272", "462 / 272" and "deliberately optimistic" with "the best known scheme costs 120 per call (any
    second line); the claim's chain figures are upper bounds";
  - label the exclusion [argued], with a margin of ×3.5.

## Files

All in `checks/wht43/anatomy/v-cx-mirror/`:
- `STATUS.txt`: predictions (13:38 and 15:33), progress and outcomes; `RUNLOG.txt`.
- `calc/v1.py` → `v1.out` and `v1.json`: the port facts, run bounds, P7, the twist-proof step and the counts.
- `calc/v2price.py` → `v2price.out`: calibration, ceiling and break-even.
- `calc/v3order.py` → `v3order.out`: the claim's chain model with lazy drops.
- `calc/v4.py` → `v4.out`: the readers, R2, random twists and the anneal.
- `calc/v5perx.py` → `v5perx.out` and `v5perx.json`: the multi-register scheme, simulated exactly.

Data read: `bound-check/k-lower/out/U8X-parsed.json.gz`, and
`anatomy/cx-mirror/job-twistall/twistall.json` (the best twist's M).
