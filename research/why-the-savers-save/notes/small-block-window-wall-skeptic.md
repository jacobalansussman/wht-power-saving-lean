> Working note written by an AI agent during the round of 15 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# v-an2-smallnet: a skeptic's check of an2-smallnet (the small-block crossing network, proposal 3)

Anatomy round 2, 15 October 2026, 11:47 to 12:15 EDT (Mac clock). I tried to refute an2-smallnet's claim:
"excluded for the 1e-1 bar, with named gaps". Every program here is my own, written without running an2-smallnet's
code, and run on the Mac through a private memory-capping wrapper (each well under 1 GB and 2 minutes). No Colab was used, no third-party
code was run, nothing was posted.

Labels: **[computed]** a program in `calc/` and its output; **[checked]** I went through the argument line by line;
**[argued]** my own reasoning, not checked by anyone else.

## Verdict: not refuted

Every number I rebuilt matches, and every step marked [proved] in an2-smallnet's NOTE holds under the conditions it
states. The gaps it names are real, and I found no error that changes the verdict. I found four limits of scope that
the claim states too broadly or does not list among its gaps (section 3). None of them opens a route to 0.1.

## 1. What I rebuilt, and whether it matches

| item | my program | result |
|---|---|---|
| design table: n, t, k, junk per port, σ, ω, tightness, t(k−1)/n for E8, P112, E7, E6, D4–D8, A3–A7 | `calc/vdesign.py` → `vdesign.json` (my own E8 built from D8 roots plus half-integer vectors; P112 = E8 minus one orthogonal frame; E7 = lines orthogonal to one root) | **all 14 match**. σ = k − 1 for every design |
| early/late counts f and the bound n·f/(s−1) (n·f/(s−2) with twin labels) | `calc/vdesign.py`, my own enumeration of the hyperplanes of S(a) | **match**: E8 32 of 56 → 69.8; P112 32 of 52 → 70.3; E7 16 of 32 → 32.5; D5 16.0; D6 17.1; A4 8.0; A5 8.6 |
| label-class floor of NOTE 3.2 | `calc/vclass.py` → `vclass.json` (search for the largest induced union of cliques) | **match**: D4 12 of 12, A3 6 of 6 (both 1.00 per port), D5 12, A4 6, A5 6, D6 ≥ 12 |
| consistency of that floor with the hub's Lean-checked word | same program on P112, E8, E7 | P112 ≥ 14, E8 ≥ 14, E7 ≥ 12 (greedy, so lower bounds): far below the hub's 64 per call, so the hub's word does not contradict the floor |
| bank-exchange lemma (NOTE 2.5) | `calc/vshears.py` → `vshears_pm2.out`, coefficients ±1 and ±2 (wider than the claim's ±1) | **match**: no exchange with at most 2 calls per bank for 2, 3 or 4 banks; with 3 calls, 2 banks can and 3 cannot. The written proof is correct for any coefficients [checked] |
| frame calculus on raw arrays | `calc/vraw.py` → `vraw.json`, `vraw.out` (numpy from our own venv). My Φ_U is a *product of rank-one kernels* over an orthonormal basis; an2-smallnet builds one diagonal | **pass** at blocks 1 to 4: Φ_0 = I, Φ_whole = C^⊗m, Φ_u² = R_u, Φ_{u⊥}Φ_u⁻¹ = C^⊗m R_u for every odd u whose complement has an orthonormal basis; Φ independent of the basis; Φ_{U+V} = Φ_U Φ_V for U ⊥ V (blocks 2 to 4; block 1 has no such pair) |
| lone pair at 2m | `calc/vraw.py`, my own network (the data array descends once, then a free three-gate swap at frame 0) | **pass**: 2m exactly, outputs exact, blocks 3 to 5 |
| two pairs, one shared bridge, junk removed by side information | `calc/vraw.py`, my own network (every gate asserts one common frame, every move asserts nesting) | **pass**: Y_a = C^⊗m R_ua x_a + Φ(u_a^⊥) y_a exactly, the same for b, scratch ends at zero, X keeps x_a, and **Y_a's frame holds the partner's label u_b** (lemma 2.2 on arrays), blocks 3 to 5 |
| ceilings from the moment equation Σ c (r/m)^(1−a) = 2 | `calc/vceil.py` → `vceil.json` | **match**: calibration 0.0963 (closed form 1 − ln 3.5/ln 4); E7 at block 13: 0.0336 / 0.0857 / 0.1400; E7 at 19: 0.0154 / 0.0389 / 0.0631; P112 at 17: 0.0339 / 0.0745 / 0.1164; P112 at 25: 0.0152 / 0.0332 / 0.0516 |

The numpy build prints spurious "divide by zero / overflow in matmul" warnings (a known quirk of numpy with the
macOS Accelerate library); every comparison that matters returned True and none returned NaN.

## 2. The proofs, line by line [checked]

- **2.1** (side information leaves a source only while it is orthogonal to the target's port). Correct in a fixed
  gauge. Along a route with no descent the frames are nested in time order, and Y_a's frame lies in u_a^⊥ for its
  whole life.
- **2.2** (every target takes up its junk partners' labels). Correct as stated, for calls y += x.
- **2.3** (some source takes up the whole window). The counting is right. A target's early arrivals all have labels
  in Y_a's frame ∩ S(a) at the last early arrival, which is a proper subspace, so there are at most f of them. A source
  with at most one late copy has at least s − 1 early copies. I also checked that the count survives side information
  sent as *shared sums* of several sources, the hub's helper style: the bound counts (source, target) incidences, not
  copies. It rests on condition A and on the natural code, as the claim says.
- **2.4** (block ≥ 1 + 2σ). Correct, given 2.3. In fact the proof gives slightly more: V_s′ is orthogonal to the whole
  window W_s, not only to V_s.
- **2.5** (some bank is a source twice). Correct, and confirmed by search.
- **2.6** (a backward call and the next forward call cannot share a window). Correct under condition A in both
  directions. One wording point: "each call-3 bridge serves one port" assumes a bridge descends only once. A bridge that
  descends again between its targets pays at least one more level per extra port, so the call is still not amortised.
- **3.2** (the class floor). Correct for bridges that each gather, descend and then scatter (rank-one pieces).
- **Mirror symmetry for x −= y calls.** Reversing time and complementing frames does map legal networks to legal
  networks. But in a backward call the constraint is carried by the *erasure*, not the creation, of a copy. So every
  mirrored statement needs condition A. The claim states condition A for 2.6 and 4; I read its general sentence "every
  statement holds for x −= y calls" as including it.

## 3. Limits of scope the claim does not list among its gaps

1. **"D4 and A3 save nothing under *any* centre"** is proved only for single-layer bridges (NOTE 3.2's rank-one
   pieces: one array gathers, descends, then scatters). The descent-counting argument of NOTE 3.1 does not carry over
   by itself to networks with several layers, or to bridges that go up and down more than once, for two reasons:
   - One descent can take an array out of several label classes at once. For D4, the frame ⟨111⟩ of the window is
     left by one descent for all three classes.
   - One descent can both make one class exit and merge two others. In D4's window, ⟨110, 011⟩ → ⟨101⟩ does this.

   A rank-flow count does seem to restore the floor [argued]. Each descent lowers the rank of the restricted label
   functionals by at most one. The span of class content that reaches each lower-rank stage grows by at most one per
   descent into it. And a frame on which two classes agree can be reached only by descents, never by climbs or gates.
   I found no D4 network below 12 levels per call. But the statement for every network is not proved, in the NOTE or
   here.
2. **σ and the floor 1 + 2σ are for each design's natural code** (every partner at ±1/2). NOTE 2.3 says that a code
   with fewer partners per port is "not explored", but sections 6 and 9 and the verdict call 1 + 2σ "the proved floor"
   without that condition.
3. **The floor 1 + 2σ needs the twice-source bank's two calls both to be amortised.** In a two-bank swap this is
   forced: a plain call costs about 1 per port, and 1 + 2c ≥ 2 for every design here (c ≥ 0.533). In a long chain, the
   regime the ceilings use (2 calls per pair), one plain call is cheap, so 2.4 alone does not force the floor. A second
   argument does seem to close it [argued]. By the mirror of 2.2 under condition A, an X array that took up W_s and is
   then the target of an amortised backward call forces that call's partner labels to be orthogonal to W_s. That again
   gives m ≥ 1 + 2σ. The claim does not make this argument; its own multi-bank remark is marked [argued].
4. **The ceilings depend on the grouping assumption, and the claim calls it "ideal".** The model puts every data
   array's levels into children of rank σ. If an array's levels came in a single child, the same equation gives 0.218
   for E7 at block 13, 0.293 for P112 at 17 and 0.247 for P112 at 25, all one-sided (`vceil.json`, sensitivity lines).
   So σ-blocks are a restriction, not ideal grouping. I argue that, at the floor block, the restriction holds for the
   source arrays that took up their window in call 1, which by NOTE 2.3's count is at least about half of them
   [argued]. Take such an X array that sends side information in its second source call. Its receivers'
   ports span σ dimensions and lie in the complement of its frame at its first send there. So that frame has
   dimension at most m − σ = 1 + σ. Its call-1 erasures (condition A) came earlier at a frame containing W_1, which has
   dimension 1 + σ. So it has gates at dimension exactly 1 + σ, and no child can exceed σ. For the other sources, for
   targets and for idle levels, the restriction is not proved.

None of these is a counterexample. Item 4 matters most: of the claim's assumptions, it is the one that the "below 0.09"
statement most depends on.

## 4. Predictions (STATUS.txt, 11:50) against outcomes

| prediction | outcome |
|---|---|
| Q1 (80 %): the raw-array checks reproduce | **held** |
| Q2 (75 %): the design table reproduces | **held** (14 of 14) |
| Q3 (60 %): the f counts reproduce | **held** |
| Q4 (55 %): a unit or break-even problem in the class floor, or its scope overstated | **half held.** The units are consistent: a dropped level costs at least one move, so the floor is conservative, and the hub's 64 is above my P112 floor. The scope is overstated (section 3, item 1) |
| Q5 (65 %): an unnamed gap that weakens a [proved] label | **held, but no step fails.** Items 2 and 3 of section 3 (the natural code; both calls amortised). The steps themselves hold under their stated conditions |
| Q6 (70 %): the ceiling table reproduces | **held** |
| Q7 (45 %): I refute the claim as stated | **did not happen** |

## Files

All in `checks/wht43/anatomy/v-an2-smallnet/`: `STATUS.txt`, `RUNLOG.txt`,
this `NOTE.md`, and in `calc/`: `vdesign.py` / `vdesign.json`, `vclass.py` / `vclass.json`, `vshears.py` /
`vshears_pm2.out` / `vshears_pm2.json`, `vraw.py` / `vraw.out` / `vraw.json`, `vceil.py` / `vceil.json`.
