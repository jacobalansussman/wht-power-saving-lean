<!-- header: frame-packing-two-arrays -->
# The skeptic's step-by-step check of the slack round

**About this file.** This is the working proof file as an AI agent wrote it during the round of 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent. Its opening lines say that no skeptic had read it; that was true when it was written. What happened afterwards is below, and where a later skeptic's verdict differs from this file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean. File names inside it refer to the programs and outputs now in `../programs/` and `../data/`; [README.md](README.md) maps the working names.

**Status now.** This is the skeptic's own record. Its verdict: both tracks upheld, with the corrections listed in [ledger-for-every-unit.md](ledger-for-every-unit.md) and [common-frame-accounting.md](common-frame-accounting.md).

---

# v-slack: step-by-step check of sl-algebra and sl-combinatorics (15 Oct 2026)

Skeptic's record. Labels: [checked] (I re-derived the step), [computed] (my own script in this folder, exact
arithmetic, none of their code imported), [read] (from a named file), [argued].

## 1. Definitions: well defined for every unit?

- **Terms and frames.** Fact 0 (sl-alg) and the "formal identity" model (sl-comb) give each content a unique
  expansion over group elements, so terms, frames (g^-1 L0) and positions (cosets Bg) are defined for every unit
  [checked]. One point to keep straight: a scalar *letter* is a coefficient, while a scalar *matrix* that arises as
  a product of letters (e.g. -I = XZXZ) is a group element acting as (-1)^f on f columns. Both notes treat scalars
  as coefficients and quotient by them in labels; nothing below depends on it.
- **Frame facts** [computed, `vframes.py`, from actual matrices at m = 1, 2, 3]: free letters fix L0; dir z acts
  as the transvection (a|b) -> (a + (z.b) z | b); every paid letter moves every Lagrangian by at most 1 and L0 by
  exactly 1; d(L0, Kbar L0) = m; 3, 15, 135 Lagrangians. The dual polar graph has triangles (a1 = 1), so lateral
  steps exist and path-length parity is not fixed.
- **Labels, gauges, coherence** (sl-comb Def. 1.2, 3.2; sl-alg Def. 2): well defined. sl-alg's "twists unique
  only up to one common right factor" is slightly too strong (each connected component of the co-occurrence graph
  of inputs has its own factor); it does not affect any result.

## 2. sl-algebra

- **Theorem 1** [checked]. Per-array path lengths sum to M; Delta = sum overlaps - sum detours; the reference form
  M = Wm - 2N_o + 2L_o + E_o with E_o >= 0, equality iff o is on a geodesic from every a_s to a*_s. Correct. Note
  2L_o = 2(steps toward o) + (lateral steps), so L_o can be a half-integer; this is why parity is not forced
  inside the class.
- **Prop. 2 (median bound)** [checked]. Delta <= 2 N_o for every o. The sink version also holds directly:
  dim(a_r cap e_r) <= d(o, a*_r) + d(o, e_r), summed, gives Delta <= 2 sum_s d(o, a*_s).
- **Prop. 3, Cor. 1, Cor. 2** [checked]. A non-monomial gate in a common-frame network joins rows with one common
  element (otherwise an output row has two elements), so re-matching at it is legitimate; with W = 2 the pairing
  a_s -> a*_s is always available; M >= 2m.
- **Cor. 2' (the step sl-alg asked a skeptic to check)** [checked]. W = 2 arrays in total (no scratch). The state
  matrix is invertible at all times (product of invertibles), so no row is zero and no column is zero. Left
  multiplication preserves each row's twist g1^-1 g2; a non-trivial combination of two rows must have equal
  elements on every input both rows hold (two distinct group elements never cancel), so every row holding both
  inputs inside one joint episode has the same twist tau. With twists (e, tau^-1) the episode is a common-frame
  network; the opening and closing gates are at a common twisted frame; each array travels from F0 to F1, so the
  episode costs >= d1 + d2 (raw distances travelled by the two inputs, tau an isometry). Separate periods are
  trivial; triangle inequality gives M >= 2m. **Correct.** It does not extend to scratch arrays (a scratch row
  could carry a second twist), and sl-alg does not claim it does.
- **Theorem 5** [checked]. Some history has product exactly pT (group elements are a basis); its frame moves <= 1
  per paid move; Sum mu_j >= Wm. Correct (and, as stated, it excludes nothing by itself).
- **Prop. 4 and the cut numbers** [computed, `vframes.py`, `vspectral.py`]. c_2 = 3/2 by brute force over all
  16,383 cuts; the spectral ceiling over all cuts equals the point cuts, 2(1 - 2^-m): 3/2, 7/4, 15/8 at m = 2, 3, 4.
- **Section 2a (widened game)** [checked]: row support at most doubles per kernel and is preserved by monomials,
  so d(B', B'T) = m; the four facts used carry over.
- **Stand-in explanation** (calibration table). The conclusion "not an engine unit" is right, but the stated reason
  ("free mixing is non-scalar; any accounting needs scalar mixing") is imprecise: the engine's own mixes
  sum mu_s b_s c_s are not scalar either. The real reason is that the stand-in is correct in one representation
  only. My check [computed, `vstandin.py`]: the 2 in 3 is exact in V(x)V; all its moves act on the two inputs
  alone (depth 1), so each output coefficient lies in span(B) + span(B R B), whose group elements all sit at frame
  distance <= 1 from L0, while K = D1 D0 sits at distance 2. So **no lift of it to the engine game is a formal
  identity**, independent of Theorem 2W.

## 3. sl-combinatorics

- **Lemma 1.1** [checked; (d) and (e) also computed]. The stabiliser argument is right: the phase letters give
  every symmetric S (rank-one zz^T span all symmetric matrices over F_2), address maps give GL(m,2).
- **Theorem 2.1, Cor. 2.2, 2.3** [checked]. Pure bookkeeping; correct (and tautological, as the note itself says).
- **Cor. 2.4** [checked], with one omission: the parity conclusion M = Wm (mod 2) also needs I0 = 0 (no move on an
  empty register); the displayed identity drops I0.
- **Lemma 2.5** [checked]: geodesic points between opposite L0, L1 are F = (F cap L0) + (F cap L1); gives F_C.
- **Prop. 3.3, Example 3.4** [computed, `vctl.py`]: the order-3 control word is an exact formal identity at
  blocks 3 and 4 with M = Wm; s is not in Paulis*scalars, so no gauge is coherent.
- **Lemma 4.1, 4.2** [checked]. Classes never cancel each other; a mix result in class lambda lies in the
  A-submodule of the same-position summands; a move on a single-class live register keeps v and moves one step;
  liveness is inherited backwards along the data flow; length counting works for any finite-length ring A (so
  whether A is read as C[Pauli group] or as M_n(C) does not matter). Correct.
- **Theorem 4.3** [checked]. The edges of H are distance-1 position pairs; Hu's two-commodity theorem (cut
  condition suffices for a fractional flow, undirected) is cited correctly; each flow path has length >= m.
  **Correct, with any number of scratch arrays.** It and sl-alg's Cor. 2' are incomparable: 4.3 needs one global
  gauge but allows multi-term contents and scratch; 2' allows a new twist in each episode but needs single-term
  contents and no scratch. So sl-comb's remark that 4.3 "generalises MAP 1.3's single-term contents" is not quite
  right: family S with twists that change between episodes is covered by Cor. 2', not by 4.3.
- **Cor. 4.4, Thm 4.5** [checked]: Papernov's classification and the O(log k) flow-cut gap are cited correctly.
- **Remark 4.6**: "c_3 >= 1.75, nothing better found" can be sharpened to c_3 = 7/4 exactly (spectral ceiling).

## 4. Where the claims disagree, and who is right

1. **Are our Lean savers coherent?** sl-comb: "E8, T74 and U8X0HQ ... coherent by construction of the hub's frame
   templates [read]; they sit under Corollary 2.4". sl-alg: our five-stage units (and the E8 unit) and our U8 Lean
   words are not common-frame in their stated frames. **sl-alg is right.**
   - The E8 unit is not a hub template. It is our fourth result, which "keeps the bridged five-stage word" with
     "One helper set for all stages. A re-framing lemma lets an invocation start with helpers that still hold old
     values" [read: the release README, section 6]. sl-alg's "[inferred]" is now [read].
   - The U8 route: "`gdirtyE` reframes per slot" [read: the U8 Lean route's `Work/Gauge/Inv.lean`,
     line 13].
   - The five-stage units also carry scratch copies of the totals that move and are consumed (cst moves per
     invocation), which is branching (B != 0), so (H2) fails for them in any gauge.
   - Whether some other global gauge makes them coherent was settled by nobody. My expectation [argued]: no. A
     re-framed helper meets bank arrays in several stages with factors that differ by k directions, and the bank
     arrays are tied together by bank-to-bank additions (7,260 per invocation [read: proveout NOTE]), so the
     relative-gauge cocycle has holonomy. By sl-comb's own definition that would make our Lean savers "saving
     units incoherent in every gauge", i.e. a "second mechanism" in name, although the deficit they carry is the
     bank's crossing slack minus copy costs.
2. **Parity.** sl-alg says nothing about parity (the brief's ADDED item). sl-comb proves parity only under
   (H1)-(H3) (and needs I0 = 0), and checks it only on the stand-ins. **Both missed a known saver with odd
   deficit:** our first published unit (m = 80, h = 16, v = 560, R = 7,122, W = 9,362) has
   D = 4v - 5 cst = 2,240 - 1,205 = 1,035 and M = 747,925 of Wm = 748,960 [computed from the README table; the
   table's own numbers agree]. An odd deficit is impossible for any unit satisfying M = Wm - 2N + 2L. So
   sl-comb's "No known unit that saves breaks the identity" and sl-alg's "No known engine unit breaks the
   identity's deficit form" are false as stated. The odd part comes from a cost (15 copies of 15 moves and one of
   16), not from a new way to save, so the claims' substantive point (no second *saving* mechanism is known)
   survives.

## 5. Parity, stated exactly (my statement)

- If a unit has a gauge with (H1) one label per register at all times, (H2) B = 0, (H3) X = Z = 0, and no idle
  moves, then Wm - M = 2(N - D) is even [sl-comb Cor. 2.4 plus I0].
- Drop any one hypothesis and odd deficits occur: a lateral step (Z = 1) gives the exact 1-array 2-move unit at
  block 1, b (dir S dir) b' = dir [computed, `vctl.py`]; sl-comb's packing gadget (M = Wm + 1) and Pauli twin
  (M = m + 1) do it with B; our first published unit does it with its copies (D = 1,035).
- Inside the common-frame class alone (H1), Theorem 1 gives Delta = sum overlaps - sum detours, with no parity
  constraint (the graph has triangles).
- **Consequence for the brief's ADDED question.** Parity excludes "save one move" units (W in Wm - 1) and the
  full-sharing U_m (odd deficit when Wm is even) **only inside the narrow (H1)-(H3) class**. A unit that saves two
  moves plus one lateral detour saves one; so parity cannot prune the tier-2 search or U_m in general.

## 6. U_m under the identity (the brief's last question; neither note answers it)

U_m has W = 2^(m-2), M = 2^m - 1, so Delta/W = m - 4 + 2^(2-m). Under (H1)-(H3) the average source level must be at
least (m - 4 + 2^(2-m))/2 (and parity forbids U_m exactly there, section 5). In the wider common-frame class
(Theorem 1) the average overlap dim(a_r cap e_r) must be at least m - 4 + 2^(2-m): each array's path is about 4
steps long, so arrays act as relays and every input's data crosses about m/4 arrays. Nothing proved excludes that:
cuts give only c_m W < 2W against U_m's 4W; Theorem R's floor 2.9447 W is below 4W; the O(log W) flow-cut bound
has an unknown constant. It does need sources spread so that every frame o has sum_s d(o, a_s) >= W(m-4)/2
(median bound): label dimensions that grow with m, unlike the lines of every known saver. Possible or not is open.

## 7. Ceilings against known units and proved bounds [computed, `vunits.py`]

Every number quoted by either claim recomputes: Theorem R ceilings 0.0628, 0.1714, 0.3191, 0.3919 (blocks 3-6),
0.4806 at block 8; thresholds m <= 6 / 29 (one overlap), 9 / 51 (two), 16 / 90 (four); the needs tables; the
two-array targets' savings; 169/57 = 2.965 < 3. No ceiling either claim states lies below a known unit: E8's share
D/(Wm) is 0.0021 (cap 1/m = 0.022), the m = 110 unit 0.0019 (cap 0.0091), the m = 80 unit 0.0014 (cap 0.0125);
published exponents (about 7.8e-4 to 1.6e-3) are below the one-overlap move-count caps (5.9e-3 at m = 45, 1.9e-3
at m = 110, 3.2e-3 at m = 73). The stand-ins save 0.415 at block 2, above the engine cap 0, but they are not
engine units (section 2).
