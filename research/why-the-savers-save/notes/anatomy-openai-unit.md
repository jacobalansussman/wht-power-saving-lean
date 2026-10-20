> Working note written by an AI agent during the round of 15 October 2026; kept as written, with private paths
> and machine names removed. Where it disagrees with this section's README.md, the README wins: it
> follows the later skeptic verdicts.

# an-openai: why OpenAI's Fourier unit saves, and what caps it

Track an-openai of the anatomy round, 15 October 2026, about 08:10 to 08:35 EDT. I read OpenAI's three preprints
first-hand (sources at the end) and checked the parts that carry the saving with my own code (section 7). No
skeptic has read this note. Nothing was posted. Labels as usual: [proved] (a written proof exists), [computed] (a
program and its output), [read] (stated in a source, not rechecked), [argued] (a short argument, checked by me only),
[guess]. **Derived here** marks reasoning added in this note.

## In short

1. **The object is one gigantic engine-game unit.** OpenAI's explicit Fourier paper uses a single finite
   network: W = 1,873,807,244,643,542,670,000 arrays at block m = 10^6, transformed with M = Wm − 6,871,402,692,000,000
   paid moves. Its saving is 2.65e-13, and 2.11e-13 after the paper pads W up to 2^71. It is the same kind of object as
   our block-3 to block-5 targets ("W in M at block m"), in the engine game, because the 2×2 matrix it powers,
   C = ½[[1+i, 1−i], [1−i, 1+i]], equals e^{iπ/4}·S⁻¹HS⁻¹ with S = diag(1, i) a free phase [computed].
2. **The saving comes from crossing routes, not from any clever move.** Every array pays exactly the change of
   its own frame label. Every data route needs the change between its two end labels. The network swaps the data of
   two banks of arrays. So routes and arrays differ, and the two sums can differ. For each swapped pair the input
   label is a line ⟨u⟩ and the output label is its complement u^⊥. The route still gets the whole transform, because
   the square of a paid kernel along u is a free address shift (C_u² = R_u). Each of the 8.46e15 bank arrays saves
   exactly one move. Everything else in the network exists to make that swap possible.
3. **The swap needs bridges, and bridges lose.** An array whose frame always contains u can never share a gate with
   one whose frame never does. So data must cross from one bank to the other through bridge arrays (the "centre"),
   and each bridge must drop a whole frame factor once per stage. That costs 200 moves per bridge array. The bridge is
   shared by 161,700 pairs, because the transfer it carries is a low-rank matrix: (|S∩T| − 1)/2 on 3-subsets of
   100 points (Alon's 1998 triple construction).
4. **Junk from the shared bridge is cleaned by side carriers.** The low-rank bridge also delivers unwanted terms
   between pairs. They are cancelled by one private carrier array per ordered pair, and those carriers cost nothing
   extra because the two pairs' frame lines are orthogonal. The carriers break even, but they are 99.9995 % of all
   arrays. So they dilute the saving by a factor of about 220,000.
5. **This is how a real saver escapes our near misses' failure modes.** No target has a single route in: each
   target is reached through 101 bridge arrays and about 147,731 side carriers per stage. The route that delivers the
   target is low rank (a zero divisor) and the route that removes junk has no diagonal. Their sum is the identity,
   so neither has to be both invertible and a zero divisor. Junk does not travel with the target in a fixed ratio,
   because it has private carriers. Scratch arrays are not required to be clean. They are themselves arrays of the
   unit, and their arbitrary contents are cancelled by mirror-image subtraction (the Bennett / Buhrman–Cleve–Koucký–
   Loff–Speelman restoration).
6. **What caps it is structural, plus parameter choice and pricing.** The saving obeys an exact identity: twice the
   total input-label dimension minus twice the dimension lost by bridges. With line labels it can never exceed 1/m
   per array, whatever else is done [argued, derived here]. The paper's own parameters are far from the best in its
   own family. Ground set h = 25 instead of 100 gives 4.18e-10, 1,600 times more [computed]. Pricing each frame step
   as one recursive call, as the hub does, adds another factor of 80 to 360 [computed].
7. **How it was found:** by structure, not search. The parts are Alon's two-field triple identities, restoration of
   dirty scratch, and an orthogonal-frame calculus on F_2^m. The constants are round numbers (h = 100) chosen so that
   the inequality v > 3h(h+1) holds with room. The companion existence paper ("Finite tensor savings") proves that
   some finite win exists by contradiction, but gives no unit and no number.
8. **Knobs, ranked by how far they could move it.** Fewer break-even carriers per saved dimension; terminal slack on
   more arrays, including carriers; higher-rank input labels; smaller blocks; pricing whole edges as one recursive
   call; ground-set size. Each is stopped by the bridge: the loss is roughly the bridge rank times the span of the
   labels it carries, and its junk needs carriers. The hub's improvements, as far as their names show, are exactly
   turns of these knobs (section 4).

## Words used here

- **Array, block m, paid move, free step, unit "W in M", saving** as in `../../exist2/NOTE.md`:
  - The saving is ln(mW/M)/ln m.
  - A paid move applies one direction of the transform to a whole array, in every digit column, and becomes a
    recursive call.
- **Label (frame).** OpenAI attaches to every gate a binary subspace U of F_2^m, shared by all arrays the gate
  touches, and represents each array in the frame Φ_U = J·diag(i^{q_U(x)})·J. Here J is the normalised Walsh matrix
  and q_U(x) is the Hamming weight of the projection of x onto U, taken mod 4. Then:
  - Φ_0 = I, and Φ of the whole space is C^{⊗m}.
  - Φ is multiplicative on orthogonal sums.
  - Moving between nested labels U ⊂ V costs exactly dim V − dim U paid moves, one per vector of an orthonormal basis
    of the part of V orthogonal to U.
  - Gates mix arrays pointwise. They commute with a common frame, so they are free.
- **Role, route.** A role is one physical array, with its own source and sink. A route is the path the data takes
  from an input role to the output role it ends on.
- **Bank, centre, side carrier.**
  - The bank is the two families of 4.2e15 data arrays X_a, Y_a, indexed by triples of 3-subsets of {1..100}.
  - The centre is 101 bridge arrays per invocation.
  - A side carrier is one array per ordered pair of "neighbouring" 3-subsets per invocation.

## 1. The saving object

**Game.** It is our engine game, with the paid 2×2 matrix C instead of H [computed: SCS = e^{iπ/4}H, see
`calc/`]. It uses these free steps:
- pointwise mixing of any set of arrays, with coefficients in Z[i, 1/2];
- invertible F_2-linear address maps on the m-bit block, applied the same way in every digit column (each directional
  kernel C_z is turned into a coordinate kernel this way);
- address shifts (C_z⁻¹ = R_z·C_z, and terminal corrections R_u);
- Pauli Z phases on 27 coordinates, one global scalar i^{27k}, signs, and relabelling roles.

It halves: the scatter coefficient is 1/2, and the side coefficients are ±1/2.

It holds for every number f of digit columns and with any spectator bits (Proposition 2.4 of the explicit paper)
[read; my checks agree].

**Target.** C^{⊗m} on every one of the W arrays, including every scratch array. Via the free phases S this is the
Walsh–Hadamard transform H^{⊗m} on every array.

**Arrays and moves** [read; recomputed exactly in `calc/family.py`]:

| part | arrays | moves each, relative to plain m | total effect |
|---|---|---|---|
| bank X_a, Y_a (a ∈ T³, v = C(100,3) = 161,700, N = v³) | 2N = 8.456e15 | −1 | gain 2N = 8.456e15 |
| centre (101 per invocation, I = 3v² invocations) | 7.92e12 | +200 | loss 2L = 1.585e15 |
| side carriers (v·d per invocation, d = 147,731 neighbours) | 1.874e21 | 0 | none (dilution only) |
| **unit** | W = 1.8738e21 | | Δ = 6.871e15; saving 2.65e-13 |

The recursion then uses this unit on W_* = 2^71 arrays, so that recursive batches are always full. The padding
roles are plain, so the exponent is θ = log_m(m − Δ/2^71) and 1 − θ = 2.106e-13 [read; recomputed].

**Other OpenAI objects.**
- The integer-multiplication paper (23 Sept) has the same complex network. It also has a *bit* network: XOR gates,
  rational labels with the form I − J/9, and shear frames (H, D) ↦ (H + MD, D) on tape address fields. That network
  saves N − 2L, one per pair. The first edge on X_a goes from −P_U to +P_U, which costs one rank yet makes two units
  of progress. Relative saving 339/22,587,335,000,000 = 1.50e-11. It lives in the tape model, not in our game.
- The "Finite tensor savings" paper has no explicit object: it proves only that some finite win exists (section 4).

## 2. Where the saving comes from, move by move

**2.1 The exact accounting** [read; brute-force checked at h = 7, 8, 9]. Labels along each array only grow, except on
the centre's single descent per stage. So the paid moves are

  M = Σ over arrays (sink label dim − source label dim) + 2 × (total dimension lost on descents) = Wm − 2N + 2L.

Gates contribute nothing: every array at a gate has the same label before and after.

**2.2 The pair of crossing routes** [read; identity checked at m = 2..5].
- The labels are:
  - source X_a = U_a = ⟨u_a⟩, with u_a = t_{a1}⊗t_{a2}⊗t_{a3}, of weight 27 and norm 1;
  - sink X_a = whole space;
  - source Y_a = 0;
  - sink Y_a = U_a^⊥.
- The scalar network sends X_a's data to Y_a and Y_a's data to −X_a.
- The route X_a → Y_a therefore gets Φ(U_a^⊥)·Φ(U_a)⁻¹ = C^{⊗m}·Φ(U_a)⁻² = C^{⊗m}·R_u: the full transform up to a free
  shift, because C_u² = R_u.
- The route Y_a → X_a gets 0 → whole space: the full transform.
- The arrays pay U → whole (m − 1) and 0 → U^⊥ (m − 1). Two full transforms cost 2m − 2.

That is the only place where any move is saved.

Two facts make this the *only* possible source in this model [argued, derived here]:
- **Sink labels are forced.** For a route's frame ratio to be the full transform up to free steps, the sink label must
  be the complement of the source label. Matching the quadratic forms mod 4 forces P_V = P_{U^⊥}. So the saving is
  exactly 2·Σ(source label dims) − 2·(descents).
- **A role whose data stays on its own array gains nothing.** Its array must go from U to U^⊥, and those are not
  nested, so it must descend by at least dim U.

Hence: **saving needs data to change arrays.** Relabelling roles is free but does not count: the crossing has to
happen through gates.

**2.3 Why a bridge is needed** [argued, derived here].
- X_a's labels always contain u_a (they start there and never descend).
- Y_a's labels always lie inside u_a^⊥.
- So X_a and Y_a can never meet at a gate, which needs one common label.
- With only the two arrays (W = 2), every arrangement of meeting gates costs at least 2m. This is the lattice
  distance argument: d(U, Z) + d(Z, whole) + d(0, Z) + d(Z, U^⊥) ≥ 2m.
- So the data of X_a must reach Y_a through a third array. That array must collect from a gate containing u_a and
  deliver to a gate avoiding it, so it must descend.

**2.4 The centre amortises the bridge** [read].
- One invocation fixes two of the three triple coordinates and handles v = 161,700 pairs at once.
- The 101 centre arrays gather G: c_j = Σ_{T∋j} x_T and c_* = Σ_T x_T, then scatter R: (Σ_{j∈S} c_j − c_*)/2.
- The product RG has entries (|S∩T| − 1)/2:
  - 1 on the diagonal;
  - 0 when the triples meet in one point (the only non-orthogonal pairs, t_S·t_T = 1);
  - ±1/2 on the rest.
- Rank 101 serves 161,700 diagonal transfers. Each centre array drops the current frame factor P⊗F, of dimension
  h = 100, once per stage: it goes up 1→3, down 3→4 and up again 4→6.
- Per invocation the gain is 2v/3 and the loss 2h(h+1). So the network saves iff v > 3h(h+1), which needs h ≥ 22
  [computed].

**2.5 Side carriers clean the junk on orthogonal pairs** [read; brute-force checked].
- The off-diagonal ±1/2 entries of RG fall on pairs with even intersection.
- Each such ordered pair (S, T) gets a carrier A_ST:
  - It copies x_T at a gate with label (B⊗F) ⊥ (P⊗⟨t_T⟩).
  - It injects −(|S∩T| − 1)/2·A_ST into y_S at a gate with label (B⊗F) ⊥ (P⊗t_S^⊥).
- Because t_T ⊥ t_S, the first label lies inside the second. The carrier only climbs, so it costs exactly m, like a
  plain array: no loss.
- RG + JV = I.
- That orthogonality is the whole reason the pairs needing correction are exactly the cheap ones: Alon's F_2 and Q
  structures are complementary.

**2.6 Dirty scratch** [read].
- Rows 0–1 subtract the effect of the old centre and carrier contents from the target.
- Rows 4–5 add the effect of the new contents.
- Rows 6–7 restore them.
- So every scratch array ends with its own C^{⊗m} for arbitrary initial values. This is what lets all W roles be
  arbitrary recursive fibres (Remark 2.7 of the explicit paper).
- A restoring centre must touch the two banks alternately (target, source, target, source, or the reverse). That
  forces at least one bank-X → bank-Y step, and hence one descent, per centre array per stage. The signed swap needs
  three shear stages. Together these fix the loss at 2 × 101 × 100 moves per invocation, in each of the three stages
  [argued].

**2.7 How it escapes the failure modes of our near misses** [argued, derived here]:

| our failure mode | what OpenAI's unit does instead |
|---|---|
| unit / zero-divisor dilemma (absorb4: one route must be invertible and kill leftovers) | the target is the sum of two routes, each singular: low-rank centre RG (rank 101 on 161,700) plus off-diagonal side routes JV. Leftovers are cancelled by mirror subtraction, not by a zero-divisor factor |
| pattern identity (a2-nag: junk travels with the target in a fixed ratio) | junk does travel with the target through the centre in fixed ratios (|S∩T|−1)/2, but every junk term has its own private carrier that never touches the target's data |
| one route per target | each y_S is touched by 2 centre gates (101 arrays) and 2 side gates (about 147,731 carriers) per stage |
| cleanliness of depth-3 carriers (o3) | carriers are dirty at full depth. They need not be clean, because they are themselves arrays of the unit and are restored exactly |
| needs halving, address maps or phases (MAP §1.3) | uses all three |

The price of this escape is **array count**: every saved dimension is accompanied by about 220,000 break-even
arrays. At W = 2 or 3 there is no room for a bridge at all (2.3).

## 3. What caps it

**3.1 Structural ceiling** [argued, derived here].
- Δ = 2·Σ(source label dims) − 2·(descents), and every array of the unit pays m.
- So the relative saving Δ/(Wm) is at most 2·(average source-label dimension)/m. With line labels, at most one
  array in two can carry one, so Δ/(Wm) ≤ 1/m.
- To reach the bar (saving 0.1) the average source-label dimension per array must be at least m(1 − m^{−0.1})/2:
  0.75 at m = 8, 1.9 at m = 16, 18 at m = 100, and 374,000 at OpenAI's m = 10^6. OpenAI's average is 4.5e-6.

**3.2 Where the paper's figure goes** [computed, `calc/caps.py`]:

| h | m | bank / W | loss / gain | saving | without side carriers | ceiling 2N/(Wm) |
|---|---|---|---|---|---|---|
| 22 | 10,648 | 6.5e-4 | 0.986 | 9.4e-11 | 1.4e-7 | 1.0e-5 |
| 25 | 15,625 | 4.2e-4 | 0.848 | 4.18e-10 | 9.9e-7 | 6.6e-6 |
| 100 (paper) | 10^6 | 4.5e-6 | 0.187 | 2.65e-13 | 5.9e-8 | 7.2e-8 |

At h = 100, the side carriers cost a factor of 220,000 and the bridge loss 19 %. The block size m = 10^6 costs
everything else: per-array slack of 1 in 10^6, then divided by ln m = 13.8.

**3.3 The parameter choice.** In its own family (3-subsets, 3 stages, 101 centres, ordered carriers), the best ground
set is h = 25, with saving 4.18e-10, 1,600 times the paper's [computed, `calc/family.py`]. Below h = 22 the bridge
loses more than the crossing gains. The paper evidently chose h = 100 for a comfortable inequality, not for size.

**3.4 The pricing.**
- OpenAI charges every residual dimension as a separate recursive call of size f, so that all children have equal
  size.
- Padding W to 2^71 loses a further 21 %.
- The hub instead charges each whole edge as one child of size r·f. This is legitimate because Φ of an r-dimensional
  orthogonal residual is address-conjugate to C^{⊗r}. The rate then solves Σ_edges (r_e/m)^τ = W.
- On OpenAI's own network this gives 9.55e-11 at h = 100 (360 times) and 3.26e-8 at h = 26 [computed,
  `calc/moment_price.py`]. Even then it is 50,000 times below the hub's 1.68e-3.

**3.5 No certificate is involved.** OpenAI prove no lower bound for their family, and none of our certificates is
anywhere near tight at block 10^6. The cap is the design: line labels, about 3·d full carriers per pair, the bridge
loss, and m = h³.

## 4. How it was found, and the lineage

**OpenAI's own sequence** [read]:
1. *Integer multiplication below n log n* (23 Sept). Both networks appear here: bit and complex, h = 100, triples.
   The parts it cites are Alon 1998 (the triple Gram identities over F_2 and Q, "ground-set size changed to
   h = 100") and Buhrman–Cleve–Koucký–Loff–Speelman 2014 (restoring dirty scratch). Relative savings: bit 1.50e-11,
   complex 3.67e-12, both rounded to an exponent 1 − 2^{−50}. The tape model then dilutes the result to κ = 2^{−182}.
   - The complex network saves twice per pair; the bit network saves once. But the complex network needs about 10
     times more carriers (147,731 even-intersection neighbours per triple, against 13,968 meeting in one point), so
     its relative saving is about 4 times smaller.
2. *An explicit power saving for the exact DFT* (25 Sept). It reuses the complex network unchanged and pads to 2^71.
   It adds the exact-width Fourier compiler (Newton / triangular Toeplitz), synchronisation over coprime prime axes,
   and chirp convolution for every length. These cost only (log log n)^{4−θ}, not exponent. Saving 1 − θ = 2.1e-13,
   stated as 1e-13.
3. *Finite tensor savings and exact Fourier circuits* (25 Sept). This is an independent existence route in a
   different game, the "finite win": one array of q^b entries, calls to A on q coordinates, monomials free.
   - The amplification lemma recurses over binomially weighted sectors.
   - If no win existed, a price function with tensor and direct-sum laws would exist.
   - The contradiction: a reflection priced at d·w by basis exchanges is also bounded by (3d/4)·w + O(w), via sparse
     symmetric perturbations, because a dense 2×2 pair has price at most 3/2.
   - The only concrete saving mechanism inside is that dense-pair lemma. B(u)⊗B(u) is reached from a cyclic 4-shear
     loop (closed by feedback, which amortises its initial state) plus 2 shears and 4 basis exchanges: 6 shears where
     4 dense pairs would naively cost 8 (the 3/4).
   - No size, saving or explicit unit can be extracted. Separation and compactness are non-constructive.

So the method was theory and structure. The README says each result took about three hours of an internal model
[read]. Nothing was searched or optimised: the constants are round, and h is 4 times the optimum.

**What the community changed** [read: hub README section names only; the details belong to an-hub]:
- The first hub extension used a "(23,23) bit interchange" and an "all-disjoint rational-center complex network on
  (24,24)". That is the ground-set knob, turned to just where my scan puts the optimum (h = 25 to 26).
- Later steps add:
  - three-stage Cayley covers;
  - partial source gauges;
  - paired-cube producers with coordinate-star centres;
  - shared cores reusing one auxiliary bank across three orthogonal blocks;
  - recycled registers and eliminated terminal accumulators, which cut W;
  - coordinated frames;
  - completed entrance banks;
  - "hyperplane-ended" and "dimension-4" helper endings (terminal slack on helpers);
  - pricing every positive-rank climb and descent as one recursive child (`native-complex-moment552` PROOF.md).
- On the Fourier side, Shea's Lean-checked #588 unit has W = 12,208 roles and 5,040 banks. That is 17 orders of
  magnitude fewer arrays than OpenAI's.

## 5. Knobs, and what stops each

| knob | effect in OpenAI's design | what stops it |
|---|---|---|
| fewer break-even arrays per saved dimension | the dominant factor: side carriers are 99.9995 % of W; without them h = 25 would give 9.9e-7 | the carriers are one per nonzero off-diagonal entry of the bridge matrix D, and a matrix with ones on the diagonal has rank at least the independence number of its support. So (bridge rank) × (carriers per pair) ≳ (pairs per invocation), and the bridge loss is about rank × label span. Fewer carriers means a bigger bridge loss [argued] |
| terminal slack on more arrays (carriers, centres) | each array with a source label of dim r saves 2r if its data crosses to another array | an array whose data stays home gains nothing (2.2). Carriers would have to be crossed in pairs too, and each crossing needs its own bridge |
| higher-rank source labels (dim U_a = r > 1) | slack 2r per pair | gates touching X must contain all of U_a and gates touching Y must avoid it. The bridge then drops the span of all the U's it carries, so the loss grows with r. Carriers need orthogonal *subspaces*, so fewer pairs qualify |
| smaller block m | the slack per array is r/m | the design needs v > 3h(h+1), so h ≥ 22 and m = h³ ≥ 10,648. The three stages (a signed swap is three shears) force a three-factor label space |
| pricing whole edges as children | ×80 to ×360 on the same network | none. It is a free gain once the recursion handles children of unequal size (hub moment equation) |
| ground-set size h | ×1,600 (h = 25) | d ≈ h³/6 and m = h³, so the saving falls like h^{−6}. Below 22 the bridge loses more than the crossing gains |
| nesting the unit or powering it | nothing | the recursion already reuses the unit. Powers and Kronecker products never beat the best single unit (MAP §1.1, Φ barrier) |
| slack on both banks (U_a on X, V_a ⊥ U_a on Y) | at most doubles the gain | needs a bridge in both directions every stage, so the loss doubles too [guess] |

## 6. What an-synth may want from this track

- **The one-line mechanism.** *Saving = 2 × (frame dimension carried across arrays by crossing routes) − 2 × (frame
  dimension bridges must drop)*. Everything else (carriers, restoration) is break-even and only dilutes.
- **A necessary condition for any unit of this type at small block** [argued]:
  - Most arrays must carry source labels and have their data end on another array.
  - The bridges must be shared by many pairs, through a matrix that is low rank over C and vanishes on the
    non-orthogonal pairs over F_2.
  - Junk cleaning must not need a full array per pair.

  At blocks 3 to 5 with W = 2 or 3 the bridge alone already costs more than any crossing can save (2.3). So an
  OpenAI-type saving cannot be what is missing from our block-3 to block-5 near misses. It needs hundreds to millions
  of arrays.
- **Seeds for scaled constructions** [guess; for an-synth to rank]:
  1. carriers that are themselves crossed pairs with slack;
  2. bridges whose junk is cleaned by another low-rank bridge instead of carriers (a second Alon pair);
  3. rank-r source labels on subspace designs with orthogonal-subspace neighbours.

  For each, the exact identity Δ = 2Σ dim σ − 2·(descents) gives the predicted saving from the counts alone. The test
  is the brute-force label check in `calc/labels_check.py`, rewritten for the new labels.

## 7. Checks run (all our own code; nothing third-party was run)

| file | what it checks | result |
|---|---|---|
| `calc/family.py` | the paper's integers: W, s, Δ, η_b, η_c, 1 − θ; scan of h | all reproduced exactly; best h = 25 (4.18e-10) |
| `calc/labels_check.py` (+ `labels_check_h8_h9.out`) | builds every label as an explicit subspace of F_2^{h³} for one invocation per stage, one representative per wire type and side-pair type. Checks nesting, residual dimensions against the table, nondegeneracy, a norm-one vector in every residual, stage boundaries and terminal labels, then sums over the whole network | h = 7, 8, 9: no problems; the only descent is centre 3→4; the sum equals Wm − 2N + 2L exactly |
| `calc/frames_check.py` | C_z² = R_z, C_z⁻¹ = R_zC_z, Φ(whole) = C^{⊗m}, Φ additive on orthogonal lines, one edge = one directional kernel, exceptional ratio R_u·C^{⊗m} | all hold at m = 2 to 5 |
| `calc/moment_price.py` | the same network priced with whole edges as children | ×360 at h = 100, 3.26e-8 at h = 26 |
| `calc/caps.py` | gain, loss and dilution per h; no-carrier and 1/m ceilings | table in 3.2 |
| one-liner in RUNLOG | SCS = e^{iπ/4}H, C² = X | holds |

## 8. Predictions written before reading (STATUS.txt), and outcomes

| prediction | outcome |
|---|---|
| P1: a small unit, tensor-powered | **wrong.** It is one huge finite unit (block 10^6, about 10^21 arrays), used recursively |
| P2: shared sums, no division | **half wrong.** The saving is crossing slack, not shared sums, and the complex network halves |
| P3: cap from design and size, about 1e-13, no route to 1e-1 | held (2.1e-13); the cap is structural |
| P4: bigger block is the knob, stopped by search cost | **wrong.** The knob runs the other way (smaller h is better, down to 22); there was no search |
| P5: lineage "transfer theorem first, explicit unit second" | **wrong in order.** The explicit network came first (integer multiplication, 23 Sept); the tensor-savings paper is a separate existence proof |

## Sources

- OpenAI, *An explicit power saving for the exact discrete Fourier transform*, 25 Sept 2026:
  `upstream/openai-math/preprints/An-explicit-power-saving-for-the-exact-discrete-Fourier-transform-September-25-2026/build/sections/{introduction,network}.tex`.
  Read in full: the introduction and the network section (Lemmas 2.1–2.3, Proposition 2.4, Lemma 2.5, Theorem 2.6,
  Remark 2.7).
- OpenAI, *Integer multiplication below n log n*, 23 Sept 2026:
  `.../Integer-multiplication-below-n-log-n-September-23-2026/build/sections/03-motifs.tex`, read in full.
- OpenAI, *Finite tensor savings and exact Fourier circuits*, 25 Sept 2026:
  `.../Finite-tensor-savings-and-exact-Fourier-circuits-September-25-2026/build/sections/`. Read in full: 01, 02
  (Definition 2.1, Lemma 2.3) and the dense-pair lemma in 07. Read in part: the statements of 04, 05, 06 and 08.
- OpenAI repository README and CONTENTS.md, families 109 and 130; Lean constants in
  `lean/OAI/Computability/FourierTransform/Goal.lean` (paperDelta = 6871402692000000).
- Hub README (`checks/wht29/sources/hub-main/README.md`): section names and κ values only.
- `checks/wht43/wave16/d-I3/src/research_native-complex-moment552_PROOF.md`: hub pricing (every
  positive-rank climb and descent is one recursive child).
- Shea's #588 RUNLOG and RESULT:
  `checks/wht43/bound-check/d-lean-route/shea-src/research/complex-transfer-588-596/` (W = 12,208;
  5,040 banks).
- Our failure modes: `moonshot/NOTE-wave2.md` (absorb4), `NOTE-wave3.md`, `meta2/NOTE.md` (pattern identity),
  `moonshot/o3/NOTE.md` (clean carriers); `meta/MAP.md` §1.3; `exist2/NOTE.md`.
