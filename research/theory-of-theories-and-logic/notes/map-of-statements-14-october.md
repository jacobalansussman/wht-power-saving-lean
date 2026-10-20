> **Working note, kept as written.** An AI agent wrote this during the theory-of-theories round of 14 October 2026 (the round's map of every statement, game and method, folded together with the second abstract wave the same evening). Private paths and machine names
> were removed, and so were the parts listed at the end of this header; nothing else was changed. No person has checked it. Where it disagrees with this section's
> [README](../README.md), the README wins, because it follows the later skeptic verdicts.
>
> What changed after this map was written (the README gives the details):
> - Theorem C, Lemmas D and E (one exponent) were upheld by one skeptic, who corrected the screen: a family must keep
>   a growth exponent of at least 0.9 at arbitrarily large blocks, not at every block.
> - The mixed units 1 + 4 in 10, 2 + 8 in 20, 4 + 21 in 49 and 7 + 38 in 88 at block 3 are excluded in the engine game
>   (Theorem M, two skeptics), and a construction round's Theorem T (upheld) excludes them in the widened game too.
> - Every block-3 lift proves at most 2.7986 moves per array (below the hook), and 18 <= b_4 <= 36 is all that is
>   shown at block 4; the statements about b_3 = 2 and b_4 = 38 or 35 here are superseded.
> - In the chi game the ceiling of every certificate is 169/57 (two skeptics, given an earlier reduction with no skeptic);
>   the interval-lift minimum there is between 9 (exact; a solver bound gives 12) and 18.
> - Block 2 never saves in any game whose free steps touch a bounded number of entries (Theorem 2W, one skeptic, 15
>   October), so 2 in 3 and 3 in 5 in the everything-free game are closed; and W + m paid moves are needed for W >= 2
>   arrays in every such game (Theorem T', upheld), which closes 2 in 6 at block 5.
> - "Making a relation-using approximate unit exact costs a factor 2 per level" (Theorem X's companion claim) is refuted.
> - The "all 652 symmetric block-4 skeletons are impossible" claim of the second wave is refuted as stated; 16 of the 17
>   skeletons that pass the relaxation with every coefficient free have no unit, and one is open.
> - The second wave's closing test for the near miss ("average sin^2 below 0.7 at depth 3") is necessary, not sufficient:
>   the correct test is that the depth-3 average is below 0.7 times the edge average.
> - Sylvester-Gallai counting is capped by some constant, but not at 1/3 per edge (refuted).
>
> Removed from this copy: the four task briefs of Appendix A (agent instructions with machine rules) and two
> bookkeeping items of section 7 about clock times and machine load.

# MAP: everything we know about the exponent saving, as a lattice of statements

meta-map, 14 October 2026, written about 08:45 to 09:20 EDT by the Mac clock (several notes carry later clock times;
see §7). Input: the brief (`BRIEF.md`), cycle 1 of this round (`spectrum/`, `border/`, `cohomology/`, `reps/`,
`literature/`), the existence round (`ex-dimension/`, `ex-laser/`, `ex-legit/`, `ex-v-legit-compute/`), and the notes
it builds on (`../proofs/NOTE.md`, `../wave16/b56-note/NOTE.md`, `../wave16/tb4-note/NOTE.md`, `../wave17/NOTE.md`,
`../wave17/xc3-note/NOTE.md`, `../wave17/ut-note/NOTE.md`, `wht35/wider/w-bounds.md`).
No skeptic has read this map. Nothing was posted. Labels as everywhere: [proved] (a PROOF.md exists), [computed],
[argued], [guess], [read]; **upheld** = two skeptics did not refute it; **unchecked** = no skeptic yet;
**derived here** = arithmetic or a short argument done for this map (checked by me only).

Saving of a unit "W in M" at block m: s = ln(mW/M)/ln m. Standing bar: s = 1e-1. Values quoted as 1.07e-1 or 0.107.

---

## 0. The answer to the user's question, in ten lines

1. **Yes, there is a theory of theories, and it already predicts.** Every lower bound we own is an additive monotone of
   one per-block monoid (Theorem R). Five invariants decide which monotones can exist and how far they go: the
   *forced kernel* of the game, the *frame geometry*, the *multiplicities* of the free group, the *characteristic*
   of the rank inequalities a ceiling proof may use, and the *θ-shadow* (the only part of a certificate the exponent
   sees). Section 1 is the stack; section 4 lists 30 predictions it made and what the calculations said.
2. **The games form a lattice, and certificates die exactly where the free group grows.** Engine (E) and unimodular
   chain maps (XQ): any Clifford representation. Running XOR (Xchi): Pauli-blind only, on a Hermitian geometry over
   F_4. Everything free (XW: chain map + Toffoli or controlled-S): no product-form certificate at all (reps Theorem 2).
   Units move the other way: a unit in a small game is a unit in every bigger one.
3. **The exponent is one number** (spectrum Theorem C), so certificates and searches answer the same question. A
   certificate family can settle the 1e-1 bar only if its θ-shadow stays ≥ 0.9 at every block. Hooks, Venn sets,
   coordinate sets and projections provably cannot. Lifts are the only candidate left (§3, §5).
4. **The construction side has its own theory now.** A saving must be one finite unit (no powering, ex-laser Theorem
   3), exact or formally polynomial (meta-border Theorem X), in special position (ex-dimension), dividing by 2,
   non-conformal, using address maps or phases, with multi-term contents (§1.3). Any complex coefficients count.
5. **Open and above the bar**, by game: block-3 mixed units 1 + 4 in 10 (0.119) and 7 + 38 in 88 (0.125), exactly
   tight on the hook (E, Wd); block 4: 3 in 10 (0.132), 4 in 13 (0.150), grouped 2,2,2,1 (0.1675, tight) (E, Wd);
   block 5: 2 in 6 (0.317, E; the bound 2.99172 misses 3 by 0.008), 3 in 7 (0.474, Wd); chi game block 3: 3 in 8
   (0.107), 2 in 5 (0.166), 3 in 7 (0.229); everything-free block 2: 2 in 3 (0.415), 3 in 5 (0.263).
6. **Four decisive pairs** (each an abstraction with the calculation that tests it), ranked by how much they change
   what we do: (1) **lifts**: local-to-global bound on lift excess plus a gluing search, deciding the only
   θ-consistent candidate and 5-in-14 / 2-in-7 / 3-in-10 / block-5 W-in-3W for lifts; (2) **hermitian**: the Fano
   LP on the subspace lattice of F_4^3, deciding whether any certificate can ever exclude 3 in 8 in the chi game;
   (3) **xw**: the first theory-constrained exact search with controlled-S on the menu, in the only game no proof can
   enter (2 in 3, 0.415); (4) **finite-group**: every certificate at block 3 as an Sp(6,2) computation, deciding
   whether the hook is optimal for 5-in-14 and for the above-bar mixed units, with a first screen of Sp(8,2).
   (The four task briefs of the original Appendix A are agent instructions and are not reproduced here.)

---

## 1. The stack: calculations, methods, theories of methods, and the meta level

### 1.1 Layer 3 (meta): one exponent, two sides

- **Theorem C** [proved, unchecked; `spectrum/PROOF.md` §2]: every monotone homomorphism of the composition semiring
  ([m]·[k] = [mk]) is n ↦ n^θ or ≤ 1 everywhere. One exponent for all block sizes; the saving is the supremum over
  single units.
- **Lemma D (θ-shadow)** [proved, unchecked]: a certificate with profile ψ_1..ψ_m forces z ≥ θ(ψ) =
  min_j ln(ψ_m/ψ_j)/ln(m/j) for every unit at that block. **Lemma E**: mixtures of certificates are never better than
  the best one (also pf-duality Lemma 1.2, ladder "sums and mixtures give exactly Theorem R" [upheld]).
- **Proposition G** [proved, unchecked]: for words without division by 2, n ↦ n (θ = 1) is a spectral point at every
  block (the mod-2 frame distance). It is the only θ-consistent point known. Lifts are its characteristic-0
  deformations.
- **B1** [proved, upheld; `../proofs/ladder/PROOF.md` §4]: every rank method also bounds approximate (border) units.
  **Theorem X** [proved, unchecked; `border/EXACTIFY.md`]: only formally polynomial border units are worth anything
  (O(log log N) overhead); relation-using ones lose a factor ≥ 2 per level under every known repair.
- **Φ barrier** [proved in the job model, unchecked; `ex-laser/PROOF.md` §3]: no combination, power, Kronecker product
  or restriction of rules beats the best single rule; Kronecker powers with any hashing tend to exponent 1.
  **τ-transfer** (Theorem 1 there): a mixed unit with rank-1 moves is worth its τ_c, Σ_k W_k k^τ_c = M.

So the two sides meet in one number. The certificate side can only bound it through θ-consistent families; the
construction side can only reach it through one finite unit.

### 1.2 Layer 2 on the certificate side: five invariants and the barriers between method levels

| Invariant | What it decides | Results | Status |
|---|---|---|---|
| Forced kernel | which representations a game allows | E: none; XQ: none at the Pauli level (cohomology class vanishes, Theorem E+); Xchi: every column-uniform Pauli and i^f (ut-direct Prop C) → Pauli-blind; XW: the whole one-column Clifford group (reps Theorem 2) → nothing of product form | Prop C upheld; E+ and Theorem 2 unchecked |
| Frame geometry | certificates as subspace families on frames | E/XQ: dual polar graph DSp(2m,2) (15, 135, 2295, 75735 frames, m = 2..5); Pauli-blindness costs nothing in E (reps Theorem 1); Xchi at alpha: generators of H(2m−1,4) (27, 891) | Averaging Lemma upheld; Theorem 1 unchecked |
| Multiplicities | cut vs label certificates | P_0-multiplicity-free module → every certificate is a cut, finitely many, spectrally capped (pf-duality Lemma 1.5, Theorem S); only repeated pieces carry labels that beat cuts; two-factor modules Λ^a ⊗ Λ^b are multiplicity-free for the widened free group (reps Lemma 4.1), so labels need ≥ 3 tensor factors there | proved; Lemma 4.1 unchecked |
| Characteristic | what a ceiling proof may use | the mod-2 family is perfect, so every proof of a ceiling below m must use a row false over GF(2) (Corollary C, B5); the DFZ inequality (arXiv:1311.4601, Thm 8.6) is the one used; every saving it allows is bought by Fano gadgets (Lemma G) | B5 upheld; D3/D4 unchecked |
| θ-shadow | what the exponent sees | hooks, Venn, coordinate, projections: θ → 0 [proved via ceilings]; lifts: 0.856, 0.711, ≈0.568, ≈0.429 at m = 3..6 (start lifts) [computed, m ≥ 5 sampled] | unchecked |

Ceilings of method levels (most moves per array a whole level can ever certify; plain is m):

| Level | m = 3 | m = 4 | m = 5 | m = 6 | Source |
|---|---|---|---|---|---|
| colourings / coordinate frame families (engine) | 1.75 | 1.875 | 1.9375 | 1.969 | reps Prop 3 [proved, unchecked] |
| same on Hermitian frames (chi game) | 21/16 = 1.3125 | 1.328 | | | reps Prop 3 |
| hooks V ⊗ Λ^k | 2.8 (reached) | 10/3 (k = 8 stops at 41/13 w.h.p.) | 62/17 = 3.647 | 42/11 = 3.818 | b56 [proved] |
| coordinate sets in tensor products of exterior powers (widened) | 3.5 | 3.75 | 3.875 | 3.9375 | b56 [proved m = 3..10] |
| lifts (engine) | 2.917 (Fano count); 105/37 = 2.838 if b_3 = 2 (§5, derived here) | 3.918 (Prop F); ≈3.453 if the interval optimum 38 is the true minimum | — | — | w-bounds Thm F; pf-duality Prop F, §4 |
| Frobenius / spectral duals | | cannot certify below 3.75 | | | pf-duality Theorem S remark |
| all Theorem-R certificates | ≤ 169/57 = 2.965 | ≤ 3.9426 | ≤ m − 0.057 | | pf-duality D3/D4 [unchecked] |
| filtered rung RF / depth | adds nothing at the tight cases of blocks 3, 4 | | closes 1-array units, 2 in 5 (Wd) | closes 2 in 6, widened 2 in 5, 3 in 7 | ladder [upheld] |
| induction from small blocks | | | ≥ 2.99172 (E) | ≥ 2.97290 (E); ≥ 2.9447 for every m ≥ 4 | pf-induction [upheld] |
| product-form certificates in XW | 0 | 0 | 0 | 0 | reps Theorem 2 [unchecked] |
| 2-adic lattice (2-integral mixing, right at f = 2) | 1.875 | 2.1875 | → √m | | ut-direct §6 [upheld]; the only bound in XW |

### 1.3 Layer 2 on the construction side: what any saving unit must be

| Necessary condition | Status | Source |
|---|---|---|
| one finite (possibly mixed) unit; no powering or mixing of non-saving rules | [proved, job model, unchecked] | ex-laser Theorem 3, Prop 4 |
| exact, or formally polynomial if approximate | [proved for known repairs] | border Theorem X |
| any complex coefficients are acceptable for the challenge | upheld refutation of ex-legit's field caveat | `ex-v-legit-compute`, HANDOFF 10:30 |
| sits in special position; dimension counts cannot prove it exists | [proved] | ex-dimension Lemmas B, C, D |
| divides by 2 | [proved at block 3 with Clifford free steps; argued for all m]; not established with Toffoli free | w-bounds Thm F; proofs NOTE |
| uses address relabellings or phases | [certificate-checked] (commutative slice: none at block 3, ≥ 15W/4 at block 4, ≥ 31W/8 at block 5) | pf-logging |
| non-conformal mixing step (not a scalar times a unitary on all arrays) | [argued] | literature T3 (Morgenstern 1973 from memory; Ailon arXiv:1305.4745 abstract) |
| total rank ≥ W + m − 1 | [proved, upheld] in E and Wd; unknown in Xchi/XW (§7 item 9) | ladder Theorem D |
| multi-term contents (single-term contents: 2 arrays need 2m moves) | [certificate-checked] | b4-depth4 family S |
| tight targets are rigid: 5 in 14 all moves rank 1 with +80 each, no move fed by a single free image; grouped 2,2,2,1 no waste, multi-term filters on every move, at least one shared move | t514 and tb4-tight [upheld] | |
| Fano gadgets with a halving step over a shared base (block 3: at the start; block 4: start, after one shared move, mirrored) | [argued] | pf-duality §7 (the dual of the LP) |
| catalysts never beat Theorem R | [proved, 3 lines, unchecked] | literature T4 |

### 1.4 The game lattice

```
            XW  (everything free: chain maps + Toffoli/CS)
           /  \
   Wd (Toffoli,   Xchi (running XOR chi; P free)
   all phases)      |
          |       XQ (unimodular chain maps row_a ^= row_b << 1, CZ chains, boundary Cliffords)
           \      /
            E (engine: Paulis, affine maps, Clifford Z4 phases)
```
Lower bounds flow downward (proved in a bigger game ⇒ holds in every smaller one); units flow upward. Variants that
change the target, not the game: spectator rows (rank-k targets), mixed targets (τ-transfer; needs a new recursion in
our engine), period-2 chain maps (= block 2m with paired directions), OpenAI's finite-win game (one array, 2^b entries).

---

## 2. The predictors (abstract results), with identifiers used below

| Id | Statement | Status | Where |
|---|---|---|---|
| A1 | Theorem R: every certificate gives M ≥ W·dem/cap | Lean-checked for the engine | `../proofs/lean/` |
| A2 | Averaging Lemma; reps Theorem 1 (engine certificates = Pauli-blind = uniform frame families on DSp(2m,2)); pf-duality §9 (the exact finite-group form) | Averaging upheld; Thm 1 unchecked; §9 argued | w-bounds §3; `reps/PROOF.md` §1 |
| A3 | Lemma 1.5 / multiplicity principle; Theorem S spectral ceiling for cuts; Prop 3 | proved; unchecked | pf-duality §1–2; reps §3 |
| A4 | Corollary C, B5: ceilings below m need char-0 rows; Theorems D3/D4 (169/57, 3.9426); Lemma G (Fano gadgets) | B5 upheld; D3/D4 unchecked | pf-duality §6–7; ladder |
| A5 | Theorem C, Lemma D, Lemma E, Prop G (spectrum) | unchecked | `spectrum/PROOF.md` |
| A6 | B1; Theorem X; relation-using units lose a factor 2 | B1 upheld; X unchecked | ladder §4; `border/EXACTIFY.md` |
| A7 | Theorem D (W + m − 1); RF adds nothing at tight blocks-3/4 cases | upheld | ladder §3 |
| A8 | Family ceilings (hooks, coordinate, code products, two-factor Venn ≤ hook) | proved; upheld | b56-note |
| A9 | Induction: units embed; ≥ 2.9447 for m ≥ 4 (engine); projections multiply by ≤ 1.1; dichotomy λ* | upheld | pf-induction |
| A10 | Cohomology: the H^2 class vanishes; Theorem E+ (Q_loc letters: 2W, 14W/5, 41W/13); obstruction words for non-banded letters; ring mechanism fails for Toffoli | unchecked (Toffoli part computed in floats) | `cohomology/PROOF.md` |
| A11 | Prop C (chi forces Pauli-blindness); Theorem R-infinity; Theorem B (column-uniform letters incl. non-Clifford + end columns) | upheld | ut-direct, ut-normal |
| A12 | Readings at alpha / 1: 27-line perfect certificate (block 2), 14W/5 for letters with a reading at 1 (block 3); no 27-line-type certificate at F_16 or F_8 (Hirzebruch; Klein) | upheld | xc3-s1s2, xc-round2 |
| A13 | reps Theorem 2 (XW: product-form certificates prove nothing) | unchecked; group facts computed at m = 2, 3 | `reps/PROOF.md` §2 |
| A14 | 2-adic lattice potential; halving necessity (Theorem F); Morgenstern–Ailon non-conformality | upheld; proved (Clifford); argued | ut-direct §6; w-bounds §4; literature T3 |
| A15 | τ-transfer, Φ barrier, truncation Lemma 6 | proved in the job model, unchecked | `ex-laser/PROOF.md` |
| A16 | Counting cannot prove existence; orbit test weaker than Theorem R; XW collapses at every fixed f (Burnside) | proved, unchecked | `ex-dimension/PROOF.md` |
| A17 | Theorem V: a shadow scheme exists for every order of moves with slack ≥ 0 in the f = 9 block-4 hook | upheld | tb4-tight |
| A18 | literature T1: Schönhage's τ-theorem → mixed room; 89/47 rank-2 certificate (excludes 7 in 13, 9 in 17) | computed, two codes, unchecked | `literature/NOTE.md` |
| A19 | Gossip bound (finite-win game H on 4, 8 entries) | proved, unchecked | wave 17 §4.3 |

---

## 3. The lattice: open case × game × method

Columns: **status** (what is proved), **predicted by** (abstractions that say something about the case), **tested by**
(calculations, with outcome), **gap** (what abstraction says that nobody computed, or what calculation would decide
an abstract question). "Owner" says which side can still settle the case: P = some proof method has a ceiling above
the target; C = no known method can, so only a construction (or a new kind of method) can.

### 3.1 Block 2 (plain 2 per array)

| Case (s) | Game | Status | Predicted by | Tested by | Gap / owner |
|---|---|---|---|---|---|
| everything (2 in 3 = 0.415, 3 in 5 = 0.263) | E | closed, 2W | A1 | Lean | — |
| same | Wd, XQ | closed, 2W | A8 (hooks), A10, A12 (reading at 1) | t1–t4 (cohomology): all pass | — |
| same | Xchi, letters with a reading at alpha (or all at 1) | closed, 2W [upheld] | A12 (27 lines, Fermat cubic) | lines27.py exact; ut-hunt / xc-search 0 hits (consistent) | — |
| 2 in 3 with mixing maps sharing no reading ({chi, R3}: 85 frames of Sp(4,4); F_16: 325; F_8: 4745) | Xchi residue | open | A12: no perfect certificate of 27-line kind at 325 / 4745 (Hirzebruch, Klein); ratio in (1.5, 2) not excluded | ut-hunt bd3 LS: 0 hits (residual 0.27, evidence only) | P/C: a non-perfect certificate with ratio > 1.5 on 85 or 325 frames would close it; none tried |
| 5 in 9 with period-2 chain maps (0.152) | = block-4 engine with paired directions | open; best 1.7826 < 1.8 | A2 at block 4 | xc-algebra numbers | P: any block-4 certificate with rank-2 ratio > 1.8 |
| 2 in 3 (0.415), 3 in 5 (0.263) | **XW** (chi or E with P, plus CS) | **open; nothing excludes them** | A13 (no product-form certificate), A16 (fixed-f collapse), A14 (lattice: 2 in 3 exactly at 1.5, so it must halve or be lattice-tight), A10 §5 (the obstruction words), T3 (non-conformal) | **no search has ever had CS on its menu** (xc-search, ut-hunt, xc3-search3: Clifford letters only) | **C. Task 3 (xw).** |

### 3.2 Block 3 (plain 3)

| Case (s) | Game | Status | Predicted by | Tested by | Gap / owner |
|---|---|---|---|---|---|
| ≤ 4 arrays (3 in 8, 2 in 5, 3 in 7, ...) | E, Wd | closed by 14W/5 | A1, A8 | Lean (except dem 224, certificate-checked) | — |
| same | XQ | closed by 14W/5 (three routes) | A12 at 1 + hook frames [upheld]; reps Cor 1.2; cohomology E+ | hook_shadow / engine_shadow: 80/144/224 over all 135 frames (reps); t1–t4 | — |
| 5 in 14 (0.063), 6 in 17 (0.052) | E, Wd, XQ | open, exactly on / just above the hook; BEST_3 ∈ [2.8, 169/57] | A4 (LP ceiling), A6 (B1: if an approximate unit exists, no rank method excludes it; Theorem X: it would be relation-using and worthless), A3 (only labels can beat cuts) | t514-cert: no certificate > 14/5 in any family tried; hookmult: engine labels in the hook module inert; V^⊗3 labels ≤ 7/3; lifts best E = 162 (2.561); pair-anchored lifts have E ≥ 150 [w-bounds §5.2] | P/C: lifts need E ≤ 67 (5 in 14) or ≤ 55 (6 in 17), and E ≥ 54 if b_3 = 2 (§5); the finite-group form is uncomputed. **Tasks 1 and 4.** |
| **mixed 1 + 4 in 10 (0.119), 2 + 8 in 20 (0.119), 7 + 38 in 88 (0.125)** | E, Wd (needs ex-laser's new recursion in our engine) | open; **exactly tight on the hook** (2.8 W1 + 1.8 W2 = M; derived here) | A15 (worth τ_c), A18 (mixed room 0.1255), A6, rigidity of t514 carries over | literature: hook, Λ²⊗Λ³ J={1}, Λ³⊗Λ⁴ J={0,2} all leave them open; ex-dimension: blind relaxation solvable | P/C: a certificate with ψ3 + 4ψ2 > 10ψ1 excludes 1 + 4 in 10; nobody has optimised one. **Task 4** (step 0: the LP for this objective). |
| rank-2 targets, 1 spare row: 7 in 13 (0.107), 9 in 17 | Wd ⇒ E | **excluded by 89/47** [computed, two codes; unchecked] | A18 | indep_venn_check.py on Colab | needs a skeptic (exact rank of the 1960 × 1960 block; keptness) |
| rank-2, 1 spare row: 10 in 19 (0.074), 11 in 21 (0.067) | Wd | open, below the bar | | | — |
| 3 in 8 (0.107), 2 in 5 (0.166), 3 in 7 (0.229) | **Xchi** (chi with Clifford letters) | open; best 231/131 = 1.763 (Fermat + Teichmüller lifts over Q(ω), exact) | A11 (Pauli-blind forced), A2 at alpha (891 frames), A3 (colourings ≤ 21/16), Fano count (lifts ≤ 2.917), Fermat-anchored lifts ≤ 2.287, Weil hooks ≤ 2.18 | xc3-cert3: every natural family stops below 7/3; anneal891 found no improvement; xc3-search3: 0 hits (low power) | **The engine analogue of an all-certificate ceiling (169/57) has never been computed at alpha. One LP decides P or C. Task 2.** |
| same | XW (chi + Toffoli/CS) | open; only the lattice bound 1.875 | A13, A10 §5.2 (ring mechanism fails for Toffoli) | none | C. Task 3, stage E. |

### 3.3 Block 4 (plain 4)

| Case (s) | Game | Status | Predicted by | Tested by | Gap / owner |
|---|---|---|---|---|---|
| 2 in 6, 3 in 9, 4 in 12; 3 in 10 with ≥ 3 rank-2 blocks | E, Wd, XQ | excluded by 41W/13 [certificate-checked; XQ via E+ / Cor 1.2] | A8 | hookG, b4-depth4, b4-cert floats agree | — |
| 2 in 7 (0.096, below the bar) | E, Wd | open | hooks ≤ 10/3 cannot; Venn two-factor cannot; spectral duals cannot decide (≥ 3.75); block-3 knowledge cannot (lp_cond: c + 1); commutative slice and family S exclude it | tb4-beam 0 hits; b4-cert: no family ≥ 3.5 | P (coordinate ≤ 3.75, lifts, labels) / C. Lifts: interval optimum 38 found, and b_4 ≥ 35 would prove lifts never close it (§5). **Task 1.** |
| 3 in 10 (0.132), 4 in 13 (0.150), 5 in 16 (0.161), 6 in 19 (0.169) | E, Wd | open | hooks cannot (ceiling 10/3, k = 8 stops at 41/13); coordinate sets could (3.75); lifts could only if the whole-graph excess per edge < 0.2 / 0.231 (local optimum ≤ 0.158 allows it) | tb4-beam 0 hits (4 in 13: 77% of final tests undecided) | P via lifts (**Task 1**) or labels in Sp(8,2) modules (**Task 4 stage 6**); C otherwise |
| grouped 2,2,2,1 (0.1675), exactly tight; 2,2,1,1,1 (0.134); 3 arrays 2,2 + six (0.163) | E, Wd | open | A17 (f = 9 hook cannot exclude any slack ≥ 0 schedule), A6, tightness rigidity | tb4-tight: ~480 + 76 single terms, 1,150 two-term filters, 11 averaging filters all fail the zero-waste test | cheap calculation pending: make "tight" exact (exact demand; exact rank-2 cap ≤ 48576). Not in this round's four (no abstraction to test), but worth one agent-hour any time |
| spectator rank-2 at 4 rows: 5 in 9 (0.152); rank-3: 2 in 5 (0.166) | Wd | open; 1.7692 < 1.8 and 2.4883 < 2.5 (missed by 0.012) | A3 (needs ≥ 3 factors in Wd), A18 (no Schönhage gap at block 4 with known certificates) | spect venn runs | P: a three-factor sweep at block 4 (literature next 2) |
| everything | Xchi, XW | not studied | | | — |

### 3.4 Blocks 5, 6 and beyond

| Case (s) | Game | Status | Predicted by | Tested by | Gap / owner |
|---|---|---|---|---|---|
| 1 in 3, 1 in 4; 2 in 5 | E, Wd | closed by Theorem D (2 in 5 also in Wd) | A7 | out_ladder_check | — |
| 2 in 6 (0.317), 3 in 9 (0.317), 2 in 7 (0.222), ... | E (block 5) | open; 2.99172 per array [upheld] misses 3 by 0.008; engine-priced cap ≈ 0.43 (lift; d ≥ 2 sampled) | A9, A5 (lifts price grouped blocks fairly), A3 | liftprof (m = 5 exact on edges, sampled d ≥ 2); exact chunks prepared, not run | P: a lift with E < 1,565,190 (now 2,056,380) closes every W in 3W (derived here); a filtered certificate with dem_1/f > 2.5 also would (ladder). **Task 1 stage D.** |
| 2 in 6 (0.317), 3 in 7 (0.474), 4 in 9 (0.496) | Wd (block 5) | open; 31/14 [computed], 62/29 [certificate-checked] | A3 (≥ 3 factors needed for labels), coordinate ceiling 3.875 | b56 two-factor sweeps complete | P in principle (three-factor coordinate sets); nobody has a candidate |
| 3 in 11 (0.193; beyond every hook), 2 in 8 (0.139; beyond every coordinate certificate) | E, Wd (block 5) | open | A8 (only label certificates could) | | C for every family we know |
| 2 in 7 (0.301), 3 in 9 (0.387) | E (block 6) | open; 2.97290 [upheld]; 2 in 6 closed by Theorem D | A9, A7 | | P/C |
| 2 in 7, 3 in 8 (0.453), 4 in 9 (0.547) | Wd (block 6) | open; 63/31 | | | |
| the bar itself (s ≥ 0.1 anywhere) | all | open | A5: a certificate family settles it only if θ ≥ 0.9 at every block; none known; lifts the only candidate (θ 0.856, 0.711, 0.568, 0.429 for start lifts) | spectrum's lift profiles | **Task 1 decides whether lifts can be θ-consistent.** |

### 3.5 Other models

- **OpenAI's finite-win game**: H on 16 entries, 28 to 31 calls open (≤ 3.9e-2, below the bar; the gossip bound's
  room peaks near 1.09e-1 around 12 bits) [unchecked]; zeta kernel on 8 entries closed at 12 [computed]. No
  abstraction from this round touches it.
- **Hub pricing (constants, not the exponent)**: ex-laser's truncation Lemma 6 gives a what-if of +8.7 % on hub
  #1042 (1.674989e-3 → 1.821393e-3) if idle children are per-array tails and the τ-transfer extends to grouped
  blocks; the recursion-lever round (`../bound-check/rl-*`) is checking this. Not part of this map's tasks.

---

## 4. Ledger: predictions from abstractions, and what the calculations said

| # | Abstraction → prediction | Calculation | Outcome |
|---|---|---|---|
| 1 | Theorem F / mod-2 family → no characteristic-0 perfect certificate; a saving at block 3 must halve | every search for certificates (ratio never 3) | confirmed; proved |
| 2 | B1 → nothing below the in-representation bound, border units included (meta-border P1) | bsearch, block 2, V⊗Λ²V | confirmed |
| 3 | meta-border P2: Theorem R is attained exactly in multiplicity-free representations | block 2 V⊗V: ceil(1.5W) exact for W = 1..4 | confirmed at block 2; block 3 unsettled (load) |
| 4 | meta-border P3: border-only schemes are relation-using | formal test on d-border's coset sector: deviation ~ δ/ε | confirmed on the only known family |
| 5 | Theorem X: formally polynomial units compose exactly; relation-using ones fail from depth 2 | exact_toy.py to depth 6 | confirmed |
| 6 | reps Theorem 1: the hook's frame shadow is uniform (one δ per distance) | hook_shadow, engine_shadow: 80/144/224 on all 135 frames | confirmed |
| 7 | the engine's 3-fold multiplicity in the hook module might beat 2.8 | hookmult, hookmult2 | refuted (labels inert) |
| 8 | Lemma 4.1: two-factor modules are multiplicity-free (widened) → Venn sets never beat cuts | explains b56's sweeps; V^⊗3 has labels | confirmed retroactively; labels found but ≤ 7/3 |
| 9 | Prop 3: colourings on 891 Hermitian frames ≤ 21/16 | xc3-cert3 fam22 wall at 1.3125 | confirmed exactly |
| 10 | Theorem 2 group facts | nogo_check m = 2, 3 | confirmed |
| 11 | cohomology: γ = 0, collapse multiplicative, ring restriction depends on L mod 4, trace condition | t1–t4 | all confirmed; ring identity fails for Toffoli |
| 12 | spectrum: Strassen axioms hold except lower boundedness; one direction = 2 kernel calls | archimedes.py (exact in Q(i)) | confirmed |
| 13 | spectrum: lifts price high-rank blocks fairly | liftprof m = 4, 5, liftsample m = 6 | confirmed in direction (block-5 engine-priced cap ≈ 0.43, not 0.83) |
| 14 | spectrum: block-6 lift θ 0.45–0.5, cap about 0.5 | liftsample m = 6 | slightly optimistic (0.429, 0.57) |
| 15 | lifts might be θ-consistent | start lifts m = 3..6 | refuted for start lifts (θ falls ~0.14 per block); open for optimised lifts (Task 1) |
| 16 | Lemma E: mixing families adds nothing | mixcap m = 4, 5, 6 | confirmed |
| 17 | literature T1: a Schönhage gap at block 3 gives mixed units more room | venn_sweep_m3 + independent code | confirmed (room 0.1255); side result 89/47 |
| 18 | τ-analog is what a schedule achieves | tau_schedule_lp, frontier_flow, lockstep_sim | confirmed numerically; then proved in the job model |
| 19 | literature T1 at block 4: no gap with known certificates | profile comparison | computed (hook k = 8 dominates ranks 2, 3, 4) |
| 20 | literature T2: approximate units where Theorem R is tight | none at block 3 | untested; meta-border predicts they would be relation-using |
| 21 | Φ barrier: Kronecker powers tend to exponent 1 | kron_check | confirmed |
| 22 | dimension counts could prove existence | exact_counts, dimcount | refuted structurally (plain schedules have negative expected dimension); savings that exist (V⊗V 2 in 3) sit on a non-dominant image |
| 23 | pf-duality: DFZ rows give a ceiling below 3 | LP | confirmed (169/57) |
| 24 | extra frames off the interval add rows | dpg.py on 135 frames | refuted (same 169/57) |
| 25 | block-3 knowledge can decide 2 in 7 | lp_cond | refuted (gives c + 1) |
| 26 | b56: the hook reaches 4(n−1)/(n+2) at k = n/2 | k = 8 at block 4 | refuted w.h.p. (stops at 41/13) |
| 27 | ladder: c2 ≤ 1.8 c1 for every certificate | Λ²⊗Λ³ | refuted (1.848) |
| 28 | the f = 1 shortcut transfers bounds to chi; class Q closed under products | P identity; E2 counterexample | refuted; repaired as Q_loc (cohomology) |
| 29 | readings: perfect 27-line certificate at alpha | lines27, 729 pairs exact | confirmed |
| 30 | searches in games later proved closed should find nothing | xc-search, ut-hunt (Clifford chain menus, block 2) | consistent (0 hits; then xc3-s1s2 proved 2W) |

What the ledger says about the method: abstractions that came with a calculation at the smallest size were right
about the *direction* almost every time; the refutations (26, 27, 28) were all extrapolations made without a small
check. The two biggest open abstract claims (lifts' θ-consistency, all-certificate ceilings off the engine game) have
no calculation behind them yet. That is where the next round goes.

---

## 5. Where abstraction now tells us something nobody has computed

1. **The lift ceiling is local, and computable** (derived here from pf-duality Lemma I2 and edge-transitivity; short
   proof, unchecked). For any lift on the whole frame graph, every edge lies in the same number of intervals (opposite
   pairs and the frames between them), and a lift restricted to an interval is an interval lift. So
   E/#edges ≥ b_m/e_m, where b_m is the least total excess of a lift of one interval (e_3 = 35, e_4 = 240,
   e_5 = 2077 covering edges), and **every lift has value ≤ m/(1 + b_m/e_m)**. pf-duality's interval annealing gives
   b_3 ≤ 2 and b_4 ≤ 38 (A-parts: 18, "never 17"); Prop F gives b_3 ≥ 1, b_4 ≥ 5. Consequences nobody has drawn:
   - if b_3 = 2, lifts are capped at 105/37 = 2.8378 at block 3 (E ≥ 54 on the 945 edges); 5 in 14 needs E ≤ 67,
     6 in 17 needs E ≤ 55: possible only with near-perfect gluing (gap ≤ 1.24 or 1.02);
   - if b_4 ≥ 35, **no lift can ever close 2 in 7**; closing 3 in 10 needs b_4 < 48 (true locally: 38) and a
     whole-graph excess per edge below 0.2;
   - θ ≥ 0.9 at block 4 needs b_4 ≤ 35.7. If b_4 ≥ 36, lifts are not θ-consistent, and with spectrum's screen no
     certificate family we know can decide the bar.
   → **Task 1.**
2. **An all-certificate ceiling for the chi game is one LP away.** Theorem I's proof carries over to the alpha
   reading (free letters fill P = Q:GL(3,4); the symbol group is distance-transitive on the 891 frames, computed). The
   interval is the subspace lattice of F_4^3 (44 frames, 147 covering edges), with 360 Fano (Baer) subplanes in place of
   one. No one has run pf-duality's LP there. For lifts the same local bound applies: b_α > 18.4 would mean no lift ever
   excludes 3 in 8. → **Task 2.**
3. **XW is unobstructed, and the theory says what a unit there must look like** (both a chain map and CS doing
   essential work; halving or lattice-tight; non-conformal; invisible at f = 1; using a free square root of a global
   Pauli string). Nobody has searched with CS on the menu. → **Task 3.**
4. **Every engine certificate is a finite-group object** (reps Theorem 1 + pf-duality §9): a module of Sp(2m,2) and a
   subspace kept by P_0. Cuts in multiplicity-free modules are finitely many and spectrally capped (Theorem S); labels
   need repeated P_0-pieces. Nobody has run through the 30 irreducibles of Sp(6,2). The above-bar mixed units at block 3
   are exactly tight on the hook, so this computation decides whether any certificate can touch them. → **Task 4.**
5. **Literature T4(c) is moot as a ceiling question**: by reps Theorem 1 a family on signed stabiliser states can never
   exceed the best frame family (they share the engine's supremum). Only useful as a different parametrisation.
6. **Block 5**: ladder's remark that a filtered certificate with dem_1/f > 2.5 would close 2 in 6 (engine), and §3.4's
   E < 1,565,190 for lifts, both target the 0.008 gap. The exact block-5 lift profile is prepared in chunks
   (`spectrum/run_m5_exact_chunks.sh`) and not run.
7. **Theorem D in XW** (see §7 item 9): if it holds, 2 in 3 is exactly depth-tight there (3 = W + m − 1), a strong
   structural constraint for Task 3.

## 6. Where one calculation would decide an abstract question

| Question | Calculation | Size | Task |
|---|---|---|---|
| Can lifts be θ-consistent? Can they close 2 in 7, 3 in 10, 5 in 14, block-5 W in 3W? | lower bound on b_3, b_4 (MILP with char-0 rank rows) + whole-graph gluing search | 16 / 67-frame intervals; 135 / 2295 frames | 1 |
| Is the chi game at block 3 certificate territory or construction territory? | Fano/DFZ LP on the 44-frame F_4 interval; interval-lift minimum b_α | small | 2 |
| Does removing all product-form certificates (XW) open sharing? | matched LS landscapes with/without CS; exact templates at f = 2, 3, 4; automaton | block 2 | 3 |
| Is the hook optimal at block 3, for the kernel and for mixed targets? | Sp(6,2) irreducibles, P_0-cuts, labels; LP with mixed objective | 30 irreducibles, degree ≤ 512 | 4 |
| Is the tight grouped 2,2,2,1 at block 4 really tight? | exact demand and exact rank-2 cap | 6435 × 6435 blocks | not assigned (cheap; no abstraction to test) |
| Is 89/47 right? | exact rank, keptness | 1960 × 1960 | skeptic |
| Does P2 (no border phenomenon in multiplicity-free reps) hold at block 3? | border runs V⊗V 4 in 7 | ~1 h quiet machine | not assigned |

---

## 7. Inconsistencies and corrections found while mapping

1. `reps/TABLE.md`, block 5 row "2 in 5 (0.431), 3 in 7, 4 in 9: Wd OPEN": **2 in 5 at block 5 is closed in both games
   by Theorem D** (W + m − 1 = 6 > 5; ladder §3, upheld). 3 in 7 and 4 in 9 stay open in Wd.
2. `reps/TABLE.md`, block 6 rows: **2 in 6 at block 6 is closed in both games by Theorem D** (7 > 6), not open; and
   **2 in 5 and 3 in 7 at block 6 are closed in Wd by Theorem D**, not open.
3. `reps/TABLE.md`, block 3 spectator row "7 in 13: E closed (2W, wave 17)": wave 17's 2W theorem is for OpenAI's
   network game (address maps that return the active bits in some order), smaller than the engine game (the skeptic
   noted that maps mixing active bits escape). The case is closed anyway in both games by the 89/47 certificate if
   that is upheld.
4. `reps/TABLE.md` and proofs NOTE: "lifts ≤ 3.918 at block 4" is the Fano-count ceiling; the interval data suggest
   ≤ 3.453 for lifts if b_4 = 38 is optimal (§5 item 1; unproved).
5. proofs NOTE / reps TABLE "lifts can close 5 in 14 with E ≤ 67": true in principle, but pair-anchored lifts have
   E ≥ 150 (w-bounds §5.2), the free annealer stops at 183–192, and E ≥ 54 if b_3 = 2.
6. wave17 NOTE §4.2, ut-direct §4.4 and §7.4: "best Pauli-blind bound at block 3 is 2.561" is wrong: ≥ 14/5
   (reps Theorem 1; xc3-cert3, upheld). Already flagged by reps.
7. b56-note: "3 arrays with four rank-2 blocks at block 6 (0.738), not excluded": excluded at engine prices by lifts
   and projected families (spectrum). Already flagged by spectrum.
8. ex-legit's field requirement (constants must lie in Q(ζ_d)) was refuted by both skeptics: any complex unit counts
   (Nullstellensatz; restriction of scalars costs a constant factor). The ex-legit NOTE still says otherwise.
9. **Theorem D outside E and Wd.** Its proof uses the coset module Q[G/F] and a formal identity in Q[G]. In Xchi and XW,
   "right at every f" yields identities such as P = (1+i)/2 + (1−i)/2 Z_str that are not formal, and Q[G/F] does not
   satisfy them. So W + m − 1 is not known in those games (Task 3 checks it).

---

## 8. The next round: four tasks, ranked by how much they would change what we do

| Rank | Task (folder) | Abstraction | Calculation | Open cases (s) | What changes |
|---|---|---|---|---|---|
| 1 | **lifts** (`meta/lifts/`) | θ-screen: lifts are the only θ-consistent candidate; local-to-global bound; Fano/DFZ obstruction; gluing | prove lower bounds on b_3, b_4; search whole-graph lifts at m = 3, 4 (5) against exact thresholds | 5 in 14, 6 in 17 (block 3); 2 in 7, 3 in 10, 4 in 13 (block 4); 2 in 6 (block 5, 0.317); the bar | either the certificate programme has a live family (and new exclusions), or lifts are provably not θ-consistent and every certificate method we know is out for the bar: machines move to constructions |
| 2 | **hermitian** (`meta/hermitian/`) | Theorem I at alpha + Corollary C + DFZ on 360 Baer subplanes | LP on the 44-frame F_4 interval (+ lift minimum b_α, gadget decoding) | 3 in 8 (0.107), 2 in 5 (0.166), 3 in 7 (0.229), chi game block 3 | decides whether the chi game at block 3 is certificate or construction territory; one LP, reusing tested code |
| 3 | **xw** (`meta/xw/`) | Theorem 2 + Theorem B + readings + lattice potential + obstruction words + Morgenstern–Ailon | landscape test, then exact templates at f = 2, 3, 4 and the automaton, with negative controls that are theorems | 2 in 3 (0.415), 3 in 5 (0.263); then 3 in 8 with Toffoli | a hit is the biggest result on the board; a clean null with power shown, or a landscape no better than the Clifford one, tells the theory there is an obstruction outside product form |
| 4 | **finite-group** (`meta/finite-group/`) | reps Theorem 1 + pf-duality §9 + Lemma 1.5 + Theorem S | the mixed-objective LP; all irreducibles of Sp(6,2) (and 2.Sp(6,2)): P_0-cuts and labels; small modules of Sp(8,2) | 5 in 14, 6 in 17; mixed 1 + 4 in 10 (0.119), 7 + 38 in 88 (0.125); screens for 3 in 10, 2 in 7 | either a certificate beats the hook (closing block 3 in the plain games), or the above-bar mixed units become proof-immune construction targets |

Why this order. Task 1 is the one where both outcomes change the whole plan: it decides whether any known certificate
family can ever speak to the bar, and its lower-bound half is cheap and may prove "no lift closes 2 in 7" outright.
Task 2 is cheap and decides a cluster of three above-bar cases in one game. Task 3 has the largest prize and is the
only place theory says calculation must lead, but a null result changes less. Task 4 refines block 3 (below the bar
except for the mixed units, which also need a new recursion) and builds the tool for block 4.

Not chosen, with reasons: the block-4 exactness check (cheap, but a pure calculation); block-3 border P2 runs (needs a
quiet machine, weak link to the bar); literature T5 (tensorised block-2 family; predicted to fail); three-factor sweeps
in the widened game at block 4 (worth doing after Task 4 shows which modules carry labels).

Skeptic queue the orchestrator may want in parallel (claims the map leans on that nobody has checked): reps Theorems 1
and 2, Prop 3, Lemma 4.1; cohomology Theorem E+; spectrum Theorem C and Lemmas D, E; meta-border Theorem X; pf-duality
D3/D4 (the 18-row certificate can be checked by hand); 89/47 (literature); ex-laser Theorems 1 and 3; ex-dimension
Lemmas C, D.

---

## Sources

Paths are folders of the private research archive; they are not part of this repository.
- This round: `wht43/meta/BRIEF.md`; `wht43/meta/spectrum/` (NOTE.md, PROOF.md, STATUS.txt, liftprof.c, out_*);
  `wht43/meta/border/` (NOTE.md, EXACTIFY.md, STATUS.txt); `wht43/meta/cohomology/PROOF.md`; `wht43/meta/reps/`
  (PROOF.md, TABLE.md); `wht43/meta/literature/NOTE.md`; existence round `wht43/meta/ex-dimension/NOTE.md`,
  `wht43/meta/ex-laser/PROOF.md` and STATUS.txt, `wht43/meta/ex-legit/NOTE.md`, `wht43/meta/ex-v-legit-compute/`.
- Earlier: `wht43/proofs/NOTE.md`; `wht43/proofs/duality/PROOF.md` (§1–4, 6–9, 11) and STATUS.txt;
  `wht43/proofs/ladder/PROOF.md` (§3, §4); `wht43/wave16/b56-note/NOTE.md`; `wht43/wave16/b56-methods/STATUS.txt`;
  `wht43/wave16/tb4-note/NOTE.md`; `wht43/wave16/b4-cert/STATUS.txt`; `wht43/wave17/NOTE.md`; `wht43/wave17/xc3-note/
  NOTE.md`; `wht43/wave17/xc3-cert3/PROOF.md`; `wht43/wave17/xc3-search3/PROOF.md`; `wht43/wave17/ut-note/NOTE.md`;
  `wht43/wave17/ut-hunt/STATUS.txt` and bdlib.py; `wht43/wave17/d-G3-SPECTATOR/more-spectators/RUNLOG.txt`;
  `wht35/wider/w-bounds.md` (§3–5); `wht33/HANDOFF.md` (STATE blocks of 13–14 October).
- Arithmetic in this map (savings, lift thresholds, mixed tightness, τ_c values): one-off python3 runs on the Mac,
  under one second each; the formulas are s = ln(mW/M)/ln m, value of a lift = m/(1 + E/#edges), and the τ_c root of
  Σ W_k k^τ = M.
- Literature named here is cited as the source notes cite it; I read none of it today beyond those notes. DFZ:
  R. Dougherty, C. Freiling, K. Zeger, arXiv:1311.4601 (local copy `wht35/parallel/p-field/lit/arxiv_1311.4601/`).


---

## Abstract wave 2 (meta2, 14 October 2026): five more tracks, folded in

a2-map. Inputs: `../meta2/BRIEF.md`; the track folders `../meta2/a2-sg-lifts/`, `a2-vi-stability/`, `a2-gct/`,
`a2-zx/`, `a2-nag/`; the skeptic verdicts (folders `../meta2/v-a2-*` and the verdicts returned to the workflow). Sections
0–8 and Appendix A above are unchanged; read this section as a delta to them. Labels and status words as at the top of
this map. Counting rule here: a claim counts only if no skeptic refuted it; **upheld** = two skeptics tried and failed;
**unchecked** = no skeptic ran (the workflow sends only "proved"/"excluded" outcomes to skeptics, so a2-gct and a2-nag,
both "narrowed", have none). Unchecked claims are reported, not counted as established. Nobody has checked this section.
Nothing was posted.

### W2.0 The two decisive questions after this wave

1. **Proof side: can any known-type proof settle the 0.1 bar?** Still undecided, but reduced to one measurable quantity.
   - **Angle law** (A23, upheld): for any frame-family certificate, ratio at depth j = spectral factor × angle factor.
     The spectral factor is ≤ 2 − 2^(1−j) (the exact maximum for m ≤ 40; < 8/3 for every m). The angle factor is
     s_1/s_j, where s_d = mean sin² per nonzero principal angle at distance d. So settling the bar at block m ≤ 40
     needs angle factor ≥ m^0.9/2 (2.88 at m = 7, 3.97 at m = 10; derived here). Every computed family stays below 2:
     block-3 cuts ≤ 1.85 over 11,048 cuts; hooks ≤ 2n/(n+2) [proved]; start lifts fall from 1.66 to 1.45 by m = 10.
   - **Lifts** stay the only named candidate.
     - Lower bounds on their excess per edge are constants: ≥ 0.072 from m = 3, ≥ 0.075 from m = 4 (A20, upheld).
     - Sylvester–Gallai-type counting is capped by a bounded constant, so it can never show that lifts fail the bar
       (A22). Prop C's value of 1/3 is refuted; the bounded cap survives.
     - Constructions grow faster than the threshold [computed; evidence only].
   - **Other kinds of proof.**
     - GCT-type obstructions prove exactly the border complexity (A25, unchecked). They gain nothing in any single
       irreducible decided at block 3 (A26).
     - ZX proves nothing beyond the support count (A27, core upheld after a repair).
   - **Net.** No known-type proof is shown able, or unable, to settle the bar. Appendix A.1's trigger for moving proof
     effort to constructions (b_4 ≥ 36 for every lift) is not reached: 18 ≤ b_4 ≤ 36.
2. **Construction side: special position.**
   - In mw-sym4's symmetric block-4 locus the special position is not in the coefficients (A28, unchecked).
     652 of 652 target skeletons have Gröbner basis {1} over Q, and each obstruction is a sparse pattern certificate.
     A unit there would need discrete coincidences: letters that put another generator on the private junk frames.
   - GCT's relaxations show the same thing numerically (A26). Known solutions sit on lower-rank strata, and generic
     descent misses them: it fails on the plain 1 in 3, which provably exists.

### W2.1 Tracks, outcomes and skeptic status

| Track | Outcome | Skeptics | What counts |
|---|---|---|---|
| a2-sg-lifts | proved (in part) | logic: **refuted** (Prop C false); compute: not refuted (Prop C must be corrected) | **upheld by both**: Lemma E3 (tight), Lemma E21, Theorem G3 (E ≥ 68), Cor G3′ (b_4 ≥ 18), b_5 ≥ 141, link Lemma L, Lemma M, Theorem SG. **Refuted**: Prop C (ceiling 1/3 per edge; "from block 18"). Reproduced but not certified: b_5^A ≥ 86 (HiGHS branch and bound over an LP of 80) |
| a2-vi-stability | proved | logic and compute: not refuted | **upheld**: Theorem F, Lemma 2, angle law, Theorem H, near-miss arithmetic. The skeptics strengthened it: s_1 = 1 for every hook at every m, k; 68,320 exact over Q; point-additivity of hooks for m ≤ 40 |
| a2-gct | narrowed | none run | unchecked: Props 1, 2, Lemma 3; the computed relaxed schemes |
| a2-zx | excluded | logic: refuted as written; compute: not refuted | core (sum-free exclusion) holds after the f-uniform repair both skeptics supplied; Lemma 2 holds; Prop 3 holds in its restricted model; several statements corrected (W2.6) |
| a2-nag | narrowed | none run | unchecked: Lemma 2 (pattern certificates), 652 exact Gröbner decisions, homotopy calibrations |

### W2.2 New predictors (continuing §2's numbering)

| Id | Statement | Status | Where |
|---|---|---|---|
| A20 | Lift gluing. **Lemma E3**: an interval lift with a bad end edge has E ≥ 3 (tight: explicit exact lift, 16 subspaces). **Lemma E21**: one more pattern is impossible. **Theorem G3**: every whole-graph block-3 lift has E ≥ 68 (64 end and 96 middle incidences per edge, 128E ≥ 8640). **Cor G3′**: b_4 ≥ 18 for all lifts (ILP optimal; LP dual 120/7). Also b_5 ≥ 141 (LP dual 5611/40, using K4 = 18 rows) and b_5^dec ≥ 172 | upheld; b_5^dec rests on branch and bound | `meta2/a2-sg-lifts/PROOF.md` §0–3; `v-a2-sg-lifts-compute/` (counts3.py, check_e3.py, dualcert.py) |
| A21 | **Theorem SG**: decomposed lifts pay ≥ (2/3)κ_m − o(1) per edge; κ_m → 1/2 (even m), 0.4716 (odd m), i.e. 1/3 and 0.314. Via DSW Cor. 1.11 (Dvir, Saraf, Wigderson, arXiv:1211.0330) with binary rank. Effective only from m ≈ 110 | upheld | same, §4; `sgl-vc/out_vc.txt`; `v-a2-sg-lifts-compute/kappa.py` |
| A22 | **Counting ceiling (corrected).** Arguments that charge one failing incidence per bad line are capped by a bounded constant per edge, so they cannot show lifts fail the bar at all large m. Prop C's 1/3 and "block 18" are refuted. Lemma L's LP is capped at (2/3)κ_m: 0.400, 0.583, 0.299, 0.387 at m = 3..6, tending to 1/3, and below m^0.1 − 1 from m = 17. Integrality-aware versions of the same counting (per-link hit bound 3x − 2x²) prove ≥ 0.382 / 0.357 per edge asymptotically. The only exhibited integral ceiling is ≈ 1 per edge (alternating levels), which m^0.1 − 1 passes only at m = 1024 | bounded cap: argued by both skeptics; the values are the skeptics' own computations, unchecked | `v-a2-sg-lifts-logic/hdp.py` (not run); `v-a2-sg-lifts-compute/out_kappa2.txt` |
| A23 | **Theorem F + angle law.** For any frame family of equal-dimension subspaces on DSp(2m,2): φ_j/φ_1 ≤ T_m(j). T_m(j) = 2 − 2^(1−j), from the eigenvalue 2^m − 3; it is exact for j ≤ m ≤ 40 with a gap ≥ 1/4 for j ≥ 2, and T_m(j) < 8/3 for every m. Ratio_j = (φ_j/φ_1)(s_1/s_j). Coordinate families have angle factor 1, which re-proves κ_m < 2 for m ≤ 40 | upheld | `a2-vi-stability/PROOF.md` §1–3; `v-a2-vi-stability-compute/omega_exact.py`, `frames_check.py` |
| A24 | **Theorem H.** For the hook N in V ⊗ Λ^k V and any unitary g: ‖(1−P)GP‖_F² = C(n−2,k−1)(n − Σ\|g_xy\|⁴). Hooks saturate Theorem F at every depth, block and k. Their whole value is the angle factor, ≤ 2n/(n+2) < 2. Every hook has s_1 = 1 (90° edges) | upheld (s_1 = 1 added by both skeptics) | same, §6; `v-a2-vi-stability-compute/hook_cells.py` |
| A25 | **GCT Prop 1**: closed obstructions (GCT occurrence and multiplicity obstructions included) prove exactly the border complexity. One can be taken in a single isotypic component of GL_W(S) × GL_W(S), and B1 caps them. **Prop 2** (relaxation duality): an obstruction that sees only the irreducibles Σ (any copies, so any labels) is capped by the relaxed game in End(V_Σ). Irreducibles couple only through shared P_0-types | unchecked | `a2-gct/PROOF.md` |
| A26 | **GCT computations.** At block 3, relaxed complexity = best cut in 15, 21(−11), 27, 35(−5) and 56: exact schemes at the tight M, plateaus one move below, two seeds. In V_56: an exact 5 in 14 (independent least squares 5.6e-13; bounded coefficients) and an exact mixed 1 + 4 in 10 (3.4e-29). 56+105(−35) at 5 in 14: 5.4e-9 after 1,500 iterations and still falling, \|x\| flat at 242 (provisional). Undecided, because the solver fails its own special control (the plain 1 in 3): 35(+15), 56+70, 56+105(+25)+7. Block 2 (XW, and the closed games xq and thmB): the commutant of the free letters is the scalars at f = 1, 2, 3 | computed, unchecked | `a2-gct/STATUS.txt` (ANSWER); `a2gct-mac/log_56+105a_5_14_s2.txt` (END line) |
| A27 | **ZX Theorem 1**: a sum-free circuit or diagram for W copies needs ≥ mWf non-monomial one-qubit gates. With boundary Clifford letters (XQ, XW) the repair is p·f + 2mWL ≥ mWf for every f, so p ≥ mW. **Lemma 2**: gadget fusion with slot weights is cost-neutral. **Prop 3**: commuting slot weights split the scheme, proved only for slot-uniform data and rank-1 idempotent weights | core upheld after repair; Prop 3 restricted | `a2-zx/PROOF.md`; `v-a2-zx-compute/` |
| A28 | **NAG Lemma 2 (pattern certificates).** If row_K = Σ c_i(α) row_{J_i} identically, the frame-sum relaxation is infeasible for every α, and min residual² ≥ 1/(1 + Σ\|c_i\|²). Constant cofactors give the solver floors 1/2, 1/3, 1/4. Polynomial cofactors give the "borders": the residual → 0 only as α → ∞. In all but 2 of about 436 classified certificates the junk rows are private (touched only by generators that reach the target) | proved, unchecked; decisions computed-exact (Singular std over Q) | `a2-nag/PROOF.md`, `STATUS.txt`, `totals.txt` |

### W2.3 Changes to the ceilings (§1.2) and the lattice (§3)

Ceilings (derived here from A20 and c2-lifts Corollary M, λ_m ≥ 0.075 for m ≥ 4):

| Level | m = 3 | m = 4 | m = 5 | Source |
|---|---|---|---|---|
| lifts, every whole-graph lift | ≤ 2835/1013 = 2.7986 (**below the hook's 2.8**) | ≤ 4/1.075 = 3.721 | ≤ 5/1.075 = 4.651 | A20 |
| θ-shadow cap for lifts | 0.937 | 0.948 | 0.955 | a2-sg-lifts |
| spectral factor of any frame certificate, depth j | 2 − 2^(1−j): 1.5, 1.75, 1.875 at j = 2, 3, 4 | same | same | A23 |
| hooks, angle factor | ≤ 2n/(n+2) | | | A24 |

Lattice:

| Case (s) | Game | Change | Status |
|---|---|---|---|
| block 3, 5 in 14 (0.063), 6 in 17 (0.052) | E, Wd, XQ | **No lift excludes either**: they need E ≤ 67 and E ≤ 55, and Theorem G3 gives E ≥ 68. No argument confined to V_56 excludes 5 in 14 (exact relaxed scheme). First coupling 56+105(−35): provisional "no gain". 56+70 and 56+105(+25)+7: undecided | G3 upheld; GCT unchecked |
| block 3, mixed 1 + 4 in 10 (0.119) | E | already excluded by Theorem M (in 105(−35)). GCT: an exact relaxed scheme exists in V_56, so every exclusion must use other irreducibles | unchecked |
| block 4, 2 in 7 (0.096), private-end branch | E, Wd | k = 8 hook slack 640 = entirely depth-3 angle spread (s_3 = 48048/68320 = 0.70328; **every frame certificate that closes the branch needs s_3 < 0.7**). 68,320 is exact over Q. Decomposed lifts cannot exclude 2 in 7 (36/240 = 0.150 > 1/7). General lifts need b_4 ≤ 34; proved 18 ≤ b_4 ≤ 36. **562 symmetric skeletons ending privately: Gröbner {1}**; mw-sym4's border skeleton has a 2-row certificate | vi, sg upheld; NAG unchecked |
| block 4, 3 in 10 (0.132), 4 in 13 (0.150), 7 in 23 (0.142) | E, Wd | symmetric skeletons 81, 4 and 5: all Gröbner {1}; one more 7-in-23 skeleton undecided (Singular halted on 410,165 rows) | unchecked |
| block 5, commutative slice | E | `moonshot/mw2-b5beyond/comm/out_cslice_5_3_12.txt` has finished: 3 in 12 (0.139), 0 of 152,430,767 κ-in-span multisets pass. So the slice is closed above the bar for W ≤ 3 at block 5. W ≥ 4 is open (4 in 16, 0.139); there is no slice LP for m ≥ 6 | computed, unchecked (both zx skeptics read the file before it finished) |
| block 5, 2 in 6 (0.317) | E, Wd | §0 item 5 and §3.4 list it as open ("misses 3 by 0.008"). **Superseded**: mw-big56's two-layer bound M ≥ W + m excludes it | proved, no skeptic yet (moonshot NOTE-wave2) |
| block 2, everything-free and the closed games xq, thmB | XW, XQ | finite-dimensional GCT proves nothing at f = 1..3; obstructions must live on all-f representations | unchecked |
| the bar | all | A23: needs angle factor ≥ m^0.9/2; every computed family < 2; lifts undecided (A20–A22) | |

### W2.4 Ledger: predictions written before computing, and outcomes

| Track | Prediction | Outcome |
|---|---|---|
| sg-lifts P1 | SG constant ≈ 1/6–1/5 asymptotically, for D ≳ 40; nothing at m = 3–6 | held in shape: (2/3)κ_m → 1/3 or 0.314, effective from m ≈ 110 |
| sg-lifts P2 | SG/design-matrix bounds stay bounded, cannot show λ → ∞ (80%) | held (bounded cap); the specific ceiling 1/3 refuted |
| sg-lifts P3 | n_5 ∈ [45, 75], search finds ≈ 60 | held: 55 |
| sg-lifts P4 | decomposed b_5/e_5 ≥ 0.2; constructions grow (65%) | not proved (≥ 0.083); constructions consistent (0.279, 0.457) |
| sg-lifts P5 | block-3 gluing optimum > 54, guess 60–80 (70%) | held: ≥ 68 proved; true optimum in [68, 162] |
| sg-lifts P6 | end-edge obstruction at block 4 for decomposed optima | held qualitatively; no quantitative gap |
| sg-lifts P7 | provable θ-shadow of lifts stays > 0.9 for m ≥ 5 (75%) | held (0.937, 0.948, 0.955) |
| vi V1 | spectral ceiling 2 − 2^(1−m), attained at i = 1 (80%) | held; theorem (exact for m ≤ 40) |
| vi V2 | hook spectral ratio in [1.3, 1.75] | held: exactly 7/4 |
| vi V4 | stable families have bounded ratio (75%) | held (the range's lower end was already known: weak test) |
| vi V5 | front certificates separate the top-rank piece, always (60%) | held for both front certificates; "always" refuted (3 of 22) |
| vi V6 | start-lift edge excess at m = 7: 1.9 ± 0.3 | held (1.78, 1.93) |
| vi P7 | start-lift ratio stable in [2.35, 2.65] at m = 8–10 | **refuted**: 2.26, 2.18, 2.02; edge excess 2.55, 3.13, 3.94 |
| vi P8 | no cut breaks the ceiling; 7/4 reached in most modules | first half held (0 of 11,048); second refuted (12 of 26). Its data were collected before it was written |
| gct P1, P4 | GCT obstructions are closed; relaxation duality | proved, unchecked |
| gct P5 | 5 in 14 exists in V_56 (0.7); 5 in 13 fails (0.97) | held |
| gct P6 | small irreducibles have schemes at the cut bound | held in 15, 21(−11), 27, 35(−5); 35(+15) undecided |
| gct P7 | block-2 free algebra full at f = 1..3 (4,368) (0.85) | held |
| gct P8 | no gain in single irreducibles; any gain only in large sums (0.65) | first half held in 5 of 6; second open (couplings undecided; 56+105(−35) provisional "no gain") |
| zx P1 | sum-free cannot save, in any game | held after the f-uniform repair |
| zx P2 | PyZX calibration: invariant = mWf; scrambled median ratio ≈ 1.5 | held; median 1.167, max 1.833 |
| zx P3 | f = 1 degenerate; f = 2, 3 jointly not | held (residuals 0.612, 0.822) |
| zx P4 | commuting slot weights split; fusion under-counts; 0 hits | held (P4(a) only in the restricted model) |
| nag P1 | the border skeleton's relaxation has no finite solution (75%) | held: 2-row certificate |
| nag P2, P3 | calibrations solvable / empty | held (17 of 18 plain-cost feasible; negative controls empty) |
| nag P4 | at most 2 relaxation solutions missed by LM in ≥ 40 skeletons (80%) | held (0 of 652); the clause "LM floors are true minima" is refuted for polynomial cofactors |
| nag P5 | finite components fail the full system | untested (no finite component found) |
| nag P6 | no unit (99%) | held where tested; mixed and XW targets not built |
| nag P7 | no hidden special components in feasible relaxations | held (9 of 9) |

### W2.5 The block-4 near miss across the tracks

- **Certificate side** (vi, upheld). The hook's depth-3 ratio is 68320/27456 = 2.488 = 7/4 (spectral; the maximum for
  any certificate) × 1.4219 (angle). Its edges are already at s_1 = 1. So the 640 of room is exactly 320 directions
  per array short of 68,640, at the same spectral mass. Any frame certificate that closes the private-end branch has
  s_3 < 0.7.
- **Block-3 precedent** (gct, unchecked). Hook-tight targets (5 in 14, 1 + 4 in 10) have exact schemes inside the
  hook's own irreducible. 1 + 4 in 10 fell only to a certificate in a different irreducible (Theorem M, 105(−35)).
  By analogy, the block-4 branch would close, if at all, through other irreducibles of Sp(8,2). vi's own next step
  names the same place (hook relatives; cuts in the 85, 238, 510 and 595); Theorem S′ allows a depth-3 ratio of 2.5
  for cuts within 20 points of balance.
- **Lifts** (sg, upheld). They cannot settle 2 in 7 yet: decomposed lifts sit at 0.150 > 1/7, and general lifts sit in
  18 ≤ b_4 ≤ 36 against a need of ≤ 34. At block 3 every lift is below the hook.
- **Construction side** (nag, unchecked). All 562 symmetric 2-in-7 skeletons ending privately are infeasible for
  every complex α. mw-sym4's border is a polynomial-cofactor certificate, not a hidden solution. Open: asymmetric
  programs (w3-asym), discrete coincidences, and whether the full system carries a border there (needs a closure in C).
  If such a border is formally polynomial, px-border Theorem FE (two skeptics) gives an exact unit with the same
  number of moves, possibly of another shape.
- **Reading** [argued]. The evidence favours "the hook is slightly loose at block 4". A unit would have to be
  asymmetric and tight to 0.47 % against the hook throughout. 2 in 7 is below the bar either way. The same split
  applies to 3 in 10 and 4 in 13, which are above it.

### W2.6 Corrections found while folding

1. a2-sg-lifts Proposition C is false as written (both skeptics).
   - The constraints of Lemma L are not satisfied by φ ≡ 1/6: φ_0 = φ_{m−1} = 0 forces φ_1 ≥ ν_m/3.
   - "(1/3)κ_m ≤ 1/6" fails at every even m and at m = 3, 7.
   - The class of arguments is not capped at 1/3 (A22).
   - The headline "SG-type methods can never decide the bar, from block 18" becomes "cannot decide it at all large
     blocks; the crossing block is unknown (between about 17 and 1024 depending on the version)".
2. a2-sg-lifts: "the block-3 optimum is in [68, 101]" means the relaxed rules only. The true whole-graph optimum is
   only known to lie in [68, 162]. STATUS "17 subspaces" should read 16. The link-ILP ceiling "123" is ≤ 111 (the
   skeptic's incumbent).
3. a2-vi-stability:
   - "Spectral factor below 2 at every block" is proved for m ≤ 40; for larger m only < 8/3 is proved.
   - Statement (3) puts point-additivity in a [proved] bullet; PROOF.md has it [argued]. It is now proved for m ≤ 40
     by the compute skeptic (the image of the incidence map is E_0 + E_1).
   - "A better certificate is likely" compares depth 3 at block 3 (the opposite pair) with depth 3 at block 4 (not
     opposite). It also mixes two block-3 hooks (0.5625 in the 56; 0.625 for V ⊗ Λ⁴). It is weak support.
4. ms-fano's odd-k hook values (the "16 × χ = 0" rule) are wrong. At k = 7 the correct depth-1 and depth-3 ranks
   are 24,024 and 59,078, not 23,744 and 58,784. This was found independently by o3-ratio and by
   v-a2-vi-stability-compute. The conclusion is unchanged: no k = 5..10 hook has a depth-3 ratio above 2.488.
5. a2-zx:
   - The premise "every free letter is monomial" is false in XQ and XW (boundary Clifford letters). The proof needs
     the f-uniform repair (A27).
   - "378 circuits" means 162 circuits and 378 pass results.
   - The cost is the number of odd eigenvalues of the slot weight, not its rank.
   - "A unit needs non-commuting slot weights with phases or relabelling" is not established beyond Prop 3's
     restricted model and blocks 3–4.
   - F1 lives in `moonshot/theory-guide/RULES.md`, not in pf-logging.
   - Its "2 in 6 at block 5 misses by 0.008" near miss is superseded (W2.3).
6. Both zx skeptics reported `cslice 5 3 12` as unfinished. It has finished, with 0 passing (W2.3).
7. a2-nag P4: the LM floors are true infima only for constant-cofactor certificates. LM's other floors (0.176,
   0.096, …) are plateaus.
8. a2-gct: the pre-restart relax.py failed its gradient check. It and its outputs are void (`a2-gct/void-14h26/`).

### W2.7 Next, ranked by what it would change

1. **Block-4 certificate with s_3 < 0.7 at full spectral budget** (vi next 1, gct next 4). Look among hook relatives
   (closed forms via Theorem H's Cauchy–Binet) and among cuts in the Sp(8,2) irreducibles 85, 238, 510 and 595. A hit
   closes the private-end branch of 2 in 7 for that family. It is cheap, and two tracks point at it independently.
2. **Lifts at large m.** One direction: a rank-deficiency Sylvester–Gallai (a 2-flat version; DSW Theorem 1.14), or
   b_4 = 36 for all lifts plus a block-4 Lemma E3; either would kill lifts. The other direction: lifts optimised for
   small far-pair angles rather than edge excess (vi next 3).
3. **Discrete special position** (nag next 1). Build symmetric 2-in-7 and 3-in-10 skeletons whose letters put another
   generator on the junk frames of the pattern certificate. Screen each with the exact relaxation test (under a second
   per skeleton) and send only survivors to the full system. Alongside: asymmetric second-route programs (w3-asym,
   running).
4. **Border 2 in 7 in the full system** (nag next 2). This needs a closure in C. If a border exists and is formally
   polynomial, px-border Theorem FE turns it into an exact unit.
5. **Stratum-aware relaxation solver** (gct next 1–2). Validate it on the plain 1 in 3 in 35(+15) and 56+70, then
   decide 5 in 14 in 56+70 and 56+105(+25)+7. An infeasible coupled relaxation whose controls pass would be a closed
   obstruction beyond the hook, excluding 5 in 14.
6. **Skeptics.** GCT Props 1–2 and Lemma 3; NAG Lemma 2 and the validity of the frame-sum relaxation in Wd; a re-check
   of a2-zx after its edits.
7. **Not reached this wave.** Block-3 mixed 1 + 5 in 12, widened 2 + 9 in 22, and XW 3 in 5: mw-sym4's builder needs a
   per-array target.

### W2.8 Sources for this section

Paths are folders of the private research archive; they are not part of this repository.
- Tracks: `meta2/a2-sg-lifts/` (NOTE.md, PROOF.md, STATUS.txt); `meta2/a2-vi-stability/` (NOTE.md, PROOF.md,
  STATUS.txt); `meta2/a2-gct/` (STATUS.txt ANSWER block, PROOF.md, `a2gct-mac/log_56+105a_5_14_s2.txt`);
  `meta2/a2-zx/` (NOTE.md, PROOF.md, STATUS.txt); `meta2/a2-nag/` (STATUS.txt, PROOF.md, totals.txt).
- Skeptics: `meta2/v-a2-sg-lifts-compute/` (STATUS.txt, out_kappa2.txt, out_dualcert_m4/m5.txt);
  `meta2/v-a2-sg-lifts-logic/hdp.py`; `meta2/v-a2-vi-stability-logic/`, `meta2/v-a2-vi-stability-compute/`
  (STATUS.txt, out_omega_exact.txt, out_hook_cells.txt); `meta2/v-a2-zx-compute/STATUS.txt`. The verdict texts as
  returned to the abstract-wave2 workflow (the logic skeptics of sg-lifts and zx saved no report file).
- Read for context: `moonshot/o3/NOTE.md` (the near miss, odd-k correction); `moonshot/NOTE-wave2.md` (mw-big56
  two-layer bound, absorb4); `px/NOTE.md` (Theorem FE); `moonshot/mw2-b5beyond/PROOF.md` and
  `comm/out_cslice_5_3_12.txt`; `wht33/HANDOFF.md` (state of the parallel rounds).
- Arithmetic marked "derived here": one-off python3 runs on the Mac (lift ceilings 2835/1013 and 4/1.075, θ caps,
  thresholds m^0.1 − 1, m^0.9/2, savings s = ln(mW/M)/ln m).
