> **Working note, kept as written.** An AI agent wrote this during the theory-of-theories round of 14 October 2026 (the literature track). Private paths and machine names
> were removed and nothing else was changed. No person has checked it. Where it disagrees with this section's
> [README](../README.md), the README wins, because it follows the later skeptic verdicts. One skeptic later recomputed the 89/47 certificate (the cap 752 exactly over the rationals) and upheld it.

# Transplants from barrier theory (meta/literature, 14 October 2026)

One agent, about 06:10 to 08:40 EDT. No skeptic has read this. Nothing posted; no Lean. Labels: [proved], [computed],
[argued], [guess], [read]. "Given Theorem R" as everywhere in wht43. Files in this folder; heavy jobs went through
a memory-capped job runner on the Mac or ran on a Colab machine (two single-core jobs under 1 GB).

## What I read, and how much

| Source | Read |
|---|---|
| Christandl, Vrana, Zuiddam, "Barriers for fast matrix multiplication from irreversibility", arXiv:1812.06952 (Theory of Computing 17, 2021) | full text sections 1-2.3 (pdftotext); definitions quoted below |
| Christandl, Vrana, Zuiddam, "Universal points in the asymptotic spectrum of tensors", arXiv:1709.07851 (JAMS 36, 2023) | abstract |
| Alman, "Limits on the universal method for matrix multiplication", arXiv:1812.08731 (CCC 2019) | abstract (+ CVZ's account of it: same barriers, asymptotic slice rank) |
| Alman, Vassilevska Williams, "Limits on all known (and some unknown) approaches to matrix multiplication", arXiv:1810.08671 (FOCS 2018) | abstract |
| Ambainis, Filmus, Le Gall, "Fast matrix multiplication: limitations of the Coppersmith-Winograd method", arXiv:1411.5414 (STOC 2015) | abstract (2.3725 exact approach, 2.3078 for a class of variants) |
| Bläser, Christandl, Zuiddam, "The border support rank of two-by-two matrix multiplication is seven", arXiv:1705.09652 (CJTCS 2018) | abstract |
| Efremenko, Garg, Oliveira, Wigderson, "Barriers for rank methods in arithmetic complexity", arXiv:1710.09502 (ITCS 2018) | abstract (no super-linear / no > 8n tensor-rank bounds for 3-tensors from rank methods) |
| Beverland, Campbell, Howard, Kliuchnikov, "Lower bounds on the non-Clifford resources for quantum computations", arXiv:1904.01124 (QST 2020) | abstract (stabiliser nullity, dyadic monotone; catalysis often necessary) |
| Bravyi, Smith, Smolin, "Trading classical and quantum computational resources", arXiv:1506.01396 (PRX 2016) | abstract (low-rank stabiliser decompositions of tensor powers of magic states) |
| Peleg, Shpilka, Volk, "Lower bounds on stabilizer rank", arXiv:2106.03214 (Quantum 6, 2022) | abstract (Ω(n) exact; via quadratic forms over F_2 and directional derivatives) |
| Heimendahl, Heinrich, Gross, "The axiomatic and the operational approaches to resource theories of magic do not coincide", arXiv:2011.11651 (JMP 2022) | abstract |
| Howard, Campbell, robustness of magic, arXiv:1609.07488 (PRL 2017); Veitch, Mousavian, Gottesman, Emerson, arXiv:1307.7171 (NJP 2014) | abstracts |
| Ailon, "A lower bound for Fourier transform computation in a linear model over 2x2 unitary gates using matrix entropy", arXiv:1305.4745 | abstract (and its quotation of Morgenstern 1973) |
| Morgenstern 1973 (JACM 20), Schönhage 1981 (SIAM J. Comput. 10), Coppersmith-Winograd 1990 (JSC 9), Strassen 1988 (Crelle 384) | NOT read; used from memory and from in-house notes (wht33/research-2026-10-09/readers/lit-lower-bounds.md). Statements used: Schönhage's τ-theorem Σ(k_i m_i n_i)^{ω/3} ≤ r and his direct-sum example; CW's use of a non-matmul helper tensor; Strassen's duality. Marked [read: memory] where used. |

## The dictionary (what corresponds to what)

| Matrix multiplication / magic | Our game |
|---|---|
| tensor rank R(T) / T-count | paid moves M |
| direct sum T ⊕ T' | more arrays (W) |
| tensor product (n,n,n)⊗(n',n',n') | substitution: block m inside block m' (block mm', rates multiply) |
| ω = log_n R̃(⟨n,n,n⟩) | τ = 1 − saving = lim log_m(M/W) |
| ⟨a,b,c⟩ with volume abc, τ-theorem Σ vol^{ω/3} ≤ r | rank-k target (k rows transformed, the rest spare) has "volume" k; mixed unit gives Σ W_k k^τ ≤ M (T1) |
| flattening / Koszul rank method | Theorem R certificate (f = a rank) |
| border rank | approximate (border) units; Theorem R bounds them too (pf-ladder B1) |
| spectral point (monotone, additive, multiplicative, normalised) | a certificate that is additive over arrays AND multiplicative under substitution: none known in characteristic 0 |
| Clifford gates free, T gates paid | monomial Cliffords (the stabiliser P of the Z-frame L_0) free, Hadamard directions paid |
| stabiliser groups / stabiliser states | frames (135 at block 3) / signed frames (1080 at block 3) |
| stabiliser nullity (Beverland et al.) | frame distance to L_0 = the mod-2 perfect certificate |
| stabiliser rank (linear combinations of stabiliser states) | cross-array mixing (linear combinations of arrays) |

## The six transplants, most direct first

### T1. Schönhage's τ-theorem and asymmetric direct sums → MIXED-TARGET UNITS. Applies directly; tested.

(a) **Abstract statement.** A unit at block m that gives W_k arrays a rank-k target (k rows transformed, m − k spare rows)
with M paid moves yields every exponent τ with Σ_k W_k k^τ ≤ M (the analogue of Σ vol^{ω/3} ≤ r) [argued: the
demand-flow schedule below; the simple "conversion" bound τ ≤ log_m((M + Σ W_k (m−k))/ΣW_k) is [proved], since a
rank-k demand plus m − k single moves is a kernel]. Theorem R is additive over targets, so for one certificate c the
bound is M ≥ Σ_k W_k r_k(c) with r_k(c) = f_c(rank-k target)/cap_c, and a direct sum of certificates never beats its
best summand for a fixed target mix (pf-duality Lemma 1.2). **Schönhage's lesson:** if different certificates are best
for different targets, a mix can be less excluded than every pure target, and the τ-theorem converts that slack into
exponent room. That is the asymmetric-direct-sum effect behind his ω < 2.55 [read: memory].

(b) **Prediction.** At block 3 (3 rows, widened game) the Theorem-R method leaves more room for mixed units
(kernel + rank-2 targets) than for either pure target. A new construction target: "1 kernel + 4 rank-2 targets in 10
moves".

(c) **Calculations.**
1. *Does the τ-analog hold?* (`tau_schedule_lp.py`, `out_tau_schedule_lp.txt`). This is an LP over per-entry demand
   flows. A k-bit demand may be split plainly (k → a + (k−a)) or fed to the mixed unit, which consumes W1 demands of
   size 3s and W2 of size 2s and returns M of size s. Lower bound for every schedule: n^{τ_c} (the potential Σ mass·k^τ
   cannot fall) [proved]. Results:
   - 1 kernel + 6 rank-2 in 14 moves: τ_c = 0.9108. LP slopes 0.9086, 0.9085, 0.9088 for n = 60 to 480, and
     cost/n^{τ_c} = 1.074, 1.072, 1.071, 1.069. So the exponent is achieved, against 0.9556 for the conversion
     bound [computed].
   - 2 + 12 in 28: identical [computed].
   - Pure control 5-in-14: cost/n^0.937 stays between 1.06 and 1.11 [computed].
   - Toy 1 + 1 in 4 (τ_c = 0.7605, conversion 0.834): slopes 0.788, 0.794, 0.800, still rising at n = 480, so not
     settled. A run with restricted splits (`tau_schedule_lp2.py`) was too constrained to say anything.
   - The two block-3 units ranked below (`tau_schedule_lp3.py`, `out_tau_schedule_lp3.txt`, run on a Colab machine):
     1 + 4 in 10: τ_c = 0.8811; slopes 0.8759, 0.8776, 0.8819; cost/n^{τ_c} = 1.096, 1.092, 1.089, 1.090.
     7 + 38 in 88: τ_c = 0.8750; slopes 0.8714, 0.8715, 0.8723; cost/n^{τ_c} = 1.105, 1.102, 1.099, 1.097.
     The conversion bound for both is 0.9372. So in the regime that matters (about 5 rank-2 arrays per kernel array)
     the τ-analog is achieved [computed, n ≤ 480; this is evidence for the schedule, not a proof].
2. *Certificate profiles at block 3* (`venn_sweep_m3.py`, `out_venn_sweep_m3.jsonl`, `mixed_room.py`). I swept all 180
   two-factor Venn certificates Λ^a ⊗ Λ^b, J at m = 3, with values on rank-2 targets with one spare row and on the
   kernel (mod-p ranks: lower bounds on f(target); exact caps). Controls reproduce the hook (80, 144, 224) and wave 17's
   Λ²⊗Λ³ J={1} (264, 488, 664). The Pareto front of (kernel ratio, rank-2 ratio):
   (2.8, 1.8) hook; (2.5152, 1.8485) Λ²⊗Λ³ J={1}; **(2.2766, 1.8936) Λ³⊗Λ⁴ J={0,2}: cap 752, f(rank-2) 1424, f(K) 1712.**
3. **Side result: a better rank-2 bound at 3 rows.** Rank-2 targets with one spare row need at least 1424/752 = 89/47
   = 1.8936 moves per array in the widened ("full") game [computed]. It is certified mod p by two codes: wave 17's
   venn_mask.py (p = 999983) and my own `indep_venn_check.py` on a Colab machine (unnormalised Hadamards, compound matrices,
   two primes near 2^31: cap 752 and 752, matching the closed form; f(rank-2) 1424 and 1424; f(K) 1712 and 1712;
   `out_3_4_J02.txt`). The same independent code reproduces the two known controls exactly: hook 80/144/224
   (`out_1_4_J1.txt`) and Λ²⊗Λ³ J={1} 264/488/664 (`out_2_3_J1.txt`). This
   **excludes 7 in 13** (13/7 = 1.857), which wave 17 listed as open (worth 1.07e-1), **and 9 in 17**. The pure rank-2
   room at 3 rows falls from 0.114 to 0.0789. Smallest open pure rank-2 cases: 10 in 19 (0.074), 11 in 21 (0.067).
   It needs a skeptic.
4. **Mixed room** (`mixed_room.py`). The smallest τ not excluded by the three front certificates is **τ = 0.8745
   (saving 0.1255)**, at W2/W1 ≈ 5.6. The pure kernel allows 0.0628 and pure rank-2 now 0.0789. Best open mixed units:
   7 kernel + 38 rank-2 in 88 moves (0.1250); 4 + 21 in 49 (0.1243). With at most 14 arrays: **1 kernel + 4 rank-2 in
   10 moves (0.1189)** and 2 + 8 in 20.
   - "1 + 4 in 10" sits exactly on the hook (2.8 + 4 × 1.8 = 10), like 5-in-14 (5 × 2.8 = 14), so the t514 tightness
     lemmas (all moves rank 1, each +80) should carry over.
   - It implies 5-in-14 by conversion [proved], so it is at least as hard. If the τ-analog holds, it is worth 0.119
     instead of 0.063.

(d) **What this means.** At block 3, the cheapest place to look for a unit is a mixed one. The cheapest proof that
closes it needs ONE certificate that is good on both targets at once: kernel ratio ≥ 2.8 and rank-2 ratio ≥ 1.85 on the
same N. At block 4 the known certificates show no gap: hook k = 8 dominates on rank 2, 3 and 4 (1.769, 2.488, 3.154). So
there the mixed room equals the best pure room, rank 2 with 2 spare rows (0.177). A block-4 sweep like calculation 2
might move that number, as it did at block 3.

### T2. Barriers for rank methods (Efremenko-Garg-Oliveira-Wigderson 2017) and border rank (Bini, Schönhage) → a universal ceiling for Theorem R. Applies directly; the test is specified, not run.

(a) Theorem R is a rank method: f = dim(N + gN) − dim N is a matrix rank, and ranks are lower semicontinuous. So
Theorem R also bounds approximate (border) units; this is in-house barrier B1 (pf-ladder). EGOW prove that rank methods
cannot give super-linear (or > 8n) tensor-rank bounds for 3-tensors [read: abstract].

(b) Predictions:
1. BEST_m (the best Theorem-R bound per array, in the sense of pf-duality) is at most the border cost at block m.
2. As for flattenings, there is a universal constant ceiling, consistent with every family we know stopping below 4
   [guess]. If so, Theorem R can never exclude a saving of 0.1 from block 5 on (that needs > 4.26).
3. Wherever Theorem R is stuck exactly on a target, expect an APPROXIMATE unit there: t514-cert found no certificate
   above 14/5, so predict a border 5-in-14, and a border "1 + 4 in 10" [guess].

(c) Smallest tests:
- Re-use xc-search's least-squares machinery with ε-scaled coefficients on 1 + 4 in 10 at block 3 (5 arrays). The
  signature of a border unit is a residual that falls as the allowed coefficient size grows.
- On the proof side, an upper bound on BEST_3 itself: the char-0 LP over interval families (pf-duality's interval
  relaxation plus Ingleton/DFZ char-0 inequalities; dfz_check.py exists there).

### T3. Morgenstern 1973 and Ailon 2013 (potential = determinant or matrix entropy) → "adelic" necessary conditions. Proof sketch here; no computation needed.

(a) **Statement** [argued; Ailon's argument adapted to arrays, with f → ∞]. Suppose every step of a correct word is
unitary up to ONE global scalar: paid C^(⊗f), monomial letters, and cross-array mixing by a unitary on all arrays. Then
M ≥ mW.
- *Proof.* Let T be the map from inputs to all arrays. It stays an isometry, so every row has norm ≤ 1 and each array's
  block has Frobenius mass ≤ 2^{mf}. Φ(T) = −Σ|T_ij|² log|T_ij|². A paid move is f layers of 2×2 unitaries inside one
  array, so it raises Φ by ≤ f·2^{mf}·log 2. Monomials leave Φ unchanged. A unitary mixing of A arrays raises Φ by
  ≤ W·2^{mf}·log A, which does not grow with f. The delivered kernel has Φ = W·2^{mf}·mf·log 2 and Φ starts at 0. So
  M ≥ mW − (F_mix·log₂A)/f for every f; the unit is fixed, so M ≥ mW. The same holds for approximate units (Φ is
  continuous). ∎
- The determinant version (Morgenstern) says nothing: det C = i, so the kernel needs no volume. The published saving
  escapes the n log n bound with ONE f-dependent final scalar 2^{mf/2} (the in-house notes say the same).

(b) **Prediction.** Every saving unit has (i) a non-conformal step: some mixing that is not a scalar times a unitary
on all arrays, e.g. a copy into scratch, a shear, or a butterfly on 2 of more than 2 arrays; and (ii) a division by 2
(the mod-2 perfect certificate, w-bounds §4; proved for Clifford free steps). These are conditions at the real place
and at the place 2. Neither implies the other: a Pythagorean mixing (1/5)[[3,4],[4,−3]] on all arrays is conformal
with no halving.
- Consequence for methods: norm and entropy monotones (Ailon, robustness of magic, mana) can never certify anything
  beyond "the unit is non-unitary". Only rank-type (T2) or valuation-type certificates survive arbitrary mixing.
- A p-adic Morgenstern at p ≠ 2 gives nothing: Theorem F forbids perfect certificates whenever 2 is invertible.

(c) Test: a necessary-condition filter for construction searches. Reject candidate words whose mixing is all
conformal, or all 2-integral. One can also confirm on the published unit's certificate file (gcert1-p11) that both
conditions occur, which they must.

### T4. Magic-state resource theory (stabiliser nullity, rank, extent; catalysis) → certificates on signed frames; catalysts cannot help. Good fit for the dictionary; one result derived.

(a) Our free group (Paulis, affine maps, diagonal Z4 phases) is the Clifford stabiliser of the computational frame,
and a direction is a Hadamard. So our game is Hadamard-counting with monomial Cliffords free and linear combinations
of arrays free.
- The frame distance to L_0 (the mod-2 perfect certificate) is Beverland et al.'s stabiliser nullity relative to the
  Z-frame.
- Theorem F (no char-0 perfect certificate) is the reason the analogue of stabiliser rank is hard.
- Frame families (Averaging Lemma) see only stabiliser GROUPS, which is exactly Pauli-blind. The hook (2.8) beats
  every frame family found so far (best 2.561; lifts are capped at 2.917), and it sees signs.
- Peleg-Shpilka-Volk prove stabiliser-rank lower bounds with quadratic forms over F_2, the same tool as wave 17's
  spectator certificate (comm_cert.py). Even so, those bounds are only linear while the truth is believed exponential.
  The magic literature therefore predicts that rank-type lower bounds stall far below the truth. That matches our
  ceilings.
- Heimendahl-Heinrich-Gross and the known non-multiplicativity of stabiliser extent say not to expect per-array
  certificates to tensorise.

(b) Derived here [proved, three lines]: **catalysts do not beat Theorem R.**
- A catalyst is an array that starts holding g·y_aux and must end holding the same thing.
- The potential of Theorem R's proof starts at f(g) instead of 0 and ends at W·f(K) + f(g), because the end map is
  block-diagonal. So Σ f(paid) ≥ W·f(K) still holds.
- This matters because OpenAI's network "restores auxiliary values" (network.tex Remark, quoted in
  research-2026-10-09/readers/paper-explicit-saving.md). Catalysis is allowed but buys nothing against any
  certificate.

(c) Construction method we have not tried: certificate families indexed by SIGNED stabiliser states (1080 at block 3)
instead of frames (135), averaged as in the Averaging Lemma. This sits between Pauli-blind lifts (≤ 2.917) and the
hook. Smallest test: block 2 (60 signed states). Check that the averaging lemma's transitivity holds on pairs of
signed states for the free group, then anneal a lift as w-bounds §5 did. Under chain maps (cross-column) this family
is forced back to Pauli-blind, so it helps only in the plain games.

### T5. Asymptotic spectrum, irreversibility, universal-method barriers (Strassen; CVZ 2018/2019; Alman 2019; Alman-Vassilevska Williams 2018; AFLG 2015) → which barriers can exist at all. Framework; partly vacuous.

(a) CVZ define the relative exponent ω(t, s) = lim (1/n) min{m : t^⊗m ≥ s^⊗n}. It satisfies the triangle inequality
ω(s,t)·ω(t,u) ≥ ω(s,u). Irreversibility is i(t) = ω(⟨2⟩,t)·ω(t,⟨2⟩) = log R̃(t)/log Q̃(t). Any method through an
intermediate t gives a bound ≥ 2·i(t), because ⟨n,n,n⟩ → numbers has rate 1/2 (Strassen's asymptotic subrank n²)
[read: arXiv:1812.06952 §1.3, Def. 4-8].
- Our analogue: rates under substitution, with K_1 (one move) as the unit resource and τ = rate(K_1 → K_m).
- The barrier needs a large "subrank" of the target. For us that is the number of independent single moves one
  kernel call can stand in for. A kernel call delivers one array's transform, and leftover directions cannot be
  undone for free, so the subrank-analogue is O(1) [argued]. The CVZ barrier for the kernel itself is therefore
  vacuous.

(b) What does transfer:
1. **Barriers for HELPER methods** (OpenAI's banks/climbs, CW-like helpers). For a helper task, a method's exponent is
   capped by any monotone evaluated on the helper. T1's mixed-room computation is exactly this, with additive Theorem-R
   monotones and rank-k helpers. Result: τ ≥ 0.8745 at block 3 for every helper built from rank-2 and rank-3 targets on
   3 rows. That is an AFLG-style number for that class.
2. **Asymptotic barriers need monotones that are multiplicative under substitution.** Ours are not; that is why every
   family has a bounded ceiling (b56-note). The only spectral point we know is the mod-2 frame distance, which is
   exactly multiplicative (rank k ⊠ rank k' = rank kk'), but valid only for 2-integral words.

(c) Smallest test: does any Theorem-R certificate tensorise along the substitution 2 ⊠ 2 = 4? Take the block-2
perfect frame family and its product on block-4 frames that split as 2 ⊗ 2, average over the block-4 free group, and
compare with 4 = L_2². Prediction: it fails (< 41/13), because the block-4 free group mixes the two factors. A pass
would mean "block 4 cannot save", a striking result.

### T6. Coppersmith-Winograd / laser method and support rank (Bläser-Christandl-Zuiddam) → construction by helpers; support relaxations useless for H. Weakest fit.

- CW start from a cheap helper that is not matrix multiplication and extract many disjoint matrix-multiplication
  pieces from its powers. For us, the helper is a mixed unit and the extraction is the multinomial class of T1's
  schedule (powering). The tool already exists: T1's LP is a laser-method calculator for helpers with rank-k targets.
  New helper candidates to rate the same way: units with rank-1 "taxis", units for the zeta kernel feeding the WHT.
- BCZ's support rank keeps only the zero pattern. The support of H^⊗f is everything, so a support relaxation of our
  target is trivial (any dense matrix, for instance the all-ones one, costs O(N)). It cannot give lower bounds for H
  [argued]. It is meaningful only for sparse kernels: in OpenAI's finite-win game for the zeta kernel, a support
  version would test whether Yates's 12 calls on 8 entries survive arbitrary weights. That is a cheap combinatorial
  check if anyone wants a short proof of the 12.

## Ranked summary for the orchestrator

1. **T1 (mixed targets, τ-theorem)**: new targets 1 + 4 in 10 (0.119) and 7 + 38 in 88 (0.125) at block 3; the
   τ-analog was achieved to 0.002 in the relevant regime; side result 89/47 for rank-2 at 3 rows (7-in-13 and 9-in-17
   now excluded; needs a skeptic).
2. **T2 (rank-method barrier)**: predicts a universal ceiling for Theorem R and border units where Theorem R is tight.
3. **T3 (Morgenstern/Ailon)**: unitary words cannot save [argued proof above]; conditions at the real place and at 2.
4. **T4 (magic)**: dictionary; catalysts cannot help [proved]; signed-stabiliser-state certificates untried.
5. **T5 (spectrum/irreversibility)**: barriers for helper methods only; multiplicative monotones are what is missing.
6. **T6 (CW, support rank)**: T1's LP is the laser calculator; support rank is useless for H.
