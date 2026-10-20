<!-- header: existence-routes -->
# The skeptic's record: price reformulation

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** This is an AI skeptic's record of its check, as it wrote it. Its verdict is the status this revision uses for the track it checked, unless a later verdict is named in the section README.

---

# Skeptic's check of ex2-reform (changing the formulation)

14 October 2026, evening. I checked the reformulation track's claim with my own code and against the OpenAI source.
Nothing was posted. Labels: [checked] = rebuilt or re-derived by me; [read] = checked against the source text;
[unconfirmed] = not written out in enough detail to check.

## Verdict in short

- **Calibration: holds only in the weak form the brief allows this track.** The dictionary is sound. Lemma 2 is a
  correct re-derivation, and the plain FFT and OpenAI's s > 0 do look natural in the price language. But the method
  certifies no unit and refuses no impossible case. The "14 recorded unit values reproduced" are evaluations of the
  defining formula, and 6 of the 14 shapes are proved impossible.
- **Proposition 3 and the quantitative obstruction (O2) are refuted as stated.** Their hypothesis H_θ cannot hold for
  any θ < 1: two of its laws give one matrix two different prices. So the "ceiling of 0.036 for OpenAI's mechanism",
  the "dense-pair price 1.04, 35% cheaper" target and the "gap factor about 1.94 at d of about 120" follow from no
  possible θ-price. The same goes for next step 2's "about 0.036 if the pivots transfer".
- **One numeric statement is wrong.** At θ = 0.9 with OpenAI's pair price, the best d is 2 (short by 7.47), not about
  120. At d = 120 the shortfall is 34.7.
- **Proposition 1 (the duality) is neither refuted nor confirmed.** Its realisation step for several kernels is not
  written out. It stays [argued].
- **Lemma 3a, the small-d sparsity check and the d > 30 calibration are right.** The sparsity result holds in a
  stronger form.
- **No existence claim was made, and none is affected.**

## 1. Calibration

**Part A (the 14 values) [checked].** With one kernel the price cancels, and the condition is W m^z = Σ r^z. That is
the formula the recorded values were computed from (Theorem C; w-engine P2). Re-evaluating it reproduces the numbers,
which only checks arithmetic. I recomputed all 14 (`out_check_prop3.txt`, section 4) and looked up what each shape is:

| shape(s) | saving the formula gives | status |
|---|---|---|
| plain 4 in 20 at block 5 | 0 | a real unit |
| 2 in 5 at block 3, groupings 1,1,1,1,1 / 2,1,1,1 / 2,2,1 | 0.166 / 0.219 / 0.327 | **proved impossible** (MAP: "closed by 14W/5") |
| 3 in 8 at block 3, eight rank 1 / 2,2,2,2 | 0.107 / 0.291 | **proved impossible** (same line) |
| mixed 1 + 4 in 10 at block 3 | 0.119 | **ruled out in the engine game** (meta2, upheld) |
| 2 in 3 and 3 in 5 at block 2 | 0.415, 0.263 | stand-ins from another game, not engine units |
| 5 in 14; 7 + 38 in 88; 2 in 7, 3 in 10, grouped 2,2,2,1 at block 4 | 0.063 … 0.168 | open; nobody has them |

So the price language gives a positive saving to six shapes that are known to be impossible. That shows it is blind
to existence. It does not calibrate anything. The track says the language is blind to per-block impossibility, but
it still calls these shapes "recorded unit values".

**What does hold [checked / read].**
- Lemma 2 is right. "Right for every f" equals "right in every irreducible representation of the letter group", by
  Burnside-Brauer and Wedderburn. Lemma 2 needs the letters to lie in a finite group. They do, paid blocks included:
  w-engine section 3.1 says "the matrices form a finite group".
- The plain FFT is the exact tensor law, p(H^(x)m) = m 2^(m−1), the butterfly count.
- OpenAI's s > 0 is "no 1-price exists" (their Theorem 4.1 and section 8).
- OpenAI's amplification exponent is the same formula as Part A. Their Lemma 2.2 condition, g ρ^α < 1 with
  ρ = 1/(b q^(b−1)), gives α = ln g / ln(b q^(b−1)).
- Part C does reproduce "incompatible when d > 30": the first contradiction comes at d = 31.

That is the weak calibration the brief allows the reformulation track (a dictionary plus one known fact re-derived).
The strong calibration (certify units, refuse impossible cases) is not met, and cannot be met in this language.

## 2. Proposition 3 rests on an inconsistent hypothesis [checked]

H_θ assumes three things:
- every law in OpenAI's Theorem 4.1 holds for some θ-price;
- the tensor law is replaced by homogeneity: p(X^(x)k)/|X^(x)k| = k^θ p(X)/|X|;
- exchanges carry tensor-power prices.

Theorem 4.1 includes the **exact direct-sum law** p(A ⊕ B) = p(A) + p(B), and p vanishes on monomials
(`04-packing-prices.tex`, lines 12-36). Take K = U ⊕ I_2, of width 4.
- By the direct-sum law, p(K) = p(U) = 1.
- K ⊗ K is a permutation conjugate of (U⊗U) ⊕ (U⊗I_2) ⊕ (I_2⊗U) ⊕ I_4. I checked this as an explicit 16 x 16
  permutation.
- Direct sums plus homogeneity on U give p(K⊗K) = 2^(1+θ) + 4.
- Homogeneity on K gives p(K⊗K) = 4·2^θ.
- These agree only when 2^θ = 2, that is, θ = 1. At θ = 0.9 they give 7.732 and 7.464.

The same holds for any non-monomial A in place of U. So for every θ < 1, no function satisfies H_θ.

What this does to the claim:
- "Even granting every OpenAI law at exponent 0.9, their gap pair reaches only 0.964" is not meaningful. Granting
  those laws refutes 0.9-prices outright, but a real 0.9-price cannot satisfy them anyway.
- The numbers 0.036, 1.04 (35% cheaper) and about 1.94 are artifacts of the inconsistent law set. They are neither
  ceilings nor targets for any valid version of OpenAI's argument at θ < 1.
- The direct-sum law is used on both sides of OpenAI's reflection comparison:
  - the cheap side: one dense pair ⊕ I_(w−2) per supported pair;
  - the pivot and pair-map costs 4w and 2w;
  - p(M) = 2p(B).

  So a θ < 1 version must replace this law as well as the tensor law.
- The right replacement is visible in OpenAI's own amplification lemma. A sparse call A ⊕ I on q of Q coordinates is
  charged by dilution, (q/Q)^θ per coordinate, not by the sum law's q/Q.
  - For a finite win (2, 2, 3) the sum law would read a saving of 0.415, while amplification yields 0.208
    (`out_check_prop3.txt`, section 5).
  - At width 2^20 the two charges differ by a factor of 3.7. At OpenAI's widths (about 4^120) the factor is about
    10^7.
- Lemma 3a is a correct substitution [checked against `07-pivots.tex` lines 106-146]. It is a symptom of this; the
  sum law is a second, independent break.
- Next step 2 promised "about 0.036 if scalar pivots transfer and the other laws follow". The other laws include the
  sum law, which provably cannot follow. As stated, that promise is void.

**The arithmetic itself is right [checked, own integer-d scan].**
- θ = 1: first contradiction at d = 31.
- θ_c = 0.96392, with best d = 262.
- The needed pair price at θ = 0.9 is 1.03884, at d ≈ 122.

**One misstatement.** The claim says "at θ = 0.9 the best d (about 120) falls short by 7.5". With OpenAI's pair price
1.6077, the best d is 2 (short by 7.474). At d = 120 the shortfall is 34.7. The d of about 122 is the optimum only for
the lowered price 1.0388. Also, "gap factor 1.94" is not reproduced. The closest reading, 2/β_needed, gives 1.925.

## 3. Proposition 1 (the θ-price duality) [unconfirmed]

- (a), the potential argument, is right as written.
- (b), Motzkin's theorem with B = I and C = −A, is right. Note that it gives a mixture that is strict in only one
  type. OpenAI's packing lemma needs strictness in every type. With conversions in both directions, the strictness
  can be spread by a small perturbation, so this is not a gap.
- (c), the realisation for several kernels, is the substance, and it is not written out. ex-laser's Theorem 1 covers
  one kernel with mixed sizes, in its job model. With several kernels, a program must hold pending tasks of the right
  types and sizes at the same time, and must re-use pass-through tasks. I found no counterexample and no proof.
- (d), Tychonoff, is right given (c).
- Proposition 1's θ-prices are per-coordinate prices on tasks K^(x)n. They are forced to be homogeneous, and they obey
  no direct-sum law. So, by section 2, they are not OpenAI's prices with θ inserted. The identification "at θ = 1
  this is OpenAI's price theorem" holds only loosely: OpenAI's comparisons are fixed words with sparse gate calls.

## 4. The small-d sparsity check [checked, stronger]

I rebuilt OpenAI's Γ, J, L and A0 = (J⊗J)L from `08-reflection.tex` and checked their structural facts:
- Γ² = I;
- nnz(Γ) = (1 + d/2)N;
- A0 is symmetric, with nnz = (d + 2)w and h = (d/2 + 1)w;
- ker L has dimension w/2, rank A0 = w/2, and A0 E = 0.

Then I computed, in exact arithmetic mod two primes, the space of real symmetric matrices that are supported on A0's
off-diagonal pairs (diagonal free) and kill E.
- Its dimension is 1 at d = 1 and at d = 2. A mod-p dimension bounds the rational one, so A0 is the only such matrix
  up to scale.
- So removing any pair forces A = 0, and the support is removal-minimal (24 and 128 pairs), as claimed.
- The track's random greedy orders were redundant. A sparser matrix with a different support is still not excluded.

## 5. The other languages

These obstruction claims are judgments, not theorems.
- The quadratic-form point is correct: y^T H x with H invertible is a change of variables away from Σ x_i y_i.
- "Only dimension counts apply" in the group-ring language overstates it. Fixed-shape decision (Nullstellensatz) is
  also available, and the note itself says so.
- I found nothing to refute in the polynomial-multiplication or semiring remarks.

## Files

In `the archive's checks/wht43/exist2/v-ex2-reform/`:
- `STATUS.txt`: predictions, written before checking, and results;
- `RUNLOG.txt`;
- `check_prop3.py` and `out_check_prop3.txt`: Proposition 3 arithmetic, the H_θ consistency test, Part A with the
  status of each shape, and a finite win priced two ways;
- `sparse_check.py`, `out_sparse_d1.txt`, `out_sparse_d2.txt`.

## Sources

- The claim: `exist2/ex2-reform/` (NOTE.md, PROOF.md, STATUS.txt, RUNLOG.txt, calib.py, out_calib.txt,
  sparse/greedy.py, sparse/out_greedy.txt).
- OpenAI, "Finite tensor savings and exact Fourier circuits", LaTeX in
  `OpenAI's openai-math repository, preprints/Finite-tensor-savings-and-exact-Fourier-circuits-September-25-2026/build/sections/`:
  - `04-packing-prices.tex` lines 1-135 (Theorem 4.1 laws; comparison (4.2); packing lemma);
  - `02-generation.tex` lines 210-306 (tensor amplification and the sector decomposition);
  - `07-pivots.tex` lines 106-146 (one-coordinate lemma) and 275-330 (dense pair);
  - `08-reflection.tex` (all).
- In-house:
  - `wht35/wider/w-engine.md` sections 3.1 and 4.2;
  - `wht43/meta/MAP.md` lines 34-36, 179-193;
  - `wht43/meta2/NOTE.md` lines 13, 67, 80, 191;
  - `wht43/meta/ex-note/NOTE.md` lines 119-123;
  - `wht43/meta/ex-laser/PROOF.md` section 2 (statement only).
