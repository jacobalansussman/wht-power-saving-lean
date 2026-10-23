# The two skeptics' reports on "Induction across block sizes"

> **Reports kept as written.** Two AI agents that had no part in writing
> [induction-across-block-sizes.md](induction-across-block-sizes.md) were each asked, on 13 October 2026, to refute
> it: one by checking the logic step by step, one by recomputing every number with its own code. Private paths were
> replaced by public ones and nothing else was changed. No person has checked these reports. Their programs and
> outputs are in `../programs/engine-bounds-blocks-5-and-6/skeptic-checks/` and
> `../data/engine-bounds-blocks-5-and-6/skeptic-checks/`. Notes they name: "w-bounds" section 1 is where correctness
> of a unit is shown to be formal (by a Vandermonde argument), section 3 is the Averaging Lemma (here: Lemma 2 of
> [two-arrays.md](../../../notes/network-coding/two-arrays.md)); "w-engine" 4.1 is where the engine's price exponent
> z < 1 is read off; both are in my private archive.

## The first skeptic: the logic, step by step

*Verdict: not refuted.*

Verdict: the main results hold. One side claim marked [proved] is false as written, and the "no growth" wording claims more than was proved. Both need correcting before anything is quoted.

**Confirmed, line by line**
- **Fact A.** The free letters give the whole stabiliser P: the address maps give GL(m,2), and the phase transvections along (0|z) give every symmetric matrix. Transitivity together with P inside the image gives the image Sp(2m,2), which is distance-transitive. These are exactly the hypotheses of the Averaging Lemma (w-bounds section 3), which allows any family, including ones whose dimension varies.
- **Fact B.** The letter determinants are powers of i, so the only scalars in Γ are ±1 and ±i, and these are free letters (i = phaseMat 0 1). The Paulis are free too, so B is the full preimage of P, and every certificate gives a family whose δ depends only on distance.
- **Lemma 2.** I re-derived all four cases (i)-(iv) and the totals sum to 4^m − 1.
- **Claims 1 and 2.** Both hold for any base family; distance-determined is not needed.
- **Theorem 4 (the coupling).**
  - C_r ≤ 2^(−j−1)·P_k and A_r ≥ (1 − 2^(−j))·P_k both check by algebra.
  - The dominating chain is monotone (p_up(j) + p_down(j+1) = 1), and the step-by-step stochastic domination is valid.
  - A chain started at its bottom state increases stochastically and stays below π*.
  - R_(M−j) is nonincreasing in j, and S_M/S_m ≥ 1.
  - c* is a certified lower bound because the normaliser is bounded above.
- **Theorem 5.** The telescoping gives S_∞ = 2(2^m0 − 1)/(2^(m0+1) + 1).
- **Theorem 6.** The pullback cut argument and the point-cut count check.
- **Proposition 7(a) and (b).** Direct sums and tensors of distance-determined families: the concavity proof is correct.
- **Proposition 8.** It works at engine prices because the engine has z < 1 (w-engine 4.1), so (j/m)^z ≥ j/m.
- **Lemma 1(a)-(d).** All hold, including the dichotomy and both unit embeddings. The lifted letters are the same kind of letter, the ranks stay below the block size, and fresh scratch arrays are used.
- **Arithmetic.** By hand, T_5 = (558·86592 + 465·68320)/(975·27456) = 834241/278850. The excluded units and the savings 0.319 and 0.392 check.
- **Independent check of the base values.** A rank-j block acts as G ⊗ I, so f splits over the untouched bits into small local ranks of at most 280 × 280. My own code (not b56-compute's), mod two primes, gives f_1 = 27456, f_2 ≥ 48576 and f_3 ≥ 68320, exactly the values used. A sanity check reproduces the block-3 values 80 and 144. Files: `../programs/engine-bounds-blocks-5-and-6/skeptic-checks/fcheck.py` and `../data/engine-bounds-blocks-5-and-6/skeptic-checks/out_fcheck.txt`.

**Problems and required corrections**
1. **Proposition 7(c) is false as stated.** The claim is "General tensor products: at most 2 × the best". The proof silently assumes each family has the same dimension at every frame. The Averaging Lemma allows families whose dimension varies, and then the claim fails.
   - Counterexample at m = 2, with cut families F_U = span{e} for U in C and 0 otherwise. Then F_C ⊗ F_D = F_(C∩D).
   - Two 9-element cuts with T = 39/56 meet in a point cut with T = 3/2. That is 28/13 = 2.154 times the larger factor, more than the claimed 2.
   - Files: `../programs/engine-bounds-blocks-5-and-6/skeptic-checks/tensor_cx.py` and `../data/engine-bounds-blocks-5-and-6/skeptic-checks/out_tensor_cx.txt`.
   - Fix: add the hypothesis "families of constant dimension". The proof is then correct, and it still covers every certificate family and every projection, pullback or sum of one.
   - The counterexample has T < 1. I did not search for one where both factors have T ≥ 1.
2. **"No growth by induction" claims more than was proved.**
   - "None of the inductions can grow" and "no theorem 'block-m bound ⇒ larger block-(m+1) bound' exists in this toolbox" are contradicted by the proof's own `out_tail.txt`:
     - projection takes the perfect block-2 family from 2 to 2.0588 at block 3, rising to 2.2;
     - the block-3 hook goes from 2.636 at m = 4 to 2.653 at m = 5.
   - A pullback may gain up to a factor κ_m, close to 2.
   - What was actually proved: each single operation stays within a fixed factor of its base. So a block-m0 computation fed through projections, pullbacks or sums certifies at most about 2·m0 at every block size.
   - Not covered: combinations such as k-fold tensors of pullbacks that were not symmetrised first. The general tensor bound weakens to k × the best factor.
   - The wording should be restated to say exactly this.
3. **A missing sentence in the Theorem 4 corollary.** The corollary checks that the profile built from the lower-bound f values is nondecreasing, not the true family's profile. It is still valid, because T_m and the R_r increase with each f_s while ā_1 depends only on the exact f_1. The proof file should say so.
4. **κ_m < 2 for every m rests on computation.** It is exact only up to m = 12 (40 in the compute skeptic's check) and argued beyond. The statement does hedge this, but "(3) [proved]" covers the bullet that relies on it.

**What survives unchanged**
- Items (1), (2), (4), (5) and (6) of the statement.
- In (3), the bullets for projection, pullback, direct sums and certificate tensor products.
- The new bounds: 2.99172 at block 5, 2.97290 at block 6, and at least 2.9447 at every m ≥ 4, in the engine game.

Nothing was posted. Every run went through a memory guard and took under 10 s.

## The second skeptic: every number again, with its own code

*Verdict: not refuted.*

I could not refute the result. Every key number came out the same in my own independent code, and each proof step I checked holds. My scripts are in `../programs/engine-bounds-blocks-5-and-6/skeptic-checks/` and my outputs in `../data/engine-bounds-blocks-5-and-6/skeptic-checks/`. All runs went through a memory guard, each under 70 s. No Colab or Lean was used, because nothing in this result is Lean.

**1. Lemma 2, the closed form for N_r(s)** (sk_nrs.py → out_sk_nrs.txt)
- I wrote a fresh brute force, separate from nrs.py. For each case it builds a pair of Lagrangians from a random symplectic basis and counts every nonzero p.
- It agrees with the closed form in all 28 cases (m = 1..7, every r), with 0 mismatches.
- A corrupted closed form (one count moved from s = r−1 to s = r−2) was rejected in 6 of 6 cases.
- I also re-derived the proof by hand: cases (i)–(iv), and the count from the perfect pairing L/H × M/H. It holds.

**2. Claim 1 of Theorem 3, tested on explicitly built families** (sk_claim1.py → out_sk_claim1.txt)
- I wrote my own symplectic identifications ι_p: p^⊥/p → F_2^4. With them I projected two block-2 families that have no symmetry to block 3:
  - a random coordinate family;
  - a random subspace family, with ranks taken mod 10007.
- I then computed the exact ā^(3)_r over all 135² ordered pairs. They equal Σ_s N_r(s)·ā^(2)_s exactly for both families: 697/3, 16163/60, 8981/30 and 493/5, 551/5, 609/5.
- So the recursion does not need the base family to be distance-determined.

**3. The headline numbers** (sk_tail.py → out_sk_tail.txt; my own recursion)
- From the base (27456, 48576, 68320, 86592) I get T_5 = 834241/278850, T_6 = 1188222/399685 and T_7 = 2361819/799370, exactly as claimed.
- The exact T_m decrease monotonically up to m = 200, with minimum 2.944709234121 ≥ 2.9447. The margin is about 9e-6, so the bound holds but only barely.
- The profiles at blocks 5, 6, 7, 30 and 60 are nondecreasing, satisfy R_r ≤ r and are subadditive. The block-200 profile is nondecreasing.
- The domination inequalities have 0 violations for k ≤ 200. A corrupted dominating chain (up-probability 2^−(j+2)) gives 741 violations.
- c*_60 = 2.944709234121 and c*_30 = 2.944709230.
- S_∞ = 2/3, 14/17, 10/11, 0.95385 for m0 = 2..5, and the ceilings are 3, 3.6429, 4.4, 5.2419.
- The excluded units check out. At blocks 5 and 6: 3 in 8, 4 in 11 and 5 in 14. "W arrays in 3W−1 moves" is excluded up to W = 120 at block 5 and W = 36 at block 6.

**4. The base data**
- In `../data/block4-hook/`, out_rank2_blocks.txt lines 3 and 8 and out_m4_hook8_G.txt line 6 give f ≥ 48576, 68320, 86592. They are certified mod p, with exact dim N = 102960.
- f_1 = 27456 is the established exact cap.
- The rounding goes the safe way: ā^(m)_1 depends only on f_1 (exact), and every N is ≥ 0. So using lower bounds for f_2..f_4 gives valid lower bounds on T_m.

**5. Re-running the author's scripts**
- All six of the author's scripts (nrs, induct 60, tail, cutceil, checks, prices) reproduce their out_*.txt byte for byte.

**6. Proof steps checked by hand**
- **Theorem 4.** The algebra for C_r and A_r is right. The dominating chain J* is monotone (p_up(j) + p_down(j+1) = 1). The step-by-step domination is right. Since J* starts at 0 it increases stochastically, so J_M ≤_st π*. The conditioning on J_M uses the Markov property correctly, and S_M/S_m ≥ 1. The monotonicity it needs is that of the computed profile, which is a lower bound term by term.
- **Theorem 5.** The telescoping products are right.
- **Theorem 6.**
  - The intersection array matches the standard dual polar graph, and the eigenvalues are θ_j = 2(2^(m−j)−1) − (2^j−1).
  - The point-cut count is right.
  - The pullback argument needs f'_1 ≤ f'_s, which the proof states.
  - I extended the exact κ_m = 2 − 2^(1−m) check, at θ_1, to every m ≤ 40 (sk_kappa.py). It is still only argued for all m, as the result labels it.
- **Proposition 7(b).** The concavity argument holds, using (A⊗B) ∩ (A'⊗B') = (A∩A') ⊗ (B∩B').
- **Fact A.** dir z conjugates Z^b to iZ^bX^z, so it is a transvection along (z|0). The free letters generate P = GL(m) ⋉ Sym(m).
- **Fact B.** The scalars in Γ are ±1 and ±i, and Z^z = (phaseMat z 0)². So B is the full preimage of P.
- **Lemma 1.** In (a)–(c) the residue isometry and the Sp-invariance argument for uniformity hold. In (d), dir(z,0) = dir z ⊗ I, the kernel factors as claimed, and the second word can use fresh scratch.
- **Proposition 8.** It relies on correctness being formal, which w-bounds section 1 derives with a Vandermonde argument. Given that, the term-history argument holds.

**Caveats, none of which is a refutation**
- The all-m constant 2.9447 is only about 9e-6 below the limit 2.9447092341.
- κ_m for m > 12 is labelled [argued]. My exact check now covers m ≤ 40, but there is still no general proof that θ_1 maximises.
- The engine-price caps are, as the result says, indicative only.
- All bounds hold in the engine game only, and depend on Theorem R.
