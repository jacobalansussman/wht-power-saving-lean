> **Working note, kept as written.** An AI agent wrote this during the logic round of 14-15 October 2026 (the skeptic of the ml-fooling track). Private paths and machine names
> were removed and nothing else was changed. No person has checked it. Where it disagrees with this section's
> [README](../README.md), the README wins, because it follows the later skeptic verdicts.

# v-ml-fooling: skeptic's verdict on meta3/ml-fooling

15 October 2026, 00:07 to about 00:40 EDT by the Mac clock. I read BRIEF.md, ml-fooling's NOTE.md, STATUS.txt and
RUNLOG.txt, its checker, its certificates and the code and logs behind the claims, plus ml-logic's NOTE, v-ml-logic's
VERDICT and ex2-modp's NOTE. I read four sources first-hand: Kollár 1988 (Theorems 1.1, 1.5, Remark 1.6, Corollary 1.9),
Laurent and Mourrain arXiv:0812.2563 (Section 1.3, Theorems 1.4 and 1.5), and Sombra alg-geom/9610006 (main theorem,
and page 6, which quotes Heintz's Bézout inequality as Heintz 1983, Theorem 1). Nothing was posted. Everything ran on
the Mac in pure Python, in seconds.

**Verdict: REFUTED in part.** The mathematics holds: the two-certificate theorem, the r bound, the six certificates,
the counts at infinity, the sweep counts and every degree I recomputed. What does not hold is the reading of what this
means. Three of the headline sentences about what logic and algebra allow are false as stated. Two of them are
contradicted by the project's own files, and one by a computation here. Separately, the sweep's "refuted" verdicts are
modulo one prime only, yet the summary calls them refutations and calls the tool an "exact refuter".

## Upheld (checked)

- **Theorem 2.5** ((R) or (M), exactly one exists). The proof is correct. The rational cofactors in (R) follow from
  linear-algebra descent (Kollár, Remark 1.6(i), read first-hand), or directly from step 2. The converse, from a trace
  functional L(p) = tr p(M)/r, is also right in characteristic 0.
- **The bound r ≤ d_1⋯d_k.** The argument is correct: a generic rational slice, then the Bézout inequality
  deg(V ∩ W) ≤ deg V · deg W, then a Galois orbit. The inequality is now confirmed as Heintz 1983, Theorem 1, quoted
  on page 6 of Sombra's paper, which I read first-hand. It is no longer from memory.
- **Literature quotes.** They match the sources: Kollár's Theorem 1.5 (with the condition d_i ≠ 2), Brownawell as
  quoted in Kollár's Theorem 1.1, Sombra's min{n,s}²(d+3n)δ and δ ≤ (d+1)^(min{s,n}−1), and Laurent and Mourrain's
  Theorems 1.4 and 1.5 (any field K, C connected to 1). The corollary in 2.3 with the relaxed condition
  deg f_i ≤ s + 1 is proved correctly: deg(f_i·q_C) ≤ 2s, so L̃ agrees with L there.
- **The six certificates.** I rechecked them with code of my own, which imports nothing from ml-fooling
  (`indep_check.py`, `out_indep_check.txt`):
  - All six are exact: the matrices commute, every slice equation vanishes, and every **original** equation vanishes
    at the matrices built from the substitution.
  - The originals are the right systems. I rebuilt the Karatsuba and W-state rank-3 equations from first principles.
    For the two nag-rerun systems I lifted the msolve exports myself (x*x read as x², symmetric lift modulo the first
    prime, checked against the second). All four match the certificates exactly.
  - Eigen-points give a residual of at most 2.2·10^−12 on the originals. The rank-3 decompositions rebuilt from them
    reproduce both tensors. The Fp2b points are not real (coordinates 0.5 ± 0.5i).
  - `check_cert.py` reports INVALID for a 10^−9 change to one entry, and for a changed constant in the original
    system.
- **Counts at infinity** (`infinity_verify.py`). The Bézout numbers are 1,728 and 144. Substituting w1 = −w0 turns the
  top-degree forms into (v0−v1)w0, (u0−u1)w0 and (u0v0−u1v1)w0. So the locus at infinity is exactly
  {w0 = w1 = 0} ∪ {u0 = u1, v0 = v1, w1 = −w0}. Its projective dimensions are 5 and 6 for Karatsuba, and 3 and 3 for
  the W-state, as the note says.
- **Degrees, with my own NS code** (`ns_mine.py`; primes 10^9+7 and 998244353, and exact over Q):
  - All 58 saved skeletons: 15 at degree 2 and 43 at degree 3. That agrees with ml-fooling, skeleton by skeleton, and
    holds exactly over Q for all 58, including the 37- and 64-unknown ones.
  - Masser–Philippon chain: 4, 7, 12 (exact over Q for n = 2, 3).
  - W-state rank 2 and Karatsuba rank 2: NS degree exactly 9 (none at 8).
  - The worst-case degree table: log10 D_eff = 17.7, 38.5, 3.8, 5.7, 617 and 8,794. The sizes 610, 2,377, 37, 69,
    1.0·10^6 and 1.6·10^8 agree.
- **Sweep counts.** 45 + 28 + 3 = 76 refuted and 57 still over the cap, matching the logs. The cores add up:
  418 + 62 and 279 + 71 make w3-asym's 830. The random calibration was 20 planted, 20 square and 20 impossible: 41
  certified and 19 refused, as claimed.
- **Fixed-shape logic.** It is right. For each shape the answer is decidable, both answers have finite witnesses, and
  neither can be independent of weak arithmetic. "Some unit saves more than 0.1" is Σ1, consistent with ml-logic's
  Theorem A, which v-ml-logic upheld.
- **Timeline.** File birth times fit the two sessions. Sombra's PDF was downloaded at 20:12. `commcert.py` was
  created at 23:06:48, after the Q predictions, and `out_commcert.txt` at 23:08. STATUS.txt and RUNLOG.txt were
  rewritten at 20:57 and 21:00, so "never edited" cannot be checked from metadata, but nothing contradicts it.

## Refuted or overclaimed

1. **"A proof that an existence proof exists is never cheaper than the existence proof itself."** No theorem in the
   note supports this, and it is false as stated.
   - *The note's own section 7, item 2.* Suppose a structural theorem bounds the PC refutation degree of a family by
     c. Then a degree-c pseudo-solution proves existence. It is not flat, it is not the unit, and it is far smaller.
   - *Computed here* (`hensel_demo.py`). For a square slice, one prime suffices. Take a point x mod p with
     f(x) ≡ 0 and det J(x) ≢ 0. Hensel's lemma lifts it to a zero in Z_p, so the ideal over Q is proper and a complex
     zero exists. For the Karatsuba slice (12 × 12):
     - p = 10007 works, with det J = 1250, and the lift to p^16 checks.
     - The certificate is 12 residues of 14 bits. It does not say what field the unit's numbers live in.
     - The W-state generic slice is also square: 8 × 8, dimension 0, vdim 48 = r, so all its points are simple. There
       a one-prime certificate of 8 residues should replace (M), which needed 48 × 48 matrices with 796-bit entries,
       54 primes and about 20 minutes to check [argued, not computed].
   - *Proof theory.* Speed-up theorems (Gödel 1936; Parikh 1973 [memory]) say that proofs of provability can be much
     shorter than the proofs themselves. So "never" is not something logic gives us for free.
2. **"The only non-constructive alternative, counting"** (claim, "meaning"), and "the one classical way to prove
   existence without exhibiting a solution is to count solutions" (NOTE 3.1). This is false. Other routes:
   - the note's own worst-case pseudo-solution (section 2.1);
   - ex2-modp's Proposition 5, which says that zeros modulo enough primes prove a zero "without saying where the
     solution is" (ml-fooling cites this as its Q3);
   - Hensel lifting, as in item 1;
   - certified numerics (Smale's α-test, interval Newton) for square slices [memory];
   - refined Bézout counting (item 3).

   Most of these fail on overdetermined unit systems, because surplus equations need exact vanishing. That is a
   different obstruction from "solutions at infinity", and the claim does not name it.
3. **"Counting needs no solutions at infinity, which is itself an impossibility statement … so counting proves only
   border units."** This is false in principle.
   - Refined Bézout counting does not need the boundary to be empty, only to be measured. This is excess intersection
     (Fulton–MacPherson; Fulton, *Intersection Theory*, chapters 12 and 13 [memory]).
   - The total 1,728 is the sum of positive contributions from the distinguished varieties.
   - For Karatsuba, an upper bound below 1,728 on the contribution of the locus at infinity would prove that a finite
     solution exists. That is a Segre-class computation, not a refutation. Here the true split is 4 + 1,724.
   - The infinity numbers are right (see Upheld); the inference "blocked, border units only" is not. What is shown is
     that *plain* Bézout is blocked.
4. **NOTE 2.4: "modulo p … is not a proof over C"**, stated generally. This holds for overdetermined or singular cases,
   which include ex2-modp's decoy. It fails for square slices with a nonsingular root (item 1). Three of the
   calibration certificates came from such slices: K3_coord, W3_coord and the r = 26 and r = 48 generic slices.
5. **The sweep's "refuted" verdicts are modulo one prime only.**
   - `sweep_ns*.py` calls `ns_degree` and `ns_partial` with the single default prime 2147483629. There is no second
     prime and no check over Q.
   - A refutation modulo p is not a proof over C. For example, {p·x − 1} is refuted mod p, yet it has the solution
     x = 1/p.
   - NOTE section 5 labels the results "[computed, mod p]". But section 0 item 10, the claim's "result" 3 and
     section 9.1 ("Singular-free **exact** refuter") drop that label, and the script's docstring says "an exact
     identity … mod p".
   - The 76 are very likely right, but they have weaker standing than w3-asym's Singular-over-Q verdicts. (v-ml-logic
     made the same kind of point about w3-asym's "certificates".)
   - "None fooled degree 3" says nothing about the 57 skeletons over the cap: they were never tested, because the
     script returns CAP before any degree test.
6. **Calibration of the exact pipeline on our own units.** The brief requires a known unit of our problem to be
   certified: the plain schedule, the block-2 stand-in 2 in 3, or the 8-move control word.
   - The commuting-matrix pipeline never certified any of them. Its calibration is on stand-ins: Karatsuba, the
     W-state, and nag-rerun's control and relaxation hit.
   - The control word was dismissed as "not a nonlinear test". But an r = 1 certificate of it would have been a cheap
     calibration on a real unit.
   - The r = 26, r = 48 and spu r = 5 certificates exist only as log lines. They were not saved and cannot be
     rechecked. The r = 48 run crashed on a missing folder.
7. **Minor points.**
   - The product form of the (R) bound for quadratic systems (our skeletons) is not Kollár's: his Theorem 1.5 excludes
     d_i = 2. It needs Jelonek or D'Andrea–Krick–Sombra.
   - The note paraphrases Kollár's Corollary 1.9 as "meet in one point … outside the hyperplane at infinity". In
     Kollár, the point is on that hyperplane; it is the tangent space that is not contained in it.
   - (M) over Q assumes rational coefficients. Games with Z4 phases have coefficients in Q(i), which needs i adjoined
     as an unknown with z² + 1 = 0 (a trivial fix).
   - "Coordinate slicing refutes it" should read "misses it": a sliced system with no zero refutes only the slice.
   - The Sources list says "four exact certificates"; there are six.

## What the claim should say instead

- For one fixed shape, exactly one of (R) and (M) exists. Both have size bounds, so the fixed-shape question is decided
  in weak arithmetic.
- The flat (standard) pseudo-solutions are the finite models.
- Existence proofs that do not exhibit (M) do exist:
  - worst-case-degree pseudo-solutions;
  - many-primes counts;
  - Hensel lifting and certified numerics on square slices;
  - refined Bézout counting;
  - a structural bound on refutation degree combined with a low-degree pseudo-solution.
- On our systems, each of these is blocked or untested, for named reasons: size, surplus equations, an unmeasured
  boundary, or no structural theorem.

## Files (all in this folder)

- `indep_check.py` → `out_indep_check.txt` (independent check of all six certificates)
- `infinity_verify.py` → `out_infinity_verify.txt` (Bézout numbers and the locus at infinity)
- `hensel_demo.py` → `out_hensel_demo.txt` (one-prime existence proof for the square Karatsuba slice)
- `ns_mine.py` → `out_ns_skeletons.txt`, `out_ns_mp.txt`, `out_ns_W2.txt`, `out_ns_K2.txt` (my own NS degrees)
- `out_deff_check.txt` (the worst-case degree table)
- `RUNLOG.txt`

Downloads (Kollár, Laurent–Mourrain, Sombra) were read with pdftotext only, in a temporary folder.
