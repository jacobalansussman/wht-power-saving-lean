<!-- header: existence-routes -->
# The skeptic's record: star lemma and slice

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** This is an AI skeptic's record of its check, as it wrote it. Its verdict is the status this revision uses for the track it checked, unless a later verdict is named in the section README.

---

# Skeptic's check of e3-slice6: 13 in 63 is impossible in the commutative slice

exist3, skeptic v-e3-slice6, 15 October 2026, 05:54 to 06:20 EDT. Nothing was posted. Predictions were written
before any reading or computing (`STATUS.txt`) and scored at the end.

## Verdict

**Not refuted. The calibration holds.**

- The headline exclusion stands: the uniform family U_6 (one move in each of the 63 directions) cannot serve 13
  arrays in the commutative slice. With 14 to 16 already excluded and 12 below the bar, U_6 has no member above the
  bar.
- Theorems A, C, D and E all survive. I checked each proof line by line and rebuilt every computation the claims rest
  on with my own code. My numbers agree with the track's everywhere.
- I found one factual slip, and one claim that was stated more narrowly than it could have been. Neither changes a
  result:
  1. **The slip.** The note says the uncorrected 2-in-7 program is "wrong at the other 448 characters". It is wrong
     at 256 and right at 256: it is right exactly when the linear part λ has even weight. The track's script reported
     only a maximum error, not a count, so the figure 448 was an inference. Theorem A is unaffected.
  2. **The narrow claim.** The note proves the "3/4 floor" of the local method only at even n and at n = 5, yet it is
     quoted as a floor in general. It does hold at every n. At odd n, take the bent matrix on n - 1 coordinates
     tensored with the 2 x 2 swap matrix; this gives L(n, 2^(n-2)). I checked it exactly at n = 3 and 5, and the
     general proof is three lines (section 4).

## 1. The exclusion (Theorem C and its use)

The chain is: Theorem 10 (b) and (c), then Lemma 12', then the star lemma, then Theorem C.

- **Theorem 10.** I re-derived it: two Schur complements of one bordered matrix. The earlier skeptic
  (v-ex2-families) had also checked it. On the genuine U_3 unit:
  - E on every even support is nonsingular (smallest singular value 0.36);
  - on odd supports, the rank beyond the slack is 0 (largest singular value there 2.6e-16).
- **Lemma 12'.** I checked the hand proof again. I also rebuilt all 2^21 quadratic forms at block 6 myself:
  - F = {z1 = z2 = 1} is an even support of 16 directions;
  - each of its 30 affine hyperplanes lies in two or three odd supports of 16 directions;
  - so at W = 13 every hyperplane block of E has rank at most 3.
- **Star lemma.** Every step of the proof is right:
  - the Schur complement T is nonsingular;
  - rank on Ψ ∪ Q_j equals k plus the rank of a block of T, because R and C lie inside;
  - the row count and the column count each use "subspaces sharing a common subspace";
  - both are closed by the inequality rank [[A,B],[C,D]] ≥ rank [A B] + rank [A; C] − rank A, which I re-proved.

  I tested it exactly over GF(2), GF(3), GF(10007) and GF(1000003): 1,432 instances, on random, low-rank-structured
  and bent matrices, with both flat partitions and arbitrary ones. There were 0 violations, and many cases were tight
  on the bent matrix. The lemma holds over any field, so the complex case is covered.
- **Theorem C.**
  - Planes: with g = 3, a nonsingular 2 x 2 block inside a plane would need 4 ≤ 3 − a, which is false. So the rank
    on every plane is at most 1.
  - Lines: with g = 7, a nonzero entry would need 2 ≤ −5a, which is false. So E = 0.

  I also checked a second route that avoids the line step. Any 3 points lie in a plane, so every principal block of
  order 3 or less has rank at most 1, and Lemma 11 with s = 1 makes E singular. The geometric fact the lemma uses
  (the (d+1)-flats through a d-flat split the rest of the space) I re-checked by brute force for n = 2 to 5.

## 2. The same argument at every block (Theorem D)

- My own recursion reproduces the track's thresholds exactly for m = 4 to 14: W ≥ 13, 26, 52, 105, 209, 419, 838,
  1677, 3354.
- It also reproduces the open ranges above the bar, for example 23 to 25 in 127 at block 7.
- Lemma 12' holds for every m ≥ 4 (I checked the hand proof), so the recursion's starting bound is sound at every
  block.
- The rule "stop when s_d ≤ d − 1" is a correct use of Lemma 11: any s_d + 2 points lie in a d-flat.

## 3. Near-uniform multisets at block 6

My own recount gives every number in the track's table:

| Arrays W | Extra moves | Excluded |
|---|---|---|
| 13 | 1 | 62 of 63 |
| 13 | 2 | 1,922 of 2,016 |
| 14 | 2 | 2,015 of 2,016 |
| 14 | 3 | 43,586 of 43,680 |
| 15 | 3 | 43,679 of 43,680 |
| 16 | 1 to 3 | all |

- Every case left open involves the direction 111111.
- I rewrote the mixed-bound refinement myself. It closes the two extra moves 111111 at W = 14, through a plane whose
  three hyperplanes all have slack 2.
- It closes none of the 95 open W = 13 cases, as the track found.
- Calibration: it excludes U_6 at W = 13 and leaves W = 12 open, as it should.

## 4. The floor of the local method

I checked these with exact rational ranks:
- the bent matrix at n = 4: rank 16, every hyperplane rank 4;
- the annealed circulant at n = 5 from `out_circ5.txt`: rank 32, hyperplane ranks at most 8. The result is the same
  under either bit order of the truth table;
- the Singer permutation at n = 3: rank 8, hyperplane ranks 2;
- the Singer permutation of F_16 at n = 4: rank 16, hyperplane ranks exactly 4. So L(4,4) holds over every field.

The bent case for general even n follows from the Fourier argument. A hyperplane block of E has rank equal to the
number of cosets {u, u+l} on which the dual sign agrees, and that number is 2^(n-2) because bent autocorrelations
vanish.

**Odd n (my addition).** Let E = B ⊗ S, where B is the bent matrix on n − 1 coordinates and S = [[0,1],[1,0]].
- Hyperplanes containing the last coordinate's line have rank 2^(n-3) · 2.
- The two hyperplanes "last coordinate fixed" have rank 0.
- The mixed ones have rank 2 times the rank of an off-coset block of B, which is 2^(n-3).
- So every hyperplane rank is at most 2^(n-2). Checked exactly at n = 3 and 5.

So the local route cannot pass 3/4 of full sharing at any block, not only at even blocks.

The remark "circulants cannot reach s = 7 at n = 5" holds for ±1 spectra by the autocorrelation sum, as written. It
also holds for any nonzero spectrum: ordered anti-equal pairs can number at most 2 · 16 · 16 = 512 in all, but 31 · 18 = 558 are needed.

## 5. Theorem E (block 5: M ≥ 4.65 W)

- My own program rebuilt all 2^15 forms at block 5 and checked the meaning of all 116 certificate rows:
  - each plane bound is "a plane inside an odd support";
  - each star row is "a line, plus the three planes through it, whose union is an even 3-flat support".
- Adding the rows in exact fractions gives Σ c ≥ 93/20, using only c ≥ 0 and w ≥ 0. VERIFIED.
- Validity calibration: the plain 1 in 5 satisfies all 6,300 star rows at block 5, with minimum slack 0, and the
  plain 1 in 4 satisfies all 196 point-and-lines rows at block 4.
- The uniform point c = 3/20 should be exactly feasible: the star rows are tight there, with each plane bound
  8c − 1 = 1/5 [argued]. This agrees with the track's floating-point optimum. The theorem needs only the lower
  bound, which is exact.
- Projection (Lemma 3) carries the bound to every block from 5 on.

## 6. Theorem A (the order-2 system is exact)

- The algebra is right:
  - the corrected move t_z u_z takes the value (−1)^Q(z) at every character;
  - v takes the value i^|λ|;
  - the target K takes the value i^|λ| (−1)^tr Q.
- The bordered-rank strengthening in section 3 of the track's note, and its block-2 example (Γ_21 = 2/k, Theorem 10
  failing only at 1/k), are right.
- My own register-level build of the skeptic's 2-in-7 order-2 solution:
  - with the corrections, it is right at all 512 characters (error 1.1e-15), and on fresh random complex arrays with
    f = 1 to 6 digit columns (relative error ≤ 1.6e-15);
  - without them, it is wrong at 256 characters, not 448 (the slip above).
- The exact plain programs at W = 1 and 2 are right at all 512 characters.

## 7. What I did not do

- No exhaustive search for L(4,3). My GF(2) annealer found L(3,2) but failed to find the known L(4,4), so it was not
  calibrated at n = 4. I did not use it as evidence and did not run it on L(4,3). The proof carries Theorem C.
- I did not re-run Theorem 17. Neither did the track; it said so.
- Leftover jobs: I could reach only two of the six Colabs (ssh4 and ssh5), and neither had any e3-slice6 job left. The
  other four timed out.
- No msolve inputs exist in e3-slice6. Its Singular runs (the superseded Singer ansatz) are unaffected by the x*x bug.

## 8. Mistakes of mine

My first geometry script ran over 2 minutes on the Mac, because it enumerated bases at n = 6. I killed it at once and
rewrote it for n ≤ 5. The mixed-bound recount took 114 s on the Mac, which is close to the limit.

## Sources

All under `the archive's checks/wht43/exist3/v-e3-slice6/`. Each program has an output file of the same name with
the prefix `out_`.

- **Predictions and log:** `STATUS.txt`, `RUNLOG.txt`.
- **Theorem C and its chain:**
  - `v_forms6.py`: all block-6 forms; the hyperplanes of F and their odd covers;
  - `v_star.py`: exact tests of the star lemma;
  - `v_geometry.py`: the flat-partition fact;
  - `v_thm10_u3.py`: Theorem 10 on the U_3 unit.
- **Theorem D:** `v_recursion.py`.
- **Near-uniform multisets:** `v_nearU6.py` and `v_mixed.py`.
- **The floor:** `v_exact_local.py` and `v_kron_odd.py`.
- **Theorem E:**
  - `v_cert5.py`, run on the track's `starlp2/remote-out/star_cert_b5.json` and `star_vars_b5.json`;
  - `v_plain4.py`, the block-4 validity calibration.
- **Theorem A:**
  - `v_thmA.py`, with outputs `out_v_thmA.txt` and `out_v_thmA_plain.txt`;
  - `v_thmA_lit.py`: the count of characters by λ for the uncorrected program;
  - `v_thmA_raw.py`: fresh random arrays.
- **Uncalibrated, unused:** `anneal/anneal_gf2.py`.
- **Read:**
  - `exist3/BRIEF.md`;
  - `exist3/e3-slice6/NOTE.md` (all of it);
  - `exist2/ex2-families/PROOF.md`: Lemma 1, Lemmas 2 and 3, Theorem 10, Lemma 11, Theorem 12, Lemma 12', and
    Corollaries 13 and 15;
  - `exist2/v-ex2-families/NOTE.md`, the verdict lines.
