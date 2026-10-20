<!-- header: existence-routes -->
# The commutative slice at block 6: Theorems A, C, D, E and the star lemma

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** Argued by hand, with exact checks by computer (re-run 19 October 2026); one AI skeptic checked every proof line by line and rebuilt every computation. Corrections: the uncorrected 2 in 7 is wrong at 256 characters, not 448; the 3/4 floor holds at every block (the skeptic's odd-block construction). A second skeptic was asked for and has not run.

---

# The commutative slice at block 6: 13 in 63 is impossible, and the order-2 test is exact

exist3, track e3-slice6, 15 October 2026, 04:32 to 05:50 EDT. Nothing was posted. No skeptic has read this.
Labels: [proved] (proof in section 7 of this note), [computed] (program and output named in the sources),
[argued], [guess]. Predictions were written before computing (`STATUS.txt`) and are scored in section 6.

## In short

1. **13 in 63 is impossible in the commutative slice** [proved].
   - This was the brief's first target: the uniform family U_6, one move in each of the 63 directions, shared by 13
     arrays (saving 0.119).
   - So the uniform family has no member above the bar at block 6:
     - 14, 15 and 16 in 63 were already excluded (ex2-families);
     - 12 in 63 and below save less than 0.1.
   - The proof is a non-linear use of the matrix E of Theorem 10, as the brief asked. It looks only inside one even
     support, the 16 directions with z1 = z2 = 1, which form a copy of the 4-dimensional affine space over F_2.
   - The new ingredient is a "star lemma" about nonsingular matrices on such a space. Applied twice, it shows that
     no nonsingular 16 x 16 matrix has rank at most 3 on all 30 of its hyperplane blocks:
     - rank at most 3 on hyperplanes forces rank at most 1 on planes;
     - rank at most 1 on planes forces every entry to vanish.

     Theorem 10 forces exactly that pattern on a 13-in-63 unit.
   - The same criterion also excludes almost every "U_6 plus one to three extra moves" above the bar [computed]: for
     example 62 of the 63 possible 13 in 64.
2. **The same argument at every block** [proved; numbers computed].
   - U_m is now impossible in the slice once W reaches about 0.81 to 0.82 of full sharing (13/16 at block 6). The
     old threshold was 0.9.
   - Block 7: 26, 27 and 28 in 127 are newly excluded (savings 0.185 to 0.222). 23, 24 and 25 in 127 stay open
     (0.122 to 0.165). This local argument cannot reach 23 or 24 at all (see the floor below).
   - Block 8: everything from 52 in 255 up is excluded. 40 to 51 stay open.
   - The method has a floor: this local argument can never go below 3/4 of full sharing. At every even block, a
     matrix built from a bent function meets the pattern exactly, and at block 7 a similar circulant does. So it
     cannot touch U_m with a constant fraction c of sharing below 3/4, and that family would still save close to 1.
3. **The order-2 test is exact** [proved; checked on U_3].
   - A slice program exists exactly when the order-2 system of Lemma 1 has a solution. That system looks only at the
     2^21 characters where every shift acts trivially. The full group has 2^27 characters at block 6, and the other
     ones add nothing.
   - The reason:
     - free shifts and copies make every combination of shifts free;
     - the group algebra then splits into 2^m blocks, one for each linear part λ of a character;
     - every block is a copy of the order-2 block, through t_z -> (-i)^[λ·z] t_z.
   - Earlier slice existence statements were solutions of the order-2 system: the 2 in 7 at block 3, and the
     solver's "found" verdicts. ex2-families noted that turning them into programs needed g-play's reduction, which
     had no written proof. Theorem A supplies it.
   - Checked on the skeptic's 2-in-7 solution:
     - realised literally, it is right at the 64 order-2 characters and wrong at the other 448;
     - with the free corrections the proof prescribes, it is right at all 512 (error 1.4e-15).
4. **The converse question has a clean answer** [proved].
   - Theorem 10's conditions on E do not determine a program. For a given program they say only that, on an even form,
     its matrix X_Q has no eigenvalue 1. A program needs X_Q = 0.
   - The exact strengthening is the order-2 system itself. So searching for "an E at block 6" is the same problem as
     searching for a program, and that is far out of reach there: about 3,600 unknowns and 3.5 x 10^8 equations.
   - Numerical search for E is unreliable, and I calibrated that on cases with known exact solutions. One early run
     (before I fixed the scale) found a nonsingular solution of the 8-point case. The other 19 runs there, and all 3
     at s = 5 on 16 points, drifted towards singular matrices: special position again.
   - So there was no block-6 candidate to send through the finite-field pipeline.
5. **The slice bound at block 5 rises from 4.306 to 4.65 moves per array** [proved; exact certificate, re-checked by
   an independent program] (Theorem E).
   - Every commutative-slice unit at block 5 has M >= (93/20) W.
   - The proof is the star lemma used as a linear inequality: one row for each line inside an even 3-flat, together
     with the three planes through that line. These rows replace the ten-piece cuts of Theorem 18 and beat them.
   - By projection, every slice unit at every block from 5 on needs at least 4.65 moves per array. The best slice
     saving at block 5 is now at most 0.045 (it was 0.093).
6. **The general slice question at block 6 stays open, but narrower.**
   - Moves per array needed now: at least 4.65 (was 4.306). The bar allows 5.02.
   - This round closes the natural family at block 6, not the whole slice.
   - LP rows of this kind cannot reach the bar there [argued]. The uniform point stays feasible at about 4.73, so a
     closing argument needs the integer rounding that Theorem C uses, or something new.

## Words used here

- **Block m, array, move, unit "W in M", saving, bar.** As in the brief. A unit serves W arrays with M paid moves. Its
  saving is ln(mW/M)/ln m, and the bar is 0.1. At block 6 the bar allows 5.02 moves per array. 13 in 63 uses 4.85
  and 12 in 63 uses 5.25.
- **Commutative slice.** The engine sub-game with paid moves t_z, free shifts X^a and free linear mixing of registers.
  It has no address maps and no phase letters. Its letters commute, so a program is an identity in the group algebra
  C[A] of an abelian group A.
- **Characters, linear part λ, order-2 characters.**
  - A character sends each move t_z to a fourth root of unity.
  - Its *linear part* λ in F_2^m records which moves go to ±i: exactly those z with λ·z = 1.
  - The *order-2 characters* (λ = 0) send every move to ±1. They are the binary quadratic forms Q of ex2-families, and
    every shift acts trivially on them. There are 2^(m(m+1)/2) of them; the full group has 2^m times as many.
- **The order-2 system (\*)** (Lemma 1 of ex2-families). Every program has a "transfer triple" (Ã, B̃, Γ). The program
  is right at every order-2 character exactly when X_Q := Ã_D (I + Γ_D)^{-1} B̃_D equals I_W for every odd form Q and 0
  for every even form, where D is the set of moves acting as -1 at Q.
- **The matrix E and its pattern** (Theorem 10 of ex2-families). E = I + Γ - K, with K = B̃Ã.
  - Even forms make E's block on their support nonsingular.
  - Odd forms make the rank of that block equal to the slack: the size of the support minus W.
- **AG(n,2), flat, hyperplane.** AG(n,2) is the n-dimensional affine space over F_2, with 2^n points. A d-flat is an
  affine subspace with 2^d points; a hyperplane is an (n-1)-flat. The support F = {z : z1 = z2 = 1} at block m is a
  copy of AG(m-2, 2).
- **The local problem L(n, s).** Is there a nonsingular 2^n x 2^n matrix, indexed by the points of AG(n,2), whose
  block on every affine hyperplane has rank at most s?
- **Uniform family U_m.** One move in each of the 2^m - 1 directions, with W arrays. Full sharing is W = 2^(m-2).

## 1. Calibration (before any claim about an open shape)

| Case | Known truth | What this track's tools say |
|---|---|---|
| plain schedule at block 3 (1 and 2 arrays, exact triples) | exists | right at all 512 characters (error 0) |
| U_3 = 2 in 7 (skeptic's solution of (\*)) | exists, if (\*) suffices | right at all 512 characters with the free corrections; this is Theorem A |
| local problem L(3,1) (block 5, 7 in 31) | impossible (Corollary 13) | refused by the star recursion; Singer ansatz refused modulo 29 |
| L(3,2) (block 5, 6 in 31) | feasible as a matrix problem (Singer permutation, exact ranks) | not refused; Singer ansatz has solutions (dimension 0) modulo 29 |
| L(4,2) (block 6, 14 in 63) | impossible (Corollary 13) | refused by the star recursion |
| L(4,4) (block 6, 12 in 63) | feasible (bent circulant, exact ranks 4 on all 30 hyperplanes, rank 16) | not refused |
| U_4 = 4 in 15, U_5 = 7 in 31 | impossible (Corollaries 5 and 13) | refused (recursion rows m = 4, 5) |
| slice at block 4 below 4W (Theorem 17) | impossible | not re-run here (Theorem 17 uses the same pattern); at block 5 the star rows keep the plain 1 in 5 feasible |
| uniform point c = 1/7 at block 5 (W = 7) | impossible | block-5 LP with star rows: infeasible; 1/6 feasible; boundary at 3/20 |
| plain 1 in 5 inside the block-5 star LP | exists | feasible |
| engine 2 in 5 at block 3 (slice version) | impossible | refused by the deletion rows (rank cannot be negative), as in ex2-families |
| E8 unit, block-2 stand-ins, control word, absorb4 | units / not units outside the slice | blind: they use address maps or phases |
| numerical search for E (Levenberg–Marquardt, `loc.py`) | exact solutions exist at L(3,2), L(4,4), L(4,5) | unreliable: 1 of 21 runs at L(3,2) found one (condition number 2e2, the first version without scale fixing); 0 of 3 at L(4,5); the rest drift to singular E |

The exact solutions named in the table (Singer permutation, bent circulant) were checked with exact rational ranks.

## 2. Theorem A: the order-2 system is exact

**Statement** [proved]. A multiset of moves, in a time order, serves W arrays in the commutative slice exactly when
the order-2 system (\*) has a solution (Ã, B̃, Γ). The paid moves stay the same, and the free steps may use scratch
registers.

**Idea of the proof** (full proof in section 7.1).
- With a scratch copy, a shift and an addition, any combination of shifts can be applied to a register for free.
- The idempotents of the shift group split the group algebra into 2^m blocks, one for each λ.
- In block λ, each move t_z squares to (-1)^(λ·z). After a free scalar (-i)^[λ·z] it becomes an involution s_z, and
  the s_z satisfy exactly the relations of the order-2 block.
- So one order-2 program, with a fixed free correction after each move, is right in every block at once.
- Concretely, the correction after a move in direction z is u_z = ((1 - i) + (1 + i) X^z)/2. The outputs are then
  multiplied by the product over j of ((1 + i) + (1 - i) X^(e_j))/2.

**Check** [computed]. I took the skeptic's numerical 2-in-7 solution of (\*) and evaluated the program at all 512
characters of the block-3 slice group.
- Realised literally (no shifts), it fails at the 448 characters off the order-2 ones (error 1.99).
- With the corrections, it is right at all 512 (error 1.4e-15).
- The exact plain programs are right at all 512 either way.

**What it changes.**
- Every slice question is now exactly a question about (\*). Exclusions were already sound, since they use only
  necessary conditions.
- Existence statements in the slice no longer need g-play's reduction. They need only an exact solution of (\*).
- One of my own runs illustrates the trap. A per-λ solver that forced the empty-word term to be I at λ ≠ 0 "failed"
  at odd λ. That failure was an artefact: there is no trivial character in those blocks.

## 3. The converse question: does Theorem 10's pattern build a program?

**For a given program, no** [proved]. Theorem 10's bordered-matrix identity gives
rank(I_W - X_Q) = W + rank E_D - |D|.
- On odd forms, rank E_D = |D| - W is equivalent to X_Q = I_W. That part is exact.
- On even forms, "E_D nonsingular" says only that I_W - X_Q is invertible, that is, X_Q has no eigenvalue 1. The
  program needs X_Q = 0.
- Smallest example: block 2, one array, moves e1 then e2. Write K_12 = k.
  - Theorem 10 holds for every Γ_21 except 1/k.
  - (\*) needs Γ_21 = 2/k exactly.

**The exact strengthening** [proved]. On an even form, X_Q = 0 is equivalent to the bordered matrix
[[E_D + K_D, K[D, :]], [K[:, D], 0]] having rank |D|, given that Ã and B̃ have full rank W, which every unit forces.
- This involves K, not E alone.
- Together with E = (unit lower triangular) - K and rank K <= W, it is equivalent to (\*), hence (Theorem A) to a
  program.
- So "search for E, then build the program" is not a shortcut: the search space is the program's own.

**At the existence level** (does some E with the bare pattern imply some program?).
- No counterexample is known.
- At blocks 3 and 4 the bare pattern already reproduces the truth for the cost: no saving (deletion, Theorem 17).
- On the small cases of ex2-families the solver found a program whenever the pattern tests passed. Those were
  solutions of (\*), which Theorem A now makes genuine programs.
- I found no case of "pattern yes, program no", and did not settle the question.

**Why no E-search at block 6.**
- For 13 in 63, (\*) has about 3,600 complex unknowns and 2^21 forms, about 3.5 x 10^8 scalar equations.
- The local E-problem alone is in special position. Levenberg–Marquardt found a known exact solution once in 21 runs
  at 8 x 8 and never at 16 x 16; the other runs drifted to singular matrices (section 1).
- Section 4 then showed that 13 in 63 has no E at all, so there is nothing to search for there.

## 4. The exclusion: 13 in 63 and the uniform family

**Step 1: the local problem** (Lemma 12' and Theorem 10 of ex2-families, read again).
- In U_m with W arrays, E's block on F = {z1 = z2 = 1} is nonsingular, because z1 z2 is an even form.
- Every affine hyperplane of F lies in an odd (m-2)-flat avoiding 0, which carries 2^(m-2) moves and so has slack
  s = 2^(m-2) - W.
- So E's block on every hyperplane of F has rank at most s.
- If L(m-2, s) has no solution, U_m with W arrays does not exist.
- For 13 in 63 this is L(4, 3).

**Step 2: the star lemma** [proved, section 7.2].
- Let E be nonsingular on AG(n,2). Take a d-flat Ψ, and a k x k nonsingular block E[R, C] with R and C inside Ψ.
- The g = 2^(n-d) - 1 flats of dimension d+1 through Ψ split the rest of the space into disjoint pieces.
- The Schur complement on E[R, C] then has rank at most (rank bound on (d+1)-flats) - k on each "star" made of Ψ's
  remaining rows and one piece.
- Counting the rank that Ψ's remaining rows and columns need gives
  (g - 2) k <= g s_(d+1) - 2^(d+1),
  where s_(d+1) bounds the rank on (d+1)-flats.
- So the rank on d-flats is at most s_d = floor((g s_(d+1) - 2^(d+1)) / (g - 2)).

**Step 3: L(4, 3) has no solution** [proved].
- Planes of AG(4,2) (d = 2, g = 3): s_2 = 3·3 - 8 = 1. So E has rank at most 1 on every affine plane.
- Lines (d = 1, g = 7 planes through each line): s_1 = floor((7·1 - 4)/5) = 0. So E vanishes on every pair of
  points, hence everywhere. Contradiction.
- Alternatively: any 3 points lie in a plane, and Lemma 11 of ex2-families (with s = 1) makes E singular.

So **13 in 63 does not exist in the commutative slice** (Theorem C).

**Near-uniform multisets at block 6** [computed with the proved criterion; `nearU6.py`]. The same argument works for
U_6 plus extra moves Z. It needs an even 4-flat F that avoids Z, each of whose 30 hyperplanes lies in an odd 4-flat
with slack 16 + #(Z in it) - W of at most 3. Above the bar this excludes:

| W | extra moves t | moves per array | excluded |
|---|---|---|---|
| 13 | 1 | 4.92 | 62 of 63 (all but an extra move in direction 111111) |
| 13 | 2 | 5.00 | 1,922 of 2,016 |
| 14 | 1 / 2 / 3 | 4.57 / 4.64 / 4.71 | all / all (the last one, two extra moves 111111, by a mixed-bound refinement) / 43,586 of 43,680 |
| 15 | up to 2 / 3 | 4.27 to 4.40 | all / all but one |
| 16 | up to 3 | 3.94 to 4.13 | all |

I checked only t <= 3.
- The uncovered cases all involve the direction 111111: either it is the extra move, or the two extra moves add up
  to it.
- A refinement with mixed bounds and a perfect-matching test (`nearU6b.py`) did not close any of the 95 open W = 13
  cases.
- They need another even support or a new idea.

**Step 4: every block** (Theorem D; `recursion.py`). Iterate the bound s_d from hyperplanes down, and stop when
Lemma 11 applies (s_d <= d - 1) or a bound goes negative.

| block m | full sharing 2^(m-2) | now excluded | before (Cor. 13/15) | open above the bar | savings of the open members |
|---|---|---|---|---|---|
| 5 | 8 | W >= 7 | W >= 7 | none (Theorem 18) | — |
| 6 | 16 | W >= 13 | W >= 14 | **none** | — |
| 7 | 32 | W >= 26 | W >= 29 | 23 to 25 | 0.122 to 0.165 |
| 8 | 64 | W >= 52 | W >= 58 | 40 to 51 | 0.109 to 0.226 |
| 10 | 256 | W >= 209 | W >= 231 | 129 to 208 | 0.101 to 0.308 |
| 14 | 4096 | W >= 3354 | W >= 3687 | 1524 to 3353 | 0.100 to 0.399 |

The excluded fraction settles near 0.82 of full sharing; the old one was 0.9.

**The floor of this method** [proved; computed at n = 4].
- Let n be even and φ a bent function on F_2^n. The matrix E_xy = (-1)^φ(x+y) is nonsingular, and its rank on every
  affine hyperplane is exactly 2^(n-2) (for n = 4, checked with exact ranks).
- So L(n, 2^(n-2)) has solutions, and the local argument can never exclude W <= (3/4) 2^(m-2) when m is even.
- In particular it can never reach 12 in 63, and it cannot settle the family with a constant fraction below 3/4.
- At n = 5 (block 7) there is no bent function. Annealing found a Boolean function ψ on F_2^5 whose nonzero
  autocorrelations are all <= 0. The circulant H diag((-1)^ψ) H / 32 has exact rank 32 and rank at most 8 on all 62
  hyperplanes (exact over Q; `circ5.py`). So L(5, 8) has solutions too.
- Hence the local route cannot exclude 24 or 23 in 127. Only 25 in 127 (s = 7) is undecided locally. Circulants
  cannot reach s = 7: the autocorrelations would sum to at most -124, but their sum is at least -32.
- The true threshold of the local problem lies between 0.75 and about 0.82.

**Ansatz results these superseded.** Before finding Step 2, I had ruled out several special forms of E for L(4,3).
- Unitary matrices: a trace count gives 112 against at most 90 [proved].
- Circulant matrices: a Turán-type count of sign changes, 64 against the 75 needed [proved].
- Permutation matrices: average 3.73 per hyperplane [proved].
- Block anti-diagonal matrices: the ranks would need 16 against at most 12 [proved].
- Kronecker products [proved].
- The Singer-symmetric ansatz: no solution modulo 31 or 61 [computed].

These pointed the way and are now consequences of Theorem C.

## 5. Star inequalities in the slice LP (block 5)

The star lemma with k = 0 is a linear inequality on the ranks of principal blocks. Let F be an even support and split
it as Ψ plus pieces Q_1, ..., Q_g. Then

2|Ψ| + (g - 2) rank E[Ψ, Ψ] <= Σ_j rank E[Ψ ∪ Q_j, Ψ ∪ Q_j].

Each right-hand term is bounded by the slack of any odd form whose support contains it. This is the same kind of row
as the rank-cover rows of Theorems 17 and 18, but it is not implied by them.
- At block 5 I added two families of rows to ex2-families' LP, for every even 3-flat F:
  - "line + its 3 planes";
  - "point + its 7 lines".
- Calibration: the plain 1 in 5 stays feasible, and the uniform point needs c >= 3/20, which is W <= 6.67 r.
- The ten-piece cuts alone allowed W <= 7.2 r.

**Result** [proved: exact certificate; re-checked].
- The LP with the star rows, the plane bounds and the deletion rows (36,009 rows) has optimum exactly 93/20 = 4.65,
  at the uniform point c = 3/20.
- An exact rational dual certificate uses 116 rows: 29 star rows and 87 plane bounds. It proves M >= (93/20) W for
  every commutative-slice unit at block 5 (Theorem E).
- A separate program of mine (`verify_star_cert.py`) checked it from scratch. It rebuilds the 32,768 quadratic forms
  at block 5 itself, without code from ex2-families, and checks the meaning of every row:
  - each star row is a line plus the three planes through it, inside an even form's 3-flat support;
  - each bound is an odd form's support containing the plane.

  It then adds the rows in exact fractions to get Σ c >= 93/20.
- Adding ex2-families' rank-cover rows cannot lower the optimum, and cannot raise it either: the uniform point
  c = 3/20 satisfies every one of them (the ten-piece cuts need only c >= 5/36) [argued; the full run was stopped
  after 17 minutes as unnecessary].

**What it means.**
- At block 5 the bar needed at most 4.257 moves per array. The slice now needs at least 4.65, so its best saving at
  block 5 is at most ln(5/4.65)/ln 5 = 0.045.
- By projection (Lemma 3 of ex2-families) every slice unit at every block from 5 on needs at least 4.65 moves per
  array. At block 6 the gap is now 4.65 to 5.02.
- Rows of this kind alone cannot close block 6 [argued]. The analogous uniform point there needs only about 4.73 moves
  per array, and LP fractions lose the integer rounding (rank 1 to rank 0 on lines) that decides Theorem C.

## 6. Predictions, scored

| # | Prediction (written before computing) | Outcome |
|---|---|---|
| P1 | E-test reproduces every block 3-5 answer (0.8) | held, where it applies: the calibrations of section 1 pass. It is blind outside the slice |
| P2 | Theorem 10's conditions are not sufficient as stated; an extra condition is needed (0.6) | held: X_Q = 0 on even forms; the exact strengthening is (\*) |
| P3 | a non-linear use of E's rank pattern excludes 13 in 63 (0.45) | held (Theorem C) |
| P4 | no 13..16 in 63 slice program found and certified (0.85) | held (13 in 63 is now proved impossible) |
| P5 | outcome "excluded" or "calibrated-no-hit" (0.85) | held: excluded |
| P6 | L(4,3) is feasible, so the local argument cannot exclude 13 in 63 (0.55) | **wrong**: L(4,3) is infeasible |
| P7 | L(4,4) is feasible (0.8) | held (bent circulant) |
| P8 | the bare pattern admits E for some known-impossible small shape, while (\*) agrees with the truth (0.6) | half: (\*) is now proved exact (Theorem A); no "pattern yes, program no" case was found |

## 7. Proofs

### 7.1 Theorem A (the order-2 system is exact)

*Setting.* The slice group A is generated by the moves t_z, with z in F_2^m \ 0. Its characters are the functions
ψ = Σ_c [z·y_c] mod 4 (Z_4-valued) for points y = (y_1, ..., y_f), with χ_ψ(t_z) = i^ψ(z). The shifts X^a = t_a^2 form
Sh ≅ F_2^m, and χ_ψ(X^a) = (-1)^(λ·a), where λ = ψ mod 2 is linear (λ = Σ y_c). The even-valued ψ are exactly
2Q for binary quadratic forms Q. So every ψ with linear part λ is [λ·z] + 2Q(z) for a unique Q, which gives
⌊ψ(z)/2⌋ = Q(z).

*Free maps.* Copy a register into scratch, apply X^a there, and add back with a scalar. Repeating this applies any
element of C[Sh] to a register. So all C[Sh]-linear maps are free.

*Correction after each move.* Put u_z = ((1 - i) + (1 + i) X^z)/2 in C[Sh], and t'_z = t_z u_z. At ψ:
- if ψ(z) = 0, 1, 2, 3, then t'_z(ψ) = 1, 1, -1, -1;
- that is, t'_z(ψ) = (-1)^⌊ψ(z)/2⌋ = (-1)^Q(z), with Q the quadratic part of ψ.

Likewise v = Π_j ((1 + i) + (1 - i) X^(e_j))/2 takes the value i^|λ| at ψ.

*The program.* Start from a solution of (\*) and build its scalar program, with W + M registers (converse of
Lemma 1). Insert u_z, on the same register, after each paid move t_z. Multiply the outputs by v.
- At any character ψ the result is v(ψ) times the order-2 program evaluated at the quadratic part Q of ψ, which is
  i^|λ| (-1)^(tr Q) I_W.
- The target is χ_ψ(K) = i^(Σ_j ψ(e_j)) = i^(Σ_j [λ_j]) (-1)^(Σ_j Q(e_j)) = i^|λ| (-1)^(tr Q).
- They agree at every character, so the word is a formal identity in C[A], right for every number of digit columns.
- The paid moves are unchanged.

Conversely, (\*) is necessary (Lemma 1 and the point evaluations of ex2-families, checked by its skeptic). □

### 7.2 Star lemma

*Two general facts.* For any block matrix [[A, B], [C, D]],
rank [[A, B], [C, D]] >= rank [A B] + rank [A; C] - rank A.
(Let U be the column space of [A; C] and V that of [B; D]. The projection onto the top rows maps U + V onto the column
space of [A B] and U onto that of A, so dim(U + V) - dim U >= rank [A B] - rank A.) For column spaces U_1, ..., U_g that all contain U_0,
dim(U_1 + ... + U_g) <= Σ dim U_j - (g - 1) dim U_0.
Both were checked on 4,000 random integer instances (`check_star_ineq.py`).

*Statement.* Let E be a nonsingular matrix on a ground set X. Split X into Ψ and pieces Q_1, ..., Q_g. Let R and C
be subsets of Ψ with |R| = |C| = k and E[R, C] nonsingular, and suppose rank E[Ψ ∪ Q_j, Ψ ∪ Q_j] <= σ_j. Then
2(|Ψ| - k) <= Σ_j (σ_j - k) - (g - 2) a,   where a = rank of the Schur complement block on (Ψ \ R) x (Ψ \ C).

*Proof.*
- Let T = E / E[R, C], the Schur complement, with rows X \ R and columns X \ C. It is nonsingular.
- Since R and C lie inside Ψ ∪ Q_j, we have rank E[Ψ ∪ Q_j, Ψ ∪ Q_j] = k + rank T[(Ψ \ R) ∪ Q_j, (Ψ \ C) ∪ Q_j].
  So the star S_j = T[(Ψ \ R) ∪ Q_j, (Ψ \ C) ∪ Q_j] has rank at most σ_j - k.
- Write A = T[Ψ \ R, Ψ \ C], B_j = T[Ψ \ R, Q_j] and C_j = T[Q_j, Ψ \ C].
- The |Ψ| - k rows Ψ \ R of T are independent. Their columns are (Ψ \ C) plus all the Q_j, so
  |Ψ| - k <= Σ_j rank [A B_j] - (g - 1) a.
- By columns, likewise, |Ψ| - k <= Σ_j rank [A; C_j] - (g - 1) a.
- Add the two, and use rank [A B_j] + rank [A; C_j] <= rank S_j + a <= σ_j - k + a. □

*On flats.* In AG(n,2), the (d+1)-flats through a d-flat Ψ number g = 2^(n-d) - 1, and they split AG(n,2) \ Ψ into
disjoint pieces of 2^d points each. This was checked by brute force for n = 3, 4, 5 (`check_geometry.py`). Taking the
largest k gives the recursion of section 4.

### 7.3 Theorem C (L(4,3) has no solution) and the corollary

Let E be nonsingular on AG(4,2) with rank at most 3 on all 30 hyperplanes.
- Planes Ψ (d = 2): g = 3 and σ_j <= 3, so 2(4 - k) <= 3(3 - k), which forces k <= 1. Hence every 2 x 2 submatrix
  E[R, C] with R and C inside one plane is singular, and the rank on every plane is at most 1.
- Lines Ψ = {x, y} (d = 1): the 7 planes through the line split the other 14 points into pairs, and σ_j <= 1. With
  k = 1 the lemma would give 2(2 - 1) <= 7(1 - 1) = 0, which is false. So every entry E_xy is 0, whether x = y or
  not, since every pair of points (and every point) lies on a line. Then E = 0. Contradiction.
- (A second route from the plane step: every set of at most 3 points lies in a plane, so all principal submatrices of
  order at most 3 have rank at most 1, and Lemma 11 of ex2-families with s = 1 makes E singular.)

*Corollary.* A U_6 slice unit with 13 arrays would make E[F, F] such a matrix (Step 1 of section 4, using Theorem 10(b)
and (c) and Lemma 12'). So none exists.

Only Theorem 10(b) and (c) and Lemma 12' of ex2-families are used (plus, on the second route, Lemma 11), and
v-ex2-families checked all of them.
The time order, the triangular structure of Γ and the bound rank K <= W are not used.

## 8. Next

1. **A skeptic for Theorems A, C and E.** All three are short. Theorem C changes the map: the uniform family has no
   member above the bar at block 6. Theorem E raises the slice floor at every block from 5 on to 4.65 moves per
   array.
2. **The block-7 members 23, 24 and 25 in 127** (savings 0.122 to 0.165) are now the smallest open members above the
   bar.
   - The local problem there is L(5, s) with s = 9, 8, 7.
   - L(5, 8) has a solution (section 4), so 23 and 24 in 127 need the global E: other even supports, and the
     triangular structure that Theorem C does not use.
   - 25 in 127 (L(5, 7)) might still close locally with a sharper star count.
3. **Block 6 for multisets that are not uniform.** The gap is now 4.65 to 5.02. LP rows lose the integer rounding
   that decides Theorem C, so the natural next tool is an integer version: a fixed W, with integer rank variables on
   the flats of each even 4-flat. Theorem E itself needs a skeptic first.
4. **Below 3/4 of full sharing** the local problem has solutions (the bent matrices), so the constant-fraction family
   needs a global argument or a construction. A construction would be an order-2 solution for U_m with about
   (3/4) 2^(m-2) arrays. By Theorem A, any exact one is a genuine slice unit.

## Sources

All paths are under `the archive's checks/wht43/` unless absolute. This folder: `exist3/e3-slice6/`.

- **Brief and reading:**
  - `exist3/BRIEF.md`;
  - `exist2/NOTE.md`;
  - `exist2/ex2-families/NOTE.md` and `PROOF.md` (Lemma 1, Theorem 10, Lemma 11, Lemma 12', Corollaries 13 and 15,
    Theorems 17 and 18);
  - `exist2/v-ex2-families/NOTE.md`;
  - `exist2/ex2-modp/NOTE.md`;
  - `meta3/NOTE.md` ("In short");
  - `moonshot/BRIEF2.md`;
  - `meta/MAP.md` section 1.3;
  - `proofs/logging/PROOF.md` section 5 (g-play's reduction).
- **Programs:**
  - `loc.py`: Levenberg–Marquardt on L(n, s), the unreliable numerical E-search;
  - `perm_check.py`: permutation matrices on AG(3,2);
  - `exact_checks.py`: exact ranks of the Singer permutation and the bent circulant;
  - `singer_ansatz.py` and `singer_gen.py`: the Singer ansatz, run in Singular over F_p;
  - `allchars.py` and `allchars2.py`: all characters, without and with the corrections of Theorem A;
  - `lamsolve.py`: the per-λ solver, whose odd-λ verdicts are void (section 2);
  - `recursion.py`: the star recursion and the table of section 4;
  - `nearU6.py` and `nearU6b.py`: the near-uniform multisets at block 6;
  - `circ5.py`: the block-7 circulant with hyperplane ranks <= 8;
  - `check_star_ineq.py` and `check_geometry.py`;
  - `starlp/rc_star5.py`: ex2-families' `rc_lp5.py` plus the star rows, with `refl.py` and `exact_cert.py` copied
    from ex2-families; `starlp2/` holds the star-only run that produced the certificate;
  - `verify_star_cert.py`: the independent check of that certificate.
- **Outputs:**
  - `out_exact_checks.txt`;
  - `out_allchars2_u3.txt` and `out_plain_triples.txt`;
  - `out_recursion.txt`, `out_nearU6.txt`, `out_nearU6b.txt` and `out_circ5.txt`;
  - `out_check_star_ineq.txt` and `out_check_geometry.txt`;
  - `sing1/`, `sing2/remote-out/ssh4/` (Singular);
  - `lm4b-partial/` (Levenberg–Marquardt tails);
  - `starlp/remote-out/out_calib_star5.txt` (LP calibration);
  - `starlp2/remote-out/out_star5_noold.txt`, `star_cert_b5.json`, `star_vars_b5.json` (Theorem E's LP and
    certificate);
  - `out_verify_star_cert.txt`.
- **The 2-in-7 solution used for the check:** `exist2/v-ex2-families/u3_solution.npz`, copied to `u3/`.
