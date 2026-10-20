<!-- header: existence-routes -->
# The skeptic's record: families commutative slice

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** This is an AI skeptic's record of its check, as it wrote it. Its verdict is the status this revision uses for the track it checked, unless a later verdict is named in the section README.

---

# v-ex2-families: skeptic's check of ex2-families (sessions 1 and 2)

15 October 2026, 00:24 to about 00:45 EDT. Everything ran on the Mac, on one core with small memory. No Colab was used and
nothing was posted. Dated predictions are in `STATUS.txt` and every run is in `RUNLOG.txt`. All programs here are my
own. The only thing I read from the track's folder was the variable list and certificate of Theorem 18, and only to
re-check that certificate.

## Verdict

**Not refuted.** Every exclusion I checked holds, given the model of the commutative slice.
- Theorem 17 (block 4: M >= 4W) and Theorem 18 (block 5: M >= 155W/36) both come out again from my own rebuild of the
  linear program, and both bounds are exact.
- The hand proof of Theorem 17 checks exactly.
- Theorem 10 holds on genuine units that share moves, found by my own solver.
- Lemma 11, Theorem 12, Lemma 12' and Corollaries 13, 15 and 19 hold up line by line.

**Calibration holds, in the limited sense the track states.** Inside the slice the method excludes no genuine unit I
built. The plain unit and the 1-in-4 unit sit exactly on the block-4 optimum. Every known-impossible slice case is
refused. Outside the slice (E8, address maps, phase letters) the method is blind, and the track says so.

## What I checked, and how

1. **The model (the likeliest weak point).**
   - Every binary quadratic form is a point evaluation of the f-column program. I checked this with exact mod-4
     arithmetic on 900 random forms at blocks 3 to 5. At those points every shift acts trivially and the kernel K acts
     by (-1)^(tr Q).
   - The transfer form of Lemma 1 matches the explicit product to 2e-14 on 200 random programs. These programs have
     rectangular mixing, scratch registers and random move orders.
   - The algebra of Lemma 1 needs only that each paid move is "free map + (chi(t_z) - 1) x rank one". So it also covers
     moves written as additions, such as c_s += λ t_z c_r.
   - The results use only g-play's Hadamard picture (each direction acts as diag(i^[z.y])). They do not use g-play's
     orbit-game reduction, so the earlier proviso "given g-play's reduction" does not apply here.
   - Hypothesis to keep in view: the theorems are about units that are right for every number of columns f. The forms
     with cross terms need many columns; my gadget uses 7 per cross term.
2. **Theorem 10** (the matrix E). I re-derived it by hand: two Schur complements of one bordered matrix.
   - Numerically I tested it on three genuine units that share moves, all found by my own solver: U_3 = 2 in 7 at block 3
     (residual 2.5e-15 on all 64 forms), U_3 + e4 at block 4 (2 in 9), and U_3 + e4 + e5 at block 5 (2 in 11).
   - The unit identity (*) holds to 1e-15 on all 64, 1,024 and 32,768 forms.
   - Even supports give nonsingular blocks (smallest singular value at least 0.148). Odd supports have rank equal to the
     slack, with 0 mismatches. E vanishes on the tight supports.
   - Every rank-cover piece bound holds on these units (945 and 15,585 checks).
3. **Theorem 17.**
   - My own LP, with direction cells on the 49 even planes and slack bounds from the 56 odd planes, has optimum exactly 4.
     My exact rational certificate has 41 rows and its bound is exactly 4.
   - The optimum is attained at the plain unit. With the ten-piece rows and deletion rows added, it is attained at
     1 in 4 = (7, 11, 13, 14).
   - The hand proof (weights 18 and 17 on the 7 even planes of the even-weight space) gives coefficient 492 on every
     direction and constant 1968, exactly as stated.
4. **Theorem 18.**
   - My own LP (225 even 3-flats, cells plus all 21 ten-piece cuts each, 240 odd 3-flats) has optimum 4.305555555556,
     at the uniform point 5/36.
   - My exact rational bound is 155/36 - 1e-15, above the bar's 5^0.9 = 4.2567.
   - I also re-read the track's 194-row certificate with my own semantics and geometry. Every row is a valid
     inequality: each slack row's 8 directions form an odd form's support that contains the piece, and each ten-piece
     row tiles F x F exactly over an even form. The rows sum to exactly sum c >= 155/36.
5. **Lemma 11.** I checked the proof line by line. Both cases of the induction are fine: a Schur complement when a small
   nonsingular principal block exists, and otherwise a chordless shortest cycle. A 4-cycle permutation matrix shows the
   order bound s + 2 cannot drop to s + 1.
6. **Theorem 12 and Lemma 12'.**
   - Theorem 12 follows from Theorem 10 and Lemma 11.
   - I checked the hand proof of Lemma 12'.
   - Brute force at blocks 4 to 6 (10, 50 and 306 subflats): every subflat of F of dimension at most m - 3 lies in an
     odd (m-2)-flat avoiding 0 with exactly one basis vector. F itself does not.
7. **Corollaries 13, 15, 16 and 19.**
   - The ten pieces of Corollary 15 tile F x F, and each piece's union is a hyperplane of F. My LP asserts both at
     block 5.
   - The arithmetic of the excluded members is right: 14 to 16 in 63, 29 to 32 in 127, 231 and up in 1023, every
     above-bar U_5^(r), every saving U_4^(r), and V_m above 0.9 * 2^(m-3).
   - Projection (Lemma 3) keeps W and drops moves that project to 0, so Corollary 19 follows.
8. **Corollary 10' (the even test).** My own even test excludes the four 4-in-17 multisets (U_4 + {e_j, complement})
   with 72 failing even supports each. It excludes the block-4 survivor 3U_4 + {7, 8} with 192. It passes 4 x (1 in 4)
   + 15, U_3 + e4 e4 and plain 3 in 12.

## Consistency with earlier results

The new bounds agree with every earlier slice result:
- block 3: no saving (wave 16);
- block 4: 2 in 7 and 3 in 10 fail (now below 4W);
- block 5: 2 in 8 and 3 in 11 fail (below 155W/36); 2 in 9 passes Theorem 18 and is excluded only by the exhaustive
  deletion test;
- meta2's open block-5 case 4 in 16 is now closed (16 < 17.2).

No known slice unit falls below 4W at block 4 or below 155W/36 at block 5.

## Small points (none changes a result)

- Corollary 19's "needs > 4.306 moves per array" is a rounding. The bound is 155/36 = 4.30556, which is below 4.306. It
  should read "at least 155/36 = 4.3056".
- Calibration (a) is weaker than "1,400+ units" suggests. The 880, 404 and 122 random units are direct sums of
  single-array units with idle moves, which share nothing. The genuine sharing units above close that gap.
- Calibration (b) for the new tests rests only on solver stalls, as the track says. The impossible slice cases with
  independent proofs (2 in 5 at block 3; the exhaustive block-4 and block-5 tables) are all refused.
- Not re-run by me: the 10-piece ILP optimality, the weight-12 remark (a side note, not used), the F_2 exhaustive check
  of Lemma 11 (I checked the proof instead), and the solver-completeness runs. None of these carries a claimed theorem.

## Files

- Programs:
  - `geo.py` (flats, forms);
  - `my_rclp.py` (rank-cover LP plus exact bound);
  - `my_hand_b4.py`;
  - `my_lemma1.py`;
  - `my_solve_u3.py` and `my_polish.py`;
  - `my_thm10_check.py`;
  - `my_lemma12.py`;
  - `my_even_test.py`;
  - `my_verify_their_b5.py`.
- Outputs:
  - `out_my_*.txt`;
  - `my_cert_m4_cells.json`, `my_cert_m5_cells-ten.json` (my own exact certificates);
  - `u3_solution.npz`.
