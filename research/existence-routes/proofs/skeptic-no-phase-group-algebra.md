<!-- header: existence-routes -->
# The skeptic's record: no phase group algebra

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** This is an AI skeptic's record of its check, as it wrote it. Its verdict is the status this revision uses for the track it checked, unless a later verdict is named in the section README.

---

# v-e3-rmbuild: skeptic's check of e3-rmbuild

15 October 2026, 07:01 to about 07:40 EDT (Mac clock). I am the skeptic of track e3-rmbuild of exist3
(`../BRIEF.md`). I wrote predictions first (`STATUS.txt`). Nothing was posted.

## Verdict

**Refuted in part. The calibration holds.**

**What stands.**
- No unit was found, in any game.
- Theorem A and Lemmas S, R, C, V and N are correct. I checked each proof line by line, and Lemma R and the lift on
  raw arrays with my own code.
- Every exact number e3-rmbuild reports for its own address groups comes out the same in my own code.

**What is false.** e3-rmbuild generalises from its three address groups to all of them, and that generalisation is
wrong:
- "With any cyclic address group ... the tight cliques drop to or below the bound Wn" (obstruction field) is false.
- So are "the slice's obstruction disappears once address maps are allowed" (result field) and "the known exact tools
  are blind to U_4 and U_5 once address maps are allowed".

**Why it fails.** The tests of Lemma V depend on the address group *and* on how it sits relative to the coordinate
frame. The target is the identity matrix I, and conjugating H turns that target into a different symmetric matrix.
For many cyclic address groups the old tight-clique test still excludes 4 in 15, with the same shape (one move per
direction) that e3-rmbuild tested:

| Cyclic address group H (block 4) | Groups tested | Excluded by tight cliques | Largest tight clique vs bound Wn |
|---|---|---|---|
| order 3 without fixed vectors (all conjugates of e3-rmbuild's Singer power) | 56 | **36** | 27 or 36 vs 12 |
| order 5 (all conjugates of e3-rmbuild's Singer power) | 336 | **48** | 50 vs 20 |
| order 15 (all Singer subgroups) | 336 | 0 | no tight form at all |
| coordinate permutations: 3-cycle, 4-cycle, the two kinds of involution | 4 | **4** | 45 vs 12, 60 vs 16, 30 vs 8 |

The permutation case is clear in theory too. When every h in H is orthogonal (h h^T = I), each form is odd either on
every channel or on none. The slack of every form is then n times its slack in the slice, so the slice's obstruction
carries over unchanged. e3-rmbuild's own groups (the Singer powers in the polynomial basis of x^4 + x + 1) do pass. So
do the matching groups for x^4 + x^3 + 1, which I also ran.

**Two smaller qualifications.**
1. The calibration "Singer maps of order 7 at block 3: deletion excludes 2 in 5 (10 < 12)" holds for e3-rmbuild's
   group, and for 6 of the 8 Singer subgroups of GL(3,2). For the other 2, no form is tight and nothing is excluded,
   although engine 2 in 5 at block 3 is known to be impossible. The refusal half of this calibration therefore depends
   on which Singer group is used.
2. The coverage wording needs two corrections.
   - Rank cover was run over partitions into at most 3 blocks only for supports of at most 8 directions. Larger
     supports got at most 2 blocks, plus the all-singletons partition. The blocks are sets of directions, not of
     virtual moves.
   - At order 5, the Levenberg-Marquardt solver missed the plain-cost control 4 in 16 in its one finished restart
     (residual 509.5, and 6,327 on all characters; `out_block4.txt`, rmG_g5c). The calibration summary does not say
     this. It says only that the order-5 runs gave "no usable numerics".

**What this changes.** The no-hit outcome stands, and is if anything stronger: U_4 is now excluded exactly in the
no-phase game for many small cyclic address groups. The useful correction is to the picture of what is left open.
What escapes the exact tests is specific choices of address group relative to the frame, among them every Singer
group of order 15. Adding address maps does not by itself make the tests fail.

## 1. Theorem A, line by line

- **Characters of the slice A.** At f digit columns, in the Hadamard basis, dir(z) acts as i^s and X^z as (-1)^s,
  where s counts the columns c with z.y_c = 1.
  - The integer s mod 4 equals [z.mu] + 2 e_2(z.y_c), where mu is the sum of the columns.
  - So every character has the form chi_{mu,Q}, with Q(z) = e_2(z.y_1, ..., z.y_f), a binary quadratic form.
  - These characters are closed under products, because [a] + [b] = [a + b] + 2ab.
  - So A/X has order 2^C(m+1,2) = |Sym_m|, and dir(z) X maps to zz^T. This is the identification psi needs.
- **psi** (sector 0). e_0 is central: A is abelian, and P_h permutes the shifts. On that block c_z = 1, so
  dir(z) e_0 = d'_z e_0. A ∩ H = 1, since a permutation matrix that is diagonal in the Hadamard basis is the identity.
  Correct.
- **phi** (the lift). It is well defined, because <Q, zz^T> = Q(z) is linear in zz^T and the d'_z act by
  (-1)^{Q(z)} (Lemma R(b)). It respects the semidirect product, because P_h c_z P_h^-1 = c_{hz}. The output
  correction is free and invertible (c_z^-1 = c_z X^z). Correct.
- **Legitimacy of the letters.** c_z is a mixing with coefficients (1 ± i)/2 of a register and its shifted copy.
  moonshot/BRIEF2.md says that any complex coefficients count for the challenge. The lifted words use the same number
  of paid moves as the words they come from.
- **Lemma C, the transfer form (\*\*).** I re-derived it myself. A move acts as I - 2 Pi on the induced
  representation. Write w_t = Pi_t v_t. The trivial character forces the free part to be 1. Rescaling gives
  `Atil Pi (I + Gamma Pi)^-1 Btil = O_Q (x) I_W`, the same equation.
- **Lemma V.**
  - (a) Deletion is a rank count.
  - (b) The bordered-matrix computation is right; L = I + Gamma[D,D] is invertible.
  - (c) If a form is tight, Atil_D vanishes off the odd rows and Btil_D off the odd columns. So K[D,D] = I + Gamma[D,D].
    Gamma has zero diagonal blocks, so K restricted to the clique is unit lower triangular.
  - (d) Each rectangle is a submatrix of E[D_G, D_G].
  - All four parts are correct.
- **Lemma N.** |Sp(6,2)| = 1,451,520, and 256 x 1,451,520 = 371,589,120, as stated. N is normal and lies in the free
  part, so e_N is central. As a necessary condition the lemma is fine.

## 2. My own computations

All of this is my own code. The only things I took from e3-rmbuild are the statements in PROOF.md.

**Raw arrays (`my_raw.py`, block 2 at f = 1 to 8, block 3 at f = 1 to 8).**
- Lemma R: every direction at every f. The largest error is 5.6e-16.
- Fano: the product of the 7 corrected moves is 1.
- P_h d'_z P_h^-1 = d'_{hz}.
- The lifted U_3 = 2 in 7 is exact at every f. The same word without the correction is off by 1.2 to 1.7.
- The two-channel 2 in 6 with an address map is exact.
- The control with a non-orthogonal map is exact at f = 1 only, and off by about 1 from f = 2 on. This matches
  e3-rmbuild.
- At block 3, f = 1 to 6 ran on the Mac and f = 7 to 8 on Colab 4 (`job_raw78/remote-out/ssh4/`).

**Lemma V deletion and tight-clique tests (`my_vtests.py`).** These reproduce every number e3-rmbuild gives:

| Case | Theirs | Mine |
|---|---|---|
| order 7, block 3, 2 in 5 | excluded, 7 forms with 10 < 12 | the same, for both primitive polynomials |
| order 7, block 3, 2 in 6 | passes | passes (7 tight forms, clique 12 of 14) |
| slice, block 4, 4 in 15 | excluded, clique 15 > 4 | the same (56 tight forms) |
| order 15, 4 in 15 and 4 in 16 | no tight form | the same |
| order 5, 4 in 15 | 5 tight forms, clique 20 = 20 | the same |
| order 3, 4 in 15 | 21 tight forms, clique 12 = 12 | the same |
| slice, block 5, 8 in 31 | 240 tight forms, clique 31 > 8 | the same |
| order 31, block 5, 8 in 31, 7 in 31, 8 in 40 | no tight form | the same; also no tight form in 300 random conjugates of the order-31 group (`my_b5_conj.py`) |

**Soundness of the tests (my implementation and Lemma V itself).** Known units must never be excluded:
- the plain schedule 4 in 16 under every conjugate of order 3, 5 and 15, and 3 in 12 and 1 in 4 under every
  conjugate of order 3;
- U_3 = 2 in 7 under all 8 Singer groups and the permutation groups;
- the slice unit 2 in 9 = {e1, e2, e3, e4} + {3, 5, 6, 7, 8} under all 56 + 336 conjugates and six other groups.

None of them was excluded (`out_my_soundness*.txt`, `out_my_conj_sweep_plain.txt`). On my first attempt I passed
the tests a direction multiset that is not a known unit; that entry in `out_my_soundness.txt` says nothing.

**Conjugate sweeps (`my_conj_sweep.py`, `my_gen_test.py`).** These give the table in the verdict. The certificates
are printed and re-verified pair by pair in `my_cert_check.py` -> `out_my_cert_check.txt`. For example, take
H = <h> with h e_1..h e_4 = 0110, 1111, 1000, 1100 (order 3, no fixed vector, a conjugate of e3-rmbuild's group):
- 15 forms are tight;
- the 27 virtual moves of directions 2, 4, 5, 8, 11, 12, 13, 14, 15 (3 channels each) lie pairwise in common tight
  supports, with 13 witness forms;
- 27 > Wn = 12, so 4 in 15 is impossible in that no-phase game.

**Other checks.**
- GL(m,2) orbits on binary quadratic forms (`my_orbits.py`): 1, 7, 7, 21, 28 at m = 3, and 1, 15, 35, 105, 168, 280,
  420 at m = 4, as claimed. The m = 4 sizes also agree with a hand count by the rank and type of the form.
- `phi_check.py` re-run: 8 of 8 as expected, identical output.
- `out_block4.txt` matches the numerical claims: targets 0 of 7, with all-character residuals 3,538.7 to 12,932.3;
  control 2 of 3, at 7.2e-25 and 8.9e-23.
- No fit was claimed, so there was nothing to re-test on fresh inputs.

## 3. Not checked

- tr_sym.py: the bound of 1.21 moves per array from orbit unions.
- Proposition P, which is labelled [argued].
- Whether some finer partition of the virtual moves makes rank cover exclude e3-rmbuild's own groups of order 3 and 5.

## 4. Files

All of these are in this folder:
- predictions and log: `STATUS.txt`, `RUNLOG.txt`;
- programs: `my_raw.py`, `my_vtests.py`, `my_conj_sweep.py`, `my_gen_test.py`, `my_cert_check.py`, `my_b5_conj.py`,
  `my_orbits.py`;
- outputs: `out_my_*.txt`, `out_rerun_phi_check.txt`, and `job_raw78/remote-out/ssh4/out_my_raw_m3_f78.txt`.
