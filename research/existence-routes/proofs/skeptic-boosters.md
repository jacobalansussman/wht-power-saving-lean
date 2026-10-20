<!-- header: existence-routes -->
# The skeptic's record: boosters

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** This is an AI skeptic's record of its check, as it wrote it. Its verdict is the status this revision uses for the track it checked, unless a later verdict is named in the section README.

---

# Skeptic's check of the booster screen (exist3, track e3-boost)

15 October 2026, 05:38 to 06:12 EDT. One agent (v-e3-boost). Nothing was posted. No Lean. Mac jobs ran only
under a memory-capped runner (each under 2 minutes, well under 1 GB). Two Colab runs, both finished and copied back. A
third launch on ssh (Colab #1) could not be reached after it started; nothing of mine is known to be running there,
and its job (a rerun of `lp_combo3.py 81 243`) was repeated on ssh5.

## Verdict in short

1. **The exact core stands.**
   - Theorem B1 (no booster whose jobs all have at most 3 rows, on any register) is right: I checked the identity
     exactly, the positivity of the three coefficients, and the validity of the two certificates it uses.
   - Corollary B2 (such relations never amplify one another) is right: I redid both Wronskians by hand. They match
     the NOTE's closed forms, and so do the sign arguments.
   - The walk condition is a valid necessary condition. The NOTE's written reason for its lower bound is incomplete;
     section 3 gives the missing step.
2. **The calibration holds.**
   - My own code reproduces the hub-word roots to every printed digit.
   - The realisable relations pass my own certificate families too.
   - The made-up booster's 53.1568 at n = 243 reproduces. That is a rerun of exist2's script, not a rebuild.
3. **The survivors are what the NOTE says.**
   - S6, S7 and S8 are exactly tight on the frame count and the plain count.
   - Their windows are exact: a sign-rule count shows D has no other roots.
   - Each passes about 65,000 of my own tests on registers up to 14 rows. These include a new certificate family:
     incidence with all lines of the label space, isotropic or not.
4. **One headline exclusion is not established: "any booster needs a job of at least 6 rows".**
   - "All jobs of at most 5 rows, registers up to 12" is labelled [computed]. With rigorous certificates only, the
     LP allows such boosters on every register from 6 to 12. So the statement rests on the argued projections, like
     the 5-row statement does.
   - Registers above 12 are not covered; the NOTE's own coverage section says so. Yet the summary concludes "so any
     booster needs a job of at least 6 rows". e3-boost's certificate library itself allows boosters with every job of
     at most 4 rows, and of at most 5 rows, on registers of 17 and 20 rows. It simply has no certificates there.
   - My own families (symmetric and computable at any block) close the case of jobs of at most 4 rows on registers
     of 17, 24 and 32 rows. They do not close jobs of 5 rows.
   - So "structure forced" holds only in the restricted form: on registers of at most 16 rows, with the argued
     projections.
5. **New, against the survivors** [derived here; short proof in section 5; nobody else has checked it].
   - In the whole-array engine game, a program that is exactly tight on the frame count has a last paid block whose
     rank equals some target rank.
   - Every listed survivor has disjoint call and task ranks. So none of them exists as an exact program without
     catalysts.
   - The NOTE's "what a construction must satisfy (no waste: exactly frame-tight)" therefore describes an impossible
     object.
   - The survivors are not removed from the screen. Versions with catalysts or with slack are untouched, and they are
     what the duality uses.

## 1. What I checked, and how

| Claim | My check | Result |
|---|---|---|
| B1 identity d^θ = λ_f + λ_p d + λ_A A(d), A = 4, 6, 7 | coefficient by coefficient in (1, 2^θ, 3^θ), exact (`vB1.py`) | holds |
| λ_f, λ_p, λ_A > 0 on (0, 1) | Descartes' rule for exponential sums (two sign changes, so at most two roots) plus the exact endpoint data | holds; values as stated (0.2877, −0.3398, 0.1178) |
| A and the frame count valid on any register | Theorem R, w-bounds §2; A is scale-free (2^M − 2^(M−d) ∝ 1 − 2^−d) | holds |
| B2 Wronskians | by hand: λ_A′λ_p − λ_Aλ_p′ = −a x + b y + (a − b) x y = x y g, and λ_f′λ_p − λ_fλ_p′ = 7a x − 6b y + 4(b − a) x y = x y h; numeric error 1e-15 | holds; g < 0, h < 0, h′ root 0.3802 |
| Walk condition | proof re-derived (section 3); my own implementation with the capped recursion | valid; it kills 2[5] ≤ 2[1] + [8] and [4] + [5] ≤ 2[1] + [7], which pass all my price tests |
| Hub roots | my own 50-digit bisection (`vbig.py`) | T74E 1.4040079421607e-3, inside [root_lower, root_upper] of the price file; #1042 1.6749890159376e-3 |
| Made-up booster | rerun of `exist2/ex2-duality/lp_combo3.py 81 243` on Colab ssh5 (`cjob2/remote-out/ssh5/remote.log`) | 30.1728, 53.1568; ln 53.1568 / ln 243 = 0.7233 |
| Certificate profiles | my own A, C_k, and distance sets from the dual polar graph's intersection array; brute force over all frames at m = 2, 3 (`vcert.py`); sphere sizes add up to the frame counts for m ≤ 8 | all 2,603 of certlib's A, C_k and Dist profiles equal mine (`vcmp.py`) |
| Projection closed forms | compared row by row with `proj/N5.txt`–`N8.txt` | all match |
| My families on realisable relations | plain, A5, Kronecker, undo (37 relations), registers m to m + 4 | no failure (`out_vvalid.txt`) |
| Distance sets on real units | parity sets, {1}, {0, 1} at block 73 on T74E and #1042 | pass |
| No booster on registers ≤ 4 | my own LP with only my families (frame, A, C_k, distance sets with every S, lines with every T, blocks m to m + 3, and the scaling b = 2), with and without hooks (`vlp.py`) | confirmed, **even without the hooks** |
| Register 5 | same LP, up to block 11 (Colab, `cjob/remote-out/ssh5/remote.log`) | my families do not suffice; consistent with the NOTE's "needs the argued projections" |

## 2. The headline exclusion and its labels

- `lp3.py` uses rigorous plus argued certificates. Rerun with rigorous only (`vlp3_rig.py`, `out_vlp3_rig.txt`):
  - ranks ≤ 4 on registers 6 to 12: no booster;
  - ranks ≤ 5: boosters allowed at every grid point, on every register from 6 to 12.

  So "all jobs ≤ 5 rows, registers ≤ 12" is [computed, given argued certificates], not [computed].
- Beyond 12 (`vlp3_big.py`, `out_vlp3_big.txt`; e3-boost's own library, rigorous plus argued):
  - registers 13, 14, 16: no booster with ranks ≤ 5;
  - registers 17, 20: boosters allowed with ranks ≤ 4 (29 of 50 grid points, from 0.44 up) and ranks ≤ 5 (all grid
    points).

  The cause is that the library stops at block 16 for distance sets and projections, not any construction.
- With my families at blocks 8, 12, 17, 24 and 32 (`vlp_big.py`, `out_vlp_big.txt`):
  - ranks ≤ 4: no booster at any of them;
  - ranks ≤ 5: boosters allowed at 12 and 17. The worst direction is about
    −0.10[1] + 0.19[2] + 0.42[3] − [4] + 0.49[5], the family of the 5-row survivor [2] + 5[4] ≤ 4[3] + 2[5], which
    only the argued projections remove.
- **So the honest statement:**
  - no booster with jobs of at most 3 rows on any register [proved];
  - none on registers of at most 4 rows [computed; rigorous];
  - none with jobs of at most 4 rows on registers of 17, 24 or 32 rows [computed here; rigorous; the other large
    registers are untested];
  - none with jobs of at most 5 rows on registers of at most 16 rows [computed, given argued projections];
  - jobs of 5 rows on registers above 16: open.

## 3. The walk condition: the missing step

The NOTE derives the lower bound 2·max R − ΣR ≤ k from "the distance from L_0 moves by at most r_i at each paid step".
That alone gives only k ≤ ΣR.

The step that is needed:
- Write x_i = d(L_0, u_i L_0). The free letter f_i fixes L_0 and is an isometry, so x_i = d(L_0, P_i u_(i−1) L_0) =
  d(P_i^−1 L_0, u_(i−1) L_0).
- P_i^−1 L_0 is at distance exactly r_i from L_0.
- So |r_i − x_(i−1)| ≤ x_i ≤ r_i + x_(i−1).
- The set of distances reachable through a multiset R is then the polygon interval [max(0, 2 max R − ΣR), ΣR], capped
  at M. The cap does not lower the left end.

With that step the condition holds as coded. I reimplemented it with the exact capped recursion; it agrees on every
case I ran.

## 4. Smaller corrections

- *"B1 contains exist2's exclusion of the spreading family as a special case."* B1 covers only members whose jobs
  have at most 3 rows, that is Rs = 2[2] ≤ [1] + [3]. The rest of the family (for example 2[4] ≤ [1] + [7]) still
  needs exist2's argument with A at every block.
- *"[1] + 5[3] ≤ 4[2] + 2[4] is killed by the block-4 hook alone."* In e3-boost's own library it fails 108 tests,
  including 28 distance-set certificates at block 4, which are derived and need no mod-p values
  (`out_vcertlib_A4tight.txt`). The parity set S = {1, 3} at block 4 gives −909. So the 4-row result does not need the
  hooks at all, which my LP confirms.
- *S6 "would give a saving of at least 0.32".* The root is 0.6806, so the saving is at least 0.319.
- *Survivor counts "2 / 8 / 31 at blocks 6 / 7 / 8".* The block-7 list contains S6, a 6-row shape, and the block-8
  list contains 7-row shapes. The counts are of vectors per register, not of new shapes.
- *The window filter in `intenum.py`* uses 41 grid points on [0.9, 0.9983]. It could drop a relation that is negative
  only between grid points. That is minor: such a relation would be a very weak booster.
- *Exponent check.* For pure units the job-model root and the pipeline's root solve the same equation
  Σ n_r r^θ = W m^θ. Their agreement is a consistency check of the data, not independent evidence that the two
  exponents agree in general. The NOTE does mark the general step as argued.

## 5. New: exactly frame-tight programs end on a target's rank

**Lemma L** [derived here; nobody else has checked it]. The setting is the engine game with all arrays of M rows,
each paid block acting on one whole array, and scratch arrays allowed. Suppose a correct program serves W live arrays
with targets of ranks k_1, …, k_W ≥ 1 and uses exactly W paid blocks. Then the last paid block's rank equals some k_l.

*Proof.*
- *The potential.* Take the frame permutation module V = Q(i)^(frames) and substitute y_l = t_l·e_(L_0). The content
  of array R at frame U is a linear form φ^R_U(t). Let Σ be the span of all φ^R_U with U ≠ L_0. Theorem R's potential
  for the frame count is Φ = dim Σ.
- *The bounds.* Φ starts at 0. It ends at W, because the targets t_l sit at frames at distance k_l ≥ 1. It never
  exceeds W.
- *Each step.* A free step never raises Φ: mixing is frame by frame, and free letters permute frames and fix L_0. A
  paid block g on array R raises it by at most 1. The only form that can newly leave L_0 is R's form at L_0, call it
  α. It moves to gL_0.
- *Tightness forces growth.* With exactly W blocks, every block raises Φ by exactly 1 and no free step lowers it. So
  Σ_j = Σ_(j−1) + Qα_j with α_j outside Σ_(j−1).
- *After the last block.* The only form off L_0 that lies outside Σ_(W−1) is α_W at the frame g_W L_0, at distance
  r_W. Afterwards, free letters keep distances from L_0 and mixing works frame by frame. So at the end, every frame
  off L_0 at a distance other than r_W carries a form in Σ_(W−1).
- *Conclusion.* If no k_l equals r_W, every t_l lies in Σ_(W−1). That space has dimension W − 1, which contradicts
  the independence of the t_l. ∎

**Consequences.**
- *The survivors.* Every survivor listed in `out_enum_678_proj.txt` is frame-tight, with call ranks disjoint from task
  ranks. Examples: S7 has calls {3, 7} and tasks {2, 5, 6}; S6 has {3, 4, 6} against {2, 5}; S8 has {3, 4, 8} against
  {2, 5, 6}. So no survivor is an exact program without catalysts.
- *The NOTE's construction target is impossible.* NOTE section 5 ("No waste: it must be exactly tight on the frame
  count") and next-step 3 ("a tightness argument for S7") are answered in the negative for that exact form.
- *What it does not exclude.* The lemma removes nothing from the screen, because the duality works with:
  - *catalysts*, the same net vector with extra matching tasks and calls. S7 + [3] ≤ 3[3] + [7] passes the graded
    form of the argument: order the calls 3, 7, 3, 3 and give the targets 6, 5, 2, 3 the levels 1 to 4;
  - *slack*, (1 + ε) more calls;
  - *portions* of registers. Here the lemma fails outright. Put x and y in the two halves of one register, split by a
    spare row s. One rank-2 block on rows {1, s}, followed by a free mixing of the halves, gives H_1 x and H_1 y. That
    is the frame-tight J-relation [1] ≤ [2] (two tasks of weight 1/2 each), whose last block has rank 2, not 1.

  A construction search has to aim at those forms.

## 6. My predictions and how they fared

Written at 05:40 (the file says 05:44; the Mac clock read 05:40:58), before reading code, outputs or background proofs
(`STATUS.txt`).
- **Q1** (90 %): B1 exactly right. *Held.*
- **Q2** (75 %, 85 %): B2's closed forms right, and its conclusion. *Both held.*
- **Q3** (60 %): the walk condition's lower bound justified. *Held,* but by a step the NOTE does not write.
- **Q4** (85 %): calibration numbers reproduce. *Held* (53.1568 by rerun only).
- **Q5** (85 %): survivors' tightness and windows as stated. *Held.*
- **Q6** (75 %): my profiles agree with certlib. *Held:* 2,603 of 2,603.
- **Q7** (35 %): some coverage figure mismatches. *Held in kind.* The figures match the files, but one label is
  wrong (the ≤ 5-row claim) and the conclusion goes beyond the coverage.
- **Q8** (15 %): I kill one of S6 to S8 with a valid certificate. *Failed* for certificates. Lemma L kills their exact
  catalyst-free forms instead.
- **Q9** (70 %): not refuted on the exact core, caveats on coverage. *Held,* but the headline structural claim is
  refuted as stated.

## Sources

All paths are under `the archive's checks/wht43/` unless absolute. Everything outside `exist3/v-e3-boost/` was
only read.
- **Read:**
  - `exist3/BRIEF.md`, `exist2/NOTE.md`;
  - `exist2/ex2-duality/PROOF.md`, `exist2/v-ex2-duality/STATUS.txt`;
  - `the archive's checks/wht35/wider/w-bounds.md` §§1–3 and 5.4;
  - `moonshot/BRIEF2.md` (letters);
  - `meta/ex-laser/PROOF.md` (labels of Theorem 1 and section 1.3);
  - `wave16/b56-compute/out_pricecap.txt`;
  - all of `exist3/e3-boost/` named in its NOTE.
- **My files** (in `exist3/v-e3-boost/`):
  - `vB1.py` and `out_vB1.txt`: B1 and B2;
  - `vcert.py`: dual polar intersection numbers, distance-set profiles, brute force at m = 2 and 3;
  - `vlines.py`: the line-incidence family, a row-by-row count checked by brute force at m = 3;
  - `vsurv.py` and `out_vsurv.txt`: the survivors and the walk condition;
  - `vvalid.py`: realisable relations;
  - `vkill4.py` and `vcertlib_A4tight.py`: the block-4 A-tight relation;
  - `vbig.py`: hub roots and distance sets at block 73;
  - `vcmp.py`: profiles and projection counts against certlib;
  - `vlp.py`, `vlp3_rig.py`, `vlp3_big.py`, `vlp_big.py`: the LPs; outputs `out_*.txt`;
  - Colab runs: `cjob/remote-out/ssh5/remote.log` (survivors up to block 14; LP at registers 5 and 6) and
    `cjob2/remote-out/ssh5/remote.log` (the `lp_combo3` rerun).
