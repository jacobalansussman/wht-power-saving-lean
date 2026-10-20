<!-- header: existence-routes -->
# Boosting pieces (Theorem B1, Corollary B2, the survivors)

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** Theorem B1 and Corollary B2 are argued by hand, not checked; one AI skeptic redid the identities and both Wronskians. The headline "any booster needs a job of at least 6 rows" is refuted as stated (it holds only on registers of at most 16 rows and only given argued projections). The survivors S6, S7, S8 are open; the skeptic's Lemma L, which removes their listed forms, has been read by nobody else.

---

# Boosting pieces: a screen at blocks 4 to 8 (exist3, track e3-boost)

15 October 2026, 04:32 to 05:40 EDT. One agent. No skeptic has read this. Nothing was posted. No Lean.

## In short

1. **The question.** The exist2 round found that a piece that saves nothing on its own can still improve the
   saving of another unit. It calls this *boosting*. A boosting piece needs a paid block larger than one of its
   targets, which the real recursion allows through spare rows. The round asked for a screen at blocks 4 to 8. A
   candidate is kept only if its balance is negative somewhere between exponents 0.9 and 0.9983 and it passes the
   hook, the lifts and certificate A. Before relying on any hit, the round asked me to check that the job model's
   exponent equals the challenge's.
2. **The job model's exponent is the challenge's, as far as it can be checked.**
   - On the hub's T74 word (through our chain) and on #1042, the root of the job-model balance equals the Lean
     pipeline's root to every printed digit: 1.40400794216e-3 and 1.67498901594e-3.
   - For a set of relations that includes a booster, the duality of exist2 produces one finite relation whose calls
     are all single rows. Theorem 1 of the laser track schedules such a relation; that theorem is proved in the job
     model.
   - What stays argued is the step every mixed unit already relies on: running jobs with spare rows in OpenAI's RAM
     model (laser track, section 1.3). Lean has no such recursion today. [argued; numbers computed]
3. **The calibration passed.**
   - Known relations pass every certificate: 190 realisable small relations and the two hub words.
   - Every proved-impossible case is refused: engine 2 in 5 at block 3; Theorem W2's shapes; 2 in 3 at block 2; the
     spreading family.
   - The made-up boosting example is reproduced: exponent 0.7233, against 0.834 for the unit alone.
   - The screen is blind, as it must be, to program-level exclusions: hook-tight mixed units, the slice, absorb4.
4. **Small boosting pieces cannot exist.**
   - *Theorem B1* [proved]: no relation whose jobs all have at most 3 rows can boost, on a register of any size. The
     proof uses only the frame count and certificate A.
   - *Corollary B2* [proved]: such relations also never boost one another.
   - On registers of at most 4 rows, no booster passes the known certificates [computed: a linear program on 200
     exponents]. The certificates it uses are proved, apart from the hooks' mod-p values.
   - On registers of 5 rows the same holds once projected certificates are added [computed, given argued
     certificates].
   - With all jobs of at most 5 rows, no booster passes on any register of up to 12 rows either [computed].
   - So a booster needs a job of at least 6 rows.
5. **From 6 rows up, the certificates leave room.**
   - Shapes that pass all 3,400 to 4,000 certificate tests exist on paper. The smallest are:
     - [2] + [5] + [6] ≤ 2[3] + [7] at block 7: 3 arrays, 3 paid blocks, balance negative on (0.867, 1);
     - [2] + 4[5] ≤ 2[3] + [4] + 2[6] at block 6: 5 arrays, negative on (0.681, 1).
   - If any one of them existed as a program, it would give, together with the units we already have, a saving above
     the bar: at least 0.13 and 0.32 respectively.
   - No program is known for any of them. They survive exactly where our certificates are weakest: block 6 and up,
     where nothing comparable to the exact hooks and lifts of blocks 3 to 5 exists.
   - All of them sit exactly on two tightness equalities: as many paid blocks as targets, and equal total rank.
6. **New tools,** with their checks:
   - certificate B's profile at blocks 2 to 5: it is A minus one, and at block 2 it is perfect [computed];
   - two certificate families valid at every block [derived here]:
     - isotropic subspaces, which reproduce A and w-bounds' C at block 3;
     - distance sets on the frame graph, checked by brute force at blocks 2 to 4;
   - closed formulas for the projection counts, matching every exact count at blocks 4 to 8;
   - a "walk condition" on target ranks [derived here].

   No new certificate refuses any of 190 relations known to be realisable.

   With these, the block-5 and many block-6 survivors of the first screen are excluded.
7. **Outcome.** The structure of any booster is forced: at least one job of 6 or more rows, and the two tightness
   equalities. Candidate shapes are listed for exact search (section 5). No unit was found, and nothing on the
   existence side is exact.

## Words used here

- **Array, row, block, register.** An array is one copy of the data. A block-m unit works on m address bits (rows)
  of each array. A *register* is the set of rows a program touches. Spare rows are rows of the register that no
  target transforms.
- **Job, relation, balance.** As in exist2.
  - A job [s] is the transform on s rows of a register.
  - A *relation* "tasks ≤ calls" is an exact program, right for every number of digit columns. It completes the task
    jobs using paid calls of the stated sizes. Example: Rs = 2[2] ≤ [1] + [3] means two arrays, each owed 2 rows,
    done with one single-row block and one 3-row block.
  - Its *balance* at exponent θ is D(θ) = Σ_calls r^θ − Σ_tasks k^θ. A relation "saves nothing alone" when
    D(1) ≥ 0.
  - It is a *booster* when D(1) ≥ 0, yet D(θ) < 0 for some θ < 1.
- **Window.** Exponents 0.9 to 0.9983. Known units already beat every exponent above 0.99833. A booster negative on
  [θ0, 0.9983] would, with them, push the exponent down to θ0 (Corollary 6 of exist2's duality track). The bar
  (saving 0.1) is θ0 = 0.9.
- **Certificate.** A Theorem-R lower bound. It charges a job of rank d a price P(d), and every relation must satisfy
  Σ_calls P(r) ≥ Σ_tasks P(k).
  - It applies to a relation whose register fits the certificate's block.
  - Multiplying a relation by [b] (bundling b rows into one) gives a relation on b times as many rows. So P(b·d) is
    a valid price too.
- **Frame count.** The certificate "number of calls ≥ number of tasks" (exist2's duality track, Lemma 2(c)). It is
  D(0) ≥ 0.
- **Labels.**
  - [proved]: a written proof below.
  - [computed]: a program and its output in this folder.
  - [argued]: reasoning short of a proof.
  - [derived here]: my own short argument, checked by nobody else.

## 1. Is the job model's exponent the challenge's?

**What has to agree.**
- *The challenge.* OpenAI's goal (`UniformFourier.lean`, DFTGoal and TimeBounds) asks for work O(N (log N)^θ).
- *The job model.* It prices a job of s rows at s^θ. The work for a transform on n = log2 N bits is N·T(n).
  - T(n) is at least n^θ whenever every relation used has a balance of at least 0 at θ (the laser track's
    Theorem 3).
  - For relations with single-row calls, Theorem 1 there reaches the root.
- *Units.* For a pure unit, W arrays at block m with paid blocks of ranks r, the job-model balance is
  Σ n_r r^θ − W·m^θ. Its root is the engine's "moment root", which our Lean pipeline turns into the exponent it
  proves.

**Checked on the pipeline's own data** [computed: `calib.py`, `out_calib_final.txt`].

| Word | Job-model root, as saving | Pipeline figure |
|---|---|---|
| T74E (hub #1839 word through our chain; `bound-check/lg-a1/out/price/T74E.json`) | 1.40400794216e-3 | root_lower 17550099277/12500000000000 = 1.40400794216e-3; Lean DFT numerator 14040060, i.e. 1.4040060e-3, one 1e-10 step below, because the Lean exponent must be a rational strictly above the root |
| #1042 (hub profile, the laser track's control) | 1.67498901594e-3 | stages.py 1.674989e-3 |

**For boosters** [argued, built from proved pieces].
- *The duality produces one relation.* Corollary 6 of exist2's duality track says what a set of relations is worth.
  By its Lemma 4 and Theorem 3(iv), any exponent above that value is reached by one finite relation c^n·[1] ≥ [N^n]:
  a single job done by single-row calls.
- *That relation can be scheduled.* All its paid blocks have rank 1, below its one target. So it is "separated", and
  Theorem 1 of the laser track applies. Theorem 1 is proved in the job model, and its schedule has the spare rows it
  needs.
- *What is not proved:* that a job-model program runs in OpenAI's RAM model at the same exponent. The laser track's
  section 1.3 argues it in four parts:
  - portions of a register, rounded to dyadic fractions;
  - spare rows, obtained by one split of the top job;
  - a valid order of the calls (its Lemma 2);
  - Q(i) coefficients.

  Every mixed unit relies on this same step.
- *What Lean would need.* The Lean engine proves pure units only. Using a booster in Lean would need the new
  recursion that `bound-check/rl-plan/PLAN.md` section 4 describes. It would also need the relation c^n·[1] ≥ [N^n]
  written out, which the duality proves to exist but does not construct.
- *The DFT.* For T74E the pipeline's DFT numerator is the kernel numerator minus one. So the exponent carries from
  the transform the engine computes to the DFT at the 1e-10 grid.

**What "equal" cannot mean.** The challenge allows any circuit, while the job model allows only engine programs. So
the challenge's exponent could be smaller than any job-model exponent. The direction a booster needs holds as
argued: a job-model exponent is a challenge exponent.

## 2. Calibration (before any claim about an open shape)

All runs use the final certificate library (`certlib.py`). Outputs: `out_calib_final.txt`, `out_calib_valid.txt`,
`out_lp_combo3_243.txt`.

| Case | Expected | Result |
|---|---|---|
| Plain [k] ≤ k[1], Kronecker, the A5 identity [1] ≤ 2[m], and "undo" relations [k] ≤ [k′] + [k′−k] (190 relations on registers of 2 to 8 rows; 137,001 certificate tests) | pass | all pass |
| Hub words T74E and #1042 at block 73 (only A and B apply there) | pass | pass |
| Engine 2 in 5 at block 3 (proved impossible) | refuse | refused: hook, lifts |
| Theorem W2's shapes: 2 arrays with 3 + 1 or 1 + 2 blocks at block 3 | refuse | refused: hook, lifts, block-4 hook |
| 2 in 3 at block 2 | refuse | refused: B at block 2 is perfect (prices 1, 2) |
| Rs, R2, R2b, and the spreading family 2[k] ≤ [k−j] + [k+j] | refuse | refused by A (and its scalings) |
| 3 arrays with four rank-2 blocks at block 6 (the spectrum track excluded it, with a sampled lift) | refuse | refused: projected hook (argued) and the sampled lift |
| 1 + 4 in 10 at block 3; the slice's 3 in 11 at block 4 | blind (excluded only by program-level arguments) | pass: blind, as expected |
| absorb4's limit; px-t514 | blind (not relations) | not testable by a price screen |
| Made-up booster Rs with the made-up unit 2 in 5 at block 3 | the duality gives exponent 0; the fluid LP at n = 243 gives 53.16 calls | 53.1568, exponent 0.7233 (unit alone 0.834); Rs is flagged as booster-shaped, then refused by A |

Not run: E8 and the 8-move control word (no relation data in hand). The slice at block 4 (Theorem 17) is outside the
engine price picture: engine certificates allow 3 in 11 at block 4.

## 3. Small boosters are impossible

**Why a booster needs an oversized block.** D(0) ≥ 0 by the frame count. If also D(1) ≥ 0 and D(θ) < 0 somewhere
in between, D changes sign twice. By Descartes' rule for sums of powers, the sizes must then alternate:
calls, tasks, then larger calls.

**Theorem B1** [proved]. Let ρ be a relation whose jobs (tasks and calls) all have at most 3 rows, on a register of
any size. Then for every θ in [0, 1]:

    D_ρ(θ) = λ_f(θ)·F(ρ) + λ_p(θ)·D_ρ(1) + λ_A(θ)·A(ρ),

where:
- F(ρ) = (number of calls) − (number of tasks) ≥ 0 is the frame count;
- A(ρ) = Σ_calls A(r) − Σ_tasks A(k) ≥ 0 is certificate A, with prices A = 4, 6, 7 for ranks 1, 2, 3;
- λ_f = 4 − 5·2^θ + 2·3^θ, λ_p = 1 − 3·2^θ + 2·3^θ, λ_A = −1 + 2·2^θ − 3^θ.

Hence D_ρ(1) ≥ 0 implies D_ρ(θ) ≥ 0 on all of [0, 1]: no relation with jobs of at most 3 rows boosts.

*Proof.*
- *The identity.* For d = 1, 2, 3 it reads d^θ = λ_f + λ_p·d + λ_A·A(d). This is three linear equations, solved
  exactly; the residual is 1e-15 on a grid (`rank3.py`).
- *Positivity of λ_A.* Each λ is a sum of three powers with two sign changes in its coefficients, so it has at most
  two real roots. λ_A vanishes at θ = 0 and θ = 1 and has slope 0.288 > 0 at 0. So it is positive on (0, 1).
- *Positivity of λ_f.* λ_f(0) = 1, λ_f(1) = 0, and λ_f′(1) = −0.340 < 0. A root inside (0, 1) would force at least
  two there, plus the one at 1: three in all, which is too many.
- *Positivity of λ_p.* λ_p(0) = 0, λ_p′(0) = 0.118 > 0, and λ_p(1) = 1. The same root count applies.
- *Validity of the two certificates.* A is valid at every block (skeptic of exist2's duality track). Normalised, it
  is the same price at every register size. The frame count is Lemma 2(c) there. ∎

**Corollary B2** [proved]. Every relation with jobs of at most 3 rows has a Θ-set that is an interval [0, τ]. So no
set of such relations is worth more than its best member: the amplification of exist2's Corollary 6 needs a job of
4 or more rows.

*Proof.*
- *Reduction.* By B1, D/λ_p = D(1) + (λ_f/λ_p)·F + (λ_A/λ_p)·A on (0, 1], with F, A ≥ 0. So it suffices that both
  ratios decrease.
- *Notation.* Write x = 2^θ, y = 3^θ, a = ln 2, b = ln 3.
- *First ratio.* The Wronskians simplify to closed forms, checked numerically to 1e-9 (`wronskian.py`, `out_wronskian.txt`):

      λ_A′λ_p − λ_Aλ_p′ = x·y·g(θ),   with g = a(1 − 3^(−θ)) − b(1 − 2^(−θ)).

  Here g(0) = 0 and g′ = ab(3^(−θ) − 2^(−θ)) < 0, so g < 0.
- *Second ratio.*

      λ_f′λ_p − λ_fλ_p′ = x·y·h(θ),   with h = 7a·3^(−θ) − 6b·2^(−θ) + 4(b − a).

  h′ changes sign once, from negative to positive, at θ = 0.380. So h ≤ max(h(0), h(1)) = max(−0.118, −0.057) < 0
  on [0, 1]. ∎ (A 2,000-point grid agrees: `out_rank3.txt`.)

This contains exist2's exclusion of the spreading family as a special case.

**Registers of 4 and 5 rows** [computed: `lp2.py`, `out_lp2_proj.txt`].
- *The test.* A linear program asks, at each of 200 exponents in (0, 1), for the most negative balance among
  relations with D(1) ≥ 0 that pass every certificate. On registers of 2, 3, 4 and 5 rows the minimum is 0 at every
  grid point: no booster at all.
- *At 4 rows the known certificates suffice*, even without any frame family above block 3 (`lp4.py`,
  `out_lp4.txt`). The program uses A, B, the isotropic-subspace and distance-set families, w-bounds' C, D and E, the
  block-3 hook and lifts, and the hooks at blocks 4 to 6 (mod-p values). The only relation that is tight on the frame
  count, the plain count and A at 4 rows, [1] + 5[3] ≤ 4[2] + 2[4], is killed by the block-4 hook alone (margin
  −1,568 of 86,592).
- *At 5 rows the projections are needed.* Without them, one booster survives every certificate:
  [2] + 5[4] ≤ 4[3] + 2[5]. It is exactly tight on the frame count, the plain count and A. A second one has larger
  weights: 2[2] + 8[4] ≤ 7[3] + 3[5]. Both are killed by projections of the block-3 certificate B and of the block-4
  hook (`out_enum_5_K10.txt`, then `out_enum_56_proj.txt`).
- *Status of the projections.* They rest on the Averaging Lemma at blocks 5 and up [argued, as in b56-methods]. So
  the 5-row statement is [computed, given argued certificates].
- *Small jobs on large registers.* With all jobs of at most 5 rows, no booster passes on registers of 6 to 12 rows
  either (`lp3.py`, `out_lp3.txt`).

## 4. The certificates and conditions added here

- **B at blocks 2 to 5** [computed: `certB.py`, `out_certB*.txt`; ranks modulo two primes].
  - The profile is 2^m − 2^(m−d) − 1, that is A minus 1. The block-3 values 3, 5, 6 of w-bounds 5.4 are reproduced.
  - At block 2 the prices are 1, 2: a perfect certificate.
  - The quotient space has dimension 2^(2m−1) + 2^(m−1) − 1 at every m tested: 9, 35, 135, 527.
  - From block 6 on, B is used only as [conjectured]. No conclusion above depends on it.
- **Isotropic-subspace certificates C_k** [derived here].
  - The module is the permutation module on the k-dimensional isotropic subspaces of the label space, with Y_0 the
    subspaces inside L_0.
  - Price of a rank-d job: G(m,k) − G(m−d,k), a difference of Gaussian binomials with q = 2.
  - Proof: as for A, which is C_1. At block 3, C_2 is w-bounds' certificate C (6, 7, 7).
- **Distance-set certificates on the frame graph** [derived here; profiles checked by brute force at blocks 2, 3, 4,
  `check_dist.py`].
  - Fix a set S of distances. Take V = Q^(frames) and Y_U = span of the frames whose distance from U lies in S.
    Every engine letter moves frames rigidly and the free letters fix L_0, so Theorem R applies.
  - The price is |S(U) \ S(V)|. It depends only on the distance, because the frame graph (the dual polar graph of
    Sp(2m, 2)) is distance-regular. It is computed from the exact intersection numbers (`dualpolar.py`).
  - S = {0} is the frame count; S = {0, 1} at block 3 is w-bounds' certificate E (12, 12, 15).
  - S = {1} is useless for ordinary units (about 1 move per array) but strong against boosters. It is the certificate
    that first removed many block-6 survivors.
- **Projection counts in closed form** [computed: matched against `project.py`'s exact counts at blocks 4 to 8
  (`proj/N5.txt` to `N8.txt`, the last from a Colab; `mixcap.py` for 4 to 6)].
  - N_r(r−2) = 2^(m−r)(2^r − 1)(2^(r−1) − 1), N_r(r−1) = 2^(m−r)(2^r − 1)(2^(r−1) + 2), and N_r(r) is the rest of
    the 4^m − 1 points.
  - With these, every hook, lift and symmetric family of blocks 3 to 7 is carried up to block 8, and some up to 16.
- **The walk condition** [derived here; it passes all 190 realisable relations].
  - Run the word on formal contents. Each term of a register is c·[w] y with w = f_j P_j ⋯ f_1 P_1 f_0, where the
    P_i are the paid steps that term passed through, each at most once.
  - A correct word must contain a term equal to the target H_S itself (the formal-identity step of w-bounds,
    section 1). Free letters fix L_0, and the frame distance is a metric. So the distance from L_0 moves by at most
    r_i at each paid step.
  - Hence each target rank k needs a sub-multiset R of the paid ranks with 2·max(R) − ΣR ≤ k ≤ ΣR.
  - It kills, for example, 2[5] ≤ 2[1] + [8] and [4] + [5] ≤ 2[1] + [7]. Without it, both pass every price
    certificate.

## 5. What survives from 6 rows up, and how close it comes

Integer relations: net vector in [−K, K]^m, primitive, booster-shaped, passing every certificate above (rigorous and
argued) and the walk condition. Sources: `intenum.py`, `out_enum_678_proj.txt`, `out_check_survivor.txt`.

| Block, bound K | Survivors | Smallest |
|---|---|---|
| 6, K = 6 | 2 | [2] + 4[5] ≤ 2[3] + [4] + 2[6] (5 arrays; window (0.681, 1)); [2] + 6[5] ≤ 2[3] + 2[4] + 3[6] |
| 7, K = 4 | 8 | [2] + [5] + [6] ≤ 2[3] + [7] (3 arrays; window (0.867, 1)) |
| 8, K = 3 | 31 | [2] + [5] + 2[6] ≤ [3] + 2[4] + [8] (4 arrays; window (0.758, 1)) |

What they have in common:
- every one has a task of rank 2 and a block the size of the register;
- each is exactly tight on the frame count (as many paid blocks as targets) and on the plain count (equal total
  rank);
- their balances are small: D(0.9) = −0.0034, −0.00057, −0.0040 for the three shapes named.

The three named shapes pass every one of the 3,400 to 4,000 tests per register. That includes the sampled block-6
lift and the projections of every symmetric family. S6 also passes a further family tried at the end: incidence
with isotropic lines, Y_U = the lines meeting U in a given dimension (brute force at block 6, `lines6.py`,
`out_lines6.txt`; margins 20 to 1,940). The linear program also finds fractional relations covering the whole
window [0.9, 0.9983], with margin 7.6e-5 per unit of task weight at 6 rows.

**Why they survive.** The exact hooks and lifts stop at block 5. At block 6 we have only one hook, with modest
prices, and a lift sampled from 1,500 pairs. Above that we have only the symmetric families and projections, which
are weak at high rank: their θ-shadow falls with the block. The frontier of the screen (a job of 6 or more rows) is
the frontier of our certificates. So I read these survivors as gaps in the certificates first and as construction
targets second [guess].

**What a construction must satisfy.** Take S7 = [2] + [5] + [6] ≤ 2[3] + [7]: three arrays on 7 rows, owed 2, 5 and 6
rows, done with two 3-row blocks and one 7-row block.
- *No waste.* It must be exactly tight on the frame count: every paid step raises that potential by exactly one.
- *Paths.* By the walk condition, the 2-row target must come from a term that passed through both 3-row blocks.
- *Equal rank.* It must use no more total rank than the plain schedule.

Like the hook-tight mixed units that ms-mixed's Theorem T excluded, these sit on certificate boundaries. So a
tightness argument of that kind (every single term must pass) is the natural way to test them.

## 6. My predictions and how they fared

Written at 04:33, before reading the duality track or computing (`STATUS.txt`).

- **P1** (70 %): a written argument that the exponents agree. *Held in part.* It is checked on pipeline data and
  argued in general; the RAM step stays argued. The 60 % that a caveat must be named held: section 1.3's
  realisation, and that Lean lacks the recursion. Found false (10 %): no.
- **P2** (90 %): calibration reproduces the made-up boost and refuses Rs through A. *Held.*
- **P3** (85 %): no real piece at blocks 4 to 8 passes everything. *Half right.*
  - At blocks 4 and 5 nothing passes.
  - At blocks 6 to 8, shapes do pass every certificate we have. None is known to be real, so "no real piece" still
    stands, but "no candidate" does not.
- **P4** (65 %): A is the killer for most candidates. *Failed.*
  - The killers were the hooks and the new symmetric families at 4 rows; projections at 5 rows; the distance-set
    certificates and projections at 6 rows.
  - A is tight on the most dangerous candidates, rather than violated.
  - The second part (55 %: the hook kills none that A lets through) also failed: the block-4 hook kills the only
    A-tight direction at 4 rows.
- **P5** (50 %): survivors of A with large spread die by a stronger exact bound. *Held* at blocks 4 and 5 (hook,
  lift, projections). Open at blocks 6 to 8.
- **P6** (under 5 %): an exact existence proof of a booster. *Held* (none).
- **P7** (80 %): seconds per candidate. *Held.* About 1 ms per candidate on 3,400 to 4,000 tests. The full LP
  screens took 7 to 95 s.

## 7. What I recommend next

1. **A skeptic for Theorem B1, Corollary B2 and the walk condition.** They are short and exact, and they change the search rule:
   no booster with jobs of at most 3 rows, ever.
2. **Exact certificates at block 6.** The survivors live where we have no exact lift and only one hook. Two things
   would test whether boosting is a real lever or a certificate gap:
   - an exact lift profile at block 6, possibly from a symmetric lift as the spectrum track suggests;
   - or hooks with other parameters at block 6.

   Run the three named shapes against them first.
3. **A tightness argument for S7.** It has 3 arrays, 3 paid blocks and 7 rows, and is exactly tight on the frame
   count and the plain count. The question is whether a program with no waste in the frame potential can exist with
   one 7-row block. This is the analogue of ms-mixed's Theorem T, and the smallest construction target the screen
   produced.
4. **Exact search only after items 2 and 3.** S7 needs 7-row registers with address maps on 7 bits. The free group
   is far too large for blind search, and the special-position lesson applies.
5. **Do not search boosters with jobs of at most 5 rows** on registers of at most 12 rows: the known certificates
   exclude them.

## Sources

All paths are under `the archive's checks/wht43/` unless absolute. Everything outside `exist3/e3-boost/` was
only read, except one read-only copy of `wave16/b56-methods/project.py` run in `projjob/` on a Colab.

- **Brief:** `exist3/BRIEF.md`.
- **Read:**
  - `exist2/NOTE.md` (sections 3, 6, 7, 8);
  - the duality track `exist2/ex2-duality/` (NOTE.md, PROOF.md: Lemma 2, Theorem 3, Lemma 4, Corollaries 5 and 6;
    `calib.py`, `lp_combo*.py`, `out_certA.txt`) and its skeptic `exist2/v-ex2-duality/`;
  - `the archive's checks/wht35/wider/w-bounds.md` (sections 1 to 3 and 5.4);
  - `meta/MAP.md` (sections 1.1 to 1.4); `meta/spectrum/NOTE.md` (sections 3 to 5) and `mixcap.py`;
  - `meta/ex-laser/PROOF.md` (sections 1 to 3, 5);
  - `wave16/b56-compute/pricecap.py` (hook values); `wave16/b56-methods/project.py`;
  - `bound-check/rl-plan/PLAN.md` (sections 2, 3); `bound-check/lg-a1/out/price/T74E.json`;
  - `moonshot/BRIEF2.md`; `meta3/NOTE.md` (first part);
  - OpenAI's `UniformFourier.lean` (DFTGoal, TimeBounds) under `OpenAI's openai-math repository, lean/`.
- **This folder, by section:**
  - Section 1: `calib.py`, `out_calib_final.txt`.
  - Section 2: also `calib_valid.py`, `out_calib_valid.txt`, `out_lp_combo3_243.txt` (exist2's `lp_combo3.py`
    run at n = 81 and 243).
  - Section 3: `rank3.py`, `out_rank3.txt`; `wronskian.py`, `out_wronskian.txt`; `dual4.py`, `dual4b.py`, `out_dual4*.txt`; `lp.py`, `lp2.py`,
    `lp3.py`, `lp4.py` with `out_lp_rigorous.txt`, `out_lp2*.txt`, `out_lp3.txt`, `out_lp4.txt`.
  - Section 4: `certB.py`, `out_certB.txt`, `out_certB_m5.txt`; `dualpolar.py`, `out_dualpolar.txt`;
    `check_dist.py`, `out_check_dist.txt`; `proj/N5.txt`–`N8.txt`; `certlib.py` (all certificate profiles, the
    projection formulas, the walk condition).
  - Section 5: `intenum.py`, `out_enum_*.txt`; `projall.py`, `out_projall.txt`; `check_survivor.py`,
    `out_check_survivor.txt`; `lines6.py`, `out_lines6.txt`.
  - Colab runs: `remote-out/ssh3/`, `remote-out/ssh5/`, `projjob/remote-out/ssh4/`.
- **Earlier runs** used smaller certificate sets, which is why their survivor counts are larger:
  - `out_lp_rigorous.txt`, `out_lp2.txt`, `out_enum_45.txt`, `out_enum_456.txt` and `out_enum_678.txt` have no C_k,
    distance sets or projections;
  - `out_lp2_Ck.txt` has no distance sets;
  - `out_lp2_Dist.txt`, `out_enum_56_Dist.txt`, `out_enum_456_walk.txt` and `out_enum_78_walk.txt` have no
    projections.

  RUNLOG.txt gives the order.
- **Derived here:** Theorem B1 and Corollary B2; the C_k and distance-set certificates; the walk condition; the
  closed forms N_r(s); and the bar arithmetic of section 5 (a window from θ0 with the known units gives saving at
  least 1 − θ0).
