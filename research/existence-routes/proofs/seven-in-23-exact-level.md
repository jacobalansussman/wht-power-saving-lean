<!-- header: existence-routes -->
# The block-4 survivor 7 in 23 at the exact level

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** The computations reproduce (one AI skeptic), but the label "excluded" is refuted: random contents fail the same test in the same way, so only one generic lift is excluded. The f = 2 blindness holds (checked by computer, re-run 19 October 2026). A later wave closed this construction in both move orders in the configuration built (moonshot-block4-block5 section).

---

# e3-bpexact: the 7 in 23 survivor at the exact level

15 October 2026, 04:33 to about 07:10 EDT (Mac clock). I am e3-bpexact (exist3, wave 6). Nothing was posted. Labels:
[proved] (an argument written here), [computed] (a script and its output in this folder), [argued], [guess]. "Derived
here" means no skeptic has seen it. RUNLOG.txt is the working record; STATUS.txt holds the dated predictions and the short
result.

The object is w5-blueprint's survivor: 7 arrays in 23 paid moves at block 4 (saving 0.142), symmetric under Z7, with orbit
order O O F O F. Its middle shared move F1 reads one term of the second orbit through a combination Y of 184 letters, and
it passes the frame-sum relaxation (the cheap test that only looks at cells). This round asked for the exact screen.

## 0. In short

- **No unit, and the survivor as built is excluded at the exact level** [computed in floating point with a wide gap:
  junk 1e-14 or less, kernel singular values 1e-15 against an O(1) scale; the F2 step is argued].
  - The survivor lifts to actual 16 x 16 Clifford matrices (section 1).
  - Its distance-4 condition is now decomposed by character of Z7 and by irreducible representation of the letter group.
    In each piece it is a linear-algebra test (section 3).
  - The lifted survivor fails that test in six representations: the natural 16-dimensional one, its conjugate, and four
    "address-type" ones, in which the Paulis act through a commutative algebra (gl14, gl15, aff16, sq16). It fails in
    every character of every one of them.
  - The terminal shared move F2 cannot repair it. One shared move can supply at most d missing ranks in total, where d is
    the representation's dimension, and the survivor misses 86 to 112.
- **The brief's f = 2 test is blind, so I did not run it as an ALS** [computed]. With the move contents as unknowns, the
  f = 2 widened test accepts ONE array served by TWO moves, which is impossible at block 4 (section 2).
- **The failure in the natural representation is a property of the lift, not of the family** [computed]. Lift each
  content letter together with Pauli multiples (a "flat" lift). Then the top orbit alone reaches rank 16 - m_a in
  character a, where m_a is how often the a-th eigenvalue occurs in the order-7 twist p (4 for a = 0, 2 for the others).
  An F2 twisted by the same p would supply exactly those m_a ranks: the condition is met with equality in all seven
  characters, under the relaxation that F2's row is free.
- **A better order: O O F F O** [computed, derived here].
  - Here two group-algebra shared moves, with different twists and characters, are both read by the top.
  - It passes the frame-sum relaxation in all seven characters with no terminal move.
  - With flat lifts, the two moves together pass the exact distance-4 test in the natural representation: 16 of 16 in
    every character, each move covering the other's own character.
  - F O O F O, built once with a uniform start, fails in its own character.
- **Where it stops.** The address-type representations (here gl14, the 14-dimensional piece of the linear address
  action on the 15 nonzero addresses, on which the Paulis act trivially):
  - there the Pauli trick is unavailable: the Paulis do not span the matrix algebra;
  - with the letters as drawn, both families reach only rank 2 of 14;
  - with kappa free, the exact criterion (section 4.4) shows that NO lift gives the top alone the full kernel block in gl14
    (one constituent needs rank 2 and gets 1). ALS converges at rank 8 (numerically), so 8 <= r^max <= 13 at a = 1,
    where F2 twisted by p needs 12. That last gap is the open point.
  - These representations are frame sums over refined cells (a cell plus a flag) [argued]. The next exact screen is
    therefore purely symplectic.
- **Predictions** (STATUS.txt, written before computing): scored in section 6.

## 1. The lift [computed: lift.py, run1.py]

- **Configuration.** The same rule as w5-blueprint/s2/t10_build.py; the cells came out the same: O1 in cell 1, O2 in 16,
  F1 in 137, kernel 731.
- **Matrices.** Every symplectic object is replaced by a 16 x 16 matrix whose symplectic image is checked.
  - The Z7 lift c was rescaled by a power of i to have exact order 7.
  - q is a letter fixing F1's move cell, with exact order 7. This works because its symplectic order is odd and the
    kernel to Sp(8,2) is a 2-group.
  - From these, p_j = m q^j m^-1 and beta_h = q^h c_h^-1.
  - c_h, c'_h = K^-1 c_h^-1 K, q, p_j and beta_h are all Z4-monomial matrices: genuine letters.
- **F1.** F1 = m . sum_h psi(h)^-1 beta_h tau_h(sum_b y_b b O2) over 400 letters b of B_m (products of 151 sampled ones;
  their slot permutations span the full 184-dimensional algebra A).
  - It is equal to the closed formula (difference 0).
  - It is absorbed exactly: tau_g(F1) = psi(g) c_g p_g^-1 F1 as formal sums, with defect 1e-15.
  - All of its array-0 terms lie in cell F.
- **Y on the lift.** Y was re-solved on the lift (residual 2e-13), with a random null-space vector added so that the
  lift is generic.
- **Relaxation check.** relax.Eng.decide, run on the lifted terms with the top content free, passes all six characters
  (gap 0.198 vs 3e-14). The control Y = 1 fails.

## 2. The f = 2 widened test with the contents as unknowns is blind [computed: f2blind.py, out_f2blind.txt]

- **Why** [argued]:
  - g (x) g preserves Alt^2 V, so the antisymmetric block is passed by the input alone.
  - On Sym^2 V = D16 + S120 one move D(z) has off-diagonal blocks of rank 8, while K (x) K has rank 15 there.
  - With free widened contents, each move moves 8 arbitrary rows between the two blocks.
- **Construction.** One array, two moves D(e1), D(e2), both reading the input, with explicit block contents and readouts.
  - rep2-style residual 6e-30, against 0.496 with one of the moves dropped.
  - Fresh random inputs: error 1e-15.
- **Consequence.**
  - A one-array program needs four moves at block 4, so this test accepts an impossible program.
  - The survivor's top orbit has seven members reading the inputs, so it passes trivially; running the ALS would show
    nothing.
  - The same rank argument goes through for engine-letter contents: their algebra is still full on S120 and A120.
    [argued]

## 3. The exact screen: distance 4, character by character, irrep by irrep

**Derivation** [argued, derived here; every ingredient checked numerically].
- **What reaches distance 4.** The kernel term [K] x_0 sits at input distance 4. In O O F O F only two things put terms
  there:
  - the top orbit's part that reads F1 (distance 3);
  - F2.
  The inputs, O1 and O2 are at distance 0 to 2.
- **The top's F1 part.** Let A in C[B] be the top representative's coefficient on F1. Since
  tau_h(F1) = psi(h) c_h p_h^-1 F1, top member h contributes psi(h) c_h D3 A p_h^-1 F1.
- **Fourier over the arrays.** In each character a of Z7 (after the psi relabelling) the distance-4 part of output 0
  must satisfy
      X_a . D3 A P_a kappa f0 |_d4  +  (F2's part)  =  K,
  where:
  - P_a = sum_j a(j) p_j;
  - kappa = m (sum y_b b) m^-1;
  - f0 = m g2;
  - X_a in C[B] is the readout, free per character.
- **The slot module.** The elements involved live in sum_d C[B] g_d f0. Here g_d = D3 b_d runs over the 30 one-move
  cosets, and slot d is the cell of g_d f0, a neighbour of F. Right multiplication by a letter permutes the cosets with a
  cocycle; in an irrep sigma of B this is a block-monomial matrix R_sigma(k) of size 30 d_sigma.
- **Decoupling.** A's components A e_a in different characters are independent, so the characters decouple. By
  Wedderburn, the irreps decouple too.
- **The test.** In each sigma and a, with w = E_d1 R_sigma(A) (the d1-row block):
      w R_sigma(e_a kappa) Pi_x = 0   and   w R_sigma(e_a kappa) Pi_d0 invertible,
  where x is the set of slots nothing else can reach. Two versions:
  - 'all30': every slot except d0 and the 4 that O1/O2 reach;
  - 'd4': the 15 other distance-4 slots, which needs no reachability argument.
  Both gave the same ranks everywhere.
- **Linear algebra.** For fixed kappa the test is linear in w. I report r_a, the largest kernel-block rank a junk-free w
  reaches.
- **The trivial irrep.** It is the distance-4 part of the frame-sum relaxation.

**F2** [argued].
- An absorbed F2 has F2[g] = pi_g F2[0] c'_g for a homomorphism pi_g = beta'_g^-1 c_g. In character a its block is
  Y2 sigma(E_a(pi)) phi, with one readout Y2 and one row phi.
- Its rank is at most m_a(pi), the multiplicity of the a-eigenvalue of sigma(pi_1). These sum to d_sigma when pi_1 has
  exact order 7.
- A row-space count gives a necessary condition for character a: r_a + m_a(pi) >= d_sigma. This holds with the top's row
  block taken arbitrary, which is an upper bound.
- Summed over a, a necessary condition for every F2: **sum_a (d_sigma - r_a) <= d_sigma**.
- When the top's rows are free (sigma = V) and phi is treated as free, the condition is also sufficient.

**Calibration** [computed, run1b.py, flat.py].
- The trivial irrep reproduces the relaxation's verdicts on the lift:
  - the designed Y passes (r_a = 1 in all six characters);
  - a = 0 alone fails (r = 0), as in 7 in 22, which needs F2;
  - a random y on the same letters fails all six.
- A conjugation bug in the first null-space code (the transpose instead of the conjugate transpose) was found and fixed
  at 04:58, and its "passes" were voided. So was a relative rank tolerance that counted noise as rank (RUNLOG 05:13).
  Every number below comes from the fixed code, with an absolute threshold and the junk printed (1e-14 or less).
- Calibration list of the brief:
  - the plain schedule, E8 and the control word were not rerun: no system of this shape exists for them;
  - the block-2 stand-ins, engine 2 in 5 and Theorem 17 belong to other pipelines.
  - What this test was calibrated on is above: pass on the designed object, fail on the random one, fail on 7 in 22's
    own character. This is a declared limit of the calibration.

## 4. Results

### 4.1 The survivor as lifted is excluded [computed: fullrow.py, d4gen.py; F2 bound argued]

r_a for a = 0..6, generic lift, upper-bound version:

| irrep sigma | d | r_a | m_a(p) | total deficit | verdict with one F2 |
|---|---|---|---|---|---|
| trivial (frame sum) | 1 | 0 1 1 1 1 1 1 | 1 0 0 0 0 0 0 | 1 | passes (F2 serves a = 0) |
| V (natural) | 16 | 0 0 0 0 0 0 0 | 4 2 2 2 2 2 2 | 112 | excluded |
| conj V | 16 | 0 for a = 1..6 (run1b; a = 0 not tested) | as V | at least 96 | excluded |
| gl14 (address action, 14-dim) | 14 | 0 2 2 2 2 2 2 | 2 (all) | 86 | excluded |
| gl15 = 1 + gl14 | 15 | 0 3 3 3 3 3 3 | 3 2 2 2 2 2 2 | 87 | excluded |
| sq16 (b o b, the D16 block of V(x)V) | 16 | 0 (all) | 4 2 2 2 2 2 2 | 112 | excluded |
| aff16 (affine permutation) | 16 | 0 1 1 1 1 1 1 | 4 2 2 2 2 2 2 | 106 | excluded |

- The structured top (A in C[B], d4gen.py on Colab #2) gives the same r_a as the upper bound in gl14, gl15, sq16 and aff16.
- A second lift has the same table in V, gl14 and aff16 (seed 31, psi = 3, other letters and twist; out_run1_s31.txt,
  out_fullrow_s31.txt). r_1 = ... = r_6 always, as they must be: a field automorphism fixing Q(i) permutes the
  nontrivial characters [proved: all data except the characters lie over Q(i)].

### 4.2 In V the failure comes from the lift [computed: flat.py, kfree.py]

- **The tensor form.** V extends to the whole group, so R_V(k) = Delta (Pi_k (x) k) Delta^-1 (checked).
- **The flat lift.** Replacing each content letter k by sum_pi c_pi (pi k), with sum_pi c_pi pi = k^-1, gives
  R_V(kappa) = Delta (Y (x) I) Delta^-1. The identity was checked on the actual Pauli-multiplied terms (64 per letter).
  The relaxation can be kept [argued]: kappa's images in the trivial and in the Pauli-nontrivial irreps of B'_m are
  independent, and a correction by the Pauli average e_P (which vanishes in V) restores Y.
- **Result.** With the flat lift, r_a = 16 - m_a(p) exactly: 12 at a = 0 and 14 at a = 1..6.
- **Why** [argued]. In the eigenblock mu of V(p), the condition becomes the frame-level condition for the character
  a - mu. The block mu = a meets F1's own character, which fails structurally (dist 0.456 at the frame level).
  - Call this the **inheritance law**: in V every character inherits F1's own-character obstruction on an m_a-dimensional
    block.
  - An F2 twisted by p supplies exactly m_a there: the criterion holds with equality.
  - So V does not exclude the family. It does say how F2 has to be built: twisted by p (or by a pi with the same
    eigenvalue multiplicities), and filling the inherited blocks exactly.
- **Kappa free.** kfree.py: 9000 letters of B'_m, a generic row block, a = 1.
  - The top alone cannot reach the full kernel block in V: relative residual 0.153 (d4) and 0.161 (all30). The matrix
    has numerically zero singular values.
  - The trivial irrep is solvable on the same letters (5e-15).
  - This is consistent with the inheritance law holding for EVERY kappa [conjecture]. If it does, F2 is necessary in
    O O F O F, and it is tight.

### 4.3 Other orders with the same construction [computed: ooffo.py, ooffo2.py, foofo.py]

- **O O F F O** (7 + 7 + 1 + 1 + 7 = 23). F1 and F2 are group-algebra moves on the same O2 term, with twists from two
  different cyclic groups and psi1 = 1, psi2 = 3.
  - **Frame sum: passes in all 7 characters with no terminal move.** Each move alone fails exactly its own character.
  - **Exact distance 4 in V, flat lifts: 16 of 16 in every character from the two together** (14 or 12 from each alone).
    Generic lifts: 0.
  - Address-type representations, generic lifts: gl14 2 of 14 (the two moves give the SAME 2-dimensional kernel rows), sq16 0,
    aff16 1.
  - Repeated with another lift (seed 31, psi1 = 2, psi2 = 5; remote-out/ssh4/out_ooffo_s31_25.txt): the same verdicts.
  - So O O F F O improves on O O F O F wherever Paulis act (no F2 tightness needed), and it meets the same wall where they
    do not.
- **F O O F O.** A uniform psi-twisted start F1 feeds O1 (two routes); F2 is a group-algebra move on O2. A route from F1
  to a distance-3 cell next to the kernel was found.
  - The frame sum fails exactly chi = psi, with or without the F1 route: in this skeleton the start does not serve F2's
    own character.
  - Only one skeleton was built: not an exclusion of the order.
- **Not built:** O F O F O, O F F O O, the non-chain posets.

### 4.4 Kappa free in the address-type representations [computed: kmax.py, ubad.py, kals.py]

**The exact criterion with kappa free** [proved, derived here; short].
- The image of C[B'_m] in R_sigma is a semisimple algebra, a sum of full matrix algebras M_(d_tau) over the
  B'_m-constituents tau of R_sigma, where tau occurs m_tau times.
- The top's rows are w' = w R(e_a). Since p lies in B'_m, they lie in a submodule whose tau-part has dimension
  n_(tau,a) x m_tau, where n_(tau,a) is the multiplicity of the a-eigenvalue of p on tau.
- Write T for the target block [Z on d0 | 0 on the constrained slots | free on the reachable ones]. For one character,
  with the top's rows free, T = w' R(kappa) for some w' and kappa **iff rank(T_tau) <= n_(tau,a) for every tau**. Here T_tau is the tau-component of T's rows,
  read as a d_tau x (rows x m_tau) matrix.
- In the trivial irrep this is exactly the frame-level structural criterion ("e_K in A.S_a + F").
- So the family question in each irrep is a rank-constrained linear problem in (Z, free blocks). It is decidable once R_sigma
  is decomposed over B'_m. I did not finish that decomposition.

**What was computed with kappa free.**
- **Generic top rows, kappa free (kmax.py, Colab #3, 4000 letters with distinct symplectic images).** In gl14 at a = 1
  the reachable kernel rank is 6 of 14 (the designed lift gave 2).
  - This is a LOWER bound on r^max: the optimum must align the rows with the target, and a generic w does not.
  - F2 twisted by p would need 12.
- **The multiplicity-zero part of the criterion (ubad.py)** gives a kappa-INDEPENDENT upper bound: T must vanish on the
  isotypic parts with no a-eigenvector of p.

  | irrep | r^up at a = 0 (d4 / all30) | r^up at a != 0 |
  |---|---|---|
  | trivial | 0 / 0 | 1 (= the frame-level obstruction: calibration) |
  | gl14 | 12 / 10 | 14 |
  | gl15 | 12 / 10 | 15 |
  | aff16 | 15 / 15 | 16 |

  - Only F1's own character is bounded, with a deficit of at least m_0(p) (2 in gl14, 3 in gl15): the inheritance law
    again.
  - This is NOT an exclusion of the family.
- **ALS over (rows, kappa) (kals.py)** aims at special rows. Results for gl14 at a = 1:
  - rank-8 target: residual 0.0056 after 40 iterations and still falling geometrically. That is numerical evidence that
    rank 8 is reachable, not a certificate;
  - rank-4 target: likewise reachable (0.0096 at iteration 20, and falling);
  - rank-12 target (what F2 twisted by p needs): 0.0185 at iteration 25, falling about as fast as the rank-8 run did. It
    is still running on Colab #5 (<colab work folder>/remote/e3-bpexact/out_kals_gl14_r12.txt, 45 iterations): no verdict.
  - The V calibrations (rank 14 must be reachable, rank 16 not) were lost when Colabs #1-#3 went down (05:59). One was
    restarted and is unfinished.
- **The exact criterion, decided in one case (crit.py, Colab #4).** gl14, a = 1, d4 mode, kappa free, rows free.
  - The commutant of R_gl14(B'_m) has dimension 79, found by propagation over the 13 orbitals of slot pairs. It matches
    the decomposition: 31 irreducible copies in 17 types, sum d m = 420, sum m^2 = 79.
  - With no help from the free slots, the full kernel block violates the criterion in 5 types.
  - In ONE type (dimension 8, multiplicity 1, n_(tau,1) = 1) the target needs rank 2 whatever the free slots add.
  - **So for every kappa and every top content, the top alone cannot give the full kernel block in gl14** [computed, exact
    criterion; floating point with clear rank gaps]: r_1^max <= 13, and 8 <= r_1^max on the ALS evidence.
  - At a = 0 (crit, the same decomposition) no type is decisive by this crude free-slot bound. The sharper ubad bound
    gives r_0 <= 12 there.
  - F2 twisted by p (m_1 = 2) could still cover a deficit of 1 or 2. So gl14 does not exclude the family at this
    resolution. It does confirm, kappa-independently, that the top alone never suffices: the O O F O F family needs F2 in
    gl14 as it does in V.

## 5. What the exact screen says about the family, and what to do next

1. **Finish the gl14 decision; it is one rank-constrained problem.** Is there a rank-12 Z, with free-slot blocks F, such
   that rank(T_tau(Z, F)) <= n_(tau,a) in all 17 types? crit.py already gives every T_tau; what remains is the minimum-rank
   completion over F. Then:
   - add the coupling that kappa is ONE element for all seven characters;
   - run the same check for O O F F O, whose two twists may cover each other's deficits here as they do in V.
2. **The wall: the address-type representations.**
   - **What they are.** gl14, gl15, aff16 and sq16 factor through the address action. For any rep induced from a
     subgroup B' of B, the sigma-test is the frame-sum relaxation with cells B' g in place of B g [argued].
     - gl15 is induced from the stabiliser of a hyperplane of L_Z, equivalently of a pair of neighbours of L_Z.
     - Its test is a **flagged frame sum**: cells are (Lagrangian, 3-subspace of it). There are 34425 of them, and the
       kernel cell has 15 flags, of which 14 must cancel.
   - **What it needs.**
     - It is purely symplectic: it needs Q = the symplectic image of B'_m (order 688128 = 64 x 10752, a 64-fold
       refinement of the slot group G), not Clifford matrices.
     - It is the natural next exact screen for O O F O F (with F2 twisted by p) and for O O F F O.
     - kappa ranges over C[Q], so the screen is decidable as a finite bilinear problem.
3. **If both families fail it.** The group-algebra construction is closed at 7 in 23. The obstruction would read: one O2
   term cannot carry enough flag information through one shared move.
4. **If one passes.**
   - Build kappa flat in every Pauli-nontrivial irrep, using the Paulis, and right in the flagged frame sums, using the
     order-64 kernel N.
   - Then run the distance-3 and lower conditions, which involve the free orbits.
   - Then the finite-field pipeline on the whole program, and raw arrays at f = 1..6.
   - No unit candidate reached the brief's step 4 (Singular, raw arrays), so those steps were not run.
5. **Housekeeping for later screens.**
   - The f = 2 widened test must never be used with free contents.
   - The per-irrep distance-4 test (d4test.py, d4gen.py, fullrow.py) applies to any program whose kernel-reaching moves
     read one shared move. Its trivial-irrep part reproduces the relaxation.

## 6. Predictions (STATUS.txt, 04:33, before computing), scored

- **P1** (calibration passes, 90 %): partly testable, and passed where tested.
  - The trivial-irrep calibration passed.
  - The brief's list (plain schedule, E8, stand-ins, control word, engine 2 in 5, absorb4, Theorem 17) was not run with
    these tools: no system of this shape exists for them.
- **P2** (the lift exists, 85 %): true.
- **P3** (f = 2 ALS reaches 0 with probability 12 %; stalls at a positive floor 80 %): moot. The test is blind and passes
  trivially.
- **P4** (the obstruction shows in the per-irrep conditions, likely one character pair, 60 %): true in a stronger form.
  Every character fails in six irreps for the lift as built, and the inheritance law explains the shape in V.
- **P5** (exact existence, 8 %): no.
- **P6** (another open order passes the frame-sum relaxation, 30 %): true. O O F F O passes in all 7 characters. That
  none of them gives an exact unit (4 %) is not decided.
- **P7** (near-miss or excluded with a named per-irrep obstruction): true. The lift as built is excluded; the family is
  open at the address-type representations.

## 7. Files (this folder)

- **Construction and screens:** lift.py (the lift), d4test.py (per-irrep distance-4 test, full rows), d4gen.py
  (structured top; gl14, gl15, sq16, aff16), fullrow.py (upper bound and F2 count), flat.py (flat lift; F2 criterion),
  kfree.py and kmax.py (kappa free), ubad.py (kappa-independent bound), kals.py (ALS over rows and kappa), crit.py
  (the exact kappa-free criterion via the commutant).
- **Calibration and orders:** run1.py and run1b.py (relaxation on the lift; all-30-slot test and controls),
  f2blind.py (blindness of f = 2), ooffo.py and ooffo2.py (O O F F O), foofo.py (F O O F O).
- **Outputs:** out_*.txt. Colab outputs: remote-out/ (fetched copies). ys_seed21_psi1.npy holds the designed content
  coefficients.
- **Library:** w5lib/ (unchanged copies of w5-blueprint's sp8.py and relax.py and of vlib/, for the Colab runs).
- **Runner:** py.sh (the Mac runner through a memory-capped runner).
