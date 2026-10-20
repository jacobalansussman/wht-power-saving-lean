<!-- header: existence-routes -->
# The skeptic's record: seven in 23 exact level

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** This is an AI skeptic's record of its check, as it wrote it. Its verdict is the status this revision uses for the track it checked, unless a later verdict is named in the section README.

---

# v-e3-bpexact: skeptic's check of e3-bpexact (the 7 in 23 survivor at the exact level)

15 October 2026, 06:49 to 07:23 EDT. I am v-e3-bpexact. Nothing was posted. Predictions are in STATUS.txt (written
before computing). RUNLOG.txt has every run. Labels: [computed] (a script and its output here), [argued], [their code] (a
rerun of e3-bpexact's own script, which checks that the result reproduces but is not independent).

## In short

- **Every computation I rebuilt reproduces.**
  - The lift is legitimate and absorbed exactly.
  - The generic lift fails the distance-4 test in V, gl15 and aff16 in every character, with the same deficits.
  - The flat lift gives exactly 16 - m_a in V.
  - The trivial-irrep calibration and the f = 2 blindness reproduce.
  - O O F F O passes the frame sum in all 7 characters and V (flat) 16 of 16.
  - The exact kappa-free criterion reproduces on gl14, and on a second lift.
- **But the outcome label "excluded" does not hold for the survivor.** Only one generic lift is excluded. Two of my results
  show that this exclusion says nothing about the w5-blueprint survivor:
  - With a random y on the same letters, the V test fails in exactly the same way (rank 0 in every character).
  - The rank count explains it [argued]. R(P_a kappa) has rank at most about 480/7 (P_a is a projector of a Z7 action), so
    any kappa with no special structure gives rank 0 against the 240 constraint columns. That is a property of generic
    contents, not of the survivor.
  - e3-bpexact's own flat-lift result keeps the survivor alive in V, and so does my rebuild of it: with an F2 twisted by p
    the criterion holds with equality. The question is still open at gl14.
  - So the correct outcome is **open (near-miss)**: generic lifts are excluded, and the flat-lift family is open at the
    address-type representations.
- **One new exact fact.** I ran the commutant criterion (their crit.py) in V. For every kappa, the top alone gives at most
  rank 14 at a = 1 in V, and the flat lift reaches it. So their conjectured "inheritance law" holds at a = 1 as an exact
  upper bound (in the full-row relaxation), and F2 is necessary in O O F O F.

## 1. The derivation, checked

I re-derived the distance-4 condition myself.
- Top member h is psi(h) c_h D3 A p_h^-1 F1 + ...
- F1's terms on array g are psi(g)^-1 p_g F1[0] c'_g.
- The psi factors combine to psi(j)^-1, so after a Fourier transform over the arrays the readout is X(chi) and the effective
  character is a = chi - psi. This matches d4test.py.
- For O O F F O the readout X(chi) is common to both moves, with a_i = chi - psi_i. That is what ooffo.py stacks.

The F2 bound (one absorbed F2 supplies at most m_a(pi) rank in character a, and these sum to at most d) is also right. If
pi is only a homomorphism on F2[0] (pi_(g+h) F2[0] = pi_g pi_h F2[0]), then sigma(pi_g) Phi = T^g Phi with T^7 = 1 on the
cyclic span. The character blocks are projections of one d-dimensional space, so their ranks still sum to at most d.

Every exclusion here is made in the 'd4' mode, which constrains only the 15 other distance-4 slots and so needs no
reachability argument.

## 2. What I rebuilt with my own code

- **Legitimacy (vd4V.py).** Each of these is Z4-monomial, Clifford (sigma asserts it) and fixes L_Z:
  - c_h, c'_h, p_j, beta_h and q;
  - the 2800 coefficient letters q^h b c_h^-1 that F1 actually applies to O2_h;
  - the 400 letters k_b and the products c_h c'_h.

  0 bad. m = D_5 x letter. q^7 = c^7 = I, and p is a homomorphism.
- **Absorption.** With my own formal sums, tau_g(F1) = psi(g) c_g p_g^-1 F1 holds with defect 7.9e-17 (seed 21) and 8.7e-16
  (seed 31). For F2 of O O F F O the defect is 1.0e-15.
- **The V test without the slot module.** In V the image of a cell's formal sum is just the actual matrix sum of its terms. So
  I built the test from actual products g_d p_j k_b f0, using my own coset representatives. The slots come from my own
  symplectic composition, spot-checked against the cells of the product matrices.

  | lift | r_a (a = 0..6), need 16 | m_a(p) | verdict |
  |---|---|---|---|
  | generic, seed 21 (D3 = D_1 and D_3) | 0 0 0 0 0 0 0 | 4 2 2 2 2 2 2 | deficit 112, excluded with any F2 |
  | generic, seed 31 / psi 3 | 0 0 0 0 0 0 0 | same | same |
  | random y (control) | 0 0 0 0 0 0 0 | same | same |
  | flat, seeds 21 and 31 | 12 14 14 14 14 14 14 | same | r_a + m_a = 16 everywhere |

  The singular-value gaps are clear: flat 3.3e-2 or 4.7e-2 against 5e-16; generic, top singular value 1e-15 on a scale of 0.9.
- **Address-type representations (vd4addr.py).** sigma(beta) comes from the address permutation of beta = M rep_c^-1.
  - gl15: r = 0, 3, 3, 3, 3, 3, 3 (m = 3, 2, ...; deficit 87).
  - aff16: r = 0, 1, ..., 1 (deficit 106).

  Both are e3's numbers. gl14 = gl15 - trivial is consistent with their 0, 2, ..., 2.
- **Calibration, trivial irrep.**
  - The designed y gives r = 0 at a = 0 and 1 at a = 1..6.
  - A random y gives 0 everywhere.

  The same holds on both seeds. This is e3's calibration, reproduced.
- **Positive control (new).** The plain chain (kappa = 1, no twist) reaches full rank at a = 0 in the trivial irrep, V, gl15
  and aff16. So the test does not reject a working chain. This is a control e3-bpexact could have run from the brief's
  list: the plain schedule does have "a system of this shape".
- **f = 2 blindness (vf2blind.py).** With my own D16 / S120 / A120 projectors:
  - every content and readout of their 1-array, 2-move program is block-diagonal, so it lies in the image of C[B'];
  - the identity K(x)K = R_x + sum_i R_i (D_i (x) D_i) A_i holds to 6.5e-16.

  The test is blind, as claimed. Their argument that the engine letters' image at f = 2 is also full M16 + M120 + M120 checks
  out [argued]: the diagonal phases i^(a.y + 2Q(y)) separate all 120 pairs, and AGL(4,2) is 2-transitive.
- **O O F F O (voofo.py, designs taken from my rerun of their ooffo.py).**
  - Frame level, d4 part: F1 alone fails chi 1, F2 alone fails chi 3, both together pass all 7.
  - V generic: 0.
  - V flat: 14/12 alone and 16/16 together in every chi.

  Their relax.Eng.decide verdict (all 7 pass; each move alone fails only its own character) reproduced [their code].

## 3. The exact kappa-free criterion (crit.py) [logic checked; their code]

- **The logic is right.**
  - The image of C[B'_m] is a sum of independent full matrix algebras.
  - The rows w R(e_a) lie in the image of e_a, of dimension n_(tau,a) m_tau. So T_tau = K_tau W_tau with col(W_tau) inside
    an n_(tau,a)-dimensional space, and the condition is iff.
  - The "lower bound with F" is valid: projecting out every free-slot coordinate can only lower rank(T_tau + A).
  - The commutant is computed from 13 generators and verified on 4 other letters, and its dimensions are consistent.
- **gl14, seed 21** (their output): a decisive type (8, 1, 1) needs rank 2 and gets 1, so r_1 <= 13.
- **V, seed 21** (my rerun): 5 types; type (96, m 1, n 14) needs rank 16 after any free-slot help. So r_1 <= 14 for every
  kappa, and the flat lift attains it.
- **gl14, seed 31 / psi 3.** The first attempt crashed with a LAPACK "SVD did not converge" (a numerical fragility of
  crit.py). I relaunched it with a gesvd fallback (critrun2/, Colab #5, 231 s). It gives the same 17 types and the same
  decisive type (8, 1, 1), which needs rank 2, so r_1 <= 13 reproduces on a second lift.

## 4. What I did not check

- gl14, sq16 and conj V directly (gl14 is consistent through gl15 - trivial).
- ubad.py's bounds.
- The ALS runs (numerical and labelled as such).
- F O O F O (one skeleton; not claimed as an exclusion).
- The Colab-only kfree and kmax runs.

## 5. Small inaccuracies in the claim

- **conj V.** The result says conj V is "0 in every character", but a = 0 was not tested there (their own NOTE table says so).
- **Calibration.** "No system of this shape" for the brief's calibration list is too strong. The plain schedule fits the d4
  test as a positive control (section 2).
- **The label.** "outcome: excluded" contradicts "The family is not excluded" in the same claim (see "In short").
