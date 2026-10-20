> **Working note, kept as written.** An AI agent wrote this during the logic round of 14-15 October 2026 (the ml-fooling track). Private paths and machine names
> were removed and nothing else was changed. No person has checked it. Where it disagrees with this section's
> [README](../README.md), the README wins, because it follows the later skeptic verdicts. Its skeptic's verdict is [two-certificates-skeptic.md](two-certificates-skeptic.md): the mathematics holds, and four
> readings are refuted ("a proof that an existence proof exists is never cheaper", "counting is the only
> non-constructive route", "counting can prove only border units", "mod p is never a proof over C").

# Existence by fooling every refutation (ml-fooling, meta3)

14 October 2026, two sessions: 20:08 to about 21:00, and 23:01 to about 00:10 EDT by the Mac clock (the first session stopped before finishing its write-up; the second added sections 2.5 and 3.1, finished section 5, and checked everything new against known cases). Nothing was posted. No skeptic has read this note.

Labels: [proved] means a proof is written here (sections 2 and 3). [computed] means a program and its output in this folder. [read: first-hand] means I read the statement in the paper itself. [memory] means I did not re-read it today. [argued] means reasoning short of proof standard. [guess] is a guess. "mod p" means computed modulo the prime 2147483629 (section 2.4 says what that does and does not prove).

## 0. In short

The question for this track: can we prove a unit exists by showing that every possible refutation of it fails?

1. **Yes in principle, and exactly** [proved, section 2.1].
   - For one fixed shape, "a unit exists" is the same statement as "no refutation of degree up to D_eff exists". A refutation is an algebraic identity 1 = Σ g_i f_i built from the shape's unit equations f_i. D_eff comes from the effective Nullstellensatz (Kollár 1988).
   - The failure of all refutations up to a degree is a linear-algebra fact. Its witness is a "pseudo-solution": a linear functional that is 0 on everything such a refutation could use and 1 on the constant 1.
2. **The worst-case degree is hopeless for us**, by hundreds of orders of magnitude [computed, section 2.2].
   - Smallest real system: a 37-unknown skeleton from w3-asym. D_eff is about 10^17.7, so the linear algebra would have about 10^610 unknowns.
   - The exact equations of the open block-4 shape 2 in 7: 10^(8.7·10^23) unknowns.
   - Our machines handle about 10^7.
   - Sparsity, symmetry and multihomogeneous bounds lower the exponent but stay hopeless (10^610 becomes 10^357 at best).
3. **The usable form is a "flat" pseudo-solution** [proved, section 2.3].
   - Flat means a moment matrix stops growing in rank. A flat pseudo-solution proves existence at any degree, however small. This rests on Laurent and Mourrain 2009, read first-hand, plus a sharper degree condition proved here.
   - It works best inside the Polynomial Calculus: close the degree-D space of consequences under multiplication by variables before testing. This is the proof system that Gröbner-basis computations follow.
   - In that form it certified the cubic tensor test systems at degree 4 instead of 7. The degree matches Singular's Gröbner basis, and the number of points matches as well [computed].
4. **Calibration passed both ways wherever the method could run** [computed, section 4]. It certified every known solvable system I could set up. The certificates were computed mod p and cross-checked exactly or numerically; the second session turned them into exact proofs (item 8):
   - the plain controls;
   - Karatsuba, the smallest saving unit for polynomial multiplication (our problem reshaped, iso/NOTE.md);
   - the W-state tensor at rank 3;
   - 57 random systems with a planted point.

   It refused every known-impossible one:
   - 58 Singular-inconsistent skeletons from w3-asym, refuted exactly over Q;
   - the W-state and Karatsuba at rank 2;
   - 73 random impossible systems.

   It never certified a false case. It is blind, for size reasons, on the block-2 stand-in 2 in 3 (1,672 unknowns) and on block-3 2 in 5.
5. **What it is, honestly.** A flat certificate is the same kind of computation as a Gröbner or border basis on the same fixed system. It is not a new source of existence proofs. Its merits:
   - the "yes" certificate is short to check;
   - the solutions (the unit's numbers) can be read off it;
   - it needs no Singular.
6. **What decides the degree: behaviour at infinity, but only for the weakest proof system** [proved one way, section 3].
   - Suppose an impossible shape has approximate solutions running off to infinity, with residual falling like (size)^(−L). Then every single-identity (Nullstellensatz) refutation has degree above L. So approximate units force long identities.
   - The Masser-Philippon chain, written in the same two-group form as our skeleton systems, needs identities of degree 4, 7, 12 at 3, 5, 7 unknowns.
   - Polynomial Calculus refutes the same chains at degree 3, 3, 4, 4, 5. So the "infinity" obstruction belongs to the identity format, not to existence certification as such.
7. **By-product on the search side** [computed, section 5].
   - All 58 saved impossible skeletons of w3-asym have short exact refutations: degree 2 or 3 as identities, degree 2 in Polynomial Calculus. The 74 impossible a2-nag systems need degree 2 to 5.
   - Re-running w3-asym's exact sweep with the same seeds reproduces its verdicts. It also decides 45 of the 62 skeletons of the (2,2) sweep that w3-asym had left undecided, all impossible. Item 10 gives the final count, with the sigma sweep.
   - ex2-modp independently built exact refutation identities for nag-rerun's skeletons. This is the same kind of object for a different family.
8. **A short exact proof format for "a unit exists"** [proved and computed, section 2.5; second session].
   - For one fixed shape, exactly one of two finite objects exists. One is a refutation (an identity 1 = Σ g_i f_i). The other is a set of commuting rational matrices at which every unit equation vanishes, a "finite model". Its existence is equivalent to a complex unit, and checking it is plain rational matrix arithmetic.
   - Flat pseudo-solutions are exactly these models in disguise. So "fooling every refutation" in the strong (flat) sense means exhibiting a unit over a number field, written as matrices.
   - Built and checked exactly over Q (no Gröbner basis or field guess is trusted):
     - Karatsuba (4 × 4 matrices, or 26 × 26 under a generic slice);
     - the W-state at rank 3;
     - nag-rerun's special-position control, where it finds a4 = 3 with no hints, as against 0 of 2,000 for random search;
     - one of nag-rerun's relaxation hits whose solutions need i (4 × 4 matrices, without being told the field);
     - random systems with up to 9 irrational solutions.
   - It refused every impossible case, including a decoy that has a clean solution modulo one prime and none over C.
   - Six certificates are on file with a standalone checker (`check_cert.py`).
   - It did not reach the second relaxation hit (25 unknowns), Toom-3 (40) or Strassen (59). There the Gröbner search modulo p is the limit, not the certificate.
9. **Existence by counting is blocked by solutions at infinity** [computed and argued, section 3.1].
   - The one classical way to prove existence without exhibiting a solution is to count solutions (Bézout's theorem). Counting needs the system to have no solutions "at infinity", and that requirement is itself a refutation. So it is the literal form of "existence proved by an impossibility proof".
   - For the Karatsuba test system, the count is 1,728 but only 4 solutions are finite. The other solutions sit at infinity, where the degenerate and border decompositions live. The route would prove border units, not units.
10. **The skeleton sweep is finished** [computed, section 5]. Of the 133 skeletons w3-asym left undecided, 76 are now refuted at degree 2 or 3: 45 + 28 at closure cap 150, and 3 more at cap 600, one of them with 1,156 unknowns. The other 57 have closures over 600 columns. None fooled degree 3.

**Bottom line.** Fooling all refutations characterises existence exactly but is not a cheaper route to it. Its practical form is Gröbner-strength linear algebra on one fixed small shape, which is what ex-note already concluded. The route proves no new unit. It leaves five things:
- an exact, checkable certificate format for "yes" (commuting rational matrices), with a working pipeline from modular search to exact check, ready to be stated in Lean;
- a fast, Singular-free exact refuter for skeletons;
- a lemma tying approximate units to long refutations;
- the observation that the only non-constructive existence route (counting) needs a refutation at infinity, which our systems fail because of degenerate solutions there;
- a precise list of what a structural theorem would need (section 7).

**For the question that started meta3** ("can one prove that a proof of the kind we want exists?"):
- On the existence side the answer is yes, in a strong and concrete form. For one shape, a unit exists exactly when a finite model exists, and no refutation exists exactly when a unit does. Both objects have a priori size bounds (products of the equation degrees). Neither side can be independent of ordinary arithmetic.
- But a proof that an existence proof exists is never cheaper than the existence proof itself. The best such proof is a flat pseudo-solution, and that is the unit, written as matrices.

## 1. The idea and the words

- **Unit system.** Fix a shape: the wiring, the letters, and which coefficients are unknown. "This shape is a unit" becomes polynomial equations f_1 = … = f_k = 0 in n unknowns over the complex numbers. In w3-asym's form the unknowns are program coefficients α and readout coefficients β, and the equations read A(α)β = e. "Correct for every f" is built in through the closures.
- **Nullstellensatz (NS) refutation of degree D.** An identity 1 = Σ g_i f_i with every deg(g_i f_i) ≤ D. Finding one is linear algebra on the space M_D spanned by the multiples m·f_i of degree ≤ D (m a monomial).
- **Polynomial Calculus (PC) refutation of degree D.** A derivation of 1 by adding lines and multiplying a line by a variable, with every line of degree ≤ D. Its degree-D space is M_D closed under "multiply the degree ≤ D−1 part by a variable". This closure is what a Gröbner computation truncated at degree D produces (Lazard). PC degree ≤ NS degree.
- **Pseudo-solution of degree D.** This is the "design" of proof complexity: a linear functional L on polynomials of degree ≤ D, equal to 1 at 1 and 0 on the degree-D space (NS or PC). A true solution x gives one: L(p) = p(x).
- **Flat at s.** The moment matrix H_s = (L(m·m′)) over monomials of degree ≤ s has the same rank as H_{s−1}.
- **The metamathematical reading.**
  - The Nullstellensatz is the completeness theorem of these little proof systems: if no refutation exists at any degree, a solution exists. Kollár's bound makes it a bounded completeness theorem.
  - Low-degree pseudo-solutions of impossible systems are their nonstandard models.
  - Flatness is the test that a nonstandard model is standard, that is, a finite combination of genuine points.

## 2. The statements

### 2.1 Pseudo-solutions are exactly the failure of refutations [proved; standard linear algebra]

A degree-D refutation fails to exist iff 1 is not in the degree-D space, iff some linear functional is 0 on that space and 1 at 1 (a subspace and its annihilator). Every solution gives such a functional at every degree.

So if D ≥ D_eff, a degree-D pseudo-solution proves that a common complex zero exists. Here D_eff is a degree that every refutation of an impossible system of this size is guaranteed to fit in.

### 2.2 The worst-case degree, and why it is hopeless here

**Bounds** [read: first-hand].
- **Kollár 1988** (J. AMS 1, 963–975), Theorem 1.5 and Corollary 1.9. Suppose f_1, …, f_k have no common zero, with degrees d_1 ≥ … ≥ d_k, all different from 2. Then an NS refutation exists with deg(g_i f_i) ≤ d_1⋯d_k if k ≤ n, and ≤ d_1⋯d_{n−1}·d_k if k > n. In all cases the bound is ≤ max(3, d)^n (as quoted by Sombra).
- **Brownawell 1987**, as quoted in Kollár's Theorem 1.1: deg g_i ≤ n·min(k, n)·d^min(k,n) + min(k, n)·d.
- **Sombra**, alg-geom/9610006: deg(g_i f_i) ≤ min(n, s)^2 (d + 3n) δ, where δ is the "geometric degree" of the system. In general δ ≤ (d + 1)^(min(s,n) − 1). The bound is polynomial when δ is small.
- **D'Andrea, Krick and Sombra**, arXiv:1103.4561, Theorem 0.1: the same products of degrees, plus bounds on the size of the coefficients.

**Sizes of our systems** (`deff.py`, `out_deff.txt`) [computed]. The columns are: unknowns n, equations k, degree d, log10 D_eff by Kollár, log10 of the number of monomials of degree ≤ D_eff (the size of the linear algebra), and the degree actually needed.

| system | n | k | d | log10 D_eff | log10 size | degree actually needed |
|---|---|---|---|---|---|---|
| w3-asym skeleton sk711_0_cumulative_4 (impossible) | 37 | 160 | 2 | 17.7 | 610 | NS 3, PC 2 |
| w3-asym skeleton sig2_901_1_cum2 (impossible) | 64 | 542 | 4 | 38.5 | 2,377 | NS 3, PC 2 |
| W-state, rank 2 (impossible; border rank 2) | 12 | 8 | 3 | 3.8 | 37 | NS 9, PC 5 |
| Karatsuba, rank 2 (impossible) | 14 | 12 | 3 | 5.7 | 69 | NS 9, PC 5 |
| Karatsuba, rank 3, gauge-fixed, 3 unknowns specialised (solvable) | 12 | 12 | 3 | 5.7 | 60 | flat at NS 7, PC 4 |
| block-2 VV stand-in, 2 in 3 (solvable) | 1,672 | 1,024 | 4 | 617 | 1.0·10^6 | not attempted |
| exact equations, block 2 plain 1 in 2 | 18,432 | 46,080 | 3 | 8,794 | 1.6·10^8 | — |
| exact equations, block 3 mixed 1+4 in 10 | 4.7·10^8 | 9.3·10^9 | 4 | 2.8·10^8 | 1.3·10^17 | — |
| exact equations, block 4, 2 in 7 | 1.1·10^12 | 1.9·10^14 | 5 | 7.8·10^11 | 8.7·10^23 | — |

The exact-equation sizes are ex-dimension's counts times |B|, the order of the free-letter group. The VV size is from `out_indep.txt`. In the last four rows the size column is itself an exponent: for 2 in 7 the linear algebra has 10^(8.7·10^23) unknowns.

**Structure lowers this, but not enough.**
- *Sombra's geometric degree.* For the 37-unknown skeleton, a two-group Bézout count C(37, 5) bounds δ [derived here, argued]. That gives D ≤ 6.7·10^10 and size 10^357 instead of 10^610.
- *Symmetry.* Averaging a pseudo-solution over a group G that permutes the equations gives an invariant one (Reynolds). It divides the size by at most |G|, which is 10^3 to 10^10 for our free-letter groups.
- *Eliminating the linear readout.* Solving A(α)β = e for β leaves the 1 to 21 program unknowns, but turns the equations into minors of degree up to about 33. With 5 unknowns Kollár gives 33^5 ≈ 4·10^7 and size about 10^36 [derived here].

**Verdict.** A pseudo-solution at the worst-case degree is a complete existence test that can never be run on a system of ours.

### 2.3 Flat pseudo-solutions prove existence at any degree [proved]

**Theorem** [read: first-hand]. Laurent and Mourrain, "A generalized flat extension theorem for moment matrices", Arch. Math. 93 (2009), arXiv:0812.2563, Theorem 1.5. It holds over any field and needs no positivity.
- Let C be a set of monomials closed under division, and C⁺ = C ∪ x_1C ∪ … ∪ x_nC.
- Let L be defined on products of pairs from C⁺.
- If the moment matrix on C⁺ has the same rank r as on C, then L extends uniquely to a functional L̃ on all polynomials whose moment matrix still has rank r.
- The kernel of that moment matrix is then an ideal whose quotient has dimension r.

**Corollary (existence certificate)** [proved here]. Let L be a pseudo-solution of degree D ≥ 2s, NS or PC. Suppose L is flat at s (take C = monomials of degree ≤ s − 1), and every deg f_i ≤ s + 1. Then the f_i have a common zero. rank H_s = r bounds the number of zeros found, counted with multiplicity.

*Proof.*
1. Let J be the kernel of L̃'s moment matrix. J is an ideal, and L̃(1) = 1, so J is not the whole ring.
2. The rank is r both on C and on all polynomials. So every polynomial q is congruent modulo J to some q_C of degree ≤ s − 1.
3. For each i: L̃(f_i q) = L̃(f_i q_C), because f_i(q − q_C) ∈ J.
4. L̃(f_i q_C) = L(f_i q_C), because deg(f_i q_C) ≤ (s + 1) + (s − 1) = 2s and L̃ agrees with L there.
5. L(f_i q_C) = 0, because f_i q_C lies in M_D, which is inside the PC space.
6. So every f_i is in J, and the ideal (f_1, …, f_k) is contained in J, which is not the whole ring.
7. By the weak Nullstellensatz the f_i have a common zero over the algebraic closure. Over Q that means over C.
∎

**Remarks.**
- Laurent and Mourrain's applications ask for s ≥ deg f_i. The weaker s ≥ deg f_i − 1 is what lets cubic systems certify at s = 2.
- The condition cannot be dropped. Karatsuba at rank 2 is impossible, yet at degree 8 its pseudo-solution had rank H_1 = rank H_0 = 1. That is flat at s = 1, below deg f − 1 = 2 [computed].
- The test was flat at s = 2, D = 4 for the sliced Karatsuba and W-state systems once the PC closure was used. Without the closure it needed D = 7 [computed].

**Specialising is allowed for existence only.** A certificate may fix a gauge, set unknowns to random numbers or add random linear equations, because a zero of the smaller system is a zero of the original. A refutation may do none of this. Positive-dimensional unit varieties (gauge freedom, free coefficients) need such slicing before a random pseudo-solution can be flat [computed: unsliced or under-sliced systems are never flat].

**Same strength as Gröbner** [computed on the calibration systems, argued in general]. Singular finds degrevlex Gröbner bases of the sliced Karatsuba and W-state systems with top degree 2 and exactly 4 and 2 solutions (`gbcmp/out_singular_gb.txt`). The PC-flat test certifies at D = 4 = 2·2 with ranks 4 and 2. Newton's method from 400 random complex starts finds exactly 4 and 2 distinct nonsingular solutions (`out_crosscheck.txt`). Flatness at s is a border basis at degree s, which is the object a Gröbner or border-basis algorithm computes. So the certificate shows up roughly where Gröbner finishes, not earlier.

### 2.4 What "mod p" proves

The flat-extension theorem is field-free. Run over Q, a flat certificate is a full proof over C. Run modulo p, it proves a zero over the algebraic closure of GF(p). That is not a proof over C: ex2-modp's Proposition 1 gives systems with a clean zero mod p and none over C. A rigorous version over C needs exact arithmetic.
- My exact refutations over Q (the 58 skeletons, the chains) are fine.
- My exact flat test with plain fractions ran over ten minutes without finishing on the 8-unknown W-state system, and I stopped it.
- A multimodular version with rational reconstruction is the fix [guess: 100 to 200 lines].

The calibration's mod-p certificates are cross-checked independently: Singular's exact solution counts, Newton's method, and for the random planted systems the recovered point checked exactly over Q.

Section 2.5 closes this gap with a different exact format, which avoids exact Gröbner or flat computations over Q altogether.

### 2.5 An exact certificate for "yes": commuting rational matrices (second session)

**Theorem (two finite certificates, exactly one exists)** [proved here; each ingredient is textbook algebra]. Let f_1, …, f_k be polynomials with rational coefficients in n unknowns. Exactly one of the following exists:
- **(R) a refutation:** polynomials g_i with rational coefficients and Σ g_i f_i = 1;
- **(M) a finite model:** for some size r ≥ 1, rational r × r matrices M_1, …, M_n that commute with each other and satisfy f_j(M_1, …, M_n) = 0 (the zero matrix) for every j.

Moreover, (M) exists exactly when the f_j have a common complex zero.

*Proof.*
1. *(M) gives a zero.* Commuting complex matrices have a common eigenvector v. Write M_i v = λ_i v. Then 0 = f_j(M)v = f_j(λ)v, so λ = (λ_1, …, λ_n) is a common zero.
2. *A zero gives (M).* If a zero exists, the ideal generated by the f_j in Q[x] is proper, so it lies in a maximal ideal m. By Zariski's lemma, K = Q[x]/m is a finite field extension of Q. Let r = [K : Q]. Multiplication by x_i on K, written in a basis of K over Q, is a rational r × r matrix M_i. These matrices commute, and f_j(M) is multiplication by f_j(x̄) = 0.
3. *Not both.* Substituting the M_i into Σ g_i f_i = 1 gives 0 = the identity matrix.
4. *At least one.* If there is no zero, the weak Nullstellensatz gives (R).
∎

**What it means.**
- *Checking a "yes" is trivial.* Checking (M) needs only exact rational matrix products: no Gröbner basis, no flat-extension theorem, no choice of number field. A unit whose numbers live in Q(i), or in a field of degree 26, is certified in the same format. This is the format to state in Lean. The one lemma it needs is the common eigenvector of commuting matrices; the rest is computation with rational numbers.
- *It is what a flat pseudo-solution really is.*
  - A flat L of rank r yields (M): the matrices of multiplication by x_i on the quotient, M_i = G^(−1) H^(i), with G = (L(b·b′)) and H^(i) = (L(b·x_i·b′)) over a monomial basis b, b′.
  - Conversely, (M) yields L(p) = (1/r)·trace p(M), which is a pseudo-solution at every degree.
  - In the words of section 1: the flat (standard) pseudo-solutions are the finite models. A pseudo-solution that never becomes flat corresponds to no finite model. Below its refutation degree, an impossible system has only pseudo-solutions of that kind.
- *Both sides have a priori bounds, both products of the degrees.*
  - (R) can be found with degree ≤ d_1⋯d_{n−1}·d_k (Kollár, section 2.2).
  - (M) can be found with r ≤ d_1⋯d_k [argued; the inequality is from memory]:
    1. Cut the solution set by a generic rational affine space whose codimension is the largest dimension of a component. The cut is finite and nonempty.
    2. By Heintz's Bézout inequality (1983), the cut has at most d_1⋯d_k points.
    3. The cut is defined over Q, so it has a point of degree at most that number over Q. Its residue field gives (M).
  - So for one fixed shape, the proof that a unit exists and the proof that none exists are both finite objects with a priori size bounds. There is no case where neither exists. This is the concrete form of the decidability that ml-logic states.
- *It gives the uniform checker that a Σ1 statement needs* [argued]. "Some shape has a unit with saving above 0.1" becomes: there are a shape and rational commuting matrices that satisfy that shape's equations. Every part of this is checkable by a fixed program, so the sentence is Σ1 (ml-logic's question (b)). For one fixed shape both answers have finite witnesses, (M) or (R). So neither the existence nor the non-existence of a unit of one fixed shape can be independent of ordinary arithmetic (Peano arithmetic, or much less).

**How the matrices are found** (`commcert.py`, `commcert_sing.py`, `commcert_ab.py`). Only the last step is trusted.
1. Search modulo several 31-bit primes: either a PC-flat pseudo-solution (Python, as in section 2.3), or a Gröbner basis and the normal forms of x_i·b (Singular).
2. Read off the multiplication matrices in a monomial basis. The basis must be the same at every prime; a prime whose basis differs is skipped.
3. Chinese remaindering and rational reconstruction of every entry. Stop when two successive reconstructions agree.
4. Exact check of (M) over Q with Python fractions.

A wrong step 1 to 3 can make the method fail. It cannot make it lie, because step 4 is the whole proof.

**Certificates on file.** `certs/*.json` hold six certificates:
- `K3_coord` (4 × 4)
- `W3_coord` (2 × 2)
- `Fp2b_ab_seed2` and `Fp2b_abscan_seed1` (4 × 4 each, the relaxation hit whose solutions need i)
- `spu_ab_seed2` and `spu_abscan_seed1` (1 × 1, the special-position control)

Each file holds the sliced system, the matrices, the original system, and the affine substitution that turns the original into the slice.

`check_cert.py` is a standalone checker that imports nothing from this folder. It checks the matrices, as above, and also that the slice really is the original system under that substitution. So a zero of the slice gives a zero of the original.
- It reports all six VALID.
- It reports INVALID when one entry of the Karatsuba matrices is changed by 10^−9.
- `add_provenance.py` rebuilt every slice from the original system with the same code and seeds, and all six matched.

A skeptic needs only `check_cert.py`, the certificate files and the theorem above. One assumption remains for the two nag-rerun systems: their "original" is nag-rerun's export lifted to small integers from one prime and checked against the second prime's export. ex2-modp made the same assumption.

**Slicing.** Positive-dimensional solution sets are cut to finitely many points first. This is allowed for existence only (section 2.3). There are three ways to do it, with different trade-offs [computed]:

| slicing | matrices | special position |
|---|---|---|
| set some unknowns to random integers | small: Karatsuba r = 4, 20-bit entries | missed: nag-rerun's control (solvable only on a4 = 3) is refuted, like ex2-modp's random search (0 of 2,000) |
| generic hyperplanes in all unknowns (h = dimension) | large: r is the degree of the solution set (Karatsuba 26, 207-bit entries; W-state 48, 796-bit entries, 54 primes); the slice makes every equation dense: for g27_Fp2b_m3_s2 (15 unknowns after slicing, degree 5) Singular did not finish the first prime in 600 s, against 4 to 9 s for the whole certificate with alpha/beta slicing | found: the control certified with r = 5 and a4 = 3 |
| generic hyperplanes within the alphas only, and within the betas only (`commcert_ab.py`) | small, and the equations keep their (alpha-polynomial) × (beta-linear) form | found: the control certified with r = 1 (two seeds); the relaxation hit g27_Fp2b_m3_s2 certified with r = 4 (two seeds; 318- and 438-bit entries; points involve i, as ex2-modp found) |

The third way is the one to use on our systems. A generic affine subspace of alpha-space meets every piece of the alphas' solution set of the right dimension, however special its position. In small-integer draws the subspace can hit a degenerate spot. `special_slice.py` (coefficients in −9..9) refused seeds 1 and 2 and certified seeds 3 and 4. I checked seed 1 by hand: its line meets a4 = 3 exactly where a1 + a2 = 0. `commcert_ab.py` draws from −30..30. A degenerate draw is refused, never wrongly certified.

**Calibration** [computed; `out_commcert.txt`, `out_special_graphA_h5.txt`, `remote-out/`]:
- *Certified exactly when solvable:*
  - Karatsuba rank 3: r = 4 under coordinate slices, 3 primes; r = 26 under generic slices, 15 primes, 207-bit entries. The eigenvectors of the 4 × 4 matrices give 4 numerical solutions with residual ≤ 3.4·10^−12. The 26 were not checked numerically; the exact check does not need that.
  - W-state rank 3: r = 2, 4 primes; r = 48 under generic slices, 54 primes, 796-bit entries.
  - nag-rerun's special-position control: r = 1 or 5.
  - g27_Fp2b_m3_s2: r = 4.
  - 41 random systems:
    - 20 with a planted point (r = 1);
    - 20 square systems, with r between 1 and 9 and never above the Bézout number (4, 8 or 9); their solutions are generically irrational;
    - 1 "impossible" draw with a genuine rational zero.
- *Refused when impossible:* W-state and Karatsuba at rank 2 (refuted modulo the first prime, PC degree 5); 19 random impossible draws; and a decoy.
- *The decoy* is {x − 1, x^2 − 1 − p0}, which has the clean solution x = 1 modulo p0 = 2^31 − 1 and no complex solution.
  - With several primes it is refused, because the second prime refutes it.
  - With p0 alone, the reconstructed matrix [1] fails the exact check (−p0 ≠ 0).
  - This is ex2-modp's single-prime trap, caught by the exact step.

**Limits.**
- Heights grow with r: the W-state's 48 × 48 matrices under generic slicing needed 54 primes (796-bit entries), and the exact check of 48 × 48 matrices took about 20 minutes in plain Python fractions.
- The search is only as strong as the Gröbner or flat step [computed, `remote-out/cs/`]:
  - certificates came in seconds for 8 to 15 unknowns (counted after slicing);
  - at 25 to 59 unknowns the modular Gröbner basis did not finish one prime within 300 to 600 s, for each of:
    - g27_Fany4_s12 (25 unknowns);
    - Toom-3, the rank-5 unit for multiplying two quadratics (40 unknowns after gauge and 5 slices, both with generic slices and with the evaluation points fixed);
    - Strassen's rank-7 2 × 2 matrix multiplication (59 unknowns after 11 slices, with generic slices and with 11 coordinates u_q1, v_q1 fixed).
- So this does not move the size frontier of section 2.2. It turns every "yes" that the search can reach into a short exact proof, whatever the field of the unit's numbers. The frontier is set by the search, at a few dozen unknowns for systems of this kind.

## 3. What decides the degree: approximate solutions at infinity

**Curve lemma** [proved here]. Suppose the f_i have no common zero, and points x(t) with ‖x(t)‖ → ∞ satisfy max_i |f_i(x(t))| ≤ C‖x(t)‖^(−L). Let I be the equations that are not exactly zero along x(t). Then every NS refutation 1 = Σ g_i f_i has deg g_i ≥ L for some i ∈ I, so its degree is at least L + min{deg f_i : i ∈ I}.

*Proof.* Evaluate the identity on x(t). Only the equations in I contribute, so 1 ≤ Σ_{i∈I} |g_i(x(t))|·|f_i(x(t))| ≤ C′‖x(t)‖^(max_{i∈I} deg g_i − L). If every deg g_i with i ∈ I were below L, the right side would tend to 0. ∎

In unit language: an approximate (border) unit whose coefficients grow like ε^(−b) while its residual falls like ε^a has L = a/b. It forces every NS refutation of the exact shape above degree a/b + min deg f. Below that degree, pseudo-solutions exist without any unit, by section 2.1.

**The other direction** [read: first-hand, as reported by Kollár; constants from memory].
- Brownawell's proof combines a lower bound max_i |f_i(x)| ≥ c‖x‖^(−λ) at infinity (a Łojasiewicz inequality) with Skoda's division theorem to bound the refutation degree.
- Kollár's Proposition 1.10 bounds that exponent by his degree bound N.
- His Corollary 1.9 locates the worst case: with k = n, refutations below degree N fail exactly when the closures of the hypersurfaces meet in one point, with one-dimensional tangent space, outside the hyperplane at infinity.

**Inside our own class of systems** [computed]. Consider the Masser-Philippon-Brownawell chain:

x_1·y_1, x_{k−1} − x_k·y_k, y_k − x_k, 1 − x_{n−1}·x_n.

The two groups are β = {x_1, …, x_{n−1}} and α = {y_1, …, y_{n−1}, x_n}. Every equation is linear in β with coefficients affine in α, so the system reads A(α)β = e(α). That is w3-asym's form, except that the right side may depend on α.
- Along x_n = t, x_{n−1} = 1/t, x_{n−2} = 1/t^2, … only x_1·y_1 fails, by t^(−2^(n−1)). So the lemma gives NS degree ≥ 2^(n−1) + 2 = 4, 6, 10 for n = 2, 3, 4.
- The exact NS degrees are 4, 7, 12 [computed; exact over Q for n = 2, 3].
- But the PC degrees are 3, 3, 4, 4, 5 for n = 2, …, 6 (`out_pc_chains.txt`).
- Chains whose contradiction sits at finite distance have linear NS degree n + 1 (`out_chains.txt`).

So:
- For the identity format no degree bound below exponential holds on the class "A(α)β = e(α)" (P3). Approximate solutions at infinity are exactly what drives the identity degree up.
- That obstruction does not carry over to Polynomial Calculus, because a PC derivation multiplies by variables step by step and does not keep the decay. A cheap existence test should therefore use the PC closure, as in section 2.3.
- Whether some other structure at infinity bounds PC degree from below is open here [guess: it does for some systems, as known PC lower bounds for combinatorial principles suggest, but I have no example in our class].

**Calibration of the lemma** [computed].

| system | lemma's lower bound | NS degree | PC degree |
|---|---|---|---|
| W-state, rank 2 (border rank 2; the approximation ((e_1 + εe_2)^⊗3 − e_1^⊗3)/ε gives L = 3) | 6 | 9 | 5 |
| Karatsuba, rank 2 (no approximate solutions: the flattening has rank 3) | nothing | 9 (equal to the flattening identity's degree) | 5 |

On our own systems [computed, `out_census_nag.txt`] I took the 74 impossible a2-nag systems (Singular basis {1}; 9 to 29 unknowns). Their NS degrees are 2 (1 system), 3 (8), 4 (39) and 5 (26), and their PC degrees are the same except in one case (5 → 4).
- Three of them have approximate solutions at infinity: a2-nag's Levenberg-Marquardt floor falls to 7.6·10^−11, 2.2·10^−7 and 5.2·10^−8 at |α| ≈ 10^3 to 10^5. All three sit at the top of the range, NS degree 5.
- But 23 of the other 71, whose residual floors stay at 0.1 to 0.6, also reach 5.

So approximate solutions push the identity degree up, as the lemma says, but they are not the only cause. The lemma is a lower bound, not a predictor.

**For the impossibility side** (ml-proofcx's question of which degree a refutation needs).
- Every approximate unit of order L gives the exact shape an NS degree lower bound of L + min deg f, and NS pseudo-units one degree below that, for free [proved].
- Conversely, an NS refutation of degree D rules out approximate solutions of order L > D − min deg f. For example, all 58 saved w3-asym skeletons (NS degree ≤ 3, min deg f = 1) have no approximate solutions whose residual falls faster than (size)^(−2) [derived here from computed].
- None of this bounds PC or Gröbner degree.

### 3.1 Existence by counting: an existence proof made of an impossibility proof (second session)

Section 2.5 says every existence proof can be turned into a finite model, but a proof need not exhibit one. The classical way to prove that solutions exist without finding any is to count them.

**The counting route** [textbook; argued for our setting].
- Take n equations in n unknowns, of degrees d_1, …, d_n.
- Suppose their top-degree parts have no common zero other than 0, that is, no solution "at infinity".
- Then Bézout's theorem says the system has exactly d_1⋯d_n solutions, counted with multiplicity. In particular it has at least one.
- The hypothesis is an impossibility statement: a homogeneous system with no nonzero solution. It has a Nullstellensatz-type certificate (every x_j^N is a combination of the top-degree parts) of degree at most Σ(d_i − 1) + 1 [memory: Macaulay's bound]. That is linear in n, not exponential.
- So this is existence proved by an impossibility proof, in the most literal sense, and with a short certificate. It is also non-constructive: it says nothing about where the solutions are or what field their numbers live in.

**Why it does not apply to our systems** [computed, `infinity_check.py`, `remote-out/`].
1. *Our systems have surplus equations.* Unit systems have more equations than unknowns. Counting gives solutions of a square subsystem, which must then satisfy the rest. That needs an extra identity putting the surplus equations in the subsystem's ideal (compare ex2-modp's Proposition 2).
2. *The calibration systems that are square are full of solutions at infinity.* The coordinate-sliced Karatsuba rank-3 system and the W-state slice happen to be square, so the route can be tested on them directly:

| system | Bézout count | finite solutions | solutions at infinity |
|---|---|---|---|
| Karatsuba rank 3, gauge + 3 slices (12 × 12) | 1,728 | 4 | a family of (projective) dimension 6 |
| W-state rank 3, gauge + 4 slices (8 × 8) | 144 | 2 | a family of dimension 3 |

   The solutions at infinity are directions along which the cubic parts cancel [argued]. One example is a term whose w-part is zero while its u- and v-parts are arbitrary. Border decompositions, the W-state's kind of approximation (section 3), run off to infinity along such directions.
3. *So counting proves existence only in a compactification.* There the answer includes border units and degenerate limits. Turning that into a real unit needs a refutation of every boundary piece. That fails here, because the boundary pieces exist.

This matches MAP's theta-shadow picture from the existence side. Counting arguments see the closure of the set of units, approximate units included, and only exact algebra separates the two. It also matches the curve lemma (section 3): the same solutions at infinity that force long refutations of impossible shapes block counting proofs for solvable ones.

**Why the boundary is never empty for decompositions** [argued]. Any system that writes a target as a sum of at least two terms has solutions at infinity. Two equal and opposite terms (u_1 = u_2, v_1 = v_2, w_1 = −w_2) can grow without bound while their sum stays zero, so their top-degree parts cancel. The same holds for our shapes whenever two moves can cancel each other.

**What would revive it** [guess]. It would take a compactification adapted to units: one that blows up or excludes such cancelling pairs, keeps the alpha and beta groups apart, and tracks the frame structure. In it the boundary would have to be provably empty (refutable at low degree) and the count nonzero. I know of no such compactification for the bilinear-algorithm setting.

## 4. Calibration (both directions)

| case | expected | result | file |
|---|---|---|---|
| plain 2 × 4 control (w3-asym sing_ctl/plain) | unit | certified: NS D = 2 (flat s = 1, rank 1); PC D = 2 | out_census_sing.txt, out_census_pc.txt |
| plain + free extra summand (plain_extra) | unit | certified: NS D = 5 (at D = 4 the top monomial a0^4 is truncation noise); PC D = 2 | same |
| plain minus a move (plain_minus) | impossible | refuted at degree 0 (it contains the equation −1 = 0) | out_census_Q.txt |
| 58 Singular-inconsistent skeletons (sing/ 25, sing_sig2/ 7, sing_l5/ 26; part of w3-asym's 755) | impossible | all refuted. NS exactly over Q: 15 at degree 2, 43 at degree 3 (the same with two primes). PC: all 58 at degree 2 | out_census_Q.txt, out_census_sing.txt, out_census_pc.txt |
| Karatsuba rank 3, gauge + 3 specialised | unit | certified mod p, rank 4: NS D = 7, PC D = 4 | out_tensors_K3.txt, out_pc.txt |
| W-state rank 3, gauge + 4 specialised | unit | certified mod p, rank 2: NS D = 7, PC D = 4 | out_tensors_W3_k4.txt, out_pc.txt |
| independent checks of those two | | Singular: Gröbner bases with 4 and 2 solutions; Newton: exactly 4 and 2 distinct nonsingular solutions | gbcmp/out_singular_gb.txt, out_crosscheck.txt |
| W-state rank 2 (border rank 2) | impossible | never flat; NS refutation at exactly 9 (none at 8), PC at 5 | out_tensors_W2.txt, remote-out, out_pc.txt |
| Karatsuba rank 2 | impossible | never flat (at s ≥ 2); NS refutation at exactly 9 (none at 8), PC at 5 | out_tensors_K2.txt, remote-out, out_pc.txt |
| 74 impossible a2-nag systems (meta2/a2-nag sys1, sys2: Singular basis {1}) | impossible | all refuted mod p: NS degree 2 to 5, PC degree 2 to 5; never certified | out_census_nag.txt |
| 200 random small systems (n = 2, 3; d = 2, 3; seed 11; a first run of 60, seed 7, agreed) | mixed | all ended certified or refuted by D ≤ 14 (NS version); see the breakdown below | out_soundness.txt |
| block-2 VV stand-in 2 in 3 (known solvable) | unit | **blind (too big)**: 1,672 unknowns, degree 4; the smallest admissible flat test (s = 3, D = 6) needs about 3·10^16 monomials | — |
| engine 2 in 5 at block 3 (proved impossible) | impossible | whole shape **blind (too big)**: exact equations with about 34·|B| = 9.4·10^7 unknowns. Skeleton level **refused**: a2-nag's six calibration skeletons cal_b3_2in5_s0, s2, s4 to s7 (9 to 10 unknowns) are refuted at NS degree 2 to 4 (PC the same) | out_census_nag.txt |
| 8-move control word (moonshot/o3/o3-build) | unit | not a nonlinear test: once its letters are fixed its free numbers enter linearly, and the linear solve is a degree-1 certificate | — |
| *second session, exact commuting-matrix certificate (section 2.5)* | | | |
| Karatsuba rank 3 (coordinate slices / generic slices) | unit | **certified exactly over Q**: r = 4, 3 primes, 20-bit entries / r = 26, 15 primes, 207-bit entries | out_commcert.txt, remote-out/ |
| W-state rank 3 (coordinate slices / generic slices) | unit | certified exactly: r = 2, 4 primes, 37 bits / r = 48, 54 primes, 796 bits (1,257 s, mostly the exact check; certificate not saved: the run's certs/ folder was missing) | out_commcert.txt, remote-out/cs/ (one Colab machine's folder) |
| nag-rerun's special-position control spu_shift.u (solvable only on a4 = 3) | solvable, special | coordinate slices: refuted (missed, as expected); generic slices: **certified exactly**, r = 5; alpha/beta slices: certified, r = 1 (2 of 2 seeds; 2 of 4 seeds with the narrower coefficient range); a4 = 3 found with no hints | out_special_graphA_h5.txt, remote-out/ |
| relaxation hit g27_Fp2b_m3_s2.u (solutions need i) | solvable | alpha/beta slices: **certified exactly**, r = 4 (2 of 2 seeds; 318 and 438 bits); alpha-only slices: r = 4 (203 bits); field not supplied. Generic slices in all unknowns: Singular did not finish one prime in 600 s | certs/Fp2b_ab_seed2.json, remote-out/ |
| relaxation hit g27_Fany4_s12.u (solutions need a cube root of unity; 90 equations, 25 unknowns) | solvable | **not reached** (search too slow): alpha-space slices with h_a = 16 down to 9 are all empty modulo p (2 seeds; h_a = 9 took 65 to 73 s), and the Gröbner basis at h_a = 8 did not finish in 300 s. Singular's dimension of the unsliced system did not finish in 900 s (twice), so generic slicing could not even be set up. ex2-modp also found this one hard (Singular over Q > 300 s, msolve 115 s) | remote-out/cs/ |
| 60 random systems (seed 3) | mixed | 41 certified exactly (20 planted, 20 square, 1 impossible-looking draw with a genuine zero); 19 impossible refused | out_commcert.txt |
| W-state, Karatsuba at rank 2; decoy {x − 1, x^2 − 1 − p0} | impossible | refused; the decoy's single-prime reconstruction fails the exact check | out_commcert.txt |

The 200 random systems broke down as follows:
- 57 planted systems were certified (56 with rank 1), and the point was recovered and checked exactly over Q.
- 69 square systems were certified with rank at most the Bézout number.
- 73 impossible draws were refuted at degree 3 to 7.
- 2 "impossible" draws were certified. Each has a genuine rational zero, recovered with residual exactly 0.
- No false certificate.
- The PC version (`soundness_pc.py`, `out_soundness_pc.txt`) gives the same verdicts, ranks and points, with refutations at degree 2 to 4 instead of 3 to 7.

So the method passes calibration wherever it can run. On the stand-ins that carry the real structure it is blind, for size reasons.

## 5. A use on the search side: the skeletons w3-asym could not decide

w3-asym decided 697 sweep skeletons with Singular and left 133 undecided ("closures over 40 columns": 62 in the (2,2) sweep, 71 in the sigma sweep). `sweep/sweep_ns.py` re-runs its exact sweep with the same seeds (801, 802), the same pools and the same subsets. Two things change:
- The closure cap is raised to 150.
- Singular is replaced by an NS refutation test at degree ≤ 3. Where the full degree-3 space is too big, a partial one is tried, with multipliers in the program coefficients only. A hit there is still an exact identity.

Results [computed, mod p]:
- **Core 0 at cap 40** reproduces w3-asym's tally exactly: 45 refuted and 15 capped, as with Singular.
- **(2,2) sweep at cap 150.** 45 of the 62 previously undecided skeletons are now refuted: 41 at degree 2 or 3, and 4 through the partial space at degree 3, including one with 244 unknowns. The other 17 still have closures over 150 columns. None survived degree 3.

- **Sigma sweep at cap 150.** 28 of the 71 previously undecided skeletons are now refuted: 27 at degree 2 or 3 in the full space (up to 244 unknowns), and 1 through the partial space at degree 3 (229 unknowns). The other 43 still have closures over 150 columns, 26 of them from one core. None survived degree 3. (`remote-out/`, the sweep folder of one Colab machine)
- **Both sweeps at cap 600** (second session, `sweep/job/sweep/sweep_ns_big.py`, only the 60 skeletons still over 150). 3 more are refuted, all through the partial space:
  - one at degree 2 with 1,156 unknowns (in 5.8 s);
  - two at degree 3 with 564 and 484 unknowns.

  The other 57 have closures over 600 columns. (`remote-out/ssh-cap600/`)

Totals: of the 133 skeletons w3-asym left undecided, 76 are refuted (45 + 28 + 3) and 57 remain undecided because of closure size, not degree. I had predicted that cap 600 would decide at least half of the 60. That was wrong: the closure, not the refutation degree, is what limits this family.

So every skeleton of this family that has been decided is impossible, and impossible for a short algebraic reason (degree 2 or 3). That fits w3-asym's Theorem R picture: the junk terms cannot cancel, and that is a low-degree fact.

## 6. Pseudo-calibration and other abstract constructions [argued]

- **Pseudo-calibration** (Barak, Hopkins, Kelner, Kothari, Moitra and Potechin 2016) [memory] builds low-degree pseudo-expectations for random instances from a planted distribution.
  - It proves that refutations need high degree, which is the impossibility side's barrier question (ml-proofcx).
  - It needs a random instance and a planted family. We have one fixed system per shape.
  - By section 2 its output proves existence only at D_eff, or if it happens to be flat.
  - The nearest thing to a planted family is an approximate unit. By section 3 those make misleading NS pseudo-solutions, not certificates. Our numerical near misses do not even approach residual 0 (w3-asym seed 301 stops at 0.2333, mw2-absorb4 at 2/3), so they plant nothing.
- **Symmetry-reduced linear algebra** is exact: Reynolds averaging keeps L(1) = 1, the vanishing, and finite rank. It saves a factor |G|. That can matter for flat tests on mid-size shapes, but not against D_eff.
- **Any abstract construction** must output a functional that is flat or of degree ≥ D_eff. A flat functional has finite rank, so it is a finite combination of genuine points (Laurent-Mourrain). So an abstract flat pseudo-solution is an abstract description of actual units. This route cannot prove existence without, in effect, locating units.

## 7. What would make it not hopeless

1. **Finitely many solutions, none at infinity.**
   - Gauge-fix and slice the shape until its solution set is finite, and check that the leading forms have no common zero at infinity (itself a Nullstellensatz question of degree about n(d−1)+1).
   - Then a refutation, if one exists, has degree at most Σ(d_i − 1) + 1 [memory: Macaulay's bound, Lazard 1983], and the flat test is complete at that degree. That replaces D_eff ≈ d^n by a degree linear in n.
   - The size is still about C(n + nd, n), so this is feasible up to roughly 15 to 20 unknowns [derived here]. For the 12-unknown Karatsuba system the bound is 25; the actual degree was 4 (PC).
   - Section 3.1 tested the hypothesis: the sliced Karatsuba and W-state systems both have large solution sets at infinity, so this item does not apply to them as written. It would need a compactification in which the boundary is empty.
2. **A structural theorem for skeleton systems.**
   - Example: "every impossible skeleton of a given family has a PC refutation of degree ≤ c". Then a degree-c PC pseudo-solution of a skeleton (plus slicing and flatness) would prove a unit exists.
   - Our data fit c = 3 for w3-asym's family (all 58 saved and all 76 newly decided skeletons), and c = 5 for a2-nag's (74 of 74). A single constant across families would have to be at least 5.
   - Section 3 shows that the identity-format version cannot hold for general systems of the form A(α)β = e(α). So a proof would have to use the closure and frame structure (w3-asym's Theorem R).
3. **Fewer unknowns.** For our skeletons the degree is not the bottleneck (2 to 3); the number of unknowns is. Every exact reduction that brings a shape down to tens or a few hundred unknowns (closures, forced readouts as in Theorem R, symmetric ansätze) makes both the refuter and the flat test practical. That is the regime where Gröbner bases already work.

## 8. Predictions (STATUS.txt, written before computing) and outcomes

| | prediction | outcome |
|---|---|---|
| P1 | D_eff astronomically large; log10 ≥ 100 for open shapes | right: log10 D_eff from 17.7 (37-unknown skeleton) to 7.8·10^11 (block-4 2 in 7); sizes 10^610 and up |
| P2 | sparse or multihomogeneous bounds stay hopeless | right: 10^610 → 10^357 for one skeleton |
| P3 | the worst case is tight; a class-wide bound must use the specific system | right for NS, now inside our class (chain degrees 4, 7, 12); **not** for PC (3, 3, 4, 4, 5) |
| P4 | flat pseudo-solutions are the useful object, of Gröbner strength | right; computed on calibration systems (PC-flat at D = 4 where Singular's basis has degree 2, with the same solution counts) |
| P5 | calibration both ways is feasible (0.6) | right for plain, skeletons, tensors and random systems; blind for the VV stand-in and block-3 2 in 5 |
| P6 | refutation degrees of the skeletons are low (2 to 6) | right: w3-asym NS 2 or 3 and PC 2 (all 58); a2-nag 2 to 5 (all 74) |
| P7 | positive-dimensional solution sets need slicing | right |
| P8 | pseudo-calibration does not transfer; symmetry gains only a factor |G| | argued, as predicted |
| P9 | hopeless for open shapes; no new unit (probability < 0.02) | right: no new unit |
| P10 | practical only for zero-dimensional, low-regularity systems with ≤ 40 unknowns | right in spirit; the sharper form is section 7 item 1. The refuter alone reaches 244 unknowns |

Not predicted:
- The curve lemma: approximate units force long identity refutations.
- The size of the NS versus PC gap. On the chains it is exponential versus linear.

Second session (predictions Q1 to Q5 in STATUS.txt, written at 23:06 before any new computing or reading):

| | prediction | outcome |
|---|---|---|
| Q1 | exact commuting-matrix certificates work on Karatsuba (r = 4), the W-state (r = 2) and planted systems with ≤ 12 primes (0.6) | right: 3, 4 and 2 primes. Generic slicing needs more: up to 30 for g27_Fp2b (438-bit entries), 54 for the W-state at r = 48 (796 bits) |
| Q2 | the single-prime decoy is refused (0.95) | right, at the second prime; with one prime alone the exact check rejects it |
| Q3 | a "many primes" count from the arithmetic Nullstellensatz is practical only for n ≤ 4 (0.7) | not computed here: ex2-modp's Proposition 5 already did it (2,636 primes for 3 unknowns, about 500,000 for 5); consistent with Q3 |
| Q4 | cap 600 decides at least half of the 60 remaining skeletons (0.5) | **wrong**: 3 of 60. The closures exceed 600 columns |
| Q5 | no new unit (0.98) | right |

Not predicted in the second session:
- Structured slicing (alphas and betas separately) finds special position while keeping r small.
- The counting route's obstruction, solutions at infinity, is the same object as the curve lemma's.

## 9. Next steps, ranked

1. **Use `nsfool.py` (PC closure plus partial multipliers) as a Singular-free exact refuter in skeleton searches.** It reproduced Singular everywhere it was tested, takes milliseconds at w3-asym's sizes, and reaches 244 unknowns. A skeleton that survives PC degree 3 or 4 deserves a full look. It is either a unit or an impossible system with long refutations; impossible systems of that kind exist, for example Karatsuba at rank 2.
2. **Make `commcert_ab.py` the existence end of every skeleton pipeline.** Any skeleton that survives the refuter goes to structured slicing plus the commuting-matrix certificate. A "yes" then arrives as an exact proof, whatever field the unit's numbers need. (Done in the second session, in place of the planned "multimodular flat certificate".)
3. **State the certificate in Lean** (on a Colab): "commuting matrices M_i with f_j(M) = 0 imply a common zero". It needs one lemma, the common eigenvector of a commuting family over C. Every future existence claim can then be checked as rational matrix arithmetic [guess: a short file, given Mathlib's eigenvector API].
4. **Skeptic.** Check the corollary in 2.3 (the relaxed degree condition, and its use with PC spaces), the curve lemma, the sweep reproduction, the theorem in 2.5 (including the r ≤ d^n remark, from memory), and the infinity computation in 3.1.
5. **Do not spend machines** on worst-case-degree pseudo-solutions, on NS-degree (rather than PC) tests, or on pseudo-calibration.

## Sources

Literature:
- J. Kollár, "Sharp effective Nullstellensatz", J. AMS 1 (1988) 963–975: Theorems 1.1 and 1.2 (Brownawell's, as quoted), Theorem 1.5, Corollary 1.9, Proposition 1.10 and Example 2.3 [read: first-hand]
- M. Sombra, alg-geom/9610006, introduction and main theorem [read: first-hand]
- C. D'Andrea, T. Krick and M. Sombra, arXiv:1103.4561, Theorem 0.1 and the sharpness example [read: first-hand]
- M. Laurent and B. Mourrain, arXiv:0812.2563 (Arch. Math. 93, 2009), Theorems 1.4 to 1.6 and 3.2 [read: first-hand]
- W. D. Brownawell, Ann. of Math. 126 (1987) [read: as quoted by Kollár and in search snippets]
- Barak, Hopkins, Kelner, Kothari, Moitra and Potechin (2016); Macaulay's bound; Lazard (1983); Polynomial Calculus (Clegg, Edmonds and Impagliazzo 1996) [memory]
- Second session, all [memory]: Zariski's lemma; common eigenvectors of commuting matrices; J. Heintz, "Definability and fast quantifier elimination in algebraically closed fields", Theoret. Comput. Sci. 24 (1983), for the Bézout inequality; Macaulay's bound for homogeneous systems with only the trivial zero. No new literature was read first-hand in the second session.

Downloads were read with pdftotext only. They are in a temporary folder.

Related project work:
- `../../exist2/ex2-modp/NOTE.md`: exact refutation identities for nag-rerun's skeletons; Proposition 1 on single primes; Proposition 5, the "many primes" count (the reason this track does not recompute it); the special-position control and the two irrational relaxation hits used in section 2.5.
- `../../meta/ex-note/NOTE.md`: the earlier verdict that existence means exact algebra on one fixed shape.

Files in this folder:
- STATUS.txt (predictions) and RUNLOG.txt.
- `nsfool.py`: NS echelon mod p and over Q, the PC closure, pseudo-solutions, the flat scan, `certify`, `pc_certify` and `ns_partial`.
- `census.py`, `census_Q.py`, `census_pc.py` and `census_nag.py`, with `out_census_sing.txt`, `out_census_Q.txt`, `out_census_pc.txt` and `out_census_nag.txt`.
- `tensors.py` with `out_tensors_*.txt`, and `pc_run.py` with `out_pc.txt`.
- `chains.py` with `out_chains.txt`, and `pc_chains.py` with `out_pc_chains.txt`.
- `crosscheck_K3.py` with `out_crosscheck.txt`; `export_sing.py` with `gbcmp/`.
- `deff.py` with `out_deff.txt`; `soundness.py` with `out_soundness.txt`.
- `gen_ctl.py` and `sys/` (trivial controls).
- `sweep/sweep_ns.py`, `sweep/job/` (the Colab package) and `remote-out/` (Colab outputs).
- Second session:
  - `commcert.py` (the exact certificate: flat search mod p, multiplication matrices, reconstruction, `exact_check`), with `commcert_run.py` → `out_commcert.txt`;
  - `special_slice.py` → `out_special_graphA_h5.txt`;
  - `commcert_sing.py` (Singular search, generic slicing) and `commcert_ab.py` (alpha/beta slicing);
  - `infinity_check.py` (section 3.1);
  - `check_cert.py` (standalone certificate checker) and `certs/` (four exact certificates);
  - `cs_job/` (the Colab package);
  - `sys/spu/` and `sys/relax/` (copies of nag-rerun's msolve exports, read with exponents accumulated);
  - `sweep/job/sweep/sweep_ns_big.py` (cap 600);
  - `remote-out/` (the Colab outputs of the sweeps at caps 150 and 600, and of `cs/`).

Inputs (read only):
- `../../moonshot/w3-asym/` (sing/, sing_l5/, sing_sig2/, sing_ctl/, exact_sweep.py and its library, NOTE.md);
- `../../meta/ex-dimension/out_exact_counts.txt`;
- `../../meta/ex-v-ex-dimension-logic/out_indep.txt`;
- `../../meta2/a2-nag/` (sys1/, sys2/, agg_sys1.txt, agg_sys2.txt).
- `../../meta2/nag-rerun/ms/` (spu_shift.u, g27_Fp2b_m3_s2.u and g27_Fany4_s12.u exports, via ex2-modp's data/msctrl/ for the first).
