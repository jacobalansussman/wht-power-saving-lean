> **Working note, kept as written.** An AI agent wrote this during the logic round of 14-15 October 2026 (the ml-proofcx track). Private paths and machine names
> were removed and nothing else was changed. No person has checked it. Where it disagrees with this section's
> [README](../README.md), the README wins, because it follows the later skeptic verdicts. Its skeptic upheld Proposition D and the three exact border curves, and refuted the step to "no Theorem-R
> certificate can ever exclude 2 in 7 at block 4" as incomplete. It also corrected "rank-type" to "closure-type".

# ml-proofcx: proofs (14 October 2026)

Nobody else has read this. Nothing posted. Labels: [proved] (proof here), [argued], [read: ...], [memory], [guess].
The plain-English account is NOTE.md; this file holds the statements and their short proofs.

## 0. Setting

Fix a game, a block size m, a shape (W arrays, M paid moves) and a **skeleton** S: the letters and the wiring are
fixed, the complex coefficients are unknown. Put all unknowns (content coefficients and readout coefficients) in one
vector c in C^n. Fix where correctness is tested: one representation (f digit columns) or the group ring (all f at
once; ex-note, "a formal identity in a group ring"). The program's outputs are then a polynomial map
phi_S : C^n -> C^N, and the target (the W transformed arrays) is a point T in C^N.

- exact unit with skeleton S  <=>  T in phi_S(C^n)  <=>  the system F_S(c) := phi_S(c) - T = 0 has a complex root;
- border (approximate) unit  <=>  T lies in the Euclidean closure of phi_S(C^n). The image is constructible
  (Chevalley), and for a constructible set the Euclidean and Zariski closures agree [memory: Mumford, Red Book I.10].
  Call this closure Z_S.

A shape has finitely many skeletons in every game with a finite letter group, so "the shape has an exact (border)
unit" is a finite union of these conditions.

**Refutation systems** (standard): a Nullstellensatz (NS) refutation of F = (f_1..f_k) is 1 = sum g_i f_i; its degree
is max deg(g_i f_i). A Polynomial Calculus (PC) refutation of degree d derives 1 using only polynomials of degree <= d,
by linear combinations and multiplication by a variable. An Ideal Proof System (IPS) refutation is a polynomial
C(c, y) with C(c, 0) = 0 and C(c, F(c)) = 1, measured by algebraic circuit size (Grochow-Pitassi). A degree-d NS
pseudo-solution is a linear functional L on polynomials of degree <= d with L(1) = 1 and L(m f_i) = 0 whenever
deg(m f_i) <= d; by linear duality it exists iff no NS refutation of degree <= d exists.

## 1. Lemma A (closure refutations, and exactly when they exist) [proved]

(i) If Q in C[y_1..y_N] satisfies Q(phi_S(c)) = 0 identically in c and Q(T) = 1, then F_S has an NS refutation of
degree <= deg Q * deg phi_S.

(ii) Such a Q exists iff T is not in Z_S, i.e. iff skeleton S has no border unit.

(iii) Every certificate whose bound also holds for border units (pf-ladder B1: all Theorem-R certificates, hence
hooks, Venn sets, coordinate sets, lifts) refutes only shapes with no border unit, and then (i)-(ii) apply to each of
the shape's skeletons.

*Proof.* (i) Write Q(y) - Q(T) = sum_i (y_i - T_i) H_i(y) with deg H_i <= deg Q - 1 (expand Q around T). Put
y = phi_S(c): -1 = sum_i (phi_i(c) - T_i) H_i(phi_S(c)), each term of degree <= deg phi_S * deg Q. (ii) The vanishing
ideal of Z_S is ker(phi_S^*) = {Q : Q o phi_S = 0}, and T is outside the zero set Z_S of that ideal iff some member is
non-zero at T; scale it. (iii) B1 (upheld) says the bound holds for border units, so a shape it excludes has T outside
every Z_S. QED.

So the whole class of "closure" refutations (those that factor through a polynomial on output space) is complete
for "no border unit" and blind exactly when a border unit exists. Rank methods are a subclass of it. The converse (ii)
is standard (session-2 correction, section 7: Grochow arXiv:1304.6333 section 3.2); only its use for our skeletons is ours.

**Theorem R as an IPS / NS proof** [argued]. In module form Theorem R bounds the rank of a matrix A(c) built from the
program by T f(D) for every c (so every (T f(D)+1)-minor of A is the zero polynomial), while the unit equations imply,
by a linear derivation, that a constant matrix B of rank W dem > T f(D) equals Lambda(c) A(c). Taking a non-singular
(T f(D)+1)-submatrix of B and Cauchy-Binet gives 1 = det(B_IJ)^-1 (det B_IJ - det(Lambda A)_IJ), a member of the
ideal. This is an IPS refutation whose circuit is a determinant (polynomial size in the certificate's dimension) and an
NS refutation of degree about (T f(D)+1)(word depth + 1). For the block-4 hook at f = 9 the ranks are near 7 x 10^4
(meta/MAP.md W2.5), so the direct NS degree is of order 10^5: Theorem R is a short proof in IPS and a huge-degree proof
in NS. Whether a much lower NS degree exists for the whole shape is open.

## 2. Lemma B (border decay forces Nullstellensatz degree; a Lojasiewicz bound) [proved]

If 1 = sum_i g_i f_i with deg g_i <= D, then for every x in C^n

    |F(x)| >= (sum_i ||g_i||_1^2)^(-1/2) * (1 + |x|_inf)^(-D),

where ||g||_1 is the sum of the absolute values of g's coefficients. Hence, if |F(x_k)| <= C |x_k|^(-a) along some
x_k -> infinity, every NS refutation has max deg g_i >= a, and NS degree >= a + min_i deg f_i.

*Proof.* |g_i(x)| <= ||g_i||_1 (1 + |x|_inf)^(deg g_i), and 1 = |sum g_i(x) f_i(x)| <= (sum |g_i(x)|^2)^(1/2) |F(x)|
(Cauchy-Schwarz). QED.

This is a2-nag's Lemma 2(b) for general certificates; the converse direction (bounds on NS degrees from Lojasiewicz
exponents) is the classical one [read: abstract, Kollar arXiv:math/9805091 lists "Lojasiewicz-type inequalities" among
the applications of the effective Nullstellensatz].

## 3. Lemma C (pseudo-units from border curves) [proved]

Let x(eps) in C((eps))^n have pole order <= b (eps^b x(eps) has no negative powers) and F(x(eps)) = O(eps^a), a >= 1.
Then L(p) := (coefficient of eps^0 in p(x(eps))) is a degree-d NS pseudo-solution for every
d <= min_i deg f_i + ceil(a/b) - 1. So NS degree >= min_i deg f_i + ceil(a/b).

*Proof.* L is linear and L(1) = 1. For a monomial m, m(x(eps)) has pole order <= b deg m, so m(x(eps)) f_i(x(eps)) =
O(eps^(a - b deg m)), whose eps^0 coefficient is 0 when b deg m < a. If deg(m f_i) <= d then
deg m <= d - min_j deg f_j <= ceil(a/b) - 1 < a/b. QED.

An exact solution x0 gives the evaluation functional p -> p(x0), a pseudo-solution of every degree.

**Polynomial Calculus escapes Lemmas B and C** [computed, nsdeg2.py, mod p and checked by a second implementation at
k = 4]. The chain sq_k: x1 y = 1, x_(j+1) = x_j^2 (j < k), x_k = 0 has the border curve x1 = eps, y = 1/eps,
x_j = eps^(2^(j-1)) with a = 2^(k-1), b = 1. NS degrees 4, 7, 12 for k = 2, 3, 4 (Lemma C: >= 3, 5, 9) and > 9 for
k = 5, 6 (Lemma C: >= 17, 33); PC degrees 3, 3, 4, 4, 5. So border decay bounds NS degree but not PC degree.

## 4. Lemma H (halving is a Nullstellensatz statement over F_2) [proved]

Suppose F has coefficients in Z[i] (true for engine and widened letters with the unnormalised Hadamard). If the
reduction of F modulo the prime (1+i) (residue field F_2) has no root over the algebraic closure of F_2 -- equivalently
1 lies in the ideal of F over F_2, which NS or PC over F_2 can certify -- then every root of F over the algebraic
numbers has a coordinate of negative 2-adic valuation, for every extension of the 2-adic valuation: the unit "divides
by 2" (possibly by sqrt 2 or 1 + i, which have valuations 1/2).

*Proof.* If all coordinates have valuation >= 0 they lie in the valuation ring; reducing modulo its maximal ideal gives
a root of F mod (1+i) over an algebraic extension of F_2. QED.

This is the frame of Prop G and Theorem F (MAP 1.1, 1.3) in proof-complexity form: those results are mod-2
refutations. They are not closure refutations (Lemma A): they say nothing about approximate units, and halving units
escape them the way border units escape rank methods. The two escapes are the archimedean and the 2-adic "points at
infinity" of the same system, matching literature T3's "conditions at the real place and at the place 2".

## 5. Pattern certificates (a2-nag Lemma 2) as NS with beta-free cofactors [proved; definitions]

For a relaxation A(alpha) beta = t (linear in the readouts beta), a pattern certificate of alpha-degree e is
y(alpha) with deg y_r <= e, y^T A = 0 identically and y^T t = 1. It is an NS refutation with cofactors free of beta.
If e = 0 the residual is bounded below by 1/|y|^2 for every alpha (a closure refutation of the relaxation: no border).
By Lemma B, a relaxation residual decaying like |alpha|^(-a) forces e >= a. betalin.py computes the least e.

Example (mw-sym4's border skeleton t27_mid3_m2_s1, chart n): a2-nag's certificate uses (1+a1)(1+a2)(1+a3), degree 3,
but the least e is 1: the rows u, a1 u, a2 u, a1 a2 u (u = b8 + a4 b9) sum with constant coefficients to
(1+a1)(1+a2) u, so row27 = (1+a3)(that sum) - 1 and y has cofactor (1+a3). The curve a1 = a4 = b9 = b5 = b7 = 0,
a2 = 1, b8 = eps, 1 + a3 = 1/(2 eps), b6 = -a3 eps makes every row 0 except the two rows equal to u, which are eps:
residual ~ 1/|alpha|, a = 1 = e. Lemma B is tight there [derived here; e computed].

## 6. Lemma A in practice: closure polynomials for relaxations [proved for the reduction; computed mod p]

For a relaxation A(alpha) beta = t the image {A(alpha) beta} is a cone (scale beta), so a closure polynomial may be
taken homogeneous of some degree k in the row values z: Q(A(alpha) beta) = 0 identically and Q(t) = 1. Degree k = 1
is exactly a pattern certificate with constant cofactors (e = 0). Because the relaxation is the image of the unit
equations under the C-linear frame-sum map Phi (a2-nag Lemma 1), Q o Phi vanishes on everything the skeleton outputs
and is 1 at the target: a closure refutation of the skeleton, which therefore has no border unit (Lemma A(ii)).
closure.py solves the linear system for the coefficients of Q (mod p), extracts one solution and checks
Q(A(alpha) beta) = 0 at 20 random points and Q(t) = 1. Calibration: the two known borders (t27_mid3_m2_s1 chart n:
explicit curve of section 5; cal_plain4F_s1 chart n: the u-chart solution with first coefficient 0 is a point at
infinity of the normalised chart) give no Q, as Lemma A requires.

## 7. Session 2 (23:01 EDT onwards): border curves, and border existence is decidable per skeleton

**Correction to section 1.** Lemma A(ii) is not new. It is the standard fact behind geometric complexity theory:
a lower bound proved by test polynomials is a lower bound for the Zariski closure ("border-C"), and the closure of a
constructible set consists of its limits [read: first-hand, Grochow, arXiv:1304.6333v1 (Computational Complexity
2015), Definition 2.10 (separating module) and section 3.2, citing Mumford Thm 2.33 and BCS97 section 20.6]. What
is ours is only the dictionary: our certificates are separating polynomials for the skeleton images.

**Lemma T (border curves)** [proved]. Let A(alpha) beta = t be a relaxation (or any system phi(c) = T). If Laurent
polynomials alpha(eps), beta(eps) satisfy A(alpha(eps)) beta(eps) = t + O(eps), then t is in the closure of the
image, so no closure polynomial Q of any degree exists (Lemma A). In particular, for pole orders (0, 1): if
A(a0) b0 = 0 and A'(a0)[d] b0 + A(a0) b1 = t, the curve alpha = a0 + eps d, beta = b0/eps + b1 is a border.
*Proof.* Expand A(a0 + eps d)(b0/eps + b1) = A(a0) b0/eps + (A'(a0)[d] b0 + A(a0) b1) + O(eps). For a general
curve, Q(t) = lim Q(A(alpha(eps)) beta(eps)) = 0 by continuity. QED.

**Proposition D (border existence is decidable for each skeleton, with a certificate either way)** [proved from
standard facts cited from memory]. For a skeleton with polynomial map phi and target T, exactly one holds:
(a) some Q with Q o phi = 0 and Q(T) = 1 (a closure refutation; finite linear algebra once deg Q is fixed);
(b) some Laurent-polynomial curve c(eps) with phi(c(eps)) = T + O(eps) (a finite polynomial system once the pole
orders are fixed, decided by a Groebner basis).
Searching (a) by increasing degree and (b) by increasing pole order, one search must stop. So "can a proof of the
closure kind (every Theorem-R certificate is one) exclude this skeleton?" is decidable, and either answer has a
checkable certificate. A shape has finitely many skeletons, and a border unit of the shape is a border unit of one of
them (pigeonhole on a convergent sequence), so the same holds for each shape.
*Proof.* (a) and (b) exclude each other by Lemma T. If (a) fails, T lies in the Zariski closure of the image (Lemma
A(ii)). Let G be the graph of phi and G' its closure in P^n x C^N; G is dense and open in the irreducible G', and some
point P of G' lies over T. Through P there is an irreducible curve in G' meeting G [memory: Mumford, Red Book I.6 /
Abelian Varieties section 6 lemma]. A branch of it at P has a Puiseux parametrisation; restricted to G and projected
to C^n it is c(eps) in C((eps^(1/N)))^n with phi(c(eps)) -> T. Put eps = s^N. The coefficients of phi(c(s)) up to s^0
use finitely many coefficients of c, so truncating c gives a Laurent polynomial with phi(c(s)) = T + O(s). QED.

So, for a fixed shape, "a rank-type impossibility proof exists" is a decidable statement with explicit witnesses on
both sides: the polynomials Q, or one border curve. It cannot be independent of any theory that proves basic algebra.
The cost of the decision is the cost of Groebner bases, which is the only obstacle.

**Lemma T' (borders pass down to relaxations, not up)** [proved]. Let Psi be a C-linear map with Psi(full unit
equations) = relaxation, as for a2-nag's frame-sum map and for each of nag-rerun's eight character maps (charrelax.py),
the readouts of the relaxation being linear images of the full readouts. If (alpha(eps), B(eps)) is a border curve of
the full system with B of pole order <= q, then (alpha(eps), Psi B(eps)) is a border curve of every such relaxation,
all at once (shared alpha), with readout pole order <= q. *Proof.* Apply Psi to A_full(alpha) B - T = O(eps). QED.
So "no first-order border curve of the combined relaxations" excludes first-order border curves of the full system;
a border of one relaxation says nothing about the full system.

**Computations (session 2)** [computed mod p at two primes unless marked]: first-order (pa, pb) = (0, 1) border
systems (bcurve.py) are consistent for the every-coefficient-free relaxations of g27_Fany4_s5, g27_Fany4_s6 and
g27_Fp2b_m3_s5, so by Lemma T no closure polynomial of any degree refutes those relaxations; f27_FAmid3_m3_s16: the
same, numerically only (tangent_num.py, residual^2 1e-29). Calibration: consistent on a toy tangent border, GB = [1] on
a toy without border and on all five relaxations that have a closure polynomial (as Lemma T forces).
