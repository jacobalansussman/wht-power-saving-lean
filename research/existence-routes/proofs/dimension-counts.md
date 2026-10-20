<!-- header: existence-routes -->
# Counting dimensions (Lemmas A to F)

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** Argued by hand, not checked; two AI skeptics re-derived it and reproduced the numbers, one with exact ranks modulo two primes. Their correction: what is proved is that counting in the whole space cannot prove existence; whether a count can exclude a shape is open.

---

# ex-dimension: proofs of the [proved] claims in NOTE.md

14 October 2026. Agent ex-dimension (existence round of the meta work). Nobody has read this. No Lean. Nothing posted.
Notation and the normal form are those of pf-algebraic (`../../proofs/algebraic/PROOF.md`, Lemma 6, Lemma 7, Theorem 8,
Corollary 9, Lemma 10), which I use as established and do not re-prove.

## 0. Setting

- G = G_m is the finite engine group (Lemma 6: scalars {+-1, +-i}, |G| = 4 * 4^m * |Sp(2m,2)|), B its free subgroup
  (|B| = 4 * 4^m * |P|, P the Siegel parabolic), L = [G:B] = number of Lagrangians = prod_{i=1..m} (2^i + 1)
  (15, 135, 2295 at m = 2, 3, 4). n_k = number of Lagrangians at distance k from L0 = 2^(k(k+1)/2) [m choose k]_2.
- A shape S is (W; step ranks r_1..r_n; target ranks k_1..k_W). Its parameter space is
  Par_S = C[B]^(slots), slots = sum_{j=1..n} (W + j - 1) + W (W + n), coordinates a_{j,z}, f_{t,z} (Theorem 8).
- The program map Phi_S : Par_S -> M_W(C[G]) sends parameters to the transfer matrix (T_{tl}): output t equals
  sum_l T_{tl} e_l. It is polynomial. The target is T* = diag(W_{k_1}, ..., W_{k_W}) (W_m = K the kernel). The unit
  variety is X_S = Phi_S^{-1}(T*). By Theorem 8, X_S(Q-bar) != empty iff a unit of shape S exists with algebraic
  mixing coefficients; Theorem 8 itself is stated over Q, and Lemma F below says why the field does not matter for
  emptiness over C versus Q-bar.
- For a representation rho of G (dimension d) and a B-invariant subspace U (dimension u; U = V allowed),
  Corollary 9 / Lemma 10 give the relaxation Phi_{S,rho} : A^(slots) -> (C_G iota_U)^(W x W), A = rho(C[B]),
  C_G = rho(C[G]) (the algebra spanned by all letters), iota_U the inclusion of U. `dimcount.py` implements it.

## Lemma A (the image is a two-sided module; orbit test) [proved]

(i) For alpha in M_W(A) and X in M_W(A_U) (A_U = {iota_U^* a iota_U : a in A}, which acts on U because U is
B-invariant), alpha . Phi_{S,rho}(Par) . X is contained in Phi_{S,rho}(Par).
(ii) Let H = GL_W(A) x GL_W(A_U) act on (C_G iota_U)^(W x W) by (alpha, X).T = alpha T X^{-1}. If T* lies in the
Zariski closure of the image of Phi_{S,rho}, then
    o := dim span{ alpha T* + T* X : alpha in M_W(A), X in M_W(A_U) }  <=  r := generic rank of the Jacobian of Phi_{S,rho}.
So "o > r" excludes the shape S in rho, hence (Corollary 9) excludes every unit of shape S, exact or border.

*Proof.* (i) Post-mixing: sum_s alpha_{ts} O_s = sum_z (sum_s alpha_{ts} f_{s,z}) z, again a choice of final
coefficients in A. Pre-mixing: replace the input e_l by e'_l = sum_k a_{lk} e_k (a_{lk} in A; free). Because a maps U
into U, a iota_U = iota_U (iota_U^* a iota_U), so e' = e . X with X the block matrix (iota_U^* a_{lk} iota_U) acting
on the right. Every step of the program is a left multiplication or a sum, which commute with right multiplication,
so running the same coefficients on e' produces every content times X; the coefficients of e'_l in the first step
are A-combinations of the e_k, so this is again a normal-form program, with transfer matrix T X.
(ii) The closure Z of the image is irreducible (image of an irreducible variety) and, by (i) and continuity, stable
under H. If T* is in Z, the orbit H.T* lies in Z, so dim(H.T*) <= dim Z. dim Z equals the generic rank of the
Jacobian of a polynomial map from an affine space (generic smoothness in characteristic 0). H.T* is a smooth
variety whose tangent space at T* is the image of the differential of (alpha, X) -> alpha T* X^{-1} at the identity,
which is {alpha T* - T* X}, of dimension o. []

Remark. Lemma A is a rank method (it uses only closed conditions), so by barrier B1 it cannot exclude a shape that
has a border scheme in rho. In every case computed it was weaker than Theorem R (NOTE.md section 3).

## Lemma B (what dimension counting can prove about existence) [proved; standard algebraic geometry]

Let Phi : C^P -> Y be a polynomial map into an irreducible affine variety Y containing its image, and T* in Y.
(a) If Phi is dominant (generic Jacobian rank = dim Y), then T* is in the closure of the image: a border scheme.
(b) If moreover the group H of Lemma A acts on Y with H.T* open in Y, then T* is in the image: an exact solution.
(c) If the Jacobian of Phi has rank dim Y at a point p with Phi(p) = T*, then the fibre X = Phi^{-1}(T*) is smooth
    of dimension P - dim Y at p and Phi is dominant. A smooth point of X alone (without the rank condition) does
    not make Phi dominant.
*Proof.* (a) Chevalley: the image is constructible and dense, so its closure is Y. (b) The image contains a dense
open subset O1 of Y (Chevalley). H.T* is open and dense (Y irreducible). Two dense open subsets meet, so some
h.T* is in the image; the image is H-stable (Lemma A (i), with h invertible), so T* is in the image.
(c) Implicit function theorem / submersion at p; dominance because the image contains a neighbourhood of T* in Y. []

Without (b) or another argument, a count can at best prove (a), and a border scheme is worth its full value only if
it is formally polynomial (meta-border EXACTIFY.md Theorem X), which is a property of the free-symbol expansion and
not of any finite-dimensional representation.

## Lemma C (the exact count, in units of |B|) [proved]

In the faithful setting (rho = the regular representation, i.e. Theorem 8 itself, U = everything):
- P = slots * |B|; the transfer matrix lives in M_W(C[G]), of dimension E = W^2 |G| = W^2 L |B|.
- Every step has a gauge of dimension at least |B| / n_{r_j}: replacing c_j by h c_j with h in C[B cap W_r^{-1} B W_r]
  invertible changes u_j into (W_r h W_r^{-1}) u_j, an invertible element of C[B] times u_j, which later
  coefficients absorb. |B cap W_r^{-1} B W_r| = |B| / n_r because P acts transitively on the n_r Lagrangians at
  distance r from L0 (Witt; pf-algebraic Lemma 1) and B = phi^{-1}(P) contains ker(phi) cap G (Lemma 6).
  So r <= P - sum_j |B|/n_{r_j}.
- For the target T* = diag(W_{k_t}): o = sum_{t,s} (2 |B| - |B W_{k_s} cap W_{k_t} B|), and
  |B W_k cap W_k B| = |(B cap W_k B W_k^{-1}) W_k| = |B| / n_k, while for k_s != k_t the two cosets lie in different
  (B,B) double cosets (Lemma 7: the rank is a double-coset invariant), so the intersection is empty.
All numbers in `out_exact_counts.txt` follow from these formulas (`exact_counts.py`). []

**Corollary C1 (the naive count fails calibration).** For every shape in the table, including the plain schedules,
P - E < 0 (for example plain 1 in 2 at block 2: P - E = (6 - 15)|B| = -9|B|). Plain schedules exist, so X_S is
not empty although its "expected dimension" is negative: X_S is an excess intersection, the W^2 |G| equations are
far from independent on it, and "expected dimension < 0" can never be used to predict that a shape is empty.

**Corollary C2 (no dominance near plain cost).** Phi_S can be dominant onto M_W(C[G]) only if slots >= W^2 L. With
n = cW steps, slots = W^2 (1 + 2c + c^2/2) - cW/2, so dominance needs c >= sqrt(2 + 2L) - 2: c >= 3.66 at block 2,
14.49 at block 3, 65.8 at block 4, against c = m for the plain schedule. So at the exact level no unit at or near
plain cost (let alone below it) is dominant, and Lemma B (a) never applies.

**Corollary C3 (no open orbit).** o <= 2 W^2 |B| < W^2 L |B| = E, since L >= 15. So H.T* is never open in
M_W(C[G]), and Lemma B (b) never applies at the exact level, whatever the shape.

Corollaries C1 to C3 together: in the faithful setting, a dimension count can neither exclude nor establish any
unit. The only output of counting that survives is Lemma A, which is a (weak) rank method.

## Lemma D (the everything-free game collapses at each fixed f) [proved, given an upheld input]

Input (wave 17 g-hypotheses, upheld; reps track): for the free group F_XW generated by the engine's free letters,
the chi-type chain map row ^= row<<1 and Toffoli, V^(x)f is irreducible for every f >= 2 tested.
Claim: then for each such f, the relaxation of any shape to the single representation V^(x)f is solvable with ZERO
paid moves; the count there is dominant at cost 0 and carries no information.
*Proof.* Burnside: the algebra spanned by an irreducible matrix group over C is the full matrix algebra, so
A = End(V^(x)f) contains K^(x)f. []
Consequence: in the XW game every obstruction, and every count, must live at the formal level (coefficients attached
to words, the same for every f). There the parameter space is C[F_XW] for an infinite group F_XW (chain maps on f
columns for all f), so "expected dimension" is undefined unless the word supports are fixed in advance; with fixed
supports it becomes a finite count (number of slots against number of distinct formal products), which depends on
the chosen support and predicts nothing about the shape as such.

## Lemma F (fields) [proved; standard]

X_S is defined by polynomials with coefficients in Q(i) (the letters have entries in Q(i), Lemma 6). If X_S has a
C-point it has a Q-bar-point (Hilbert's Nullstellensatz: a Q(i)-ideal that is proper stays proper over Q-bar).
A Q-bar-point is a unit with algebraic coefficients, not yet a Lean-admissible one: the engine mixes with rational
coefficients and has i as a free letter, so Q(i)-points are what it can run; descending from Q-bar to Q(i) is an
arithmetic question that no dimension count addresses: a positive-dimensional variety over Q(i) with Q-bar-points can
have no Q(i)-point (for instance x^2 = 3 in the (x, y)-plane, since sqrt(3) is not in Q(i)).
