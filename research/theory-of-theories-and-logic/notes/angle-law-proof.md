> **Working note, kept as written.** An AI agent wrote this during the second abstract wave of 14 October 2026 (the stability track). Private paths and machine names
> were removed and nothing else was changed. No person has checked it. Where it disagrees with this section's
> [README](../README.md), the README wins, because it follows the later skeptic verdicts. Two skeptics later upheld Theorem F, the angle law and Theorem H, and showed more: every hook has right-angle edges,
> and the depth-3 count 68,320 of the near miss is exact over the rationals. "Below 2 at every block" holds only
> up to block 40 (exact computation); beyond that only "below 8/3" is argued.

# a2-vi-stability: the Frobenius ceiling and the angle law (14 October 2026)

One agent; no skeptic has read this. Nothing posted. Labels: [proved] (proof below), [computed] (script and output
named), [argued], [guess], [read]. Notation: DSp(2m,2) is the frame graph (dual polar graph), N frames, degree
k = 2(2^m − 1), intersection numbers b_j = 2^(j+1)(2^(m−j) − 1), c_j = 2^j − 1, eigenvalues θ_i = 2(2^(m−i) − 1) − (2^i − 1)
(i = 0..m), cosine sequences ω_i(j) (ω_i(0) = 1, ω_i(1) = θ_i/k, three-term recurrence). The eigenvalue list is the
classical-parameter formula for dual polar graphs (Brouwer, Cohen and Neumaier, *Distance-Regular Graphs*, 1989,
Corollary 8.4.2; cited from memory, not re-read today); `frob_ceiling.py` checks it for m ≤ 40 by the closing equation
of the recurrence (exact rationals), together with integrality of every multiplicity.

## 1. Theorem F (Frobenius ceiling for every frame family) [proved]

Let U ↦ Y_U be ANY assignment of an r-dimensional subspace of R^D or C^D to every frame U (no symmetry assumed). For a
pair (U, V) put δ(U,V) = dim(Y_U + Y_V) − r (the rank excess, what Theorem R certificates count) and
φ(U,V) = Σ sin²(principal angles between Y_U and Y_V) = r − ‖Q_U* Q_V‖_F² (Q orthonormal bases). Let ā_d, φ̄_d be the
means over all ordered pairs at distance d. Then

    φ̄_d = Σ_{i=1..m} c_i (1 − ω_i(d))   with every c_i ≥ 0,

hence φ̄_j / φ̄_1 ≤ T_m(j) := max_{i ≥ 1} (1 − ω_i(j)) / (1 − ω_i(1)) for every j.

*Proof.* K(U,V) = ‖Q_U* Q_V‖_F² = tr(P_U P_V) is the Gram matrix of the orthogonal projections P_U in the real inner
product space of Hermitian matrices, so K is positive semidefinite, and K(U,U) = r. Let A_d be the distance-d matrix
and E_0..E_m the primitive idempotents of the Bose–Mesner algebra; A_d = Σ_i k_d ω_i(d) E_i. Then
K̄_d = ⟨K, A_d⟩/(N k_d) = Σ_i ω_i(d) tr(K E_i)/N and tr(K E_i) ≥ 0 (both PSD). Put c_i = tr(K E_i)/N; since ω_i(0) = 1,
φ̄_d = K̄_0 − K̄_d = Σ_{i≥1} c_i (1 − ω_i(d)), and 1 − ω_i(1) > 0 for i ≥ 1. A ratio of nonnegative combinations is at most
the largest termwise ratio. ∎

The same argument applies verbatim to mean ‖f(U) − f(V)‖² for any vector-valued f on frames (Delsarte's linear-programming
method; pf-induction Theorem 6 used it for 0/1 functions, i.e. cuts).

## 2. The i = 1 term, for every m [proved]

**Lemma 2.** ω_1(j) = ((2^m + 1)/2^j − 2)/(2^m − 1) for 0 ≤ j ≤ m, for every m. Consequently
(1 − ω_1(j))/(1 − ω_1(1)) = 2 − 2^(1−j) exactly, independent of m, and ω_1(m) = −2^(−m).

*Proof.* Write x = 2^j, y = 2^m, u = (y + 1)/x, so w(j) := ω_1(j) = (u − 2)/(y − 1), w(j − 1) = (2u − 2)/(y − 1),
w(j + 1) = (u/2 − 2)/(y − 1). Here k = 2y − 2, θ_1 = y − 3, b_j = 2y − 2x, c_j = x − 1, a_j = k − b_j − c_j. Multiplying the
recurrence c_j w(j−1) + a_j w(j) + b_j w(j+1) = θ_1 w(j) by (y − 1), the constant terms give −2(c_j + a_j + b_j) = −2k on the
left and −2θ_1 on the right, and the u-terms give u(2c_j + a_j + b_j/2) = u(k + c_j − b_j/2) on the left and uθ_1 on the
right. So the identity is u(k − θ_1 + c_j − b_j/2) = 2(k − θ_1). Now k − θ_1 = y + 1 and c_j − b_j/2 = 2x − 1 − y, so the
left side is u·2x = 2(y + 1), which is the right side. The same computation covers j = m (where b_m = 0 is already the
value of 2y − 2x), and w(0) = 1, w(1) = (y − 3)/(2y − 2) = θ_1/k. Finally
(1 − w(j))/(1 − w(1)) = (y + 1 − u)·(2/(y + 1)) = 2(1 − 1/x). ∎ (Machine cross-check on a rational grid:
`omega1_symbolic.py`.)

So the eigenvalue 2^m − 3 alone already gives T_m(j) ≥ 2 − 2^(1−j) at every m.

## 3. The other eigenvalues

- [computed, exact rationals, m = 1..40, `frob_ceiling.py` → `out_frob_ceiling.txt`] T_m(j) = 2 − 2^(1−j) for every j,
  attained only at i = 1; ω_i(m) = (−1)^i 2^(−i(m+1−i)) for every i; the termwise top ratios (1 − ω_i(m))/(1 − ω_i(1))
  tend to 2^i/(2^i − 1) (2, 4/3, 8/7, 16/15, ...) as m grows.
- [proved, every m, given the eigenvalue list] T_m(j) < 8/3. For i ≥ 2: |ω_i(j)| ≤ 1 (ω_i(j) = ⟨E_i e_x, E_i e_y⟩/‖E_i e_x‖²
  for any x, y at distance j, and ‖E_i e_x‖ = ‖E_i e_y‖, so Cauchy–Schwarz), so 1 − ω_i(j) ≤ 2; and θ_i ≤ θ_2 = 2^(m−1) − 5, so
  1 − ω_i(1) ≥ 1 − θ_2/k = 3(2^(m−1) + 1)/(2^(m+1) − 2) > 3/4. Hence those ratios are < 8/3, and the i = 1 ratio is < 2.
- [argued] T_m(j) = 2 − 2^(1−j) for every m: it needs only ω_i(j) ≥ 1 − (2 − 2^(1−j))(1 − θ_i/k) for i ≥ 2, which the exact
  table confirms with a wide margin for m ≤ 40 and which the closed form (b) would give in general; (b) is the standard
  q-Krawtchouk evaluation at the last distance (from memory, not re-derived). pf-induction listed exactly this as the
  missing step for its κ_m at m > 12; the exact check now reaches m = 40.

## 4. Corollary A (the angle law) [proved]

For pairs at distance d let s_d = φ̄_d/ā_d, the mean of sin² over the nonzero principal angles (weighted by count). Then
for every frame family, hence for every engine-game Theorem-R certificate (reps Theorem 1 and the Averaging Lemma make
every certificate a frame family with profile ψ_d = ā_d):

    ψ_j/ψ_1 = (φ̄_j/φ̄_1) · (s_1/s_j) ≤ T_m(j) · s_1/s_j ≤ T_m(j)/s_j.

The first factor is spectral and capped by T_m(j) = 2 − 2^(1−j) < 2 (m ≤ 40; < 8/3 for every m). **Every move per array
that a certificate proves above 2 comes from its angle ratio s_1/s_j**: adjacent frames' subspaces must meet at larger
angles than distant frames' subspaces. Consequences:
- (i) Coordinate families (all angles 0 or 90 degrees, s ≡ 1) are capped by T_m < 2 (pf-induction Theorem 6(i); now for
  m ≤ 40 exactly and < 8/3 for every m).
- (ii) A certificate with θ-shadow ≥ 0.9 at block m needs ψ_m/ψ_1 ≥ m^0.9, so s_m ≤ T_m s_1 m^(−0.9) < 2 m^(−0.9): the
  nonzero angles between opposite frames' subspaces must be small on average (sin² below 0.35 at m = 7, 0.25 at m = 10,
  0.14 at m = 20), although there are as many of them as the rank says.
- (iii) Scale invariance (linear form). If every Y_U is the graph of a map A_U over one common r-dimensional Y_0, then
  δ(U,V) = rank(A_U − A_V) exactly, and replacing A_U by tA_U changes no rank. The vector form of Theorem F bounds
  mean ‖A_U − A_V‖_F² at distance j by T_m(j) times its value on edges. So the angle law says, scale-free: **the mean
  squared nonzero singular value of A_U − A_V must be at least m^0.9/T_m times smaller on opposite pairs than on
  edges.** Opposite pairs must differ by matrices of full rank but small Frobenius mass per rank unit. Rank counts tiny
  singular values as whole units, and only that lets a certificate escape the factor 2. This is where
  characteristic 0 differs from the mod-2 family (Proposition G), where there are no small singular values. It is a
  border-rank-like phenomenon on the certificate side [argued: the reading; the identity itself is proved].

## 5. Remarks

- Lifts are not obstructed by the angle law alone [argued]: the rational subspaces reducing to a given F_2-subspace are
  dense in the real Grassmannian (the condition is 2-adic), so real angles can be chosen almost freely without
  changing the mod-2 data. The law says where a θ-consistent lift must live (many small singular values on opposite
  pairs). It does not say whether one exists.
- Equality in Theorem F at j = m forces all of the spectral weight onto the eigenvalue 2^m − 3 (c_i = 0 for i ≥ 2). The
  block-3 hook does exactly this: its Frobenius ratio is 7/4 = T_3 to machine precision (`out_cut_frob.txt`, module 8).

## 6. Theorem H (every hook saturates Theorem F) [proved]

Let V = C^n, n = 2^m, 1 ≤ k < n, X = V ⊗ Λ^k V with orthonormal basis e_x ⊗ e_T (T a k-subset), and the hook
N = span{e_x ⊗ e_T : x ∈ T} (dim N = k·C(n,k)). For every unitary g on V, with G = g ⊗ Λ^k g:

    ‖(1 − Π_N) G Π_N‖_F² = C(n−2, k−1) · (n − Σ_{x,y} |g_xy|⁴).

*Proof.* G is unitary, so ‖(1 − Π)GΠ‖² = dim N − ‖ΠGΠ‖², and
‖ΠGΠ‖² = Σ_{x∈T, y∈S} |g_xy|² |det g[T,S]|² = Σ_{x,y} |g_xy|² Σ_{T∋x, S∋y} |det g[T,S]|².
For a fixed k-set T the rows g_T are orthonormal, so by Cauchy–Binet Σ_S |det g[T,S]|² = det(g_T g_T*) = 1, and deleting
column y gives Σ_{S∌y} |det g[T,S]|² = det(I − c c*) = 1 − ‖c‖² with c = g_T[:, y]; hence Σ_{S∋y} |det g[T,S]|² =
Σ_{x'∈T} |g_x'y|². Summing over T ∋ x, and using Σ_{x'} |g_x'y|² = 1 (columns of a unitary),
Σ_{T∋x, S∋y} |det g[T,S]|² = C(n−1,k−1)|g_xy|² + C(n−2,k−2)(1 − |g_xy|²) = C(n−2,k−2) + C(n−2,k−1)|g_xy|².
So ‖ΠGΠ‖² = n·C(n−2,k−2) + C(n−2,k−1)·Σ|g_xy|⁴, while dim N = k·C(n,k) = n·C(n−1,k−1) = n·(C(n−2,k−2) + C(n−2,k−1)). ∎

For the engine's depth-j element g_j = dir(e_1)···dir(e_j), dir z = ((1+i)/2) I + ((1−i)/2) X^z (proofs/induction Fact A),
g_j is a tensor product of j copies of a 2 × 2 matrix whose four entries have |.|² = 1/2 with the identity on the other
m − j bits, so Σ|g_xy|⁴ = 1^j · 2^(m−j). Hence

    φ_j = C(n−2, k−1) · 2^m · (1 − 2^(−j)),   φ_j/φ_1 = 2 − 2^(1−j)   for every m and every k:

equality in Theorem F at every depth. By the strict gap of section 3 (m ≤ 40), equality forces all spectral weight onto
the eigenvalue 2^m − 3, so the Gram kernel of the projections P_U lies in span(E_0, E_1). That span is the image of the
point–frame incidence map C[points] → C[frames] [argued for every m: the common constituents of the two permutation
modules; checked by dimensions at m = 3 (63 = 1 + 27 + 35 against 135 = 1 + 35 + 84 + 15) and m = 4 (255 = 1 + 119 + 135
against 2295 = 1 + 51 + 135 + 918 + 1190)]. Therefore every hook family is **point-additive**:
P_U = c·I + Σ_{p ∈ U} M_p for fixed Hermitian matrices M_p, one per point of PG(2m−1, 2).

Checks [computed, `hookfrob/hookfrob.py` → `hookfrob/out_hookfrob_small.txt`, dense SVD]: (m, k) = (2, 1), (3, 2), (3, 4),
(4, 2): φ equals the closed form to 1e-6 at every depth; ranks 80, 144, 224 at (3, 4) reproduce the hook's cap and demand.
At every computed case ψ_1 = φ_1 (all angles at edges are 90 degrees); at block 4, k = 8 the closed form gives φ_1 = 27456,
equal to the certified cap, so the same holds there.

**Consequence for the hook family.** Its value is (2 − 2^(1−m)) · s_1/s_m with s_1 = 1 in every computed case, and b56's
ceiling 4(n−1)/(n+2) becomes: the hook's angle ratio is at most 2n/(n+2) < 2. The block-4 k = 8 hook, with certified ranks
27456, 48576, 68320, 86592 and φ = 27456 × (1, 3/2, 7/4, 15/8), has s = 1, 0.84783, 0.70328, 0.59451.

**The block-4 near miss (ms-fano).** A 2 in 7 whose arrays both end privately must deliver two depth-3 terms in five moves,
which a certificate excludes when ψ_3/ψ_1 > 5/2. With a full depth-3 budget (φ_3 = (7/4) φ_1, the most Theorem F allows)
and 90-degree edges, this reads s_3 < 7/10. The k = 8 hook has s_3 = 48048/68320 = 0.703279...: it would need rank 68640
at depth 3 for the same Frobenius mass and has 68320, i.e. 320 per array, which is ms-fano's slack of 640 for two arrays.
Its budget and its edges have no slack left (proved above). Any certificate that closes this branch must spread its
depth-3 Frobenius mass over more directions (smaller mean sin²) than the hook does.
