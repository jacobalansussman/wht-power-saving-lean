<!-- header: existence-routes -->
# Limits in the number of columns and of address bits

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** Argued by hand, not checked; one AI skeptic checked the theorems step by step and found two wrong calibration reasons (the limit sees one residue of f modulo 4, not every irreducible at once; px-t514's program is not a unit at any f) and one fixable gap (a global 8th-root phase). The conclusions stand.

---

# ex2-limits: limit objects and transfer — statements and proofs

ex2-limits (track D of `../BRIEF.md`), 14 October 2026, Mac clock 19:37 onwards. One agent. **No skeptic has read
this.** Nothing posted. Labels: [proved] = every step is here or is a cited standard fact; [computed] = program and
output named; [argued]; [guess]; [read] = read today at the place cited.

## 0. Setting (as in `../../moonshot/mw2-xwlayers/PROOF.md` §0)

Block m, f digit columns, H_f = C^N, N = N_f = 2^(mf). A *letter* is a family (b_f)_f of N × N matrices. Paid letters
D_z = dir(z)^(⊗f), dir(z) = a·I + b·X^z, a = (1+i)/2, b = (1−i)/2. Target K = (dir(e_1)···dir(e_m))^(⊗f). A unit
"W in M" is a straight-line program u_k = D_k c_k, c_k = Σ_(i<k) B_ki u_i + γ_k, with outputs
K e_l = Σ_k A^(l)_k u_k + A^(l)_0, where every B, γ, A is a linear combination, **with coefficients independent of
f**, of products of free letters; it must be right at every f. Its *residual* is the finite list of family-algebra
elements R^(l) = (output) − K e_l (one per output and input slot); it is a unit iff every R^(l)_f = 0.

Two facts used throughout [proved, elementary]:
- (F1) dir(z) is the Fourier multiplier ξ ↦ i^[ξ·z] ([t] ∈ {0,1} is t mod 2): on the character χ_ξ(x) = (−1)^(ξ·x),
  X^z χ_ξ = (−1)^(ξ·z) χ_ξ, so dir(z) χ_ξ = (a + b(−1)^(ξ·z)) χ_ξ and a + b = 1, a − b = i. Hence all paid letters
  commute, and K is the multiplier i^|ξ| (|ξ| = Hamming weight), per column.
- (F2) Normalised Hilbert–Schmidt residual: ρ_x(f) := tr(x_f^* x_f)/N_f for a family x.

## 1. The limit structures, and what a "limit unit" is in each

| | Limit | Limit unit | Transfer back |
|---|---|---|---|
| L1 | f → ∞, tracial ultralimit τ_U(y) = lim_U tr(y_f)/N_f on the family algebra 𝒜 (mw2-xwlayers §2), von Neumann algebra M_U | a program with C-coefficients whose identities hold in L²(M_U), i.e. lim_U ρ_R(f) = 0 | §2 (E, Wd: exact on a residue class of f); §3 (XQ, Xchi: exact on everything, Theorem 9; XW: exact through one paid letter, Theorems 8, 8′; CS with two or more paid letters open) |
| L1′ | L1 with coefficients in M_U or its affiliated operators | soft coefficients | E: none (the free algebra is finite-dimensional, §2.4) |
| L2 | algebraic ultraproduct ∏_U M_(N_f)(C), coefficients in C^U | internal unit | Łoś: f-dependent units, which do not count; Noetherian compactness §5.1 |
| L3 | m → ∞: the Cantor group Z_2^N in each column | a program on L²(Z_2^N)^(⊗f) | §4: restriction is a homomorphism; no finite-move limit unit; infinite ones = coherent families |
| L4 | ultraproduct over primes p of F̄_p | a unit over F̄_p for infinitely many p | Lefschetz principle: exact (ex2-modp's track) |
| L5 | ε → 0: coefficients in C((ε)) (border units) | residual O(ε) | Theorem X (formally polynomial: full value), v-meta-border (window cost) |
| L6 | Hilbert–Schmidt / sofic: f-dependent bounded-degree free letters, error o(1) in ρ | HS-approximate unit | no stability theorem applies; violates "one fixed recipe" and exactness (§5.3) |

## 2. Column-uniform games (E, Wd): the f-limit is the regular representation, and it is rigid

**Hypothesis (CU).** Every letter, free or paid, is a family h^(⊗f) with h a unitary 2^m × 2^m matrix. Let 𝔊 ⊂ U(2^m)
be the group these h generate, Z = 𝔊 ∩ C^× I its scalars. In E, 𝔊 lies in the m-qubit Clifford group (Paulis, affine
maps of F_2^m, Z4 phases and dir(z) = √X^z are all Clifford), so 𝔊 and Z are finite; all entries lie in Q(i),
so Z ⊆ {±1, ±i} and its exponent q divides 4. In Wd, Z is finite whenever the
phase patterns used have bounded order (assumed below; see Remark 2.3).

A word w has family h_w^(⊗f). A program's residual is a formal element x = Σ_i r_i h_i of C[𝔊] (distinct h_i ∈ 𝔊),
with family x_f = Σ_i r_i h_i^(⊗f). Call h_i, h_j *proportional* if h_j = μ h_i, μ ∈ Z; group the indices into
proportionality classes C with representatives h_C, h_i = μ_i h_C (μ_i ∈ Z).

**Theorem 1** [proved]. Under (CU):
- (a) ρ_x(f) = Σ_(i,j) conj(r_i) r_j λ_ij^f with λ_ij = tr(h_i^* h_j)/2^m: an exponential polynomial in f.
- (b) |λ_ij| ≤ 1, with equality iff h_i, h_j are proportional.
- (c) ρ_x(f) = Σ_C |s_C(f)|² + Σ_(i ≁ j) conj(r_i) r_j λ_ij^f with s_C(f) = Σ_(i∈C) r_i μ_i^f; the second sum tends to 0
  geometrically (rate max |λ_ij| < 1).
- (d) Let q be the exponent of Z and U a non-principal ultrafilter on N; U contains exactly one residue class
  f ≡ f_0 (mod q). Then lim_U ρ_x(f) = Σ_C |s_C(f_0)|², and lim_U ρ_x = 0 **iff x_f = 0 for every f ≡ f_0 (mod q)**.
- (e) x_f = 0 for every f iff x = 0 in C[𝔊] (distinct words are independent as families).

*Proof.* (a) tr((h_i^(⊗f))^* h_j^(⊗f)) = tr(h_i^* h_j)^f. (b) h_i^* h_j is unitary; its trace is a sum of 2^m numbers
of modulus 1, of modulus 2^m iff all are equal, i.e. iff h_i^* h_j is scalar. (c) For i, j in one class,
h_i^* h_j = conj(μ_i) μ_j I, so λ_ij^f = conj(μ_i^f) μ_j^f, and the class contributes |s_C(f)|²; the other pairs have
|λ_ij| < 1. (d) μ_i^q = 1, so s_C(f) = s_C(f_0) on the class; the decaying part vanishes in the limit. If all
s_C(f_0) = 0 then x_f = Σ_C s_C(f) h_C^(⊗f) = 0 for every f ≡ f_0. The converse is clear. (e) If x_f = 0 for all f, then
by (c) Σ_C |s_C(f)|² = −(decaying part) → 0, and s_C(f) is periodic in f, so s_C ≡ 0; the μ_i (i ∈ C) are distinct, so
the Vandermonde system Σ_(i∈C) r_i μ_i^f = 0 (f = 0..|C|−1, using periodicity) forces r_i = 0. ∎

**What the limit object is.** By (c)–(d), the GNS space of τ_U is ℓ²(𝔊/Z) twisted by the character μ ↦ μ^(f_0) of Z:
τ_U is the regular (Plancherel) trace of the twisted group algebra. In E this is the finite-dimensional algebra
C[𝔊] ≅ ⊕_π M_(d_π)(C), and τ_U gives each irreducible π the weight d_π²/|𝔊/Z| (twisted). So **the f-limit of the
engine game is the regular representation of a finite group: every irreducible, each with positive weight.** This is
px-t514's Lemma A ("tight at every f ⟺ tight in every irreducible", `../../px/px-t514/STATUS.txt` 15:20) seen from the
limit: the finite-f representations h ↦ h^(⊗f) pick up the irreducibles one by one (Burnside–Brauer), the limit has
them all at once.

**Corollary 2 (no softening in E and Wd)** [proved]. Fix a skeleton (the words of the program; unknown coefficients
c ∈ C^p). Then:
1. a program is a limit unit along U iff it is an exact unit at every f ≡ f_0 (mod q) (Theorem 1(d) applied to each
   residual coordinate);
2. if lim_U ρ_R(c) can be made arbitrarily small with ‖c‖ bounded, an exact unit (on that residue class) exists with
   the same skeleton: each s_C(f_0; c) is a polynomial in c, and a continuous function on a compact set attains its
   infimum;
3. so the only points of the f-limit that are not exact units are limits of coefficient sequences with ‖c‖ → ∞,
   i.e. **border units** (L5). The f-limit adds nothing that L5 does not already contain.

**2.3 Remark (phases of infinite order).** If Z is infinite (Wd with arbitrary constant phases e^(iθ) per column, a
letter whose family is e^(ifθ) I), then (d) fails: x = 1 − e^(iθ)-letter has lim_U ρ = 0 along any U with
e^(ifθ) → 1, but is nonzero at every f ≥ 1 when θ/π is irrational. Such "limit units" are exact only after replacing
a scalar family by its U-limit, which is not a program. They do not transfer. This is the f-limit's own version of a
border phenomenon; it disappears when phases have bounded order.

**2.4 Soft coefficients (L1′).** In E the free letters (Paulis, affine maps, Z4 phases of F_2^m) form a finite group
F ⊂ 𝔊, so the coefficient algebra C[F] is finite-dimensional and already equal to its von Neumann closure in M_U.
Allowing coefficients in M_U changes nothing [proved]. Affiliated (unbounded) coefficients: C[F] is a finite-
dimensional C*-algebra, in which every element affiliated with it already lies in it [proved: a finite-dimensional
von Neumann algebra has no unbounded affiliated operators]. Same in Wd with bounded-order phases.

**2.5 Calibration of §2** [computed, `e1/remote-out/ssh/out_ctl.txt` and `out_analyze_ctl.txt`]: on 14 column-uniform
block-2 words the exact Grams obey N·G_f = (tr h_i^* h_j)^f: G_2∘G_3 = G_5 and G_2∘G_2 = G_4 with error 0
(entrywise products).

## 3. Games that couple columns (XQ, Xchi, XW): transfer matrices, faithfulness, and the one case left open

**Proposition 3 (every element is a translation-invariant matrix product operator)** [proved]. Every letter of the
current games is an MPO with f-independent site tensors and bond dimension bounded independently of f:
column-uniform letters and D_z (bond 1); the chain map x_j ^= x_(j−1) (bond = the previous column's bits); prefix
XOR, its mirror and the parity phases P_r = i^(parity of row r) (bond 2 per row: the running parity); boundary
letters (bond 1, special end tensor). Products multiply bond dimensions and sums add them, so every x ∈ 𝒜 is such an
MPO with some bond dimension D_x. Consequently tr(x_f^* x_f) = ⟨ℓ| E^(f−2) |r⟩ with a transfer matrix E of size D_x²,
and:
- (a) ρ_x(f) is a linear recurrence sequence of order ≤ D_x²; it is bounded (the letters have bounded operator norm),
  so ρ_x(f) = Σ_k c_k β_k^f + (terms with |β| < 1, polynomial coefficients allowed), |β_k| = 1.
- (b) Every-f test from finitely many f: if ρ_x(f) = 0 for D_x² consecutive f ≥ f_1 then ρ_x ≡ 0 for f ≥ f_1
  (Cayley–Hamilton). (The in-house automaton test `../../wave17/xc-legit/every_f.py` is a version of this.)
- (c) lim_U ρ_x is the peripheral part Σ c_k lim_U β_k^f.
- (d) The f-limit is strictly softer than the exact problem iff some fixed x ∈ 𝒜 has ρ_x ≢ 0 but only decaying modes:
  an **asymptotic-only relation**, a family that is wrong at every large f on a fraction of the space that shrinks
  like N^(−δ). In E and Wd there are none (Theorem 1). Here τ is not the regular trace: P_r has τ(P_r) = a ≠ 0
  (the identity P_r = a + b·Zs_r), so Theorem 1 does not apply.

**Proposition 4 (two-term faithfulness, by Reed–Muller distance)** [proved]. Let u be a monomial word of XW (free
letters only). Its permutation part and its Z4 phase are given by F_2-polynomials in the mf input bits whose degree D
is bounded independently of f (affine, chain and prefix maps have degree 1; Toffoli and CS add a factor 2; a sum of f
column phases mod 4 has high bit e_2 of the low bits, degree ≤ 2·(degree of the low bits); the parity phase has degree
1). If u_f ≠ 1 then ρ_(1−u)(f) = 2 − 2 Re τ_f(u) ≥ 2^(1−D).
*Proof.* ρ_(1−u)(f)·N = 2·#{v : u(v) ≠ v} + Σ_(v fixed) |1 − i^(φ(v))|² ≥ 2·#{v : u acts nontrivially at v}. That set is
the union of the supports of the nonzero polynomials among (u(v)_k ⊕ v_k) and the two phase bits. A nonzero F_2-
polynomial of degree ≤ D in n variables is nonzero on at least 2^(n−D) points (minimum distance of the Reed–Muller code
RM(D, n)). ∎
So two distinct monomial families never merge in the f-limit. For combinations of three or more words the argument
breaks: the event "w_i(v) = w_j(v)" is an AND of mf conditions and can be thin (Fix(chain map) has density 4^(−(f−1))
at block 2), so a combination could in principle vanish everywhere except on thin coincidence sets. For letters with
affine address maps it does not:

**Theorem 8 (the f-limit is faithful on the free algebra of an affine monomial game)** [proved here, unchecked;
found at 20:35 after E1 had pointed to it]. Suppose every free letter is monomial, x ↦ (A x + c, i^(φ(x))), with A
affine-linear over F_2 and the two bits of the Z4 phase φ given by F_2-polynomials of degree bounded independently
of f. (This covers every free letter of `xw2lib.py` at block 2 — column-uniform Paulis, S, CZ, CS, P, C01, C10, chain
maps and their prefix inverses, parity phases, boundary versions — and every free letter of E, XQ and Xchi; it
excludes Toffoli.) Let x = Σ_(i≤n) r_i w_i be a fixed combination of free words and U a non-principal ultrafilter.
If lim_U ρ_x(f) = 0, then x_f = 0 for U-most f.

*Proof.* Column v of x_f is Σ_i r_i i^(φ_i(v)) e_(w_i(v)); it is zero iff every class of the partition
"i ~ j iff w_i(v) = w_j(v)" has zero sum Σ r_i i^(φ_i(v)). Every word has an affine address map, so each
F_ij := Fix(w_j^(−1) w_i) = {v : w_i(v) = w_j(v)} is an affine subspace (or empty), and every phase bit is a polynomial
of degree ≤ D (D bounded in f; PROOF Prop. 4). Along U, the codimension of each F_ij is either bounded (call the pair
*fat*) or tends to infinity (*thin*; the empty set counts as thin); finitely many choices, so one pattern holds on a
U-large set of f, and the union of the thin F_ij has density δ_f → 0 along U.
Fat pairs that meet at a point compose: if (i,j) and (j,k) are fat and both hold at v, then F_ik ⊇ F_ij ∩ F_jk has
bounded codimension, so (i,k) is fat. Hence for each v the set T(v) of fat pairs holding at v is an equivalence
relation; call its classes the T-blocks. The full partition at v is obtained from the T-blocks by merging along thin
pairs that hold at v. Write C_T = {v : T(v) = T} and E(T) = ∩_((i,j)∈T) F_ij ⊇ C_T, an affine subspace of bounded
codimension c_T.
Claim: for U-most f, every T-block B has zero sum β_B(v) := Σ_(i∈B) r_i i^(φ_i(v)) at every v ∈ C_T. Suppose not:
β_B(v_1) ≠ 0 at some v_1 ∈ C_T. For each fat pair p ∉ T, F_p ∩ E(T) is an affine subspace not containing v_1, so it
lies in an affine hyperplane {ℓ_p = ℓ_p(v_1) + 1} of E(T); let E′ = {v ∈ E(T) : ℓ_p(v) = ℓ_p(v_1) for all fat p ∉ T}.
Then v_1 ∈ E′ ⊆ C_T and E′ has codimension ≤ k := n² in E(T). Let Y = {v ∈ E′ : φ_i(v) = φ_i(v_1) for i ∈ B}; its
indicator on E′ ≅ F_2^d is a product of ≤ 2n conditions of degree ≤ D, nonzero at v_1, so by the Reed–Muller distance
|Y| ≥ 2^(−2nD)|E′| ≥ 2^(−2nD−k−c_T) N. Every v ∈ Y outside the thin sets has partition exactly the T-blocks and
β_B(v) = β_B(v_1) ≠ 0, so column v is nonzero with ‖column‖² ≥ c_r > 0 (the least nonzero value of |Σ_(i∈B) r_i i^(e_i)|²
over the finitely many phase vectors). Hence ρ_x(f) ≥ c_r (2^(−2nD−k−c_T) − δ_f), which stays bounded below along U:
a contradiction. So the claim holds, and at every v the full classes are unions of zero-sum T-blocks: x_f = 0. ∎

So the free algebra has no asymptotic-only relations, and the f-limit is rigid on it.

**Theorem 8′ (one paid letter)** [proved here, unchecked; 20:45]. Same hypotheses. Let x = x_0 + x_1, where x_0 is a
combination of free words and x_1 = Σ_i r_i w_i a combination of words with exactly one paid letter,
w_i = L_i · D_(z_i) · R_i (letters act left to right; as a matrix w_i = R_i D_(z_i) L_i, L_i and R_i free words). If
lim_U ρ_x(f) = 0, then x_f = 0 for U-most f.
*Proof.* (i) Layers separate: for a free word g and a word w with one paid letter, |tr(g^* w)|/N ≤ (max entry of w) =
2^(−f/2), so ρ_x = ρ_(x_0) + ρ_(x_1) + o(1); both terms are ≥ 0, so both tend to 0 along U, and Theorem 8 handles x_0.
(ii) For x_1: dir(z) e_u = a e_u + b e_(u+z) in each column, so with t ∈ F_2^f and t·z := (t_1 z, ..., t_f z),
  w_i e_v = 2^(−f/2) ω^f Σ_t i^(−|t|) i^(φ_(L_i)(v) + φ_(R_i)(L_i v + t·z_i)) e_(R_i(L_i v + t·z_i)),   ω = e^(iπ/4),
using a^(f−s) b^s = 2^(−f/2) ω^f i^(−s). The target τ_i(v, t) = R_i(L_i v + t·z_i) is affine in (v, t) and injective in
t, and the phase is a bounded-degree polynomial in (v, t) (i^(−|t|) has bits |t| mod 2 and e_2(t) mod 2). Work on
Ω = F_2^(mf) × F_2^f with uniform measure. For i, j let F_ij = {(v, t) : τ_i(v, t) ∈ τ_j(v, F_2^f)}, an affine
subspace of Ω on which the partner t′ = t′_ij(v, t) with τ_j(v, t′) = τ_i(v, t) is an affine function. Fat and thin
are defined by codimension along U as before. Fat coincidence is transitive: F_ik contains the preimage of F_jk
under (v, t) ↦ (v, t′_ij(v, t)) inside F_ij, so codim F_ik ≤ codim F_ij + codim F_jk. At each v, the entries of column v
of x_1 group into classes (pairs (i, t) with the same target), each a union of fat blocks merged along thin
coincidences, with at most one member per word. Since every class has ≤ n members,
  ρ_(x_1)(f) ≥ (c_r / n) · μ_Ω{(v, t) : the class of (i_0, t) at v has nonzero sum}   for each i_0,
with c_r > 0 the least nonzero value of |Σ_(j∈J) r_j i^(e_j)|² (the common factor 2^(−f/2) ω^f cancels against the
normalisation). From here the proof of Theorem 8 runs verbatim on Ω, with "column v" replaced by "the class of
(i_0, t) at v": a nonzero fat-block sum at one point (v_1, t_1) of a cell gives, by the hyperplane cut and the
Reed–Muller count, a set of measure ≥ 2^(−2nD−k−c) in Ω, minus thin sets, on which the class is exactly that fat block
with the same nonzero sum; so ρ_(x_1) stays bounded below, a contradiction. Hence every fat-block sum vanishes
everywhere, every class sums to 0, and x_1 = 0. ∎

What Theorems 8 and 8′ do not cover by themselves: words with **two or more paid letters**, and Toffoli-type letters
(Fix sets are then varieties, and the hyperplane step fails). With two paid letters the paths (t_1, t_2) of one word
can meet each other (D_z · D_z = X^z collapses 4^f paths onto 2^f points), so a class can hold exponentially many terms
of size 2^(−f) and the minimum-value step fails. For Clifford letters the paths recombine in a controlled way:

**Theorem 9 (Clifford games: the f-limit is faithful on the whole family algebra)** [proved here, unchecked; 20:55].
Suppose every letter, free and paid, is at each f a Clifford unitary on the mf qubits (true for every letter of E,
XQ and Xchi: Paulis, affine and chain maps, prefix XOR, S and CZ phases, parity phases P_r = √(Zs_r) up to a scalar,
boundary Cliffords, and D_z = √(X^z)^(⊗f); false for CS and Toffoli). Let x = Σ_(i≤n) r_i w_i be any fixed combination
of words, with any number of paid letters. If lim_U ρ_x(f) = 0, then x_f = 0 for U-most f.
*Proof.* (i) Columns are stabiliser states. A word w is a Clifford unitary, so its Choi state is a stabiliser state:
⟨y|w|v⟩ = 2^(−d_w/2) i^(ψ_w(v, y)) on an affine subspace of (v, y) and 0 elsewhere, with ψ_w a Z4-valued function whose
two bits are F_2-polynomials of degree ≤ 2 (the standard form i^(l·t) (−1)^(Q(t)) of a stabiliser state, with the carry
bit of a Z4 sum of bits of degree 2). Since w is invertible every column is nonzero, so the fibres have the same
dimension d_w(f), and column v is 2^(−d_w/2) Σ_(t∈F_2^(d_w)) i^(ψ_w(v,t)) e_(s_w(v) + A_w t) with s_w affine, A_w injective
and ψ_w of bounded degree in (v, t).
(ii) Layers separate. Two stabiliser columns with supports of dimensions d, d′ have overlap at most
2^(min(d,d′)) · 2^(−(d+d′)/2) = 2^(−|d−d′|/2). Along U each difference d_(w_i) − d_(w_j) either is eventually constant or
tends to ±∞; group the words into layers of bounded differences. Cross terms between layers tend to 0, so
ρ_x = Σ_layers ρ_(x_layer) + o(1), and each layer's ρ tends to 0.
(iii) Within a layer, fix i_0, put d = d_(w_(i0)) and Ω = F_2^(mf) × F_2^d. Word j contributes the coefficient
r_j 2^((d − d_j)/2) i^(ψ_j) at its points, a finite set of possible values along U, and each word meets each target at
most once (A_j injective). For (v, t) ∈ Ω let F_(i0 j) = {(v, t) : s_(i0)(v) + A_(i0) t ∈ s_j(v) + range A_j}, an affine
subspace with an affine partner map, as in Theorem 8′. The proof of Theorem 8′ (ii) now applies verbatim: transitivity
of fat coincidences, a cell, the hyperplane cut, the Reed–Muller count for the bounded-degree phases, and the bound
ρ_layer ≥ (c/n) · 2^(d_(i0) − d) · μ_Ω(bad) with c > 0 the least nonzero class sum. So every class sums to 0. ∎

So in E, XQ and Xchi the f-limit is rigid on everything, and a limit unit with fixed coefficients is an exact unit at
U-most f. In XW the rigidity is proved for every word with at most one paid letter (Theorems 8, 8′, which allow CS);
words that combine CS with two or more paid letters are the one case left. There the coset sums of a word's paths are
exponential sums of a degree-4 phase, whose values are not confined to a finite set, and the minimum-value step fails.
A unit's residual has terms with up to M paid letters, so this case does matter for XW units. **Experiment E1** checks
Theorem 8 numerically and **E1b** checks Theorem 8′.

**E1** [computed; `e1/e1_gram.py`, `e1/e1_analyze.py`; outputs `e1/remote-out/ssh/out_*.txt` and the exact Gram files
`e1_gram_<menu>_f<f>.npy`; Colab #1]. Exact Gaussian-integer Grams N·G_f(i,j) = tr(w_i^* w_j) of menus of free XW words
at block 2 (letters of `../../moonshot/mw2-xwlayers/xw2lib.py`), after merging words that agree up to an f-dependent
global phase. An asymptotic-only relation would appear as an eigenvalue of G_f, outside the exact null space, that
decays towards 0 as f grows.
- Menu "main": the 57 free words of mw2-xwlayers' relation finder, 26 generators and all their 676 products in pairs
  (759 words, 514 projective families). Rank 152, 391, 400, 402, 402, 402, 402 at f = 2..8. The smallest nonzero
  eigenvalue at f = 4..8 is 0.0129, 0.0286, 0.0345, 0.0372, 0.0386: rising and levelling off near 0.04. **No decaying
  eigenvalue: no asymptotic-only relation in this menu.** 111 relations hold with fixed coefficients at f = 5..8;
  one more holds at even f only and another at odd f only (X1 flips the parity of row 0 iff f is odd, so X1.Pr0 and
  Pr0.X1 differ by a sign (−1)^f), the residue-class effect of Theorem 1(d). On the complement of the 111, the
  smallest eigenvalues at f = 5..8 are 0.0286, 0.0345, 0.0372, 0.0386 (`out_analyze_main_3to8.txt`).
- Menu "chain3": all products of up to three of 13 chain, parity, boundary and CS letters (2,380 words, 1,178
  families). Rank 890, 975, 994, 996, 994 at f = 3..7. 182 relations hold with fixed coefficients at f = 5..7; two more
  hold at odd f only (e.g. Pr0 commutes with B.Ai iff f is odd). On the complement of the 182 the smallest eigenvalues
  are 0.0115, 0.0125, 0.0130, 0.0132 at f = 5..8 (`out_analyze_chain3_3to7.txt`; f = 8 run from `e1c/`, analysis
  `e1an/out_analyze_chain3_4to8.txt`); rank 996 again at f = 8. Again rising: no asymptotic-only relation, as Theorem 8
  requires.
- **E1b, one paid letter** (`e1p/e1_paid.py`, Colab #5; analysis `e1an/out_analyze_paid_3to6.txt`): words L.D_z.R with
  L, R among 16 free words and z = 1, 2, 3 (768 words, 722 distinct up to a global phase; a built-in control checks the
  sparse word matrices against xw2lib's apply_word). Rank 256, 711, 712, 712, 712 at f = 2..6; 10 fixed-coefficient
  relations at f = 4..6; on their complement the smallest eigenvalues are 0.0115, 0.0135, 0.0146 at f = 4, 5, 6, and the
  minimising combination's norm stays at 0.0140 ± 0.0006 from f = 3 to 6. No asymptotic-only relation with one paid
  letter, as Theorem 8′ requires (E1b was run before the proof; f ≤ 6, because a word with a paid letter has 2^f
  entries per column).

## 4. The m → ∞ limit: the Cantor group, and why it restricts back exactly

**Setting.** Replace each column F_2^m by the Cantor group Z_2^N (rows 1, 2, 3, ...), with Haar measure; one array is a
vector of H = L²(Z_2^N)^(⊗f). By Pontryagin duality H has the orthonormal basis of characters χ_ξ, ξ ∈ (⊕_N Z_2)^f
(finitely supported). Paid letters D_z, z ∈ Z_2^N \ {0} (any support), act by the multiplier i^[ξ_j·z] on each column
(F1; ξ_j·z is a finite sum).

**Lemma 5 (the limit target exists)** [proved]. K_∞ := the multiplier i^|ξ| (per column) is unitary, and
K_∞ = s-lim_n ∏_(r≤n) D_(e_r): on each χ_ξ the product is constant once n ≥ max supp ξ, and uniformly bounded
operators converging on a dense set converge strongly. So **the plain schedule is a limit unit**: the strongly
convergent infinite product of the m = ∞ directions.

**Free letters at m = ∞.** Call a letter *finitary causal* if, for some n_0, it maps L²(𝓕_n) into itself for every
n ≥ n_0, where 𝓕_n is the σ-algebra of rows 1..n in every column. Examples: affine maps x ↦ Ax + b with A lower
triangular in the row index (finite or infinite); Z4 phases depending on finitely many rows; chain maps, prefix XOR
and parity phases between columns (they act row by row); boundary letters; Toffoli and CS on finitely many rows;
anything acting on rows ≤ n_0 only.

**Lemma 6 (restriction is a homomorphism)** [proved]. Let 𝒜_n be the algebra of operators on H that map
H_n := L²(𝓕_n)^(⊗f) into itself. Then y ↦ y|_(H_n) is an algebra homomorphism 𝒜_n → B(H_n). Every D_z and K_∞ lie in
every 𝒜_n, with D_z|_(H_n) = D_(z|n) (the block-n paid letter, or the identity if z|n = 0) and K_∞|_(H_n) = K_n; a
finitary causal letter restricts to a letter of the same kind at block n (n ≥ n_0).
*Proof.* (yz)|_(H_n) = y|_(H_n) z|_(H_n) when both preserve H_n; linear combinations likewise. D_z and K_∞ are
diagonal in the character basis, and H_n is spanned by the χ_ξ with supp ξ_j ⊂ [n]; on those, ξ_j·z = ξ_j·(z|n). ∎

**Theorem 7 (m-limit: nothing new)** [proved, given Theorem D of the ladder, which mw2-xwlayers' Theorem D′ extends to
every game with bounded-support letters: a unit with W arrays at block n has M ≥ W + n − 1].
- (a) No limit unit at m = ∞ has finitely many paid moves: its restriction to H_n would be a block-n unit with at most
  M moves for every large n, against M ≥ W + n − 1.
- (b) Call an infinite program *locally finite* if for each n only finitely many of its moves have z|n ≠ 0 and only
  finitely many letters act non-trivially on H_n. Its restrictions are exact finite units U_n with M_n moves, and
  conversely a limit unit is determined by them. So a limit unit is exactly a coherent family of finite units, and
  its value is sup_n ln(nW/M_n)/ln n: no more than its best finite restriction.
- (c) Calibration. The plain schedule has M_n = nW (value 0). A finite unit U ("W in M" at block b) gives the limit unit
  U^(⊗∞) acting on consecutive blocks of b rows; M_(kb) = kM, and its restrictions are k-fold Kronecker sums with value
  ln(bW/M)/ln(kb), largest at k = 1. A non-unit at block n is never the restriction of a limit unit (Lemma 6 is a
  homomorphism). So known units appear and known non-units do not.
- (d) [argued] Non-causal letters such as the two-sided shift on Z_2^Z break Lemma 6. For a translation-invariant
  layered program of finite range on Z_2^Z with ρ moves per array per row, wrapping onto the cycle Z/(Lp) (L large
  enough for the light cone) gives a block-Lp unit with ρLp moves per array, value ln(1/ρ)/ln(Lp), decreasing in L.
  So shift-invariant limits are worth at most their smallest admissible cycle. Any saving at m = ∞ must come from a
  scale hierarchy, which is the composition [m]·[k] = [mk] of meta-spectrum, where the exponent is one number
  (Theorem C) and nothing is amplified (ex-laser Theorem 3).

## 5. The other limits

**5.1 Ultraproduct over f with nonstandard coefficients (L2)** [proved]. By Łoś's theorem an internal unit of
∏_U M_(N_f)(C) with coefficients in C^U is the same as units with f-dependent coefficients at U-most f, which do not
count (ex-note §1, item 1). One genuine transfer holds, by Noetherian compactness: fix a skeleton and let S_F ⊂ C^p be
the variety of coefficient vectors right at every f ∈ F. If S_F ≠ ∅ for every finite F, then the S_F form a
downward-directed family of Zariski-closed sets, which stabilises at some finite F_0 (Hilbert basis theorem), so
∩_(all f) S_f = S_(F_0) ≠ ∅: a unit right at every f exists. This needs a solution common to F, not one per f.

**5.2 Ultraproduct over primes (L4).** A fixed skeleton's existence statement is a first-order sentence over Z[i, 1/2];
by the Lefschetz principle (ACF_0 is complete and is the limit theory of the ACF_p) it holds over C iff it holds over
F̄_p for infinitely many p. This is exact and is ex2-modp's route; nothing to add here.

**5.3 Hilbert–Schmidt / sofic limits (L6)** [argued]. Allow f-dependent free letters with uniformly bounded row and
column supports and entries, and require only ρ_R(f) → 0. Three things block transfer: (i) the unknowns of a unit are
scalars multiplying fixed families, so stability theorems for approximate representations (which move matrices, not
scalars) do not apply; (ii) a correction found by stability need not have bounded support; (iii) f-dependent letters
break "one fixed recipe". The known lower bounds survive anyway: mw2-xwlayers' Lemma O and Theorem T′ use only uniform
support and entry bounds, so they hold for f-dependent letters and HS-approximate identities alike.

**5.4 The ε-limit (L5).** A border unit is a unit over C((ε)) with residual O(ε). Its transfer theorem exists in house:
meta-border's Theorem X (formally polynomial border units count at full value, with O(log log N) overhead) and
v-meta-border's correction (relation-using ones exactify at a cost between τ_exact and (a+1)M, sometimes factor 1)
[read: `../../meta/border/EXACTIFY.md` §0, §3; `../../moonshot/BRIEF2.md`, update of ~12:20]. This is the only limit
whose transfer is non-trivial, and it is not a softening of the f- or m-direction: by Corollary 2, item 3, the
f-limit's soft points in E and Wd are exactly border points.

**5.5 By-product: the trace prices moves** [argued]. Let p = (1 + Z_r^(⊗f))/2 (Z_r^(⊗f) is a column-uniform Pauli,
free in E) and let z have z_r = 0. After an affine relabelling with z ↦ e_1 that keeps row r, D_z acts on the
2^((m−1)f) lines "bit 1 varies in every column, all other bits fixed", one recursive call of size 2^f per line. Row r
is constant on every line, so a content p·c vanishes on half of the lines, and D_z(p·c) needs only half of the
recursive calls: its real cost is τ(p) = 1/2 of a move. In general a move whose content is supported on a union of
D_z-lines of density δ costs δ, which is the f-limit trace of the support projection. The games charge 1 per move.
Whether the move-count bounds (Theorem R, D, T′) hold with these trace-weighted costs is unchecked [guess: yes,
because a certificate kept by Z_r^(⊗f) splits into its two eigenspaces, and at odd f the free letter X_r^(⊗f) swaps
them].

## 6. The dichotomy

[Proved in the instances §2, §4, §5.1, §5.2; argued as a pattern.] Every limit examined falls into one of two kinds.
- **Rigid.** A homomorphism (L3, Lemma 6), a faithful trace (L1 in E and Wd, Theorem 1; in XQ and Xchi, Theorem 9;
  in XW through one paid letter, Theorems 8 and 8′), elementary equivalence (L4)
  or Noetherian compactness (L2) carries limit units back to finite exact units of the same shape. Then the limit adds
  no existence power: its units are the finite units, and calibration passes because the limit *is* the finite
  problem.
- **Soft.** The limit forgets an exact constraint (bounded coefficients: L5; scalar drift: Remark 2.3; exactness and
  uniformity: L6; asymptotic-only relations in XW among words that combine CS with two or more paid letters, if
  any — Theorems 8, 8′ and 9 exclude them everywhere else). Then it accepts objects that do not transfer at full value,
  and the transfer theorem must restore exactly the forgotten constraint. The one such theorem we own is
  exactification for L5.

So limits cannot replace constructions. In particular the room at large blocks (m^0.9 moves per array allowed against
≥ 2.9447 proved) cannot be claimed by a limit argument: the m-limit restricts to finite blocks exactly (Theorem 7) and
the f-limit is the regular representation of the same finite problem (Theorem 1).

## 7. Calibration against the brief's list

| Item | f-limit (L1) | m-limit (L3) |
|---|---|---|
| plain schedule | limit unit (exact) | limit unit: strong infinite product (Lemma 5) |
| E8 unit | limit unit (it is a formal E unit; Theorem 1(e)) | U^(⊗∞) (Theorem 7(c)) |
| 8-move control word (o3-build) | limit unit (verified formally there; Theorem 1(e)) | U^(⊗∞) |
| block-2 stand-ins (2 in 3, 3 in 5 in V⊗V) | not applicable: a fixed representation has no f; the f-limit (all irreducibles) correctly does not accept them as engine units (excluded in E by 2W) | not applicable |
| engine 2 in 5 at block 3 | refused: L1 = the formal problem, where 14W/5 is proved | refused: restriction to block 3 is a homomorphism |
| mw2-absorb4 near-miss limit (residual 2/3) | refused: L1's residual is the formal residual, bounded below along the family | refused (restriction) |
| px-t514's 14-move program (f = 5 passes, f = 6 fails) | refused: L1 contains every irreducible (Lemma A there), including the one that first appears at f = 6 | — |
| non-exactifying border families | refused member by member (formal residual ≠ 0); accepted only by L5, whose transfer needs Theorem X | — |

Everything that applies passes. It passes because the transferable limits coincide with the exact problem; for the
same reason they are blind as existence tools.
