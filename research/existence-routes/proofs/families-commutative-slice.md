<!-- header: existence-routes -->
# Families of units and the commutative slice (Theorems 10 to 19)

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** Argued by hand, with exact certificates checked by computer (re-run 19 October 2026); one AI skeptic rebuilt the linear programs and checked the hand proofs. Theorem 18's 155/36 is superseded by Theorem E (93/20) and the 10 % exclusion of U_m by Theorem D ([star-lemma-and-slice.md](star-lemma-and-slice.md)).

---

# ex2-families: proofs

14 October 2026, evening. Track C of exist2 (`../BRIEF.md`). No skeptic has read this. Nothing was posted.
Labels: [proved] (proof here), [computed] (program and output named), [argued], [guess].

## 0. Setting: the commutative slice

The *commutative slice* CS_m is the sub-game of the engine game used in pf-logging §5b and mw2-b5beyond: paid moves
`t_z` (one direction z in F_2^m \ 0 applied to one register), free shifts `X^a = t_a^2`, free mixing of registers by any
linear map (copies, scratch, erasing); no address maps, no phase letters. Everything commutes, so the letters generate an
abelian group A, and a word is right for every number of digit columns iff it is a formal identity in C[A].

In the X-basis `t_z = diag_x(i^{<z,x>})`. The characters of A that are trivial on the shifts are exactly
`chi_Q(t_z) = (-1)^{Q(z)}` for the binary quadratic forms `Q(z) = sum_{j<=k} q_jk z_j z_k` (Q(0) = 0), and on them the
kernel `K = t_{e_1} ... t_{e_m}` takes the value `(-1)^{tr Q}`, `tr Q = sum_j Q(e_j)`. (The vectors `(Q(z))_z` form
RM(2,m) without its constant; the full character group is the Z_4-valued quadratic forms.) A program for W arrays with
moves z_1, ..., z_M in time order must, at each chi_Q, map the W inputs to the W delivering registers by `(-1)^{tr Q} I_W`.
Write `D_Q = {i : Q(z_i) = 1}` (the moves that act as -1 at chi_Q). A form is *odd* when tr Q = 1.

## Lemma 1 (transfer form) [proved]

For any program in CS_m (any number of scratch registers, any linear mixing), there are `Ã in C^{W x M}`,
`B̃ in C^{M x W}` and a strictly lower-triangular (in time order) `Γ in C^{M x M}` such that at every chi_Q the
input-to-delivery map is

    I_W - 2 Ã_D (I + Γ_D)^{-1} B̃_D,   D = D_Q,

where `Ã_D` keeps the columns in D, `B̃_D` the rows in D, `Γ_D` the principal block. Conversely every triple
(Ã, B̃, Γ) comes from a program with W + M registers. Hence the program serves the W arrays on these characters iff

    Ã_D (I + Γ_D)^{-1} B̃_D = [Q odd] I_W   for every Q.                                   (*)

*Proof.* Let the program be `Λ_M D_M Λ_{M-1} ... D_1 Λ_0` (Λ's the mixing maps between moves, any shapes; free shifts
act as the identity at chi_Q; `D_i = I + c_i e_{r_i} e_{r_i}^T` with `c_i = chi_Q(t_{z_i}) - 1 in {0, -2}`). Expanding the
product, the term for a set `T = {t_1 < ... < t_k}` of moves is `c_T α_{t_k} γ_{t_{k-1} t_k} ... γ_{t_1 t_2} β_{t_1}`
with `α_i = E_out^T Λ_M ... Λ_i e_{r_i}`, `β_i = e_{r_i}^T Λ_{i-1} ... Λ_0 E_in`, `γ_{ij} = e_{r_j}^T Λ_{j-1} ... Λ_i e_{r_i}`;
the empty term `E_out^T Λ_M ... Λ_0 E_in` equals I_W because the trivial character must give K's value 1. Let
`Ã = [α_1 ... α_M]`, `B̃` have rows β_i, and `G_{ji} = γ_{ij}` for i < j (strictly lower triangular). The sum over
T ≠ ∅ is `Ã C (I - G C)^{-1} B̃` with `C = diag(c_i) = -2 N_D` (N_D the 0/1 diagonal of D), and the block form of
`I + 2 G N_D` with respect to (D, rest) gives `-2 Ã_D (I + Γ_D)^{-1} B̃_D` with `Γ = 2G`. Converse: use registers 1..W
(arrays) and W+1..W+M; put u_i = (α_i ; e_i), v_j = (β_j^T ; x_j) with x_j chosen so that `v_i^T u_i = 1` and
`v_j^T u_i = Γ_{ji}/2` for i < j; a mixing Λ with `Λ u_i = e_r` and `v_i^T Λ^{-1} = e_r^T` exists because
`v_i^T u_i = 1`, and `Λ^{-1} (I + c e_r e_r^T) Λ = I + c u_i v_i^T`; the ordered product of these is the program.
[The expansion was checked against explicit products of reflections: `t_refl.py`, error 9e-15.]

## Lemma 2 (two-point / deletion test) [proved; this is proofs/logging Lemma F in this language]

If a program exists, every odd Q has `|D_Q| >= W`. *Proof.* The left side of (*) factors through C^{|D_Q|}, so its rank is
at most |D_Q|; the right side has rank W. Equivalently: after deleting any W - 1 moves, `I = class(K)` stays in the
F_2-span of the classes `z z^T`.

## Theorem 1 (tight-pair clique bound) [proved]

Call an odd form Q *tight* when `|D_Q| = W`. Join two moves when some tight form contains both, and call a move *live*
when some tight form contains it. If a program exists, every set Y of live moves that are pairwise joined has `|Y| <= W`.

*Proof.* Let `K = B̃ Ã` (M x M). Its rank is at most W. For a tight Q, (*) reads `Ã_D (I + Γ_D)^{-1} B̃_D = I_W` with
square factors, so `Ã_D` and `B̃_D` are invertible and `I + Γ_D = B̃_D Ã_D = K[D, D]`. So on every tight support, K has
1 on the diagonal and 0 above it (Γ is strictly lower triangular in time order). If Y is pairwise joined and live, every
entry of K[Y, Y] above the diagonal lies in some tight support and is 0, and every diagonal entry is 1. So K[Y, Y] is
unit lower triangular, `|Y| = rank K[Y, Y] <= rank K <= W`.

*Remarks.* (i) The test needs no choice of move order and no coefficients: it reads only the move multiset.
(ii) It uses only the order-2 characters, so it holds in CS_m with free shifts and with any scratch.
(iii) Slack version (used nowhere below): for an odd Q with `|D_Q| = W + s`, `K[D, D] (I + Γ_D)^{-1}` is an idempotent of
rank W, so K[D, D] is a unit lower-triangular matrix minus a matrix of rank <= s; for an even Q,
`(K[D, D] (I + Γ_D)^{-1})^2 = 0`.

## Lemma 3 (projection) [proved]

Let `π: F_2^m -> F_2^k` keep the first k coordinates. A CS_m program with moves z_1..z_M serving W arrays becomes a CS_k
program serving W arrays with the moves π(z_i) (moves with π(z_i) = 0 are the identity and cost nothing). *Proof.* Every
(Z_4- or binary) quadratic form on F_2^k pulls back along π to one on F_2^m, with `tr(Q∘π) = tr Q + (m - k) Q(0) = tr Q`,
and `chi_{Q∘π}(t_z) = chi_Q(t_{π z})`. So the projected word is right on every character of block k.

## The uniform-direction family U_m^(r)

`U_m^(r)`: r copies of each nonzero direction of F_2^m (M = r(2^m - 1)), serving `W = r 2^(m-2)` arrays. Savings if it
existed: m = 4: 0.047; m = 5: 0.158; m = 6: 0.235; m = 7: 0.292; m = 10: 0.398; tending to 1.

**Lemma 4 (it passes every earlier test) [proved].** (a) Every nonzero binary quadratic form vanishing at 0 is a nonzero
word of RM(2,m), so its weight is at least 2^(m-2), all of it on nonzero points; hence `|D_Q| >= r 2^(m-2) = W` for every
odd Q: the deletion test passes for every m and r. (b) Its moves per array, `(2^m - 1)/2^(m-2)`, equal the one-distribution
LP bound of pf-logging at m = 4 and 5 (15/4, 31/8) [read: theory-guide RULES E1; mw2-b5beyond]. (c) In the engine game
every certificate we have allows it: hook 41/13 < 15/4 at block 4, 2.99172 < 31/8 at block 5, 2.9447 < 4 - 2^(2-m) from
block 4 on, and M >= W + m - 1. So U_m is the natural (GL(m,2)-invariant) family that no earlier test excludes at any m.

**Corollary 5 [proved, with one finite computation].** For every m >= 4 and r >= 1, U_m^(r) does not exist in the
commutative slice.

*Proof.* For U_4^(r) the minimum-weight odd forms are the indicators of the 56 affine planes of F_2^4 that avoid 0 and
contain an odd number of basis vectors; their supports hold 4r = W moves, so they are tight. `tightpairs.py` checks that
every pair of nonzero points (and every point) lies in one of these 56 planes [computed: out_tightpairs.txt, "uncovered
pairs 0"]. Then all 15r moves are live and pairwise joined (two copies of one point are joined too), and Theorem 1 gives
15r <= 4r, false. For m > 4, Lemma 3 projects U_m^(r) onto U_4^(r 2^(m-4)), which does not exist. [Direct checks: the same
covering holds at m = 5 (465 flats, 240 odd) and m = 6 (1953 flats, 992 odd); at m = 3 it fails (9 uncovered pairs), as it
must, since U_3 = 2 in 7 exists.]

## Lemma 6 (what a "family" can give: the dichotomy) [proved; elementary]

For a family with `M_m / W_m = C m^d (1 + o(1))`, the saving `s_m = ln(m W_m / M_m) / ln m` tends to `1 - d`. If W_m and
M_m are polynomials in m, or in q = 2^m (counts of orbit-natural configurations: subsets, subspaces, flags), then d is
an integer or the ratio grows like a power of q. So a natural family either has saving tending to 0 or below (d >= 1),
or saving tending to 1 (bounded moves per array). Intermediate savings need non-polynomial counts. With d = 1 and
leading constant C < 1 the saving `ln(1/C)/ln m` is largest at the smallest member, so such a family adds nothing beyond
its first block (consistent with ex-laser's "one finite unit"). Families that grow by a fixed gadget
(`c(m + k) <= c(m) + γ`) are of the d = 1 kind. *Consequence:* a family is worth building at large m only if its moves per
array stay bounded; U_m is the natural one of that kind, and Corollary 5 closes it in the commutative slice.

## Corollary 7 (everything that projects onto a tight uniform block-4 target) [proved]

Any CS_m program (m >= 4) for W arrays whose moves, projected onto some 4 coordinates (Lemma 3), give r copies of each
nonzero direction of F_2^4 with W = 4r, does not exist (Lemma 3 + Corollary 5). Example: the odd-weight family
`V_m` = all 2^(m-1) directions of odd weight, W = 2^(m-3) arrays (4 moves per array; savings 0.139, 0.226 at m = 5, 6,
tending to 1): each w in F_2^4 has exactly 2^(m-5) odd-weight preimages, so V_m projects onto U_4^(2^(m-5)) with
W = 4 * 2^(m-5). Excluded for every m >= 5. [Direct checks: `oddfamily.py`: m = 5 (W = 4, clique 16) and m = 6 (W = 8,
clique 32) EXCLUDED; m = 4 (2 in 8, the plain pair of representations) passes, as it must.]

## Calibration record (what the tests and the solver were run on)

| Case | Known truth | Deletion | Theorem 1 clique | LM solver (order-2 system) |
|---|---|---|---|---|
| plain 1 in 3, 2 in 6, 3 in 9 (block 3); 2 in 8, 4 in 16 (block 4); 3 in 15 (block 5) | exist | pass | 1, 2, 3; 2, 4; 3 (pass) | found (1 in 3, 2 in 6 interleaved) |
| 1 in 4 with directions 7, 11, 13, 14 (block 4) | exists | pass | 1 (pass) | found |
| 2 in 8 from the two 4-element representations, interleaved (block 4) | exists | pass | 2 (pass) | found, K identity checked |
| 4 in 16 = two copies of that (block 4) | exists | pass | 4 (pass) | not run to the end (too slow on the full system) |
| U_3 = 2 in 7, all directions, lex and Singer order (block 3) | exists | pass | 2 (pass) | found, K identity checked |
| 2 in 3 (block 2), 2 in 5 (block 3, two multisets), 2 in 7 (block 4), 2 in 6 non-plain and 3 in 8 (block 3) | impossible | fail | - | floors 3.00, 20.0, 20.0, 288, 18.0, 24.0 |
| U_4 = 4 in 15 | (excluded here) | pass | 15 | floors >= 540 after 120 s (3 orders) |
| U_4 plus e1 (4 in 16), U_4 plus e1, e2 (4 in 17) | (excluded here) | pass | 14, 8 | not run |
| U_4^(2), U_5, U_6, V_5, V_6 | (excluded here) | pass | 30, 31, 63, 16, 32 | not run |

The tests and the solver say "exists" on every known unit of the slice they were given and "no" on every known
impossible one; no case has the tests excluding something the solver found. Outside the slice (E8, the block-2
representation stand-ins, the o3 control word, absorb4's limit, px-t514) this route is blind: those objects use address
maps, phase letters or non-abelian letters, and the transfer form of Lemma 1 needs commuting letters.

## Corollary 8 (the slice LP bounds are never attained at blocks 4 and 5) [proved, given an LP computation]

The one-distribution LP `min sum_z p_z  s.t.  sum_z p_z Q(z) >= 1 (every odd Q), p >= 0` has value 15/4 at block 4 and
31/8 at block 5, and its optimum is UNIQUE there: p_z = 1/4 (resp. 1/8) for every z [computed, floating-point LP with
HiGHS, range of every p_z over the optimal face: `lp_unique.py` -> out_lp_unique.txt]. A slice unit with exactly
M = 15W/4 (block 4) passes the deletion test, so c/W is LP-optimal, so its multiset is U_4^(W/4); Corollary 5 excludes
it. The same at block 5 with U_5^(W/8). (At block 3 the LP value is 3 and the optimum is not unique; no strictness claim.)

## Proposition 9 (exhaustive at small W) [computed]

For a target (m, W, M), every multiset passing the deletion test has counts c_z in [W lo_z, W hi_z], where lo_z, hi_z are
the least and greatest p_z on the LP region {sum p <= M/W, constraints}. `box_enum.py` lists every integer vector in that
box with sum M, applies the exact deletion test, then the clique test. Results are in out_box_enum_b4*.txt and
out_box_enum_b5.txt. Block 4: (W, M) = (4, 15), (5, 19), (6, 23) -- the only costs between the LP bound and plain for
W = 4, 5, 6 -- have 1, 2 and 3 deletion-passing multisets, ALL excluded. So in the commutative slice at block 4 no unit for
at most 6 arrays saves anything. (The LP box uses floating-point LP bounds with a 1e-9 tolerance; an exact version would
certify each bound by a rational dual solution.)

---

# Session 2 (14 Oct 2026, 23:01 EDT onwards): the matrix E, local-rank singularity, and the near-tight uniform family

Same setting and notation as Lemma 1: `Ã` (W x M), `B̃` (M x W), `Γ` strictly lower triangular in time order,
`K = B̃ Ã` (M x M, rank <= W), and for a binary quadratic form Q its support `D = D_Q` (moves acting as -1). The
programs' identity on the order-2 characters is (*): `X_Q := Ã_D (I + Γ_D)^{-1} B̃_D = [Q odd] I_W`. The *slack* of an
odd form is `|D_Q| - W` (>= 0 by Lemma 2). (In refl.py's convention Γ is strictly upper and K = B^T A; nothing changes.)

## Theorem 10 (the matrix E) [proved]

Let `E := I_M + Γ - K`. For every binary quadratic form Q with support D:
 (b) Q even  =>  `E[D, D]` is nonsingular;
 (c) Q odd   =>  `rank E[D, D] = |D| - W` (the slack of Q);
 (a) if moves i, j (i = j allowed) lie in a common tight support, `E_ij = 0`.

*Proof.* Let L = I + Γ_D (unit lower triangular, invertible) and consider the bordered matrix
`N = [[I_W, Ã_D], [B̃_D, L]]` of order W + |D|. Its Schur complement with respect to L is `I_W - X_Q`; with respect to
I_W it is `L - B̃_D Ã_D = E[D, D]`. Hence `|D| + rank(I_W - X_Q) = rank N = W + rank E[D, D]`. Even Q: X_Q = 0, so
rank E[D,D] = |D|. Odd Q: X_Q = I_W, so rank E[D,D] = |D| - W. (a): for a tight D the factors `Ã_D`, `B̃_D` are square
and invertible, so `L = B̃_D Ã_D = K[D, D]`, i.e. `E[D, D] = 0`. []
[Numerical check on 6 solver-found units (blocks 3, 4): |E| on joined positions <= 1e-9, smallest singular value of
E_D over even forms >= 0.019, singular values beyond |D| - W on odd forms <= 1e-9: `t_theorem10.py`, out_t_theorem10.txt.]

*Corollary 10' (E-pattern tests; order-free, coefficient-free).* Let Z be the set of positions forced to 0 by (a). If for
some even Q the bipartite graph of positions of D_Q x D_Q outside Z has no perfect matching, no program exists (every
matrix with that zero pattern is singular). If for some odd Q that graph has no matching of size |D_Q| - W, no program
exists. (`even_test.py`.) Theorem 1 is the special case "a clique of joined live moves of size > W": then K is unit
triangular there; Corollary 10' instead uses the even forms.

## Lemma 11 (local rank forces singularity) [proved; exhaustive check over F_2]

Let s >= 0 and let E be an n x n matrix over a field, n >= s + 1, all of whose principal submatrices of order <= s + 2
have rank <= s. Then E is singular.

*Proof* by induction on s. Suppose E nonsingular.
(i) If E has a nonsingular principal submatrix E_A with 1 <= t = |A| <= s, let S = E/E_A (Schur complement, order
n - t >= s - t + 1, nonsingular). For B disjoint from A, `rank E[A∪B, A∪B] = t + rank S[B, B]`; for |B| <= s - t + 2
the left side is <= s, so every principal submatrix of S of order <= (s - t) + 2 has rank <= s - t. By induction S is
singular: contradiction.
(ii) Otherwise every principal submatrix of order <= s + 2 is singular; in particular E_ii = 0 and E_ij E_ji = 0 for
i != j. Let G be the digraph with an arc i -> j iff E_ij != 0: no loops, no 2-cycles. Since det E != 0, some permutation
has all E_{i σ(i)} != 0, so G has a directed cycle; let C = c_1 -> ... -> c_L -> c_1 be a shortest one (L >= 3). C has
no chord: an arc c_a -> c_b with b != a + 1 (mod L) closes the cycle c_a -> c_b -> c_{b+1} -> ... -> c_a of length
((a - b) mod L) + 1 <= L - 1. So on any set P of vertices of C the nonzero entries of E[P, P] are exactly the arcs of C
inside P, in distinct rows and columns. If L <= s + 2, E[C, C] is a monomial matrix, nonsingular: contradiction with
(ii). If L > s + 2, P = {c_1, ..., c_{s+2}} gives rank E[P, P] = s + 1 > s: contradiction. []
[Exhaustive check over F_2: no nonsingular counterexample for (n, s) = (3,0), (3,1), (4,0), (4,1), (5,1), (5,2); the
order bound is sharp: with order <= s + 1 instead there are 522 nonsingular 4 x 4 counterexamples for s = 1 and 59,794
nonsingular 5 x 5 ones for s = 2 (`lemma11/l11.c`, `l11b.c`, out_lemma11_F2.txt).]

## Theorem 12 (slack-s covering test) [proved]

If some even form Q_e with |D_{Q_e}| >= s + 1 has every set of at most s + 2 of its moves inside the support of an odd
form of slack <= s, then no commutative-slice program exists.

*Proof.* By Theorem 10(b) E[D_{Q_e}, D_{Q_e}] is nonsingular. A principal submatrix of it on <= s + 2 moves lies inside
E[D, D] for an odd D of slack <= s, so its rank is <= s by Theorem 10(c). Lemma 11 makes E[D_{Q_e}, D_{Q_e}] singular. []
(`triple_test.py` is s = 1, `cover_test.py` any s. Copies of one direction lie in the same supports, so the test runs
on direction sets.)

## Lemma 12' (covering subflats by odd flats) [proved by hand; checked by brute force for m <= 7]

Let F = {z in F_2^m : z_1 = z_2 = 1} (an (m-2)-flat avoiding 0 that holds no basis vector; the support of the even
form z_1 z_2). Every affine subflat A of F of dimension d <= m - 3 lies in an (m-2)-flat G with 0 not in G holding
exactly one basis vector; G is the support of an odd binary quadratic form with Q(0) = 0.

*Proof.* Write A = a + V, V a d-dimensional subspace of {z : z_1 = z_2 = 0}. An (m-2)-flat through a containing A is
G = a + ker ψ for a surjective linear ψ = (φ_1, φ_2): F_2^m -> F_2^2 vanishing on V; then 0 is not in G iff
ψ(a) != 0, and e_j is in G iff ψ(e_j) = ψ(a). The functionals vanishing on V form a space V^⊥ of dimension m - d >= 3
containing z_1 and z_2; pick t in V^⊥ outside span(z_1, z_2). Take φ_1 = z_1 (so φ_1(a) = 1 and only e_1 has
φ_1(e_j) = 1) and φ_2 = t or t + z_2, whichever satisfies φ_2(e_1) = φ_2(a) (exactly one does: adding z_2 keeps
φ_2(e_1) and flips φ_2(a), because a_2 = 1). Then ψ is surjective, ψ(a) = (1, φ_2(a)) != 0, and the only basis vector
with ψ(e_j) = ψ(a) is e_1. G's indicator φ_1(z) (φ_2(z) + φ_2(a) + 1) is quadratic, vanishes at 0, and has
trace = #(basis vectors in G) = 1. []
[`flatcover.py`: every subflat of F of dimension <= m - 3 is covered, for m = 4..7 (6, 28, 14, 140, 30, 620, 62 subflats
checked); F itself (d = m - 2) is not, as it must not be.]

## Corollary 13 (the near-tight uniform family, small slack) [proved]

For every m >= 4, r >= 1 and 0 <= s <= m - 4, the uniform multiset U_m^(r) (r copies of every nonzero direction) serving
W = r 2^(m-2) - s arrays does not exist in the commutative slice.

*Proof.* The odd (m-2)-flats avoiding 0 have |D| = r 2^(m-2) = W + s moves (slack s). F of Lemma 12' is the support of
an even form with r 2^(m-2) >= s + 1 moves. A set of <= s + 2 moves of F has directions spanning an affine flat of
dimension <= s + 1 <= m - 3, which lies in an odd (m-2)-flat avoiding 0 (Lemma 12'). Theorem 12 applies. []
(s = 0 re-proves Corollary 5 without the clique. Direct runs agree: cover_test excludes U_5 W = 7, U_5^(2) W = 15,
U_6 W = 15, U_6 W = 14 and passes U_5 W = 6, U_6 W = 13, U_6 W = 12, where s > m - 4: out_triple_family.txt,
out_cover_family.txt.)

## Corollary 15 (the uniform family, slack up to a tenth) [proved]

For every m >= 4 and r >= 1, U_m^(r) serving W arrays with `s = r 2^(m-2) - W < r 2^(m-2) / 10` does not exist in the
commutative slice. Equivalently: W > (9/10) r 2^(m-2) is impossible.

*Proof.* Let N = r 2^(m-2) = |D_F| (F as in Lemma 12'); E_F := E[D_F, D_F] is nonsingular (Theorem 10(b)). If R ∪ C
(sets of moves) has its directions inside an affine hyperplane of F (an (m-3)-flat), then E[R, C] is a submatrix of
E[D_G, D_G] for an odd (m-2)-flat G of slack s (Lemma 12'), so rank E[R, C] <= s. Pick a hyperplane H of F, its
complement H' = F \ H, and a second affine functional h_2 on F independent of H's; it splits H = K ∪ L and H' = K' ∪ L'.
The sets H, H', K ∪ K', K ∪ L', L ∪ K', L ∪ L' are all hyperplanes of F ({h = 0}, {h = 1}, {h_2 = 0}, {h + h_2 = 0},
{h + h_2 = 1}, {h_2 = 1}). Write E_F as the sum of its blocks on H x H, H' x H', and the eight cross blocks
K x K', K x L', L x K', L x L' and their transposes' positions K' x K, L' x K, K' x L, L' x L. Each of these ten blocks
has rank <= s, so N = rank E_F <= 10 s. []
At m = 4 this excludes every U_4^(r) target that would save (saving needs W > 15r/4 = 3.75 r > 3.6 r). For m >= 5 it
leaves (2^m - 1) r / m < W <= 0.9 r 2^(m-2) open.

## Corollary 16 (transfer by projection) [proved]

Any commutative-slice multiset at block m >= 5 whose projection onto some 5 coordinates (Lemma 3) is U_5^(r) is
impossible with W > 7.2 r arrays (Corollary 15 at block 5; projection keeps W). Such a multiset has at least 31 r moves,
so a unit of this kind needs more than 31/7.2 = 4.31 moves per array. Example: the odd-weight family V_m (all 2^(m-1)
odd-weight directions) projects onto U_5^(2^(m-6)) for m >= 6 (each w in F_2^5 has 2^(m-6) odd-weight preimages; those
of w = 0 cost nothing), so V_m is impossible with more than 0.9 * 2^(m-3) arrays. At block 5 itself every above-bar
member of U_5^(r) is excluded for every r: the bar needs W >= 31 r / 5^0.9 = 7.28 r > 7.2 r. At block 4 every saving
member of U_4^(r) is excluded (W > 15r/4 > 3.6 r).

## Corollary 14 (the block-4 sweep, completed to W = 13) [computed + proved]

The 12 multisets that survived the deletion and clique tests in the exhaustive block-4 sweep (W = 12, M = 47: 4; W = 13,
M = 51: 8; all "three copies of every direction plus extras") all fail the even test of Corollary 10' (192 even forms
each without a perfect matching; out_even_test.txt). With out_box2_b4.txt this gives: in the commutative slice at
block 4, no unit for W <= 13 arrays has fewer than 4W moves; for W = 14 the levels M = 53, 54 are done (M = 55 was still
running). Caveat as in Proposition 9: the LP box uses floating-point LP bounds with a 1e-7 rounding tolerance; an exact
version would certify each bound by a rational dual solution.

## Session-2 calibration record

| Check | Result |
|---|---|
| clique + even + odd + triple tests on 880 existing units (random unions of representations at blocks 3, 4, 5, with up to 4 idle moves added) | 0 excluded (out_t_falsify2.txt) |
| known units in the calibration lists (plain at blocks 3-5, 1 in 4 = (7,11,13,14), two-representation units, U_3 = 2 in 7, U_3^(2)) | all pass (out_even_test.txt, out_triple_calib.txt) |
| known impossible: 2 in 5 at block 3 | refused (deletion) |
| block 3, W = 2, 3, 4, every deletion-passing multiset with M <= 3W + 2 (90 multisets) | new tests exclude none; the solver finds a solution for every one tried (residual <= 1e-14) (out_cross_small_3_*.txt) |
| block 4, (W, M) = (2, 9), (3, 12), (3, 13), (4, 16): deletion passers | new tests exclude none beyond the clique test; solver finds every control tried (out_cross_box_*.txt) |
| block 4, (4, 17): 211 deletion passers | clique excludes 116, the new tests 4 more; the solver stalls on all 4 at 261-272 (3 orders each) and finds every control tried (out_cross_box_4_17.txt) |
| restricted system for U_5, W = 7 (the 241 forms the proof uses) | Colab, time-limited: residuals 187-220 with F, 189-197 without (controls): no counterexample, inconclusive (remote-out/ssh/out_t_thm12_solver.jsonl) |
| Lemma 11 | exhaustive over F_2 for n <= 5, sharp in the order bound |
| Lemma 12' | brute force m = 4..7 |

## Theorem 17 (the commutative slice at block 4 never saves) [proved; exact LP certificate, independently re-checked]

Every commutative-slice unit at block 4 has M >= 4W.

*Rank-cover inequality (block m, any even form F).* Group the moves by direction: C_x moves in direction x, and write
E[x, y] for the block of E on (moves of direction x) x (moves of direction y). For an even form F with support P (as a
set of directions), E[D_F, D_F] is nonsingular of order sum_{x in P} C_x (Theorem 10(b)) and is the sum of its blocks,
so for any partition of P x P into "pieces" R x C (R, C subsets of P),
    sum_{x in P} C_x  <=  sum over pieces of rank E[R, C].
Each piece obeys rank E[R, C] <= min(sum_R C, sum_C C) and rank E[R, C] <= |D_G| - W for every odd form G whose support
contains R ∪ C (E[R, C] is a submatrix of E[D_G, D_G], Theorem 10(c)). Dividing by W (c_x = C_x / W) gives linear
inequalities in c and in one bound variable per piece; every unit gives a feasible point (piece variable = its rank / W,
taking the larger of E[x, y] and E[y, x] for a pair of transposed cells).

*Proof of Theorem 17.* At block 4 use, for every even 2-flat F (49 of them), the partition of its 4 x 4 direction blocks
into 4 diagonal blocks and 12 cross cells. The LP "minimise sum_z c_z subject to these inequalities (plus the deletion
rows)" has optimum exactly 4 (`rc_lp.py`, `rc_lp_m.py`; the deletion rows are not even needed). An exact rational dual
certificate with 83 rows (7 even planes, 33 diagonal-slack rows, 43 cross-slack rows) is in `rc_cert_b4.json`
(`rc_lp2.py`). `verify_cert_b4.py` re-derives every row from its label with its own form tables (oddness of each
support, membership, each F an even 2-flat) and sums the rows exactly: the piece variables cancel and the result is
sum_z beta_z c_z <= -4 with every beta_z >= -1, so sum_z c_z >= 4 (out_verify_cert_b4.txt). []
*Structure.* Averaging over the 24 coordinate permutations gives a symmetric certificate (415 rows, `rc_cert_b4_sym.json`,
re-verified exactly by `sym_cert_b4.py`). Its only even forms are the 7 affine planes of the even-weight subspace
{z : z_1 + z_2 + z_3 + z_4 = 0} (two orbits, two weights), and every direction ends with coefficient exactly 1. 

*Hand proof of Theorem 17* (found from that structure; `hand_b4.py`, out_hand_b4.txt).
- Let E* be the 7 nonzero even-weight directions (a Fano plane). Its 7 affine planes avoiding 0 are the complements in
  E* of its 7 lines. 4 of them contain 1111 (call them type A) and 3 do not (type B).
- Each such plane P holds no basis vector. Its indicator is the product of two affine linear functions, so it is the
  support of an even binary quadratic form.
- For each P apply the rank-cover inequality with its 16 direction blocks:
      sum_{x in P} C_x <= sum_{x in P} rank E[x, x] + sum_{x != y in P} rank E[x, y].
- Bound rank E[x, x] by the average of |D_G| - W over all odd planes G through x. Bound rank E[x, y] and E[y, x] by the
  average over all odd planes G through x and y. (Each odd plane qualifies by Theorem 10(c).)
- Divide by W and add the inequalities with weight 18 for each type-A plane and 17 for each type-B plane.
- Counting incidences of the 56 odd planes, every c_z (z = 1..15) ends with the same coefficient, 492. Each piece
  contributes -1, so the constant is 16 (4 * 18 + 3 * 17) = 1968.
- Hence 492 sum_z c_z >= 1968, i.e. M >= 4W. []
(The weights come from equalising three numbers. Averaged this way, type-A planes give coefficients 16 on odd-weight
directions, 302/21 on weight-2 directions and 180/7 on 1111. Type-B planes give 12, 96/7 and 12/7. With weights
18 : 17 all three equal 82/3 per unit weight of type A.)

## Theorem 18 (the commutative slice at block 5 cannot reach the bar) [proved; exact LP certificate, independently re-checked]

Every commutative-slice unit at block 5 has M >= (155/36) W = 4.30556 W. The bar needs M <= 5^0.9 W = 4.25669 W, so no
slice unit at block 5 reaches it; the best conceivable slice saving at block 5 is ln(5 * 36/155)/ln 5 = 0.0929.

*Proof.* Rank-cover inequalities for every even 3-flat F (225) with the ten-piece partitions of Corollary 15 (two
complementary planes H, H' of F as diagonal pieces; four line-by-line rectangles K x K', K x L', L x K', L x L' and their
transposes, for each of the 7 plane pairs and 3 line directions), together with the 36-piece partition (singletons and
cross cells), with slack bounds from the 240 odd 3-flats. LP optimum 155/36 at the uniform point c = 1/7.2
(`rc_lp5.py`). Exact certificate: 194 rows (29 ten-piece rows, 52 plane-piece bounds, 113 line-rectangle bounds; no
deletion row), `rc_cert_b5.json`. `verify_cert_b5.py` rebuilds the forms, checks every row's meaning (odd supports
containing the piece, planes and lines are flats, each ten-piece row tiles F x F exactly, F an even 3-flat) and sums
exactly to sum_z c_z >= 155/36 (out_verify_cert_b5.txt). []

*Calibration of the inequalities.* With c fixed to an existing unit's counts / W the LP must stay feasible: 404 existing
block-4 units (random unions of representations, some with idle moves) satisfy every block-4 row with margin >= 0
(out_rc_calib.txt); 122 existing block-5 units are feasible for the block-5 LP (out_rc_calib5.txt). The excluded U_4,
U_4^(3) W = 12, U_4^(5) W = 19, U_5 W = 8 and U_5^(5) W = 37 violate them, and U_5^(5) W = 36 does not (Corollary 15's
boundary).

*Validity of the characters (no reliance on g-play's reduction).* In the f-column slice program a move t_z multiplies
the entry at (y_1, ..., y_f) by i^{n(z)}, n(z) = #{c : z.y_c = 1} (g-play's D), shifts X^a = t_a^2 and scalar mixing are
free. Every binary quadratic form Q is a point evaluation of this kind for a suitable f: two equal columns e_j give
(-1)^{z_j}, and the columns e_j, e_k, e_j + e_k plus two copies each of e_j and e_k give (-1)^{z_j z_k}; at these points
every shift acts trivially and K acts by (-1)^{tr Q} (`t_charreal.py`: 600 random forms, m = 3..5). So the order-2
conditions used in Lemma 1 and Theorems 10-18 are necessary for any slice word that is right for every f; g-play's
reduction (needed for sufficiency) is not used.

*What these do to the map.* MAP section 1.3 lists "uses address relabellings or phases" as certificate-checked, with
the slice bounds 15W/4 (block 4) and 31W/8 (block 5) from the one-distribution LP. At block 4 the old bound 15/4 already exceeded the bar's 4^0.9 = 3.48, so
the bar was closed there before; Theorem 17 adds that the slice cannot save at all. At block 5 the old 31/8 = 3.875 was
below the bar's 4.257; Theorem 18 is the first proof that no slice unit at block 5 reaches the bar.

*Limits of the method* [argued]. For the uniform family the best rank-cover partition of an even flat has ten pieces (the cross
part between complementary halves always needs eight), so on its own the method gives at most about 4/0.9 = 4.44 moves
per array at large blocks, while the bar allows m^0.9. Even weight-12 supports at block 5 have no ten-piece scheme with
6 + 6 halves (checked on 300). Beyond block 5 the method cannot certify the bar.

## Corollary 19 (every block) [proved]

Every commutative-slice unit at block m >= 4 has M >= 4W, and at block m >= 5 has M >= (155/36) W.
*Proof.* Lemma 3 projects it onto block 4 (or 5) with the same W and at most M moves; apply Theorem 17 (or 18). []
So in the slice a family can never go below 4.306 moves per array. That is a constant, so this does not stop a slice
family from saving at large blocks, where the bar allows m^0.9 moves per array.
