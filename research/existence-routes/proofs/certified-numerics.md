<!-- header: existence-routes -->
# Certified numerics, Theorem 9 and the relaxations

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** Theorem 9 is argued by hand and checked by computer (two exact checkers, re-run 19 October 2026); one AI skeptic upheld it. Lemma 10 holds under its stated condition; the summary sentence "such relaxations have depth 1" is false for one 35-dimensional piece, and "inputs-only fails with three or more pieces" overstates (the skeptic). "Every unit at a saving is unstable" is argued only.

---

# ex2-certnum: the statements behind the tool (14 Oct 2026)

Setting. F = (f_1, ..., f_N) are polynomials in x in C^n with coefficients in Q(i). The unit equations of a fixed skeleton
have this form: x holds the unknown mixing coefficients, and F(x) = P(x) - t, where P is the program map and t the exact
target. A box X is a product of complex rectangles; X(y, rho) is the box of half-width rho around y. Ball enclosures are
python-flint arb/acb (outward rounding). Nothing here is new mathematics: Lemmas 1 and 5 are standard; Lemmas 2 to 4
and Proposition 6 are short consequences, written out so that a skeptic can check exactly what the code relies on.

## Lemma 1 (Krawczyk test, complex form) [standard; Krawczyk 1969, Moore 1977, cited from memory]
Let G be square (n equations, n unknowns), X a box, y in X, C any n x n matrix. If
    K(y, X) = y - C G(y) + (I - C J_G(X)) (X - y)
is contained in the interior of X (componentwise, J_G(X) a ball enclosure of the Jacobian over X), then G has exactly
one zero z* in X, and z* lies in K(y, X).
*Remark on the complex case.* The real mean-value theorem is not used. For holomorphic G and x, y in X,
G(x) - G(y) = (integral over s in [0,1] of J_G(y + s(x - y)) ds) (x - y), and the integral lies in J_G(X) because a
product of complex rectangles is convex. With this the usual proof goes through (Brouwer for existence, the
contraction of I - C J for uniqueness).

## Lemma 2 (refusal in a box from a square subsystem) [proved here]
Let G consist of some rows of F together with slice equations S, and suppose Lemma 1 holds on X with zero z*. If for
one remaining row f_j the enclosure f_j(K(y, X)) does not contain 0, then F has no zero in X intersected with V(S).
*Proof.* A zero of F in X with S = 0 is a zero of G in X, hence equals z*, and z* lies in K(y, X); but f_j does not
vanish anywhere on K(y, X). []
With no slice (S empty) the refusal is unconditional on X. The plain box test "0 not in f_j(X) for some j" is the
special case without Krawczyk.

## Lemma 3 (exact closure by recognition) [trivial]
If a point q with Gaussian-rational coordinates satisfies F(q) = 0 exactly, F has a zero. The numerics only proposes
q (the simplest Gaussian rational inside the Krawczyk box, coordinate by coordinate); the proof is the exact
evaluation. When q lies in X and G(q) = 0, Lemma 1 also says q = z*.

## Lemma 4 (Route S: exact closure without recognising the root) [proved here]
Let z* be the unique zero of the square subsystem G in X (Lemma 1). Suppose that for a remaining row f there are
polynomials u, q_1, ..., q_n with Gaussian-rational coefficients such that
    u f = q_1 g_1 + ... + q_n g_n        identically, and   0 is not in the enclosure u(X).
Then f(z*) = 0. If this holds for every remaining row, z* is an exact zero of F, whether or not it is rational.
*Proof.* Evaluate the identity at z*: u(z*) f(z*) = 0 and u(z*) != 0. []
*Completeness at nonsingular roots.* If det J_G(z*) != 0 and f(z*) = 0, such u, q exist for some degree: the local
ring of C[x]/(G) at z* is C (multiplicity one), so f lies in (G) localised at the maximal ideal of z*; clearing
denominators gives u f in (G) with u(z*) != 0. For each degree bound the admissible (u, q) form the kernel of a matrix
with entries in Q(i); a C-basis of that kernel can be chosen in Q(i), and u(z*) is linear in the kernel vector, so if
some complex solution has u(z*) != 0, a Gaussian-rational basis vector does; shrinking X makes 0 not in u(X). So
Route S fails only by the degree searched, never in principle, at a nonsingular isolated root. []
In the code, the kernel is computed exactly (fmpz_mat.nullspace, real and imaginary parts split), and every basis
vector is tested with ball arithmetic.

## Lemma 5 (global refusal) [standard, the easy half of the Nullstellensatz]
If c_1 f_1 + ... + c_N f_N = 1 identically for some polynomials c_k, then F has no zero in C^n. Degree 0 (constant
c_k): then |F(x)| >= 1/|c| for every x (Cauchy-Schwarz), so the infimum of the residual is bounded away from 0 by an
exact number. The code finds degree-0 certificates in two ways: (i) read off a numerical near miss, y = conj F(x0),
rescaled, rounded to small denominators and checked exactly [heuristic for finding; the check is exact]; (ii) exact
rational linear algebra (rank [A | b] = rank A for the coefficient matrix A of the products m f_k with multiplier
monomials m of degree <= d). Restricting the multipliers to the alpha unknowns is allowed: it only shrinks the search.
*Why (i) works at constant-cofactor floors* [argued]: at a minimiser x0 of |F|^2 with |F(x0)| > 0, J_F(x0)^H F(x0) = 0;
if a degree-0 certificate y exists, every point attaining the bound 1/|y| has F(x0) parallel to conj(y). At borders
(infimum 0, approached only as x -> infinity) no degree-0 certificate exists, and the near miss gives no y; then (ii)
at degree >= 1 is needed.

## Proposition 6 (what numerics alone can and cannot prove) [proved here; part (c) argued]
Write the unit equations as P(x) = t. Let T be the smallest affine subspace containing the image P(C^n); it is
computable exactly (it is spanned by the coefficient vectors of P), and t in T is an exact linear check.
(a) *Open (regular) points.* If at some solution x* the Jacobian of P has rank dim T, then choosing dim T exact linear
    coordinates on T and n - dim T rational slices gives a square system with a nonsingular zero near x*, which
    Lemma 1 certifies from a good enough approximation; every other equation is an exact affine consequence of the
    chosen ones (the components of P - t normal to T are the constants of the exact check t in T). So existence is
    proved by numerics plus exact LINEAR algebra only, with no recognition of x*. Such points are stable: every target
    near t in T is still solvable.
(b) *Closed points (special position).* If the image of P has empty interior in T (for instance, the generic rank of
    J_P is below dim T), then the targets in T that have no solution are dense, and Gaussian-rational ones among them
    come arbitrarily close to t. Every step of Lemmas 1 and 2 depends on t only through enclosures, so a certificate
    built from those steps alone, valid for t, would be valid for targets within its working radius, some of which are
    unsolvable. Hence in case (b) every sound existence certificate must contain an EXACT identity involving t: a
    recognised point (Lemma 3), a Route-S identity (Lemma 4), or an equivalent piece of exact algebra. []
(c) *Where our objects sit.* By ex-dimension (meta/ex-dimension/PROOF.md, Lemmas B and C, upheld), the exact program
    map of a unit never fills a neighbourhood of its target below about 3.7, 14.5 and 66 moves per array at blocks 2, 3
    and 4; so every unit with a saving is in case (b) [argued: it identifies their ambient output space with T; if T
    were smaller than that space, "not filling the ambient space" would not exclude filling T, and the claim would need
    their count redone inside T]. The 15 distinct solvable plain-cost relaxations tested here are in case (b) by direct
    computation (opentest.py). In the V(x)V stand-in the 2-in-3 point is in case (b): the Jacobian rank at the
    solution is 520 against 544 effective equations (skeptic v-ex-dimension-logic, out_indep.txt). The GCT relaxation
    of 5 in 14 inside V_56 behaves as case (a): a generic perturbation of its target is absorbed (gct/out_tperturb.txt).

## Consequence for refusals that rest on a rank deficit [argued]
"Rank <= r" is a closed condition, "rank >= r" an open one. A certificate of the second kind is numerical (a ball
enclosure of a nonzero r x r minor, or a certified lower bound on the r-th singular value). A certificate of the first
kind needs an exact kernel. A modular rank gives only rank_Q >= rank_p (for integer matrices), so it never proves a
deficit either. px-t514's headline "fails at f = 6" is a deficit (potential 8,848 of the 8,960 needed, read from a
floating-point eigenvalue gap, px/px-t514/STATUS.txt 16:30), so it is float evidence until 112 independent exact
kernel vectors are exhibited. *Corrected in session 2 (23:24, 14 Oct):* session 1 called its modular slot count
(slot 2's f = 6 slice 560 against 448, at two primes) an excess that modular arithmetic proves. That was wrong. The slice
is a NULL SPACE (px/px-t514/f6slice.py takes F.nullspace of the stacked generator constraints), and for an integer
matrix nullity_Q <= nullity_p. So the modular count gives nullity_Q <= 560, and the failure needs nullity_Q > 448: the
same closed (rank-deficit) kind as the float count. Both are strong evidence; neither is a proof. A proof needs 560
independent kernel vectors over Q(i), checked exactly (modular kernels plus rational reconstruction).

# Session 2 (14 Oct, 23:01 on): the V(x)V stand-in, proved

Setting (ex-dimension's stand-in, as in meta/ex-v-ex-dimension-logic/indep_core.py). V = C^4 (block 2, two bits),
B = the group generated by X0, X1, the two CNOTs, S0, S1 and i (|B| = 3072, monomial matrices), rho(g) = g (x) g on
V (x) V (16-dim), A = span{rho(b) : b in B} (the free mixing), R = rho(D0) the paid move, target Kr = rho(D1 D0), with
D0 = (1+i)/2 I + (1-i)/2 X0, D1 likewise with X1. A W-in-s scheme: sources src_0..src_{W-1} = E_l (identity in block
l of a 16 x 16W matrix); move j: src_{W+j} = R (sum_i M_ji src_i), M_ji in A; output t: sum_i N_ti src_i = Kr in block t,
N_ti in A. (This is exactly the system of solve23.py / standin.py.)

## Lemma 7 (the free algebra splits into three matrix blocks) [proved; exact check in sx/construct.py]
Let U6 = span{e_a e_b + e_b e_a : a < b}, U4 = span{e_a e_a}, L2 = span{e_a e_b - e_b e_a : a < b}, so
V (x) V = U6 + U4 + L2 and Sym^2 V = U6 + U4. Then A = End(U6) + End(U4) + End(L2) (block-diagonal matrices).
*Proof.* Each b in B is monomial, b e_a = lambda_a e_{pi(a)}; so rho(b) maps e_a e_b to lambda_a lambda_b e_{pi a} e_{pi b}:
diagonal tensors to diagonal tensors, and symmetric (antisymmetric) off-diagonal ones to symmetric (antisymmetric)
off-diagonal ones. Hence A lies in the block-diagonal algebra, of dimension 36 + 16 + 36 = 88. The exact rank of the
1536 distinct matrices rho(b) is 88 (realified integer matrix, python-flint), so A is all of it. []
R and Kr commute with the swap, so they preserve Sym^2 V and L2 (checked exactly); they mix U6 and U4.

## Lemma 8 (the stand-in is a subspace game) [proved]
Write each source s as its U6-rows s_6 (6 x 16W, in the adapted basis) and U4-rows s_4 (4 x 16W), and restrict to
Sym^2 (the L2 part is treated at the end). Let P = sum_i rowspace(s_i,6) and Q = sum_i rowspace(s_i,4).
(a) The A-combinations of the sources are exactly the matrices [X; Y] whose 6 rows lie in P and 4 rows lie in Q
    (Lemma 7: the mixing acts as an arbitrary 6 x 6 block on U6-rows and an arbitrary 4 x 4 block on U4-rows).
(b) An output equals the target iff the U6-rows of the target lie in P and its U4-rows lie in Q: P >= (U6 + U6 K) (x) C^W
    and Q >= (U4 + U4 K) (x) C^W, with U6 K the row space of the U6-rows of K (the Sym^2 part of Kr), etc.
(c) A move with c = [c_6; c_4] (rows of c_6 in P, of c_4 in Q, arbitrary) adds to P the rows of R_66 c_6 + R_64 c_4,
    whose rows modulo P are those of R_64 c_4: so P grows by at most rank R_64 = 2, by exactly the span of any two
    vectors of the old Q one chooses (R_64 = a_1 b_1^T + a_2 b_2^T, and b_k^T c_4 can be any two vectors of Q).
    Likewise Q grows by at most rank R_46 = 2, by any two chosen vectors of the old P.  [ranks computed exactly]
*Proof.* (a), (b) are the row-space description of a left module over a full matrix algebra; (c) as written. []
Extra rows in P or Q never hurt (b is a containment), so nothing like "junk" exists in this game.

## Theorem 9 (exact schemes in the stand-in) [proved; exact verification for W = 2, 3, 4]
In the V(x)V stand-in at block 2, a W-in-s scheme exists if and only if s >= ceil(3W/2). In particular 2 in 3, 3 in 5
and 4 in 6 exist, with every coefficient a Gaussian rational, and 1 in 1, 3 in 4 do not.
*Proof.* Exact data: dim(U6 + U6 K) = 9 and dim(U4 + U4 K) = 7, and P starts as U6 (x) C^W, Q as U4 (x) C^W.
Only if: P must gain (9 - 6) W = 3W dimensions and gains at most 2 per move (Lemma 8c).
If: the needed gains are G_P = ((U6 + U6 K) cap U4) (x) C^W, inside the initial Q, and G_Q = ((U4 + U4 K) cap U6) (x) C^W,
inside the initial P, both of dimension 3W. Move j (on the inputs alone) transfers the next two basis vectors of G_P
into P and the next two of G_Q into Q (Lemma 8c). After ceil(3W/2) moves P and Q contain the targets' rows, so each
output is an exact linear solve (Lemma 8b). On L2 the target already lies in A (Lemma 7: A contains all of End(L2)):
take every move's L2 block 0 and output t's L2 block K|L2 on input t. []
The scheme for W = 2 is in sx/out/scheme_2in3.json in the original 16-dim coordinates. Two exact checks pass: construct.py
(sympy; outputs equal Kr in all 512 entries per output) and verify_scheme.py (own group BFS with Fractions, membership
of every mixing matrix in A by exact rank, program evaluated with Fractions). For W = 3, 4 construct.py's check passes.

## Remarks [argued]
1. Where the special position sits. Random moves fail (session 1) because the transferred vectors must lie in the fixed
   3W-dimensional spaces G_P, G_Q; that is a linear condition on the move coefficients in the adapted basis. So the
   "exact step" Proposition 6(b) demands came here from the representation, not from recognising numbers.
2. Why the stand-in shows a saving at all. The construction writes Kr = N_0 + sum_j N_j R M_j with N, M in A: the target
   has paid depth 1 in this representation, while in the faithful setting K = D1 D0 needs two paid levels. Generally,
   for a multiplicity-free free algebra S = sum End(tau) and a paid letter h, S h S contains all of Hom(tau', tau) as soon
   as the block h_{tau,tau'} is nonzero; so once h connects the types, S + S h S is everything and any target has depth
   1. Multiplicity-free relaxations are depth-blind: what is left is a capacity question (ranks of the blocks of h).

## Lemma 10 (depth 1 in multiplicity-free relaxations) [proved]
Let S = End(tau_1) + ... + End(tau_r) act block-diagonally on V = tau_1 + ... + tau_r (each piece once), and let h be any
matrix on V with blocks h_{ab}: tau_b -> tau_a. Then S h S contains all of Hom(tau_b, tau_a) whenever h_{ab} != 0.
Hence if every off-diagonal block of h is nonzero, S + S h S = End(V), and every target k has depth 1.
*Proof.* S h S contains X h_{ab} Y for all X in End(tau_a), Y in End(tau_b); their span is a nonzero
End(tau_a)-End(tau_b) sub-bimodule of Hom(tau_b, tau_a), which is simple, so it is everything. The diagonal blocks are
in S. []
[computed, gct/out_capacity.txt] For a2-gct's pieces irr 4 (7 + 14), irr 5 (6 + 21) and V_56 (7 + 21 + 28) every
off-diagonal block of h has norm at least 1.3 (float data; nonvanishing is an open condition), so k has depth 1 there.
For irr 2 (1 + 7 + 7) one block is 0 and S + ShS has dimension 211 of 225.

## Corollary 11 (two-piece multiplicity-free relaxations) [proved, given the ranks]
In the module game of a2-gct's relax2.py (outputs right iff, for every piece tau, the target's tau-rows lie in the span
of the sources' tau-rows) with two pieces tau, tau', each occurring once: let need_tau = rank k[tau, tau'] and
cap_tau = rank h[tau, tau']. Then a W-in-M relaxed scheme exists iff M >= W need_tau / cap_tau and
M >= W need_tau' / cap_tau'.
*Proof.* Exactly as Theorem 9: the tau-row space grows by at most cap_tau per move (Lemma 8c), must grow by W need_tau,
and the needed vectors lie in the tau'-coordinates of the inputs, from which cap_tau of them can be copied per move. []
With the float ranks (gaps 0.94 against 1e-14): irr 4 and irr 5 both give need 6, cap 3 in both pieces, so 1 in 2 and
2 in 4 exist and 1 in 1, 3 in 5 do not, IF the ranks of the exact data are the float ranks. For irr 5 the needs are full
row ranks (6 = number of rows), an upper bound that holds automatically, and the caps are lower bounds, an open
condition; so for irr 5 only data closeness is missing. For irr 4 the need 6 of 7 is a closed condition and needs the
exact representation.
*Three or more pieces* [argued]: the bound M >= W need_tau / cap_tau holds for every piece by the same Lemma 8c; whether
it is attained depends on copying into several pieces at once from shared sources. For V_56 (5 in 14, tight at the
28-dim piece only) the numerical solution says it is.

## Computation 12 (GCT's V_56 with exact data) [proved where stated; v56/v56.py, outputs v56/out_v56_*]
Data: c2-finite-group's exact integral model (Ind_B^G(lambda_3) on Z^2835 by signed permutations; the 56 cut out by the
class sums of 2a, 2b with membership checked exactly, as in meta/c2-finite-group/PROOF.md (c)-(d)).
(a) V_56 restricted to P_0 is 7 + 21 + 28: the A-weight spaces of the 28-orbit give the 28; those of the 7-orbit give a
    28-dimensional part whose commutant under P_0 is exactly 2-dimensional, split by its idempotents into 7 + 21. Each
    piece has commutant dimension exactly 1 (computed exactly over Q from equations selected modulo p), so the pieces are
    irreducible, and they are pairwise non-isomorphic (distinct dimensions). Hence the relaxation's free algebra is all
    block-diagonal matrices (Burnside), and Lemma 8 applies with three pieces.
(b) Exact ranks: need of the kernel k per array (7, 21, 28), of the two-row target h2 (6, 18, 18); capacity of one
    paid move (4, 10, 10). By Lemma 8c, every relaxed scheme in V_56 with W1 kernel arrays and W2 two-row arrays needs
    M >= (28 W1 + 18 W2) / 10 moves. So 5 in 14 and 1 + 4 in 10 sit exactly on the bound. This re-derives c2's hook
    profile of the 28-piece (10, 18, 28) through the same exact model.
(c) Not proved: existence of 5 in 14 or 1 + 4 in 10 in V_56. Two exact inputs-only designs fail (the delivered row spaces
    miss the kernel array's needs by 21/21/42 and 7/9/16 dimensions), and float inputs-only optimisation stalls at an
    objective of 1.337e-2 for 5 in 14 (3000 iterations, the same value from two seeds), while a2-gct's chained solution reaches 1e-26. Inputs-only schemes
    do exist numerically in the two-piece irreducibles and in the 15-dim one (1e-30). So V_56's 5 in 14, unlike the
    stand-in, appears to need chained moves [computed, float; not a proof of necessity].

## Computation 13 (GCT's 21- and 27-dimensional pieces with exact data: relaxed "savings" at block 3) [proved, given the
module lemma (Lemma 8a/b) and Corollary 11; v56/v56.py with DEG/CHI2A/LAM, outputs v56/out_21_*, v56/out_27_*]
The irreducibles 27 (chi(2a) = 15, inside Ind(lambda_4)) and 21 (chi(2a) = -11, inside Ind(lambda_7)), each of
multiplicity one there, built exactly as in Computation 12. Under P_0 they split as 6 + 21 and 7 + 14, every piece with
exact commutant dimension 1. Exact ranks: in both, need (kernel) 6 and 6, capacity 3 and 3. By Corollary 11 a relaxed
W-in-M scheme exists iff M >= 2W. Exact constructions with inputs-only moves were also checked (containment of every
needed row by exact ranks) for 1 in 2 and 2 in 4 in both pieces, and 1 in 1 is refused by the bound.
Consequence: at block 3, where every engine unit needs at least 2.8 moves per array (the hook, Theorem R), these two
relaxations admit 2 per array exactly, with rational coefficients. So no argument confined to either piece can prove
even that one array needs 3 paid moves (a2-gct Proposition 2, relaxation duality). This is depth-blindness (Lemma 10:
k lies in S + ShS there), now with exact data.

## Computation 14 (the two 35-dimensional pieces) [proved where stated; v56/out_35*]
35(-5) (inside Ind(lambda_14)): pieces 7 + 28, exact commutants 1; need (kernel) 7 and 7, capacity 6 and 6. By
Corollary 11 a relaxed W-in-M scheme exists iff M >= 7W/6; checked by an exact inputs-only construction for 6 in 7, and
1 in 1 refused by the bound. So at block 3 this relaxation admits 6 arrays in 7 moves.
35(+15) (inside Ind(lambda_4); a2-gct's undecided piece): pieces 1 + 6 + 7 + 21, exact commutants 1 (the 7-dim weight
part splits as 6 + 1); need (kernel) 1, 6, 7, 13; capacity 1, 4, 4, 6. Lemma 8c gives the exact lower bound
M >= 13W/6 for every relaxed scheme there. The 21-piece's need 13 is one below the generic 14 (special position), which
fits a2-gct's report that searches there fail even on the plain 1 in 3. An inputs-only design cannot even produce the
plain 1 in 3 there (exact run: delivered ranks 0, 0, 0, 12 of needs 1, 6, 7, 13); the plain schedule is chained.
Whether the bound 13/6 is attained is not decided here.
