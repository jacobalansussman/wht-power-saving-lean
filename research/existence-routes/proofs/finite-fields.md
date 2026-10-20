<!-- header: existence-routes -->
# From solutions modulo a prime to complex units (Propositions 1 to 6)

**About this file.** This is a working file as an AI agent wrote it during the rounds of 14 and 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent, and "the user" is me. Its own labels ([proved], [computed], upheld, unchecked) are the agent's; the status below is the one this revision uses, and where a later verdict differs from the file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean unless the status says so. Working names of folders are mapped in [README.md](README.md).

**Status now.** Argued by hand, not checked; one AI skeptic checked Lemma 0 and Propositions 1 to 6 line by line and the Krick-Pardo-Sombra constant on the paper ([skeptic-finite-fields.md](skeptic-finite-fields.md)). The certificates are checked by computer (240 of 240 exact identities, re-run 19 October 2026). Side claim refuted: "random search never finds the a4 = 3 control" (an artefact of the search box).

---

# ex2-modp: proofs (14 October 2026, evening)

Nobody else has read this. Nothing posted. Labels: [proved] (the proof is below), [computed] (a program in this folder
and its output), [read], [argued].

## Setting

A skeleton (wiring and free steps fixed, mixing coefficients unknown) gives, for each array l, a system
f_{l,1} = ... = f_{l,N_l} = 0 in the shared content coefficients alpha and that array's coefficients beta_l
(nag-rerun's closure "cores": linear in beta_l, polynomial in alpha, integer coefficients). A unit exists for the
skeleton iff the **joint** system X (all arrays, alpha shared) has a complex point. Write X for any finite system
f_1, ..., f_N in Z[x_1, ..., x_n], usually N > n (overdetermined), and X(R) for its points over a ring R.

## Lemma 0 (characteristic-0 transfer) [proved; standard]

If X(K) is non-empty for some field K of characteristic 0 (a number field, Q_p, C_p), then X(C) is non-empty.

*Proof.* If X(C) were empty, the weak Nullstellensatz over Q-bar and descent to Q give 1 = sum g_i f_i with
g_i in Q[x]. Substituting a K-point gives 1 = 0 in K. []

## Proposition 1 (no test that reads one prime can decide existence) [proved]

For every prime p:
(a) X_1 = {x - 1 = 0, x^2 - 1 - p = 0} has X_1(C) empty, yet X_1(F_p) = {1}. At that point the Jacobian (1, 2x)^T has
    full rank 1 = n mod p, so the point is reduced, isolated and non-degenerate mod p, of the expected dimension 0.
(b) X_2 = {p x - 1 = 0} has X_2(Q) = {1/p} but X_2(F_p) empty.

*Proof.* (a) x = 1 forces -p = 0. (b) Clear. []

So nothing computed from the equations modulo one prime (points, tangent spaces, Groebner bases, dimension) can
certify a complex unit, and nothing can exclude one either. For a SQUARE system Hensel's lemma does lift a
non-degenerate point; X_1 shows exactly where overdetermination breaks it: the square subsystem {x - 1} lifts, the
extra equation is not zero at the lift. Every rigorous criterion below therefore adds either characteristic-0
information (Theorem 1, Proposition 2) or many primes (Proposition 5).

## Theorem 1 (exact point; the criterion actually used) [proved; trivial]

If x is a point with coordinates in a number field K and every f_i(x) = 0 holds exactly in K, then X(C) is non-empty.

The mod-p machinery is the SEARCH that produces x, never part of the proof:
1. find a point of X over F_p (alpha by search, guided by the alpha-only elements of a Groebner basis mod p; beta by
   linear algebra mod p);
2. lift it p-adically (Proposition 4) with some coordinates held at small integers;
3. reconstruct the lifted coordinates as rationals (rational reconstruction) or as elements of a quadratic field
   (LLL on (1, v, v^2) and (1, theta, v));
4. recompute beta EXACTLY over K and evaluate every equation of the original system exactly over K.
An error in steps 1-3 (a wrong prime, a wrong hint, a wrong reconstruction) can only lose a certificate; it cannot
create a false one, because step 4 alone decides. Programs: `modp_exist.py` (joint arrays, rational),
`special_exist.py` (hints, rational), `padic_exist.py` (hints, p-adic lift, quadratic fields).

## Proposition 2 (a one-prime criterion with a characteristic-0 certificate: smoothness over Z_p) [proved]

Let x-bar be in X(F_p). Suppose
(a) some r of the equations, f_S, have Jacobian J_S(x-bar) of rank r over F_p, and
(b) for every other equation h there are u_h, a_{h,s} in Z[x] with u_h(x-bar) != 0 in F_p and
    u_h h = sum_{s in S} a_{h,s} f_s identically.
Then X has a Z_p-point reducing to x-bar (in fact an (n - r)-parameter family of them), so X(C) is non-empty.

*Proof.* Choose r coordinates on which an r x r minor of J_S(x-bar) is invertible; hold the other n - r coordinates
at integer lifts of x-bar. The square system f_S in the r chosen coordinates has invertible Jacobian mod p, so
Hensel's lemma gives y in Z_p^n with f_S(y) = 0 and y = x-bar mod p. Then u_h(y) = u_h(x-bar) != 0 mod p is a unit of
Z_p, so h(y) = u_h(y)^(-1) sum a_{h,s} f_s(y) = 0 for every other h. Lemma 0 finishes. []

Remarks. (i) (b) is a polynomial identity over Z, checked exactly and cheaply; it is the characteristic-0 input that
Proposition 1 shows to be unavoidable. (a)+(b) amount to X being smooth over Z_p at x-bar (Jacobian criterion;
conversely smoothness gives such identities locally), which is the "right lifting criterion" for an overdetermined system. (ii) For systems linear in beta,
(b) is a statement about rows: each surplus row is a combination of r basic rows with polynomial coefficients, its
own coefficient u_h not vanishing at x-bar. The exclusions of a2-nag are the mirror image: the target row minus a
combination of junk rows is a nonzero constant. (iii) Unlike Theorem 1, Proposition 2 never needs the point's field.
It was not needed in this round: every system met here had an exact point of degree <= 2 (Theorem 1 was enough).

## Proposition 3 (the obstruction test) [proved]

Let x-bar be in X(F_p), with some coordinates F held at fixed integers, and suppose the Jacobian of the whole system
with respect to the remaining coordinates, J-bar = J(x-bar) mod p, is injective. Lift digit by digit:
x_{k+1} = x_k + p^k delta, where J-bar delta = -F(x_k)/p^k (mod p). Then x-bar has at most one lift to Z_p with the
coordinates F fixed, and if this linear system is inconsistent at some step k, x-bar has NO such lift.

*Proof.* Let y be a lift and suppose x_k = y mod p^k (true for k = 1). For k >= 1, F(x_k + p^k delta) = F(x_k) +
p^k J(x_k) delta mod p^(2k), and J(x_k) = J-bar mod p, so mod p^(k+1) the digit delta of y at position k solves the
step-k system. By injectivity it is the only solution, so the iteration follows y forever; an inconsistent step
contradicts the existence of y. []

Without injectivity the lift along one path proves nothing (another choice of kernel components may lift). This is
the rigorous form of "a spurious mod-p point is refused" (`zsys.hensel`, `spurious_scan.py`).

## Proposition 4 (lifting a forced coordinate) [proved; standard]

If a coordinate is forced by a univariate hint g(t) (coefficients in Z_p from earlier lifted coordinates) and its
residue t-bar is a simple root mod p (g'(t-bar) != 0 mod p), Newton's iteration t <- t - g(t)/g'(t) converges
p-adically to the unique root of g congruent to t-bar. (`padic_exist.lift_forced`.) Used only inside the search of
Theorem 1.

## Proposition 5 (solutions modulo many primes) [proved, given the cited theorems]

(i) If X(C) is non-empty, X has a point over some number field K, and X(F_p) is non-empty for every prime p that
splits completely in K and divides none of the point's denominators: a set of primes of Dirichlet density at least
1/[L : Q], L the Galois closure of K (Chebotarev [read, standard]). So X(F_p-bar) is non-empty for all but finitely
many p.
(ii) If X(C) is empty, Krick-Pardo-Sombra Theorem 1 [read: arXiv math/9911094, p. 2] gives a in Z \ {0} and g_i in
Z[x] with a = sum g_i f_i and ln|a| <= H := 4n(n+1)d^n (h + ln s + (n+7) ln(n+1) d), where d is the largest degree,
h the largest log|coefficient| and s the number of equations. A point of X over F_p-bar forces p | a.
(iii) Hence X(C) is non-empty iff X(F_p-bar) is non-empty for infinitely many p, iff for more than H/ln P distinct
primes >= P. The last form is a finite, rigorous test that never computes in characteristic 0 and never says where
the complex point is.

*Proof of (iii).* (i) gives "=>". For "<=": at most ln|a|/ln P <= H/ln P distinct primes >= P divide a. []

Cost. KPS's d^n is unavoidable for any bound of this type (the Masser-Philippon example on p. 3 of the same paper has
h(a) >= d^(n-1) h). The number of primes needed is 2,636 at n = 3, d = 3, about 5 * 10^5 at n = 5, d = 4, and
10^13 at n = 15 (`out/kps_table.txt`). The same paper's refinements do not rescue it [read, pp. 4-5]: Theorem 2 replaces d^n by the
"degree" delta of the system (useful only if delta is small and itself proved), and the sparse Corollary 3 bounds h(a)
by 2(n+1)^3 d V (h + ln s + 2^(2n+4) d ln(d+1)) with V the normalised volume, which still carries 2^(2n+4).
So (iii) is a proof tool for toy sizes only; at our sizes the many-primes idea
is useful as a FIELD DETECTOR: the residue classes of the primes that carry points predict the field of the point
before any reconstruction (this round: points exactly at p = 1 mod 4 -> Q(i); exactly at p = 1 mod 3 -> Q(sqrt(-3));
both then confirmed by exact points). The census is heuristic (finitely many tries per prime).

## Proposition 6 (the refutation side, for the tracks that exclude) [proved]

(a) X(F_p-bar) empty for one p (a Groebner basis {1} mod p) does not prove X(C) empty: Proposition 1(b).
(b) It does if the equations are homogeneous and X is read in projective space: Y = V(f) in P^n_Z is proper over
Spec Z, so its image is closed; if the generic fibre Y_Q were non-empty the image would be all of Spec Z and Y(F_p-bar)
would be non-empty. Our affine systems with "= 1" rows are not of this form, and their homogenizations have points at
infinity (alpha -> infinity is exactly the "border" behaviour), so (b) rarely applies.
(c) A rigorous exclusion of an affine system is an exact certificate 1 = sum g_i f_i. a2-nag's pattern identities
(row_K - sum c_i(alpha) row_J = -1, c_i integer polynomials) are such certificates over Z, so they exclude the
skeleton over EVERY field. This matches what `spurious_scan.py` saw: none of the 15 excluded systems scanned has a
point over F_p for any p tested (2 to 43). nag-rerun's "GB = 1 modulo two primes" verdicts are therefore very strong
evidence but not proofs; an exact certificate makes them proofs.
(d) The certificate can itself be found modulo p: compute cofactors g_i mod p (Singular `lift` in characteristic p),
reconstruct their coefficients as rationals, and check sum g_i f_i = 1 exactly over Q. As with Theorem 1, the modular
step is search and the exact identity is the proof; a wrong reconstruction can only fail the check
(`lift_cert.py makep` / `check`). This is Theorem 1's principle on the refutation side.
(e) For systems linear in beta (rows f_r = A_r(alpha) . beta + c_r(alpha)), a vector y(alpha) of polynomials with
sum_r y_r A_r = 0 and sum_r y_r c_r = 1 is a certificate (cofactors g_r = y_r): sum_r y_r f_r = 0 + 1. Finding y is
a membership problem for the vector (0, ..., 0, 1) in a MODULE over Q[alpha], with only the 2-4 alpha variables,
instead of an ideal problem in alpha and all the betas. This is exactly the shape of a2-nag's pattern identities. It
is not complete in principle (a certificate might need beta), but on the five cores where the full ideal lift ran for
more than 10 minutes, the module version took at most 1 s each, and every result passed the exact check
(`lift_cert.py makem`). It found no certificate for 4 of the 22 g27_Fany4_s17/s19 charts, even with the chart
relaxation's rows added; (f) handled those.
(f) nag-rerun's "extras" and `--linsub` substitutions are consequences of the chart relaxation computed modulo p (an
alpha-only polynomial e, or x - L(alpha)). Each becomes exact once the system (relaxation) + {t e - 1} is certified
empty: then e vanishes at every point of the relaxation (Rabinowitsch), hence at the alpha of every unit (a2-nag
Lemma 1), and e may be used as an exact row or substitution (`extras_cert.py`, `subst_cert.py`). With the 24 extras
certified this way, the full lift of the 4 hard charts took about 1 s each.

## Proposition 7 (soundness requirements for any "unit exists" claim from these systems) [argued]

1. **Joint arrays.** One alpha must serve every array. A per-array feasibility test (nag-rerun's decide.py labels a
   chart "UNIT CANDIDATE" when each array's system is feasible mod both primes, separately) is necessary, not
   sufficient. The certificate must be for the joint system (`zsys.joint`, used by `modp_exist.py`).
2. **Every f.** The system must encode the identity for every number of digit columns (the closure / formal
   identity), not at one f. px-t514's 14-move program passes every test at f = 5 and fails at f = 6: a system built at
   f = 5 would be certified, wrongly, by any method. The modular method has no view of f; it is only as good as the
   system it is given.
3. **The prime 2** (rests on MAP Proposition G, unchecked). A unit whose coefficients are 2-integral reduces to a program over F_2, where no saving exists
   (MAP Proposition G, "divides by 2"; unchecked there). So the complex points of a saving shape's system cannot be
   2-integral, p = 2 is bad by necessity, and Hensel at 2 is useless; use odd primes. Observed in a relaxation:
   t27_mid3_m2_s1's relaxation has exact rational points (with 2 in the denominators) and no point over F_2.
