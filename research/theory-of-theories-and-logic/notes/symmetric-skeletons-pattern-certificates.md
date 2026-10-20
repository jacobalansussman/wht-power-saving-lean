> **Working note, kept as written.** An AI agent wrote this during the second abstract wave of 14 October 2026 (the numerical algebraic geometry track). Private paths and machine names
> were removed and nothing else was changed. No person has checked it. Where it disagrees with this section's
> [README](../README.md), the README wins, because it follows the later skeptic verdicts.
>
> **Refuted as stated.** This note's claim that all 652 symmetric block-4 skeletons are impossible for every complex
> choice of coefficients did not survive its skeptic: the normalisation "first coefficient of every move = 1" leaves
> out "first coefficient = 0". With every coefficient free, 17 of the 652 pass the relaxation; 16 of those have no
> unit and one is open (README, section on homotopy continuation).

# a2-nag: proofs (14 October 2026)

Nobody else has read this. Nothing posted. Setting and notation: mw-sym4 (`../../moonshot/mw-sym4/STATUS.txt`, sections
A-E; its code is copied unchanged into `sym4/`). B = the Z4-monomial matrices (the widened free group, which contains
the engine's); a *frame* is a right coset B g; a *skeleton* is a program whose letters (free factors b and paid
directions) are fixed and whose content coefficients alpha (first summand of every move normalised to 1) are unknown.

## Lemma 1 (frame-sum relaxation; mw-sym4, restated)

If a skeleton is a unit for array 0 with some alpha, then there are scalars beta_s with
  sum_s beta_s Phi(u_s)(alpha) = e_(K-frame, array 0),
where Phi([g] x_l) = e_(Bg, l). *Proof* (mw-sym4 phitest.py docstring; repeated): the end condition is
[K] x_0 = sum_s f_s u_s with f_s in C[B]. The map [g] x_l -> e_(Bg, l) is C-linear and sends [b g] x_l to the same
vector as [g] x_l for b in B (b B g = B g), so it sends f_s u_s to eps(f_s) Phi(u_s), eps the augmentation. []

So every skeleton has a polynomial system R(alpha, beta) = 0 (one equation per (frame, array) row, linear in beta,
polynomial in alpha), and "R has no solution over C" excludes that skeleton for every choice of complex coefficients.

## Lemma 2 (pattern certificates)

Write the relaxation as A(alpha) beta = t, t = e_K. If there are polynomials c_1(alpha), ..., c_k(alpha) and rows
J_1, ..., J_k (different from the K row) such that, identically in alpha,
  row_K(A(alpha)) = sum_i c_i(alpha) row_(J_i)(A(alpha)),                                         (*)
then (a) the relaxation has no solution for any alpha in C^n (so the skeleton is excluded, whatever its coefficients),
and (b) for every alpha, min_beta ||A(alpha) beta - t||^2 >= 1 / (1 + sum_i |c_i(alpha)|^2).
*Proof.* Put y(alpha) = e_K - sum_i c_i(alpha) e_(J_i). By (*), y^T A(alpha) = 0, and y^T t = 1. (a) A beta = t would
give 0 = y^T A beta = y^T t = 1. (b) For every beta, 1 = |y^T (t - A beta)| <= ||y|| ||t - A beta|| (Cauchy-Schwarz),
and ||y||^2 = 1 + sum |c_i|^2. []

Consequences used in the STATUS file:
- If the c_i are constants, the floor is bounded away from 0 uniformly: 1/(k+1) when every c_i = 1 (1/2, 1/3, 1/4 for
  k = 1, 2, 3), which are the floors mw-sym4's solvers kept hitting.
- If some c_i is an unbounded polynomial, the bound tends to 0 as |alpha| -> infinity: the residual can go to 0
  only along alpha -> infinity. This is the "border" behaviour mw-sym4 saw (res^2 ~ 1e-7 to 1e-10 with |alpha| ~ 1e3
  to 6e4), and by (a) it is NOT an approximation of any finite solution of the relaxation.
- A Groebner basis {1} over Q (Singular) is equivalent to the existence of SOME polynomial certificate
  1 = sum_r h_r(alpha, beta) R_r(alpha, beta) (Nullstellensatz). In every infeasible system of this round, the
  minimal infeasible subsystem found by deletion had 2, 3 or 4 rows and the lifted certificate had exactly the
  linear-in-beta shape (*): h_K = -1 and h_(J_i) = c_i(alpha) (files sys1/*.mout, sys2/*.mout).

## Example: the border skeleton (mw-sym4 cand1b.py: block 4, Z2swap, P:inp,P:prev,L:prev,P:prev, menu 2, seed 66/1)

Rows 27 (the K frame of array 0) and 28 (a frame of array 1, two moves from the kernel frame) of the exported system
(`sys1/t27_mid3_m2_s1.sing`, computed exactly; the same skeleton as `sys/b4_2in7_border.json`):
  row 27:  (1 + a1)(1 + a2)(1 + a3) (b8 + a4 b9) - 1
  row 28:  b8 + a4 b9
(b8, b9 = the last orbit u and its image tau u; a1..a4 = content coefficients.) So
row 27 - (1 + a1)(1 + a2)(1 + a3) row 28 = -1: Lemma 2 with k = 1, c_1 = (1 + a1)(1 + a2)(1 + a3). The skeleton's
relaxation has no solution; LM's 4.5e-7 is the bound 1/(1 + |c_1|^2) being driven down by |alpha| -> infinity.
(Expand by hand: the 16 monomials of row 27 are exactly (1 + a1)(1 + a2)(1 + a3) times b8 and times a4 b9.)

Reading [argued]: every term of u and tau u in array 0's kernel frame comes from the shared move's content (the
coefficients are the 8 products of (1 + a1)(1 + a2)(1 + a3)), and each of u and tau u carries exactly one term of
array 1's material, from the same shared ancestor, into a frame that no other generator reaches. So target mass and
junk mass ride together in a fixed ratio; this is the frame-sum form of ms-fano's "patterns survive every later step"
(`../../moonshot/ms-fano/NOTE.md` section 3) and of mw-sym4's observation that "array 0's residual sits half on the
target and half on one junk row".
