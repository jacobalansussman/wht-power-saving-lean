# Proofs: the Singer-invariant ansatz (Propositions S1 to S9)

*Status of this file.* The proof file of the Singer-invariant track of the Reed-Muller project (15 October 2026), as written. An independent AI skeptic checked it and upheld it with four wording points, none of which changes a result: "Corollary L3 holds only in the polynomial placement" should read "is proved only there" (Theorem C3, argued by hand, later settled it for every placement); at block 5 "the polynomial placement" depends on the polynomial (deletion bound 3 for the two trinomials, 2.5 for the four pentanomials); "every long run leaves the trivial split" describes this track's solver (in the skeptic's solver 1 of 8 runs left it clearly); and the uniform test's threshold lies somewhere in (2, 2.5]. No person has read it and no proof checker has checked it. Statements inside marked "No skeptic has read this" describe the moment of writing. Paths that begin with `checks/` refer to the private working folders of the research and are not part of this repository; private machine paths were removed. Where this file and the section README disagree, the README gives the later verdict. The labels inside ([proved], [computed], [argued], [guess]) are the working note's own: "[proved]" means that a proof is written there, not that it is proved in Lean or checked by a person; the statuses this revision uses are the ones in the section README.

---

# rm-singer: proofs

15 October 2026. Track rm-singer of rmfam (`../BRIEF.md`). No skeptic has read this. Nothing was posted.
Labels: [proved] (proof here), [computed] (program and output named), [argued], [guess].

## 0. Setting

Everything rests on two proved results of exist3/e3-rmbuild (`../../exist3/e3-rmbuild/PROOF.md`):

- **Theorem A.** The no-phase engine game with a group H of address maps has a unit "W in M" exactly when the group
  algebra `C[Sym_m x| H]` has one. Moves multiply by `[zz^T]`, the free part is `C[H]` with scalar mixing, the target
  is `[I]`, and h acts by `S -> h S h^T`.
- **Lemma C.** Equivalently, for every binary quadratic form Q, in the n channels of `Ind chi_Q = C^n` (n = |H|):
  `Atil Pi_Q (I + Gamma Pi_Q)^-1 Btil = I_W (x) O_Q`. Here Atil, Btil and Gamma are made of n x n circulants
  (elements of C[H] acting left-regularly) and Gamma is strictly lower triangular in move order. The move t acts on
  channel k as `(-1)^{Q(h_k^-1 z_t)}`. `O_Q` is the 0/1 diagonal of the channels k with `<Q, h_k^-1 h_k^-T> = 1`.

H = <sigma> is the Singer cycle (multiplication by a primitive alpha of F = F_{2^m}), n = 2^m - 1, h_k = sigma^k.
Since `dir(hz) = P_h dir(z) P_h^-1` and H is transitive on nonzero directions, every move may be taken in one
direction z_0 = 1 at no cost (e3-rmbuild, section 0). So "one move per direction" is not a constraint here: only W and
M matter.

## Proposition S1 (Singer symmetry must act inside time slots) [proved]

Take a commutative-slice program, or a multi-target one (w8-rmexact Lemma 1'), whose moves are permuted by a cyclic
group C = <g> that preserves its transfer data: `Gamma_{gs, gt} = Gamma_{s,t}`, with Atil and Btil equivariant for some
permutation of the arrays. If C acts freely on an orbit O of moves, then Gamma vanishes on O x O. If all M moves form
one free orbit (one move per direction, permuted by the Singer cycle), the program serves no array at all at block
m >= 2.

*Proof.* Suppose `Gamma_{s,t} != 0` with s, t in O. Write s = g^r t with r != 0 mod |O|. Then
`Gamma_{g^i s, g^i t} != 0` for every i. Gamma is strictly lower triangular in time, so time(g^i t) < time(g^{i+r} t)
for every i. Following t, g^r t, g^{2r} t, ... returns to t, which contradicts the strict order. So Gamma = 0 on
O x O.

If O is all the moves, Gamma = 0. Lemma 1 of ex2-families then gives, at the character Q, the output
`I - 2 sum_t [Q(z_t) = 1] Atil_t Btil_t`. Array j needs its diagonal entry to be `1 - 2 [<Q, S_j> = 1]` for every Q.
The functions `Q -> [<Q, S> = 1] = (1 - chi_S(Q))/2` for distinct nonzero S are linearly independent, together with the
constant function, because the characters chi_S are distinct. The target S_j has rank m >= 2, so it is not of the form
`z z^T`, and its character has coefficient 0 on the left side. That is a contradiction. []

*Consequence.* A "Singer-invariant program" in which the Singer cycle permutes moves in time cannot exist. The only
consistent version lets the cycle act inside each time slot. There one paid move is the n virtual moves
`(t, k)`, in directions `sigma^k z_t`, all simultaneous. That is exactly the no-phase game with H = Z_n: Theorem MT
of w8-rmexact, read through Lemma C. Its virtual arrays (l, k) are the arrays "permuted by a compatible action".
**So the Singer-invariant ansatz of the brief is the no-phase game with the Singer cycle as address group.**

## Theorem S2 (the canonical Singer system; the linear-part theorem) [proved; computed m = 2..5]

Use coordinates in a self-dual basis w_1..w_m of F over F_2 (`Tr(w_i w_j) = delta_ij`; such bases exist for every
m). Then:
(a) sigma is a symmetric matrix. The virtual targets `sigma^-k sigma^-kT = sigma^-2k` are the trace forms
    `Tr(alpha^-2k x y)`.
(b) For every binary quadratic form Q and every d in F*, `sum_i Q(d w_i) = Tr(a_Q d)`. Here
    `a_Q = sum_{y in F*} Q(y) y^-1` is the coefficient of x in the univariate polynomial of Q (its *linear part*).
(c) Hence the no-phase Singer game has a W-in-M unit iff there are circulant Atil, Btil, Gamma with, for every Q,

        Atil[:, D] (I + Gamma[D, D])^-1 Btil[D, :] = I_W (x) P_J,
        D = {(t, k): Q(alpha^-k) = 1},    J = {k : Tr(a_Q alpha^-k) = 1}.

    In words: each paid move reads Q along the Singer orbit, and the target is the linear part of Q, read along the
    same orbit.

*Proof.* (a) `Tr(alpha x y) = Tr(x alpha y)`, so multiplication by alpha is self-adjoint for the trace form, which
is the dot product in these coordinates. (b) A quadratic form has the polynomial
`Q(x) = sum_a c_{2^a} x^{2^a} + sum_{a<b} c_{2^a+2^b} x^{2^a+2^b}`, with `c_{2^a} = c_1^{2^a}` and `c_1 = a_Q` (the
discrete Fourier formula over F*; signs do not matter in characteristic 2). The basis is self-dual, so
`sum_i w_i^{2^a} w_i^{2^b} = delta_ab`: the matrix `(w_i^{2^a})` is orthogonal. In particular `sum_i w_i = 1`. So
`sum_i (d w_i)^{2^a + 2^b}` vanishes for a != b, and `sum_i (d w_i)^{2^a} = d^{2^a}`. Summing gives `Tr(a_Q d)`.
(c) In these coordinates `sigma^-2k` is the Gram matrix of the vectors `alpha^-k w_i`. So
`<Q, sigma^-k sigma^-kT> = sum_i Q(alpha^-k w_i) = Tr(a_Q alpha^-k)`. Insert this into Lemma C with z_0 = 1. []

[`struct_check.py` -> `out_struct_check.txt`: m = 2, 3, 4, 5. sigma is symmetric. (a) holds. (b) holds with 0
violations over all 8, 64, 1024, 32768 forms and every channel and every d. The vectors `(Q(alpha^-k))_k` form a
cyclic code of dimension m(m+1)/2.]
[`exact_plain.py` -> `out_exact_plain.txt`: in exact rational arithmetic the plain schedule
(`Atil_t = P^{r_t}`, `Btil_t = P^{-r_t}`, `Gamma_{s,t} = 2 P^{r_t - r_s}`, `alpha^{r_t} = w_t`) satisfies every fibre
equation at m = 3 and 4, in both placements. Dropping one move, or using a non-orthonormal basis, violates every fibre.]

## Proposition S3 (where the Fourier split happens: never where the target lives) [proved; computed]

(a) If `sigma^r` fixes Q, then `a_Q = alpha^r a_Q`. So every odd form (a_Q != 0, i.e. J nonempty) has trivial
    stabiliser in Z_n.
(b) For an odd Q the fibre `Ind chi_Q` is irreducible of dimension n, by Mackey (the inertia group is trivial). In it,
    C[Z_n] acts by its regular representation, so every one of the n Fourier components of every coefficient enters
    that fibre's equation. The Z_n characters do not split the system on any fibre that carries the target.
(c) An even form fixed by a subgroup Z_s gives a fibre that splits into s blocks of dimension n/s under the Fourier
    transform. These blocks only require the program to act as the identity there.

*Proof.* (a) `a_{Q o sigma^r} = sum_y Q(alpha^r y) y^-1 = alpha^r a_Q`. (b) is standard Mackey theory for an abelian
normal subgroup (Sym_m) with trivial inertia group. (c) The fibre's channel data are (n/s)-periodic, so their Fourier
transforms are supported on the multiples of n/s. []

[`out_struct_check.txt`, orbit table. At m = 3 and 5 (n prime) every nonzero orbit is free. At m = 4 the only
non-free orbits are even: one orbit of size 3 (|K| = 10), and three of size 5 (|K| = 6, 6, 12).]

## Proposition S4 (one target set for all odd fibres) [proved; computed]

The rotation `a -> alpha^r a` is simply transitive on F*. So each odd Z_n-orbit contains exactly one form with
a_Q = 1, namely `Q = Tr(x) + q(x)` with q purely quadratic (no linear part). There are `2^(N-m)` odd orbits,
N = m(m+1)/2: 8, 64 and 1024 at m = 3, 4, 5. Every odd fibre then has the same target set
`J_0 = {k : Tr(alpha^-k) = 1}` (|J_0| = 2^(m-1)), with move set `K = J_0 xor supp(q)`. The even fibres are the
nonzero purely quadratic q themselves, with `K = supp(q)` and target 0.
[`out_struct_check.txt`: odd orbits 8, 64, 1024.]

## Proposition S5 (the placement of the Singer group matters; the symmetric one is the least constrained) [proved; computed]

A placement of the Singer group relative to the target is a choice of coordinates. Conjugating by g turns target I
and group H into target `g^-1 g^-T` and group `g^-1 H g`.
(a) In the symmetric placement (Theorem S2) deletion `M |K| >= W |J|` gives exactly M >= 2W, at every block m >= 3.
    Every odd fibre has |J| = 2^(m-1). The smallest odd support is `2^(m-2)`, attained by `Tr(beta x) Tr(gamma x)`,
    whose linear part is `(beta gamma)^(1/2) != 0`.
(b) In the polynomial-basis placement used by e3-rmbuild and w8-rmexact, deletion gives M >= 3W at m = 3, 4 and 5
    [computed, `delcheck.py` -> `out_delcheck.txt`]. w8's Corollary L3 (Singer maps at block 3 never save) is this
    deletion bound. It is placement-specific: in the symmetric placement deletion allows 2 in 5 at block 3. That case
    is still impossible, by the engine theorem, which the brief lists as known.
(c) w8's decoupled multi-target rank-cover test (uscan.py), run on the symmetric placement (`uscan_sd.py`, a copy with
    one added option), fails at M/W = 2 and passes from 2.5, at blocks 3 and 4. In the polynomial placement it passes
    from 3 and from 67/22 [computed, `out_uscan_sd_34.txt`; w8's outputs].
(d) The symmetric placement has the largest symmetry. The Frobenius x -> x^2 is orthogonal for the trace form, fixes
    the direction 1 and the target, and normalises H (`phi sigma phi^-1 = sigma^2`). So `k -> 2k` on channels is a
    symmetry of the canonical system.

So, among the exact tests we have, the canonical Singer ansatz is constrained only by engine-wide certificates: the
hook's 2.8 at block 3, 41/13 = 3.15 at block 4, and 2.99 at block 5 (ex2-families Lemma 4(c)). U_4 needs 3.75 and U_5
needs 3.875.

## Remark S6 (with phase letters) [proved: the identity; argued: the reading]

In Lemma N's block `C[Sp(2m, 2)]`, write `x(S) = [[I, S], [0, I]]` for the move side (moves `x(zz^T)`, the kernel
`x(I)`) and `y(C) = [[I, 0], [C, I]]` for the free phase side. Then `y(I) x(I) y(I) = w = [[0, I], [I, 0]]`, checked by
multiplying out. So with phase letters free, the target is the long Weyl element w up to free factors. With only the
Singer cycle as address maps, the free part is `Z_n x| {y(C)}`. The Singer ansatz with phases is: reach w with M root
elements `x(zz^T)`, on W arrays, through that free part. Like Lemma N itself, this is necessary, not sufficient; the
Pauli-nontrivial blocks are left out. Nothing was computed with phases this session.

## Proposition S7 (convex tests cannot see the Singer restriction) [proved]

Let T_n be the multi-target slice that Theorem MT of w8-rmexact attaches to the Singer game: nM virtual moves (t, k)
in directions alpha^-k (M copies of each nonzero direction), and nW virtual arrays (l, k) with targets
`sigma^-k sigma^-kT` (W copies of each). A Singer program gives a program of T_n whose coefficients are circulant, but
T_n also contains programs with arbitrary (non-circulant) coefficients.

Claim: any test of the form "a convex set C of numbers attached to the virtual moves, arrays and forms (ranks,
weights, LP variables) contains the data of every program of T_n, and C is empty" excludes a Singer program exactly
when it excludes every program of T_n, provided C is mapped to itself by the rotation k -> k + 1.

*Proof.* The address map sigma^-1 is a change of coordinates g with `g z = alpha^-1 z` and
`g sigma^-2k g^T = sigma^-2(k+1)` (sigma is symmetric in the self-dual coordinates, Theorem S2(a)). Changing
coordinates maps slice programs to slice programs, so the rotation k -> k + 1 of virtual moves and arrays maps T_n to
itself, and by hypothesis it maps C to itself. If C is not empty, averaging a point of C over the n rotations gives a
rotation-invariant point of C, which is where the data of a circulant program lie. So C is empty iff it has no
rotation-invariant point. []

*Consequence.* w8-rmexact's rank-cover LPs and the deletion and uniform tests of section S5 decide T_n, not the
Singer game. Adding "the data are rotation-invariant" to such a test changes nothing (that is the averaging above).
To exclude a Singer shape that T_n allows, a test must use facts that hold for circulant programs but not for all of
T_n (for example: at prime n, by Chebotarev's theorem on the minors of the Fourier matrix, a circulant of rank r
keeps rank exactly r on any rows R and columns C with |R|, |C| >= r), or exact algebra, or a residual floor proved
exactly.

## Observation S8 (residual floors and the trivial split) [computed; evidence only]

For a shape W in M, let f(W, M) be the smallest residual the solver reaches (sum over all fibres of |X_Q - target|^2).
If W - 1 arrays are served exactly by the plain schedule on m(W - 1) moves, the last array is left with
M - m(W - 1) moves, so f(W, M) <= f(1, M - m(W - 1)) always (the *trivial split*). At block 4, symmetric placement,
f(1, 3) = 294.5794 (5 runs agree to 1e-5). For 2 in 7, two short random runs stopped on the plateau at the split
(294.61, 294.74, after 2,600 and 3,000 steps); starts at the split with noise 0.3 crept below it (294.30, 294.43);
and all six random starts run past about 5,000 steps fell well below it, to 239.3, 244.8, 250.0, 257.5, 261.7 and
275.2, still creeping down (sections 4.1-4.2 of NOTE.md). So the trivial split is a saddle, and two arrays do share
moves in the Singer game at block 4. The proven-impossible 2 in 6 floors at 322.2 and 326.6. The floors seen are far
from 0. [`lbf.py`; jobs c4b, c4e, w4, w4b, c5f, c4g, c4h, c4i, c5h, c5g.]

## Remark S9 (the semilinear group is a legal, larger address group) [proved]

In the self-dual coordinates of Theorem S2 the Frobenius map phi(x) = x^2 is orthogonal, because
`Tr(x^2 y^2) = Tr((xy)^2) = Tr(xy)`. It fixes the direction 1 and the target I, and `phi sigma phi^-1 = sigma^2`. So
the group generated by sigma and phi (the semilinear maps x -> a x^(2^j), of order mn) is an address group in the
same placement. More address maps are more free letters, so every Singer unit is a unit for this group, and the
canonical form of Theorem S2 holds for it with channels indexed by the group (coefficients c in C[H] act by
R(c)[a, b] = c[a b^-1]). `lbfH.py` searches this game; with H = Singer it is the same system as `lbf.py`.

