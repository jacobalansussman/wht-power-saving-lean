# Proofs: 4 in 15 without phase letters (Theorems MT, O, E, Corollary L3, Propositions U, N, D)

*Status of this file.* The proof file of wave 8's exact 4-in-15 track (15 October 2026), as written. An independent AI skeptic read Theorem MT, Theorem 10', Theorem E and Proposition D line by line and found them correct; its corrections (only P15j, not P15g, is a symmetric class) are listed in the section README. No person has read it and no proof checker has checked it. Statements inside marked "No skeptic has read this" describe the moment of writing. Paths that begin with `checks/` refer to the private working folders of the research and are not part of this repository; private machine paths were removed. Where this file and the section README disagree, the README gives the later verdict. The labels inside ([proved], [computed], [argued], [guess]) are the working note's own: "[proved]" means that a proof is written there, not that it is proved in Lean or checked by a person; the statuses this revision uses are the ones in the section README.

---

# w8-rmexact: proofs

15 October 2026, wave 8 (`../BRIEF8.md`, section w8-rmexact). No skeptic has read this. Nothing was posted.
Labels: [proved] (proof here), [computed] (program and output named), [argued], [guess].

## 0. Setting

From e3-rmbuild (`../../exist3/e3-rmbuild/PROOF.md`), Theorem A: the no-phase engine game with a group H of address
maps has a W-in-M unit with move directions z_1..z_M iff the group algebra `C[Sym_m x| H]` has one: moves are left
multiplications by `[z z^T]`, the free part is `C[H]` (plus scalar mixing of registers), the target is `[I]`. Here
`Sym_m` is the additive group of symmetric m x m matrices over F_2 and `h` acts by `S -> h S h^T`. For H = 1 this is the
commutative slice in its order-2 normal form (e3-rmbuild Corollary R1). The skeptic of e3-rmbuild was still running
when this was written. Every exclusion below uses only the direction "game unit => group-algebra unit" of
Theorem A (its Lemmas S and R), and nothing else from that folder.

A **multi-target slice program** at block m: registers hold elements of `C[Sym_m]^N` (row vectors indexed by N
arrays), paid moves multiply one register by `[z z^T]`, every scalar mixing of registers is free. It *serves* array j
with target `S_j` in Sym_m when some delivering register ends as `[S_j] e_j`. With every `S_j = I` this is a
commutative-slice unit in the sense of ex2-families (`../../exist2/ex2-families/PROOF.md`).

## Theorem MT (channels: an address-map unit is a multi-target slice program) [proved]

Let H be a finite group of address maps, n = |H|. If the no-phase game with H has a W-in-M unit with moves
z_1..z_M, then there is a multi-target slice program at the same block with nW arrays and nM paid moves, in which
array (l, h) has target `h h^T` (l = 1..W, h in H) and the n moves coming from move t have directions `h z_t` (h in H).

*Proof.* `C[Sym x| H]` is a free right `C[Sym]`-module with basis `[k]`, k in H, because `[S][k] = [k][k^-1 S k^-T]`.
Left multiplication is right-`C[Sym]`-linear, so it gives an injective algebra homomorphism
`Phi: C[Sym x| H] -> M_n(C[Sym])`, with `Phi([S]) = diag_k([k^-1 S k^-T])` and `Phi([g])` the permutation matrix of
`k -> g k`. Apply Phi to every register of the unit (Theorem A form). A register holding a row vector over `C[Sym x| H]`
becomes n "virtual registers" (the rows of its image, indexed by k), each holding a row vector over `C[Sym]` indexed
by (l', k'). A free letter `[g]` becomes a permutation of virtual registers, scalar mixing stays scalar mixing, and a
paid move `[z z^T]` becomes n paid slice moves, one on each virtual register k, in direction `k^-1 z`. The input
register `e_l` becomes the virtual registers `e_(l,k)`, and the delivered `[I] e_l` becomes `[k^-1 k^-T] e_(l,k)`.
Re-index k -> k^-1. []

*Remarks.* (i) This is Lemma C / Lemma V of e3-rmbuild read at the level of programs: the virtual system there is
exactly the transfer form of this program. (ii) The converse fails in general: the image program has coefficients
that commute with the right regular action of H on channels, a structure the multi-target program forgets.

**Lemma 1', Theorem 10' (multi-target transfer form)** [proved; same proofs as ex2-families Lemma 1 and Theorem 10].
A multi-target program with N arrays and moves y_1..y_P has, at every binary quadratic form Q,
`Atil_D (I + Gamma_D)^-1 Btil_D = O_Q`, where D = D_Q = {moves with Q(y) = 1} and O_Q is the 0/1 diagonal of the arrays
with `<Q, S_j> = 1` ("odd at Q"). With `E = I + Gamma - Btil Atil`:
`rank E[D_Q, D_Q] = |D_Q| - n_odd(Q)` for every Q, where n_odd(Q) is the number of arrays odd at Q.
*Proof.* Lemma 1's expansion never uses that the targets are equal; each array's value at Q is `(-1)^{<Q, S_j>}`.
The bordered matrix `[[I_N, Atil_D],[Btil_D, I + Gamma_D]]` has Schur complements `I - O_Q` (rank N - n_odd(Q)) and
`E[D, D]`. []
**(a) Tight supports.** If slack(G) = 0, then `E[D_G, D_G] = 0`. *Proof.* Let L = I + Gamma_D (invertible) and O the
odd projector, of rank |D|. From `Atil_D L^-1 Btil_D = O`, the factor `L^-1 Btil_D` has full row rank |D|, so the rows
of Atil_D at even arrays vanish; likewise the columns of Btil_D at even arrays vanish. On the odd arrays the identity
reads `Atil_odd L^-1 Btil_odd = I_|D|` with square factors, so `L = Btil_odd Atil_odd = (Btil Atil)[D, D]`, i.e.
`E[D, D] = 0`. (This is e3-rmbuild's Lemma V(c) and ex2-families Theorem 10(a).) []
Hence every rank-cover inequality of ex2-families Theorem 17 holds with `slack(F) := |D_F| - n_odd(F)`: for any
partition of supp F x supp F into rectangles R x C, `slack(F) <= sum rank E[R, C]`, with `rank E[R, C] <= #moves in R`,
`<= #moves in C`, and `<= slack(G)` for every form G whose support contains R u C; and `slack(G) >= 0` (deletion).

## Theorem O (orthogonal address maps never help without phases) [proved]

Let H consist of address maps that preserve the dot product (`h h^T = I` for every h in H; for example coordinate
permutations, or `I + J` with J the all-ones matrix, m even). Then every W-in-M unit of the no-phase game with H gives
a commutative-slice unit nW in nM, n = |H|. Consequently, exactly as in the slice:
- block 3: M >= 3W [the LP of Proposition L with H = 1 certifies 3 exactly; `out_cert.txt`];
- block m >= 4: M >= 4W (ex2-families Theorem 17 and Corollary 19);
- block m >= 5: M >= (155/36) W = 4.306 W, above the bar's 5^0.9 W = 4.257 W at block 5 (Theorem 18, Corollary 19).

*Proof.* In Theorem MT every virtual target is `h h^T = I`. The image program is then a commutative-slice word in the
order-2 normal form serving nW arrays with nM moves; Theorems 17, 18 and Corollary 19 of ex2-families are proved from
the order-2 system (*) alone, which this program satisfies. []
[Computed agreement: the cell LP of Proposition L gives exactly 4 for H = all 24 coordinate permutations and for the
whole orthogonal group O_4(F_2) (48 elements), certified exactly (`out_cert.txt`). The skeptic of e3-rmbuild found the
cyclic coordinate-permutation groups of orders 2, 3, 4 excluded by cliques, as this predicts.]

So a saving unit of the no-phase game must use address maps outside the orthogonal group. With such maps the virtual
targets `h h^T` differ from channel to channel, and that is the only thing the relaxation sees.

## Proposition L (the multi-target rank-cover LP for cyclic and full address groups) [computed; exact certificates]

LP: variables `c_x` = (virtual moves in direction x) / (virtual arrays), one piece variable per cell {x, y} (the larger
of rank E[x, y] and rank E[y, x], over virtual arrays); rows: one cover row per form F with the cell partition, every
cell bounded by c and by every slack(G) with G containing it, and deletion for every G. Minimise `sum_x c_x`, a lower
bound on M/W for every unit of the H-game (any directions). `mtlp.py`; exact dual check with independent code
`certify.py` (rebuilds every row from its label, own Singer and Gram code; rationalised duals; all negative reduced
costs zero).

| block | H | virtual targets | LP bound on M/W | exact certificate |
|---|---|---|---|---|
| 3 | trivial | 1 | 3 | 3 (17 rows) |
| 3 | Singer Z_7 | 7 | **3** | 3 (7 rows, below) |
| 3 | GL(3,2) | 28 | 2.1786 (float) | not needed (no exclusion) |
| 4 | trivial, all coordinate permutations, O_4 | 1 | 4 | 4 (40 rows, O_4) |
| 4 | Singer power 5 (order 3) | 3 | 79/23 = 3.4348 | 79/23 (78 rows) |
| 4 | Singer power 3 (order 5) | 5 | 102/31 = 3.2903 | 102/31 (77 rows) |
| 4 | Singer cycle (order 15) | 15 | 67/22 = 3.0455 | 67/22 (96 rows) |
| 4 | GL(4,2) | 420 | 2.0346 (float) | not needed (no exclusion) |

Richer partitions (all set partitions of the 4-point supports, `mtlp2.py`, 339,246 rows) do not move the order-3 value
(3.434783) or the trivial one (4). With the order-3 group and uniform directions the optimum is exactly 15/4 (below).

**Every cyclic address group of order 3, 5 or 15 at block 4** [computed, float LP; `sweep.py`, out_sweep_ab.txt].
The cell LP depends on the group's position relative to the dot product, not only on its isomorphism type.

| order | subgroups | LP range | LP >= 3.75 (4 in 15 excluded for every move pattern) | LP = 4 (no saving at all) |
|---|---|---|---|---|
| 3 | 616 | 2.667 to 4 | 76: the 4 orthogonal ones (LP 4) and 72 others (LP 3.820 or 3.857) | 4 (exactly the orthogonal ones) |
| 5 | 336 | 3.027 to 3.556 | 0 | 0 |
| 15 | 336 | 2.087 to 3.087 | 0 | 0 |

The 72 non-orthogonal order-3 groups with LP >= 3.75 were each certified exactly by certify.py (`run_cert72.sh`,
out_cert72.txt). The values are 531/139 = 3.8201 (24 groups) and 27/7 = 3.8571 (48 groups); in one case a
rationalisation gave a longer fraction with the same decimal value. So for these 76 groups, 4 in 15 is excluded
for every move pattern [proved, given Theorem A]. The other values in the table are floats.

**Corollary L3 (Singer address maps at block 3 never save)** [proved]. In the no-phase game at block 3 with the Singer
cycle Z_7 as address group (multiplication by a root of x^3 + x + 1 in the polynomial basis, the group e3-rmbuild
used), every unit has M >= 3W. [Position matters: of the 8 cyclic subgroups of order 7 of GL(3,2), 6 give the LP
value 3 and 2 give only 16/7 = 2.2857 (float LP, `sweep3.py`, out_sweep3.txt); the skeptic of e3-rmbuild saw the
same 6 + 2 split for 2 in 5 by deletion.] *Proof (the 7-row certificate, by hand).* Put
`C_x` = number of virtual moves in direction x, N = 7W virtual arrays. The seven forms
`x1 + x1x3, x1x2 + x1x3, x2 + x1x2, x1 + x1x2 + x1x3 + x2x3, x2 + x1x2 + x2x3, x3 + x2x3, x3 + x1x3 + x2x3`
each have a 2-point support, and each is odd at 6 of the 7 virtual targets (a = 6/7). Deletion (Lemma 1') gives
`C(G) >= 6W` for each. Their supports cover every direction exactly twice (they form a 7-cycle on the 7 directions),
so `2 sum_x C_x >= 7 * 6W`, i.e. `7M = sum_x C_x >= 21 W`. [The two "cell" rows of the computed certificate are
deletion rows in disguise: a cell rank is >= 0.] So M >= 3W: no unit with these maps beats the plain cost. In particular
e3-rmbuild's lost block-3 run (5 in 14 with Singer maps), if it used this group as its other block-3 runs did, is
decided: no such unit exists. []

**Block 3, every address group** [computed, float LP; `sweep3.py`, out_sweep3.txt]. GL(3,2) has 179 subgroups. The
cell LP equals 3 (no saving possible, if the float value is exact) for: all 21 of order 2, 16 of 28 of order 3, 19 of
35 of order 4, 4 of 28 of order 6, 6 of 8 of order 7, 9 of 21 of order 8, 1 of 14 of order 12, 1 of 14 of order 24
(78 subgroups in all, the 6 orthogonal ones included, Theorem O). The lowest values: GL(3,2) itself 2.1786, the Frobenius groups of order 21
2.1905, groups of order 12 and 24 with 12 targets 2.2083. A block-3 saving without phases would need one of the
groups with LP below 3 (the engine's hook bound 2.8 applies on top).

## Proposition U (4 in 15 with cyclic address groups of order 3, 5, 15: what survives) [computed, exact]

*Lemma (fixed directions decouple the LP)* [proved]. When the counts c_x are fixed, the piece variables appear only in
the cover rows (with negative coefficients) and in upper bounds that are constants. So the LP is feasible iff every
deletion row holds and every cover row holds with every piece at its upper bound
`u(R, C) = min(c(R), c(C), min_{G contains R u C} slack(G))`. This is an exact finite check per direction pattern.

For a group acting freely on directions, a W-in-M unit has `c_x = n_o / (W n)`, where n_o counts the physical moves
in the orbit o of x. `fixc.py` checks every composition (n_o) of M = 15 over the orbits, with W = 4, using deletion,
cells of every support, and all set partitions of the 4-point supports, in exact fractions.

| address group | compositions | survive (4 in 15 not excluded) |
|---|---|---|
| Singer power 5 (order 3), 5 orbits (smallest directions 1, 2, 3, 4, 5) | 3,876 | 12: (3,3,3,2,4), (3,3,3,3,3), (3,3,4,1,4), (3,3,4,2,3), (3,3,5,1,3), (3,4,3,1,4), (3,4,3,2,3), (3,4,4,1,3), (4,3,3,1,4), (4,3,3,2,3), (4,3,4,1,3), (4,4,3,1,3) |
| Singer power 3 (order 5), 3 orbits | 136 | 17, among them the U_4 pattern (5,5,5) |
| Singer cycle (order 15), 1 orbit | 1 | 1 (no tight form at all) |

*The U_4 pattern with the order-3 group sits exactly on the bound.* With uniform directions c_x = 1/4, the 21 forms
with a 4-point support that are odd at all three virtual targets have slack exactly 0 (deletion with equality), and the
LP optimum under uniform directions is exactly 15/4 (`mtlp2.py --uniform`). These are e3-rmbuild's 21 tight forms;
its clique of 12 = Wn is the same boundary seen through Theorem 1.

## Theorem E (4 in 15 with the order-3 Singer group and evenly spread moves does not exist) [proved]

Address group H = <sigma>, sigma = multiplication by alpha^5 (alpha a root of x^4 + x + 1, polynomial basis; order 3;
e3-rmbuild's group). Suppose a W = 4, M = 15 unit of the no-phase game uses every orbit of directions exactly three
times (e3-rmbuild's pattern 1,2,4,8,3 repeated three times, or all 15 directions once). Contradiction.

*Proof.* Directions are written as integers x1 + 2 x2 + 4 x3 + 8 x4.
- *Image program.* By Theorem MT there are 12 virtual arrays, 4 in each channel. Every direction carries exactly 3
  virtual moves. The three channel targets h h^T are
  `I`, `[[0,1,1,1],[1,0,1,1],[1,1,1,1],[1,1,1,0]]` and `[[1,0,1,0],[0,1,1,0],[1,1,0,1],[0,0,1,1]]`.
- *Tight forms.* A form whose support has 4 directions and which is odd at all three targets has slack
  `3*4 - 4*3 = 0`. There are 21 such forms (e3-rmbuild's 21 tight forms). On each, E vanishes on `D_G x D_G`
  (Theorem 10'(a): a tight support makes the transfer blocks square and invertible, so `K[D,D] = I + Gamma_D`).
  Six of them:
  - `x1 + x1x4` (support 1, 3, 5, 7);
  - `x1x2 + x1x3 + x1x4 + x2x4 + x3x4 + x4` (3, 5, 8, 14);
  - `x2 + x2x3 + x2x4 + x3x4` (2, 3, 12, 13);
  - `x1x2 + x1x4 + x2x3 + x2x4 + x3x4 + x4` (3, 6, 8, 13);
  - `x1 + x1x2` (1, 5, 9, 13);
  - `x1x2 + x1x3 + x2 + x2x3` (2, 5, 10, 13).

  So E vanishes on every block E[x, y] with x, y in {3, 5, 13}, the diagonal blocks included.
- *The form that fails.* `F = x1x2 + x1x3` has support {3, 5, 11, 13}, so |D_F| = 12. It is odd at exactly one
  channel target: `S_12 + S_13` = 0, 0, 1 on the three targets. Theorem 10' then gives `rank E[D_F, D_F] = 12 - 4 = 8`.
- *Rank count.* In `E[D_F, D_F]`, the 9 rows of directions 3, 5, 13 vanish outside the 3 columns of direction 11, so
  they have rank at most 3. The 3 rows of direction 11 add at most 3. So the rank is at most 6 < 8. []

[Checked by two programs with no shared code: `epat.py` (uses forms.py) and `epat_check.py` (own GF(16), forms and
targets; finds 81 violated forms, the first being this F).]

*Calibration of the E-pattern test.*
- It passes the slice's plain 4 in 16.
- With the order-3 group it passes the plain schedules in the bases e1..e4 and {7, 11, 13, 14} (the latter has
  `sum z z^T = I`). Both give the orbit counts (4,4,0,4,4), so this is a single check.
- It refuses the slice's U_4 (known impossible): the linear form x1 needs rank 4 and has none.
- The order-5 group with the pattern (5,5,5) passes, and the order-15 group has no tight form, so the test is empty
  there.

*All direction patterns* (`epat.py --comps all`, out_epat_z3.txt, out_epat_z5.txt). With the order-3 group, deletion
plus this test leaves 7 of the 3,876 orbit splits of 15 moves. All 7 also pass the rank-cover check of fixc.py. With
the order-5 group, 16 of 136 survive both, (5,5,5) among them. With the order-15 group, the one pattern survives.

**Proposition RC (rectangle partitions with forced zeros)** [proved; `rcover.py`]. For a form F, split the blocks
E[x, y] (x, y in supp F) that are not forced to zero into rectangles R x C. Each nonzero block lies in exactly one
rectangle; forced-zero blocks may lie in several rectangles or in none. Then `E[D_F, D_F]` is the sum of the
rectangles' submatrices, so
`slack(F) <= sum_rect min(C(R), C(C), min_{G contains R u C} slack(G))`.
`rcover.py` finds the cheapest such partition exactly for every 4-point support, by dynamic programming over the
covered cells.

*Overlapping covers are not allowed.* The off-diagonal part of J - I (n x n, rank n) is covered by 2 log2 n
rectangles of rank 1. A first version of rcover.py allowed overlaps; it was wrong and was replaced before any result
was used.

*Calibration.* The test passes the slice's plain 4 in 16 and the order-3 plain 4 in 16. It refuses the slice's U_4:
the 4-point form with support {3, 7, 11, 15} has every block forced to zero but needs rank 4.

*Result.*
- Order-3 group: 6 of the 7 surviving splits are excluded. **Only (3,3,4,2,3) survives every test**, and it has no
  tight form at all.
- Order-5 group: 3 of the 16 are excluded; 13 survive, (5,5,5) among them.
- So for the order-3 group, 4 in 15 is decided for 3,875 of the 3,876 direction splits.

## Theorem R over GL-submodules of the permutation blocks (brief step 1) [computed]

For each GL-orbit O of binary quadratic forms, V_O = C[O] is a G-module (move = diag((-1)^{Q(z)}), target =
diag((-1)^{tr Q})). `trsub.py` decomposes V_O into irreducible GL-submodules (eigenspaces of a generic Hermitian
element of the orbital algebra), enumerates every sum of them, and samples mixed submodules in the isotypic
components with multiplicity (14 three times, 56 three times and 20 twice in the 420-block, 14 twice in the 105-block,
20 twice in the 280-block). f(g) = rank((1 - P) Pi_g P).

| block | best f(I)/f(zz^T) | needed to exclude |
|---|---|---|
| 3 (blocks of 7, 21, 7, 28 forms) | 1.4 (21-block) | 2.5 (to refuse 2 in 5) |
| 4 (blocks of 15, 105, 35, 420, 280, 168 forms) | 2.1667 (35-block, support 12) | 3.75 (to refuse 4 in 15) |

So this family fails its block-3 calibration and is blind for 4 in 15. Since every move is GL-conjugate, an LP over
several such Y_0 gives nothing beyond the best single ratio.

## Lemma N's block C[Sp(2m, 2)] (brief step 4): calibration fails [computed]

`sp.py` (block 3): permutation modules of Sp(6,2) on the 63 nonzero vectors, the 135 Lagrangians and the 315 isotropic
planes; constituents of dimension 1, 15, 27, 35, 84, 168 (multiplicity free). Deletion in a representation reads
`M r_t >= W r_K` (r = dimension of the -1 eigenspace). Ratios r_K / r_t: 2 (dim 27), 1.6 (35), 1.25 (168), 1 (84),
0.8 (15). Theorem R with Y_0 = unions of orbits of the free parabolic (the stabiliser of {(0|b)}): at most 1.75.
Refusing engine 2 in 5 needs more than 2.5. A general limit [proved]: K is a product of three commuting transvections
conjugate to t, so in the joint eigenspace decomposition `r_K = 3a + c` and `r_t = a + 2b + c` (a, b, c: dimensions
with exactly one, exactly two, all three signs -1), hence `r_K <= 3 r_t`, and a ratio above 2.5 needs `a > 10b + 3c`.
Not run on 4 in 15 (it failed calibration).

**Proposition N (deletion in Lemma N's block is capped by the no-phase GL game)** [proved]. In every representation
of Sp(2m, 2), the deletion ratio r_K / r_t is at most `max_O (2^m - 1) a_O / s_O` over the GL(m,2)-orbits O of nonzero
binary forms (s_O = support size, a_O = odd fraction of the orbit at the target I). This is 2 at block 3 and 2 at
block 4. So deletion in that block can refuse neither engine 2 in 5 (needs 2.5) nor 4 in 15 (needs 3.75), whatever
irreducibles one reaches; the 24 constituents not computed above cannot change this.
*Proof.* The moves t_z and the target K lie in the abelian unipotent radical U_X = {[[I, S],[0, I]]} = Sym_m of the
X-parabolic. Restricted to U_X a representation is a sum of characters chi_Q (Q a binary form), with multiplicities
m(Q) constant on GL(m,2)-orbits, because GL(m,2) inside Sp(2m,2) normalises U_X. t_z = -1 exactly on the chi_Q with
Q(z) = 1 and K = -1 exactly on those with tr Q = 1. So `r_t = sum_O m_O |O| s_O / (2^m - 1)` and
`r_K = sum_O m_O |O| a_O`, a weighted average of the orbit ratios. Orbit values (`forms.py`): block 3: 1, 2, 0.67, 0.75;
block 4: 1, 2, 0.57, 0.93, 1.21, 0.79 (the maximum 2 from the 2-point (block 3) and 4-point (block 4) supports). []

## Additions after the resume (15 October, from 09:13)

**Proposition D (rank depth: one array needs m moves, for every address group)** [proved]. In `C[Sym_m x| H]`, for any
group H of address maps, every unit serving one array uses at least m paid moves.
*Proof.* Let `F_k` be the span of the elements `[h][S]` (h in H, S in Sym_m) with S a sum of at most k matrices
`u u^T`. Free letters keep `F_k`: `[g][h][S] = [gh][S]`, and scalar mixing is linear. A paid move maps `F_k` into
`F_(k+1)`: `[z z^T][h][S] = [h][S + u u^T]` with `u = h^-1 z`. The input `[1]` lies in `F_0`, so after k paid moves
every register lies in `F_k`. The target `[I]` = `[1][I]` is not in `F_(m-1)`: the elements `[h][S]` are a basis of
the group algebra, and I is not a sum of fewer than m rank-one matrices (it has rank m). []
*Use.* 1 in 3 at block 4 is impossible with any address maps. It is therefore a negative control for every
necessary test.

**Calibration failure of the multi-target relaxation at order 15** [computed, exact; `fixc_p15.sh`, out_fixc_P15.txt].
Deletion plus rank cover (all set partitions of the 4-point supports) passes 1 in 3 for 6 of the 10 placement classes
of the order-15 Singer groups (P15a, b, c, g, i, j) and refuses it only for P15d, e, f, h (e3-rmbuild's group is in
P15e). For every fixed-point-free class of order 3 and 5 it refuses 1 in 3 (out_fixc_classes.txt). So in most
order-15 placements the relaxation cannot see Proposition D, and "passes" there carries no information.

**Monotonicity** [proved, immediate]. If `H' <= H`, a unit with maps in H' is a unit with maps in H. So an exclusion for H
excludes every subgroup of H, and an existence for H' gives existence for every overgroup. Exclusions for small cyclic
groups do not propagate upward.

**Placement classes** [computed; `placements.py`, out_placements.txt]. Conjugating by g in GL(4,2) maps the H-game with
target I to the `gHg^-1`-game with target `g g^T`. So the game depends on the pair (H, I), that is on H up to
conjugation by the orthogonal group O_4 (48 elements, the stabiliser of I). Classes: order 3: 616 subgroups in 28
classes (5 fixed-point-free classes, of sizes 24, 12, 8, 6, 6); order 5: 336 in 10; order 15: 336 in 10 (named
P15a..P15j in `clm/placements15.txt`, in the order of out_placements.txt). e3-rmbuild's groups lie in the classes
P15e (order 15), [3,9,14,12] (order 5) and the class of size 8 (order 3).

**Symmetric placements** [proved, immediate; `symclass.py`]. If H is generated by a symmetric matrix g (g = g^T), then
every element of H is symmetric and every virtual target is `h h^T = h^2`, an element of H. The Singer cycle written
in a self-dual basis of F_16 over F_2 is such a group (multiplication by alpha is self-adjoint for the trace form).
Classes with a symmetric generator: order 3: the fixed-point-free class of size 6 [4,8,5,10] (and two classes with
fixed vectors); order 5: [4,12,7,2]; order 15: P15g and P15j.

**4 in 15, 2 in 7 and 3 in 11 for every fixed-point-free class of order 3 and 5** [computed, exact; `chain_classes.py`:
fixc.py, then epat.py on its survivors, then rcover.py on theirs; every split of the moves over the direction orbits].
Rows give the splits left after all three tests; "excluded" means none is left, so the shape is impossible with those
address maps, whatever the directions and their order.

| order | class (size, generator rows) | 4 in 15 (3876 or 136 splits) | 2 in 7 (330 or 36 splits) |
|---|---|---|---|
| 3 | 24, [6,8,15,10] | excluded | excluded |
| 3 | 12, [2,3,13,6] | excluded | excluded |
| 3 | 8, [6,13,10,3] (e3-rmbuild's) | 1 left: (3,2,3,4,3) | excluded |
| 3 | 6, [12,7,14,3] | excluded | excluded |
| 3 | 6, [4,8,5,10] (symmetric) | 115 left, all-once (3,3,3,3,3) among them | 2 left |
| 5 | 48, [3,8,13,15] | 3 left | excluded |
| 5 | 48, [3,9,14,12] (e3-rmbuild's) | 13 left, (5,5,5) among them | excluded |
| 5 | 48, [2,10,14,15] | 15 left | excluded |
| 5 | 48, [3,10,7,5] | 1 left | excluded |
| 5 | 48, [4,10,11,6] | 13 left | excluded |
| 5 | 24, [5,15,14,3] | excluded | excluded |
| 5 | 24, [3,6,10,13] | 5 left | excluded |
| 5 | 24, [5,10,11,13] | 16 left | 1 left |
| 5 | 12, [2,6,10,3] | 16 left | 2 left |
| 5 | 12, [4,12,7,2] (symmetric) | 20 left | 1 left |

(The e3-rmbuild rows agree with the earlier runs on its own generators: 1 of 3,876 and 13 of 136.)

*3 in 11 (saving 0.063), same pipeline* (out_chain_311.txt): excluded for every split in the order-3 classes of
sizes 24, 12, 8 and 6 ([12,7,14,3]); the symmetric class keeps 17 of 1,365. Order 5: excluded in [3,8,13,15],
[3,10,7,5] and [5,15,14,3]; the others keep 2 to 8 of 78.

*Soundness on known units* (`soundness_classes.py`, out_soundness_classes.txt): in every fixed-point-free class of
order 3 and 5 the three tests pass 1 in 4, 2 in 8 and 4 in 16 with directions 1, 2, 4, 8, and 4 in 16 with
directions 7, 11, 13, 14 (60 of 60).

