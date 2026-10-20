# Proofs: block 3 under symmetry (S1 to S10)

*Status of this file.* The proof file of the tier-2 block-3 track (15 October 2026), as written. An independent AI skeptic refuted it in part: the multiple-of-7 exclusions hold only under exact symmetry, not under gauge symmetry. It did not check S6 and S9 or the six-shape enumeration. No person has read it and no proof checker has checked it. Statements inside marked "No skeptic has read this" describe the moment of writing. Paths that begin with `checks/` refer to the private working folders of the research and are not part of this repository; private machine paths were removed. Where this file and the section README disagree, the README gives the later verdict. The labels inside ([proved], [computed], [argued], [guess]) are the working note's own: "[proved]" means that a proof is written there, not that it is proved in Lean or checked by a person; the statuses this revision uses are the ones in the section README.

---

# t2-block3: proofs behind the [proved] labels (15 Oct 2026; no skeptic has read this)

Setting (engine game, block m, here m = 3 unless said): registers hold formal contents sum_l x_l e_l with x_l in C[G],
G the letter group; free steps are left multiplication by letters of B (the stabiliser of the input frame L_Z:
Paulis, affine address maps, Clifford phase letters) and scalar mixing; a paid move multiplies one register on the
left by a letter of B dir(z) B. Array l is delivered when some register holds K e_l exactly (K = dir1 dir2 dir4).
Cells: a group element g lies in the cell g B, i.e. the Lagrangian sigma(g) L_Z; its distance is the frame distance
from L_Z (0..3); left multiplication by a free letter keeps the distance, a move changes it by at most 1, and [K] is
at distance 3. L = B cap K B K^-1 ("the Levi"); its image in Sp(6,2) is the stabiliser of the opposite pair
(L_Z, L_Y = sigma(K) L_Z), which is GL(3,2).

**Symmetric program.** H a finite group, c: H -> B a homomorphism, c'_h := K^-1 c_h^-1 K, an action of H on the arrays.
tau_h maps a content g e_l to c_h g c'_h e_{h.l} and a move letter m to c_h m c_h^-1. A program in normal form (every
move writes a new register, a combination of earlier registers and inputs, multiplied by one move letter) is
H-symmetric if tau maps its set of steps onto itself. An array is regular (orbit of size |H|) or fixed. For a fixed
array tau_h must fix its target: c_h K c'_h = K, i.e. c'_h must be a letter, i.e. c_h in L. For regular orbits c'_h
may be any homomorphism into B (targets c_h K c'_h are in B K B, which costs nothing). Orbit members of one move orbit
never read each other (the read relation is a DAG and tau-invariant).

## S1. Which groups can act, and trivial actions [proved]

(a) The image of B in Sp(6,2) is P_0 = 2^6 : GL(3,2), of order 10752 = 2^9 . 3 . 7, and the kernel B -> P_0 is the
Pauli group times scalars, a 2-group. So every element of B has order dividing 2^k . 3 . 7. **B has no element of order
5.** A group of order 5 (or containing Z_5: S_3 x Z_5, D_5, Z_15, ...) can act on the letters only through a
homomorphism whose kernel contains its Sylow-5 subgroup.

(b) **Trivial action decouples.** Suppose H_0 (order n) acts on letters trivially (c_h = c'_h = 1). By S2 with c = 1
the symmetric program is a direct sum of a program for the trivial character (the a orbit sums and the b fixed arrays,
with j + s moves) and n - 1 programs for the non-trivial characters (a arrays each, j moves each; fixed registers hold
only trivial-character content, and fixed arrays have P_chi = 0 there). So a trivially acting symmetry gives a unit only
if one of these smaller programs already is one. With the arrays in regular orbits (b = 0): 5 in 14 with Z_5 needs one
array in floor(14/5) = 2 moves; 10 in 29 needs 2 in 5; 15 in 44 (Z_5 or S_3 x Z_5) needs 3 in 8; all impossible (one
array needs 3; Theorem R gives 5.6 and 8.4). **The brief's Z_5 and S_3 x Z_5 ansatzes are empty:** either the arrays
form regular orbits and the target is excluded, or the program is a direct sum of smaller programs (for instance a Z_5
on 10 in 29 with a = 1, b = 5, j = 3, s = 14 is a 6-in-17 program plus four plain arrays).

(c) A fixed array needs c(H) inside L (above). Every subgroup of order 7 of B is conjugate in B to a subgroup of the
Levi (Sylow in P_0, then Schur-Zassenhaus in the 2-group kernel).

## S2. Free move orbits split the problem by characters [proved; calibrated numerically, out_s3_decomp_check.txt]

H abelian of order n, c a homomorphism into B, and **every move orbit free**. For a register orbit (r, h) put
y_{r,h} = c_h^-1 x_{r,h} and Y_{r,chi} = sum_h chi(h)^-1 y_{r,h}. Then:
- a move orbit (letter m on (r, 0)) acts as Y_{r,chi} := m Y_{r,chi} for every chi;
- a mix x_{r',0} := sum_{s,k} mu_{s,k} b_{s,k} x_{s,k} (and its tau-images) acts as
  Y_{r',chi} := sum_{s,k} mu_{s,k} chi(k) b_{s,k} c_k Y_{s,chi}, a free C[B]-combination;
  proof: x_{r',h} = sum mu (c_h b c_h^-1) x_{s,h+k}; multiply by chi(h)^-1 c_h^-1, sum over h, substitute h' = h + k and
  c_{h'-k}^-1 = c_k c_{h'}^-1;
- a regular array orbit enters as yhat_{chi} = sum_h chi(h)^-1 c'_h y_h (a free invertible change of inputs) and leaves as
  K yhat_{chi}; a fixed array f enters (through its orbit of input copies c_h c'_h y_f) as n P_chi y_f with
  P_chi = (1/n) sum_h chi(h)^-1 c'_h, and leaves as K P_chi y_f because Q'_chi K = K P_chi, Q'_chi := (1/n) sum chi(h)^-1 c_h^-1
  (termwise c_h^-1 K = K c'_h).
So character chi carries a standalone program with j = M/n moves for a full arrays and b projected pieces P_chi y_f.
Conversely, any family of such programs (one per chi, each with j moves) assembles into a symmetric program: solve
sum_k chi(k) alpha_{s,k} c_k = beta_{s,chi} by the inverse transform; different move letters per character differ by
free letters on both sides, absorbed into the mixes. Fixed registers are replaced by orbits of free copies.
**Corollary S2.1.** If also b = 0 (all arrays regular), the symmetric W-in-M program exists iff an (W/n)-in-(M/n) unit
exists. A free symmetric ansatz never produces a unit that a smaller unit does not already give.

## S3. Block 3, order 7: no fixed moves [proved + computed]

Let s be an element of order 7 of B; WLOG (S1c) in L, sigma(s) = A (+) A^-T on L_Z (+) L_Y with A a Singer cycle. The minimal
polynomials of A and A^-T are x^3 + x + 1 and its reciprocal x^3 + x^2 + 1, different, so L_Z and L_Y are
non-isomorphic irreducible F_2[s]-modules; the centraliser of sigma(s) in Sp(6,2) keeps both and is F_8^* = Z_7, inside P_0.
A move letter m commuting with s would have sigma(m) in P_0, i.e. m in B: impossible, as m lies at frame distance 1.
Computed (out_s2_split.txt): the Levi Singer lift fixes exactly 2 of the 135 cells (L_Z, L_Y) and none of the 945
edges of the dual polar graph; all 64 Pauli twists of exact order 7 likewise have no fixed move (out_s4_lifts.txt).
**Consequences.** Every move orbit of a Z_7 acting non-trivially is free, so M = 7j; with a trivial Z_7 action S1(b)
applies (direct sums of smaller programs). 21 in 62: 62 is not divisible by 7. **No 21 in 62 has a symmetry group
containing a non-trivially acting Z_7 (Z_7, F_21, GL(3,2)); the same holds for 5 in 14, 7 in 20, 10 in 28, 10 in 29,
12 in 35 and 15 in 44 by S4-S5 (table S7).**
F_21 and GL(3,2) contain Z_7, so everything proved for Z_7 binds them.

**S3' (general form, any block m) [proved; table computed, out_s8_hyperplanes.txt].** If a free letter s commutes with
a move letter m (up to a scalar), then sigma(s) fixes L_Z and sigma(m) L_Z, hence their intersection, a hyperplane of L_Z.
So **an element whose linear part on L_Z fixes no hyperplane (equivalently, no non-zero vector) fixes no move**, every
move orbit of the cyclic group it generates is free (when its order is prime), and M must be a multiple of that order.
Elements with no invariant hyperplane: in GL(3,2) exactly the 48 elements of order 7; in GL(4,2) all elements of
order 5 (1344) and 15 (2688), the 112 fixed-point-free elements of order 3 and the 1680 fixed-point-free elements of
order 6. Block-4 consequences (for t2-block4): a non-trivial Z_5 needs M divisible by 5, so 5 in 19 has no Z_5 symmetry
at all, and Z_5 acting on regular orbits only is equivalent to a unit with a fifth of the arrays and moves (S2.1);
3 in 11, 6 in 23 or any M not divisible by 3 with a Z_3 symmetry must use the order-3 class that fixes a 2-space.

## S4. The hook bound character by character [proved; splits computed]

Restrict the inputs of a symmetric program to the chi-isotypic subspace: regular orbit inputs y_{l,h} = chi(h) c'_h^-1 v_l
(v_l arbitrary), fixed array inputs y_f in U_chi = the chi-eigenspace of rho(c'). On these inputs the register values
of an orbit are free images of one another (c_h keeps the certificate Y), so in Theorem R (module M = V (x) Lambda^4 V,
Y = N) a free move orbit raises the potential by at most 80, and a fixed move m (c m c^-1 = m) on a fixed register,
whose value lies in the c-eigenspace for chi(1)^-1, by at most sigma_{-k}(m) = multiplicity of that eigenvalue of rho(c)
on (N + rho(m)N)/N. Demands: 224 per regular orbit, d_k = multiplicity of chi(1) in rho(c') on N/E_N per fixed array.
So for every character k:  80 j + sum_{fixed moves} sigma_{-k} >= 224 a + d_k b.
Computed splits (out_s2_split.txt, out_s4_lifts.txt; every Pauli-twisted lift):
| element | d (fixed arrays) | sigma (fixed moves) |
|---|---|---|
| Singer (order 7), every lift | 32 x 7 | no fixed move |
| coordinate 3-cycle | (76,74,74) or (74,74,76) | (28,26,26) or (26,28,26) |
| coordinate transposition | (104,120), (112,112), (120,104) | (36,44), (40,40), (44,36) |
| Pauli X^v | (112,112) | (40,40) for every dir z |
Tight shapes force counts: for the 3-cycle (lift A) a tight shape needs s = b; for the transposition (lift A) s = 2b.

## S5. Depth in the trivial character [proved]

Theorem D's proof (ladder PROOF.md section 3) uses the coset module Q[G/B], Y = the line [B], V_j = cells at distance
<= j. With inputs restricted to the trivial character: c' fixes [B], so every fixed array's input lies in Y and is
visible; a free move orbit counts as one step of Theorem RF (its members are independent and their constraints
coincide on these inputs), a fixed move as one step. Hence **j + s >= a + b + 2**. For a non-trivial character the
fixed-array pieces vanish on [B] (P_chi [B] = 0), so frame arguments see only the a regular arrays: j + s >= a + 2.
Every piece still needs three paid letters: x P_chi has its terms in the cells of x (c' in B), and K P_chi != 0 lives
at distance 3.
**Corollary S5.1 (decomposing arrays does not pay).** With Z_7 and a = 0: 7j >= 7b + 14 > 3b, never a saving. In
general a saving with Z_7 needs 4b <= 14a - 15, so at least two regular orbits.

**S5' (orbit count, any symmetry group H acting through B, any block m) [proved].** Restrict the inputs to the
H-trivial isotypic part of the coset module Q[G/B]: every c'_h fixes the coset [B], so each array orbit contributes
exactly one copy of the line Y = [B] (inputs y_{h.l} = c'_h^-1 y_l, consistent on stabilisers), and on these inputs the
values of the registers of a move orbit are free images c_h (value) of one another, so their constraints "value in V_j"
coincide and a move orbit is one step of Theorem RF (its members never read each other). Hence
**(number of move orbits) >= (number of array orbits) + m - 1.** S5 is the case H = Z_7. The more symmetric an ansatz,
the stronger this is: large groups with few move orbits are excluded outright.

## S6. Only two distinct problems for the Singer [proved + computed, out_s7_equiv.txt]

(a) There is n in L of order 3 with n c n^-1 = c^2 exactly; n' = K^-1 n K is free and n' c' n'^-1 = c'^2, so
n' P_chi' = P_chi n' with chi'(k) = chi(2k). Given a program for character chi: feed it n' y (free) for every input,
run it, apply n^-1 to every output; since K n' = n K this computes K on the chi'-pieces and on the full arrays.
(b) Complex conjugation maps letters to letters (conj dir z = X^z dir z, conj K = X^7 K, free generators to free
letters) and conj(c') = Z^7 c' Z^7, so it maps the chi-problem to the conj(chi)-problem up to a Pauli in L.
The orbits of k -> 2k and k -> -k on Z_7 \ {0} are all six non-trivial characters. **A Z_7-symmetric unit exists iff
two problems with j moves are solvable: (a full arrays + b trivial pieces) and (a full arrays + b pieces of chi_1).**

## S7. What survives Z_7 [computed from S1-S6, out_s5_table.txt]

Z_7 (and so F_21, GL(3,2)): every listed brief target is excluded; the one-move-saving survivors with W <= 40 start at
15 in 42 (a = 2, b = 1, j = 6; tight; 0.0628), 17 in 49 (2, 3, 7; 0.0364), 22 in 63 (3, 1, 9; 0.0423), 24 in 70,
25 in 70, 26 in 77, 27 in 77, ... (all need a >= 2). Z_3, Z_2 and Pauli symmetries leave many shapes for every target,
because their fixed moves are plentiful; there the per-character screens are weak.

## S8. F21 = <c, n>: the trivial character splits again; 15 in 42 is excluded [proved + computed, out_s10_f21.txt]

n in L of order 3 with n c n^-1 = c^2 (S6). In the Z_7-trivial problem the residual Z_3 = <n> acts. Array orbits of F21:
size 21 (three Z_7-orbits permuted by n: three full arrays forming a free <n>-orbit), size 7 with stabiliser a conjugate
of <n> (its orbit sum is an n-fixed full array, split d = (76,74,74) of c'_n on N/E_N), size 3 (three trivial pieces
permuted freely), size 1 (an n-fixed trivial piece, split p = (12,10,10) of c'_n on the 32-dimensional c'-invariant part
of N/E_N). Move orbits: size 21 (one move per character of <n>) or 7 (an n-fixed move; the one found, a B-conjugate of
dir2, has cap split sigma = (28,26,26)). No move orbit has size 1 or 3 (S3).
**15 in 42.** Z_7 forces two regular orbits and one fixed array (S7), so under F21 the arrays are 7 + 7 + 1 (three
orbits). Moves: 42 = 21 r + 7 q with 3 r + q = 6. Tightness of the Z_7-trivial problem in each character of <n> (the analogue of S4) needs
80 r + q sigma_{-k} = 2 d_k + p_k = (164, 158, 158): (r, q) = (2, 0) fails at k = 0, (0, 6) needs 6 sigma = (164, 158, 158),
not integral, so (r, q) = (1, 3) (with sigma = (28,26,26), which the found move has). That is 4 move orbits for 3 array
orbits, against S5': 4 < 3 + 2. **No F21-symmetric 15 in 42 exists, whatever the fixed-move letters; hence none with a
symmetry group containing F21 (e.g. a lift of GL(3,2)).** The same screens leave 17 in 49, 22 in 63, 24 in 70, ... for
F21 (out_s11_f21_table.txt, which assumes sigma = (28,26,26) for every n-fixed move: only that letter type was found).
Note [computed, out_s9_levigroup.txt]: the 168 canonical Levi lifts P_A diag(i^q_A) multiply exactly up to Z-Paulis, so
GL(3,2) itself is realised inside B only as an extension 2^3 . GL(3,2); F21 (odd order) lifts exactly.

## S9. A chain's last register has a single reader [proved]

Let an output O deliver K y for a full array y, and suppose exactly one register R of the program contains terms of y
at distance 3. Write D for the distance-3 part of R's y-content. Distance-3 y-terms of any combination come only from R
(free letters and mixes keep distances), so O = c R + ... gives c D = [K]. Then C[B] K = C[B] c D lies in C[B] D, and
C[B] K has dimension |B| (B acts freely on G by left multiplication), while C[B] D is the image of e -> e D on C[B], of
dimension at most |B|. So e -> e D is injective: Ann(D) = 0. Every y-free combination (an output for another array, or a
register fed into a later move whose y-content must vanish) has distance-3 part c'' D = 0, hence c'' = 0. **A second
reader of an array's last register needs a second register with distance-3 terms of that array** (the "second route"
of mw2-absorb4 and an-synth, in a clean form).
Application [proved for this case]: in the per-character problem "2 full arrays + 1 piece in 6 moves" (S7's 15 in 42),
if every move touches exactly one full array, each array has its own 3-move chain whose last register is the only one
with distance-3 terms of it, so the piece's output reads neither last register. The piece must then reach distance 3
through the first two moves of the chains and y-free transfers between them (which need zero divisors of C[B] at the
transfers). The case where some move touches both full arrays, and the general decision, are open.

## S10. S3 = the qubit permutations, per irreducible representation [proved; splits computed]

The coordinate permutations commute with K, so S3 acts exactly (c = c'^-1 = the permutation matrix) and every array may
be fixed. For a non-abelian H the same restriction argument as S4 works isotypic component by isotypic component: the set
of admissible inputs is H-stable, an array orbit H/S contributes Ind_S^H of its demand space and a move orbit H/T
(letter fixed by T) at most Ind_T^H of (N + mN)/N, so for every irreducible rho
   sum over move orbits <rho|_T, (N+mN)/N>_T  >=  sum over array orbits <rho|_S, N/E_N>_S        (Frobenius reciprocity).
Characters [computed, out_s2_split.txt, out_s13_fixedmoves.txt]: on N/E_N the transposition has trace -16 and the
3-cycle trace 2, so N/E_N = 30 (trivial) + 46 (sign) + 74 (2-dim); on (N + mN)/N every fixed move letter found has
traces -8 and 2, i.e. 10 + 18 + 26. Per orbit type (trivial, sign, 2-dim):
| stabiliser | array orbit demand | move orbit cap |
|---|---|---|
| S3 (size 1) | 30, 46, 74 | 10, 18, 26 |
| Z3 (size 2) | 76, 76, 148 | 28, 28, 52 |
| Z2 (size 3) | 104, 120, 224 | 36, 44, 80 |
| 1 (size 6) | 224, 224, 448 | 80, 80, 160 |
With S5' for S3 and for its subgroups (out_s12_s3.txt; calibrated: plain 6 in 18 passes, 2 in 5 fails), **exactly six
S3-symmetric shapes of 5 in 14 survive** (arrays | moves, as orbit sizes):
2+3 | 2+3+3+6;  2+3 | 1+1+6+6;  1+2+2 | 1+2+2+3+6;  1+1+3 | 1+1+3+3+3+3;  1+1+1+2 | 1+1+1+2+3+3+3;
1+1+1+2 | 1+1+1+1+1+3+6. All are exactly tight in each of the three irreducibles. Limitation: fixed-move letters were
taken from a menu (84 letters fixed by a transposition, 6 by the 3-cycle, all with the same traces); a fixed letter with
other traces would change the caps.
Structural remark [proved, from S3']: a move fixed by the 3-cycle (hence every S3-fixed move) has its cell meeting L_Z
in the only 3-cycle-invariant hyperplane, the even-weight plane {000, 011, 101, 110}: it is a "direction 111" move
(normal vector 111). A move fixed by the transposition (01) has normal 001, 110 or 111. So in every S3 shape the fixed
moves alone never span the three directions; each array's distance-3 paths must use moves from the orbits of size 3
and 6.

