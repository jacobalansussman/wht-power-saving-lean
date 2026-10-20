# Proofs: phase letters, Theorem T, the angle barrier and phase letters at the ends

*Status of this file.* The proof file of the phase-letter track of the Reed-Muller project's second round (15 October 2026), as written. An independent AI skeptic upheld it, with one scope note: "at the ends" must mean letters on the inputs, a commutative program that holds all the mixing, and letters on the outputs. The skeptic did not check the tiling and ILP exclusions, the completeness of the 4-in-15 enumeration, the numerical runs, or the inputs taken from other rounds. No person has read it and no proof checker has checked it. Statements inside marked "No skeptic has read this" describe the moment of writing. Paths that begin with `checks/` refer to the private working folders of the research and are not part of this repository; private machine paths were removed. Where this file and the section README disagree, the README gives the later verdict. The labels inside ([proved], [computed], [argued], [guess]) are the working note's own: "[proved]" means that a proof is written there, not that it is proved in Lean or checked by a person; the statuses this revision uses are the ones in the section README.

---

# rm2-phase: proofs

15 October 2026, track rm2-phase of rmfam round 2 (`../BRIEF2.md`). No skeptic has read this. Nothing was posted.
Labels: [proved] (proof here), [computed] (program and output named, exact unless said), [derived] (a short argument
from results of other tracks, which keeps their status), [argued], [guess].

## 0. Lemma N's block as a game

Block m. Vectors of F_2^(2m) are Pauli labels u = (x|b), with the form omega(u, v) = x.b' + b.x'. The engine group G
maps onto Sp(2m, 2) with kernel the Paulis and scalars (e3-rmbuild Lemma N; pf-induction Fact B). The images are:

| engine letter | image in Sp(2m, 2) | name here |
|---|---|---|
| paid move dir z | (x\|b) -> (x + (b.z) z \| b) | t_z = x(z z^T) |
| kernel K = dir e_1 ... dir e_m | (x\|b) -> (x + b \| b) | x(I) |
| address map A | (x\|b) -> (A x \| A^-T b) | lev(A) |
| phase letter diag(i^[w.x + c]) | (x\|b) -> (x \| b + (w.x) w) | y(w w^T) |
| shifts, Z-phases, scalars | identity | (the kernel N) |

Here x(S) = [[I, S], [0, I]] and y(C) = [[I, 0], [C, I]]. The frame of the identity is L0 = {(0|b)}.

**The game G_Sp(m)** (Lemma N's block). Registers hold rows over C[Sp(2m, 2)]. Free steps: left multiplication of a
register by an element of P = Stab(L0) = {lev(A) y(C)}, and scalar mixing of registers. A paid move multiplies one
register by t_z. A unit "W in M" delivers x(I) y_l in W registers with M moves. Every engine unit maps to a unit of the
same shape (e3-rmbuild Lemma N: the idempotent e_N of the kernel is central and lies in the free part). The converse is
not known: the blocks of C[G] where the kernel acts non-trivially are left out.

**Proposition 0** [computed, `spg.py`, `t_setup.py` -> `out_t_setup.txt`, exact, m = 2, 3, 4, 5].
(a) All the images above are symplectic; K = t_e1 ... t_em = x(I); y(I) x(I) y(I) = w = [[0, I], [I, 0]].
(b) lev(A) and y(C) fix L0; every t_z moves L0 to a frame at distance 1; K and w move it to distance m.
(c) t_z = lev(A) t_e1 lev(A)^-1 whenever A e_1 = z (checked for every z, random A).
(d) The frames number 15, 135, 2295, 75735, and P has exactly m + 1 orbits on them, the distance classes (enumerated at
    m = 2, 3, 4; at m = 5 the distance classes have the sizes 2^(d(d+1)/2) [m choose d]_2).

*Consequence of (c)* [proved]. In G_Sp(m), as in the engine (e3-engineE Lemma E0), every move is a free conjugate of
t_e1, so the free steps absorb the direction. **With phases, "one move per direction" is no constraint, and the
Singer cycle plays no role: the family's member W in 2^m - 1 is the generic engine shape "W arrays, 2^m - 1 moves".**

## 1. Theorem T: Theorem R lives in Lemma N's block, without loss

A *certificate* of a game is a representation V of its group with a subspace Y_0 kept by the free part; its cap is
f(t) = dim(Y_0 + t Y_0) - dim Y_0 for a move t, its demand f(K), and Theorem R gives M f(t) >= W f(K) (all moves have
the same f by Proposition 0(c)).

**Theorem T** [proved]. The ratios f(K)/f(t) reached by certificates of the engine game are exactly those reached by
certificates of G_Sp(m). In particular every bound that Theorem R gives in the engine game holds in Lemma N's block:
every unit of G_Sp(m) has M >= 2W at block 2, M >= 14W/5 at block 3, M >= 41W/13 at block 4, M >= (834241/278850) W at
block 5, M >= (1188222/399685) W at block 6, and M >= 2.9447 W from block 4 on.

*Proof.* (i) From Sp to the engine. A representation of Sp(2m, 2) pulled back along G -> Sp(2m, 2) is a representation
of G; a P-stable Y_0 is kept by the free group B, because B maps into P. Caps and demands do not change.
(ii) From the engine to Sp. Let (V, Y_0) be an engine certificate. The free group B is the whole preimage of P: it
contains the kernel (shifts X^a, the squares Z^b of phase letters, the scalar i) and maps onto P (address maps give the
Levi factor, phase letters the y(w w^T), which span every y(C)) (pf-induction Facts A and B). For a frame U pick g in G
with g L0 = U and put Y_U = g Y_0. If g' L0 = U too, then g^-1 g' lies in the preimage of P, which is B, so
g' Y_0 = g Y_0: the family is well defined and G-equivariant (Y_(hU) = h Y_U). Hence delta(U, V) = dim Y_U -
dim(Y_U ∩ Y_V) depends only on the orbit of the pair (U, V), and Sp(2m, 2) is transitive on ordered pairs at each
distance. The pair (t L0, L0) is an edge and (K L0, L0) an opposite pair, so abar_1 = f(t) and abar_m = f(K). Now the
Averaging Lemma (wht35 w-bounds section 3), run inside G_Sp(m): in Map(Sp(2m, 2), V) with (g v)_gamma = v_(g^-1 gamma),
the subspace Yhat = {v : v_gamma in Y_(gamma^-1 L0)} is kept by P (it fixes L0), and f(g) = |Sp| abar_(d(gL0, L0)). So
it is a certificate of G_Sp(m) with ratio abar_m/abar_1 = f(K)/f(t). (iii) The listed bounds are Theorem R bounds of
the engine game (proofs/NOTE.md: 2W and 14W/5 Lean-checked, 41W/13 certificate-checked, the block-5 and block-6 values
and 2.9447 proved through frame families in pf-induction). []

*Remarks.*
- (ii) is the statement "every certificate can be put in frame form without loss" (meta2, upheld), read in Lemma N's
  block. What is new here is only the reading: the phase version is not blind. The test that sees phases is Theorem R
  itself, and it already lives entirely in C[Sp(2m, 2)].
- The mod-2 certificate (F_2^(2m), L0) is a representation of Sp(2m, 2) over F_2, so in G_Sp(m) too a unit whose
  mixing has odd denominators needs mW moves: **a saving must divide by 2 in Lemma N's block as well.**
- *Why deletion and rank tests are blind there and Theorem R is not.* Deletion (w8-rmexact Proposition N) restricts
  to the abelian group U_X = {x(S)} and never looks at P; it is capped by 2. Rank and spectrum tests read only the
  similarity class of the output in each representation, and K is conjugate in G to a free letter (e3-engineE
  Proposition E3). Theorem R uses a P-stable subspace, and x(I) moves L0 to an opposite frame while P fixes it.

**Proposition T2 (explicit certificates computed in Lemma N's block)** [computed, `t_avg.py`, `t_lifts.py`].
(a) Block 2, a code check of (ii): in Map(Sp(4, 2), Q^4) (dimension 2880) with a random family of unequal dimensions,
    P (order 48) keeps Yhat, and f(g) = 720 abar_d holds for all 720 elements g (`out_t_avg.txt`).
(b) Block 2: a lift of the 15 frames with no excess (abar_1 = 1, abar_2 = 2), found by search: every unit of G_Sp(2)
    has M >= 2W, so Lemma N's block has no saving at block 2 (`out_t_lifts.txt`, `lift_m2_perfect.json`).
(c) Block 3: w-bounds' lift E162, re-checked from its rows with my own code (mod-2 spans are the 135 Lagrangians of the
    seven-point symplectic space; exact integer ranks): abar = 41/35, 1307/630, 3, the file's values. So every unit
    of G_Sp(3) has M >= 105W/41 = 2.561 W: **engine 2 in 5 is refused in Lemma N's block by an explicit frame family**
    (5 x 41/35 = 41/7 < 6). The hook's 14/5 refuses it too, through Theorem T.

## 2. The family against Theorem R: an angle barrier

pf-induction and meta2 (a2-vi-stability) write every frame family's ratio at depth j as a spectral factor times an
angle factor: abar_j/abar_1 = (phi_j/phi_1)(s_1/s_j), where phi is the mean total sin^2 of the principal angles and s
the mean sin^2 per nonzero angle. The spectral factor is at most T_m(j) = 2 - 2^(1-j) (exact for m <= 40; below 8/3
for every m). Hooks have spectral factor exactly 2 - 2^(1-j) and angle factor at most 2n/(n+2) < 2 (a2-vi-stability
PROOF, the hook section). Lifts have abar_m = m and abar_1 = 1 + lambda (lambda = edge excess per edge), with
lambda >= 0.075 from block 4 on (meta2, upheld).

**Proposition H (the family sits on twice the spectral ceiling)** [derived; exact arithmetic `thresholds.py`].
(a) For every m, the full-sharing member's cost per array is (2^m - 1)/2^(m-2) = 2 T_m(m).
(b) Hence a frame family refuses W in 2^m - 1 (sharing fraction c = W/2^(m-2)) only if its angle ratio s_1/s_m
    exceeds 2/c (for m <= 40, where T_m(m) is exact).
(c) **No hook refuses the full-sharing member at any block, in any game:** a hook's value is at most
    (2 - 2^(1-m)) x 2n/(n+2) < 2 T_m(m).
(d) **No lift refuses 4 in 15:** that needs 4/(1 + lambda) > 15/4, i.e. lambda < 1/15 = 0.0667, below the floor 0.075.
    No hook does either (hooks are below 10/3 at block 4). So 4 in 15, with phases, can be refused by Theorem R only
    through a certificate of another kind with angle ratio above 2. The duality track's ceiling for all certificates at
    block 4 (3.94261, not yet checked by a skeptic) leaves room for one.
(e) At block 5 a lift refuses 8, 9, 10 in 31 exactly when its excess per edge on the whole frame graph is below
    9/31 = 0.2903, 14/31 = 0.4516, 19/31 = 0.6129. Known: whole-graph start lifts 0.88 (ratio 2.67, meta2);
    decomposed lifts of a single interval 0.279 (an interval is not a certificate; gluing costs more).

*Proof.* (a) (2^m - 1)/2^(m-2) = 2(2^m - 1)/2^(m-1) = 2(2 - 2^(1-m)). (b) A refusal needs ratio > M/W =
2 T_m(m)/c, and the spectral factor is at most T_m(m). (c) From (a) and the hook bound. (d), (e): a lift proves m/(1 +
lambda); solve m/(1 + lambda) > (2^m - 1)/W. []

*Reading.* The angle ratio s_1/s_m is the factor by which edge subspaces must meet at larger angles than opposite
ones. Every family computed so far has it below 2 (hooks < 2, start lifts 1.45 to 1.66, block-3 cuts at most 1.85;
a2-vi-stability). The full-sharing family is exactly the shape that needs more than 2, at every block. So it is not an
accident that no certificate comes near it.

## 3. What phases add to the commutative construction: two-sided targets

**Lemma A (per-array targets with phases)** [proved]. Let a commutative program (paid moves dir z, shifts, mixing; in
Theorem A's sector-0 form: C[Sym_m] with moves [z z^T]) deliver [S_l] y_l in register l, for nonsingular symmetric
matrices S_l over F_2, with M moves. Then the engine game, with phases, has a unit W in M.

*Proof.* e3-rmbuild's map phi turns the program into a slice word delivering g_l y_l with g_l in G and image x(S_l) in
Sp(2m, 2). x(S_l) sends L0 to {(S_l b | b)}, which meets L0 only in 0 (S_l nonsingular): it lies in the big cell
P w P, like x(I). So x(S_l) = p x(I) p' with p, p' in P, and g_l = b K b' with b, b' in B (B is the whole preimage of
P, Theorem T (ii)). Apply the free letter b'^-1 to input register l first and b^-1 to output register l last. The
program is linear in its inputs, so register l ends with b^-1 g_l b'^-1 y_l = K y_l. []

Without phases, address maps give only the targets S = h h^T, the non-alternating nonsingular ones (rm-theory Theorem
B). With phases the **alternating** nonsingular S are new. They exist only at even m: none at blocks 3 and 5, 28 at
block 4 (of 448 nonsingular symmetric 4 x 4 matrices).

**Proposition A2 (alternating targets are expensive)** [proved; computed `t_alt.py` -> `out_t_alt.txt`].
(a) An alternating nonsingular S is a sum of at least m + 1 matrices z z^T (their sum's diagonal is the sum of the z,
    which must vanish, so the z are dependent, and m independent ones are needed for rank m). So a lone array owed it
    needs m + 1 moves in the commutative program (e3-rmbuild's rank-depth filtration; computed: 3 at m = 2, 5 at m = 4).
(b) If every array is owed the same alternating S, then each form Q = (a.z)(b.z) with a^T S b = 1 is odd for every
    array (<Q, S> = a^T S b) and has support 2^(m-2); averaged over these forms, deletion gives
    M >= W (2^m - 1)/2^(m-2): exactly the family's cost. So the family with one alternating target is deletion-tight
    on every such form (computed at m = 2, 4, 6: all forms odd, support 2^(m-2), every direction covered equally).
*Consequence.* At block 5 (odd) Lemma A adds nothing. At even blocks it adds targets that cost more alone and put the
family on a tight deletion bound, where the zero-block test has teeth. Mixed target sets were not examined.

## 4. Calibration summary

| Case | Theorem R in Lemma N's block | status |
|---|---|---|
| plain schedule, every block | a unit (K = t_e1 ... t_em, Proposition 0(a)); m per array, above every bound | passes |
| E8 unit (block 45, 240 data arrays, at most 4 x 120 moves saved) | at least 43 moves per data array; its image in C[Sp(90,2)] is a unit; every bound here is at most 3 per array | passes |
| any saving at block 2 | refused by the perfect lift (2W) | refused |
| engine 2 in 5, block 3 | refused by E162 (2.561 W) and by the hook (2.8 W) | refused |
| 4 in 15, block 4 (3.75 per array) | best bound 41/13 = 3.154; no hook, no lift can ever refuse it (Proposition H) | not refused |
| 8 in 31, block 5 (3.875 per array) | best bound 2.992; no hook can ever refuse it; a whole-graph lift needs excess < 0.2903 | not refused |
| 9, 10 in 31 | lifts need excess < 0.4516, 0.6129 (start lifts 0.88) | not refused |

## 5. Phase letters at the two ends (restart session, 15 Oct, 15:21 to 17:40 EDT)

Call a unit *ends-only* when every address map and phase letter acts on an input register before the first paid
move or on an output register after the last one, and the middle is a slice program (paid moves, shifts, scalar
mixing).

**Corollary A4 (where phase letters can matter)** [(a) argued; (b), (c) proved].
(a) The ends-only units are exactly the Lemma A programs: register l is owed x(S_l) with S_l nonsingular symmetric,
    delivered by a commutative program in Theorem A's sector-0 form. One direction is Lemma A. For the other, the
    middle delivers out_l = sum over l' of c_(l,l') in_l', with c in C[G_slice], and b_l c_(l,l') b'_l' = delta K.
    So c_(l,l) = b_l^-1 K b'_l^-1 is a single group element, which therefore lies in the slice group. Its image
    x(S_l) lies in U_X and in P x(I) P = P w P, which for x(S) means that S is nonsingular.
(b) Over F_2 a nonsingular symmetric matrix is congruent to I unless its diagonal is zero, and the alternating ones
    form one congruence class, which exists only for even m (checked by orbit enumeration at m = 2 to 5,
    `tile/congr_classes.py`: classes of sizes 1 + 3, 28, 28 + 420, 13,888). So at odd blocks the ends-only units with phase letters
    are exactly rm-theory's per-array-target programs without phase letters. At block 3 Theorem C3 then forbids
    every saving. At block 5 the ends-only phase version of 8, 9 and 10 in 31 is exactly the system of rm2-solve
    and rm2-exact.
(c) At even blocks the new targets are the alternating matrices: 1 at block 2, 28 at block 4, and 13,888 at block 6
    (|GL(6,2)| / |Sp(6,2)| = 20,158,709,760 / 1,451,520).

So phase letters can matter only (i) between paid moves, or (ii) at even blocks, through alternating targets.

**Tools** (all in `tile/`, each calibrated before use; Colab jobs under `jobs/`).
- `fam_enum.py`: tests (S), (D), (C), (Z) of t_small_alt.py for one direction multiset. It fixes one target per
  orbit of the multiset's stabiliser in GL(4,2), which loses nothing (see its docstring), and deduplicates survivors
  up to the stabiliser. Cross-checks: on two test classes it gives the same classes as t_small_alt.py (6 and 24); on
  4 in 14 the orbit sizes of its alternating classes add up to 77,840 tuples, exactly altA's t_small_alt.py count
  (`tile/orbit_count.py`). Its sharded runs reproduce the unsharded one.
- `tile_alt.py`: rm2-exact's rank-cover tiling test (rcdp.c, copied unchanged; rank E[D_F, D_F] = |D_F| - r_F with
  rectangle bounds), for targets of any kind. It holds for alternating targets because rm-theory 8.5 uses a target
  only through p_l(Q) = <Q, T_l>. Every exclusion is re-verified independently. For that the verifier recomputes
  <Q, T> as sum_j Q(v_j) over a minimal decomposition T = sum v_j v_j^T; an alternating T needs m + 1 vectors.
  Controls: plain programs with two alternating targets pass at every form; impossible single-array cases are
  refused; and the six zero-block exclusions printed by altA are re-excluded with verified certificates.
- `ilp_alt.py`: rm2-exact's stronger ILP bound (rc_ilp.py, copied unchanged), with the same verifier change.
- `lm_alt.py`, `lm_batch_alt.py`: rm2-solve's Levenberg-Marquardt (glm.py, copied unchanged). Calibrated in every
  run: the plain 4 in 18 with two alternating targets is found (2.3e-27), and a zero-block-excluded 3 in 11 is
  refused (288).

**Proposition A3 (block 4, ends-only units with phase letters)** [computed with exact tests; LM is evidence only].
(a) *No ends-only unit, with or without phase letters, for 2 in 7, 3 in 10, 3 in 11 or 4 in 13.* Exhaustive over
    every direction multiset up to GL(4,2) and every multiset of nonsingular symmetric targets. Tuples with all
    targets congruent to I are rm-theory's runs (`rm-theory/out_t_free_b4_7.txt`, `out_t_small_b4_W3_M10*.txt`,
    `jobs/mtD/remote-out/ssh/out_W3_M11_nofilter.txt`, `out_W4_M13.txt`). Tuples with at least one alternating target
    are mine (`out_alt_b4_W2_M7.txt`, `out_alt_b4_W3_M10.txt`, `jobs/altA/remote-out/ssh2/out_alt_b4_W3_M11.txt`,
    `out_alt_b4_W4_M13.txt`). The free-direction skip used for 4 in 13 stays valid with alternating targets, because
    such a target needs m + 1 rank-ones (Proposition A2(a)).
(b) *Consequence.* The engine's bound 41/13 leaves exactly two above-bar shapes at block 4 with at most 4 arrays,
    3 in 10 (0.132) and 4 in 13 (0.150): W arrays need at least 41W/13 moves, and the bar allows at most 3.482 W. So,
    with address maps and phase letters used only at the two ends, **block 4 has no above-bar unit with at most 4
    arrays.**
(c) *4 in 14 (0.096).* In altA's full run, only one direction class survives with an alternating target: 14
    distinct directions, each used once, which leaves out one direction. In that class (`fam_enum`, all 448 targets,
    `jobs/altE/`):
    - 411 target classes, up to the point stabiliser of order 1,344, pass S+D+C+Z. 328 have no alternating target
      (this class was never reached by rm-theory's partial 4-in-14 run), and 83 have at least one. Every one of them
      has a tight form.
    - The tiling test excludes 393 of them, with verified certificates (`tile/out414/`). That leaves 18: 12 with no
      alternating target, and 6 with one (4 classes), two (1) or three (1).
    - The ILP bound excludes none of the 6 alternating ones (`jobs/altK*/`, all forms, 367 to 534 s each; the
      margins are the DP's).
    - LM, 3 time orders x 2 restarts on each of the 18 (`jobs/altJ*/`), stops at 494.7 or above in every case:
      two or three arrays served exactly and the rest not (the trivial-split pattern); no near-zero point.
(d) *4 in 15, each direction once (the family member, 0.047).* With at least one alternating target the survivors
    of S+D+C+Z form exactly 1,612 classes up to GL(4,2) (`tile/surv_fe_alt_415_all.txt`): 1,428 with one alternating
    target, 163 with two, 15 with three and 6 with four. 114 of them have no tight form. The run fixed one
    alternating target (one orbit) and went through all C(450, 3) = 15,086,400 triples of the other three; deletion
    removed none of them, and 1,122,378 passed. It was split 3 ways on the first free target
    (`jobs/altH/remote-out/ssh3/`, `jobs/altH2/remote-out/ssh5/`). Each third alone produced the same 1,612 classes,
    which is a consistency check, since the thirds enumerate disjoint sets of tuples.
    The classes behave like the no-phase canonical sets of rm2-exact:
    - all 114 classes with no tight form pass the tiling test at every form, with the blind margin pattern
      (0: 127 forms, 1: 128, 2: 384, 3: 384). Three of them have all four targets alternating (`jobs/altL/`,
      `jobs/altN/`);
    - of the classes with tight forms, a sample was tiled (`tile/sample415_j200.txt`, one first-level branch,
      1,170 classes, all among the 1,612): 360 of the first 641 pass as well (`jobs/altL/remote-out/ssh/`; stopped
      at the time-box);
    - LM on 9 of the no-tight classes (`jobs/altM/`) stops at 384.0, which is three arrays served exactly with 12
      moves and the fourth not, or above.
(e) *Not completed.* The negative control 5 in 15 (impossible in every game by 41/13) needs about 1.7e9 leaves
    in this enumeration, so I stopped it; whether S+D+C+Z plus tiling refuse it is unknown.

*Reading.* Phase letters at the two ends change nothing that the exact tests can see at block 4. Small shapes are
excluded whatever the targets. At 4 in 14 and 4 in 15 the alternating targets survive only where the
non-alternating ones survive too: the sets with no tight form are invisible to every test, as in the no-phase game. So the
family with phase letters at block 4 is open in the same way, and for the same reason, as without them. The
numerical search gives the same trivial-split floors.

