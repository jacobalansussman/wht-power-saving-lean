<!-- header: frame-packing-two-arrays -->
# The frame game, its lift and the position test (the build track of packing round 1)

**About this file.** This is the working proof file as an AI agent wrote it during the round of 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent. Its opening lines say that no skeptic had read it; that was true when it was written. What happened afterwards is below, and where a later skeptic's verdict differs from this file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean. File names inside it refer to the programs and outputs now in `../programs/` and `../data/`; [README.md](README.md) maps the working names.

**Status now.** Upheld by one AI skeptic for its theorems and its search tallies; one sub-claim refuted: the "closest design" fails the position test exactly, for every choice of coefficients (its near-pass was a limit in which the shared move stopped packing). The skeptic also found that the lift works at f = 3 mod 4 too, and that real landings were at most 6 or 7 per design, not 12. Its record: [frame-game-sandwich-skeptic.md](frame-game-sandwich-skeptic.md). The open part of this file (the even column counts) was settled later: see [frame-game-completion.md](frame-game-completion.md).

---

# PROOF.md (fp-build): the frame game, its lift, and the screens used to build packing-first units

fp-build, packing round, 15 October 2026. No skeptic has read this. Labels: [proved] (proof written here), [computed]
(program and output in this folder), [read] (from a named file, not re-derived), [cited], [argued].
Programs: `fpa.py` (frame-game engine, exact end solve, psi screen, constructive charts), `t_cal.py`, `t_blocks.py`,
`t_o3ctl.py`, `t_lift.py`, `t_lift_dbg.py`, `t_psi.py`, `t_psi_repair.py`, `t_psilin.py`, `t_split.py`, `t_psibil.py`,
`design.py`; outputs `out_*.txt`, Colab outputs in `remote-out/`.

## 0. Setting

The engine game at block m (sl-combinatorics PROOF.md section 0; theory-guide RULES section 0). G is the group generated
by the letters: paid `dir z = ca 1 + cb X^z`, free shifts, address maps and phase letters `diag(i^((z.x + c) mod 2))`;
B is the free subgroup. A word is right at f digit columns when every letter g acts as `g^(x)f` and the delivering
registers end with `K^(x)f` applied to the inputs. A unit is right at every f.

Symplectic images. Every g in G acts on Pauli labels (a|b) (a = x-part) through `gbar` in Sp(2m,2). `dir z` maps to
the transvection `T_z: (a|b) -> (a + (z.b) z | b)`; free letters map into the parabolic `P` = stabiliser of the raw
frame `L0 = {(0|b)}`; K maps to `Kbar = [[I, I], [0, I]]`. The level of `lambda` in Sp is the rank of its x-block on
L0 (paid length, sl-combinatorics Lemma 1.1(d)); its position is `lambda^-1 L0`, a point of the dual polar graph
Gamma_m. [computed: `t_lift.py` checks that the monomial lifts used below and every `dir z` have exactly these images.]

## 1. Blocks of the engine algebra

**Lemma 1.1** [read: e3-rmbuild Lemma N; computed at m = 2]. Let N be the Pauli matrices times the scalars {1, i, -1,
-i}. N is normal in G, N lies in B (shifts; `Z^b` is a squared phase letter; `i = i^[0.x + 1]` is a phase letter), and
`G/N = Sp(2m,2)`. At m = 2 the group G modulo scalars has 11,520 = 16 x 720 elements and every fibre over Sp(4,2) has
16 elements (`t_blocks.py`, `out_blocks.txt`).

**Lemma 1.2 (blocks, and the skeleton is shared)** [proved]. Let `s = i.1` (a free letter, central in G) and
`e_k = (1/4) sum_j i^(-jk) [s^j]` for k mod 4. These are central idempotents of C[G], they lie in C[B], and at f
columns `e_k` acts as 1 if f = k (mod 4) and as 0 otherwise. Hence:
- a word is right at every f if and only if it is right in each block `C[G] e_k`;
- the free coefficients of a word can be chosen independently in each block: if for each k a word with the same
  SKELETON (same steps `c_r := sum_s beta_s c_s` and `c_r := dir z . c_r`, same registers, same order) is right in
  block k with coefficients `beta^(k)` in `C[B] e_k`, then the word with coefficients `sum_k beta^(k)` is right at every
  f, with the same paid moves.

*Proof.* `[s]^(x)f = i^f`, so `e_k` evaluates to `(1/4) sum_j i^(j(f - k))`. Centrality: s is central. The second
statement: the e_k are orthogonal central idempotents, so a word's image in block k only sees `beta e_k`, and paid
letters commute with every e_k. []

**Proposition 1.3 (the faithful block is a matrix algebra over C[Sp])** [proved; computed at m = 2]. Block 1 is the
quotient `C[G]_mu4` in which `[i g] = i [g]`. This is exactly the convention of `slack/sl-combinatorics/ledger.py`
(group elements modulo {1, i, -1, -i}, the scalar moved into the coefficient). The map

    Phi(g) = g (x) gbar :  C[G]_mu4  ->  M_n(C) (x) C[Sp(2m,2)],   n = 2^m,

is an algebra isomorphism. Under it:
- `C[B]` goes onto `M_n(C) (x) C[P]`;
- `dir z` goes to `(dir z (x) 1)(1 (x) T_z)`, a free element times `1 (x) T_z`;
- `K` goes to `(K (x) 1)(1 (x) Kbar)`.

*Proof.* Phi respects `[i g] = i [g]` and is multiplicative because both factors are. Over each `lambda` in Sp, the
fibre of G/{scalars} is `ghat . (Paulis)`, and the n^2 Pauli matrices are linearly independent; so Phi is injective
fibre by fibre and fibres go to different `lambda`. Both sides have dimension `n^2 |Sp|` (Lemma 1.1). The fibres over
P are `phat . (Paulis)`, whose span is `M_n (x) phat-bar`. [`t_blocks.py`: at m = 2 all 720 fibres have linearly
independent matrices.] []

*Remark 1.4 (what a "formal identity" in the ledger convention certifies).* By Lemma 1.2, an identity in `C[G]_mu4`
is rightness at f = 1 (mod 4) only. Raw-array checks at even f (and f = 3) remain necessary. Every unit claim in this
round's pipeline uses raw arrays at f = 1..6 as well, so nothing is affected; the remark is for future readers of
ledger-style checks.

*Remark 1.5 (what the other blocks look like)* [argued; dimension counts checked by hand at m = 1, not computed].
- At even f the Pauli matrices commute after the f-fold action, because `(XZ)^(x)f = (ZX)^(x)f` for f even. So in
  blocks 0 and 2 the Pauli part is a commutative algebra. Its 4^m characters are free idempotents, and the paid
  moves permute them.
- By Clifford theory (Mackey), the blocks should then be:
  - block 0 (f = 0 mod 4): `C[Sp]` (this is the frame game itself: trivial character, Lemma N), plus
    `M_(4^m - 1) (x) C[Sp_v]`, where Sp_v is the stabiliser of a nonzero vector;
  - block 2 (f = 2 mod 4): `M_(N+) (x) C[O+(2m,2)]` plus `M_(N-) (x) C[O-(2m,2)]`. Here the characters are the
    quadratic forms of Arf invariant 0 and 1, with `N+- = 2^(m-1)(2^m +- 1)`; possibly twisted by cocycles.
- At m = 1 the dimensions add up to `|G/mu4| = 24` in both blocks.
- So a frame-game unit becomes a full unit exactly when the same skeleton can also be solved in the orthogonal-group
  games and in the vector-stabiliser game. Those are the next things to test for any frame-game find.

## 2. The frame game and the sandwich

**Definition 2.1.** The *frame game* at block m: registers hold W-tuples of elements of `A = C[Sp(2m,2)]`; free steps
are `C[P]`-linear combinations of registers; a paid step is `c_r := T_z c_r`; input l starts as `[1] e_l`; register r
delivers when it ends as `Kbar e_l` (times a nonzero scalar). Terms of a content are *walkers* `(l, lambda)`.

**Theorem 2.2 (sandwich)** [proved]. For every block m, every W, M and every skeleton:

    engine units (right at every f)  ⊆  frame-game programs  ⊆  engine words right at every f = 1 (mod 4),

with the same W, the same M and the same skeleton at both inclusions.

*Proof.* (Left) The Pauli average `e_N = |N|^-1 sum_{n in N} [n]` is a central idempotent in C[B] (N normal, N in B),
and `C[G] e_N = C[G/N] = C[Sp]`. Under this map free letters go into `C[P]`, `dir z` to `T_z` and K to `Kbar`. So the
image of an engine unit is a frame-game program with the same steps (e3-rmbuild Lemma N).

(Right) `Lambda = Phi^-1 (1 (x) .)` is an algebra homomorphism `C[Sp] -> C[G]_mu4`.
- It sends `C[P]` into the image of `C[B]`. Concretely, `Lambda(p) = Pexp(phat^-1) . phat`, where phat is a monomial
  lift and Pexp is the Pauli expansion, a formal combination of free letters.
- It sends `T_z` to `c_z . dir z`, one paid letter times the free `c_z = ((1-i)[1] + (1+i)[X^z])/2`. This is
  e3-rmbuild Lemma R's corrected move.

Applying Lambda step by step to a frame-game program gives a word with the same skeleton and the same paid moves. Its
block-1 image is `1 (x)` (the program), so it delivers `1 (x) Kbar`. One more free factor `Pexp(K)` turns this into
`Phi^-1(K (x) Kbar) = [K]`. By Lemma 1.2 the word is right at every f = 1 (mod 4). []

[computed, `t_lift.py`, `out_lift.txt`]:
- **Conventions.** Monomial lifts of 20 random p in P, and every `dir z`, have exactly the frame-game images.
- **Homomorphism.** Lambda is an exact homomorphism on raw arrays at f = 1 and f = 5. It fails at f = 2, 3 and 4
  (relative error about 1).
- **An actual lift.** The order-3 control's frame-game readout, found by `fpa.end_solve`, lifted and replayed on random
  complex arrays:
  - right at f = 1, 3 and 5 (errors 2e-16, 5e-16, 7e-16);
  - wrong at f = 2 and 4 (errors 1.7).

- **Lifts with phase letters** (`t_lift2.py`, `out_lift2.txt`; block 2). Two frame programs of the form
  `q p3 T p2 T p1` (random free letters with phases), lifted letter by letter:
  - exact at f = 1, 5 and 9 (errors at most 5e-16);
  - wrong at f = 2, 3, 4 and 6.

  So with phases, f = 3 mod 4 is not automatic either; only f = 1 mod 4 is guaranteed, as the theorem says.

  The control word written with its own engine letters is right at every f (o3-build). So the even-f blocks need
  coefficients of their own, and Lemma 1.2 says exactly that they may have them: the skeleton is shared, the
  coefficients are not.

**Corollary 2.3** [proved, given the hook's f]. Every lower bound proved by a certificate at a single f = 1 (mod 4)
holds in the frame game, with the same W and M. Conversely, every exclusion proved in the frame game is an exclusion
for engine units. The hook certificates live in `V (x) Lambda^k V`, inside f = k + 1 columns [argued from bm.py's
definition, not re-derived]. That gives:
- 14W/5 at block 3 (hook at f = 5);
- 41W/13 at block 4 (k = 8, f = 9).

So both floors hold in the frame game. In particular 2 in 5 at block 3 and 2 in 6 at block 4 are impossible there,
and they serve as the negative controls of section 4.

**Remark 2.4 (lanes; the Pauli part does not add packing capacity)** [proved].
- Read block 1 through Proposition 1.3 one column of `M_n` at a time. An engine word, W in M, right at f = 1 (mod 4),
  becomes a frame-game program for nW inputs. Each register carries n coordinates ("lanes"), and a paid move applies
  `T_z` to all n lanes at once, so the program has nM single applications.
- Conversely, Lambda embeds every frame-game program, N in C, as an engine word N in C (the n lanes carry identical
  copies).
- So the best moves-per-array at f = 1 (mod 4) equals the best T-applications-per-input in the frame game. The 2^m-fold
  multiplexing that the matrix factor seems to offer buys nothing: every lane is an input that needs its own
  `Kbar`.

## 3. Packing in the frame game

**Lemma 3.1 (rigid bodies)** [proved]. Every step multiplies all walkers of a register on the left by the same element.
So for two walkers `lambda_i, lambda_j` of one register:
- the relative element `lambda_j^-1 lambda_i` never changes while they stay together;
- their positions satisfy `X_j = (lambda_j^-1 lambda_i) X_i` at all times, so each walker's trajectory is the image of
  the other's under one fixed isometry of Gamma_m;
- a paid step `T_u` moves the walker at `X = lambda^-1 L0` to the neighbour `T_(lambda^-1 u) X`;
- two walkers of one register never cancel each other later: their relative element is fixed and not 1.

A packed move is therefore a rigid motion of several walkers, and every cancellation of junk is between different
stored states.

**Lemma 3.2 (landing)** [proved; constructive in `fpa.find_chart_c`].
- *Statement.* A walker `(l, g)` of a stored state, injected with a free chart mu before a stretch of letters w of some
  register, ends at a prescribed position X if and only if `d(g^-1 L0, X) = level(w)`. When it does, the chart mu can
  be computed.
- *Proof.* The walker ends at `g^-1 mu^-1 w^-1 L0`. So we need `mu(gX) = w^-1 L0`. P acts transitively on each sphere
  around L0, and `d(L0, gX) = d(g^-1 L0, X)`, while `d(L0, w^-1 L0) = level(w)`. Explicitly, `normal_P` brings any
  Lagrangian at distance d to the standard one: a Levi element straightens its x-part, and a phase removes the
  symmetric remainder. []
- *Checked* [computed]: 200 random cases at each of blocks 3 to 6, all correct.

**Lemma 3.3 (one shared move with single walkers: the junk has no partner)** [proved; a special case of sl-algebra
Corollary 2′ and t2-pair Lemma S/B, stated in walker terms because the search below starts from it].
- *Setting.* Two lineages share exactly one move (the trunk). Every stored state has one walker per input.
- *Claim.* The y2 walker that branch A carries to its end cannot be cancelled by any other stored state. The same holds
  for branch B's y1 walker.
- *Proof.* By Lemma 3.1, A's y2 walker follows the image under a fixed isometry sigma of y1's geodesic from the trunk
  to L1. So its path is a geodesic that starts with the ascending trunk step `P_a2 -> P_a2+1` of y2's own chain.
  - It cannot end on y2's prefix chain: that would need `d(P_a2, X) = d(P_a2+1, X) + 1` and
    `d(P_a2+1, X) = d(P_a2, X) + 1` at once.
  - The only branch-B state with a y2 walker at the same position is the one after the same number of moves. Its two
    walkers then share positions with A's, so the two states differ by one free factor (`q_A q_B^-1` lies in P), and
    using one to cancel the other also cancels A's target.
  - States before the trunk carry one input only. []

**Proposition 3.4 (the psi screen: positions with scalar coefficients)** [proved].
- *Definition.* Let `psi : A^W -> C[Gamma_m]^W` send a walker `(l, lambda)` to `[lambda^-1 L0] e_l`, extended
  linearly.
- *Key identity.* `psi(beta x) = eps(beta) psi(x)` for every beta in C[P], where eps is the augmentation, because a free
  factor never changes a position.
- *Screen.* A frame-game unit, and so every engine unit by Theorem 2.2, must satisfy: for each input l,
  `[L1] e_l` is a scalar combination of the psi-images of the stored states (inputs and move outputs).
- RULES D2's top augmentation is this condition summed over the level-m sphere. Psi keeps every position separately,
  so it is strictly finer.

[computed, `t_psi.py`, `out_psi.txt`]:
- the plain 2 in 2m passes at blocks 3 and 4;
- the same word cut by one move fails for the cut array;
- 300 single-trunk designs without injections at each of blocks 3 and 4 fail it for both outputs (as Lemma 3.3
  predicts).

*Variants used in the search.*
- **Linearised psi.** Every injection scalar is a symbol, and every (state, monomial) pair becomes its own column. If
  this fails, no choice of injection scalars and readout scalars passes psi.
- **Bilinear psi.** The scalars are unknowns, minimised numerically. This is evidence only.

**Lemma 3.5 (repeated branch moves close one output, never both by this route)** [proved].
- *Setting.* Branch B's first t letters are exactly A's first t letters, so before any landing `B_t = A_t`. A
  single-walker state pi of input 2 is landed into B before its move k <= t, on the frame X of A's final y2 walker,
  arriving at B_t (Lemma 3.2 needs `d(pos pi, X)` = the paid length of B's letters k..t).
- *Claim.* Then `A_final - c p (B_t - A_t) = Kbar y1` for a suitable scalar c and free p. So output 1 closes.
- *Proof.* `B_t - A_t` is exactly the landed walker. Its frame is X, so one free p and one scalar match it to A's
  junk term. []
- *The catch.* The landed walker stays in B: it is in `B_final`, on a frame of its own. So output 2 now has two junk
  terms to cancel, B's y1 walker and the landed one. Landing the mirror image into A, to cancel B's y1 walker, puts a
  new walker into `A_final` in the same way. [computed: designs in which one output is exact and the other is left
  with junk without a partner occur in the runs with the constructive chart finder: 17 of 330 at block 3, 4 of 300 at
  block 4 with two shared moves, 8 of 398 at block 5 with two shared moves, 3 of 153 at block 6. Both outputs at once:
  never. Four such one-output solutions were re-verified term by term (`t_oneyes.py`, `out_oneyes.txt`).]

## 4. The exact pipeline and its calibration

`fpa.end_solve(states, target)` decides exactly whether a target content lies in the C[P]-span of the stored states.
1. *Closure.* From the target rows, add every column (state r, free element p) that produces a row, then that column's
   other rows, and repeat.
2. *Peel.* A non-target row with a single live producer forces that column to 0.
3. *Solve.* Exact Gaussian elimination over Q (Fractions); any solution found is re-verified term by term
   (`verify_combo`).

"cap" means the closure passed 6,000 columns (or 20,000 rows) and was not decided.

| case | block | expected | result | file |
|---|---|---|---|---|
| plain W = 1, 2 | 2-5 | yes | yes | out_cal.txt |
| plain minus the last move of array 0 | 2-5 | no | no (target row unproduced) | out_cal.txt |
| phases inserted between K's directions | 2-5 | no (not K) | no | out_cal.txt |
| order-3 control C6O3 / C8O3 (r1 = y1 + y2, r2 = y1 - s y2) | 3, 4 | yes, raw-gauge P = m | yes both outputs; P = 3, 4; solutions re-verified | out_o3ctl.txt |
| order-3 control cut by one move | 3, 4 | no | no, no | out_o3ctl.txt |
| single trunk, no injections (Lemma 3.3) | 3, 4 | no | no (all designs), psi fails | out_psi.txt; st.py |
| 2 in 5 at block 3 (14W/5 holds in the frame game) | 3 | never yes, yes | never | out_design_m3_*; remote-out (cal job) |
| 2 in 6 at block 4 (41W/13 holds in the frame game) | 4 | never yes, yes | never | out_design_m4_try; remote-out (cal job) |
| 3 in 7 at block 3 (14W/5) | 3 | never all three | never (300 designs) | remote-out/out_designw_3in7_b3_s92.txt |
| 2 in 3 at block 2 (Theorem 2W, an engine theorem; not shown for the frame game: consistency check only) | 2 | never yes, yes | never (200 designs) | out_design_m2_s3.txt |
