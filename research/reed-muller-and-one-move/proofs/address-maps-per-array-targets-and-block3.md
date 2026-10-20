# Proofs: address maps as per-array targets (Theorem B) and no saving at block 3 without phase letters (Theorem C3)

*Status of this file.* Section 8 of the working note of the Reed-Muller project's general-m track (15 October 2026), as written. An independent AI skeptic re-derived Theorem B and recomputed the certificate behind Theorem C3 with its own exact code (`programs/address-maps-and-small-lift/independent-check/vk_face3.py`) and found no error. No person has read it and no proof checker has checked it. Statements inside marked "No skeptic has read this" describe the moment of writing. Paths that begin with `checks/` refer to the private working folders of the research and are not part of this repository; private machine paths were removed. Where this file and the section README disagree, the README gives the later verdict. The labels inside ([proved], [computed], [argued], [guess]) are the working note's own: "[proved]" means that a proof is written there, not that it is proved in Lean or checked by a person; the statuses this revision uses are the ones in the section README.

---

### 8.1 Theorem B

(a) Let the multi-target program take inputs `x_1..x_W` in the slice and deliver `[T_l] x_l`, with `T_l = h_l h_l^T`.
- Its output register l holds `sum_l' c_(l l') x_(l')` with `c_(l l')` in `C[Sym_m]`.
- Correctness as a formal identity means `c_(l l') = [T_l]` if l' = l, and 0 otherwise.
- In the no-phase game with H, first apply the free letter `[h_l]` to input register l, so `x_l = [h_l] e_l`.
- Run the same moves and mixings; they are letters of `C[Sym_m x| H]`. Output l is then `[T_l][h_l] e_l`.
- In `Sym_m x| H`, `[h]^-1 [S] [h] = [h^-1 S h^-T]`. So `[h_l]^-1 [T_l][h_l] = [h_l^-1 h_l h_l^T h_l^-T] = [I]`.
- Applying the free letter `[h_l]^-1` to output l therefore delivers `[I] e_l`. The moves are unchanged in number and
  direction.

(b) This is w8-rmexact's Theorem MT (`../../moonshot/w8-rmexact/PROOF.md`): left multiplication on the free right
`C[Sym]`-module `C[Sym x| H]` with basis H.

(c) By (b), a ratio r = M/W of the H-game is achieved by a multi-target program with all of T(H), equally used. By
(a), that program is an H-game unit with the same ratio. Conversely, any multi-target ratio with targets in T(H) is an
H-game ratio by (a). []

*Formal check* (`t_lift.py` -> `out_t_lift.txt`, blocks 3 and 4, exact group-algebra arithmetic in `Sym_m x| Z_n`).
- Lifted correctly: frames for I and for `sigma sigma^T`, and a 4-move non-frame decomposition of `sigma sigma^T`.
  Each delivers `[I]`.
- A wrong relabelling (`[sigma^2]` for target `sigma sigma^T`) does not.

### 8.2 Field model (F1 to F3)

- Let f_1..f_m be a self-dual basis (it exists because the trace form is non-alternating; `gf.py` finds one by
  greedy search). The coordinate map `x -> (Tr(x f_j))_j` sends the trace form to the dot product.
- **Characteristic vector.** For all x, `Tr(x^2) = Tr(x) = Tr(x · 1)`, so 1 is the characteristic vector of the trace
  form. In an orthonormal basis that vector is `sum_j f_j`.
- **Moore matrices.** Let M be the matrix `(f_j^(2^a))` (row a, column j). Self-duality says `M^T M = I`, so also
  `M M^T = I`. That gives `sum_j f_j^(2^a) f_j^(2^b) = delta_(ab)`; with a = 0 and b = i, `sum_j f_j^(1+2^i) = 0` for
  i ≠ 0.
- **The parity.**
  - `tr Q = sum_j Q(f_j) = Tr(b sum_j f_j) + sum_i Tr(c_i sum_j f_j^(1+2^i)) = Tr(b)`.
  - For copy k, `(Q o sigma^-k)` has linear coefficient `b alpha^-k`, and its parity is `Tr(b alpha^-k)`.
  - As k runs over Z_n, `b alpha^-k` runs over all nonzero field elements when b is not 0. Exactly 2^(m-1) of them
    have trace 1.

### 8.3 The placement theorem

- **(a).** Copy k of a move in direction z acts in direction `alpha^-k z`. As k runs over Z_n these are all nonzero
  points. So `|D_Q| = sum over moves of wt(Q) = M wt(Q)`.
- **(b).** Write the target as `T = g g^T`, with frame g_1..g_m. Then `n_odd(Q) = #{y != 0 : Q~(y) = 1}` with
  `Q~(y) = sum_j Q(y g_j)`, a binary quadratic form in y vanishing at 0.
  - Every nonzero quadratic form vanishing at 0 has weight between 2^(m-2) and 3·2^(m-2), since its complement is a
    nonzero word of RM(2, m) of weight at least 2^(m-2).
  - So for Q ≠ 0, `n_odd(Q) <= 3·2^(m-2) <= 3 wt(Q)`, and `M wt(Q) > 3W wt(Q) >= W n_odd(Q)` whenever M > 3W.
  - Hence no violation and no tight form.
- **(c).** By F3, `n_odd(Q) = 2^(m-1) [b != 0]`, and the forms with b ≠ 0 have weight at least 2^(m-2), with equality
  for the flats. So rho = 2.
- **(d).** In w8-rmexact's normalisation (per virtual array, c = (M/W)/n per direction), every slack satisfies
  `slack(G) >= (t 2^(m-2) - 2^(m-1))/n` with t = M/W.
  - That is at least c = t/n once `t >= 2^(m-1)/(2^(m-2) - 1)`. Then every cell bound equals c.
  - The cover row of a form F with w directions then reads `w^2 c >= slack(F)`. It holds because `slack(F) <= t w/n`.

### 8.4 The splitting

- Z_n has odd order, so F_2[Z_n] is semisimple and Alt has a Singer-invariant complement in Sym.
- **Uniqueness.** Two complements differ by a Singer map from Sym/Alt (the natural module, cyclotomic class of 1) to
  Alt. By the field model Alt is dual to the pure forms, which are the classes 1 + 2^i.
  - Those classes are distinct from the class of 1: the coset `{2^s(1 + 2^i)}` contains a number with two binary
    ones.
  - So the map is 0, and the complement is unique.
- U (F4) is a Singer-invariant complement, so it is that one. []

### 8.5 Theorem 10 and Theorem 1 with a target per array

- ex2-families Lemma 1 does not use the targets until its last line. So a multi-target program exists exactly when
  `X_D := Ã_D (I + Γ_D)^-1 B̃_D = diag_l(p_l(Q))` for every form Q, with D = D_Q and `p_l(Q) = <Q, T_l>`.
  (At Q the program must deliver `(-1)^(p_l)`, and `I - 2X = diag((-1)^(p_l))`.)
- **Ranks.** Put K = B̃ Ã and E = I + Γ - K. Then `rank E[D, D] = |D| - r_Q`, with `r_Q = sum_l p_l(Q)`.
  - `E[D, D] = (I + Γ_D)(I - C Ã_D)` with `C = (I + Γ_D)^-1 B̃_D`, and I + Γ_D is invertible (unit triangular).
  - `v -> Ã_D v` maps the kernel of `I - C Ã_D` one-to-one onto the kernel of `I - Ã_D C`, with inverse `w -> C w`.
  - `Ã_D C = X_D` is the 0/1 diagonal of the odd arrays, so that kernel has dimension r_Q.
- **Tight forms** (|D| = r_Q >= 1): `E[D, D] = 0`, that is `K[D, D] = I + Γ_D` (unit diagonal, zero above it in time
  order).
- **Cliques.** If Y is a set of live moves, any two of which lie in a common tight support, then `K[Y, Y]` is unit
  lower triangular, so `|Y| <= rank K <= W`. This is Theorem 1 word for word.

### 8.6 The free-direction lemma

**Statement.** Let y_1..y_d be the distinct directions of a multi-target program, used μ_1..μ_d times, and suppose
the matrices `y_i y_i^T` are linearly independent over F_2. Then `M >= sum_l rank T_l`, which is mW when every target
is congruent to I.

*Proof.*
- The pairing of symmetric matrices with forms is perfect, and `Q(y) = <Q, y y^T>`. So by independence there is, for
  each i, a form Q_i with `Q_i(y_j) = [i = j]`. Its support D is the μ_i copies of y_i.
- Every target lies in the span of the `y_j y_j^T`: a form vanishing at every move has D empty and X = 0, so every
  p_l is 0 there. Write `T_l = sum over j in S_l of y_j y_j^T`; S_l is unique.
- `p_l(Q_i) = [i in S_l]`. Deletion (the rank in 8.5 is at least 0) gives `μ_i >= #{l : i in S_l}`.
- Summing over i, `M >= sum_l |S_l| >= sum_l rank T_l`, since a sum of |S_l| matrices of rank one has rank at most
  |S_l|. []

**Where dependencies come from** [proved; computed by `t_free.py deps` at blocks 3, 4, 5, `out_t_free_deps.txt`].
- `sum over z in R of z z^T = 0` exactly when every form with Q(0) = 0 has an even number of ones on R.
- So the dependent sets R are the supports of the dual of the shortened RM(2, m), which is RM(m-3, m) punctured at 0.
- RM(m-3, m) has minimum weight 8, reached by the 3-flats. Puncturing at 0 gives weight 7 for the 3-dimensional
  linear subspaces, and nothing lighter.
- Computed: the dependency spaces have dimensions 1, 5 and 16 at blocks 3, 4 and 5. The minimum dependency has 7
  points, and the minimum words are exactly the 1, 15 and 155 three-dimensional subspaces.

### 8.7 The zero-block test

**Statement.** Let Y be a set of live moves, any two of which lie in a common tight support. Then `E[Y, Y] = 0`, and
every form F with Y inside D_F satisfies `|D_F| + r_F >= 2|Y|`.

*Proof.*
- Each tight support D has `E[D, D] = 0` (8.5). So every diagonal entry at a live move is 0, and so is every entry
  between two moves that share a tight support.
- A d x d matrix whose rows and columns indexed by Y meet in a zero block has rank at most `(d - |Y|) + (d - |Y|)`:
  the rows outside Y contribute at most d - |Y|, and the rows in Y are supported on the d - |Y| columns outside Y.
- With 8.5: `|D_F| - r_F = rank E[D_F, D_F] <= 2(|D_F| - |Y|)`. []

This is the case "two rectangles" of w8-rmexact's Proposition RC. `t_small.py` applies it to every form with
`|D_F| + r_F < 2c` (c the clique size), using the largest such Y inside D_F.

### 8.8 When a form is tight (one move per direction)

**Statement.** Take one move in each nonzero direction, W arrays, and the canonical placement (targets T_l of
linear part gamma_l, so that `p_l(Q) = Tr(b gamma_l)` for Q of linear part b). If `W < 3·2^(m-3)`, some form is tight
exactly when some nonzero b has `#{l : Tr(b gamma_l) = 1} = 2^(m-2)`. In particular there is no tight form when
`W < 2^(m-2)`.

*Proof.*
- A form Q is tight when `wt(Q) = |D_Q| = r_Q <= W`.
- The nonzero weights of RM(2, m) below `3·2^(m-3)` are only `2^(m-2)`, reached exactly by the (m-2)-flats. So a
  tight form is a flat with `r_Q = 2^(m-2)`.
- Every nonzero b is the linear part of a flat: `Tr(ux)Tr(vx)` has linear part `sqrt(uv)` (section 2, F2 applied to
  the diagonal terms), and u != v with `uv = b^2` exist. Its odd arrays are those with `Tr(b gamma_l) = 1`. []

Checked at block 5: all 31 linear parts occur on flats, and the criterion agreed with the computed tight-form count in
120 of 120 random target sets (`out_t_tightcode.txt`).

### 8.9 Theorem C3: no saving at block 3 without phases, for any address group

**Statement.** In the no-phase game at block 3 (paid moves dir(z); free shifts, mixing with scratch, and address maps
from any group H of GL(3,2)), every unit for W arrays uses at least 3W moves.

*Proof.*
- *The group.* Let A be the group of 8 x 8 matrices generated by the dir(z) and the shifts X^a. Every element of A is
  a polynomial in the shifts, so A is abelian and commutes with every shift. An address map x -> hx conjugates X^a to
  X^(ha), so it lies in A only if h = 1: A ∩ H = {1}. And `h^-1 dir(z) h = dir(h^-1 z)`, `h^-1 X^a h = X^(h^-1 a)`, so H
  normalises A. The letters generate the semidirect product A ⋊ H.
- *Theorem MT at this level* (w8-rmexact's proof, with A in place of Sym). C[A ⋊ H] is a free right C[A]-module with
  basis H. Left multiplication gives an injective algebra map rho into n x n matrices over C[A] (n = |H|), with
  `rho(a) = diag_h(h^-1 a h)` for a in A and rho(h') a permutation matrix.
- *The new program.* Apply rho to every step of the unit, which is a formal identity in the group algebra (w-bounds
  section 1: right for every number of columns is the same as a formal identity). Each register becomes n registers.
  A move dir(z) becomes n moves dir(h^-1 z), one on each copy. A shift becomes shifts. An address map becomes a
  permutation of registers, which is mixing. So we get a program with no address maps, for nW arrays, with nM moves.
  Array (l, h) is owed `h^-1 K h = dir(h^-1 e_1) dir(h^-1 e_2) dir(h^-1 e_3)`, the kernel in another basis.
- *The certificate.* Take f = 3 columns, the antisymmetric part `Lambda^3(C^8)` (every letter acts the same way on
  each column, so it keeps this part), and N = the span of the 24 face triples (wht38/tiny/close).
  - Shifts and all diagonal matrices keep N (wht38 sentence 3; recomputed).
  - For each direction z, `dim(N + dir(z) N) - dim N <= 8` (wht38 sentences 4-6, by hand; computed: exactly 8).
  - For every basis b, `dim(N + K_b N) - dim N = 24` [computed, `t_face3.py`, ranks modulo two primes. A rank
    modulo p is at most the true rank, so "at least 24" is rigorous, and 24 = dim N is the maximum].
- *Theorem R* (w-bounds section 2) applied to the new program: `8 nM >= 24 nW`. []

**What it covers and what it does not.**
- Every address group and every placement at block 3, hence w8-rmexact's Corollary L3, the canonical Singer case, and
  every small-lift survivor at block 3.
- With phases it survives as long as the moves, shifts and phases generate a group that contains no address map
  except 1 and is normalised by H: for example Pauli Z phases and scalars (they turn dir(z) into a scalar times its
  inverse). The engine's phase letters break it. They are `phase_mat(z, c) = diag(i^[z·x + c])` (regmachine.py,
  Letters.lean:42), so `S dir(e_1) S` is e^(i pi/4) times a Hadamard on bit 1, `S_1 S_2 S_(e1+e2)^-1` is CZ, and
  `H_2 CZ H_2` is CNOT, an address map. Then rho would produce address maps inside the new program, and those leak
  against N. So the engine game at block 3 is not covered.
- It says nothing at block 4 or above. There the transfer needs a certificate of the address-map-free game whose
  demand is the same for every basis; none with demand above 3.058 per array is known even for the standard basis.

**The same proof without MT, and why it does not contradict Theorem F** [proved].
- Put Y = ⊕ over g in H of gN inside ⊕ over g of V. Let A-elements act on every summand alike, and let an address map
  h act by `(h v)_g = h v_(h^-1 g)`. This is a representation of A ⋊ H (it is V tensored with the regular
  representation of H, pulled back along A ⋊ H -> H). Y is kept by shifts and by every address map.
- Its cap on dir(w) is `sum_g cap_N(g^-1 w)` and its demand is `sum_g dem_N(K in the basis g^-1 e)`. Take H = all
  of GL(m,2): every unit with a smaller address group is also a unit with GL. GL is transitive on directions and on
  ordered bases, so Theorem R gives: moves per array >= (average demand over all bases) / (average cap over all
  directions). So any certificate for programs without address maps transfers to every address group, with these
  two averages in place of the standard basis and the worst direction. At block 3 every basis has demand 24 and
  every direction cap 8, so nothing is lost.
- w-bounds' Theorem F says no perfect certificate exists at block 3 when address maps AND phases are free. There is
  no contradiction: with the engine's phases, the moves, shifts and phases generate CNOT, so there is no
  homomorphism from the game's group onto H that kills them, and Y is not a representation of the engine game. The
  block-3 face certificate thus shows that an engine saving at block 3 must use both kinds of letters: address maps
  (wht38: with phases alone the face certificate gives 3W) and phases (Theorem C3: with address maps alone, 3W
  again). Each kind alone is harmless; only their interaction (which generates CNOT-like address maps from moves)
  can save.

