<!-- header: frame-packing-two-arrays -->
# The skeptic's step-by-step check of the build track

**About this file.** This is the working proof file as an AI agent wrote it during the round of 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent. Its opening lines say that no skeptic had read it; that was true when it was written. What happened afterwards is below, and where a later skeptic's verdict differs from this file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean. File names inside it refer to the programs and outputs now in `../programs/` and `../data/`; [README.md](README.md) maps the working names.

**Status now.** This is the skeptic's own record. Its verdict: upheld with corrections, one sub-claim (the closest design) refuted.

---

# v-fp-build: step-by-step check of fp-build (15 Oct 2026)

Skeptic's record. Labels: [checked] (I re-derived the step), [computed] (my own script in this folder; fp-build's code is
used only where I say so, and then only to regenerate its designs), [read] (from a named file), [argued].

## 1. Blocks and the faithful block

- **Lemma 1.2 (blocks)** [checked]. `s = i.1` is a free letter (phase letter with z = 0, c = 1), central, and
  `e_k = (1/4) sum_j i^(-jk)[s^j]` evaluates to 1 at f = k (mod 4) and 0 otherwise. The blocks are independent, and
  coefficients may differ per block with the same paid moves. The delivered scalars may differ per block too. They are
  merged by the free element `sum_k u_k^-1 e_k`. Correct.
- **Prop 1.3** [checked; computed at block 2, `vgroup.py`, `out_vgroup.txt`, own 4 x 4 matrices].
  - |G| = 46,080, and the scalars in G are exactly {1, i, -1, -i}, so |G/mu4| = 11,520 = 16 x 720.
  - There are 720 symplectic images, each fibre has 16 elements modulo mu4, and every fibre's matrices are linearly
    independent.
  - The image of B is exactly P (48 elements), and the preimage of P is exactly B (3,072 elements).
  - So Phi(g) = g (x) gbar is well defined on C[G]_mu4, multiplicative, injective fibre by fibre, and equal in
    dimension to its target. It is an isomorphism, and C[B] goes onto M_n (x) C[P].
  - At block 3 the group order 256 x |Sp(6,2)| [read: e3-rmbuild Lemma N] gives the same count. Correct.
- **Remark 1.4** [read: `ledger.py` lines 15, 102, 149]. ledger.py normalises group elements modulo {1, i, -1, -i} and
  moves the scalar into the coefficient. That is the quotient [ig] = i[g], i.e. block 1, so its identities certify
  f = 1 (mod 4) only. Correct; it agrees with p6-prove's warning in the brief.
- **Remark 1.5 (even-f blocks)** [argued; dimension counts checked by hand for every m].
  - Block 0: |Sp| + (4^m - 1)^2 |Sp_v| = 4^m |Sp|.
  - Block 2: N+^2 |O+| + N-^2 |O-| = (N+ + N-) |Sp| = 4^m |Sp|, since |Sp|/|O±| = N±.
  - Both add up. The cocycle question is open, as fp-build says.
  - **Omission: block 3 (f = 3 mod 4) is automatic.** g -> conj(g) (x) gbar is an isomorphism of block 3 onto
    M_n (x) C[Sp], by the same fibre argument. Its lift Lambda_3 = Phi_3^-1(1 (x) .) sends T_z to (free) x [dir z] and
    C[P] into C[B]. [computed: Lambda_3 is an exact homomorphism at f = 3 and 7 (errors 3e-16, 4e-16) and wrong at
    f = 1, 2, 4, 5, 6, 9.]
  - With Lemma 1.2 (coefficients Lambda_1(.) e_1 + Lambda_3(.) e_3), **every frame-game program gives an engine word,
    with the same moves, that is right at every odd f.**
  - So PROOF §2's "f = 3 mod 4 is not automatic either" is true only of the f = 1 lift. The remainder after a
    frame-game find is the even-f blocks alone, as the claim's "meaning" item 2 says. Certificates proved at
    f = 3 (mod 4) also hold in the frame game.

## 2. The sandwich

- **Left inclusion** [checked]: e3-rmbuild Lemma N, read through e_N. Correct.
- **Right inclusion** [checked; computed].
  - Lambda_1 = Phi^-1(1 (x) .) is an algebra homomorphism C[Sp] -> C[G]_mu4.
  - [computed] Lambda_1(T_z) = ((1-i)/2)[dir z] + ((1+i)/2)[X^z dir z] = c_z [dir z]. That is Lemma R's corrected
    move, exactly.
  - [computed] Lambda_1(p) is supported on free letters for all 48 p in P.
  - [computed] Lambda_1 is an exact homomorphism on raw random complex arrays at f = 1, 5 and 9 (errors at most
    4.4e-16) and wrong at f = 2, 3, 4, 6 and 7 (relative errors 1.4 to 2.0).
  - In the module form (L_t = L_{t-1} + C[B] dir(z_t) v_t), the lifted word has the same W and the same paid moves,
    and it is right at every f = 1 (mod 4). Correct.
  - **Precision correction.** The lifted word's free coefficients are multi-term elements of C[B]: c_z after every
    move, Lambda(p) (16 terms at block 2), and Pexp(K). In the engine's step format (scalar mixes, one free letter per
    step; RULES §0 and sl-comb §0) these need scratch registers, or repeated-source mixes, which the brief's notice
    calls scratch words.
    - So "same skeleton" holds for W, M and the move sequence, not for the register count.
    - This changes nothing used: units allow scratch, and Theorem R certificates hold step by step.
    - It does matter for any reading of Theorem 2.2 next to "two arrays, no scratch" results (Corollary 2').
- **Corollary 2.3** [checked; read].
  - The block-4 hook is V (x) Lambda^8 V inside V^(x)9: bm.py's header defines the potential in "the hook module
    V (x) Lambda^k V".
  - MAP A17 independently names "the f = 9 block-4 hook".
  - 9 = 1 (mod 4), and the block-3 hook is at f = 5. So 14W/5 and 41W/13 hold in the frame game.
  - Theorem R's monotonicity holds for each elementary step, scratch included, so it covers the lifted words.
  - 41 x 2/13 = 6.31 > 6, so 2 in 6 is excluded at block 4. 14 x 2/5 = 5.6 > 5 and 14 x 3/5 = 8.4 > 7, so 2 in 5 and
    3 in 7 are excluded at block 3.
  - No stated ceiling falls below a known unit. The round adds no new ceiling, and the frame-game floors are the
    engine's.
- **Remark 2.4 (lanes)** [checked].
  - Rows of M_n(C[Sp]) are mixed by free left multiplication, and columns are kept. So a W-in-M engine word gives an
    nW-in-nM frame program, after one final free K^-1 (x) 1 on each delivering register.
  - Lambda embeds frame programs back.
  - The ratio is the same both ways. Correct.

## 3. Packing lemmas and the screens

- **Lemma 3.1** [checked]: left multiplication by one element per step.
- **Lemma 3.2** [checked; computed, `vland.py`, own normal form: a Levi element straightens the x-projection, then a
  symmetric phase].
  - 200 of 200 random cases are correct at each of blocks 3, 4, 5 and 6.
  - The necessity direction d(pos g, pos(w mu g)) = level(w) held on every sampled mu.
- **Lemma 3.3** [checked].
  - Branch A's y2 junk sits at X with d(P_a2, X) = d(P_a2+1, X) + 1, so it is not on y2's prefix chain.
  - The only branch-B state with a y2 walker at X is the one after the same number of moves. Both branches inherit
    the trunk's relative element, so that state's y1 walker sits at A's y1 position. The two states differ by one
    element of P, with the same coefficient ratio, so cancelling y2 cancels y1.
  - Pre-trunk states hold one input. Correct.
- **Prop 3.4 (psi)** [checked]. psi(beta x) = eps(beta) psi(x), and every content is a C[P]-combination of inputs and
  move outputs. The linearised version is a valid relaxation: span{v_(s,mono)} contains every span{psi(x_s)(c)}.
  - [computed, `vpsilin.py`: my own symbolic build and exact elimination over Q, on designs regenerated with
    design.py] It agrees bit for bit with fp-build's psi_lin_screen on 75 of 75 designs: 30 at block 4, 30 in dup
    mode at block 4, 15 at block 5.
- **Lemma 3.5** [checked].

## 4. Measured quantities

- **Tallies** [computed, `vtally.py`, own parser of the raw Colab lines]. Every production row reproduces exactly:
  - 5,271 designs, 0 with both outputs, 0 CANDIDATE lines;
  - 3,729 (0.707) dead for every scalar choice, 1,187 (0.225) open, 167 with one output exact. That is 166 two-array
    designs plus one 3 in 10.
  - The controls agree too: 3 in 7, 300 designs, 249 dead; 2 in 3, 200 designs, 70 dead.
  - 8,522 = 5,271 + 3,251, and 1,130 = 330 + 300 + 300 + 200.
- **Null repairs** [computed, `vselfinj.py`]. I reran design.py unchanged for the first 30 designs of production
  seeds 61 and 62 (dup), and classified each injection with my own arithmetic:

  | seed | injections | genuine (new walkers) | pure rescalings (the register's own current content, identity chart) | reweightings (subset of current walkers) |
  |---|---|---|---|---|
  | 61 | 174 | 99 | 50 | 25 |
  | 62 dup | 163 | 116 | 24 | 23 |

  - Genuine landings per design were at most 6 (seed 61) and 7 (seed 62).
  - A rescaling leaves the C[P]-span unchanged. If its drawn coefficient is -1 (probability 1/15), it deletes the
    register.
  - So "up to 12 landings" overstates the depth. The negative result is still a negative result for the family, but
    a weaker one.

## 5. The closest design is not a psi pass (refutes the "closest" field)

[computed, `loadhit.py` (restricted unpickler: numpy reconstruction globals only), `vfg.py`, `vhit.py`,
`out_vhit.txt`]

- **Structure.** Of its 12 injections, 11 inject the register's own current content with the identity chart, which is
  a pure rescaling.
  - The twelfth, A slot 0, lands input y2 with a nontrivial chart. It is the only genuine landing.
  - The trunk injection ('Rt', 0, 0) adds c0 y1 to a trunk that already holds y1 (a1 = 0). The trunk's y1 weight
    becomes alpha = 1 + c0.
- **fp-build's scalars.**
  - c0 = -1.00000004 - 3.5e-8 i, so alpha = -3.7e-8 - 3.5e-8 i.
  - At A slot 0, 1 + c5 + c6 + c11 = 2.3203 - 3.2971 i, which is minus the landing scalar c1 = -2.3203 + 3.2971 i.
- **Exact obstruction.**
  - The psi coordinates (y1, X1), with X1 the trunk's y1 position at level 1, and (y2, L1) are touched by exactly one
    stored state, the trunk output, with coefficients (alpha, 1). This was checked with exact support at random
    rational scalars; support can only shrink at special scalars.
  - So output 2 passes psi iff alpha = 0. Its least-squares residual is |alpha|/sqrt(1 + |alpha|^2) exactly. I
    reproduced 0.7071, 0.0995, 1e-2, ..., 1e-8 as alpha -> 0, also with random other scalars.
  - At alpha = 0, y1 enters no state except the raw input: a1 = 0, and no injection carries y1. So output 1 fails.
  - **The design is psi-dead for both outputs at every scalar choice.**
- **What the optimiser found.** The 5.2e-8 residual is the infimum approached as alpha -> 0, with readout coefficients
  of about 1/|alpha| (2.7e6 at |alpha| = 1e-8). It is not a pass.
- **The continuations** (out_hit2_*). They start from these scalars and sit at the same residual, 5.13e-8 to 5.22e-8 in
  most steps. That is the same degenerate valley.
- **Consequences.**
  - "Passes the psi screen for both outputs with special scalars, all nonzero" is false.
  - "12 landings" is one landing, one trunk reweighting and ten rescalings.
  - "It is the only design in the round that got past positions" and "it shows that the psi screen is necessary but
    not sufficient" are unsupported.
  - Corrected: none of the 341 bilinear samples passes psi; the smallest genuine residual is 0.21.
