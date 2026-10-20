<!-- header: frame-packing-two-arrays -->
# Closed cancellation patterns, the parity law and Conjecture TC (the orbit track of packing round 2)

**About this file.** This is the working proof file as an AI agent wrote it during the round of 15 October 2026, copied with private file paths removed and nothing else changed. "I" in it is that agent. Its opening lines say that no skeptic had read it; that was true when it was written. What happened afterwards is below, and where a later skeptic's verdict differs from this file, the later verdict wins. No person has read this proof and nothing in it is checked in Lean. File names inside it refer to the programs and outputs now in `../programs/` and `../data/`; [README.md](README.md) maps the working names.

**Status now.** Upheld by one AI skeptic, with six corrections that change no result: Remark 7's example is described wrongly (the point survives one repair earlier); "none delivered even one output exactly" says more than was decided (519 of 4,584 target outputs were refused by no screen); the controls did not show "the same mix of verdicts"; two figures have no output file behind them (the skeptic's reruns give a cycle-orbit rate of 8.5 to 16 % at block 4, not 20 %); the proof that orbits multiply junk covers only copies inside one array; bookkeeping slips. Conjecture TC is open.

---

# PROOF.md (pk2-orbits): closed cancellation patterns in the frame game

pk2-orbits, packing round 2, 15 October 2026. No skeptic has read this. Labels: [proved] (proof written here),
[computed] (program and output in this folder), [argued], [read]. Setting and notation as in `../fp-build/PROOF.md`:
the frame game at block m, `A = C[Sp(2m,2)]`, free algebra `C[P]` (P = stabiliser of the raw frame L0), paid letters
`T_u` (u not in L0), target `K`, L1 = K^-1 L0 the target frame. Contents are W-tuples (one entry per input); a
*walker* is one term `(l, g)`; its frame is `pos(g) = g^-1 L0` and its level is `d(L0, pos g)`. Left multiplication by
P never changes a frame (fp-build Lemma 3.1). A *state* is an input or the output of a paid move; the readout of
output l is a C[P]-combination of states equal to `K e_l` (module game, p6-prove Lemma M).

## 1. The closed pattern that is known to work: the order-3 control, generalised

**Lemma 1 (two-frame closed patterns: the parity law)** [proved; computed]. Let sigma lie in the two-frame
stabiliser `H2 = P ∩ K^-1 P K` (≅ GL(m,2): the elements `[[A, 0], [A + A^-T, A^-T]]`), and let
`sigma' = K sigma K^-1` (in P). The program

    r1 = y1 + y2,   r2 = y1 - sigma y2,   then K (m moves) on each register

(2 in 2m) delivers both outputs **if and only if the order n of sigma is odd**.

*Proof.* After the moves, `r1 = K y1 + K y2` and `r2 = K y1 - sigma' K y2` (because `K sigma = sigma' K`), and these
are the only states with walkers at level m; their walkers are all at level m. Left multiplication by P keeps levels,
so the readout's level-m part must already give `K e_l`: with a, b in C[P] and the right factor K cancelled (left
multiplication by K is injective on C[Sp]), output 1 needs `a + b = 1`, `a = b sigma'`, i.e. `b (1 + sigma') = 1`;
output 2 needs `b (1 + sigma') = -1`. In a finite-dimensional algebra a one-sided inverse is an inverse, so both
outputs exist iff `1 + sigma'` is invertible in C[P], and then `a r1 + b r2` is exactly `K e_l` (no lower walkers).
If n is odd, `(1 + s)(1 - s + s^2 - ... + s^(n-1)) = 1 + s^n = 2` for s = sigma', so `1 + sigma'` is invertible.
If n is even, `(1 + s)(1 - s + ... - s^(n-1)) = 1 - s^n = 0`, so `1 + sigma'` is a zero divisor and not invertible. []

[computed, `cal.py`, `out_cal.txt`] 41 random sigma in H2 at blocks 3, 4, 5, orders 2 to 31: the exact end solve
says YES exactly for the odd orders (3, 5, 7, 15, 21, 31; readouts of 2n terms) and NO for every even order (2, 4,
6, 8, 12, 14). The 3-cycle address map of o3-build is the case n = 3; cut by one move, it fails.

So the "closed pattern" of the order-3 control is a group-algebra inversion over a cyclic subgroup of P. It works
because both registers end with all their walkers on L1 (it is *top-closed*, below), and it costs 2m: each of its
two coefficient rows `(1, 1)`, `(1, -sigma')` is paid for by a full K.

## 2. Orbit injections

An *orbit injection* (orb.py) adds, at one site (register r, before its k-th move),

    sum_{i<n} c_i sigma^i mu0 S_src,     sigma in P of order n,

where mu0 is the landing chart that puts a walker g of S_src on the frame X of a chosen junk walker at a chosen
later output of r (fp-build Lemma 3.2). Let w be the letters of r from slot k to that output and `Y = w^-1 L0`.
Copy i of the walker ends at `w sigma^i mu0 g`, on the frame `g^-1 mu0^-1 sigma^-i Y`.

**Lemma 2 (stabiliser orbits are group-algebra coefficients at the landing output)** [proved]. If sigma fixes Y
("stab" orbits), every copy lands on X, and at that output the orbit contributes

    (sum_i c_i tau^i) . (w mu0 g),     tau = w sigma w^-1 in P,

i.e. exactly one landing with a coefficient in C[<tau>] instead of a scalar. At an earlier output of r (letters w' from slot
k, `Y' = w'^-1 L0`) copy i sits on the frame `g^-1 mu0^-1 sigma^-i Y'`; these are different frames for different i
unless the powers of sigma also fix Y'. So the copies spread over up to n frames between the injection and the landing,
and merge in frame (not in element) at the landing output.

*Proof.* `w sigma^i = (w sigma w^-1)^i w` and `w sigma w^-1` fixes `w Y = L0`, so it lies in P. []

**Lemma 3 (copies never combine)** [proved]. Two copies `sigma^i mu0 g` and `sigma^j mu0 g` (i ≠ j) of an orbit
injection are different group elements, and they stay different at every later output of their register (every
step multiplies both on the left by the same element). So copies never combine with or cancel each other: an orbit of
length n adds exactly n walkers to every later state of its register (unless a copy meets a walker that was already
there). Only the landing copy, or copies that land on other junk, do useful work; the rest is new junk.

*Consequence* [computed, `growth.py`]: the orbit repair multiplies the unpartnered junk faster than one-at-a-time
landing (section 4 of NOTE.md).

## 3. Why the order-3 trick does not share a move: one joined pair gives rank one

**Lemma 4 (one joined pair, rank one)** [proved]. Let `v = lambda0 e1 + mu0 e2` be the content of a register
right after a join (one walker per input; `rho = lambda0^-1 mu0`). Every state descending from v alone (any moves and
free factors, and any C[P]-combinations of such states) is `Phi v` with Phi in C[Sp]. If such a state is
*top-closed* (all its walkers on L1), then `Phi lambda0 = beta K` with beta in C[P], `s := K rho K^-1` lies in P, and

    Phi v = beta (K e1 + s K e2).

So the top-closed states of one joined pair have left rank one over C[P]: a readout made of them can give neither
`K e1` nor `K e2` (from `sum gamma_r beta_r = 1` one gets `sum gamma_r beta_r s = s ≠ 0`, and symmetrically). Both
outputs therefore need either a second join, or a single-input state on L1 (m moves of that input alone).

*Proof.* The first statement is linearity. `Phi lambda0` is a sum of walkers on L1, each of the form pK, so it is
`beta K`. If beta ≠ 0 then `Phi mu0 = Phi lambda0 rho = beta K rho` has its walkers on `rho^-1 L1`, which must be L1,
so rho fixes L1 and s fixes L0, i.e. s is in P; then `beta K rho = beta s K`. []

This is where every attempt to "share a move" in the order-3 pattern breaks. The shared move must serve both
registers, so both descend from one joined pair, and their top-closed states are proportional. A second coefficient
row needs a second join, whose new walker must still climb all the way to L1.

**Lemma 5 (disjoint histories)** [proved; immediate]. The *history* of a walker is the set of paid moves applied to
it and its ancestors. A walker on L1 has a history of at least m moves. So if a readout uses, with nonzero weight, a
walker of input 1 on L1 and a walker of input 2 on L1 whose histories are disjoint, then M ≥ 2m.

**Proposition 5' (where a second join puts its walker)** [argued; proved for the shapes listed]. In a top-closed,
cancellation-free two-input design whose joins are one joined pair plus injections of single-input prefix states,
Lemma 4 forces some readout state to contain an injected prefix walker nu on L1. In each of these shapes there is then
a readout walker of the other input whose history is disjoint from nu's, so Lemma 5 gives M ≥ 2m:
- nu injected into a new register (its own lineage): the pair's walker of the other input avoids that lineage;
- nu injected into the pair's register at time tau, followed by two branches A, B: nu's history is (its prefix) +
  (its branch, say B); the walker of the other input in branch A has history (its prefix) + (the pair's moves before
  tau) + (branch A), disjoint from nu's.
  The count: `M ≥ (a1 + f1' + f_A) + (a2'' + f_B) ≥ 2m`.

*Not covered:* readouts that use states with walkers below L1 (junk that must cancel), and designs with cancellation
inside registers. **Conjecture TC** [argued]: every top-closed two-input program needs 2m moves. If it holds, every
two-array saver must cancel junk below the top level, and closed patterns of the order-3 kind (group-algebra inversion
on L1) can never share a move.

**Lemma 6 (no exact landing of prefix walkers on branch junk)** [proved]. In the trunk family (prefixes a_1, a_2, one
trunk move that raises both walkers, branches), no walker of a single-input prefix state can be injected into branch
A so that it lands *exactly* (same group element, not only the same frame) on A's junk walker J_A; likewise for B.

*Proof.* Let J_A's frames along branch A be Z_0 = P_{a2+1}, Z_1, ..., where P_j are the frames of input 2's prefix
chain. By fp-build Lemma 3.3 the path P_{a2}, P_{a2+1} = Z_0, Z_1, ... is a geodesic. An exact landing of the prefix
walker of level j at slot k needs its element to equal J_A's element at slot k, so `P_j = Z_k`. Then
`d(P_{a2}, P_j) = k + 1` along the geodesic, while along the chain `d(P_{a2}, P_j) = a2 - j`, so `a2 - j = k + 1`;
but `d(Z_0, P_j) = k` and `d(P_{a2+1}, P_j) = a2 + 1 - j = k + 2`. Contradiction. []

[computed, `t_lemma6.py`, `out_lemma6.txt`: 1,000 (slot, prefix walker, junk) tests on 120 random trunk skeletons
at blocks 3, 4, 5: no exact landing; 60 frame-level landings.]

So a landing of a prefix walker always leaves a free factor q ≠ 1 between the landed walker and the junk (it lands on
the frame, never on the element); fp-build's repairs and all orbit repairs are of this kind.

## 4. Screens

**Remark 7 (the monomial-split relaxations ignore rigidity)** [proved; computed]. fp-build's linearised psi screen and
its chart-level analogue used here (`chlin`: end_solve over the monomial components of every state) are valid
necessary conditions: if `sum beta_k S_k(c) = K e_l` then `K e_l = sum_k sum_mono (beta_k c^mono) S_k^mono`. But the
monomial components of a state are not states: splitting separates an injected copy from the walkers of the register
it was injected into, as if it rode in a register of its own. Example (wave A, `two` m = 4, seed 103, trial 1): prefixes
a = (3, 3), no branch moves, two-frame orbits injected into the trunk and into input 2's prefix register. The trunk
output is the only state with walkers at level 4, and *at least W states must reach level m* (the level-m part of the
readout lies in the C[P]-span of the level-m parts of the states, which has dimension at most |P| per state, while it
must contain `C[P] K e1 + C[P] K e2`, of dimension 2|P|). So at most one output can be delivered; yet both relaxations
pass for both outputs. So a pass of psilin or chlin says little for orbit designs; a failure is exact for every choice of
the orbit scalars.

**Lemma 8 (layered screens)** [proved]. For d = m, m-1, ..., if `sum_r beta_r S_r = K e_l` with beta_r in C[P],
then also `sum_r beta_r S_r^{≥d} = K e_l`, where `S^{≥d}` keeps the walkers of level at least d (left
multiplication by P keeps levels, and K e_l has level m). So "the layer-d system is inconsistent" is an exact NO for
the scalars used. orb.py runs it from the top down whenever the full solve exceeds the cap ("L3n" = dies at layer 3).
This is the frame-game form of wave 7's layered test, without the numerical step: every layer is solved exactly.

## 5. Calibration of the pipeline [computed, `out_cal.txt`]

| case | block | expected | psi | exact | heightfilter |
|---|---|---|---|---|---|
| plain 2 in 2m | 3, 4, 5 | yes | pass | yes (1-term readouts) | height system exists (cannot save), bound 2m |
| plain cut by one move | 3, 4, 5 | output 2 fails | 10 | yes, no | — |
| order-3 control (3-cycle address map) | 3, 4, 5 | yes | pass | yes (6-term readouts) | height system exists |
| order-3 control cut | 3, 4, 5 | no | fails | no, no | — |
| two-frame sigma, odd order | 3, 4, 5 | yes | pass | yes, 2n-term readouts (22 cases) | height system exists |
| two-frame sigma, even order | 3, 4, 5 | no (Lemma 1) | pass | no (19 cases) | height system exists |
| orbit designs at 2 in 5, block 3 (14W/5) | 3 | never | | never (section 4 of NOTE) | |
| orbit designs at 2 in 6, block 4 (41W/13) | 4 | never | | never | |

The even-order rows show the psi screen passing a pattern that is dead at chart level: psi forgets the left P
factor sigma', which is exactly what decides invertibility. So psi is necessary, not sufficient, with an exact example
(fp-build's skeptic noted that no example had been given).
