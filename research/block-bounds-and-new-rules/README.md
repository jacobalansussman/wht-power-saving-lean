# Block bounds in Lean, Toffoli and column-mixing maps (13-14 October 2026)

On 13 and 14 October 2026 I asked two questions about the small units that would give a large saving. First, how
many costly steps does a unit of each block size need at the least, and can those lower bounds be checked by a
machine? Second, does OpenAI's cost model allow free steps that the usual game leaves out (non-linear relabellings
such as Toffoli, maps that mix neighbouring digit columns, address rows the unit leaves alone), and do they open a
door? This folder tells what I tried, what held, what broke and what is still open, with the programs, certificates
and Lean files that produced each result.

## In short

- **No unit was found**, in any game and at any block size. What this part of the work produced is lower bounds,
  limits of the method behind them, and a map of which new free steps are allowed and what they buy.
- **Proved in Lean** (in the engine game of this repository, with only Lean's three standard axioms): the counting
  rule (Theorem R) for the engine's own programs; block 2 never saves (every unit needs 2 paid moves per array);
  block 3 needs at least 14/5 = 2.8 paid moves per array, so a block-3 unit saves at most 6.28e-2. The block-3 proof
  has no hypothesis left: its one numeric fact is checked by Lean's kernel.
- **Checked by computer, not in Lean:** block 4 needs at least 41/13, about 3.154 moves per array, so a block-4 unit
  saves at most 0.171 (re-run for this revision); certificates for the widened game at blocks 5 and 6. **Argued by
  hand from exact computations, with two AI skeptics:** at least 2.99172 moves per array at block 5 and 2.97290 at
  block 6 in the engine game.
- **The method has a ceiling.** No hook, cut or coordinate certificate, the kinds behind the bounds above, can reach
  4 moves per array at any block size, so from block 5 on these tools can never rule out a saving of 0.1 (argued by
  hand, not checked; upheld by AI skeptics; section 2). Certificates that only the engine game keeps, such as lifts,
  are not covered by this ceiling.
- **Toffoli and every other relabelling are allowed** in OpenAI's model if the engine gets an AND table next to its
  XOR table (argued by hand, not checked), **but they buy nothing at blocks 2 and 3**: the bounds above hold with them
  free. At block 4 one array gains many more routes, but no sharing gain was found.
- **Maps that mix neighbouring digit columns** (chain maps, `row ^= row << 1`) are allowed too (checked by computer)
  and break every certificate of the usual kind. The lead was mostly closed again by new proofs: block 2 needs 2 moves
  per array in every game whose free steps touch boundedly many entries per row, and block 3 keeps 2.8 for the chain
  maps that have a reading at 1 (both argued by hand, not checked; upheld by AI skeptics). **Still open:** 3 arrays in
  8 moves at block 3 with the running XOR free (saving 0.107), where the best bound is 1.763 moves per array and more
  than 2.667 is needed.
- **Spare rows** do not help the mechanism of OpenAI's own network on rank-2 and rank-3 targets (argued by hand, not
  checked; one AI skeptic); one rank-3 case is open by a margin of 0.012 moves per array. **OpenAI's own finite-win
  game** is closed on 4 and 8 entries (argued by hand from a classical theorem, no skeptic) and open on 16 entries
  between 28 and 31 calls; for the subset-sum kernel on 8 entries no word with 11 calls exists (checked
  by computer, and re-run for this revision).
- No person has read any of the proofs in this folder. Most results were checked by one or two independent AI
  skeptics; each status below says by whom or by what.

## The question

A **unit** is a small, fixed program that transforms W arrays at once with fewer costly steps than doing them one
at a time. At **block size m** each array has 2^m entries and must be transformed along m address directions; one
**paid move** applies one direction's 2-point transform to one array (inside the recursion it becomes one recursive
call), and everything else is a **free step** that costs a fixed amount of work per entry: relabelling addresses,
multiplying entries by phases, adding multiples of one array to another. The plain schedule needs W·m paid moves. A
unit "W in M" with M < W·m gives, used recursively, a Walsh-Hadamard transform with about N (log N)^(1-s)
operations, where the **saving** is s = ln(mW/M)/ln m; the earlier revisions explain how a saving carries over to
the Fourier transform and to convolution. A single small unit would be enough: 3 arrays in 8 moves at block 3 would
save 0.107, about 60 times the records of this repository, which come from very large units.

Two kinds of work can settle such a question. A construction finds a unit. A **lower bound** proves that every unit
of some shape needs at least so many moves, which excludes whole families of targets at once. This section is about
lower bounds, and about whether the game everybody analyses is really the game OpenAI's cost model allows. The
standing bar for a new route is a saving of 0.1; since 15 October a saving of 0.01 counts as a stepping-stone tier.

## Words used here

- **Games.** They differ in what is free. The **engine game** has the free steps of the Lean engine of this
  repository: affine (XOR-linear) maps of the address bits, shifts, the phase patterns i^(quadratic form) that go
  with them, and rational mixing of arrays. The **widened game** also makes every permutation of the addresses
  (Toffoli included) and every pattern of phases 1, i, -1, -i free. A bound that holds in the widened game also holds
  in the engine game. The **commutative slice** has only shifts, paid moves and mixing.
- **Digit columns, f; rows.** Inside the recursion an address is cut into f digit columns of m bits. The recursion
  uses the unit at every size, so a unit must be right for every large f. The same address read the other way is m
  rows of f bits.
- **Certificate, Theorem R, cap, demand.** A certificate is a subspace of test vectors that every free step keeps.
  Its **cap** is the most one paid move can enlarge it, and its **demand** is how much one finished array needs.
  Theorem R says every unit needs at least W × demand / cap paid moves. The ratio demand/cap is the certificate's
  "moves per array".
- **Hook.** My best family of certificates: the span of the tensors e_x ⊗ e_T with x in T inside V ⊗ Λ^k V, with
  V = C^(2^m). At block 3 (k = 4, f = 5) its cap is 80, two moves in a row raise it by at most 144, and an array
  needs 224, so 14/5. At block 4 (k = 8, f = 9) the numbers are 27,456, 48,576 and 86,592, so 41/13.
- **Frame.** Where a piece of work stands up to free steps. One paid move takes it to a neighbouring frame; the raw
  input and the finished transform sit m steps apart. At block 3 there are 135 frames.
- **Tight.** A unit is tight on a certificate when it uses exactly the number of moves the certificate allows, so
  every move must raise the count by its full cap.
- **Toffoli.** The relabelling (a, b, c) → (a, b, c XOR (a AND b)); on 8 addresses it swaps 110 and 111. It is the
  simplest relabelling that is not linear. **Controlled-S** multiplies by i where two chosen bits are both 1.
- **Chain map.** A free step that XORs each digit column's bits into the next column (`row ^= row << 1`). Its
  inverse is the **running XOR** along each row. `row0 ^= row1 << 1` adds row 1, shifted by one column, to row 0.
- **Reading at a value.** Write each row as a polynomial in a shift S; a chain map multiplies the rows by a matrix
  G(S) of polynomials. Putting a value λ in place of S gives a small matrix that describes what the map does far from
  both ends of the address. A map has a reading at λ when that matrix is invertible. λ = 1 lives in the two-element
  field; alpha, with alpha² = alpha + 1, in the four-element field.
- **Pauli, Pauli-blind.** A sign flip or a bit flip of one address bit; a column-wide Pauli does the same flip in
  every column. A certificate is Pauli-blind if it treats every column-wide Pauli as doing nothing.
- **Spare rows.** Address bits that a unit leaves untouched; the target is then a **rank-k target** that applies
  only the first k directions.
- **OpenAI's finite-win game.** The game of OpenAI's finite-tensor paper: one array of 2^b entries; each **call**
  applies the 2 × 2 kernel to one pair of entries; shuffling and rescaling entries is free. OpenAI's own theorems turn
  a word with fewer calls than the plain b·2^(b-1) into a saving for the Fourier transform.
- **Status words.** *Proved in Lean*: Lean checked every step from the engine's own definitions, with only the
  three standard axioms (propext, Classical.choice, Quot.sound). *Checked by computer*: an exact computation, or a
  stored certificate checked by a small separate checker. *Argued by hand, not checked*: a written proof; where it
  says so, an independent AI skeptic re-derived it and found no error, but no person and no proof checker has
  checked it. *Refuted*: shown false, and by what. *Open*: nobody has settled it either way.

## 1. The lower bounds, block by block

### 1.1 The counting rule, in Lean

**Theorem R** takes any subspace that every free step keeps and turns it into a lower bound: summed over the paid
moves, the growth of the subspace is at least what the targets need. It holds even when the mixing is singular and
when scratch arrays are left holding garbage at the end. It is a general form of the counting lemma of a pull request
in the community's repository (see Credits). Lean proves it on the engine's own definitions, for any kept subspace.
Lean also proves that every free step of the engine keeps the hook, and that one paid move raises the hook by at most
(2^m/2)·C(2^m - 2, k - 1) at every block size (that this cap is attained is checked by computer at blocks 2 to 6).
**Proved in Lean.** Two AI skeptics rebuilt the four files independently, and five deliberately broken versions failed
to build, as they should.

### 1.2 Block 2: no saving

Every unit at block 2 needs at least 2W paid moves, which is the plain cost.
- In the engine game: **proved in Lean** (`block2_bound`; demand 8, cap 4 at f = 3).
- With every monomial map free: a certificate of ratio 2 (dimension 24, cap 8, demand 16), **checked by computer**
  and replayed by a second agent with its own code. Separately, 3 arrays in 5 moves, and every 3k in 5k, is
  impossible for any complex phases: an exact certificate shows that the first move can add nothing inside the
  target. That is **argued by hand, not checked**, with an exact computation; no skeptic read it. So a block-2 unit
  that could be lifted to block 4 (it would have been worth 0.1315) does not exist.
- In every game whose free steps touch a bounded number of entries per row, which includes Toffoli, controlled-S,
  chain maps and the running XOR: Theorem 2W of 15 October, **argued by hand, not checked; one AI skeptic re-derived
  it**. It is told in `../moonshot-small-blocks/`. It supersedes the partial block-2 results for chain maps in
  section 4 below.

### 1.3 Block 3: at least 2.8 moves per array

Every unit at block 3 needs at least 14W/5 paid moves.
- **In the engine game: proved in Lean, with no hypothesis left** (`block3_bound`). The bound needs one number,
  the hook's demand of at least 224, which is the rank of a 224 × 224 system. Lean's kernel is far too slow for that
  in one piece, so the proof splits it into 8 independent blocks of 28 × 28, using the transform's own symmetries
  (characters), and Lean proves that the split is valid. Each block's entries come from one closed formula proved in
  Lean, and each block is shown invertible modulo the prime 65537 by an inverse that was computed outside Lean and
  that the kernel checks entry by entry (`decide +kernel`, no `native_decide`). Changing one entry of an inverse by 1
  makes Lean reject the check (a negative control). The file `Dem3.lean` has 432 lines; the whole development, 1,662.
  Two AI skeptics checked it: one restated the theorem in its own words and Lean accepted the proof for it, and built
  a one-array program that meets the theorem's hypothesis, so the theorem is not vacuous; the other rebuilt
  everything from source in a fresh copy, replayed it with Lean's independent checker and recomputed 224 twice
  outside Lean. The 13-October version of these files proved the bound only given the demand as a hypothesis; the
  14-October version removed it. **Rebuilt again for this revision** on 19 October 2026 from the published files, in
  a fresh copy of the same Lean tree: the build passes in 11 minutes with the same three axioms, and the negative
  control is still rejected (`programs/lean-block2-block3/logs/rebuild_2026-10-19.txt`).
- **In the widened game** (every address permutation and every phase pattern free): the same bound. The argument
  is Theorem R with the hook, which every permutation and every diagonal map keeps. **Argued by hand, not checked**;
  three AI skeptics upheld it: one re-derived Theorem R line by line against the Lean engine, one checked in exact
  arithmetic that 199 free steps (Toffoli, CCZ and controlled-S among them) keep the subspace and found no violation
  in random words, one redid the numbers. The numbers themselves (80, 144, 224) are **checked by computer** by
  certificates with separate checkers. It is not in Lean, because the Lean engine's relabellings are affine only.
- **Grouped moves.** Two moves in a row raise the block-3 hook by at most 144 (**checked by computer**, a
  certificate of 136 exact vectors). So the cap of 6.28e-2 also holds when grouped blocks are priced the way the
  engine prices them.
- **Excluded at block 3:** 2 in 5, 3 in 7, 3 in 8 (the old moonshot), 4 in 11, 5 in 13, 6 in 16, and every unit with
  4 or fewer arrays. The bound is also new in the engine game, where the earlier table gave 105/41 = 2.561.
- **Left open at block 3:** only units on or above the bound. The smallest are 5 in 14 (6.28e-2, exactly on the
  bound, section 3) and 6 in 17 (5.20e-2).

### 1.4 Block 4: at least 41/13 moves per array

Every unit at block 4 needs at least 41W/13, about 3.154W, paid moves. The certificate is the hook with nine columns.
Its cap of 27,456 per move is exact (Lean proves it as an upper bound by the cap formula, and a separate exact check
confirms that it is attained in all 15 directions); its demand of at least 86,592 = 16 × 5,412 comes from a 5,412 ×
5,412 minor whose determinant is nonzero modulo the prime 2,147,483,629, verified by a checker that shares no code
with the program that found the minor. Three programs written by different agents found the same demand (one of them
in floating point), a skeptic rebuilt it with its own code from exact minors, and the certificate program above found
the same rank again. **Checked by computer** (the counting theorem and the cap are proved in Lean; the demand is a
computer-checked certificate). **Re-run for this revision:** the demand certificate's determinant was checked again on
19 October 2026 (the same value, 1797789153), and `hookG.py` found the same demand by its own route. In the widened
game the bound also rests on the keptness argument of 1.3 (argued by hand).
- **Excluded:** 2 in 6, 3 in 9, 4 in 12, and every "W arrays in 3W moves". Any block-4 saving is at most 0.171.
- **Grouped moves:** two moves in a row raise the hook by at most 48,576 (**checked by computer**, a certificate
  with a 3,399 × 3,399 minor). The value 68,320 for three moves in a row is computed modulo a prime but not
  certified, so the engine-price cap of 0.176 at block 4 holds only with high probability.
- **Open:** 2 in 7 (9.63e-2), 3 in 10 (0.132), 4 in 13 (0.150), 5 in 16 (0.161). The searches for them are told in
  `../moonshot-block4-block5/`.
- **Not in Lean.** The block-3 route would need 16 character blocks of rank 5,412 each, roughly 10^7 times the
  kernel work of block 3.

### 1.5 Block 5 and beyond

- **Engine game, every block from 4 on: at least 2.9447 moves per array**; 2.99172 at block 5 and 2.97290 at block
  6. The method carries the block-4 hook upward through projections between frame graphs, controlled by a coupling
  argument. **Argued by hand, not checked**, from exact computations; two AI skeptics reproduced every number with
  their own code. Newly excluded: 3 in 8 and 4 in 11 at blocks 5 and 6; 5 in 14 at block 6; W arrays in 3W - 1 moves
  for W up to 120 (block 5) and 36 (block 6). The largest saving left open in the engine game is 0.319 at block 5
  and 0.392 at block 6; at block m it is at most ln(m/2.9447)/ln m.
- **Widened game:** 62W/29 at block 5 and 63W/31 at block 6, **checked by computer** (full-row certificates);
  31/14 at block 5 is computed exactly but not certified.
- **Depth bound:** every unit's total rank is at least W + m - 1, in both games, over any field, and for
  approximate units too. **Argued by hand, not checked; two AI skeptics upheld it.** It closes 1 in 3 and 1 in 4 at
  block 5, and several one- to three-array shapes at block 6. A later, stronger bound (W + m for W ≥ 2) is in
  `../moonshot-small-blocks/`.
- **Approximate units.** Every bound of this kind also holds for families of programs that are only correct in a
  limit (the "border barrier"). **Argued by hand, not checked; two AI skeptics upheld it.** The reverse side: if an
  approximate 5 in 14 exists, no argument of this kind can ever exclude it.

### 1.6 Two families excluded by certificate

- **The commutative slice** (only shifts, moves and mixing): no saving at block 3; at block 4 every scheme with
  fewer than 15W/4 moves fails (2 in 7, 3 in 10 and 4 in 13 among them); at block 5 every scheme with fewer than
  31W/8 fails. One probability distribution over test functionals per block settles every W at once, and all
  20,058,300 cases of 4 in 13 at block 4 were also checked one by one. **Checked by computer**, given an earlier
  reduction of the slice to a test on direction multisets that was not re-derived. So a unit must use address maps
  or phases. (Later work showed it must use both together at block 3; see `../reed-muller-and-one-move/`.)
- **Single-term contents** (the depth lemma): if every paid move acts on a content that holds at most one term of
  each input, two arrays need at least 2m moves, at every block size and whatever the free group. **Argued by hand,
  not checked**, with Farkas certificates **checked by computer** for every m; one skeptic's own program agrees on
  the hardest case. So a 2 in 7 at block 4 must put two different terms of one input into a single content.

### The table

| block | game | best bound (moves per array) | status | largest saving left open | smallest open targets |
|---|---|---|---|---|---|
| 2 | engine | 2 | proved in Lean | 0 | none |
| 2 | every game with bounded free steps | 2 | argued by hand, one skeptic | 0 | none |
| 3 | engine | 14/5 = 2.8 | proved in Lean | 6.28e-2 | 5 in 14, 6 in 17 |
| 3 | widened | 14/5 | argued by hand, three skeptics; numbers checked by computer | 6.28e-2 | 5 in 14, 6 in 17 |
| 4 | engine and widened | 41/13 = 3.154 | checked by computer | 0.171 | 2 in 7, 3 in 10, 4 in 13 |
| 5 | engine | 2.99172 | argued by hand, two skeptics | 0.319 | 2 in 7 (0.222), 3 in 9 (0.317); 2 in 6 is closed by the two-layer bound W + m |
| 5 | widened | 62/29 = 2.138 | checked by computer | 0.528 (0.506 by the exact value 31/14, which is computed but not certified) | 2 in 7 (0.222), 3 in 8 (0.391), 4 in 9 (0.496); the two-layer bound excludes 2 in 6 and 3 in 7 (`../moonshot-small-blocks/`) |
| 6 | engine | 2.97290 | argued by hand, two skeptics | 0.392 | 2 in 8 (0.226), 3 in 9 (0.387); 2 in 7 is closed by the two-layer bound |
| 6 | widened | 63/31 = 2.032 | checked by computer | 0.604 | 2 in 8 (0.226), 3 in 9 (0.387); the two-layer bound excludes 2 in 7, 3 in 8 and 4 in 9 |
| m ≥ 4 | engine | 2.9447 | argued by hand, two skeptics | ln(m/2.9447)/ln m | |

## 2. What the certificates can never do

The bounds above stop where the method stops. These limits matter as much as the bounds: they say which questions
a better certificate could still settle and which need a construction or a new kind of argument.

- **Ceilings.** The hook family can never certify more than 4(n - 1)/(n + 2) moves per array, n = 2^m: 10/3 at block 4
  (argued by hand), and the k = 8 hook almost surely stops at 41/13 (three codes agree, computed with high
  probability). Every coordinate certificate that all address permutations keep is capped at 4(n - 1)/n (argued by
  hand for m = 3 to 10). These ceilings do not cover certificates that only the engine game keeps, such as lifts and
  affine-orbit sets, which can do better (at block 4 the affine-only ones have no ceiling below 4). Two-factor Venn
  sets never beat the hook (computed: complete for m = 3 to 5). At block 4 no family screened reaches 3.5, so no hook
  or two-factor certificate can ever exclude 2 in 7 or 3 in 10. **Argued by hand, not checked; upheld by three AI
  skeptics**, who corrected several labels (for example "ceiling exactly" became "at most").
- **From block 5 on, these certificates cannot rule out a saving of 0.1.** Ruling out 0.1 at block m needs a bound
  above m^0.9 (3.48 at block 4, 4.26 at block 5, 5.02 at block 6), and the hook, cut and coordinate kinds all stop
  below 4. The saving they must leave open is about 0.16 at block 5, 0.24 at block 6 and 0.29 at block 7. Certificates
  built from whole non-repeated pieces (cuts) are capped below 4 at every block (Theorem S′, two skeptics; told in
  `../theory-of-theories-and-logic/`). **Argued by hand, not checked; upheld.** At engine prices the room is larger:
  the hooks cannot cap the saving below 0.83 at block 5 or 0.91 at block 6.
- **The best bound Theorem R can ever give** (the duality track). Relaxing Theorem R over every rank inequality
  valid over the complex numbers, including the Dougherty-Freiling-Zeger inequality, the best possible bound lies
  between 2.8 and 169/57 = 2.9649 at block 3, and between 3.154 and 3.94261 at block 4 in the widened game; from
  block 4 on it is at most m - 0.0574. So Theorem R can never exclude a saving of 0.0107 or less at block 3. All of
  the dual's slack comes from rows that read as "Fano gadgets" ending in a halving step. **Argued by hand, not
  checked**; the exact certificates (18 rows at block 3, 173 at block 4) were checked only by the track's own
  verifier, and no skeptic has read the track.
- **The ladder above Theorem R.** Reversing a unit and transposing every step gives a unit with the same moves, so
  every first-move lemma is also a last-move lemma. Beating 2.8 at block 3 needs a certificate whose two-move growth
  exceeds 1.8 times its one-move growth; excluding 2 in 7 at block 4 needs more than 1.75 and three-move growth above
  2.5 (the block-4 hook has 1.77 and 2.49). The guess "never more than 1.8" is false: the certificate Λ² ⊗ Λ³ reaches
  1.848 (a skeptic confirmed it exactly at two primes), but its three-move growth is only 2.515. Sums, mixtures,
  restrictions, reversal, single cuts, tensor squares and any linear programme over the scheme's subspace dimensions
  all give exactly Theorem R again. **Argued by hand, not checked; two AI skeptics upheld it**, with a correction:
  the claim that excluding 5 in 14 or 2 in 7 "needs one of" a short list of tools went too far.
- **Induction from small blocks.** Projection from block 4 can raise the bound by a factor of at most 1.1, pullbacks
  by less than 2, and direct sums never beat their best part, so no chain of these operations turns a fixed
  computation into a bound that grows with m. Whether the best Theorem-R bound grows in step with m is open (the
  track guesses it does not). **Argued by hand, not checked; upheld by two skeptics, one of whom refuted the claim
  "general products are at most twice the best factor"** with a counterexample (it holds for certificates of
  constant dimension).
- **"A unit exists" as equations.** Any run of moves on one array between two mixing steps equals free steps, then
  the first r coordinate moves, then free steps; so "a unit exists" is one polynomial system over the group algebra
  of the free steps, and each representation gives a relaxation. Deciding 2 in 7 this way is out of reach (about
  3.6 × 10^10 unknowns per step). **Argued by hand, not checked; no skeptic read it.** Its programs are not published
  here.

## 3. Toffoli and the other new relabellings

### 3.1 Are they allowed?

Two generators noticed, independently, that OpenAI's cost model treats linear and non-linear relabellings the same
way. The engine uses linearity in one place only: the pass that moves addresses computes each new address bit as an
XOR, by a look-up in an XOR table. A Toffoli also needs an AND, and an AND table is built like the XOR table at the
same cost. Every permutation is then a fixed chain of such gates and every phase pattern a product of table
look-ups, so each new free step is one pass over the data, like the linear ones; OpenAI's model charges every
arithmetic step one unit and never asks whether an index map is linear. **Argued by hand, not checked** (read from
the Lean files and OpenAI's model; no skeptic was assigned). Using such a step in a Lean proof would need about 400
to 450 lines in the engine's free-step layer; they were not written, because no unit appeared. The plain-English
account of this idea, as written on 13 October with the figures brought up to date, is
[toffoli-explained.md](toffoli-explained.md).

### 3.2 They buy nothing at blocks 2 and 3

The block-2 and block-3 bounds of section 1 hold in the widened game. On 13 October the first certificates in the
widened game had reached only 7/3 moves per array at block 3, which reopened 2 in 5 (0.166), 3 in 7 (0.229) and
3 in 8 (0.107); the hook closed all three the same evening.

### 3.3 Block 4 in the widened game

- **Many more routes for one array** (**checked by computer**; a skeptic confirmed it with an independent model).
  A route is a shortest way to finish one array alone. One array at block 4 has 736,785 routes in the widened game,
  against 315 in the engine game (105 against 21 at block 3; at block 2 the games coincide). The symmetry group of
  the route set is exactly the engine's affine group; the extra routes come from splitting: after a halving the work
  left is two independent half-size transforms, and a non-linear relabelling can send each half along its own route.
- **No help with sharing so far** (**argued by hand, not checked; one AI skeptic upheld it**). Junk left by one
  array lands in a place fixed by the relative twist, whatever route the other array takes, and one shared register
  with single-term private chains never pays at block 4 (moves ≥ 3m/2 + L). A gain only the widened game could make
  would need contents with two different terms of one input, with junk cancelling between branches or registers.
  Nobody has tested that case.
- **A beam search found nothing** in 12 target runs (2 in 7, grouped 2 in 7, 3 in 10, 4 in 13), and its controls
  passed. **Refuted as written** by its skeptic: the grouped runs never tested a scheme that could be finished, so
  they say nothing about the grouped target; "the obstacle is junk" holds only for rank-1 2 in 7; 77 % of the 4-in-13
  final tests were undecided. A sampled search tried 7,764 random 2-in-7 schemes and found none; that is a sample,
  not a search (**open**).
- **The grouped schedule 2, 2, 2, 1** (two arrays, three rank-2 blocks and one move; it would save 0.1675) is exactly
  tight on the block-4 hook. The f = 9 hook alone cannot exclude it, because a "shadow" scheme that lives only in the
  certificate's space passes every test the certificate can make (Theorem V; **argued by hand, not checked; one AI
  skeptic upheld it**), and tightness forces a rigid shape: no tested single-term content of depth 1 to 3 passes the
  zero-waste test (about 400 terms, and 76 more tested by a skeptic with its own code, two primes; **checked by
  computer**; a sample, since the free steps are far too many to enumerate). On 16 October a later round combined this
  with the term-by-term filter of section 3.5 and found that no exactly tight unit exists at block 4 in the engine
  game (**argued by hand from an exhaustive computer check of the 310 positions of depth 1 and 2; one program, no
  skeptic**; the grouped schedules only with high probability, because the rank-2 cap is not certified). That is told
  in `../idea-rounds/`. The widened game is not covered by it.
- **One phase between two moves removes a loss.** At block 4 the second of two plain moves in a row gains only 24,288
  of the cap 27,456. A skeptic computed that one controlled-S phase between the moves removes the loss entirely
  (**checked by computer**). So an earlier claim that plain chains waste 17,184 of the 19,008 spare capacity of 2 in 7
  holds only for plain chains, and only 6,336 of it is established (**refuted in part**). In the widened game a
  candidate 2 in 7 should place phases between moves.

### 3.4 The flagship search for 3 in 8 at block 3 was stopped

Before 14W/5 was known, a flagship track set out to decide 3 arrays in 8 moves at block 3 exactly (worth 0.107 to
0.127). Written as equations it has about 2.3 × 10^8 rational unknowns. It produced tools that carry over: each
input may be pre-twisted by a free step that stays free after the transform, so only 64 kinds of twist matter when
an input first enters (5 when two fresh inputs meet); a quotient game 256 times smaller in which every lower bound
still holds; and free finishing moves. Every family it searched was empty (complete only for 2 shared moves followed
by private chains of up to 2 moves), and its controls passed (all 9,261 plain 3-in-9 schedules found, none for
2 in 5). **Argued by hand and checked by computer by its author; no skeptic.** It was stopped once 14W/5 excluded
its target; a planned 72-hour run could only have confirmed the bound.

### 3.5 5 in 14 at block 3: exactly on the bound

5 arrays in 14 moves (6.28e-2) is the one block-3 target that the bound leaves with no slack: 14 × 80 = 1,120 =
5 × 224. It is **open**. What is known:
- In a tight 5 in 14 every move is rank 1 and gains the full 80; no move's content is a single free image of an
  earlier output (such a move gains at most 72; **checked by computer**, a certificate with exact vectors); each input
  has at most two private moves; the fresh-entry family of the flagship contains no 5 in 14. **Argued by hand, not
  checked; upheld by three AI skeptics**, whose own code reproduced the key numbers.
- A term-by-term filter (from Theorem T of the mixed-units track): in the engine game only 30 of the 135 frames can
  hold a term, and from each of them 7 of the 14 kinds of paid step keep a term legal. For all 8 certificates on
  which 5 in 14 is tight, the target looks exactly like 7 other targets, so no certificate argument can decide
  5 in 14. **Argued by hand, not checked; no skeptic has read it.**
- A greedy 14-move program is tight on all 8 certificates at f = 5 and passes every term test there, yet fails at
  f = 6 by 112 out of 8,960 (0.175 of a move); all 112 come from one array's end condition. So checking each
  array's end condition at the next f is strictly stronger than the cheap tests. **Checked by computer** (one
  program, not repeated).
- 6 in 17 (5.20e-2) is tight on no known certificate.
- Later results (15-16 October): a 5 in 14 must use address maps and phase letters together; under the order-7
  symmetry none exists; a tight one lives in a subgame described by the group S8. These are told in
  `../reed-muller-and-one-move/` and `../idea-rounds/`.

## 4. Which assumptions does OpenAI's model really force?

On 13 October I asked whether we were still tuning, or questioning the assumptions baked into our games. A round of
agents listed every assumption behind my impossibility proofs and my engine, and asked of each one: does OpenAI's cost
model force it, and if not, what opens when it is dropped?

| assumption | forced? | what was found | status |
|---|---|---|---|
| free steps only move entries and multiply them by phases | no | cheap singular maps are useless (their terms must cancel, so they can be replaced by zero); no cheap invertible non-monomial map is known | argued by hand, not checked |
| a paid move acts on a whole array | no | at block 3 a partial move adds exactly its share of a whole move's cap (exact rank 20 for one fibre call), so partial moves open nothing there; block 4 untested | checked by computer, one program |
| the same map on every digit column | no | column swaps, a fixed number of column-specific variants and maps on the end columns are useless | argued by hand; the end-column part later argued again and upheld by AI skeptics (section 5) |
| free steps never couple neighbouring columns | no | chain maps: section 4.1 | |
| correct for every f | yes, for all large f | correct at infinitely many f is already enough for the bounds | argued by hand; later upheld by AI skeptics (section 5) |
| exact correctness | yes | approximate units obey the same bounds anyway | argued by hand, two skeptics |
| each array gets the transform of its own independent input | partly | correlated inputs gain nothing; spare rows: section 4.2 | argued by hand, not checked |
| the same unit at every level of the recursion | no | mixing units by level only averages them | argued by hand, not checked |
| arithmetic over Q(i) | no | the bounds hold over any field | argued by hand, not checked |
| free steps cost nothing, paid moves a recursive call | every step is charged | the line between free and paid is work per entry: a fixed map of one column's digit, in all columns at once, costs O(1) per entry (OpenAI's explicit paper, Lemma 2.5), while a paid move's size grows with the problem; even per-entry work growing slower than the paid part leaves the exponent unchanged | argued by hand from OpenAI's model and papers, not checked |
| the kernel is the 2-point butterfly | no | OpenAI's finite-win game: section 4.3 | |

No loophole was found: data are only added, subtracted and scaled, so every program is a charged linear
straight-line program, and packing several values into one register cannot be undone without a floor or a
comparison, neither of which exists (**argued by hand from reading the model; not checked**). Three traps are not
loopholes: least squares converging to "solutions" that exist only in a limit, checking correctness at a few small
f, and searching over the field with three elements, which cannot rule out units over the rationals.

### 4.1 Chain maps: the open door, and how far it closed

**They are allowed** (**checked by computer**; one AI skeptic re-ran it). Two independent simulators wrote chain
maps as programs in OpenAI's model and counted every operation: work per entry stays constant as f grows (block
sizes 2 to 4, up to 7 columns). A paid move seen through a chain map is still exactly one transform call per fibre
between two free passes (argued by hand for every f, checked on raw arrays in 143 cases). The Fourier result uses the
kernel program as a black box, so any faster kernel program carries over. In this repository's Lean pipeline the new
step would need about 300 to 600 lines, plus 300 to 3,500 for correctness at every f (estimates).

**Correctness at every f is a finite test** (**checked by computer**). A word whose letters are finite-state across
the columns is a matrix product operator, and two such words agree at every f exactly when two weighted automata are
equivalent, which an exact computation decides (Schützenberger 1961, Tzeng 1992). Testing a few small f proves
nothing: the map `row0 ^= row1 << e` composed with the transform agrees with the transform up to f = e and first
differs at e + 1, and the test finds it.

**They break every certificate of the usual kind** (**checked by computer** at block 2 for 2 to 5 columns and at
block 3 for 3 columns; **argued by hand** for every f at block 2 by a skeptic). With one chain map free, no subspace
other than zero and the whole space is kept by all free steps. A chain map and an ordinary shift together give a
shift of one column only, and a paid move seen through `row ^= row << 1` contains a new paid operation, a chain of CZ
gates between neighbouring columns.

**A wall: Pauli-blindness** (**argued by hand, not checked; checked by computer up to f = 8 and by the automaton
test; upheld by two AI skeptics**). With `row ^= row << 1` free, the phase i^(parity of one address row over all
columns) is a free step, and at every f it equals (1+i)/2 + (1-i)/2 × (the column-wide sign flip of that row). So
the shortcut of judging every free step at f = 1 is not valid in this game, and any certificate of the usual product
form must treat every column-wide Pauli as doing nothing. The certificates behind 2 and 2.8 are not of that kind.

**Block 2 was closed again**, in two steps. First (14 October), the frames of the reading at alpha are the 27 lines
on a cubic surface, and the Fermat cubic gives an exact copy in complex space whose counts equal the frame distance
on all 729 pairs: a perfect certificate. With two lemmas about readings, this gives at least 2 moves per array at
block 2 for every unit right at infinitely many f whose free steps are Clifford steps of the listed kinds with a
reading at alpha (or all with a reading at 1). **Argued by hand, not checked; two AI skeptics** (both re-ran the
27-line computation, one checked the proofs and narrowed the scope to the listed kinds). The next readings (the 16-
and 8-element fields) admit no certificate of the 27-line kind, by Hirzebruch's inequality on line arrangements and
by Klein's classification of finite rotation groups (**argued by hand; upheld**). Then (15 October) Theorem 2W closed
block 2 in every game with bounded free steps (section 1.2), which covers every chain map, non-Clifford steps
included. Maps that repeat every 2 columns turn block 2 into block 4 with paired directions; there a unit needs at
least 1.7826 pairs per array (**checked by computer**, a skeptic recomputed it with separate code), which excludes
2 in 3, 3 in 5 and 4 in 7; the remaining 5 in 9 (0.152) is covered by Theorem 2W, since that game is block 2 with
free steps that touch one entry per row (argued by hand, one skeptic).

**Block 3: closed for one kind of chain map, open for the other.**
- Every free step a Clifford step with a reading at 1, for example `row0 ^= row1 << 1` and its inverse: at least 2.8
  moves per array, because the block-3 hook depends only on frames (its counts are 80, 144 and 224 at frame distances
  1, 2 and 3), together with a proof that the reading at 1 behaves as needed for every word. So 3 in 8, 3 in 7 and
  2 in 5 are excluded there. **Argued by hand, not checked; both halves checked by AI skeptics** (one recomputed the
  counts exactly over the rationals, one modulo two large primes). This replaces the figure 105/41 = 2.56 that the
  fifth-result notes of the morning of 14 October quote for this case. A later round extended 14W/5 (and 2W, 41W/13)
  to every chain map with a polynomial inverse, when every other step is Clifford (`../theory-of-theories-and-logic/`,
  Theorem E+; argued by hand, one skeptic).
- With `row ^= row << 1` or its running-XOR inverse free (no reading at 1): the best certificate gives
  231/131 = 1.763 moves per array, from the 405 planes of the Fermat cubic fourfold plus lifts of the other 486 of the
  891 frames, exact over the rationals with a cube root of unity (**checked by computer**; a skeptic re-ran the
  counts). Two natural constructions are capped below what is needed: certificates from the Weil modules at most 2.18
  (computed, not checked) and certificates that keep the 405 Fermat planes at most 2.287 (reproduced by a skeptic).
  Lifts of the frames are capped only at 2.917 (argued by hand, a counting argument over Fano subplanes), which is
  above 2.667, so counting alone does not rule them out; a later round lowered this cap to 147/52 = 2.827, still
  above 2.667 (checked by computer, two AI skeptics; `../theory-of-theories-and-logic/`). **3 in 8 is open in this
  game** (more than 2.667 is needed), and so are 2 in 5 and 3 in 7. A later bound caps every certificate in this game
  at 169/57 = 2.965, which leaves room above 2.667 (an exact certificate checked by computer and upheld by two AI
  skeptics, given a figure of the duality track that no skeptic has checked; `../theory-of-theories-and-logic/`).
- **A weak bound for every cheap relabelling** (chain maps, carries, affine maps of all the address bits): if no
  mixing coefficient has an even denominator and the unit is also correct at f = 2, it needs at least 1.5 moves per
  array at block 2 and 1.875 at block 3. **Argued by hand, not checked; checked by AI skeptics**, one of whom wrote
  out a missing lemma. It excludes no target (2 in 3 sits exactly at 1.5), and it says nothing about units whose
  mixing divides by 2.

**No unit was found** (**open**: evidence, not proof). Least squares over full templates tried 2 in 3 and 3 in 5 at
block 2 and 2 in 5 at block 3; the best residual for 2 in 3 was 0.124, against about 0.18 for an engine-only control
that is impossible by proof. A skeptic rebuilt every step (99 of 99 checks) and re-ran it with a new seed. Exact
searches with fixed shapes found no 3 in 8 (869 random patterns; 0 of 256 assignments in each of two shapes) and no
2 in 5 (2,473 patterns; 0 of 3,072), but the same search found none of the existing units on 1,095 random 2-in-6
patterns, so its random rows are a sample. One apparent hit was an artefact of testing at small f and failed on fresh
inputs.

### 4.2 Spare rows

OpenAI's own network leaves some address rows untouched, so its target is a rank-k target. In the game that network
lives in, every unit that delivers a rank-2 target to W arrays needs at least 2W moves, and a rank-3 target at least
3W, for every number of spare rows. **Argued by hand, not checked**; the proof is Theorem R with a small certificate
built from quadratic forms over the two-element field, and one AI skeptic re-derived it by hand and tested 150 random
words in each of four sizes. To escape it a unit must use a phase i^(z·x) with z ≠ 0, or an address map whose active
rows are not a permutation of the active bits (the skeptic added: or a map that mixes the active bits among
themselves). With every phase and address map free the certificates weaken as
rows are added: 1.8485 (rank 2, 3 rows), 1.7692 (rank 2, 4 rows), 2.4883 (rank 3, 4 rows), with floors of 1.5 and
1.75 for any number of rows (**checked by computer**; 1.7692 reproduced by a skeptic's separate code, 2.4883 from one
code family). Later, Theorem 2W′ showed that rank-2 targets never save in any game with bounded free steps
(`../moonshot-small-blocks/`), and a Venn certificate excluded 7 in 13, 8 in 15 and 9 in 17 with one spare row
(`../theory-of-theories-and-logic/`). **Open:** the rank-3 target 2 in 5 with 4 rows (0.166), which the certificate
misses by 0.012 moves per array.

### 4.3 OpenAI's own finite-win game

- **Hadamard kernel on 4 and 8 entries: no word beats the plain count.** With free shuffles, calls spread
  information exactly as in the classical gossip problem, where everyone learning everything among N ≥ 4 people
  takes 2N - 4 calls (Tijdeman 1971; Baker and Shostak 1972; Hajnal, Milner and Szemerédi 1972). The target is dense,
  so a word needs at least 4 calls on 4 entries and 12 on 8, which are the plain counts; this holds for every dense
  2 × 2 kernel. **Argued by hand, not checked; no skeptic.**
- **On 16 entries, 28 to 31 calls are open** (they would save 3.9e-2 down to 9.2e-3). Every word with 31 or fewer
  calls avoids every "layer boundary" state, so no Cooley-Tukey-shaped word can be shortened (argued by hand); the
  case reduces to one unproved inequality, "a Kronecker layer followed by a direct sum cannot share calls". The
  28-call gossip scheme passes density and orthogonality but fails a rank condition. Annealing found no word below
  32, but with its strongest objective it also failed to re-find a 32-call word, so its failures carry no weight.
  **Open; no skeptic.**
- **Subset-sum (zeta) kernel [[1,1],[0,1]] on 8 entries: no word with 11 or fewer calls**, so Yates's 12 additions are
  optimal there (11 calls would have saved 3.5e-2). A cut bound excluded 5,528 of the 6,828 classes of output order
  with exact integer certificates; an exhaustive search covered the other 1,300 and found nothing; words without
  cancellation need 12 calls (a separate breadth-first search). **Checked by computer**; one AI skeptic re-ran it on
  another machine and wrote an independent search for the hardest classes (its files were not kept). **Re-run for this
  revision** (section 9).
- The gossip bound leaves room that peaks at about 0.109 around 12 bits, so the game as a whole could clear the
  bar in principle; 16 entries is its smallest open case.

## 5. Conjecture U: one bound for every cheap free step

Instead of closing one lead at a time, a round tried to prove one statement for all cheap free steps at once.
**Conjecture U:** making more steps free (anything OpenAI's model does for a fixed amount of work per entry) never
pushes a unit below the engine's bounds of 2, 2.8 and 41/13 moves per array at blocks 2, 3 and 4. If it were true,
block 2 would allow no saving and block 3 at most 6.28e-2. It is **open** in general.

- **Shown for every column-uniform free step plus any end-column steps.** Free steps may be any relabelling that
  is the same on every column (Toffoli included) with any phases, plus any fixed invertible operators on the first and
  last few columns (S, T or H on the last column, a CZ across the last two, even a free transform on one end column).
  Then every unit right at every large f needs 2W, 14W/5 and 41W/13 moves at blocks 2, 3 and 4. **Argued by hand, not
  checked**: two independent proofs, each checked by two AI skeptics; one tested 38 to 40 words that mixed uniform
  and end-column steps, including random non-unitary end operators. So end-column tricks cannot rescue the Toffoli
  lead.
- **Theorem R-infinity.** The counting theorem needs correctness only at infinitely many f, run at one large f, if
  the certificate treats as "do nothing" every free step whose normalised trace stays away from 0 (the trace
  condition) and the steps are unitary. **Argued by hand, not checked; upheld by AI skeptics**, with one correction:
  the phases need entries in Z[ζ8, 1/2], not just in Q(ζ8); the skeptic's counterexample is (3 + 4i)/5.
- **The wall** of section 4.1 (Pauli-blindness) is where the two leads meet: the Toffoli certificate needs to see
  Paulis, and the only Pauli-blind certificates known here work through readings, which exist only for Clifford steps.
- **Refuted side claims of the round:** "the end-column theorem is a special case of the general reduction" is
  false (the global phase i^f breaks that route); the proposed lemma for chain maps with polynomial inverses is false
  as stated; a lemma on families with banded label action is false as written (a skeptic's counterexample: two
  global-parity maps whose product leaves the class); "the f = 1 shortcut is dead however the runs are rewritten" was
  narrowed to "it is not a valid certificate in this game".
- **No counterexample.** A hunter assumed U false and searched for 2 in 3 at block 2 with chain maps and end-column
  steps (least squares at two f, re-checked at f = 3 to 6); best residual 0.115, no hit (**open**: one start per
  pattern). Its first, interrupted run left an apparent hit that was an artefact of small f.
- **What is left of U** after Theorem 2W and Theorem E+: block 3 with the running XOR (section 4.1), non-Clifford
  steps together with chain maps at block 3, carries and other finite-state rules, and at block 4 every column mixing
  that Theorem E+ does not cover, the running XOR among it (not studied).

## 6. Smaller results of 13 and 14 October

- **Hub-side levers closed** (argued by hand, not checked): folding small blocks into the community's idle children
  is dead (0 blocks qualify); a weaker halving rule costs nothing in the community's chain but its template never uses
  it; Kerdock sets and mutually unbiased bases have copy charge exactly 1, against E8's 0.6.
- **Odd-or-zero line systems beyond E8** (the structure behind the E8 unit of the fourth revision) are closed at
  dimensions 10 and 12 by a heuristic search (best about 0.21 of E8's ceiling; **open as a proof**). One open object,
  1,023 lines in 31 dimensions, would beat E8 by a constant factor (about 1.4 times at 9 stages).
- **Splitting a move's capacity** (augmentation lemma; argued by hand, no skeptic): with every free step set to 1,
  the moves that reach the final frame have W independent "landing vectors", so a shared start can serve two arrays
  only if its branches land different parts of the shared content or a register holds two separate bundles. At block
  3 a single-term second move wastes 8 of its cap of 80 (computed at two primes).
- **Line lifts and twiddles** (computed, not checked): at block 2 a "line lift" gives a perfect frame family, which
  vanishes at block 3 (the Fano obstruction); if every per-column diagonal is free (like FFT twiddles), the best
  certificates found give only 1.83 to 1.88 moves per array at block 2 (floating point).
- **A free Toffoli between two moves** makes the second move perfect for the block-3 hook (exact over Q(i);
  **checked by computer**, one program, not repeated). Counting waste move by move cannot close 7 in 13, 9 in 17 or
  5 in 14, because the potential depends only on the state (argued by hand, not checked).

## 7. What failed or was refuted

Failures are as useful as successes, so here they are in one place. Where notes disagree, the later skeptic verdict
is the one used.

- **An "approximate unit" at block 2** was a modelling bug: the model dropped a factor of each term. Refuted by its
  author; no skeptic could be spawned.
- **"3 in 8 at block 3 is open even in the plain engine"** (a chain-map round's report): refuted by 14W/5, which
  excludes it there. It is open only with the running XOR free.
- **"With chain maps only 105/41 carries over at block 3"** (the morning summary of 14 October in the fifth-result
  notes): superseded the same day by 2.8 for every Clifford chain map with a reading at 1; the open case is the
  running XOR, with best bound 1.763.
- **"2 in 3 at block 2 is open with chain maps"**: closed by Theorem 2W (argued by hand, one skeptic).
- **The block-4 beam search's claims**: refuted as written (section 3.3).
- **The block-4 waste count of plain chains**: overstated; refuted in part by one controlled-S phase (section 3.3).
- **"The grouped 2,2,2,1 at block 4 is the most valuable small target"**: excluded in the engine game two days later
  (argued by hand from an exhaustive computer check, no skeptic; with high probability; section 3.3).
- **"The hook reaches 4(n - 1)/(n + 2) at every block"**: very likely false at block 4, where the k = 8 hook stops
  at 41/13 (three codes agree).
- **"General products of certificate families are at most twice the best factor"**: refuted by a skeptic's
  counterexample when the families' dimensions vary.
- **Four side claims of the Conjecture U round**: refuted or narrowed (section 5).
- **The Toffoli explainer's 7/3**: the first widened-game certificate; superseded the same evening by 2.8.
- **"Approximation never helps"** does not carry over to the everything-free game, whose running XOR generates an
  infinite group of commuting steps; that is told in `../moonshot-small-blocks/`.

## 8. What is still open, with a first experiment for each

| question | worth | first experiment |
|---|---|---|
| 5 in 14 at block 3 (and its multiples) | 6.28e-2 | Decide the S8 sub-game of `../idea-rounds/` exactly for one sector: an exact search with one representative per symmetry orbit, every candidate judged by the per-array end condition at f = 6 (`programs/five-in-fourteen-block3/`), not by the cheap tests. |
| 3 in 8 at block 3 with the running XOR free | 0.107 | Solve for the best Pauli-blind frame family on the 891 frames of the reading at alpha, as a linear programme over orbit-sums of the frames (the 169/57 cap leaves room above 2.667); in parallel an exact search for 3 in 8 with chain-map menus, correct at f = 2, 3 and 4 together, every hit confirmed by `every_f.py`. |
| 2 in 7 and 3 in 10 at block 4 by proof | 0.096, 0.132 | Compute the exact finite-group form of the duality relaxation, first at block 3 (does any certificate beat 2.8?), then at block 4: a certificate above 3.5 would close 2 in 7. The searches are in `../moonshot-block4-block5/`. |
| rank-3 target 2 in 5 with 4 spare rows | 0.166 | Repeat the averaging proof that closed 2 in 5 with 3 rows on the 2,295 Lagrangian frames of 4 rows, or try three-factor certificates (`programs/spare-rows/`). |
| Hadamard finite win on 16 entries, 28 to 31 calls | 3.9e-2 to 9.2e-3 | Prove "a Kronecker layer followed by a direct sum cannot share calls"; or an exhaustive structural search built like the zeta search (`programs/finite-win-game/zeta-kernel-8-entries/tdfs.c`), with the cell and singleton-state prunes. Cheap add-on: the zeta kernel on 16 entries (32 calls plain). |
| Conjecture U at block 3: Toffoli and chain maps together | up to 0.107 | Look for a Pauli-blind certificate that Toffoli keeps, or an identity that forces such a certificate to be empty, first with Toffoli plus `row0 ^= row1 << 1`. |
| Column mixing at block 4 beyond Theorem E+ (the running XOR, non-Clifford steps) | up to 0.171 | Not studied: run the period-2 computation of `programs/chain-maps-certificates/period-2-maps/` at block 4 with period 2 (block 8 with paired directions) as a first calibration. |
| 41W/13 in Lean | a Lean proof of section 1.4 | Split the 16 character blocks further by the bit-swap group of `hookG.py` (four characters, blocks of about 1,400), and time the Lean kernel on one block before committing. |
| Does the best Theorem-R bound grow with m? | the room at large blocks | Compute the duality bound at block 5 restricted to the frames within distance 2 of the start and of the target. |
| A skeptic for the duality track | method limits | Check the 18-row block-3 certificate by hand and re-run `verify3.py` from scratch with independent code. |
| 1,023 odd-or-zero lines in 31 dimensions | a better constant for E8-type units | A sign search with a SAT solver on a Linux machine. |
| Toffoli in the Lean engine | needed only if a Toffoli unit appears | The AND table and the widened free-step layer, about 400 to 450 lines. |

## 9. Programs and data

Each folder under `programs/` has a README with what it shows, the command, time and memory, and the expected
output. Small certificates and outputs are in `data/`; nothing here is larger than 5 MB. `results.json` lists every
finding with its status, figures and evidence files.

| folder | what it shows | re-run for this revision |
|---|---|---|
| `programs/lean-block2-block3/` | the Lean files of Theorem R, the hook, block 2 and block 3 (11 files, 1,662 lines), the certificate generator, the skeptic's restatement | yes: rebuilt from these exact files on 19 October on a Linux machine of Google's Colab service (11 minutes, 10 GB, the three standard axioms, negative control rejected); the generator re-run gives the same certificate |
| `programs/hook-certificates-blocks-2-to-6/` | certificate generators and separate checkers: caps, demands, rank-2 blocks, the chain gain, blocks 5 and 6 | yes: blocks 2 and 3, every cap, the chain gain, both rank-2 caps, the negative controls; and the block-4 demand determinant (56 minutes, the same value) |
| `programs/widened-game-hook-block3/` | the hook in the widened game; the skeptics' recomputation and keptness check | yes |
| `programs/block4-hook-41-over-13/` | 41/13 by symmetry-reduced ranks; the skeptic's rebuild | yes (`hookG.py`, 48 seconds) |
| `programs/commutative-slice-certificates/` | the slice exclusions for every W | yes (the LP certificates and 2 in 7) |
| `programs/family-S-depth-lemma/` | Farkas certificates for single-term contents | yes |
| `programs/engine-bounds-blocks-5-and-6/` | 2.9447, 2.99172, 2.97290 | yes |
| `programs/certificate-ceilings/` | ceilings of the hook and Venn families; the arithmetic of the block table | yes (`numbers.py`, `table.py`) |
| `programs/ladder-above-theorem-R/` | growth ratios c2/c1, c3/c1; Λ² ⊗ Λ³ | no |
| `programs/duality-best-theorem-R-bound/` | the duality certificates and their verifier | no |
| `programs/five-in-fourteen-block3/` | the single-term gain table (certificate) | yes |
| `programs/toffoli-routes-block4/` | route counts in the widened game; the skeptic's census | yes (`census.py` at blocks 3 and 4) |
| `programs/chain-maps-cost-and-every-f/` | cost of chain maps in OpenAI's model; the automaton test | yes |
| `programs/chain-maps-certificates/` | the 27 lines, readings, the 891 frames, the parity identity, the refuted banded lemma, period-2 maps | yes (`lines27.py`, `relcheck.py`, `e2_counter.py`) |
| `programs/spare-rows/` | the spare-row certificates | no |
| `programs/finite-win-game/` | the zeta-kernel search; the 16-entry Hadamard tools | yes (the zeta search, about 15 minutes) |
| `programs/tightness-at-block-3/` | exact hook gains; a free Toffoli between two moves | no |

Large outputs are not published: the block-4 zero-waste spaces (three files of 50 MB each), a 105 MB pickle of the
5-in-14 depth-3 terms, and the 40 MB and 3.9 MB certificates of the 4-in-13 and 3-in-10 slice checks; the programs
above regenerate them (their READMEs give the commands). Pickles are not published, because loading one can run
code.

## Credits

I did this work with a team of AI agents (Claude) that I directed: I set the goals, chose which ideas to pursue and
decided at each stage what to do next; the agents ran the searches, wrote the proofs and programs, and tried to refute
one another's work. The engine and the engine game are those of this repository and its earlier revisions; the hook
certificates were found in this work. Theorem R is a general form, with scratch arrays and singular mixing, of the
certificate lemma (its "FO lemma") of eumemic's pull request #288 to the community's repository
CrocSwap/integer-mult-bounds, prepared with AI assistance by its own statement; that pull request also names the graph
of frames and reports, for its model, that block 2 never saves. It credits its authors' earlier work and says that
papers of Swapnil Jain prove overlapping statements. The fourth revision's note for network-coding readers gives the
same credit. The Lean files were built in a copy of Ryan Shea's Lean tree (shea256/fourier-transform- below-nlogn,
folder `lean/complex588`), which pins the same Lean and Mathlib versions as this repository and carries byte-identical
copies of the engine files they import. The cost model, the explicit construction and the finite-win game are OpenAI's
(the Lean file `UniformFourier.lean` of openai/math, and OpenAI's preprints "An explicit power saving for the exact
discrete Fourier transform" and "Finite tensor savings and exact Fourier circuits"). Classical results used: the
gossip bound (Tijdeman 1971; Baker and Shostak 1972; Hajnal, Milner and Szemerédi 1972); the Dougherty-Freiling-Zeger
inequality (arXiv:1311.4601); Hirzebruch's inequality on line arrangements (1983); Klein's classification of the
finite rotation groups (1884); the 27 lines on a cubic surface (Cayley and Salmon, 1849); weighted-automaton
equivalence (Schützenberger 1961; Tzeng 1992); Yates's algorithm for subset sums. The parallels with convolutional
codes (Forney 1970), network coding (Koetter and Médard 2003) and Clifford quantum cellular automata (Schlingemann,
Vogts and Werner 2008; Haah 2013) are noted in the chain-map notes; claims about them are not made here.

## Sources

All of this section comes from my working notes of 13 to 16 October 2026, which are kept in a private archive;
paths are relative to its `checks/wht43/` folder. Archive catalogue entries are in brackets; they are also the
`archive_entry` fields of `results.json`.

- Round notes: `proofs/NOTE.md` (the six proof tracks: lean, logging, ladder, induction, duality, algebraic),
  `px/NOTE.md` (px-lean, px-border, px-t514), `wave16/NOTE.md`, `wave16/b56-note/NOTE.md`, `wave16/tb4-note/NOTE.md`,
  `wave17/NOTE.md`, `wave17/xc-note/NOTE.md`, `wave17/xc3-note/NOTE.md`, `wave17/ut-note/NOTE.md`,
  `TOFFOLI-EXPLAINED.txt`.
- Section 1 [L1, L6, L7, L9, L10, L11, L12, X14]: `proofs/lean/` (PROOF.md; skeptics v-pf-lean-logic,
  v-pf-lean-compute), `px/px-lean/` (PROOF.md, STATUS.txt, logs; skeptics v-px-lean-logic, v-px-lean-compute) →
  `programs/lean-block2-block3/`; `proofs/logging/` (foundation, matroid, familyS, t514R1) →
  `programs/hook-certificates-blocks-2-to-6/`, `programs/commutative-slice-certificates/`,
  `programs/family-S-depth-lemma/`, `programs/five-in-fourteen-block3/`; `wave16/d-tof-cert/`, `wave16/v28-theorem/`,
  `wave16/v28-kept/`, `wave16/v28-numbers/` → `programs/widened-game-hook-block3/`; `wave16/b56-compute/`,
  `wave16/b4-cert/`, `wave16/v-b56-compute/`, `wave16/b4-depth4/` → `programs/block4-hook-41-over-13/`;
  `proofs/induction/` → `programs/engine-bounds-blocks-5-and-6/`; `wave16/d-tof-search/NOTE.md`.
- Section 2 [T4, T6, T10, T11]: `wave16/b56-ceilings/`, `wave16/b56-methods/`, `wave16/b56-v-skeptic/`,
  `wave16/v-b56-methods/` → `programs/certificate-ceilings/`; `proofs/ladder/` → `programs/ladder-above-theorem-R/`;
  `proofs/duality/` → `programs/duality-best-theorem-R-bound/`; `proofs/algebraic/` (not published);
  `proofs/v-pf-induction-logic/` (the product counterexample).
- Section 3 [M15, X13, X27, N3, C20]: `wave16/d-tof-cert/LEGIT.md`; `wave16/tb4-beam/`, `wave16/tb4-vbeam/`,
  `wave16/tb4-tight/`, `wave16/tb4-v-tight/`, `wave16/tb4-nonlinear/`, `wave16/tb4-nonlinear-skeptic/` →
  `programs/toffoli-routes-block4/`; `wave16/b4-search/`; `wave16/f-reduce/NOTE.md`; `wave16/t514-cert/`,
  `wave16/t514-search/`, `wave16/t514-tight/` and their skeptics; `px/px-t514/` (PROOF.md, greedy/);
  the 16-October recombination of the wht44 idea round (`../wht44/ideas/RECOMBINE.md` section 2.1).
- Section 4 [L16, L18, L19, L20, M15, X35, T9]: `wave17/g-rules/`, `wave17/g-hypotheses/`, `wave17/g-history/`;
  `wave17/xc-legit/`, `wave17/xc-v-legit/` → `programs/chain-maps-cost-and-every-f/`; `wave17/xc-algebra/`,
  `wave17/xc-v-algebra/`, `wave17/xc-round2/`, `wave17/v-xc-round2/`, `wave17/xc3-s1s2/`,
  `wave17/v-xc3-s1s2-compute/`, `wave17/xc3-cert3/`, `wave17/v-xc3-cert3-compute/`, `wave17/v-xc3-cert3-logic/` →
  `programs/chain-maps-certificates/`; `wave17/xc-search/`, `wave17/xc-v-search/`, `wave17/xc3-search3/`;
  `wave17/d-G3-SPECTATOR/more-spectators/` → `programs/spare-rows/`; `wave17/d-R4/T-kernel-b3/`,
  `wave17/d-R3/C-kernel-b4/` → `programs/finite-win-game/`.
- Section 5 [L17, C10]: `wave17/ut-normal/`, `wave17/ut-direct/`, `wave17/ut-hunt/` and their skeptics
  `v-ut-normal-logic`, `v-ut-normal-compute`, `v-ut-direct-logic`, `ut-v-direct-compute`.
- Section 6 [R26, R27, T16]: `wave16/d-I12/`, `wave16/d-absorb/`, `wave16/d-I3/`, `wave16/d-halving/`,
  `wave16/d-K6/`, `wave16/g-kerdock/`, `wave16/x-cross/`; `wave17/d-G3-SPECTATOR+G1/tightness-at-block-3/` →
  `programs/tightness-at-block-3/`.
- Section 7 [C8, C9, C10, C20]: `wave16/d-B1(g-play)+T2-BORDER-NUM/NOTE.txt`; `wave17/xc3-note/NOTE.md`; the fifth-
  result candidate's `notes/open-directions-3.md` section 8 (14 October, morning).
- Cross-references to later rounds: Theorem 2W and the two-layer bound [L4, L5] (`moonshot/w7-xwborder/`,
  `moonshot/mw-big56/`); Theorem C3 [L22] (`rmfam/rm-theory/`); Theorem E+ [L18] (`meta/cohomology/`); the chi-game
  cap [T9] (`meta/c2-hermitian/`); the Venn certificate with one spare row [X12] (`meta/literature/`).
- Credit for Theorem R's lemma: pull request #288 of CrocSwap/integer-mult-bounds (eumemic, opened 10 October 2026),
  as credited in this repository's `notes/for-network-coding-readers.md` and `notes/network-coding/two-arrays.md`
  (fourth revision); the first general statement in my notes is `../wht33/hub288/verify/v-certificate.md` section 3.
- The re-runs of 19 October 2026 were done on an Apple M4 Mac mini and on one Linux machine of Google's Colab
  service (8 processors, about 50 GB); their outputs are the `rerun_*` files in `data/` named in the program READMEs.
