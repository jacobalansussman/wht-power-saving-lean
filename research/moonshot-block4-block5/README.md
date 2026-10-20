# Moonshot waves 1-8 at blocks 4 and 5: the 2-in-7 near miss and the shape searches

On 14 and 15 October 2026 I ran eight "moonshot" waves and a side round on order-3 letters, looking for a small unit
that saves at block size 4 or 5. This folder tells what was tried there, what held, what broke, and what is still
open, with the programs that produced each result.

## In short

- **No unit was found** at block 4 or block 5, in any shape that was searched. Nothing in this section is proved in
  Lean, and none of the proofs below has been checked by a person; each status says who or what checked it.
- **The near miss is still open.** A 2-in-7 unit at block 4 in which each array ends on a private move is allowed by
  every bound I have, by a margin of 640 out of 137,280 (0.47 %). Order-3 letters, and every other odd-order letter,
  leave that margin unchanged. Many shapes inside this branch are now excluded; the branch as a whole is not.
- **Excluded at block 4, each under the hypotheses stated below** (targets 2 in 7, 3 in 10, 4 in 13): the symmetric
  2 in 7 with single-term shared and last moves; the Singer-cycle template; private layer-2 moves with one shared
  layer-3 move (Theorem P, single letters); for 2 in 7, Theorem A's start and shape R with one term per frame. Each
  of these is argued by hand and upheld by one AI skeptic; none is checked in Lean. **Narrowed, not excluded:** one shared move followed by two private layers (a second route there is only a change of
  basis, and 755 small skeletons of it have no solution), and the 2-in-7 searches (1,256 of their 1,261 saved prefixes
  excluded exactly against every end move).
- **The symmetric targets 7 in 23 and 8 in 26** (savings 0.142 and 0.150): the blueprint as written is impossible
  (argued by hand, upheld by a skeptic); 8 in 26 with a group-algebra middle move fails for every twist subgroup, and
  the 7-in-23 construction is excluded in both move orders in the configuration that was built (both checked by
  computer, each upheld by a skeptic). Its nearest variant, 7 in 24 (0.111), passes every distance-4 test computed and
  is open.
- **Block 5** (2 in 7, saving 0.222; 2 in 8, saving 0.139): the single-letter butterfly 2 in 7 is never a unit
  (argued by hand, upheld by a skeptic); every old 2-in-8 near miss fails exactly (checked by computer, one skeptic);
  one new 2-in-8 design, s1501, is undecided (open).
- **Two method lessons.** The forced-zero tool used in every block-5 wave from the fourth on missed forced zeros and
  kept dead designs looking alive; once it was corrected, every old near miss was decided exactly. A cheap test that
  forgets signs and phases (the frame-level test) cannot tell an invertible step from a zero divisor, and several
  "near misses" passed it only for that reason.

## The question

A **unit** is a small, fixed program that transforms W arrays at once with fewer costly steps than doing them one at
a time. At **block size m** each array must be transformed along m address directions; one **paid move** applies one
direction's 2-point transform to one array (inside the recursion it is one recursive call), and everything else in a
program is a **free step**: relabelling addresses, multiplying entries by phases, adding multiples of one array to
another. The plain schedule needs W·m paid moves. A unit "W in M" that needs only M < W·m gives, used recursively, an
algorithm for the Walsh-Hadamard transform with about N (log N)^(1-s) operations, where the **saving** is
s = ln(mW/M)/ln m (the earlier revisions explain how such savings carry over to the Fourier transform and to
convolution). A single unit is
all a saving needs. The standing bar for a new route is a saving of 0.1; since 15 October, 0.01 counts as a
stepping-stone tier.

Earlier rounds had settled the smallest blocks: block 2 never saves, and block 3 needs at least 2.8 moves per array,
so it saves at most 6.28e-2 (both proved in Lean for the engine game). Block 4 needs at least 41/13, about 3.154 moves
per array, so it saves at most 0.171 (checked by computer: four programs by different agents found the same count and
a skeptic rebuilt it from exact minors; not in Lean). These are told in `../block-bounds-and-new-rules/`. That left
these targets, which are what this section is about:

| block | target | saving | why it was open |
|---|---|---|---|
| 4 | 2 in 7 | 0.096 | below the bar, but the cheapest member of the families below, and the one near miss |
| 4 | 3 in 10 | 0.132 | allowed by every bound |
| 4 | 4 in 13 | 0.150 | allowed by every bound |
| 4 | 7 in 23, 8 in 26 | 0.142, 0.150 | symmetric programs suggested by the best certificate's dual |
| 5 | 2 in 7 | 0.222 | exactly on the bound W + m (two-layer bound) |
| 5 | 2 in 8 | 0.139 | a split-trunk 2 in 7 gives a 2 in 8, so it is easier and a target in its own right |
| 5 | 3 in 11 | 0.193 | allowed by every bound; not reached by any track |

An older block-5 near miss, 2 in 6 (0.317), which missed exclusion by 0.008 on the old certificates, is excluded by
the two-layer bound "W arrays at block m need at least W + m moves" of wave 2 (argued by hand, not checked; one AI
skeptic upheld it; see `../moonshot-small-blocks/`). So that near miss was a weak bound, not a hidden unit.
Several early notes quoted 0.317 for the block-5 2 in 7; the right figure is 0.222.

## Words used here

- **Games.** They differ in what is free. The *engine game* has the free steps of the Lean engine of this repository
  (linear address maps, shifts, Clifford phase patterns, rational mixing of arrays). The *widened game* also makes
  every address permutation and every pattern of phases 1, i, -1, -i free. Unless a line says otherwise, a result
  holds in the engine game.
- **Frame, depth, layer.** Free steps move a term between fixed positions called frames. At block 4 the engine game
  has 2,295 frames, and they form a well-known geometry (the dual polar graph of the symplectic group Sp(8,2)); a paid
  move takes a term to one of the 30 neighbouring frames. A term's depth (or layer) is its distance from the input
  frame; the finished target K x sits at depth m, on K's frame.
- **Shapes.** A *start* is an array's first, private move; a *shared move* (or *merge*) mixes several arrays; the *end*
  or *top* moves are the last moves; the *readout* is the free combination of moves that forms each output. A
  *private end* means each array finishes on a move of its own. *Leftovers* (junk) are terms that land in the wrong place
  and must cancel. An *injection* adds an extra input term to a later move.
- **Certificate, hook, cap, demand, room.** A certificate is a counting bound: a space of test functions that every
  free step keeps in place, a *cap* on how much one paid move can raise the count, and a *demand* each finished piece
  needs. The *hooks* are my best family of certificates, one for each size k from 1 to 15. The *room* is what the caps
  allow minus what the demands need; a unit needs room of at least 0. The counting rule behind every such bound
  (Theorem R) is proved in Lean for the real engine (see `../block-bounds-and-new-rules/`); the hook values themselves
  are computed, not proved in Lean.
- **Frame-level test, full test.** The frame-level test replaces every free letter by the number 1 and tracks only how
  much lands on each frame. Every unit passes it, so passing proves little. The full test checks the real identity.
- **Invertible, zero divisor.** A combination of free letters is invertible if another combination undoes it, and a
  zero divisor if it sends some nonzero input to 0. With sigma a letter of order 3, (1 + sigma)(1 - sigma + sigma^2) = 2,
  so 1 + sigma is invertible and undoing it divides by 2; with u a letter of order 2, (1 + u)(1 - u) = 0, so 1 + u is a
  zero divisor. The frame-level test turns both into 2 and cannot tell them apart.
- **Symmetric programs, orbits, characters.** A symmetry can rotate the W arrays. An *orbit* (O) is a set of W moves,
  one per array, that the symmetry permutes; a *fixed* or shared move (F) serves all arrays. Under the symmetry every
  test splits into W independent parts, one per *character*, and a unit must pass all W.
- **Skeptic.** An independent AI agent that wrote its own code and tried to refute a claim. "Upheld" means the skeptic
  tried and failed. Most claims in this section had one skeptic, not two; second skeptics for the results that matter
  most are still to do (see the open questions).

## What I tried

I ran the work as rounds of tracks. Each track had one question and a time box; when it finished, a separate skeptic
agent tried to break its claims, and I wrote up the round from both. In order:

1. **Wave 1** (14 October, morning): the Fano-gadget blueprint at block 4, suggested by the dual of the best
   certificate. It produced the near miss.
2. **Wave 2**: symmetric programs at block 4, which produced a second near miss (a swap-symmetric 2 in 7); then its
   analysis, the Singer-cycle template, and early-separation results at block 5. A skeptic of wave 1 had meanwhile found
   that an order-3 letter separates two arrays with only two pieces, so I started a side round on order-3 letters.
3. **The order-3 round**, three rounds: the hook bound with order-3 letters free, constructions with order-3 ends, an
   exact model of the block-4 engine game, a finished search of every saved 2-in-7 prefix, shape R, and shared ends for
   3 in 10 and 4 in 13.
4. **Wave 3**: the shared move early with two private layers after it (block 4); the butterfly at block 5; the
   blueprint targets 7 in 23 and 8 in 26.
5. **Wave 4**: the block-5 trunk with several branches.
6. **Wave 5**: the open orders of 7 in 23 and 8 in 26; two shared moves at layer 2 then one shared move at layer 3.
7. **Wave 7**: block 5 again (collision orbits, the 2-in-8 shape as a target).
8. **Wave 8** (15 October): the 7-in-23 construction decided in the 14-dimensional representation; the old block-5
   near misses decided exactly with a corrected forced-zero tool.

(Wave 6, and the block-2 and block-3 parts of the other waves, are told in `../existence-routes/` and
`../moonshot-small-blocks/`; the Reed-Muller family of wave 8 is in `../reed-muller-and-one-move/`.)

The rules throughout: predictions written down before computing, controls on known cases first, only exact checks
count as a verdict, every numerical fit re-checked on fresh inputs, and heavy jobs on Linux machines of Google's Colab
service (8 processors, about 50 GB each), not on my Mac.

## What I found

Each finding carries one of these statuses: **proved in Lean** (none in this section), **checked by computer** (an
exact or exhaustive computation; the line says who repeated it), **argued by hand, not checked** (a written proof,
checked neither by a person nor by a proof checker; the line says whether an AI skeptic re-derived it), **refuted** (and by
what), and **open**. Numerical searches that found nothing are evidence, not proof, and are labelled so. When two
readings disagree, the later skeptic's verdict is the one used, and the line says so.

### 1. The near miss: 2 in 7 at block 4 with private ends

If both arrays of a 2 in 7 end on private moves, cut those two moves off: the first five moves must deliver two
depth-3 pieces, one per array. Under the k = 8 hook (at f = 9 digit columns), one paid move is worth at most 27,456 and
one depth-3 piece needs 68,320. Five moves allow 5 × 27,456 = 137,280; two pieces need 2 × 68,320 = 136,640. The bound
closes this branch only if a depth-3 piece is worth more than 2.5 moves, and here it is worth 68,320/27,456 = 2.48834.
So 640 of the 137,280 are left over.

| finding | status |
|---|---|
| The k = 8 hook values 27,456 and 68,320, and the room 640. | **Checked by computer**: exact ranks modulo a large prime; an independent skeptic reproduced every number by another method, and a skeptic of a later round confirmed that 68,320 is the exact value and not only a lower bound computed modulo primes (floating-point singular values with a wide gap, and two more primes agreeing on every block). So the room is real. The move value 27,456 is the case m = 4, k = 8 of the hook's cap formula, which is proved in Lean (`../block-bounds-and-new-rules/`). I re-ran the four rank blocks behind both values on 19 October 2026 (`programs/hook-ratio-640/`). |
| No hook with k from 1 to 15 reaches 2.5 (k = 7 and 9 give 2.4591, k = 6 and 10 give 2.3766), and mixing hooks cannot beat k = 8. | **Checked by computer**; one skeptic. |
| Every free step that relabels addresses and multiplies by phases keeps the hook's test space in place, so making more such letters free (order-3 shifts, sign products, Toffoli, controlled-S, other roots of unity) can never change the 640. Every separating end of this shape asks the first five moves for exactly two depth-3 values. | **Argued by hand, not checked**; an AI skeptic re-derived it and found no error. |
| In the widened game the 640 rests on a guess: that targets outside the engine game are never cheaper at k = 8. It was computed only at k = 2, 3 and 4. | **Open**. |
| Values the first round gave for the odd hooks (k = 5, 7, 9): they multiplied one character block by 16, which is valid only for even k. | **Refuted** by the order-3 round, which ranked all 16 character blocks; the maximum is still at k = 8, so the near miss stands. |
| One closing sentence of an earlier note, which named the test a certificate must pass to close the branch. | **Refuted** by a later certificate search (the correct test is s3 < 0.7·s1, not s3 < 0.7; one orbit cut has the right spectral factor and still misses). That search, which tried about 4,000 certificates and found none better than the hook, is told in `../theory-of-theories-and-logic/`. |
| The near miss itself. | **Open.** No program for it exists, and no certificate closes it. |

There are two readings: the hook is slightly loose here, and a sharper count would close the branch; or a real unit
sits in the sliver. A unit there would have to waste less than 640 under the k = 8 hook, less than 1,964 under k = 7
and less than 3,952 under k = 6, all at once, which forces every move to be a multi-term mix that wastes almost nothing.
On present evidence I think the first reading is likelier, but nothing below decides it.

### 2. Order-3 letters: a two-piece separation that is, at the ends, only a change of basis

Let sigma be the cyclic shift of three address bits, a free letter of order 3. From the two pieces x1 + x2 and
x1 - sigma x2 one recovers K x1 exactly: Z1 (x1 + x2) + Z2 (x1 - sigma x2) = K x1, with Z2 = K (1 - sigma + sigma^2)/2
and Z1 = Z2 sigma. This works because 1 + sigma is invertible. Until then I had thought the smallest separating set
needed four pieces (a Fano gadget); that is true only when every letter relating the pieces has order 2.

| finding | status |
|---|---|
| The two-piece separation above, and its failure for order-2 letters. | **Checked by computer**: explicit matrices at blocks 3 and 4, an exact identity; found by the skeptic of wave 1. Re-run on 19 October 2026 (`programs/order-3-two-piece-separation/`). |
| An order-3 pair placed at the start or the end of a program is only a change of input basis. So the order-3 "halving" end is a private end in other coordinates and inherits the same 640. | **Argued by hand, not checked**; the skeptic of that round confirmed this part. |
| A correct 8-move program with the order-3 halving end exists. It is the plain cost: in the basis y1 = x1 + x2, y2 = x1 - sigma x2 it is two private chains. | **Checked by computer**: exact end solve, formal identity, random complex arrays at f = 1 to 5 (errors below 7e-16); the skeptic replayed it with its own code. Re-run on 19 October 2026 (`programs/order-3-halving-end-control/`). |
| Four shapes are dead for every choice of letters: one depth-2 carrier plus private helpers; single pieces with two-valued frustrated patterns; one depth-3 carrier for both targets; two depth-3 carriers that are not clean. A genuinely shared end for W arrays needs at least W + 1 pieces. | **Argued by hand, not checked**; the skeptic read each proof and found it correct within its hypotheses, and noted that the first one silently assumes the helper is shallower than the target. |
| "Sign patterns alone need four pieces, so order-3 letters do better at the end." | **Refuted** by the skeptic: with multi-term pieces, sign letters also reach three pieces for two arrays and four for three. |
| The first searches' claims to be complete. | **Refuted** by the skeptic: some end moves were never tested, and one template computed nothing although its prefixes were counted. Section 5 below finishes that search. |

### 3. The symmetric 2 in 7: not a unit, and why

A search over symmetric programs (a cyclic symmetry rotating the arrays) turned every block-4 target of the form
W in 3W + 1 into "three orbits plus one shared move", and all 216 skeletons it decided at target cost failed. Its best
object was a swap-symmetric 2 in 7 with the shared move in the middle, which met the frame-level test to 2.4e-10 while
two coefficients grew to 5.7e4. It looked like an approximate unit.

| finding | status |
|---|---|
| The growth was only rescaling. Normalised, the program is small and ordinary; it meets the frame-level test exactly and fails the full test by a fixed amount (the full residual tends to 2/3 and never reaches it; at the exact frame-level solution it is 1). | **Checked by computer**; the skeptic rebuilt the closure and the tests with its own code. |
| Why it fails (the unit/zero-divisor dilemma): each array reaches its target through one filter 1 + a·u with u a free letter. The filter must be invertible to deliver the target, yet clearing the other array's leftovers forces a = ±1, which makes it a zero divisor. This holds for every filter, every letter order and every twist, in shapes where the shared move and the last moves each have a single term. | **Argued by hand, not checked**; one AI skeptic upheld it and proved a stronger statement about the twisted gadget. The two group-algebra facts it uses are checked by computer in `programs/symmetric-2-in-7-dilemma/` (re-run on 19 October 2026). |
| Every repair failed: the 16 order-3 twists of the shared move (residual 1.0), a shared "junk absorber" (floor 1/3), and 407 + 18 twisted symmetric gadgets (all pass the frame-level test and f = 2, and fail at f = 1). | The twists and the absorber: a numerical search found nothing, which is evidence, not proof. The gadgets: **checked by computer**, and the skeptic proved that every gadget of that kind fails at f = 1 (argued by hand, not checked). |
| A second limit skeleton (seed 301): its least residual in an exact f = 2 test, over every coefficient choice and every limit, is 0.2333; a unit needs 0. | **Open as a matter of proof**: a numerical minimisation that the skeptic reproduced (0.233327) from many starts. The track's reason for trusting the minimum was wrong, so this is evidence, not proof. |

The lesson I took: a unit at block 4 needs a second route into each target, related to the first by a genuine
odd-order letter and formed after a shared move. A symmetric program supplies its second piece as a mirror image, and
then the relating letter is forced to be the wrong kind (an order-2 letter, the identity, or a central sign).

### 4. The Fano blueprint and the Singer template

| finding | status |
|---|---|
| Once two arrays share a generator, later steps can separate them only with "frustrated" signs (signs that no consistent sign rule can produce), and separation always divides by 2 (Lemma P). Two literal readings of the Fano blueprint fail for every choice of letters (Lemma S). | **Argued by hand, not checked**; the skeptic of wave 1 checked them and found that one part of Lemma P needs an invertible coefficient. |
| The ends must be shared: at block 4 a private last move wastes 9,184 of a move's 27,456. So 3 in 10 and 4 in 13 can end privately for at most one array, and an exactly tight schedule for none. | **Checked by computer** (exact ranks modulo a prime, given the counting rule); the skeptic held the arithmetic. |
| "The smallest frustrated set spread over frames is the Fano gadget." | **Refuted** by the order-3 pair of section 2. |
| Template searches: 0 of 820 sampled 2-in-7 and 0 of 146 sampled 3-in-10 programs. | Evidence only: the samplers were weak, so these zeros exclude only what was sampled. |
| The Singer-cycle template (2 in 7, 3 in 10, 4 in 13, arrays separated by an odd-order address relabelling) is excluded for every separating matrix: read backwards it has only one move that starts a full-length chain from the inputs, and a unit needs W such moves. | **Argued by hand, not checked**; the skeptic re-derived it and found no error, and its own joint solver said "yes" 0 times in 2,052 calls at 2 in 7 and in 45 calls at 3 in 10. |
| "The number of independent top routes equals the number of shared moves at depth 3." | **Open**: an empirical law in one calibrated family, not proved. |
| The smallest odd-order separators are 1 + sigma and 1 + sigma - sigma^2, both with denominator 2; an order-3 pair formed from the raw inputs is only a change of coordinates. | **Checked by computer** (exact enumeration; re-run on 19 October 2026, `programs/singer-template-separators/`); the second part argued by hand. |

### 5. Closing in on 2 in 7 by shape: Theorem A, shape R and the 1,261 prefixes

The order-3 round found an exact model of the block-4 engine game: its frames are the 2,295 Lagrangian planes of an
8-dimensional binary symplectic space, and a paid move is one step to a neighbouring plane. With it, the private-end
branch of 2 in 7 could be attacked by proof, shape by shape.

| finding | status |
|---|---|
| The model agrees with the real 16 × 16 matrices of the engine (400 of 400 random words; the skeptic, 1,500 of 1,500). | **Checked by computer**; re-run on 19 October 2026 (`programs/dual-polar-model-and-theorem-A/`). |
| Lemma 5+: if exactly two moves touch the two targets, their 2 × 2 coefficient matrix is invertible. Theorem A: if the first three moves are u1 = D b1 y1, u2 = D b2 y2 and u3 = D(r u1 + c r′ u2) with single letters, in any input basis, no 4th and 5th moves deliver both depth-3 targets (the two target moves get a rank-1 matrix). It holds in the engine and widened games. | **Argued by hand, not checked**; the skeptic rebuilt the model by brute force, read every proof line by line, and found no error. |
| "No letter advances one term of a shared move while holding the other back" (an observation of the first order-3 search). | **Refuted**: such a selective advance exists exactly when two particular planes meet, which happens for 43,008 of the 100,352 single-letter shared moves (43 %). The first search's "0 of 7,477" came from how it built its first moves. Checked by computer, recounted by the skeptic, and re-run on 19 October 2026. |
| Theorem R of shape R (Theorem A's start, but the third move also carries raw inputs or extra letters): with one term per frame of the third move, no 2 in 7 exists in the engine game, for every choice of letters, every rational coefficient and any number of extra terms. The exact solve reduces to a finite system with 30 unknowns per generator (Theorem Q). A second letter on u2 can never absorb a selective letter's junk. | **Argued by hand, not checked**; the skeptic read every proof line by line, rebuilt the test from scratch, and tested 1,588 prefixes and 26,869 target pairs of its own: none passed. |
| The shape-R families: 6,820 prefixes and 70,228 target pairs, none passing. | **Checked by computer**, with controls (two private chains pass, Theorem A's prefixes fail); the controls and two of the families re-run on 19 October 2026 (`programs/shape-R-exact-test/`). |
| Shape R with two or more terms on one frame of the third move, and shape R in the widened game. | **Open**: searched by the track and the skeptic, nothing passed, no proof. |
| The 1,261 saved 2-in-7 prefixes of all the order-3 searches: 1,256 are excluded exactly over the rationals, for every end move and every invertible free mixing of the two end moves, so in every input basis. 802 fall to a proved criterion that needs no linear solve (Criterion 0); certificates that the other ingredients are dead exclude 454 more (median certificate 10 junk rows, at most 100). | **Checked by computer**; the skeptic rebuilt both criteria in its own code and reached the same verdict on all 1,261, with the same 5 left over. I re-ran Criterion 0 on all 1,261 on 19 October 2026 (`programs/two-in-seven-prefix-end-tests/`): every value matches. |
| Of the 5 left over, 3 are decided against every end move inside a smaller group of 24,576 letters (bit maps, bit flips, sign flips and the scalar i); 2 (named anyb G3S s61 #11 and w3 STG2 #17 in the data) are excluded only against the "pool" end moves. | The 3: **checked by computer** (one of them by the skeptic). The 2: **open**; their runs against every end move were left running at the end of the round and never collected. |
| The second round's claim that 1,251 of the 1,261 prefixes were ruled out for every end move. | **Refuted** on coverage by its skeptic: on the prefixes that already had a term on K's frame the test tried only the "pool" end moves, and 22 of the claimed exclusions rest on those alone. The corrected count was 1,229 for every end move, 22 against the pool only and 10 open; superseded by the row above. |
| The 3-in-10 analogue of shape R: dimension 0 in all 1,200 natural extensions (not even one array can finish). | **Checked by computer** by one program, not repeated by anyone; no proof for three arrays. |

These verdicts decide given prefixes, not whole shapes. The criterion behind them (a 2 in 7 needs two end-move classes
whose coefficient sums on K's frame are independent) cannot by itself exclude a shape, because a real 2 in 7 passes it
trivially in the basis where its end is private. Still open by proof inside the private-end branch: shape R with two
or more terms on a frame; first moves that mix the arrays or carry several terms; a chain that reaches depth 3 by the
third move; three target moves.

### 6. One or two shared moves early: the one-shared-move shape and Theorem P

Wave 2's lesson asked for a second route into each target, formed after a shared move. Wave 3 put one shared move
early, with two private layers after it (W private starts, one shared move, then an end and a top per array: 3W + 1
moves, so 2 in 7, 3 in 10, 4 in 13). Wave 5 tried two shared moves at layer 2 and one at layer 3.

| finding | status |
|---|---|
| In the one-shared-move shape with single-letter starts and merge, only the main chain reaches the target layer, each output's use of the top moves is forced (Theorem R of that shape), and every two-move path from the merge into a target is a free multiple of one fixed element. So two pieces related by sigma, or by anything else, are only a change of basis of the top moves. Necessary conditions Lemmas 2 to 5 (with single-letter ends, the leftover must sit exactly two layers out, a special position 6 % of random merges have). | **Argued by hand, not checked**; the skeptic checked the proofs line by line, twice. |
| A cheap exact test at f = 2: the widened free letters split the data space into three pieces (16, 120 and 120 dimensions), and in each piece "the program reaches the target" is a row-space test. It sees the zero-divisor failure that the frame-level test misses (the wave-2 near miss fails it at 0.267) and runs in milliseconds. | **Argued by hand, not checked** (a valid necessary condition); upheld by the skeptic. Program: `programs/one-shared-move-exact-decisions/rep2.py`. |
| 755 small skeletons built in the special positions (bare, with the landing injections the lemmas ask for, singly, cumulatively and in pairs, and with sigma in the merge using both kernel routes) have no solution for any complex coefficients: Groebner bases over the rationals in Singular return {1}. Controls behave: the plain 8-move program is solvable, the plain program minus one move is not. | **Checked by computer**; the skeptic re-ran every exact decision, also modulo a prime, and rebuilt the pipeline in its own code. On 19 October 2026 I regenerated the 25 systems of the first batch on the Mac; all are byte-for-byte identical to the decided ones (Singular itself not re-run). |
| 144 further skeletons were too large and are undecided (133 counted by the track, 11 more found by the skeptic). | **Open.** A later round claimed to exclude 76 of them with certificates of degree 2 or 3; its skeptic found that those certificates hold modulo one prime only, so the 76 remain undecided (see `../theory-of-theories-and-logic/`). |
| Theorem P: private starts, private layer-2 moves, one shared layer-3 move and private tops never make a unit, for any number of arrays W ≥ 2, with single letters and any injections, in the engine game. So no 2 in 7, 3 in 10 or 4 in 13 of this form exists. | **Argued by hand, not checked**; the skeptic read it line by line and checked its four geometric facts and Lemma D4 exhaustively over all 2,295 frames with its own code, under weaker hypotheses than the proof: 0 survivors (1,344 with one condition dropped, so the check has teeth). That check re-run on 19 October 2026 (`programs/theorem-P-frame-check/`). |
| The variant with shared layer-2 moves (including a sigma-pair): all 68 "cross" cores and 34 sigma-pair cores tested are excluded at the frame level over every engine letter, and 395 + 598 + 18 + 18 skeletons have no solution by Groebner bases. There is no theorem for this variant. | **Checked by computer** for the cores and skeletons tested; the skeptic upheld the core verdicts. |
| Lemma D′ as stated (the first-round lemma in the cross form). | **Refuted** by the skeptic: it misses two ways to cancel some leftovers, and 49,152 letters per such leftover pass a condition that the proof says none can. The core verdicts do not depend on it. |
| "The 99 + 9 skeletons that Singular left undecided are excluded too." | **Open**: not confirmed. Some of them contain two-letter injections outside the model the exclusion used; the skeptic's job that would say which was left running and never collected. |
| Extending the first round's private-ending count to the endings with leftovers that Theorem R forces (which would close 3 in 10 and 4 in 13 in the one-shared-move shape). | **Open**: argued not to extend by itself; the computation needed (exact ranks of dense matrices up to 205,920 × 49,104) is out of reach as posed. |

### 7. Shared ends for 3 in 10 and 4 in 13

3 in 10 and 4 in 13 must end on genuinely shared moves (section 4), and their ends must be "rationed" (no two pieces
extend to a basis).

| finding | status |
|---|---|
| Over an abelian group of ratio letters, a rationed shared end for W arrays needs at least W fields in the rational group algebra. So the order-3 shift alone, or any letter of prime order, never works; the order-4 bit cycle works for 3 arrays only; 4 arrays need two signs, order 6, order 8 or more. If every piece uses the same paid end move, W generic mixtures of the pieces are W private ends in some basis. | **Argued by hand, not checked** for abelian groups; the class counts are **checked by computer** and the skeptic reproduced them with its own code (re-run on 19 October 2026, `programs/shared-ends-3-in-10-4-in-13/`). Non-abelian ratio groups were not classified. |
| "In the coordinate world all end moves form one class", which would have closed it for every abelian ratio group. | **Refuted** by the skeptic: with sign letters the four coordinate end moves form 4 classes, with the order-3 group 2. So the remaining window is wider than first stated (re-run on 19 October 2026). |
| Searches: 1,447 prefixes for 3 in 10 and 678 for 4 in 13, across 10 templates; none passes even the relaxed level. The best 3-array prefix carries only what two private chains carry; for 4 arrays only the coordinate world does better (three arrays' worth out of four). | Evidence only (one template hit its size cap on up to 25 of 57 prefixes in one job). |

### 8. The symmetric targets 7 in 23 and 8 in 26, and 7 in 24

Each is three orbits of W moves plus two shared moves (3 × 7 + 2 = 23, 3 × 8 + 2 = 26), with a symmetry of order 7 (from
a Singer cycle) or a group of order 8. The orders below list the moves from first to last.

| finding | status |
|---|---|
| The blueprint as written (both shared moves first), and both shared moves last, are impossible: read backwards a unit needs W independent routes out of the inputs, and only 2 of the 7 exist. The test splits into W independent characters. A single-term shared start (F O O O F) fails, and so does a single-term shared move right after the first orbit (O F O O F). | **Argued by hand, not checked**; the skeptic upheld these and checked O F O O F exhaustively in the engine game. |
| "Theorem R′": the same holds one step deeper (order O O F O F). | **Refuted** as written by the skeptic, which found 3,024 of 286,720 configurations escaping it. Wave 5 then showed that every one of these escapes is spurious (the term the middle shared move reads must be exactly 2 moves from the target, and in all 3,024 it is 3 or 4): with that constraint 0 of 35,840 configurations escape. |
| O O F O F with a single-term middle shared move is excluded for 7 in 23 and 8 in 26, for every way of building that shared move (a Sylow argument for 7 in 23; all 11,391 relevant subgroups listed for 8 in 26); 1,216 skeletons of the other open orders fail every character tested; two single-term shared moves on one frame serve no character (0 of 28,311,552). | **Checked by computer** by one program, not repeated by anyone (the Sylow step for 7 in 23 is argued by hand): the skeptic confirmed only the file counts. Of the 1,216 skeletons, 520 used random coefficients, which is a screen only; for pairs of shared moves with different characters 6 of 8 parts were collected. |
| A 7-in-23 program of order O O F O F whose middle shared move reads one term through 184 letters at once (a "group-algebra" content) passes the frame-level test in all 6 characters outside its own. | **Checked by computer** and reproduced by the skeptic with independent code, but the skeptic showed the pass is generic (180 equations, 208 unknowns, full row rank), so it carries little evidence towards a unit. Re-run on 19 October 2026 (`programs/seven-in-23-and-8-in-26-blueprints/skeptic-wave5/own.py`, mode z7). |
| 8 in 26 with a group-algebra middle shared move is excluded for every twist subgroup: whatever its parts, the product the top orbit receives stays in a fixed space, and the target lies outside it in at least 4 of the 7 characters tested, where at most 2 can be set aside. Calibration: 7 in 22, which the hook already excludes, fails structurally at a fixed distance 0.456. | **Checked by computer**; the track tested 3,264 of the 11,391 twist subgroups without saying why, and the skeptic tested all 11,391 (they all fail). The later verdict is the one used. I re-ran the skeptic's test over all 11,391 subgroups on 19 October 2026 (8 minutes on the Mac): every one fails in at least 4 characters, exactly as recorded. |
| The 7-in-23 construction is excluded in both orders, in its configuration, by rank counts in the 14-dimensional representation gl14 of the address maps. O O F O F: in the first shared move's own character the top reaches rank at most 11 where 12 is needed (the constituent (14,2)). O O F F O: once only the four slots that other terms can reach are free, each shared move's own character caps the rank at 13 where 14 is needed (the constituent (6,5)). Requiring one content for all seven characters adds no constraint. | **Checked by computer**: exact bounds from linear algebra in floating point, with rank gaps of 13 orders of magnitude or more; ranks 11 and 13 are reached. Not refuted: the skeptic rebuilt every space on both lifts with its own code, checked the decomposition on 199 fresh letters, and closed a route the track left open (all order-7 elements are conjugate, so the O O F F O exclusion holds for every pair of order-7 twists with different labels). I re-ran both bounds on 19 October 2026 (`programs/seven-in-23-gl14-rank-bounds/`). |
| An earlier round's claim that this construction was already excluded (order O O F O F). | **Refuted** in that round: the exclusion applied only to a generic lift, random contents fail in the same way, and a corrected lift kept the family alive; told in `../existence-routes/`. Wave 8's exclusion above replaces it. |
| Scope of the exclusion: the configuration as built, in which both shared moves read only the second orbit's term. Not covered: other cells, the orders O F O F O and O F F O O, orders that are not chains, and a variant in which the second shared move reads the first. | **Open.** |
| The nearest variant, O O F F O F (7 in 24, saving 0.111): O O F F O misses by one rank, only in gl14 and only in its two own characters; in the other representations the top already gives the full block. A terminal shared move twisted by an order-7 element adds up to 2 ranks, so 7 in 24 passes every distance-4 test computed. | **Open.** The top's part is computed; the terminal move's part is argued by hand, not checked. Nothing about distance 3 and below has been computed. |

### 9. Block 5: the butterfly 2 in 7 and the 2-in-8 shape

At block 5 a 2 in 7 (0.222) sits exactly on the two-layer bound W + m. Wave 2 showed that both inputs must enter within
the first two shared moves and that the top needs a special merge. The shapes tried: the *butterfly* (private starts, a
shared trunk of three moves, private ends), a trunk whose third move has several branches, and the *split trunk* in
which two branches serve different arrays. A split-trunk 2 in 7 gives a 2 in 8 (0.139) of the shape "two starts, two
shared moves, two branch moves, two ends", so the 2-in-8 shape became a target in its own right.

| finding | status |
|---|---|
| Corollary B++: subtracting one array's part from the other's leaves a one-array program that shared moves cannot touch, so both inputs must be injected within the first two shared moves; the ends split cleanly only for special merges. | **Argued by hand, not checked**; the skeptic checked it. |
| "Leftovers find no partner, so the chain of cancellations stalls." | **Refuted** by the skeptic as a sampling artefact: by the exact rule, every leftover at layers 3 and 4 has a partner. |
| Theorem W3: the butterfly is never a unit when the starts, the merge and the trunk use single letters, whatever the input injections into the first two trunk moves, the side terms into the last one and the end contents. The proof rests on an exact census (393,216 configurations per trunk path; 3,360 paths) and the producer rule: a leftover cancels only against a term at exactly the same group element, and the only move that produces it is the one that also carries that array's target. | **Argued by hand, not checked**; upheld by the skeptic, which redid the census over all 3,360 paths with code written from scratch (0 violations of the fact the proof uses over 1,720,320 pairs). Both censuses re-run on 19 October 2026 (`programs/block5-butterfly-census/`). Not covered: multi-term trunk, merge or starts; side terms from the starts into the second trunk move; the other shapes E5, E7, E12 and H. |
| The trunk with several branches: every column that delivers one array's target also carries the other's leftover with the same factor unless a "companion" term exists; companions can sit only within 3 steps of the other array's anchor; each target touches exactly three branch frames and only one of them is close; a split-trunk 2 in 7 gives a 2 in 8; Lemma Q for injections into the first trunk move. | **Argued by hand, not checked**; upheld by the skeptic, which turned the three-frames count into a proof and checked it over every anchor pair (re-run on 19 October 2026, `programs/block5-two-branch-census/`). |
| Lemma Z′ for three or more branch letters, and "a second branch in the same frame changes nothing". | **Refuted** by the skeptic, which built configurations that escape (none of them a unit). Lemma Q for injections into the second trunk move is not proved. |
| Wave 7: collisions between injection streams depend only on their relative position (Lemma E1-σ); in the most separated layout they are forced (Lemma T); if every producer of a pollution row lies in one class family, the readout is not invertible (Corollary O). | **Argued by hand, not checked**; confirmed line by line by the skeptic. |
| "Only orbits that swap the two branch frames can help", and one bullet of Lemma S in the 2-in-8. | **Refuted** by the skeptic: cycle orbits help too (one reached a layer-3 residual of 0.047), and in the 2-in-8 a two-letter start whose copy moves the target frame passes layers 5 and 4 exactly. The track's best closure (0.432) cannot be reproduced, because the forced-zero tool's result depended on Python's hash order. |
| The forced-zero tool used from wave 4 on had two gaps: a column forced to zero before it was built stayed listed as a producer, so later forced zeros were missed; and a row met while all its built producers were dead was dropped for good. It never produced a false exclusion, but it kept dead designs looking alive. A corrected tool fixes both and audits every closure. | **Checked by computer**; the skeptic upheld it and verified all 19 forced-zero chains of wave 8 with its own brute-force code. Programs: `programs/block5-sparse-designs/peel2.py` (and the old `peel.py`). |
| Along each orbit the layer-4 coefficients must form a geometric sequence whose ratio is a root of unity, so layer-3 tests at random coefficients mean nothing. | **Checked by computer** (exact computation of the layer-4 solution set). |
| All 20 rebuildable layer-3 near misses of wave 7 and its skeptic (old residuals 0.047 to 2.0) fail layer 3 exactly, for every choice of injection coefficients with every injection on: 19 by forced-zero chains, the best (s1302) by an exact rank test (168 against 169). The five designs whose closures had outgrown the solver (s801 to s805) fail within two steps of the seed: a ball of 22 to 42 rows is already inconsistent for every value of the free parameters. The new union designs s1401 and s1402 fail by Groebner bases (msolve and Singular, over the rationals and two primes). | **Checked by computer**; one skeptic, which reproduced s1401 and s1402 at a new prime and over the rationals, checked by brute force 161 and 128 kill chains, and ran its own ball test on all of s801 to s805. On 19 October 2026 I regenerated the 18 radius-2 systems of s801, s1401 and s1402 from the skeptic's data; all are byte-for-byte identical to the decided ones (`programs/block5-ball-test/`). Not re-checked by anyone: s1302's rank test, s1503, and s1407 to s1409. |
| Lemma F: in the 2-in-8 with the most separated layout, a column through a far branch dies on its own layer-4 row unless some stream has a branch term at one prescribed frame. | **Argued by hand, not checked**; no skeptic has read it. |
| The survivor s1501 (2 in 8, saving 0.139; a swap orbit plus two extra stream orbits, 32 injections) passes layer 4 exactly. Its two-step ball at layer 3 is consistent numerically with every coefficient nonzero (residual 5e-29); at three steps every numerical run that gets near zero switches an injection off. Groebner programs gave no verdict in 20 minutes. | **Open.** Layer 4 exact (the skeptic); everything else numerical only. |

The two-step and three-step ball systems of s1501 are published in `programs/block5-ball-test/systems/` so that anyone
can attempt the exact verdict.

## What failed, in one list

These are the claims of this section that broke, with what broke them. Each is told in context above.

- The first round's hook values at odd k: refuted by ranking all 16 character blocks (the conclusion stands).
- One closing sentence of the near-miss note: refuted by a later certificate search.
- "The smallest frustration is a Fano gadget": refuted by the order-3 two-piece separation.
- "Sign patterns need four pieces": refuted (multi-term pieces reach three).
- The first order-3 searches' completeness, and the second round's coverage claim: refuted on coverage; superseded by
  the exact decision of 1,256 of the 1,261 prefixes.
- "No selective advance": refuted (43 % of single-letter shared moves have one).
- The symmetric 2 in 7 as an approximate unit: it is a small exact program that fails the full test, for a reason
  argued by hand and upheld by a skeptic (the unit/zero-divisor dilemma).
- "The shared ends' coordinate end moves form one class": refuted (4 classes with signs, 2 with order 3).
- Theorem R′ one step deeper (O O F O F): refuted as written, then its 3,024 escapes shown spurious.
- Lemma D′ as stated: refuted (two cancellation routes missed).
- "Leftovers find no partner" at block 5: refuted (a sampling artefact).
- Lemma Z′ for three or more branch letters: refuted.
- "Only swap orbits help" and a bullet of Lemma S: refuted (cycle orbits; a two-letter start in the 2-in-8).
- The old forced-zero tool: two gaps, which kept dead block-5 designs looking alive; wave 8 decided every old near
  miss again, exactly.
- Wave 7's best layer-3 closure (0.432): not reproducible (it depended on Python's hash order).

## What is still open, and the first experiment for each

| open question | worth | first experiment |
|---|---|---|
| The 640 room of the private-end 2 in 7 at block 4. | 0.096 | Search for a certificate with the full spectral factor 7/4 at depth 3 whose angle factors satisfy s3 < 0.7·s1 (s3 and s1 are the average sin² per angle at depths 3 and 1; this is the corrected closing test), starting from the hook's relatives and the orbit cuts in the pieces of dimension 85, 238, 510 and 595. Or confine an exact search to what a later round says a 2 in 7 must use: address maps in the middle of the program, or phase letters between moves (`../reed-muller-and-one-move/`). |
| The 2 prefixes still open against end moves outside the pool. | (closes the sampled search) | `poolsuff.py p32.jsonl out.jsonl anyb_G3S_s61:11,w3_STG2:17` in `programs/two-in-seven-prefix-end-tests/`, run to the end; the pass criterion is a line "pool-suffices" per prefix. |
| Shape R with two or more terms on one frame; shape R in the widened game. | 0.096 | Extend Theorem R's type-by-type argument to rows that take two contributions (where the zero divisors of section 3 can enter), and test it with `theoremR_check.py` on the families COS and MIX, which already contain such prefixes. |
| The other shapes of the private-end branch: first moves that mix the arrays or carry several terms; a chain at depth 3 by the third move; three target moves. | 0.096 | Apply shape R's small system and Lemma G to each; begin with u3 and u4 as the target moves. |
| The 144 undecided one-shared-move skeletons, and the 99 + 9 undecided two-shared skeletons with two-letter injections. | 0.096 to 0.150 | Decide them exactly with Singular on a high-memory machine, one at a time, starting with the smallest; for the two-shared ones, first list which skeletons contain two-letter injections. |
| 3 in 10 and 4 in 13 with shared ends. | 0.132, 0.150 | Settle the coordinate world exactly with sign letters and several end classes (a finite group of order 6,144); then classify non-abelian ratio groups. |
| 7 in 24 at block 4 (O O F F O F). | 0.111 | The structured top: a joint search over the top content and the shared move's content in gl14, then the conditions at distance 3 and below, then the finite-field pipeline. Screen the other orders and configurations of 7 in 23 with the exact rank bounds of `pair.py` and `one65.py` first; each takes minutes at most. |
| s1501, the one undecided block-5 2-in-8 design. | 0.139 | Split the layer-4 solution set by which of the 8 Fourier characters each cycle uses, work modulo a prime p ≡ 1 (mod 8), and run the one-shot test at three steps on each piece with Singular (the systems are in `programs/block5-ball-test/systems/`). |
| Block-5 2 in 7 beyond the butterfly and the split trunk: multi-term trunks, merges or starts; the shapes E5, E7, E12 and H; the 2-in-8 cycle orbits' readout invertibility; the z-cover designs. | 0.222, 0.139 | Run the z-cover list (`z.list`) through `batch.sh` with the corrected forced-zero tool, then test the readout of every surviving cycle orbit for invertibility. |
| Results relied on that no skeptic has checked: the first wave-5 session's single-term exclusion and 1,216 skeletons; Lemma F; s1302's rank test; s1503 and s1407 to s1409. | - | A first skeptic for each. They are cheap, and later results rest on them. |

## How this fits the rest of the research

Later rounds narrowed the 2-in-7 question from other directions. With address maps and phase letters used only at the
ends, there is no 2 in 7, 3 in 10, 3 in 11 or 4 in 13 at block 4 (checked by computer, exhaustive tables); without
phases 2 in 7 is also excluded for most small cyclic address groups (checked by computer, one skeptic); there is no
single-term 2 in 7 in the block-4 cube (checked by computer, exhaustively, one skeptic); and any two-array saver must
split one input and pack in every gauge (argued by hand, not checked, one skeptic). A proof of Proposition 6 would
close every two-array target, 2 in 7 at blocks 4 and 5 and 2 in 8 at block 5 included; a proof of Conjecture T would
close their single-term branch. Both are open. These are told in `../block-bounds-and-new-rules/`,
`../frame-packing-two-arrays/`, `../reed-muller-and-one-move/` and `../theory-of-theories-and-logic/`.

## Programs and data

Every folder below has its own README with the command, the time and memory it took, and the expected output.
Python 3.9 or later with numpy unless the README says otherwise.

| folder | what it holds | re-run for this revision |
|---|---|---|
| `programs/hook-ratio-640/` | the hook values 27,456 and 68,320 (exact ranks modulo a prime) | yes |
| `programs/order-3-two-piece-separation/` | the two-piece order-3 separation, at blocks 3 and 4 | yes |
| `programs/order-3-halving-end-control/` | the 8-move program with the order-3 halving end | yes |
| `programs/symmetric-2-in-7-dilemma/` | the two group-algebra facts behind the dilemma | yes |
| `programs/singer-template-separators/` | the smallest odd-order separators | yes |
| `programs/dual-polar-model-and-theorem-A/` | the exact block-4 model, its check against the matrices, the selective-advance census | yes |
| `programs/shape-R-exact-test/` | shape R's exact test, the families, the Theorem R check | controls and two families |
| `programs/two-in-seven-prefix-end-tests/` | the 1,261 saved prefixes and the end-move criteria | Criterion 0 on all 1,261 |
| `programs/shared-ends-3-in-10-4-in-13/` | the classification of shared ends over abelian ratio groups; the skeptic's class count | yes |
| `programs/one-shared-move-exact-decisions/` | the f = 2 test, the skeleton builders, the Singular systems and verdicts | systems regenerated; Singular not re-run |
| `programs/theorem-P-frame-check/` | the skeptic's exhaustive check of Theorem P's geometry | yes |
| `programs/seven-in-23-and-8-in-26-blueprints/` | the blueprint searches of waves 3 and 5 and both skeptics' code | the skeptic's z7 and z2 modes |
| `programs/seven-in-23-gl14-rank-bounds/` | the gl14 decompositions and the rank bounds | the two exclusions |
| `programs/block5-butterfly-census/` | the census behind Theorem W3 and the skeptic's geometry check | yes |
| `programs/block5-two-branch-census/` | the three-middle-frames census | yes |
| `programs/block5-sparse-designs/` | the block-5 design code, the corrected forced-zero tool, summaries of the decisive runs | no (needs python-flint, Singular, msolve) |
| `programs/block5-ball-test/` | the decisive radius-2 ball systems, their verdicts, and the open s1501 systems | systems regenerated; solvers not re-run |
| `programs/track-code/` | the code of the remaining tracks, as it ran | no |

`results.json` lists every finding of this section with its status, its figures and its evidence files.
`data/README.md` lists the data files with their sha256, and the unpublished large files with how to rebuild them.

Two warnings apply to every program here that writes systems for a Groebner-basis program: msolve 0.6.5 reads a
product written `x*x` as `x` (write powers as `x^k`, and cross-check every msolve verdict with Singular), and the
helper `sing2ms.py` corrupts any entry that is a sum of two terms. Neither affected a verdict of this section: no
system here contains a power or such an entry, and every msolve verdict was cross-checked with Singular.

Large outputs are not published: the order-3 round's folders (79 MB in all), the gl14 track's run folders (63 MB) and
the last block-5 skeptic's job folders (58 MB) are mostly search logs and solver dumps, and the programs above
regenerate what matters.
The layered block-5 systems are Python pickles and are not published, because loading a pickle can run code; the
systems that were decided are published as plain text instead.

## Credits

I did this work with a team of AI agents (Claude) that I directed: I set the goals, chose which ideas to pursue and
decided at each stage what to do next; the agents ran the searches, wrote the proofs and programs, and tried to refute
one another's work. The exact decisions use Singular (W. Decker, G.-M. Greuel, G. Pfister and H. Schönemann) and msolve
(J. Berthomieu, C. Eder and M. Safey El Din); the finite-field arithmetic of the block-5 code uses python-flint (Python bindings to the FLINT library); the
geometry relies on standard facts about dual polar graphs and on Sylow's theorems. No outside construction is used in
this section; the engine, the games and the certificates are those of this repository and its earlier revisions. The
cost model in which the paid moves are counted is OpenAI's (the Lean file `UniformFourier.lean` of openai/math and
OpenAI's two preprints).

## Sources

All of this section comes from my working notes of 14 and 15 October 2026, which are kept in a private archive;
paths below are relative to its `checks/wht43/moonshot/` folder. Archive catalogue entries are given in brackets.

- Round notes: `NOTE.md` (wave 1), `NOTE-wave2.md`, `NOTE-wave3.md`, `NOTE-wave4.md`, `NOTE-wave5.md`,
  `NOTE-wave7.md` (block-5 part), `NOTE-wave8.md` (w8-b5sparse and w8-bpgl14 parts); `o3/NOTE.md` (rounds 1 to 3).
- Section 1 [N1, C19]: `ms-fano/` and `v-ms-fano/`; `o3/o3-ratio/` (PROOF.md Lemmas 1 to 5; chirank.py, blkk.py,
  rankp.c → `programs/hook-ratio-640/`) and `o3/v-o3-ratio/`; the closing-sentence correction and the certificate search
  are in `../meta2/` (follow-ups, cert-07); the exact value of 68,320 in `../meta2/v-a2-vi-stability-compute/`
  (STATUS.txt, Q5). The block-4 bound 41/13 and its four programs: catalogue entry [L9].
- Section 2 [N10, X21]: `v-ms-fano/cex.py` (→ `programs/order-3-two-piece-separation/`); `o3/o3-build/` (PROOF.md
  Lemmas 1 to 6; verify_ctl.py → `programs/order-3-halving-end-control/`) and `o3/v-o3-build/RUNLOG.txt`.
- Section 3 [X16]: `mw-sym4/` (→ `programs/track-code/symmetric-block4-skeletons/`), `mw2-absorb4/` (NOTE.md sections 4
  to 9; → `programs/track-code/symmetric-2-in-7-repairs/`, cyclelemma.py → `programs/symmetric-2-in-7-dilemma/`, whose
  output was first saved by the re-run of 19 October 2026),
  `v-mw2-absorb4/`; the seed-301 limit in `w3-asym/NOTE.md` section 4 and `v-w3-asym/`.
- Section 4 [X15, X17]: `ms-fano/PROOF.md`, `ms-fano/NOTE.md` (→ `programs/track-code/fano-blueprint/`);
  `mw2-singer4/` (PROOF.md, sep3.py → `programs/singer-template-separators/`; the rest →
  `programs/track-code/singer-template/`) and `v-mw2-singer4/`.
- Section 5 [X18, X19, X20, C18]: `o3/o3-close/` (PROOF.md: the model, Proposition S, Lemmas 4′ and 5+, Theorem A) and
  `o3/v-o3-close/` (→ `programs/dual-polar-model-and-theorem-A/`); `o3/o3-shapeR/` and `o3/v-o3-shapeR/` (→
  `programs/shape-R-exact-test/`); `o3/o3-finish/`, `o3/v-o3-finish/`, `o3/o3-finish2/` (work/, PROOF.md: Lemmas A, S,
  J, P, Criterion 0, Lemmas D1 to D4) and `o3/v-o3-finish2/` (→ `programs/two-in-seven-prefix-end-tests/`).
- Section 6 [X23, X24]: `w3-asym/` (NOTE.md sections 1 to 4; → `programs/one-shared-move-exact-decisions/`) and
  `v-w3-asym/RUNLOG.txt`; `w5-twoshared/` (NOTE.md section 3; → `programs/track-code/two-shared-layer-2-moves/`) and
  `v-w5-twoshared/` (myframes.py → `programs/theorem-P-frame-check/`); the 76 skeletons whose low-degree
  certificates held modulo one prime only in `../meta3/ml-fooling/` and `../meta3/v-ml-fooling/`.
- Section 7 [X22]: `o3/o3-shared/` and `o3/v-o3-shared/` (→ `programs/shared-ends-3-in-10-4-in-13/`).
- Section 8 [X28, X29, X30, X40, N4]: `w3-blueprint/` and `v-w3-blueprint/`; `w5-blueprint/` (NOTE.md sections 2 to 9)
  and `v-w5-blueprint/` (→ `programs/seven-in-23-and-8-in-26-blueprints/`); `w8-bpgl14/` and `v-w8-bpgl14/` (→
  `programs/seven-in-23-gl14-rank-bounds/`); the earlier exclusion claim in `../exist3/e3-bpexact/` and its skeptic.
- Section 9 [X31, X32, X33, X41, N9, T12]: `mw2-b5beyond/` and `v-mw2-b5beyond/` (→
  `programs/track-code/block5-early-separation/`); `w3-b5/` (PROOF.md, Theorem W3) and `v-w3-b5/` (→
  `programs/block5-butterfly-census/`); `b5-twobranch/` and `v-b5-twobranch/` (→ `programs/block5-two-branch-census/` and
  `programs/track-code/block5-two-branch-trunk/`); `w7-b5/` and `v-w7-b5/` (→
  `programs/track-code/block5-split-trunk-orbits/`); `w8-b5sparse/` (→ `programs/block5-sparse-designs/`) and
  `v-w8-b5sparse/` (work/vball.py and the ball systems → `programs/block5-ball-test/`).
- Cross-references: archive catalogue entries L4, L24, L26, L27, X39, X43, N12, N13 and the rounds tier2, rmfam,
  prop6, packing.
- The re-runs of 19 October 2026 were done on an Apple M4 Mac mini; their outputs are the `rerun_*` files named in the
  program READMEs.
