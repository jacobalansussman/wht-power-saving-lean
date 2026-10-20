# What our methods can and cannot do: the abstract rounds

This section covers three research rounds of 14 and 15 October 2026. They did not hunt for a unit. They asked a
question one level up: which kinds of lower-bound proof, and which kinds of construction, can exist at all, in which
game, and how far can each kind ever go? The first round I called a "theory of theories": every abstract claim had to
predict a specific calculation, written down before the calculation ran. A second abstract wave and its follow-ups
took up spectrum duality, the Sylvester–Gallai theorem, representation stability, geometric complexity theory (GCT),
ZX rewriting and numerical algebraic geometry (homotopy continuation and exact Gröbner bases). A third round asked
what mathematical logic says: can one prove that the proof I want exists?

I did this work with a team of AI agents (Claude) that I directed. I set the goals, chose which ideas to pursue and
decided at each stage what to do next. The agents ran the searches, wrote the proofs and programs, and tried to
refute one another's work. Nothing in this section is checked in Lean. No person has read any of these proofs. Where
I say "a skeptic", I mean one AI agent that tried to break a result with its own reasoning or its own code. It is not
a human referee. Most results here had one skeptic; a few had two, which is what the rounds' own rule asked for.

**In one paragraph.** The best possible saving is a single number, whatever block sizes are used, and every
lower-bound method of our kind meets it through one "growth exponent" per block. That gives a screen: to rule out a
saving of 0.1, a family of certificates must keep proving more than m^0.9 moves per array at arbitrarily large
blocks m. Cuts, colourings and hooks cannot (hand proofs that AI skeptics re-derived), and coordinate sets fall short
wherever they have been computed.
Lifts are the one named family that might,
and whether they can is open; the evidence leans against them. Every certificate's value is a spectral factor, which
is below 2 at every block up to 40, times an angle factor, so everything depends on how far the angle factor can grow. On the way, a
finite-group computation found a certificate in a 105-dimensional piece of the group Sp(6,2) that rules out four
block-3 "mixed" units I had thought were my best targets above the bar. On the construction side, 16 of the 17
symmetric block-4 programs that pass a cheap test have no unit, and one is open. Logic gives no shortcut: "some unit
saves more than 0.1" is the kind of statement that is provable exactly when it is true, so any existence proof must
work on the unit equations themselves. At one fixed shape exactly one of two finite certificates exists, a refutation
or a matrix model of the unit, so no fixed shape can be independent of ordinary arithmetic. Independence could live
only in the unbounded ceilings, and only if the best saving is not a computable number. No unit was found in any of
these rounds, and none was expected from them.

## Contents

- [The question](#the-question)
- [Words used here](#words-used-here)
- [What I tried](#what-i-tried)
- [What I found](#what-i-found)
  1. [The saving is one number, and every method meets it through one exponent](#1-the-saving-is-one-number-and-every-method-meets-it-through-one-exponent)
  2. [The ladder of games: what each game lets a certificate be](#2-the-ladder-of-games-what-each-game-lets-a-certificate-be)
  3. [How far each kind of certificate can ever go](#3-how-far-each-kind-of-certificate-can-ever-go)
  4. [A certificate from the group Sp(6,2) rules out four block-3 mixed
  units](#4-a-certificate-from-the-group-sp62-rules-out-four-block-3-mixed-units)
  5. [A Venn certificate with a spare row](#5-a-venn-certificate-with-a-spare-row)
  6. [Lifts: dead at block 3, undecided at large blocks](#6-lifts-dead-at-block-3-undecided-at-large-blocks)
  7. [The angle law, and the block-4 near miss as one number](#7-the-angle-law-and-the-block-4-near-miss-as-one-number)
  8. [The chi game: every certificate stops at 169/57](#8-the-chi-game-every-certificate-stops-at-16957)
  9. [Cohomology: the hook bounds carry over to Clifford chain
  maps](#9-cohomology-the-hook-bounds-carry-over-to-clifford-chain-maps)
  10. [Approximate units](#10-approximate-units)
  11. [One array needs m moves in every game; ZX rewriting](#11-one-array-needs-m-moves-in-every-game-zx-rewriting)
  12. [Geometric complexity theory](#12-geometric-complexity-theory)
  13. [Exact algebra and homotopy continuation on symmetric block-4
  programs](#13-exact-algebra-and-homotopy-continuation-on-symmetric-block-4-programs)
  14. [Logic, part 1: what kind of statement the question is](#14-logic-part-1-what-kind-of-statement-the-question-is)
  15. [Logic, part 2: when can an impossibility proof of our kind
  exist?](#15-logic-part-2-when-can-an-impossibility-proof-of-our-kind-exist)
  16. [Logic, part 3: proving existence by showing that every refutation
  fails](#16-logic-part-3-proving-existence-by-showing-that-every-refutation-fails)
  17. [The general pattern: how to show that no proof of a given kind
  exists](#17-the-general-pattern-how-to-show-that-no-proof-of-a-given-kind-exists)
  18. [Claims refuted along the way](#18-claims-refuted-along-the-way)
  19. [Statements of these rounds that later work changed](#19-statements-of-these-rounds-that-later-work-changed)
- [What is still open, with a first experiment for each](#what-is-still-open-with-a-first-experiment-for-each)
- [Programs, data and notes](#programs-data-and-notes)
- [Credits](#credits)
- [Sources](#sources)

## The question

A *unit* "W in M" at *block* m is a fixed program that applies the m-direction transform to W arrays using M paid
moves. The plain schedule needs W·m. A unit with fewer moves gives the saving ln(W·m/M)/ln m to the whole algorithm.
The bar I set for any new route is a saving of 0.1; at block m that means at most m^0.9 moves per array (2.69 at block
3, 3.48 at block 4, 4.26 at block 5). Since 15 October a saving of 0.01 also counts, as a stepping stone.

Our lower bounds all work the same way. A *certificate* is a quantity that the free steps never change and that each
paid move raises by at most a fixed *cap*. If each array needs a *demand* of it, then M ≥ W·demand/cap (Theorem R of
the earlier rounds). Our best family, the *hook*, proves 14/5 = 2.8 moves per array at block 3 and 41/13 = 3.154 at
block 4.

For matrix multiplication the path to today's knowledge went: small calculations (Strassen's seven multiplications,
Bini's approximate schemes), then abstraction (Schönhage's τ-theorem, Strassen's asymptotic spectra), then barrier
theorems that rule out whole methods (Ambainis–Filmus–Le Gall, Alman, Christandl–Vrana–Zuiddam), then calculations
redirected by the barriers. Our saving is also a rate that improves when units are combined, so I asked the agents to
build the same machinery for our problem and to let it predict. The concrete questions were:

1. Is the best saving one number, or can different block sizes give different answers? What does each certificate
   see of it?
2. Which certificates can exist in which game, and what is the most each kind can ever prove?
3. Can any known kind of certificate settle the bar at all?
4. What must any saving unit look like?
5. Does spectrum duality, Sylvester–Gallai, representation stability, GCT, ZX rewriting, or numerical algebraic
   geometry open a route or close one?
6. What does logic say? Can existence of a unit, or of a proof, be shown without finding the object?

## Words used here

- **Array, block m, digit columns.** An array is one vector being transformed. A unit at block m works on m address
  bits of each array. Inside the recursion an address splits into f digit columns, and a unit must be right for every
  f.
- **Paid move, free step.** A paid move applies one direction's 2-point transform to one array; it becomes a recursive
  call, and it is what I count. A free step costs a fixed amount per entry: adding arrays, phases, relabelling
  addresses.
- **Games.** Each game says which steps are free.
  - The *engine game* has the free steps of our Lean engine: Pauli flips, affine address maps, Clifford phases, and
    rational mixing of arrays.
  - *Unimodular chain maps* add XORs from one digit column into the next that are invertible as polynomials.
  - The *chi game* adds the running XOR (each column XORed with the one before it).
  - The *widened game* adds every address reshuffle, Toffoli included, and every pattern of i-phases.
  - The *everything-free game* has chain maps together with Toffoli or controlled-S (a phase i on the input 11).
- **Moves per array, ratio, ceiling.** A certificate's *ratio* is the number of moves per array it proves. The
  *ceiling* of a kind of certificate is the most any member can ever prove.
- **Frame, distance, edge.** A frame is a maximal set of compatible Pauli flips. There are 135 frames at block 3 and
  2,295 at block 4. Two frames are at distance d when they differ in d dimensions; frames at distance 1 form an edge.
  An *interval* is a pair of opposite frames with every frame between them.
- **Irreducible piece, cut, label.** Under the free steps a certificate's space splits into smallest pieces that the
  symmetry group cannot break further. A *cut* is built from whole pieces that do not repeat. A *label* combines
  repeated copies of a piece in a new way. Only labels can beat every cut.
- **Lift, excess.** Modulo 2 there is a perfect certificate made from the frames. A lift copies it over the rationals:
  each frame becomes a rational subspace that reduces to it modulo 2. Neighbouring frames then overlap less than they
  should, and the total shortfall over all edges is the *excess* E. A lift at block m proves m/(1 + E/#edges) moves
  per array. It is the Fano plane, which cannot be drawn with straight lines over the rationals, that forces some
  excess.
- **Growth exponent.** If a certificate charges ψ_j for j directions and ψ_m for the whole transform, its growth
  exponent at block m is the least of ln(ψ_m/ψ_j)/ln(m/j). A certificate with growth exponent θ rules out savings above
  1 − θ, and no more.
- **Spectral factor, angle factor.** Two subspaces meet at a list of angles. A certificate's ratio at depth j is a
  spectral factor (how much total turning there is between far frames' subspaces compared with neighbours') times an
  angle factor (how thinly that turning is spread over many small angles). Section 7 makes this exact.
- **Mixed unit.** A unit that finishes some arrays only partly. "1 + 4 in 10" at block 3 does one full transform and
  four two-row transforms with 10 moves. Its value assumes a recursion with mixed sizes that our engine does not have
  yet.
- **Skeleton, relaxation.** A skeleton is a program whose steps are fixed and whose mixing numbers are unknowns. A
  relaxation is a smaller test that every unit of a skeleton must pass.
- **Approximate (border) unit.** A family of programs that come as close to a unit as you like, usually as some
  numbers grow without bound.
- **Σ_1 statement.** "There is a finite object that passes a mechanical check." "Some unit saves more than 0.1" is one.
- **s\*.** The best saving that units can reach, over all blocks.

## What I tried

Three rounds, all on 14 and 15 October 2026.

1. **The theory-of-theories round** (14 October, morning to early afternoon). A first cycle of five tracks: the
   algebra of combining units (spectrum), approximate units (border), a cohomology question left open by an earlier
   round, representation theory of certificates, and transplants from the literature on matrix-multiplication
   barriers. Then a map of everything known, as a lattice of open case × game × method, with 30 earlier predictions
   checked against their outcomes. Then four second-cycle tracks chosen from the map: lifts, the chi game's geometry,
   the everything-free game, and every block-3 certificate as a computation in one finite group. Skeptics checked
   nine of the tracks.
2. **The second abstract wave and its follow-ups** (14 October, afternoon to night). Five tracks, each writing its
   predictions before computing: Sylvester–Gallai and gluing for lifts, representation stability, GCT, ZX rewriting,
   and numerical algebraic geometry. Then four follow-ups: building programs by choosing letters, a search for a
   sharper certificate at the block-4 near miss, the skeptic that the algebraic-geometry track never had, and a full
   re-run of that track's 652 programs with every coefficient free.
3. **The logic round** (14 and 15 October). Three tracks: what kind of logical statement the question is; how hard
   impossibility proofs are and when one of our kind can exist; and whether existence can be proved by showing that
   every refutation fails. Each had one skeptic.

The searches ran on an Apple M4 Mac mini and on Linux machines of Google's Colab service (8 processors, about 50 GB).

## What I found

### 1. The saving is one number, and every method meets it through one exponent

Units combine: running two side by side, or nesting a block-k unit inside the paid moves of a block-m unit to get a
block-mk unit. The spectrum track wrote these rules as an algebra of the kind Strassen and Fritz use for matrix
multiplication, and showed by a hand argument that the algebra has only one dimension (Theorem C): every quantity that respects the
rules behaves like n^θ for a single θ, or is degenerate. So the best saving s\* is one number, the same whether one
uses blocks 2^N, 3^N or anything else. For matrix multiplication there are many such quantities, and finding new ones
drove progress. Here there are none to find. One detail is worth knowing: our algebra is not "Strassen preordered"
(an oracle for a large block cannot act on a narrower array), so Strassen's and Fritz's theorems do not apply as
stated, and Theorem C is argued directly.
*Status: argued by hand, not checked. One AI skeptic re-derived it and found no error.*

Every lower-bound method then meets the saving through one number per block, its growth exponent (Lemma D), and
mixing several certificates never beats the best single one (Lemma E, by Farkas' lemma). This gives a screen:

- To rule out a saving of 0.1, a family needs a growth exponent of at least 0.9 at *arbitrarily large* blocks. The
  first version said "at every block"; the skeptic corrected it, because a unit at block n can also be caught at
  block n², n³ or nk.
- A family that never proves more than 4 moves per array cannot rule out 0.1 at any block from 5 on, since such a unit
  needs at most 5^0.9 = 4.26 moves per array.

*Status of Lemmas D and E and the screen: argued by hand, not checked; one AI skeptic re-derived them.*

Two small facts came with it. One direction costs at most two calls of the whole block transform, by an exact
identity in Q(i) (*checked by computer* at blocks 1 to 3 and argued for every block). And units whose mixing never
divides by 2 never save, at any block, in any composition, because the mod-2 frame distance is then a perfect
certificate (Proposition G; *argued by hand, not checked*; it assembles two results of earlier rounds).

The growth exponents of the lifts I had ("start lifts") fall with the block: 0.85, 0.71, 0.57 and 0.44 at blocks 3
to 6. At the engine's own prices for grouped blocks, lifts cap the saving near 0.43 at block 5 and near 0.56 at block
6, where hooks allow about 0.83 and 0.91.
*Status: checked by computer (exact at blocks 3 and 4, with exact edges and sampled pairs at block 5, sampled at
block 6); one AI skeptic recomputed the block-4 pair sums exactly over 2,632,365 pairs and enlarged the block-6
sample, which moved its growth exponent to about 0.437.*

### 2. The ladder of games: what each game lets a certificate be

```
                 everything-free
                /               \
          widened             chi game
             |                    |
             |          unimodular chain maps
              \                  /
                    engine
```

Going up the ladder, more steps are free. Lower bounds flow down (a bound proved in a game holds in every game
below it), and units flow up. The round added a third rule: **ceilings flow up**. A certificate for a game with more
free steps is also a certificate, with the same numbers, for every game below it, so a ceiling computed in the engine
game holds in every game above it. *Status: argued by hand, not checked; two AI skeptics found no error.*

What each game allows:

- **Engine game.** A certificate loses nothing by ignoring Pauli flips (Theorem 1 of the representation track). So
  every engine certificate is a pattern of subspaces attached to the frames, and its ratio is read from averages over
  pairs of frames at each distance. The hook's frame pattern alone gives 14/5 at block 3: on all 18,090 ordered
  pairs of block-3 frames it takes one value per distance, 80, 144 and 224 at distances 1, 2 and 3, and 224/80 = 14/5.
  *Status: argued by hand, with the frame-pair values checked by computer; one AI skeptic.*
- **Chi game.** Every certificate must ignore Pauli flips (an earlier result), and the frames become 891 subspaces over
  the four-element field F_4.
- **Everything-free game.** Every certificate of our usual product form must ignore everything, so it proves nothing
  (Theorem 2 of the representation track). It rests on two identities, Tof·X0·Tof⁻¹ = X0·CNOT(1→2) and
  CS·X0·CS⁻¹ = X0·S1·CZ. *Status: argued by hand, not checked; one AI skeptic.* So a lower bound there must be
  arithmetic, or built from all column counts f at once. That is what later happened: the counting bounds of 15
  October (Theorem 2W and the two-layer bound, in [../block-bounds-and-new-rules/](../block-bounds-and-new-rules/))
  work in every game whose free steps touch a bounded number of entries, the everything-free game included.

### 3. How far each kind of certificate can ever go

The most moves per array each kind of certificate can ever prove. The plain cost is m; the bar needs more than
2.69 at block 3, 3.48 at block 4 and m^0.9 in general.

| kind of certificate | block 3 | block 4 | large blocks | status |
|---|---|---|---|---|
| colourings and coordinate frame families, engine game | 1.75 | 1.875 | below 2 | argued by hand, not checked; one skeptic |
| the same, chi game | 21/16 = 1.31 | | below 4/3 | argued by hand, not checked; one skeptic |
| cuts (no repeated pieces), engine game | 28/9 = 3.11 | 60/17 = 3.53 | below 4 | argued by hand, not checked; two skeptics |
| cuts, chi game | 28/11 = 2.545 | | | argued by hand, not checked; no skeptic |
| hooks | 2.8 (reached) | 10/3 (41/13 reached) | below 4 | earlier rounds (in [../block-bounds-and-new-rules/](../block-bounds-and-new-rules/)) |
| lifts, engine game | 2835/1013 = 2.7986 | 3.72; 80/23 = 3.478 for decomposed lifts | depends on how the excess grows | block 3: two skeptics; block 4: see section 6 |
| lifts, chi game | 147/52 = 2.827 (2.774 if a solver bound is right) | | | checked by computer; two skeptics |
| every certificate, engine game | 169/57 = 2.965 | 3.9426 | at most m − 0.057 | two programs of an earlier round; no skeptic; at block 4 the engine-game figure also needs one step that is only argued |
| every certificate, chi game | 169/57 = 2.965 | | | checked by computer; two skeptics, given the reduction of the earlier round |
| product-form certificates, everything-free game | nothing | nothing | nothing | argued by hand, not checked; one skeptic |
| the 2-adic bound (units that never divide by 2, right at f = 2) | 1.875 | 2.1875 | √m | earlier rounds, upheld by their skeptics |

What follows from the table (each step argued by hand, not checked; one skeptic or two as marked above):

- **No certificate whose pieces do not repeat can settle the bar.** All of them stay below 4 moves per array, which
  from block 5 on is below m^0.9 (at blocks 3 and 4 the caps 3.11 and 3.53 still lie above 2.69 and 3.48, so cuts can
  matter there). Hooks and colourings also stay below 4 at every block, and the coordinate sets of the earlier
  revisions stay below 4 at every block where they were computed (3 to 10). Only labels, lifts among them, can settle
  the bar.
- **Mixing certificates never helps, and there are no hidden spectral quantities to look for** (section 1).
- **Counting cannot prove a unit exists**, and a search residual says nothing about a game unless a planted control on
  the same menu shows that the search can succeed. A game in which a hand proof rules out every unit reached lower residuals
  than an open one on its larger menu.

The cut cap 28/11 for the chi game was derived in the round's summary from the frame graph's eigenvalues (degree 42,
second eigenvalue 9, so 2/(1 − 9/42)). The eigenvalue arithmetic is *checked by computer*
(`programs/frame-graph-eigenvalues/`), and it reproduces 5/14 at block 3 and 13/30 at block 4; whether the cut theorem's
proof carries over to the chi game is *argued by hand, not checked*. If it holds, no cut can ever exclude 3 in 8 in
the chi game.

### 4. A certificate from the group Sp(6,2) rules out four block-3 mixed units

In the engine game every block-3 certificate is an object of one finite group, Sp(6,2) (order 1,451,520): a
representation and a subspace kept by the stabiliser of one frame. That turned an open-ended search into a finite
list, the 30 irreducible representations of the group. The prediction, written before the calculation, was that no
certificate would touch the block-3 mixed units, which sit exactly on the hook bound. It was wrong.

**Theorem M.** In the 105-dimensional irreducible representation of Sp(6,2) whose value on transvections is −35, the
sum N of the weight spaces of one 28-element orbit has dimension 56 and charges 20 per paid move, 38 per two-row target
and 49 per full transform. The hook charges, in the same units, 10, 18 and 28. Theorem R then gives
20·M ≥ 49·W1 + 38·W2 for a unit with W1 full arrays and W2 two-row arrays, which fails for

| mixed unit | 20·M | 49·W1 + 38·W2 | saving it would give |
|---|---|---|---|
| 1 + 4 in 10 | 200 | 201 | 0.119 |
| 2 + 8 in 20 | 400 | 402 | 0.119 |
| 4 + 21 in 49 | 980 | 994 | |
| 7 + 38 in 88 | 1760 | 1787 | 0.125 |

So these four units do not exist in the engine game, and the engine game's room for mixed units falls from 0.1255 to
0.1160. The ranks were computed exactly over the rationals by fraction-free elimination and again by plain Gaussian
elimination. *Status: checked by computer; two AI skeptics each rebuilt the certificate independently, in a
different model (through the Weyl group of E7), with code that shares nothing with the original.* I re-ran both checks for
this revision: `verify_exact.py` (5 seconds) and the skeptic's `skeptic_e7_check.py` (8 seconds) both print the
profile (20, 38, 49).

The certificate does not apply in the widened game, whose free steps do not act through Sp(6,2). The same day a
construction round showed, by a different hand argument, that no hook-tight unit can give any array a two-row target in the
engine or the widened game (Theorem T, upheld by its skeptic); it excludes these four units in both games, and more.
That work is in [../moonshot-small-blocks/](../moonshot-small-blocks/).

**One spectral cap for every cut (Theorem S′).** The same computation produced a pattern: the eigenvalue ratio 5/14
appeared in 29 of the 30 irreducibles. It became a theorem for every block: every cut at block m proves fewer than
4(2^m − 1)/(2^m + 1) moves per array, which is 12/5, 28/9 = 3.11, 60/17 = 3.53, 124/33 and 252/65 at blocks 2 to 6,
and below 4 always. The proof compares the certificate with the frame graph's second eigenvalue, 2^m − 3 out of a
degree of 2(2^m − 1). A cut can beat a target only if it is balanced to within a few percent of half its
representation: at block 4, beating 10/3 needs balance within 2.9 points and beating 3.5 within 0.42 points.
*Status: argued by hand, not checked; two AI skeptics.* All 24,609 cuts of the 30 irreducibles were computed, and none
beats the hook's 14/5, which is itself a cut in the 56-dimensional irreducible (*checked by computer; one program, not
repeated by anyone*).

### 5. A Venn certificate with a spare row

The literature track transplanted Schönhage's τ-theorem: a mixed unit is worth the root τ of Σ W_k·k^τ = M, and mixed
units have more room than pure ones. A sweep of 180 certificates found room 0.1255 for mixed units at block 3, against
0.0628 for a pure kernel and 0.079 for pure rank-2 targets (*checked by computer; one skeptic*). As a side result it
found a new certificate for rank-2 targets with one spare row, in the space Λ³V ⊗ Λ⁴V with intersection pattern
{0, 2}: cap 752 (exact over the rationals), rank-2 demand 1,424 and full demand 1,712 (modulo two primes, which can
only understate them). So a rank-2 target with one spare row needs at least 1424/752 = 89/47 = 1.8936 moves. That
excludes 7 arrays in 13 moves (which would need 13.255), 8 in 15 and 9 in 17 (17.043), in the widened game and so in
the engine game. The smallest open case is 10 in 19 (saving 0.074). *Status: checked by computer; one AI skeptic, who
recomputed the cap exactly with its own code.* Re-run for this revision: the skeptic's exact cap (half a second) and
the modular check (37 seconds, about 1.5 GB of memory) both reproduce 752, 1,424 and 1,712.

### 6. Lifts: dead at block 3, undecided at large blocks

Lifts are the one family the screen of section 1 leaves alive: a lift family settles the bar only if its excess per
edge stays at or below m^0.1 − 1 at arbitrarily large blocks. That threshold is 0.116, 0.149, 0.175 and 0.196 at blocks
3 to 6, and it keeps growing slowly.

What the second wave showed (*argued by hand, with exact linear-programming certificates; two AI skeptics found no
error in the parts listed here*):

- **Block 3.** An interval lift with a bad edge at either end costs at least 3, not 2, and an explicit exact lift
  attains 3. Every edge is an end edge in 64 of the 160 intervals through it. Counting with these facts, every lift of
  the whole block-3 frame graph has excess at least 68 on its 945 edges (it was 54). So **every block-3 lift proves at
  most 2835/1013 = 2.7986 moves per array, below the hook's 2.8**: lifts can never beat the hook at block 3, and no lift
  can rule out 5 in 14 or 6 in 17 (they would need excess at most 67 and 55). The true whole-graph optimum lies
  between 68 and 162, the best lift known.
- **Block 4 and 5.** Every block-4 interval lift has excess at least 18 (an exact dual certificate of value 120/7,
  *checked by computer*), so 18 ≤ b_4 ≤ 36; block 5 has at least 141 (dual value 5611/40). The excess per edge is at
  least 0.072 from block 3 and 0.075 from block 4. Decomposed lifts, the kind every good lift known so far is, need
  exactly 36 on a block-4 interval, so they prove at most 80/23 = 3.478 at block 4, a hair under 4^0.9 = 3.482. (The
  36 comes from an exhaustive search of the first-cycle lifts track over 32 million sets of lines; the second wave's
  skeptics used it and did not refute it, but none re-derived it.)
- **Sylvester–Gallai transfer (Theorem SG).** Through the rank bound of Dvir, Saraf and Wigderson for design
  matrices, decomposed lifts pay at least about 1/3 per edge at large even blocks and 0.314 at large odd ones. It
  takes effect only from about block 110.

What failed. The track claimed that every counting argument of this Sylvester–Gallai kind is capped at exactly 1/3 per
edge, so that the bar's threshold overtakes it from block 18 on. **Refuted** by both skeptics. What survives: such
counting is capped by some constant, so it can never prove that lifts fail at all large blocks. The linear version is
capped at (2/3)κ_m = 0.400, 0.583, 0.299 and 0.387 at blocks 3 to 6; a sharper version proves at least 0.38 per edge;
the only cap shown for the whole class is about 1 per edge, which the threshold passes only at block 1024. So where the
threshold overtakes these methods is unknown.

The constructions point the other way. The best decomposed lifts found have 0.057, 0.150, 0.279 and 0.457 excess per
edge at blocks 3 to 6, growing faster than the threshold, and start lifts get worse with the block (ratios 2.26, 2.18
and 2.02 at blocks 8 to 10, against a predicted 2.35 to 2.65). *Status: checked by computer; this is evidence about
the lifts I can build, not a proof about all lifts.* So whether lifts can settle the bar is **open**, and the evidence
leans towards no.

For readers of the fourth revision's note for network-coding readers: its question 3 asked whether the
seven-subspace inequalities give a positive lower bound on how far a labelling of the block-3 frame graph must be from
exact. For lifts, the labellings that reduce to the frames modulo 2, there is now such a bound: total excess at least
68 on the 945 edges. It comes from counting the bad edges that Fano planes force, not from the two inequalities that
question names, so the question as asked stays open.

### 7. The angle law, and the block-4 near miss as one number

**Theorem F and the angle law** (stability track). For any certificate made of equal-dimension subspaces on the frames,
at every block up to 40, the spectral factor from edges to depth j is at most 2 − 2^(1−j): 1.5, 1.75 and 1.875 at
depths 2, 3 and 4, and this is the exact maximum there. At larger blocks only a cap of 8/3 is shown. A certificate's ratio at depth j is its
spectral factor times its angle factor, where the angle factor is s_1/s_j and s_d is the mean sin² per nonzero angle
between subspaces of frames at distance d. Each nonzero angle counts as a whole unit, however small, so spreading the
same turning over more, smaller angles raises the count.

**Theorem H.** Hooks use the full spectral factor at every depth, block and k, so a hook's whole value is its angle
factor, which stays below 2. The skeptics added that every hook has right-angle edges (s_1 = 1).

*Status of Theorem F, the angle law and Theorem H: argued by hand, not checked; two AI skeptics, who strengthened
them.* Consequence: settling the bar at block m (up to 40) needs an angle factor of at least m^0.9/2, which is 2.88 at
block 7 and 3.97 at block 10. Every family computed stays below 2: block-3 cuts at most 1.85 (11,048 cuts), hooks below
2, start lifts falling from 1.66 to 1.45 by block 10. Representation stability cannot decide blocks 5 and up.

**The block-4 near miss.** Doing 2 arrays in 7 moves at block 4 would save 0.096, just under the bar. In the branch
where each array ends on a private move, the first five moves must build a depth-3 piece for each array. Under the best
hook (k = 8) one move is worth at most 27,456, one depth-3 piece needs 68,320, and five moves supply 137,280 against the
136,640 needed: 640 of room, under half a percent. The angle law shows where that room sits. The hook's depth-3 ratio
68,320/27,456 = 2.488 is the full spectral factor 7/4 times an angle factor of 1.4219; its edges are already at right
angles; its average sin² per angle at depth 3 is 48,048/68,320 = 0.7033, and closing the branch needs an average below
0.7. A skeptic recomputed 68,320 exactly over the rationals, so the shortfall is real. *Status: checked by computer;
two AI skeptics.*

**The search for a sharper certificate.** A follow-up track searched about 4,000 certificates: relatives of the hook, cuts
in the irreducibles of Sp(8,2) of dimensions 85, 238, 510 and 595, all 2,047 orbit cuts, and overlap cuts. No hit. The
best is still the hook at 2135/858 = 2.488; the next are 2.4868, 2.4858 and 2.4855; the four named irreducibles reach
only 13/7, 9/4, 11/6 and 7/3. The block-3 precedent did not carry over: at block 4 the hook ties several irreducible
pieces together, and each alone reaches at most 1.79. Building block-4 certificates from block-3 ones would need a
block-3 certificate above 2.849, and none is known. *Status: checked by computer; one AI skeptic confirmed the main
result with its own code (the hook's numbers modulo two primes, the orbit cuts, the overlap cuts by a different
reduction) and computed the 22 orbit cuts the search had missed (best 2.304); the hook-relative runners-up were not
rechecked.*

The skeptic also **refuted** one sentence, and with it the closing test as first stated. Orbit cut 1434 has the full
spectral factor 7/4 and a depth-3 average of 0.699, below 0.7, yet it misses (2.483) because its edges are not all
right angles (s_1 = 858/865 = 0.992). The correct test: with the full spectral factor, a certificate closes the branch
exactly when its depth-3 average is below 0.7 times its edge average. Re-run for this revision: the skeptic's
`vc_1434.py` reproduces cut 1434 and its mirror 854 (14 seconds, about 3 GB of memory).

The near miss stands. Later rounds narrowed 2 in 7 from the construction side; that is in
[../moonshot-block4-block5/](../moonshot-block4-block5/).

### 8. The chi game: every certificate stops at 169/57

With the running XOR free, three block-3 targets were open: 3 in 8 (saving 0.107), 2 in 5 (0.166) and 3 in 7 (0.229).
They need more than 8/3 = 2.667, 5/2 and 7/3 moves per array. The best certificate known in that game proves
231/131 = 1.763. The chi game's frame geometry over F_4 contains 360 copies of the Fano plane where the engine game has
one, and the prediction was that they would push the ceiling below 8/3. It did not.

- **The ceiling is 169/57 = 2.965 for every certificate of the chi game**, the same as in the engine game. One proof is
  the "ceilings flow up" rule of section 2. The Hermitian track also showed it natively, with an exact 13-row
  certificate of rank inequalities on the 44-frame interval over F_4 (submodularity, monotonicity, Ingleton's
  inequality and one Dougherty–Freiling–Zeger inequality, the one that fails in characteristic 2). A single Fano
  gadget, with weight 2/57, does all the work; the other 359 buy nothing. *Status: checked by computer; two AI
  skeptics rebuilt the geometry and the rows with their own code. Both routes rest on a reduction from an earlier
  round (from every certificate to rank inequalities on one interval) that no skeptic has checked.* Re-run for this
  revision: `verify_q.py` (16 seconds) and the skeptic's `skeptic_mycert.py` (13 seconds) both verify 169/57.
- **So rank inequalities on one interval cannot decide any of the three cases.** All three thresholds lie in the window
  between 1.763 and 2.965.
- **Lifts survive locally, by a hair.** The least excess b of a lift of one 44-frame interval lies between 9 (an exact
  integer-programming optimum, which a skeptic confirmed by showing that 8 is impossible) and 18 (an exact lift found
  by annealing, rechecked by both skeptics). Its 9 bad edges are exactly the 9 tangent flags of the Hermitian curve.
  With extra rows from Hirzebruch's inequality a solver bound gives 12. With b = 18 the local ceiling is
  147/55 = 2.6727, above 8/3 by 0.006; had b been 19, no lift could ever have excluded 3 in 8. A lift of the whole
  891-frame geometry excludes 3 in 8 only with total excess at most 2,338 on its 18,711 edges (3,742 for 2 in 5, 5,345
  for 3 in 7); the best known has 13,122. *Status: checked by computer; two AI skeptics.*

The chi-game cluster is **open**. Only labels, lifts among them, a bound from the representation theory of the
unitary group acting on these frames, or a construction in which the running XOR does essential work, could settle
it.

### 9. Cohomology: the hook bounds carry over to Clifford chain maps

An earlier round had left one question: when the free steps include chain maps (XOR of one digit column into the
next), does a certain sign defect obstruct carrying the hook certificates over? The cohomology track showed that the
defect vanishes identically, that the naive way of collapsing a Clifford cellular automaton to one column is
multiplicative with signs included, and so that the difference of two cohomology classes is zero. **Theorem E+:** when
every letter is a locality-preserving Clifford map with a translation-invariant bulk (one-column Cliffords, the
directions, the block transform, every unimodular chain map such as row0 ^= row1 << 1 and its mirror, CZ and i-phase
couplings, cyclic column rotations, and every Clifford acting near the boundary), and every free letter is monomial or
acts near the boundary, then every certificate kept by monomial Clifford tensor powers gives the same bound as without
chain maps: at least 2W moves at block 2, 14W/5 at block 3 and 41W/13 at block 4.
*Status: argued by hand, with randomized exact checks; one AI skeptic re-checked it with its own code on 1,440 random
Clifford automata and 720 pairs (all multiplicative, defect zero on 8,640 pairs of vectors), with a control showing
that the test has power.* Re-run for this revision: the skeptic's three checks reproduce (about 2 seconds in all).
The skeptic noted that 41W/13 excludes neither 2 in 7 nor 3 in 10 at block 4. The same mechanism cannot carry
Toffoli: a necessary trace condition fails for Clifford + Toffoli at every odd ring length from 3 to 41 (*checked by
computer in floating point; one program, not repeated by anyone*). So the mixed case of chain maps with non-Clifford
letters is open.

### 10. Approximate units

For matrix multiplication, approximate (border) schemes cost nothing in the exponent. For us they need care, because
our recursion substitutes approximate sub-calls into each other instead of taking tensor powers.

- **Theorem X.** Call an approximate unit *formally polynomial* when its cancellations of negative powers of ε hold
  with the paid moves treated as unknown symbols. Such a unit runs exactly at a cost factor of order log log N, so its
  saving is unchanged. *Status: argued by hand, not checked; one AI skeptic found no error.*
- **Refuted: "making any other approximate unit exact costs at least a factor 2 per level."** The skeptic built a
  counterexample from the track's own toy unit: copy-by-copy exactification needs 3 coefficient-level paid calls for
  M = 3, and 1 after pruning, so the factor can be 1 (`programs/exactification/skeptic_exactify_count.py`, re-run for
  this revision). The correct statement: the cost lies between the exact cost and (a + 1)·M, and only moves whose
  inputs carry several ε-coefficients pay extra. So approximate units are worth whatever their exact version costs,
  and they are back on the table at blocks 3 and 4.
- **Refuted: "a family of approximate units gives a real unit only if its error is exactly zero in the limit"** (from
  the logic round). The skeptic's counterexample is to scale by 1 + 1/k. Exact units form a closed set, so every
  family with bounded numbers whose error tends to 0 has an exact unit as its limit; only families whose numbers grow
  without bound can fail.
- If a target has an approximate scheme, no certificate of our kind can exclude it (an earlier round's result,
  upheld by its skeptics). So for a target whose only schemes are approximate, rank methods are blind.

Making formally polynomial approximate units exact with the same number of moves is in
[../block-bounds-and-new-rules/](../block-bounds-and-new-rules/) and [../existence-routes/](../existence-routes/).

### 11. One array needs m moves in every game; ZX rewriting

**One array needs m paid moves, in every game** whose free steps do a bounded amount of work per output entry,
including the chi game and the everything-free game. The proof counts nonzeros along paths. This settles W = 1 in every
game. *Status: argued by hand, not checked; upheld by AI skeptics in two rounds.*

**ZX rewriting.** The ZX calculus is a diagram language for quantum circuits, with a standard optimiser, PyZX. The
track asked whether its rewriting can reduce the number of paid moves.

- Without sums of arrays, no diagram goes below m paid moves per array. PyZX agreed on 162 circuits (378 optimiser
  runs): none went below. The proof first assumed that every free step only permutes and rephases entries, which is
  false in games with Hadamards on the end columns. Both skeptics supplied the same repair: p·f + 2mWL ≥ mWf at every
  f forces p ≥ mW once f > 2mWL. *Status: argued by hand, not checked; one skeptic refuted it as first written, both
  skeptics agree that the repaired core holds.*
- PyZX's costs do not match ours: free steps show up as Hadamard edges, and PyZX counts a move applied to two arrays
  as one. Fusing moves never saves. With sums, ZX's strongest fragment is a slice I had already searched. Three
  apparent hits each broke a rule of the game: the same program at every column count, cost per array, and no
  post-selection.
- **A bug in PyZX 0.10.7.** Its `basic_optimization` changes the linear map of circuits that contain XPhase gates. The
  skeptic shrank it to three gates on three qubits: CNOT(2,1), CNOT(1,2), XPhase(1, phase 1/2). I re-ran the
  reproduction for this revision against PyZX 0.10.7 from PyPI, and it still fails the same way. I have not reported
  it to the PyZX maintainers yet.

ZX gives no new search family and no new proof. Its lasting use is a list of rules that any diagram search must obey.

### 12. Geometric complexity theory

GCT proves lower bounds with polynomial equations that every M-move program satisfies but the target does not,
organised by symmetry.

- The best GCT-type obstruction proves exactly the border complexity (what approximate programs cannot reach), and the
  barrier that caps rank methods caps it too. An obstruction that sees only some irreducible pieces is capped by a
  smaller game played inside those pieces, and that game can be computed. Pieces interact only through shared
  sub-pieces. *Status: argued by hand, not checked; no skeptic ran.*
- At block 3, in every single irreducible the track could decide (dimensions 15, 21, 27, 35 and 56), the smaller game
  costs exactly what the best cut proves. Inside the hook's own 56-dimensional piece, an exact 5 in 14 and an exact
  mixed 1 + 4 in 10 both exist (least-squares residuals 5.6e-13 and 3.4e-29). So **no argument confined to the hook's
  piece, GCT or any other, can rule out 5 in 14**, and the exclusions of 1 + 4 in 10 (Theorem M in section 4, and
  Theorem T) must use other pieces. *Status: open as a theorem; this is numerical evidence (a residual near zero is
  not an exact solution), no skeptic.*
- The first coupled pair of pieces (56 + 105) descended steadily towards an exact 5 in 14: residual 5.4e-9 after 1,500
  iterations, still falling, with bounded numbers. *Provisional; open.* Two other combinations stayed undecided because
  the solver failed even on a program known to exist, which sits on a lower-rank stratum that ordinary descent misses.
- At block 2 in the everything-free game, the free steps generate the whole matrix algebra at 1, 2 and 3 digit
  columns, so finite-dimensional GCT proves nothing there.

GCT is not a separate route to the bar: whatever it could add lies in couplings between pieces, where labels and lifts
already live. The GCT programs are not in this revision (see [Programs, data and notes](#programs-data-and-notes)).

### 13. Exact algebra and homotopy continuation on symmetric block-4 programs

The algebraic-geometry track asked whether units hide in "special position", coefficient choices that random sampling
misses but exact algebra finds. It took 652 symmetric block-4 skeletons from a construction round's builder (562 aimed
at 2 in 7, 81 at 3 in 10, 4 at 4 in 13 and 5 at 7 in 23), decided each exactly over the rationals with Gröbner bases,
and reproduced the answers by homotopy continuation.

- **Refuted: "all 652 are impossible for every complex choice of coefficients."** The track had fixed the first
  coefficient of every move to 1, which leaves out "first coefficient = 0". Its skeptic found that with every
  coefficient free, 4 of 93 sampled skeletons pass the relaxation exactly, the track's own "border" skeleton among
  them; the "border" was an exact solution at coefficient 0, not numbers running off to infinity. 157 of the 562
  "2 in 7" skeletons really cost 5 or 6 moves, and the 652 are a random sample from one builder, not every symmetric
  program.
- **What survives, with every coefficient free.** A full re-run found that 17 of the 652 pass the relaxation, all of
  them 2 in 7 at a true cost of 7 (17 of 405, about 4 %); none of the 3-in-10, 4-in-13 or 7-in-23 skeletons pass.
  **16 of the 17 have no unit for any complex coefficients**: 15 excluded by exact computation over the rationals, and the 16th with
  its first steps (which coefficients must vanish) done modulo two primes. **One, g27_Fany4_s7, is open**: 49 of its 64
  charts are excluded exactly over the rationals, 15 are open, and in each open chart the full system fails at random
  points. *Status: checked by computer; one AI skeptic rebuilt every system with its own code, found zero
  disagreements on all 677 relaxations it decided (651 skeletons and 26 calibrations), and independently re-decided 9 of
  the 15 exact exclusions over the rationals; five of the other six were only screened at one exact point each, and one
  (g27_mid3_m3_s13) was not rechecked.* (The round's notes give 45 of 64, a count made
  during the run; the final count, 49, is the later one and is used here.)
- **The pattern behind each failure** (*argued by hand, not checked*): the row for the target equals a multiple of a
  few "junk" rows that only the target's own moves touch. When the multiple is a constant, solvers stall at 1/2, 1/3
  or 1/4, the floors earlier searches kept hitting. When it is a polynomial, the residual shrinks as numbers grow but
  never reaches zero. It held in 81 of 82 checked certificates, and only where the first coefficient is nonzero.
- **Choosing letters instead of numbers** (*checked by computer, one program, not repeated by anyone*). A follow-up
  changed letters so that some other move lands on the junk. On controls it found real special position: 4 of 1,116
  single letter replacements pass, and all 4 are exact programs. It produced the first 71 two-in-seven programs that
  pass the relaxation without dropping a summand (40 symmetric, 31 asymmetric). All 71 fail the full test for one
  reason: the only coefficients that pass turn a two-letter step into 1 + u with u·u = 1, and (1 + u)(1 − u) = 0, so
  the step destroys part of the target it carries. About 9,400 single and 136 paired letter changes all kept this
  zero divisor. The one aligned 3-in-10 skeleton fails its relaxation exactly. Caveat: this follow-up inherited the
  normalisation hole, so the cases with a first coefficient of 0 are unchecked.
- A further check on these 17 hits is in [../existence-routes/](../existence-routes/): every one of them is feasible
  over the complex numbers at the relaxation level, so that relaxation can exclude none of them.

**A warning for anyone who uses msolve.** msolve 0.6.5 reads a monomial written `x*x` as `x`; exponents are not
accumulated. The system {x·x·y·z − x·x, x·x·y·z} written with `*` gives the basis [x], and written with `^` gives
[x²]. Always write powers as `x^k`, and cross-check every msolve verdict with Singular. It invalidated several
full-system verdicts of these rounds, which were all redone; the relaxation systems never repeat a variable and were
not affected. The one msolve output from before the bug that is published is the relaxation table
(`data/symmetric-block4-skeletons/relaxation_table.txt`, two primes): the bug cannot touch it, because relaxation
systems never repeat a variable, and a skeptic re-decided 677 of its 678 systems with Singular over the rationals,
with no disagreement (the 678th it decided modulo two primes, with the same answer).

### 14. Logic, part 1: what kind of statement the question is

All three logic tracks had one skeptic each, and all three verdicts were "refuted in part": in each case the
mathematics held and some readings of it were too broad. Everything here uses the corrected versions.

- **Every f reduces to finitely many f (Theorem F of the logic round, not the one of section 7).** Some free letters act on the digit columns like a small machine
  reading them one at a time: column-uniform letters, chain maps, carries, neighbour phases. For those, a bound N(S)
  computable from the shape alone suffices: right at f = 1 to N(S) means right at every f. So each shape is one finite
  polynomial system over the complex numbers, and Gröbner bases decide it. *Argued by hand, not checked; one skeptic
  found no error.*
- **In the engine game the bound depends only on the block.** The engine's letters and paid moves generate a finite
  group: 96 elements at block 1 and 46,080 at block 2. Being right at every f, the formal identity in the group algebra,
  and being right at f = 1 to r′ are equivalent, where r′ is the number of distinct nonzero traces: 12 at block 1, 20 at
  block 2, at least 28 at block 3 (sampled). At block 1 the shortest range that always suffices is f = 1 to 8: some
  combination of group elements vanishes at f = 1 to 7 and not at 8. A chain shift by e looks right up to f = e and
  first fails at f = e + 1, so chain-map games have no block-wide bound, and with controlled-S the block-2 group is
  infinite. *Checked by computer; the skeptic reproduced every number in exact arithmetic.* Re-run for this revision:
  `fbound.py` (3 seconds) and the skeptic's `engine_group_exact.py` (90 seconds) reproduce them. This does not replace
  the missing Lean lemma of the earlier rounds, as the track first said: its route to "right at every f" passes
  through that lemma (**refuted** reading).
- **The logical form (Theorems A and I).** "Some unit saves more than q" is Σ_1; "no unit saves more than q" is the
  opposite form, Π_1. s\* is computable from below: list the shapes, decide each, keep the running best. A true
  Σ_1 statement is provable in a very weak theory (Robinson arithmetic) by exhibiting the unit, so **the existence side
  can never be independent**. No bounded question ("no unit with at most B steps saves more than q") can be either,
  since it is decidable. **Independence could live only in the unbounded ceilings, exactly as far as s\* fails to be
  computable**: if some sound theory proved ceilings arbitrarily close to s\*, listing its proofs would compute s\*;
  conversely, if s\* is computable, weak arithmetic plus one true axiom proves every true ceiling. Whether s\* is
  computable is open, as it is for the Shannon capacity of graphs and for asymptotic tensor rank. *Argued by hand, not
  checked; one skeptic found no error. One further item, that the full circuit model drops to Σ_1 exactly if every
  circuit saving comes from a finite unit, is argued in one direction only.*
- **Self-reference gives nothing for free (Theorem P).** For a Σ_1 statement, being true, being provable in weak
  arithmetic, being provable in ZFC and being provably provable all coincide. Löb's theorem, reflection principles,
  consistency axioms and large cardinals add no new statements of this kind, only shorter proofs; reasoning about
  provability alone, without the unit equations, would prove 0 = 1. A unit in a nonstandard model shows only that
  impossibility cannot be proved. The Lefschetz principle: a fixed shape has a complex unit exactly when it has
  solutions over the algebraic closure of the field with p elements for infinitely many primes p. *Argued by hand, not
  checked; one skeptic.*
- **The one loophole is proof length.** A short abstract proof of a unit too large to write down is logically
  possible, and it would count for OpenAI's challenge. Logic gives no way to find one. A search that finds no unit up
  to size B rules out explicit constructions of about that size, and nothing abstract.
- **The record of this track is not genuine.** The skeptic found that every file of the track was written within about
  ten minutes, before most of the times its files give for their predictions, runs and reading. So "predictions written
  before any computation" cannot be verified for this track. The mathematics does not depend on the record, but the
  discipline was not kept, and I say so here.

**What it means.** There is no logical shortcut to existence. Any existence proof must work on the unit equations:
algebra, structure, or a contradiction argument like OpenAI's prices.

### 15. Logic, part 2: when can an impossibility proof of our kind exist?

- **Our certificates are all closure-type**: each finds a polynomial that vanishes on everything a skeleton can
  output, even approximately, but not on the target. Theorem R, hooks and lifts are of this kind.
- **Proposition D.** For each skeleton exactly one of two things exists: such a closure polynomial, or a *border
  curve* (an approximate unit), which proves that no closure-type proof can work. Searching for both by size, one
  search must stop. *Argued by hand, not checked; one skeptic found the proof sound.* The skeptic's correction:
  Theorem R, hooks and lifts form narrower classes than closure-type, and for them only the "no" side is decided. A
  border curve proves that no Theorem-R certificate works; no curve does not produce one. The fact that a closure-type
  proof exists exactly when no approximate unit sits on the skeleton is a standard fact of GCT (Grochow).
- **Three 2-in-7 skeletons at block 4 have exact border curves on their relaxations**: g27_Fany4_s5, g27_Fany4_s6 and
  g27_Fp2b_m3_s5. Their relaxations have no exact solution, yet no closure-type proof can refute those relaxations. The
  skeptic confirmed the curves at two new primes, Singular agreed on two, and one curve was refined to 80 digits:
  along it the relaxation's error is 0.447·ε while its numbers grow like 0.57/ε. *Checked by computer; one skeptic.*
  A fourth, f27_FAmid3_m3_s16, has a curve found numerically only and not reproduced. Re-run for this revision: the
  skeptic's 80-digit refinement `refine.py` (12 seconds) reproduces the curve. **Refuted:** the step from these curves
  to "no Theorem-R certificate can ever exclude 2 in 7 at block 4"; it needs the border to be right at every f with
  legitimate letters, and the skeleton to lie in the branch the certificate search targets. The full systems of these
  skeletons are not decided, so whether approximate 2-in-7 units exist there is **open**.
- **Refutations of our impossible systems are cheap.** 1,304 infeasible relaxations at blocks 3 and 4 and 58 skeletons
  decided by Singular refute at degree 2 to 6, mostly at block size + 1. The degree tracks the block size, not the
  number of moves; the track had predicted the opposite. At block 4, 2 in 7 and 7 in 23 both refute at degree 5. The
  skeptic redid the 58 exactly over the rationals: 15 at degree 2 and 43 at degree 3. *Checked by computer for the 58
  (one skeptic); the 1,304 were computed modulo a prime only, not repeated.*
- **Halving shows up at block 4.** 15 of the 17 two-in-seven relaxations that have complex solutions have no solution
  over the two-element field, so every solution of those 15 must divide by 2. *Checked by computer; one program, not
  repeated by anyone.*
- Approximate units force long Nullstellensatz identities but not long Gröbner proofs: on a test chain the identity
  degrees are 4, 7 and 12 while Polynomial Calculus needs only 3 to 5. The famous barriers mostly do not reach our
  games: the rank barrier of Efremenko, Garg, Oliveira and Wigderson applies to Theorem R directly (our ceilings 169/57
  and 3.9426 are concrete cases of it), while the algebraic natural-proofs and bounded-arithmetic results concern
  growing families in the full circuit model and say nothing about one fixed shape.

### 16. Logic, part 3: proving existence by showing that every refutation fails

- **The principle is exact and the degrees are hopeless.** A solution exists exactly when no refutation of degree up
  to the effective-Nullstellensatz bound exists (Kollár). On our systems that bound means linear algebra with about
  10^610 unknowns for the smallest real 37-unknown skeleton, and 10^(8.7·10^23) for the exact equations of 2 in 7 at
  block 4, against machines that handle about 10^7. Structure brings 10^610 down to 10^357 at best. *Checked by
  computer; one skeptic.*
- **Two certificates, exactly one exists.** At a fixed shape, either a refutation identity 1 = Σ g_i·f_i exists, or a
  set of commuting rational matrices at which every equation vanishes: a matrix model of the unit, whose common
  eigenvector gives a solution. The matrix size is at most the product of the equation degrees (Heintz's Bézout
  inequality). So existence at a fixed shape cannot be independent of ordinary arithmetic. *Argued by hand, not
  checked; one skeptic confirmed the theorem and its source.*
- **A working pipeline and six certificates.** The pipeline searches modulo several primes, rebuilds rational
  matrices and keeps only what passes an exact check. Six certificates are on file with a standalone checker
  (Karatsuba multiplication, the W-state tensor, a special-position control, and relaxation hits whose solutions need
  i); a tampered certificate is rejected. *Checked by computer; the skeptic rechecked all six with independent code,
  including that each belongs to the right original system.* Re-run for this revision: `check_cert.py` prints VALID for
  all six (4 seconds). None of the six is a unit of our problem; the pipeline was never run on one, and it stalls at a
  few dozen unknowns (it did not reach Toom-3 with 40 or Strassen with 59).
- **Plain Bézout counting is blocked.** For a square Karatsuba slice the count is 1,728, but only 4 solutions are
  finite; the rest lie at infinity. *Checked by computer; one skeptic.*
- **Refuted readings.** "A proof that an existence proof exists is never cheaper than the existence proof", "counting
  is the only non-constructive route", "counting can prove only border units", and "modulo p is never a proof over the
  complex numbers" are all false. The skeptic computed a counterexample to the first and last: for a square slice with
  a nonsingular root, one prime plus Hensel's lemma proves existence; for the 12 × 12 Karatsuba slice p = 10007 works,
  and the certificate is 12 residues (re-run for this revision). The real obstacle on our systems is a different one:
  they have more equations than unknowns.
- **A sweep extension, modulo one prime.** Of 133 skeletons an earlier search had left undecided, 76 are refuted at
  degree 2 or 3 modulo one prime. That is not a proof over the complex numbers (p·x = 1 is refuted modulo p and has the
  solution x = 1/p); the 76 are very likely right but stand below Singular's verdicts. *Open until rechecked.*

### 17. The general pattern: how to show that no proof of a given kind exists

All three logic tracks met the same pattern: to show that no proof of some kind exists, exhibit the object that
fools that kind. Only a real unit fools every impossibility proof.

| kind of proof | what shows that no proof of this kind exists | decidable at one fixed skeleton? |
|---|---|---|
| a proof that a unit exists | a refutation (existence is provable exactly when it is true) | yes: a matrix model or a refutation |
| Nullstellensatz identity of degree d | a pseudo-solution of degree d | yes, by linear algebra |
| Polynomial Calculus of degree d | a pseudo-solution closed under its derivation steps | yes, by linear algebra |
| closure-type (polynomials vanishing on all outputs, approximate ones included) | a border curve | yes (Proposition D), given finitely many letters and a finite every-f system |
| Theorem R, hooks, lifts | a border curve shows none can work; no curve does not give one | only the "no" side |
| mod-2 arguments (the halving results) | a solution over the algebraic closure of the two-element field | yes, by Gröbner bases over that field |
| all Theorem-R certificates at one block | the ceilings on Theorem R (169/57 at block 3; 3.9426 at block 4, no skeptic) | ceilings computed by an earlier round |
| every sound impossibility proof | a real unit | at a fixed shape, yes; for the unbounded ceiling, possibly independent, exactly when s\* is not computable |

### 18. Claims refuted along the way

| claim | refuted by |
|---|---|
| A family needs growth exponent 0.9 "at every block" to settle the bar | the spectrum skeptic: "at arbitrarily large blocks" (a unit at block n is also caught at n², nk) |
| Making a relation-using approximate unit exact costs at least a factor 2 per level | a skeptic's counterexample with factor 1 (section 10) |
| No certificate touches the block-3 mixed units | Theorem M (section 4) |
| The old spectral bound caps most irreducibles below 2.8 | the same eigenvalue ratio 5/14 appeared in 29 of 30 (it became Theorem S′) |
| The chi game's 360 Fano planes push its ceiling below 169/57 | the native certificate: exactly 169/57 (section 8) |
| A chi-game interval lift needs excess about 30 | an exact lift with excess 18 |
| A block-4 interval lift needs excess between 35 and 38 | the upper end: 36 suffices; only 18 is shown for every lift |
| With controlled-S free, searches for 2 in 3 at block 2 do far better | no 2-in-3 found, and controlled-S never beat CZ in the same places (two skeptics); later the case was closed by Theorem 2W |
| Sylvester–Gallai counting is capped at 1/3 per edge, overtaken by the bar from block 18 | both lifts skeptics (section 6) |
| Every free step is a permutation with phases (ZX track) | false with Hadamards on end columns; proof repaired (section 11) |
| All 652 symmetric block-4 skeletons are impossible | the normalisation hole (section 13); 16 of 17 hits excluded, one open |
| "The construction side's near miss was not a near miss" because numbers ran off to infinity | the border was an exact relaxation point at coefficient 0 (section 13) |
| A closing certificate for the near miss needs only full spectral factor and depth-3 average below 0.7 | orbit cut 1434 (section 7); the test is depth-3 average below 0.7 times edge average |
| Approximate families give real units only if their error is exactly zero in the limit | scaling by 1 + 1/k (section 10) |
| A Lean-free finite check (Corollary F1) replaces the missing Lean lemma | its route passes through that lemma (section 14) |
| Four broad readings of the existence-by-refutation track | Hensel lifting and the other routes (section 16) |
| No Theorem-R certificate can ever exclude 2 in 7 at block 4 (from the border curves) | the step is incomplete (section 15) |
| Start lifts keep a ratio of 2.35 to 2.65 at blocks 8 to 10 | measured 2.26, 2.18, 2.02 |
| The refutation degree grows with the number of moves | it tracks the block size (section 15) |

### 19. Statements of these rounds that later work changed

Where the notes of 14 October say something that later results changed, the later verdict wins.

| statement on 14 October | what holds now | where |
|---|---|---|
| 2 in 3 and 3 in 5 at block 2 in the everything-free game: "nothing excludes them" | block 2 never saves in any game whose free steps touch a bounded number of entries (Theorem 2W, 15 October, one skeptic) | [../block-bounds-and-new-rules/](../block-bounds-and-new-rules/), [../moonshot-small-blocks/](../moonshot-small-blocks/) |
| the proof that W arrays need W + m − 1 moves "cannot carry over" to the everything-free game | the two-layer bound W + m holds there and in every bounded game (Theorem T′, upheld) | [../block-bounds-and-new-rules/](../block-bounds-and-new-rules/) |
| 2 in 6 at block 5 "misses by 0.008" | excluded by W + m | same |
| mixed 1 + 4 in 10 and 7 + 38 in 88 "open in the widened game" | excluded by Theorem T (upheld) | [../moonshot-small-blocks/](../moonshot-small-blocks/) |
| lifts: "b_3 = 2", "b_4 between 35 and 38" | every block-3 lift ≤ 2.7986; 18 ≤ b_4 ≤ 36 (section 6) | this section |
| "every symmetric program of this shape is exactly impossible" | 651 of the 652 excluded, one open; the sample is not every symmetric program | this section |
| a saving unit "relabels addresses or uses phases" | at block 3 it must use address maps and phase letters together; either kind alone gives exactly the plain cost (Theorem C3, 15 October, one skeptic) | [../reed-muller-and-one-move/](../reed-muller-and-one-move/) |

The working map of 14 October is published as written in
[notes/map-of-statements-14-october.md](notes/map-of-statements-14-october.md),
with a header that lists these changes.

## What is still open, with a first experiment for each

1. **Can lifts keep up with the bar?** This is the last named certificate route to a saving of 0.1. Shown (by hand,
   with exact LP certificates; two skeptics): excess per edge at least 0.075 from block 4, and 18 ≤ b_4 ≤ 36. If every block-4 lift needs 36, lifts miss the threshold at
   block 4, and if their excess per edge grows faster than m^0.1 − 1, no certificate family I know can decide the
   bar, and proof effort on the bar should move to constructions.
   *First experiment:* show that "defect collisions" never pay at block 4, which looks like a finite check; it would
   give b_4 = 36 for every lift. In parallel, a rank-deficiency version of Sylvester–Gallai (about 2-flats rather than
   lines; Theorem 1.14 of Dvir, Saraf and Wigderson), tested first at block 5.
2. **Can any certificate family make its angle factor grow like m^0.9/2?** Every family computed stays below 2.
   *First experiment:* build block-5 lifts optimised for small angles between far frames (the depth-3 average sin²),
   not for small edge excess, and measure their angle factor with `programs/angle-law/` against the start lifts.
3. **The block-4 near miss, certificate side** (2 in 7, saving 0.096, a gauge for 3 in 10 at 0.132 and 4 in 13 at
   0.150). No certificate searched has depth-3 average below 0.7 times its edge average at full spectral factor.
   *First experiment:* compute the cuts of the Sp(8,2) irreducibles of dimensions 1071, 1275, 1360 and 1785 with
   `programs/near-miss-certificate-search/sp8lib.py`, and score each with the corrected test; then overlap cuts whose
   first set has 4 or more elements, and whole-graph lifts built for depth 3 (a mean edge excess below 0.2 would
   suffice).
4. **The chi-game cluster: 3 in 8 (0.107), 2 in 5 (0.166), 3 in 7 (0.229).** *First experiment:* have a skeptic check
   the derived cut cap 28/11; then build the irreducible pieces of the unitary group that acts on the 891 frames and
   list their cuts, the way section 4 did for Sp(6,2). In parallel, try to prove that an interval lift with o bad edges
   at its two ends has excess at least 18 + o/6, which would mean no lift ever excludes 3 in 8.
5. **The last symmetric block-4 skeleton, g27_Fany4_s7.** 15 charts open. *First experiment:* run the exact chart test
   on each open chart with Singular over the rationals (5 to 30 minutes per chart on a Colab machine); the chart list
   is `data/symmetric-block4-skeletons/g27_Fany4_s7_open_charts.txt`.
6. **3 in 10 at block 4 by letter choice** (0.132, above the bar). *First experiment:* make the junk count include the
   other arrays' target rows, rerun the descent from the one aligned 3-in-10 skeleton and its neighbours, and screen
   each aligned skeleton with every coefficient free; also redo the letter-choice verdicts with first coefficients
   allowed to be 0, which takes seconds per program.
7. **Approximate 2-in-7 units on the four border skeletons.** Their full systems are undecided. *First experiment:*
   find the numerical curve of f27_FAmid3_m3_s16 again and save it, use it to seed a search in the full system, and run
   the other three on all eight characters with hours of msolve time (powers written as x^k) or a sparse closure.
8. **The 76 refutations modulo one prime.** *First experiment:* recheck them at a second prime, or exactly over the
   rationals with Singular, before calling them refutations.
9. **The existence pipeline on a real unit of ours.** *First experiment:* run it on the 8-move control word of the
   earlier rounds, which should give a 1 × 1 certificate cheaply, and on the plain schedule; then state the
   certificate lemma (commuting matrices have a common eigenvector) in Lean on a Colab machine.
10. **Is s\* computable?** This is the precise open logical question on the lower-bound side. No machine work is
    suggested for it. *First step:* check whether the arguments known for the same question about the Shannon
    capacity of graphs (Boche and Deppe) and asymptotic tensor rank (Christandl, Hoeberechts, Nieuwboer, Vrana and
    Zuiddam) carry over to our algebra of combining units.
11. **GCT couplings.** *First experiment:* a solver that prescribes the rank drop of the solution stratum, validated
    on the plain 1 in 3 inside the pieces 35 and 56 + 70 where ordinary descent fails; then decide 5 in 14 in 56 + 70.
12. **Smaller loose ends.** Find r′(3) exactly by enumerating the block-3 engine group on a Colab machine, and test
    whether f = 1 to 8 already suffices at block 2. Get second skeptics for Theorem C and Lemmas D and E, Theorems F
    and I of the logic round, Theorem E+, the 89/47 certificate, the GCT propositions and the pattern identity of
    section 13.

## Programs, data and notes

Everything here runs with Python 3.9 or later unless marked. "numpy" means the program needs numpy; I ran those with
`uv run --no-project --with numpy python PROGRAM` (scipy, mpmath and pyzx were added the same way where marked). Run
each command from the program's folder. Some programs write files into the current folder, so run them from a copy if
you want to keep the shipped outputs. [data/](data/) holds the outputs of the original runs under the same names, and,
for a few programs, a new output from 19 October (names ending in `_rerun.txt`). [notes/](notes/) holds the working
proofs and notes, kept as written with a header giving their status now. [results.json](results.json) is the
machine-readable summary: every finding with its status in the words used here, its numbers, evidence and programs;
the open questions with their first experiments; and every program with its command, time and memory.

**Re-run for this revision.** On 19 October I re-ran 38 programs (41 runs) from a clean copy of these folders on the
Mac, with the commands given below. Every one ran without error. Where an output is shipped, the new output equals it
line for line, apart from timing lines, the file paths that some programs print in their first line, the last two
digits of one floating-point LP value (the certified bound 120/7 is the same), and for `vc_1434.py` the choice of
cuts (the shipped log has cuts 1434 and 1432, the new run 1434 and its mirror 854; the numbers for 1434 agree). The
programs marked "not re-run" were not re-run for this
revision, because they need Singular, GAP or the Colab machines, or produced the large searches.

| folder | what it shows | command (from that folder) | time, memory | re-run |
|---|---|---|---|---|
| [frame-graph-eigenvalues](programs/frame-graph-eigenvalues/) | exact eigenvalues of the frame graphs from intersection numbers; cut caps 2/(1 − θ1/k): 5/14 and 13/30 in the engine game, 28/11 in the chi game | `python3 dpg_eig_check.py` | under 1 s | yes |
| [one-exponent](programs/one-exponent/) | the exact identity "one direction = two block calls" (`archimedes.py`); lift profiles and growth exponents (`liftprof.c` on `data/one-exponent/lifts/g3c.txt` and `g4c.txt`; `liftsample.py` samples start lifts) | `python3 archimedes.py`; `cc -O2 -o liftprof liftprof.c -lm -lpthread && ./liftprof ../../data/one-exponent/lifts/g4c.txt 4`; `python3 liftsample.py 4 3000 1`; `python3 gen.py 5 g5c.txt` writes the block-5 lift | each under 2 s; `gen.py 5` 33 s, 85 MB | yes (blocks 3 and 4; `gen.py` reproduces `g3c.txt`, `g4c.txt` and the unshipped `g5c.txt` byte for byte; the block-5 and block-6 profile runs not re-run) |
| [sp62-mixed-unit-certificate](programs/sp62-mixed-unit-certificate/) | Theorem M: exact ranks of the 105-dimensional certificate (`verify_exact.py` on the saved matrices); the skeptic's independent rebuild in the Weyl group of E7 (numpy); the producers (`exact_cert.py`, `fgcore.py`, `fg_build.py`, `fg_cuts.py`, GAP scripts in `gap/`) | `python3 verify_exact.py ../../data/sp62-mixed-unit-certificate/exact_mats_lam3.json`; `python3 skeptic_e7_check.py 1` | 5 s; 8 s, 50 MB | yes; producers not re-run (GAP; `fg_build.py` writes a 228 MB cache of the irreducibles, not shipped) |
| [venn-certificate-spare-row](programs/venn-certificate-spare-row/) | the 89/47 certificate: exact cap (skeptic) and modular demands; the mixed-room computation | `python3 skeptic_exact_cap.py 3 4 0,2`; `python3 indep_venn_check.py 3 4 0,2` (numpy); `python3 mixed_room.py ../../data/venn-certificate-spare-row/out_venn_sweep_m3.jsonl` | 1 s; 37 s and 1.5 GB; under 1 s | yes |
| [chi-game-ceiling](programs/chi-game-ceiling/) | the 13-row certificate for 169/57 over F_4 (two independent checkers) and the lift of excess 18 | `python3 verify_q.py ../../data/chi-game-ceiling/certificate_q4_m3_K5_u_1D.json` (numpy); `python3 skeptic_mycert.py` with the same file (numpy); `python3 skeptic_mycheck_lift.py ../../data/chi-game-ceiling/lift_alpha_s2.txt` | 16 s; 13 s; 5 s | yes |
| [exactification](programs/exactification/) | Theorem X's toy unit composed exactly to depth 6; the skeptic's counterexample to the factor 2 | `python3 exact_toy.py`; `python3 skeptic_exactify_count.py` (numpy) | each about 1 s | yes |
| [clifford-cohomology](programs/clifford-cohomology/) | Theorem E+: the skeptic's independent checks (random Clifford automata, controls, ring restriction) and the author's library | `python3 indep.py 1 2 3 4`; `python3 ctrl.py`; `python3 ring.py` | about 2 s in all | yes (`cqca.py`, `t1_bulk.py` not re-run) |
| [lift-gluing](programs/lift-gluing/) | the counting behind E ≥ 68, the constants κ_m, the exact lift with excess 3, and the LP dual certificates for b_4 ≥ 18 and b_5 ≥ 141 | `python3 counts3.py`; `python3 kappa.py`; from `data/lift-gluing/`: `python3 ../../programs/lift-gluing/check_e3.py`; `python3 dualcert.py 4 1` (numpy, scipy) | each under 7 s | yes (the block-5 dual not re-run) |
| [angle-law](programs/angle-law/) | Theorem F exactly up to block 40; Theorem H on hooks; frame-family checks; the spectral ceiling table | `python3 omega_exact.py 40`; `python3 frob_ceiling.py`; `python3 hook_cells.py` (numpy); `python3 frames_check.py` (numpy) | under 1 s; under 1 s; 1 s; 14 s, 640 MB | yes |
| [near-miss-certificate-search](programs/near-miss-certificate-search/) | the follow-up searches (hook relatives, Sp(8,2) cuts, projections) and the skeptic's independent library | `python3 vc_1434.py 1434,854` (numpy); `python3 proj_counts.py` | 14 s, about 3 GB; under 1 s | the two named; the searches (`hookrel.py`, `mcuts.py`, `vc_run.py`, `vc_ab.py`) not re-run |
| [zx-sum-free](programs/zx-sum-free/) | the ZX track's PyZX tests and the skeptic's checks, including the PyZX bug | `python3 skeptic_pyzx_bug_repro.py` (numpy, pyzx 0.10.7) | 3 s | the bug reproduction only; `zx1`–`zx4` ran on a Colab machine and were not re-run |
| [symmetric-block4-skeletons](programs/symmetric-block4-skeletons/) | the skeptic's own rebuild of the 652 symmetric skeletons with every coefficient free, its Singular driver and chart decision; recipes of all 678 systems | `tar -xzf ../../data/symmetric-block4-skeletons/recipes_all678.tar.gz`, then `python3 vnr_run.py relax NAME`, or `python3 vnr_run.py fz NAME` followed by `python3 chart.py NAME 1800` (numpy; needs Singular on the path) | seconds to 30 minutes per skeleton | not re-run (Singular) |
| [every-f-and-logical-form](programs/every-f-and-logical-form/) | the finite engine groups and r′ (author and skeptic), why twiddle factors are not finite-state, the widened check, the transfer counterexample and the infinite controlled-S element, block-3 traces | `python3 fbound.py` (numpy); `python3 engine_group_exact.py`; `python3 transfer_and_cs.py` (numpy); `python3 hankel.py`; `python3 block3_traces.py`; `python3 widened_check.py` (all numpy except the second) | 3 s; 90 s, 200 MB; each other under 2 s | yes (`widened_check.py` not re-run) |
| [two-certificate-theorem](programs/two-certificate-theorem/) | the standalone checker of commuting-matrix certificates; the skeptic's one-prime Hensel certificate and points-at-infinity check | `python3 check_cert.py ../../data/two-certificate-theorem/certs/*.json`; `python3 skeptic_hensel_demo.py`; `python3 skeptic_infinity_verify.py` | 4 s; each under 1 s | yes |
| [border-curves](programs/border-curves/) | the skeptic's 80-digit refinement of the g27_Fp2b_m3_s5 border curve; its parser and search | `python3 refine.py` (numpy, mpmath) | 12 s | yes |

**Not in this revision.** The programs of the GCT track, of the letter-choice follow-up, and of the original
algebraic-geometry track and its re-run are long multi-folder pipelines that import the private builders of the
construction rounds; they could not be made runnable here in the time I had, so the findings that rest on them are
marked above with their status, and the skeptic's self-contained rebuild stands in for the skeleton verdicts. The
sweeps of the logic round (1,304 relaxations, the 133-skeleton extension) and the existence pipeline's search are not
included either; their certificates and verdict tables are. The 228 MB cache of Sp(6,2) irreducibles and a 107 MB
border-track cache are regenerable and not shipped, nor is the 36.6 MB block-5 lift `g5c.txt` behind the block-5
profile; `programs/one-exponent/gen.py 5 g5c.txt` writes it again in about 33 seconds (SHA-256 in the Sources). The only msolve output from before the bug of section 13 is the relaxation table,
which the bug cannot affect. Raw agent logs,
machine job folders and the round briefs are not published.

## Credits

- The frame-graph eigenvalues are the classical formula for dual polar graphs (Brouwer, Cohen and Neumaier,
  *Distance-Regular Graphs*, 1989), cited from memory and checked by computer up to block 40.
- The character table of Sp(6,2) and its irreducibles were computed with GAP 4.12.1 and T. Breuer's character table
  library; the ATLAS naming (105a) is the library's.
- The algebra of section 1 follows T. Fritz, "A generalization of Strassen's Positivstellensatz", arXiv:1810.08667
  (Definitions 2.4 and 2.8, Theorem 2.12), and "Resource convertibility and ordered commutative monoids",
  arXiv:1504.03661 (abstract only), and V. Strassen's asymptotic spectra (1988, through Fritz).
- The inequality behind the ceilings is Theorem 8.6, inequality (65), of R. Dougherty, C. Freiling and K. Zeger,
  arXiv:1311.4601; Ingleton's inequality is classical. Hirzebruch's inequality was used as stated by T. Bauer et al.,
  arXiv:1407.2966, Theorem 3.4.
- The Sylvester–Gallai transfer uses the rank bound for design matrices of Z. Dvir, S. Saraf and A. Wigderson,
  arXiv:1211.0330.
- Approximate units: M. Bläser, "Fast Matrix Multiplication", Theory of Computing Graduate Surveys 5 (2013), section
  6, for Bini's and Schönhage's methods (Bini, Capovani, Lotti and Romani 1979; Schönhage 1981, cited from memory).
- Barrier literature read as abstracts or first-hand by the tracks: Christandl, Vrana and Zuiddam (arXiv:1812.06952,
  arXiv:1709.07851); Alman (arXiv:1812.08731); Alman and Vassilevska Williams (arXiv:1810.08671); Ambainis, Filmus and
  Le Gall (arXiv:1411.5414); Efremenko, Garg, Oliveira and Wigderson (arXiv:1710.09502); Ailon (arXiv:1305.4745,
  quoting Morgenstern 1973).
- ZX: PyZX 0.10.7 (A. Kissinger, J. van de Wetering and contributors); the ZX papers whose abstracts the track read are
  arXiv:1902.03178, 1903.10477, 2004.05164, 2302.12135, 1705.11151 and 1905.00041.
- GCT: J. Grochow, arXiv:1304.6333, section 3.2, for the fact that a closure-type proof exists exactly when no
  approximate unit sits on the skeleton.
- Logic: J. Kollár, "Sharp effective Nullstellensatz", J. AMS 1 (1988); Heintz's Bézout inequality as quoted by
  M. Sombra, alg-geom/9610006; D'Andrea, Krick and Sombra, arXiv:1103.4561; Laurent and Mourrain, arXiv:0812.2563; the
  Stanford Encyclopedia of Philosophy entry "Provability logic"; Löb's and Solovay's theorems and the Lefschetz
  principle from memory; Christandl, Hoeberechts, Nieuwboer, Vrana and Zuiddam, arXiv:2411.15789, and Boche and Deppe,
  arXiv:2010.06873, for the parallel open questions. The barriers of section 15 that do not reach our games, read as
  abstracts: algebraic natural proofs (Grochow, Kumar, Saks and Saraf, arXiv:1701.01717; Forbes, Shpilka and Volk,
  arXiv:1701.05328) and unprovability in bounded arithmetic (Razborov, Izvestiya Math. 59, 1995; Pich and Santhanam,
  STOC 2021).
- Tools: Singular, msolve 0.6.5 (with the parsing bug of section 13), HiGHS through scipy, mpmath.
- OpenAI's challenge statement and preprint define the finite win and the cost model these rounds read against.
- The certificate method rests on a pull request by eumemic to the community's hub repository (10 October 2026, prepared
  with AI assistance by its own statement; listed in the Sources). It has the certificate lemma (its "FO lemma") of
  which Theorem R is a general form, the reading of schedules as codes on the frame graph (the dual polar graph), and the
  statement that gates without a (1 + i) denominator never save, which Proposition G of section 1 extends to every
  composition.
- Theorem R in its general form, the hook, the frame certificates, the 2-adic bound and the two-layer bound are my
  earlier work; some of it is in the earlier revisions of this repository and the rest is in the sections linked
  above.

## Sources

The hub pull request of the Credits: #288 of `CrocSwap/integer-mult-bounds`, opened by eumemic on 10 October 2026.

Working names and dates. The rounds ran on 14 and 15 October 2026 (Eastern time). The names below are the working
folders of the rounds and tracks; the programs' comments and the notes refer to each other by them. The full map from
every published file to its source folder is kept privately.

| folder here | round and track (working name) |
|---|---|
| `programs/frame-graph-eigenvalues/` | theory-of-theories round, the summary's own check (`meta/note-derived`) |
| `programs/one-exponent/`, `notes/one-exponent-proof.md` | theory-of-theories round, spectrum track (`meta/spectrum`); lifts `g3c`, `g4c` and their generator `gen.py` from wave 16 (`wave16/b56-methods`) |
| `programs/sp62-mixed-unit-certificate/`, `notes/sp62-certificate-and-cut-cap-proof.md` | second cycle, finite-group track (`meta/c2-finite-group`, skeptic `skeptic-logic/`, compute skeptic `meta/v-meta2-finite-group-compute`) |
| `notes/certificates-by-game-proof.md` | first cycle, representation track (`meta/reps`; skeptic `meta/v-meta-reps`) |
| `programs/venn-certificate-spare-row/`, `notes/transplants-and-venn-certificate.md` | first cycle, literature track (`meta/literature`; skeptic `meta/v-meta-literature`) |
| `programs/chi-game-ceiling/`, `notes/chi-game-ceiling-proof.md` | second cycle, Hermitian track (`meta/c2-hermitian`; skeptics `meta/v-meta2-hermitian-logic`, `meta/v-meta2-hermitian-compute`) |
| `programs/exactification/`, `notes/exactification.md` | first cycle, border track (`meta/border`; skeptic `meta/v-meta-border`) |
| `programs/clifford-cohomology/`, `notes/clifford-cohomology-proof.md` | first cycle, cohomology track (`meta/cohomology`; skeptic `meta/v-meta-cohomology`) |
| `programs/lift-gluing/`, `notes/lift-gluing-and-sylvester-gallai-proof.md` | second wave, lifts track (`meta2/a2-sg-lifts`; skeptic `meta2/v-a2-sg-lifts-compute`); first-cycle lifts track `meta/c2-lifts` |
| `programs/angle-law/`, `notes/angle-law-proof.md` | second wave, stability track (`meta2/a2-vi-stability`; skeptic `meta2/v-a2-vi-stability-compute`) |
| `programs/near-miss-certificate-search/` | follow-up `meta2/cert-07`; its skeptic `meta2/v-cert-07` |
| `notes/gct-proof.md` | second wave, GCT track (`meta2/a2-gct`) |
| `programs/zx-sum-free/`, `notes/zx-proof.md` | second wave, ZX track (`meta2/a2-zx`; skeptic `meta2/v-a2-zx-compute`) |
| `programs/symmetric-block4-skeletons/`, `notes/symmetric-skeletons-pattern-certificates.md` | second wave, algebraic-geometry track (`meta2/a2-nag`); follow-ups `meta2/v-nag`, `meta2/nag-rerun` (verdict tables, recipes) and its skeptic `meta2/v-nag-rerun` (programs); letter choice `meta2/dsp-build` |
| `programs/every-f-and-logical-form/`, `notes/logic-classification-*.md` | logic round, `meta3/ml-logic` and skeptic `meta3/v-ml-logic` |
| `programs/border-curves/`, `notes/proof-complexity-proof.md` | logic round, `meta3/ml-proofcx`; programs from its skeptic `meta3/v-ml-proofcx` |
| `programs/two-certificate-theorem/`, `notes/two-certificates-*.md` | logic round, `meta3/ml-fooling` and skeptic `meta3/v-ml-fooling` |
| `notes/map-of-statements-14-october.md` | the round's map (`meta/MAP.md`), Appendix A and two bookkeeping items removed |

The round summaries behind this README are the three round notes (`meta/NOTE.md`, `meta2/NOTE.md` with its follow-ups
section, `meta3/NOTE.md`) and the follow-up notices (`meta2/FOLLOWUP.md`: the normalisation hole and the msolve bug).
Entries of my private results catalogue that this section covers: M1 (one dimension), T5 (growth exponents), T6 (cuts
and colourings), T8 (lifts), T9 (angle law, chi-game ceilings), T13 (relaxation blind spots), T14 (cheap
refutations), L3 (one array needs m moves), the cohomology part of L18 (Theorem E+), the Theorem X part of L14, X8
(Theorem M), X12 (the Venn certificate), X25 (symmetric skeletons), X26 (letter choice), N2 (border curves), M9 to M12
(logic) and the refutations C11, C13, C14, C15 and C16. The final chart count of g27_Fany4_s7 (49 of 64) is from the
re-run's last log entry of 15 October, 01:37; the catalogue's 45 of 64 is the earlier count.

Not shipped, by SHA-256: the block-5 lift `g5c.txt` (36.6 MB) used for the block-5 lift profile, regenerable with
`programs/one-exponent/gen.py 5 g5c.txt` (checked on 19 October):
`3d13df4a99c6f92b876aa7d4703872b8f2631e4ba0f548e33e692cb1d1311a5c`; the Sp(6,2) irreducibles cache `irreps.npz`
(228 MB), regenerable with `programs/sp62-mixed-unit-certificate/fg_build.py`:
`1f2197d205c57d6f67db0141c689960e88cfc6840eab9faa5193ea4bea82c1bf`; the border track's cache `algcache_m3_VV.npy`
(107 MB): `506e627c3ad5c0e8c8ea57db8862f617a4cbdc0621b7f33a54a1004007867151`.
