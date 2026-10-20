# How savings happen: crossing slack, frame packing, the frame game, and Proposition 6

This section covers four research rounds of 15 October 2026. Together they asked how a unit saves a move at all.
Every known unit that saves anything does it by *crossing slack*: OpenAI's unit, the community's units and mine. I
wanted to know whether that is the only way. If it were, could that rule exclude large parts of the open
targets? The rounds also asked how much one paid move can carry, whether a much smaller "frame game" holds
the whole problem, and what a unit with only two arrays would have to do to save.

I did this work with a team of AI agents (Claude) that I directed. I set the goals, chose which ideas to pursue and
decided at each stage what to do next. The agents ran the searches, wrote the proofs and programs, and tried to
refute one another's work. Nothing in this section is checked in Lean. No person has read any of these proofs. Each
result was checked by at most one independent AI skeptic, not the two that the rounds' own rule asked for. Where I
say "a skeptic", I mean one agent that tried to break the result with its own reasoning and its own code. It is not a
human referee.

**In one paragraph.** Every saving of every unit is *frame packing*, paid moves that carry data sitting at several
places at once, as long as that packing is not paid back later. Crossing slack is the same thing seen in other
coordinates. No second mechanism turned up among the known units. Whether every saver can be written in the crossing
form is open, and that open question is exactly where the moonshot targets live. Packing itself is cheap and
plentiful. The record units at blocks 45 and 73 carry an estimated 5 to 7 extra labels per paid move. Almost all of
it is paid back, and what survives is the saving itself, about 0.005 moves per paid move at best. One move can carry any amount. What limits it is that
the terms of one array move together. The frame game, the engine with its Pauli part dropped, turned out to be the
whole problem: every program there completes to a genuine unit. A two-array saver would have to do three awkward
things at once. Does any two-array unit save? That is Proposition 6 for two arrays, and it is open. I built 10,529
designs aimed at two-array and three-array savers, and none of them delivered its outputs. No target shape was newly
excluded at any block.

## Contents

- [The question](#the-question)
- [Words used here](#words-used-here)
- [What I tried](#what-i-tried)
- [What I found](#what-i-found)
  1. [Every saving is frame packing that is not paid back](#1-every-saving-is-frame-packing-that-is-not-paid-back)
  2. [Parity prunes nothing](#2-parity-prunes-nothing)
  3. [Common-frame designs: exact accounting and what it caps](#3-common-frame-designs-exact-accounting-and-what-it-caps)
  4. [The refund identity: packing is plentiful and almost all of it is paid back](#4-the-refund-identity-packing-is-plentiful-and-almost-all-of-it-is-paid-back)
  5. [How much one paid move can carry](#5-how-much-one-paid-move-can-carry)
  6. [The frame game is the whole problem](#6-the-frame-game-is-the-whole-problem)
  7. [Theorem C: no crossing, no saving](#7-theorem-c-no-crossing-no-saving)
  8. [What a two-array saver must do](#8-what-a-two-array-saver-must-do)
  9. [Constructions that failed](#9-constructions-that-failed)
  10. [Claims refuted along the way](#10-claims-refuted-along-the-way)
  11. [Warnings for anyone who searches](#11-warnings-for-anyone-who-searches)
- [What is still open, with a first experiment for each](#what-is-still-open-with-a-first-experiment-for-each)
- [Programs, data and proofs](#programs-data-and-proofs)
- [Credits](#credits)
- [Sources](#sources)

## The question

A *unit* "W in M" at *block* m is a fixed program that applies the m-direction transform to W arrays using M paid
moves. The plain schedule needs W·m. A unit with fewer moves gives the saving ln(W·m/M)/ln m to the whole algorithm.
"2 in 7 at block 4", for example, would transform 2 arrays with 7 paid moves instead of 8. That is worth 0.096.

The anatomy round of 15 October (in [../why-the-savers-save/](../why-the-savers-save/)) found that OpenAI's unit, the
community's T74 and my U8X0HQ all save in the same way. Each array starts at one input's "source frame" and ends
holding a different input's transform. When the two frames overlap, the array gets there in fewer than m moves. That
round wrote the count as M = W·m − 2N + 2L. Here N is the sum of the source labels' dimensions and L counts steps
back down. This section asks four questions about that picture:

1. **Is crossing slack the only way to save a move?** If it were, could the identity rule out large parts of the
   moonshot targets? Could it settle the hope that the number of saved moves is always even?
2. **How much can one paid move carry?** Does the answer give a ceiling at each block?
3. **Is the frame game the whole problem?** The frame game is the engine with the Pauli part dropped, and it is far
   smaller. Can every search run there?
4. **Must two-array savers pack?** What exactly must a two-array unit do to save? Can it be proved that none does?
   That last question is Proposition 6 for two arrays.

The open two-array targets are 2 in 7 at block 4 (worth 0.096), 2 in 7 and 2 in 8 at block 5 (0.222 and 0.139),
and 2 in 8, 2 in 9 and 2 in 10 at block 6 (0.226, 0.161 and 0.102). The pair targets of the stepping-stone search
are 2 in 9 at block 5 (0.065) and 2 in 11 at block 6 (0.049). Smaller two-array shapes are excluded: 2 in 6 at
block 4 by my block-4 certificate (41/13 moves per array), and 2 in 6 at block 5 and 2 in 7 at block 6 by the
two-layer bound (every unit with two or more arrays needs at least W + m moves). Both are described in
[../block-bounds-and-new-rules/](../block-bounds-and-new-rules/).

## Words used here

- **Array, paid move, free step.** An array is one vector being transformed. A paid move applies one direction's
  2-point transform to one array and is what we count. Free steps cost nothing in the count: relabelling addresses,
  phases, bit flips and sign flips, and mixing arrays with number coefficients.
- **Term.** An array's content is a sum of terms. Each term is one input passed through one fixed product of steps,
  times a number. A unit must be right at every number of digit columns, so it must be right as such a formal sum,
  and its terms are well defined.
- **Digit columns f.** A unit must be correct for every number f of digit columns. The proofs split this by f mod 4.
- **Frame, frame graph, level.** Each term sits at a point (a *frame*) of a finite graph, the frame graph. Its points
  are the Lagrangian subspaces of a 2m-dimensional space over the two-element field. There are 15, 135 and 2,295 of
  them at blocks 2, 3 and 4. Every input starts at the *raw frame*, and its transform sits at distance m from it.
  Free steps move nothing. A paid move takes every term of its array one step: up, down (a *descent*) or sideways.
  A term's *level* is its distance from the raw frame.
- **Coordinates (gauge).** A bookkeeping choice of one reference point per input. Frames are then read relative to
  it. OpenAI's source labels are such a choice. In *raw coordinates* each input is read from its own raw start.
  Changing coordinates changes every bookkeeping number except M.
- **Common-frame design.** A design in which, in some coordinates, every array holds data at only one frame at every
  moment, so additions join only arrays at the same frame. OpenAI's unit and the community's units are common-frame
  as designed.
- **Crossing slack.** In a common-frame design, the overlap between an array's start frame and its end frame. That
  overlap is what such a design saves.
- **Frame packing.** A paid move on an array whose terms sit at two or more frames, so one paid move advances all of
  them. *Live* packing counts only terms whose data reach an output. *Irreducible* packing P* (for two arrays) is the
  least live packing over all choices of coordinates.
- **Refund.** Packing that is paid back later, by additions that cancel terms, by descents or by sideways steps.
- **Same-input split.** One array holds two different terms of the same input at the same moment.
- **Scratch array.** Any array beyond the W data arrays.
- **Height system.** One function per input on the frame graph that is 0 at the raw frame and m at the transformed
  frame, and changes by at most one per step. The level is one example.
- **Frame game.** The same game with the Pauli part (bit flips, sign flips, the scalar i) dropped. Each term is then a
  2m × 2m binary matrix instead of a 2^m × 2^m complex one, and the free group is about 4^m times smaller.
- **Certificates (Theorem R).** My lower bounds on moves per array, from the earlier rounds: 2 and 2.8 at blocks 2
  and 3 (proved in Lean), 41/13 = 3.154 at block 4 (checked by computer), and 2.992 and 2.973 at blocks 5 and 6 and
  2.9447 from block 8 on (argued by hand, not checked; two AI skeptics). They cap the saving at 0, 0.063, 0.171,
  0.319 and 0.392 at blocks 2 to 6, and at 0.48 at block 8. They are described in
  [../block-bounds-and-new-rules/](../block-bounds-and-new-rules/).
- **Skeptic.** An independent AI agent that tried to refute a result with its own reasoning and its own code.

## What I tried

All four rounds ran on 15 October 2026. Each one paired builders or provers with skeptics.

- **The slack round.** Two agents attacked the first question in parallel, one with group algebra and one
  combinatorially. One skeptic then checked both step by step with its own exact code. I also gave the round the
  parity question: is the number of saved moves always even?
- **The packing round, first part.** Three tracks. One counted the packing in every known unit. One worked out what
  packing can buy and how much one move can carry. One tried to build small two-array units designed around
  packing. Each had its own skeptic.
- **The packing round, second part.** Two tracks. One tried closed cancellation patterns built from orbits of a
  free symmetry. The other asked whether every frame-game program can be completed to a genuine unit. Each had its
  own skeptic.
- **The Proposition 6 round.** The slack round ended with a conjecture. No unit costs less than the cheapest
  common-frame design with the same number of arrays. One agent tried to prove it and another tried to refute it,
  each with its own skeptic.

The computations are exact: rational and Gaussian-rational arithmetic, exact group elements, exhaustive state
searches. The exceptions are the searches marked as numerical. Most ran on a Mac in seconds. The largest exhaustive
searches ran on Linux machines of Google's Colab service.

## What I found

Each finding carries its status in words. "Argued by hand, not checked" means a written proof that an independent
AI skeptic re-derived without finding an error, and that no person and no proof checker has checked. "Checked by
computer" means an exact computation. I say whether a skeptic repeated it with its own code.

### 1. Every saving is frame packing that is not paid back

**The ledger.** Take any engine unit and any coordinates. Its moves satisfy

  M = W·m − 2N + X + 2D + Z + B − P   (plus one for each move on an empty array).

N is the sum of the source levels. X measures how far the sources sit off the straight path from raw to
transformed. D counts descents, Z sideways steps, B branching (level destroyed when additions cancel terms, net of
copies) and P packing. In raw coordinates N = X = 0, so

  moves saved = P − B − 2D − Z.

So every saving of every unit is packing that is not refunded. A change of coordinates moves saving between the 2N
term and the P term. Crossing slack is packing that some change of coordinates removes. They are one mechanism seen in
two coordinate systems.

- Status: **argued by hand, not checked** (one skeptic re-derived it and found no error). It is also **checked by
  computer** with no exception on 34 named units, 630 random units, 1,280 random words, every word up to length 4
  at blocks 2 and 3 (about 2.9 million runs), and a census of small two-move units at block 1. That check is one
  program. For this revision I re-ran its block-2 part, 692,402 runs: see
  [Programs, data and proofs](#programs-data-and-proofs).
- The anatomy's form M = W·m − 2N + 2L is the special case X = Z = B = P = 0. It is exact only for common-frame
  designs whose reference frame lies between every input's start and end frames (true for OpenAI's and the
  community's designs). Outside that class it is bookkeeping and excludes nothing. Its N is not even intrinsic: the
  plain schedule can be written with N = 0 or with N = W·m. **The claim that the anatomy's identity holds for every
  unit is refuted** by this ledger and by the units below.
- **No second way to save among the known units** (the skeptic agrees). Several known units break the anatomy's form,
  though:
  - My first published unit (block 80) uses 747,925 moves instead of 748,960. It saves 1,035, an odd number, which
    the anatomy's form cannot produce. The odd part is a cost of its copies.
  - My five-stage Lean units (the block-110 unit and the E8 unit at block 45) re-frame their helpers between stages.
    At block 110, 11,550 of the 51,781 additions in each invocation join arrays at different frames. In Lean's own
    coordinates they are therefore outside the common-frame class. Whether some other coordinates make them
    common-frame is unsettled. The skeptic expects not, but that is argued, not proved. Either way, the deficit they
    carry is still the data arrays' crossing slack minus copy costs.
  - The order-3 control word of the moonshot rounds is common-frame in no coordinates. It saves nothing, because
    its packing is refunded exactly at its readout.
- **What this leaves open:** must every saver be common-frame? If yes, its packing can always be rewritten as
  crossing slack. If no, there is a second mechanism. The moonshot targets live exactly in that gap.

### 2. Parity prunes nothing

The hope was that the number of saved moves is always even. If it were, every "save one move" target would be
impossible, and so would the full-sharing family U_m whenever W·m is even (4 in 15, 8 in 31, 16 in 63).

- **Refuted** by my own block-80 unit, which saves 1,035 moves, and by small exact units. The count is even only
  when all four of these hold: the unit is common-frame in some coordinates; there is no branching; there are no
  sideways steps and every source sits on the straight path; and no move acts on an empty array. Drop any one and
  odd counts occur. A one-array unit at block 1 with one sideways step uses 2 moves instead of 1.
- Even inside the common-frame class, parity is not forced, because the frame graph has triangles. A unit can save
  two moves and spend one on a sideways detour, saving one in total.
- Status of the four-condition statement: **argued by hand, not checked** (one skeptic, who also noted the fourth
  condition, which the track had missed). The odd counts are **checked by computer**.
- Consequence: the stepping-stone searches need not favour "save two moves" targets. U_m stays in play.

### 3. Common-frame designs: exact accounting and what it caps

These results are about common-frame designs only. They bound every unit only if Proposition 6 holds (see
[What is still open](#what-is-still-open-with-a-first-experiment-for-each)).

- **Exact accounting.** In a common-frame design, moves saved = the sum over arrays of the overlap of start and end
  frames, minus detours. An array that ends holding its own input saves nothing, so pure routing saves nothing.
  **Argued by hand, not checked** (one skeptic). Checked on every subspace frame up to block 4, and on OpenAI's family
  by re-running its brute-force check.
- **A general bound for every unit.** Moves saved ≤ the sum over paid moves of (number of inputs whose data that
  move carries − 1). It holds in the engine and widened games, but by itself it excludes nothing. **Argued by hand,
  not checked** (one skeptic).
- **Median bound.** The deficit is at most twice the summed distance from the best single frame to the inputs'
  start frames (or to their end frames). Designs whose sources cluster save little, and line labels, the known
  savers' choice, cluster one step from the raw frame. **Argued by hand, not checked** (one skeptic).
- **Clustered sources or sinks never save.** A common-frame design whose inputs all start at one frame, or all end
  at one frame, needs at least W·m moves, for any number of arrays. The proof is max-flow min-cut. **Argued by hand,
  not checked** (one skeptic).
- **Cuts are weak beyond two arrays.** The cut condition gives at best 1.5, 1.75 and 1.875 moves per array at blocks
  2, 3 and 4. The truth at block 2 is 2, and my certificates give 2.8 at block 3. So cut and flow arguments cannot
  replace the certificates. **Checked by computer** (exhaustive over all 16,383 frame sets at block 2; the skeptic
  recomputed 3/2, 7/4 and 15/8 with its own code).
- **Any number of arrays.** A common-frame design needs at least W·m/(C·log 2W) moves, from a published flow–cut
  gap. The constant C is unknown, so this is harmless at our blocks. **Argued by hand, not checked** (one skeptic;
  the published theorem it uses was cited from memory).
- **Network coding.** A common-frame unit is exactly a scalar linear network code for W unicast sessions on the
  frame graph, each session joining two opposite frames. A saving is a coding gain over routing. My certificates
  bound that gain: none on the block-2 graph, at most 15/14 on the block-3 graph. **Argued by hand, not checked**:
  this is a remark that the skeptic did not check separately. I make no claim about any conjecture in that field.
  The network-coding round of the same day is in [../other-fields/](../other-fields/).
- **Designs built like the known savers cannot reach 0.1.** These designs have one frame per array, equal-frame
  additions and line labels. Even with free bridges, a design with overlap one per array (line labels, every known
  design) reaches 0.1 only up to block 6 and 0.01 only up to block 29. Overlap two reaches them up to blocks 9 and
  51, and overlap four up to blocks 16 and 90. These are move counts. **Argued by hand, not checked**; the skeptic
  recomputed every figure. With the anatomy's argued walls (index-code centres need windows of 7 or 9, stages need
  block 25 or more), 0.1 is out for line designs, and 0.01 needs a block between 25 and 29 with nearly free bridges.
  Those walls are argued, not proved.
- **The full-sharing family U_m** (2^(m−2) arrays in 2^m − 1 moves, one move per direction shared by every array,
  about 4 moves per array) would need sources spread out over the frame graph, with label dimensions that grow with
  m, unlike the lines of every known saver. **Nothing proved excludes it.** That is open.

### 4. The refund identity: packing is plentiful and almost all of it is paid back

**The refund identity.** Follow each label (a term up to bit flips, sign flips and scalars) from its creation to its
destruction or to the end. Give it the paid moves it rides, its descents and sideways steps, and its levels at
creation and at the end. The ledger then splits label by label. A label that is created at an addition and never
delivered contributes exactly zero to P − B − 2D − Z, however many moves it rides. A delivered label contributes its
final level. So packed junk is always paid back in full, and the amount of packing means nothing by itself. Status:
**argued by hand, not checked** (the skeptic checked it line by line against the ledger's definitions).

**How much the real savers carry** (in raw coordinates, measured with exact fractions):

| unit | block | extra labels per paid move | surviving: moves saved per paid move |
|---|---|---|---|
| E8 unit (mine) | 45 | 5.1 to 5.7 | 0.0021 |
| U8X0HQ (mine, on the community's U8 circuit) | 73 | 6.2 to 6.8 | 0.0047 |
| T74, per replica (the community's) | 73 | 6.0 to 6.6 | 0.0049 |
| block-110 unit (mine) | 110 | about 55 | 0.0019 |
| OpenAI's unit, at the size in its paper | 10^6 | 0.0034 | 3.7e-12 |

- The third column is an estimate, except for OpenAI's unit, where it is exact. The units' global words were not
  replayed. Each call was taken as one shear of a three-shear swap, in equal shares. Its range runs from helpers
  that start empty to helpers that hold old contents. The fourth column is exact: it is each unit's own saving per
  move.
- Per call the counts are **checked by computer**: exact replays of one call of each local word. The skeptic
  repeated every one with its own code and got the same numbers to the unit. One paid move inside a call of the
  community's words and mine advances 8.6 to 12.3 labels besides its own, and about 102 in the circuit of my
  block-110 unit. Additions in the same call destroy 79 to 95 % of it again, and descents take most of the rest. In
  OpenAI's whole unit the additions destroy 82 % at the smallest size measured and 99.96 % at the size in its paper.
- Gross packing exceeds surviving packing about 1,300-fold on the community's U8 line, 2,500-fold for the E8 unit
  and a billion-fold for OpenAI's unit. More packing did not buy more saving: the block-110 unit packs about 10 times
  more per move than the E8 unit and saves less.
- The most on one move is 161,701 labels, on the centre wire of OpenAI's network at the size in its paper. Among the
  community's savers and mine, the most is 886 on the block-110 targets and 77 on the U8 and T74 targets.
- About 85 % of the packed labels sit on distinct frames (the skeptic got 86, 85 and 84 %).
- I re-ran the E8 call for this revision. The ledger of all its call types and helper models came out identical to
  the shipped output.
- Corrections by the skeptic, none of which changes a headline figure:
  1. "The third shear packs only the crossing" is false on the community's U8 and T74 lines. There the targets read
     each other for a while, so the third shear packs 1,680 (U8) and 1,968 (T74), not 896.
  2. The models with old helper contents are not consistent with the frame model's own rule for labels. With old
     contents at twist 0, their extra packing is 0.57 per move, not 0.7.
  3. The helpers carry most of the backward packing only in the second shear.
  4. "A ceiling on gross packing cannot cap the saving" is too strong. In every measured case B ≥ 0, so the saving is
     at most P, but B ≥ 0 is not proved in general. Any such ceiling would have to allow at least 5 to 100 per move at
     the known savers' blocks, which makes it useless there.

### 5. How much one paid move can carry

- **What packing buys, for every unit** (any number of arrays, any scratch): moves saved ≤ the sum over paid moves
  of (live frames the move raises − 1). Units in which no paid move raises two live frames never save. No paid move
  lifts two live frames off the raw frame, and none finishes two at once. **Argued by hand, not checked** (one
  skeptic).
- **Packing versus coding gain.** In any coordinates, moves saved ≤ live packing + coding gain of the network those
  coordinates define. In raw coordinates the coding gain is 0 and everything is packing. In a common-frame design's
  own coordinates the packing is 0 and everything is coding gain. **Argued by hand, not checked** (one skeptic).
- **Two arrays: M ≥ 2m − P\***, where P* is the irreducible packing. With no scratch, a two-array unit can save only
  if some array, at some moment, holds one input on two different labels. Two terms that differ only by bit and sign
  flips are not enough (Proposition D). **Argued by hand, not checked.** The skeptic found the reason given for one
  step wrong and supplied a correct one: every step is an invertible 2 × 2 matrix.
- **Capacity is unbounded.** One paid move can serve any finite set of steps with suitable contents. Three frames
  one step below the top can block every move, but two never can. So capacity gives no ceiling at all. **Argued by
  hand, not checked.**
- **Rigidity is the only limit found.** The terms of one array move together, and their distances never change while
  they share the array. The average a packed pair can add to the saving over its time together is **checked by
  computer**:

  | block | relative positions | average | positions that can carry a full m |
  |---|---|---|---|
  | 2 | 120 (all) | 1.18 | 1 |
  | 3 | 8,640 (all) | 1.77 | 1 |
  | 4 | 300 (a sample) | 2.31 | none in the sample |

  The one position that can carry a full m fixes both the raw and the transformed frame. The order-3 controls use it
  and pay all of it back.
- **Every ceiling is a refund ceiling.** This step is **argued by hand, not checked**; even the track that made it
  calls it an argument, not a proof (the skeptic recomputed the figures). Read
  through my certificates, at most this many net saved moves per paid move can survive. A unit with saving s needs
  at least m^s − 1 of them:

  | block | 2 | 3 | 4 | 5 | 6 | 8 | 16 | 73 |
  |---|---|---|---|---|---|---|---|---|
  | most that can survive (from the certificates) | 0 | 0.071 | 0.268 | 0.671 | 1.018 | 1.717 | 4.43 | 23.8 |
  | needed for 0.1 | excluded | excluded | 0.149 | 0.175 | 0.196 | 0.231 | 0.320 | 0.536 |
  | needed for 0.01 | excluded | 0.011 | 0.014 | 0.016 | 0.018 | 0.021 | 0.028 | 0.044 |

  The best known object has 0.0049 per move, at block 73. That is about 9 times short of 0.01 and about 110 times
  short of 0.1. These are move counts. They leave out the community's grouped pricing.
- Corrections by the skeptic: counting packing by image frames is not always the sharpest count (an exact block-2
  unit counts more sharply by labels); one parenthetical of a corollary holds only in raw coordinates. The skeptic
  also said the two-array table had left out 2 in 7 at block 6. **That correction is itself mistaken, and I am
  correcting it here.** 2 in 7 at block 6 is already excluded by the two-layer bound (at least W + m = 8 moves),
  which had been upheld by 14 October. The skeptic took the target from an older map. So 2 in 7 at block 6 is not among
  the open targets below.

### 6. The frame game is the whole problem

- **Both directions.** Every engine unit gives a frame-game program, which was known before. Every frame-game program
  completes to an engine unit: one explicit rule turns it, step by step, into an engine word. That word has the same
  arrays, the same paid moves in the same order and the same data flow, and it is right at every number of digit
  columns at every block, with no extra condition. Before each paid move the completed word applies one fixed free
  step of six terms. Free letters get corrections of 1 to 8 terms when address maps are written as elementary
  steps, and up to 2·4^m terms for a general address map. Nothing paid is added. Status: **argued by hand, not
  checked**, with one skeptic for each half. Both halves are also **checked by computer**:
  - the track checked 34 program runs on raw arrays (f = 1 to 8 at blocks 2 and 3, 1 to 6 at block 4) and exact
    identities at block 2;
  - the skeptic built its own engine and its own lift, checked exact identities at blocks 2 and 3, and checked all
    10,800 products of the 720 frame-game elements with 15 generators at block 2. At block 2 that settles the
    rule's key property by exhaustive computation.
- **The caveat.** The completed word's multi-term free steps need one scratch array. Units may use scratch, so this
  costs nothing. But results about two arrays *without* scratch do not transfer to the frame game. Whether every
  frame-game program also has a completion without scratch is open.
- **What it changes.** Every construction search can run in the frame game and be certified there. `lift.py` writes
  out the engine word, and replaying it on raw arrays then checks the code, not existence. Engine theorems hold in
  the frame game too (2 in 3 is impossible at block 2 there). No further reduction of this kind exists from block 3
  on, because the frame game's group is simple there (a standard fact).
- Along the way, the build track found the lift back to the engine for f = 1 (mod 4), its skeptic found a second
  lift for f = 3 (mod 4), and the completion track then closed the even column counts.

### 7. Theorem C: no crossing, no saving

**Theorem C (heights).** Take any number of arrays, any scratch and any contents. Suppose some height system puts all
terms of every moved array at one height right after each move. Then the unit needs at least W·m moves. The proof
splits everything the unit can hold by height and adds up each input's coefficients. A move adds at most one
direction at one height. It needs only the f ≡ 0 (mod 4) column counts, so it holds for every unit.

- So every saver needs, for every height system, a paid move whose output array straddles two heights. This is the
  general form of "no crossing, no saving": some form of crossing is necessary for every saving. It does not say
  crossing is the only way to save. It rules out no target shape, but it filters designs.
- Quantitative forms bound M from below by sums of shortest-path distances. `heightfilter.py` finds a height system,
  or the best such bound, for any concrete word in milliseconds.
- Status: **argued by hand, not checked**. One skeptic checked it line by line and rebuilt the test units in its own
  exact engine (raw arrays f = 1 to 8 up to block 3, f = 1 to 6 at block 4). All twelve of its predictions held.

### 8. What a two-array saver must do

For two arrays the cheapest common-frame design costs 2m moves at every block. So a two-array saver must escape
every one of these results:

| result | setting | conclusion | status |
|---|---|---|---|
| Corollary 2 (slack round) | two arrays in all, scratch included, common-frame | at least 2m | argued by hand, not checked (one skeptic) |
| Corollary 2′ (slack round) | two arrays, no scratch, at most one term per input in every content | at least 2m (while two inputs are mixed, their twist stays fixed) | argued by hand, not checked (one skeptic) |
| Theorem 4.3 (slack round) | one or two data arrays, any scratch | if some coordinates have no frame packing, at least W·m (via Hu's two-commodity flow theorem) | argued by hand, not checked (one skeptic) |
| Proposition D (packing round) | two arrays, no scratch | saving needs one input on two different labels in one array | argued by hand, not checked (one skeptic) |
| Theorem A (Proposition 6 round) | two arrays, no scratch, plain mixes, one label per entry at every moment (bit and sign flips allowed) | at least 2m | argued by hand, not checked (one skeptic) |
| Theorem B and Corollary T₀ (Proposition 6 round) | two inputs, any scratch, every moved array holds at most one term of each input | with free twists, at least 2m; with any twists, a saver needs at least two different twists and some combination of them (possibly a single twist) that carries the raw frame closer than m to the transformed frame | argued by hand, not checked (one skeptic) |
| Theorem C | every unit | needs a move whose output straddles two heights, for every height system | argued by hand, not checked (one skeptic) |

**Together.** A two-array saver must, at some moment, hold two different terms of one input in one array. It must
pack in every choice of coordinates. Its packing must also survive the readout. The order-3 control word shows how
that last part fails: its readout destroys every packed term at full level, so it saves nothing. In numbers, a
two-array unit that saves k moves needs irreducible packing P* ≥ k: at least 1 for 2 in 7 at block 4; 3 and 2 for
2 in 7 and 2 in 8 at block 5; 4, 3 and 2 for 2 in 8, 2 in 9 and 2 in 10 at block 6; and 1 for each pair target.

**Lemma M** (argued by hand, not checked; one skeptic). With unlimited scratch, a unit is exactly a growing module. Start from the
free mixtures of the inputs. Each paid move adds the free mixtures of one moved vector, and the unit is done when
every transformed input lies in the module.

**Rigidity is essential** (argued by hand, not checked; one skeptic). If a pair's two halves could move
independently, two inputs would finish in 3m/2 moves. So any proof that two-array units never save must use the fact
that a pair's twist stays fixed.

### 9. Constructions that failed

No construction in this section produced a unit. Every result below is a negative result for one family of designs,
not for the targets.

**Packing-first designs (packing round, first part).** Both inputs share one paid move (two for the 2 in 2m − 2
targets; one shared by all three inputs for 3 in 10). Each input climbs alone before it and finishes alone after it.
Without extra terms such a design can never cancel its junk (argued by hand). The repair was to land stored arrays
on the junk one at a time and solve again exactly.

- 8,522 designs at 2 in 7 (block 4), 2 in 8 and 2 in 9 (block 5), 2 in 10 and 2 in 11 (block 6) and 3 in 10
  (block 4): none delivered all its outputs. None of 1,130 designs at the negative controls did either, as must be
  (2 in 5 at block 3 and 2 in 6 at block 4 are excluded by the block-3 and block-4 floors, which hold in the frame
  game). **Checked by computer**; the skeptic recounted every tally from the raw outputs.
- Each landing fixes one junk term and creates more. In one typical run the unpartnered junk went 2, 4, 7, 9, 13,
  20, 48, 64 over eight landings. Landed terms stay in every later state of their array.
- The *position test* keeps only each term's frame and asks whether the target can still be produced from the
  stored arrays with plain coefficients. With the landing coefficients left free, it refuses 71 % of the 5,271
  production designs exactly, for every choice of coefficients. A numerical search over the coefficients of 341 of
  the rest found none that passes. That search is evidence, not proof.
- **Refuted: the "closest design".** It was reported to pass the position test with special coefficients. The
  skeptic showed that it fails exactly, for every choice. The reported near-pass was a limit in which the shared
  move stopped packing and the readout coefficients grew without bound. The best genuine attempt stayed far off (a
  residual of 0.21 against 0 needed). The skeptic also found the repairs shallower than stated: 15 to 29 % of repair
  steps only rescaled an array, and real landings were at most 6 or 7 per design, not 12.

**Orbit designs (packing round, second part).** Instead of one partner at a time, each repair added a whole orbit
of a free symmetry, hoping the copies would cancel each other's junk in a closed pattern. Five kinds of symmetry
were tried.

| verdict on 2,007 designs (6,244 orbits) at 2 in 7 (block 4), 2 in 8 and 2 in 9 (block 5), 3 in 10 (block 4) | designs | share |
|---|---|---|
| refused by the linearised position test, so dead for every choice of scalars | 1,723 | 86 % |
| refused at the drawn scalars | 147 | 7 % |
| undecided (over the exact solver's budget) | 137 | 7 % |
| delivered all outputs | 0 | |

- **Checked by computer**, one skeptic, who recounted every tally with its own frame game and solver. The accurate
  statement is that no exact success was observed. At the end, 519 of the 4,584 target outputs were refused by no
  screen, and the solver ran out of budget on 3,476 output solves.
- Orbits made the junk grow faster than single landings did. After two repairs the stored arrays held about 66 to 99
  terms, against 17 to 20. Why, in part (argued): a free step never changes the frames of an array's terms, so an
  orbit's copies start on one frame. The paid moves that follow spread them out, and inside one array they stay
  distinct for good. The proof covers only copies inside one array.
- **The parity law** (argued by hand, and checked by computer on 41 cases at blocks 3 to 5 with symmetry orders 2 to
  31; the skeptic checked 40 with its own code, with no exception). The order-3 control is a *top-closed* pattern:
  both arrays end with every term on the target frame. The same pattern works with a free symmetry σ that fixes
  both the raw and the target frame exactly when the order of σ is odd. The even-order cases pass the position test
  and are still impossible. So the position test is necessary but not sufficient.
- **One shared pair gives one ratio** (argued by hand). Once two inputs start sharing an array, every top-closed array
  descended from that moment holds them in one fixed ratio. No readout of such arrays can then separate the two
  outputs. A second ratio needs a second point where sharing starts, and in the two shapes analysed that costs the
  saving back. The skeptic agrees for those two shapes.

**The single-term game (Proposition 6 round, refute track).** Take the class "two inputs, any scratch, every moved
array holds at most one term of each input". That class is an exact game on frames. Its solver was checked against
actual modules in the group algebra: the skeptic made 2,732 checks in nine games and found nothing missed and nothing
extra, and two deliberately broken versions of the game fail that test. The cube (the 2^m frames of one "apartment")
and the Hamming game sit exactly inside the frame game.

- **The block-4 cube has no single-term 2 in 7.** **Checked by computer**, exhaustively. The skeptic re-ran it
  independently on all 107.0 million six-move states. 1.44 million of them could still finish, and from none does a
  seventh move deliver both inputs. So the cube's cheapest single-term two-input program costs exactly 8.
- **No single-term 2 in 5 in 19 small games** with the target three steps away (the first size the two-layer bound
  leaves open there). **Checked by computer**, repeated by the skeptic with its own decider.
- **The key inequality held everywhere searched.** The highest level reached by input 1 plus the highest level
  reached by input 2 never exceeded the number of moves. It held in the block-3 cube and Hamming game to 6 moves (80
  million states), the block-4 cube to 6 moves (91 million) and the block-4 Hamming game to 5 moves (111 million).
  **Checked by computer.** A beam search by the skeptic in the block-4 Hamming game to 7 moves found no violation;
  that search is not exhaustive.
- **Why the obvious crossing twist does not save.** A pair whose twist swaps the raw and the transformed frame,
  branched, delivers both inputs in m moves. But creating that pair costs at least m, so the total is 2m again.
  Locked potential is cheap to build, and unlocking it costs exactly what it would save.
- **What this does not cover.** The full block-4 frame game (2,295 frames) was not searched, so a single-term 2 in 7
  at block 4 is narrowed, not excluded. Multi-term contents are not covered.

### 10. Claims refuted along the way

Failures are as useful as successes, so here they are plainly:

- **"The anatomy's identity M = W·m − 2N + 2L holds for every unit."** Refuted by the general ledger (section 1). It
  holds only for common-frame designs, under extra conditions.
- **"Saved moves are always even."** Refuted by my own block-80 unit (1,035 saved moves) and by small exact units
  (section 2).
- **"No known saving unit breaks the anatomy's identity", and "the E8 unit, U8X0HQ and T74 are all common-frame by
  construction of the community's frame templates".** Both were claimed in the slack round and refuted by its
  skeptic: my two units, the E8 unit and U8X0HQ, re-frame their helpers, and the block-80 unit's count is odd.
- **"A change of input basis makes every unit common-frame" (Proposition S).** Refuted by the prove track itself.
  The exact unit SPLIT1 (one input, m + 1 moves) holds one input at two frames in every basis and every choice of
  coordinates. It wastes a move, so it does not touch Proposition 6's cost form. It does show that Proposition 6
  cannot be proved by a change of basis alone.
- **"The single-term argument extends to three arrays" (Proposition H3).** Refuted by the prove track itself. The
  exact unit HOLO3 (three arrays, 3m moves) keeps one term per input throughout and still has holonomy. Theorem C
  still covers it.
- **"For common-frame designs, height systems exist exactly when every start frame is opposite every end frame."**
  False as the prove track's summary stated it, as its skeptic showed: the plain schedule with twists is common-frame
  and has a crossing, yet its levels form a height system. It holds for heights built from one function through the
  coordinates, which is how the proof file states it.
- **"Corollary T₀ needs at least two twists, combined."** Overstated, as the same skeptic found. What is argued is
  that at least two different twists occur and that some combination of them crosses, and that combination may be a
  single twist.
- **The packing round's "closest design"** fails the position test for every choice of coefficients (section 9).
- **"None of the orbit designs delivered even one output exactly"** says more than was decided (section 9).
- **"The negative controls show the same mix of verdicts as the targets"** is false. Dead for every choice of scalars:
  60 % and 74 % at the controls, against 86 % at the targets.
- **"A multi-term 2 in 7 in the block-4 cube would be an engine unit."** Refuted by the refute track's skeptic. The
  cube's symmetries do not lift to the engine's letters: no engine element acting as an apartment swap squares to
  the identity (all 16 lifts checked for each of three kinds of swap). Such a unit would be a frame-game program,
  not directly an engine unit. With the completion rule of section 6 it should then lift to an engine unit that uses
  scratch. That last step combines two results, and no skeptic has checked the combination.
- **"Counting packing by image frames is the sharpest count."** Refuted by the capacity skeptic with an exact block-2
  unit.
- **"The third shear packs only the crossing"**, on the community's U8 and T74 lines (section 4).
- **The capacity skeptic's own correction** that 2 in 7 at block 6 is an open two-array target. I am correcting it
  here: the two-layer bound already excludes it (section 5).
- Two figures of the orbit track have no output file behind them. It gave the cycle-orbit rate at block 4 as 20 %,
  but the skeptic's reruns gave 8.5 to 16 %. And of the seven designs it searched for special scalars, one is
  recorded only in its run log.

### 11. Warnings for anyone who searches

1. **`ledger.py` certifies a word only when the number of digit columns is 1 more than a multiple of 4.** It treats
   scalar matrices made of letters (for example −I = XZXZ) as plain numbers, but in the real model they act as
   (−1)^f. A three-line word that `ledger.py` accepts is not a unit: it fails at every even f (`t_scalar.py`; the
   skeptic confirmed it). Every search whose last check is `ledger.py` must add a check in all four column-count
   blocks, or on raw arrays at f = 1 to 6 or 8, or certify in the frame game instead (`lift.py`).
2. **The "one term per input" filter holds only for plain mixes without scratch.** A two-array word in `ledger.py`'s
   format that reads one array twice is really a scratch word. The order-3 control written that way keeps one term
   per input throughout and still has holonomy.
3. **Run `heightfilter.py` on every candidate word before solving for coefficients.** If it finds a height system,
   the word cannot save, at any block and with any number of arrays.
4. **On designs with injected orbit copies, the linearised position test is blind to rigidity.** It treats each copy
   as if it rode in an array of its own, so its passes there mean nothing. Its refusals stand.
5. **Screens for any candidate, cheapest first.** Height filter; packing counted by live frames, at least the target
   saving and irreducible in every choice of coordinates; for two arrays without scratch, one input on two labels in
   one array; the position test at the drawn coefficients; the top-state count (at least W stored arrays must reach
   the target level); the layered screens (the part of the readout from level d or higher must already be right);
   then the exact solve in the frame game, and `lift.py`.

## What is still open, with a first experiment for each

1. **Proposition 6.** No unit costs less than the cheapest common-frame design with the same number of arrays. For
   two arrays it says exactly "no two-array unit saves", so a proof would close every open two-array target
   (2 in 7 at block 4; 2 in 7 and 2 in 8 at block 5; 2 in 8 to 2 in 10 at block 6; the pair targets). A
   counterexample would be a second saving mechanism. For three or more arrays even the cheapest common-frame cost is
   unknown. Together with the network-coding round's bounds for common-frame designs with three arrays, a proof would
   probably also close 3 in 10 and 3 in 11 at block 4 (see [../other-fields/](../other-fields/)). *First
   experiment:* a two-height version of Theorem C. Combine two height functions by their maximum, which is the
   network-coding round's lead. Test it against the calibration units in `programs/heights/checks.py` before trying
   to prove it.
2. **Conjecture T** (the single-term branch). Two inputs, scratch and single-term moves need 2m moves. It holds at
   blocks 2 and 3, where my certificates (proved in Lean there) already exclude every two-array saver. For free
   twists and for a single twist it is argued by hand (not checked; one skeptic). It is open from block 4. Two candidate invariants held on
   every state tested. *Conjecture S* says that two frames made available by two different objects have levels
   adding up to at most the number of moves. *Conjecture Φ* says that a dimension count over levels never exceeds the
   moves, and it would prove T at every block. Both come down to one question: when a release splits one lineage
   between two objects, has their shared prefix already been paid for? *First experiment:* write a candidate
   invariant and test it with `programs/single-term-game/t_inv.py`, which checks it on every reachable state of the
   small games in seconds. In parallel, write a compiled version of `stg.py` for the block-4 Hamming game at 6 to 7
   moves and the block-5 cube (2 in 8, 0.139), since Python stops near 10^8 states.
3. **The smallest open configuration.** Two inputs, scratch, single-term moves, and two kinds of twist, one trivial
   (both terms at the same group element) and one crossing. Every theorem here leaves it open. *First experiment:*
   search single-term 2 in 7 at block 4 restricted to these two twist classes, using composed pairs and releases
   through single terms, in the full block-4 frame game rather than the cube.
4. **Conjecture TC.** Every top-closed two-array program needs 2m moves, cancellations inside arrays included. If it
   holds, every two-array saver must cancel junk *below* the target frame, against arrays that have cost fewer moves,
   and the order-3 trick never does that. *First experiment:* build closed patterns from exact relations among paid
   moves rather than free symmetries. Moves along orthogonal directions commute, and by the braid relation three
   moves can equal one. Aim them at junk below the target frame, and run each design in the frame game through the
   position test at the drawn scalars, the top-state count and the layered screens (`programs/orbit-designs/orb.py`
   has all three) before the exact solve.
5. **The refund question for two arrays.** Is irreducible live packing always paid back? A proof would settle every
   open two-array target negatively, and it would amount to a coding-gain bound for networks whose arrays move in
   rigid bundles. *First experiment:* compute P* and the refund for every exact two-array unit at blocks 2 to 4 that
   `ledger.py` and `fpa.py` can enumerate, including the order-3 controls, and look for a quantity that pays P* back
   in every case.
6. **Must every saver be common-frame? Are my Lean savers common-frame in some coordinates?** The skeptic expects not,
   which would make them savers outside the class in name only. *First experiment:* replay one call of the E8 word
   (the certificate is in the repository's `tools/certificate/`) with `programs/packing-measurement/fm_ledger.py`,
   which already tracks every register's frame classes. Then test whether one choice of reference per input puts
   every register on a single class at every step. That test needs new code: `ledger.py`'s `find_coherent_gauge` does
   it only at blocks 1 to 4.
7. **Can re-framing raise the deficit, not only cut W?** My record words already re-frame their helpers, which cuts
   the number of arrays but not the deficit. The general bound of section 3 allows more. *First experiment:* take
   the smallest re-framed exact unit `ledger.py` can hold at block 3 or 4, and measure whether re-framing changes P −
   B − 2D − Z.
8. **A completion without scratch.** The completion rule always uses one scratch array, so the no-scratch theorems
   (Corollary 2′, Theorem A, Proposition D) do not transfer to the frame game. *First experiment:* at block 2,
   enumerate small frame-game programs with `programs/frame-game/fpa.py` and ask, for each, whether a scratch-free
   completion exists, by solving for the free coefficients directly.
9. **A multi-term 2 in 7 in the block-4 cube, with scratch.** This is open. A find would be a frame-game program.
   `lift.py` would then write it out as an engine word, to be checked in all four column-count blocks. *First experiment:* run the module form of `fpa.py` restricted to the
   cube's frames, with the screens of warning 5 first.
10. **Common-frame coding gain with three or more arrays.** Can a common-frame unit with three or more arrays save at
    blocks 3 to 8? Cuts cannot decide it, and my certificates cap the gain on the block-3 graph at 15/14. *First
    experiment:* a small exact search on the block-3 frame graph (135 frames, 3 to 6 sessions), starting from
    `programs/common-frame-accounting/frames.py`.
11. **The full-sharing family U_m** (4 in 15, 8 in 31, 16 in 63). Nothing proved excludes it. It would need source
    labels whose dimension grows with m. *First experiment:* at block 4, test whether any choice of plane labels
    (dimension 2) for 4 in 15 passes the median bound and the cut condition, using `programs/ledger/cutconst.py`'s
    cut machinery.
12. **"2 in m + 2 is impossible"** would close 2 in 7 at block 5 and 2 in 8 at block 6. The refute track's first
    session (numerical frame-game searches; no skeptic has checked that session) found that in two small games the
    top two levels alone admit a 2 in 5, and only the third level rules it out. *First experiment:* write the
    third-layer argument for those two small games (D_7, and S_4 with a free subgroup of order 3) and check it against
    the exact single-term solver.
13. **Measurement.** Replay one whole unit exactly, to replace the whole-unit estimates of section 4. The E8 unit is
    the smallest (1,263 arrays at block 45, 5 calls). Count *useful* packing, meaning moves that carry two or more
    labels that are each later delivered, starting with OpenAI's unit, whose contents are already exact (`oa_ledger.py`).

A later idea round (16 October, in [../idea-rounds/](../idea-rounds/)) noted that a streamed, unrolled unit escapes
Proposition 6, Conjecture T and Conjecture TC formally, because its number of arrays grows without bound. Its
smallest cells are two-array programs with split targets, though, so the theory here returns in that form. That
remark is argued by hand, not checked, and no skeptic has read it.

## Programs, data and proofs

Everything runs with Python 3.9 or later. numpy is needed by most programs and scipy by a few (`cutlp.py`, `vflow.py`,
`t_phi.py`, `t_bilin.py`, `t_psibil.py`, `t_hit.py`, `t_hit2.py`, `vhu.py`). Each folder in [programs/](programs/) has its own
README with the commands. Some programs write their JSON output next to themselves or into `out/` in the current
folder, so run them from a copy if you want to keep the shipped outputs. [data/](data/) holds the outputs of the
original runs, with the same file names. [proofs/](proofs/) holds the working proofs, each with a header giving its
status now.

| folder | what it shows | findings |
|---|---|---|
| [programs/ledger/](programs/ledger/) | exact engine at blocks 1 to 4 with the full ledger, coordinate tests and the exhaustive identity checks; cut constants | sections 1, 2, 3; warning 1 |
| [programs/common-frame-accounting/](programs/common-frame-accounting/) | frame geometry from actual matrices; the 0.1 and 0.01 needs | section 3 |
| [programs/heights/](programs/heights/) | `heightfilter.py` (Theorem C as an exact filter), the named test units (SPLIT1, HOLO3, the order-3 control) in the true model, the `ledger.py` scalar warning | sections 7, 8, 10; warnings 1 to 3 |
| [programs/packing-measurement/](programs/packing-measurement/) | packing ledger of one call of each local word and of OpenAI's network | section 4 |
| [programs/packing-capacity/](programs/packing-capacity/) | Theorem A on actual units, the rigidity tables | section 5 |
| [programs/frame-game/](programs/frame-game/) | the frame game, its exact solver, the position test, the packing-first designs | sections 6, 9 |
| [programs/frame-game-completion/](programs/frame-game-completion/) | the completion rule and `lift.py` (frame-game program to engine word), with float and modular evaluators | section 6 |
| [programs/orbit-designs/](programs/orbit-designs/) | orbit designs, the screens, the parity-law calibration | section 9 |
| [programs/single-term-game/](programs/single-term-game/) | the exact single-term game solver and the invariant tests | section 9; open items 2 and 3 |
| [programs/skeptic-checks/](programs/skeptic-checks/) | the skeptics' independent programs (own engines, own solvers) | every section |

**Re-run for this revision.** On 19 October I re-ran 40 programs (43 runs) from a clean copy of these folders on a
Mac. Every one ran without error. Where an output file is shipped, the new output equals it line for line, apart from
timing lines, the paths printed in warnings, one summary line that the shipped output of `t_strict.py` lacks, and
one line that the shipped output of the skeptic's `vhit.py` repeats.
Among them are the block-2 exhaustive ledger check (692,402 runs, 0 violations, 6 minutes), the E8 call of the
packing measurement and the block-2 exhaustive check of the completion rule (10,800 products). The list of what was
re-run and what it was compared with is in [programs/README.md](programs/README.md).
The large searches were not re-run for this revision: the 107-million-state cube run, the 111-million-state Hamming
run, the design waves of the packing rounds and the block-3 exhaustive ledger check. Their outputs are in `data/` as
they came out.

**Not included.** Three of the words measured in section 4 are my conversions of community words into the checker's
format: the U8 circuit (U8X), Jeremy Kothe's X4 word (X4X) and a T74-line word (T74g). They are not copied here. Their
SHA-256 hashes are listed in the Sources so that anyone holding them can check them. The E8 and block-110
certificates that the same programs read are already in this repository under `tools/certificate/`. Raw agent logs,
machine job folders and the round briefs are not published.

## Credits

- OpenAI's unit and its source labels are from OpenAI's "Exact Fourier transforms below n log n"
  ([openai/math](https://github.com/openai/math)). Its network was measured here from the preprint's rows and frame
  table.
- The community words measured in section 4, from the CrocSwap/integer-mult-bounds leaderboard:
  - the U8 circuit: pull request #869 by Jeremy Kothe, on his #799 word, which carries Chafik Boukhalfa's #741, and on
    icekylinx's #645, with the lineage those pull requests credit;
  - Kothe's #799 X4 word;
  - the T74 words: Rohan Arun (#1218, #1297, #1319 and the T74 exact-fill owner construction, with OpenAI Codex
    assistance), Marius-Paul Dumitrean (#1625), and Jeremy Kothe and LeverRank (#1672), on Kothe's k8 words with Thomas
    DiFiore, Jeremy Diviney and Chafik Boukhalfa. The T74-line word measured here is the variant that the
    community's registry names after the gauge44 construction, which the community's notices credit to icekylinx
    (#645).

  The full lineages are in this repository's README and NOTICE.
- Published theorems used: T. C. Hu, "Multi-commodity network flows", Operations Research 11 (1963), the
  two-commodity flow theorem behind Theorem 4.3; the multicommodity flow–cut gap of N. Linial, E. London and
  Y. Rabinovich (Combinatorica 15, 1995) and Y. Aumann and Y. Rabani (SIAM Journal on Computing 27, 1998), behind the
  W·m/(C·log 2W) bound; the Papernov–Lomonosov classification of demand graphs, as surveyed in A. Schrijver,
  *Combinatorial Optimization*; the simplicity of Sp(2m, 2) for m ≥ 3; and Clifford–Mackey theory, for one remark of
  the completion track. The agents cited these from memory and did not re-read them.
- The anatomy round's identity, the certificates of Theorem R, the two-layer bound and the order-3 control words are
  my earlier work, described in the sections linked above.

## Sources

Working names and dates. All four rounds ran on 15 October 2026 (Eastern time). The names below are the working
folders of the rounds and tracks. The programs' comments and the proof files refer to each other by them.

| folder here | round and track (working name) |
|---|---|
| `programs/ledger/`, `proofs/ledger-for-every-unit.md` | slack round, combinatorics track (`slack/sl-combinatorics`) |
| `programs/common-frame-accounting/`, `proofs/common-frame-accounting.md` | slack round, algebra track (`slack/sl-algebra`) |
| `programs/skeptic-checks/slack/`, `proofs/slack-skeptic.md` | slack round skeptic (`slack/v-slack`) |
| `programs/packing-measurement/` | packing round 1, measurement track (`packing/fp-measure`); `labels_check.py` from the anatomy round (`anatomy/an-openai`) |
| `programs/skeptic-checks/packing-measurement/` | its skeptic (`packing/v-fp-measure`) |
| `programs/packing-capacity/`, `proofs/packing-capacity-and-rigidity.md` | packing round 1, capacity track (`packing/fp-capacity`) |
| `programs/skeptic-checks/packing-capacity/` | its skeptic (`packing/v-fp-capacity`) |
| `programs/frame-game/`, `proofs/frame-game-sandwich.md` | packing round 1, build track (`packing/fp-build`) |
| `programs/skeptic-checks/frame-game/`, `proofs/frame-game-sandwich-skeptic.md` | its skeptic (`packing/v-fp-build`) |
| `programs/frame-game-completion/`, `proofs/frame-game-completion.md` | packing round 2, completion track (`packing/pk2-complete`) |
| `programs/skeptic-checks/frame-game-completion/` | its skeptic (`packing/v-pk2-complete`) |
| `programs/orbit-designs/`, `proofs/closed-patterns-and-conjecture-tc.md` | packing round 2, orbit track (`packing/pk2-orbits`) |
| `programs/skeptic-checks/closed-patterns/` | its skeptic (`packing/v-pk2-orbits`) |
| `programs/heights/`, `proofs/heights-theorem-c.md` | Proposition 6 round, prove track (`prop6/p6-prove`) |
| `programs/skeptic-checks/heights/` | its skeptic (`prop6/v-p6-prove`) |
| `programs/single-term-game/`, `proofs/single-term-game.md` | Proposition 6 round, refute track (`prop6/p6-refute`); `data/single-term-game/runN/` are its machine runs (job folders 10 to 19) |
| `programs/skeptic-checks/single-term-game/`, `proofs/single-term-game-skeptic.md` | its skeptic (`prop6/v-p6-refute`) |

Entries of my private results catalogue that this section covers: M23 (every saving is net frame packing), M24 (the
refund identity), M25 (capacity and rigidity), M26 (parity), M27 (the frame game is the whole problem), L26
(Theorem C), L27 (what a two-array saver must do), L28 (common-frame accounting), X43 (the single-term game), X44
(packing-first and orbit designs), T17 (screens), T18 (the `ledger.py` warning), N12 (Proposition 6 and Conjecture
T), N13 (Conjecture TC), C21, C22, C25 and C26 (refuted claims). The two-layer bound used in sections 5 and 10 is
entry L4 (upheld as Theorem T′ by 14 October). The network-coding track (`packing/nc-gain`) and its skeptic,
which ran in the same round, are covered in [../other-fields/](../other-fields/).

Words measured in section 4 and not shipped, by SHA-256 of the gzipped file in the checker's format:
- U8X (the U8 circuit with exit frames): `62fbff5c7e789c896deaa4b65d22d0ebb9374dd0f2326cdf9c332c34dba495d0`
- X4X (Kothe's #799 X4 word): `f046713d7ab279d38d4d0aa083b931ccd23d8a197fcbc18a7e557de7a847c78d`
- T74g (the T74-line word, gauge44-566): `1981cb2f93ff98aa9ce76e1ff16a7ba5c0c941bf489df25c99e0d8322ef91233`
  (its lowering record is `data/packing-measurement/lower_x-T74g.json`)
- Read from this repository: `tools/certificate/gcert1-e8-r783.json.gz`
  (`4b92f00fc7454b9b15e6d71a3b05062792eacca0a0ca4efee46f87827d57c8aa`) and `tools/certificate/gcert1-p11-pr193.json.gz`
  (`15251deb0c11e0504d7b297c9786818888a93f08b92b03c6c122784b70488633`).

Machines: the exact computations ran on an Apple M4 Mac mini, and the exhaustive searches of section 9 and some
skeptic runs on Linux machines of Google's Colab service (8 processors, about 50 GB). The re-runs of 19 October ran on
the Mac. A private provenance list maps every file here to its working file.
