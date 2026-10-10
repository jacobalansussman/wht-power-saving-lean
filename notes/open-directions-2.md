# Open directions, second note

Author: Jacob Sussman. Text of 2026-10-10, written with the fourth revision of this repository.

The subject of this repository: the Walsh-Hadamard and Fourier transforms of length n can be computed in
O(n (log n)^z) operations with an exponent z below 1, and the proof assistant Lean checks the proofs. The
saving is 1 - z, how far the exponent is below 1. Larger is better. This revision, the fourth, has a new
small unit (the fixed program that the method repeats at every scale) with a saving of 8.7625e-4.

This note lists what has not been tried and may be promising, after the work behind this revision. There are
eight directions. I ordered them by promise, the most promising first, by my own judgement. Each direction
has three parts: what is known, a concrete first experiment, and an honest guess of the odds.

The main conclusion: the two directions with fair odds are small ones. Tuning the program that chose the
design may give about 1 % more saving, and more only if new sums pay. Designing for fewer stops (the halts
in an array's climb, explained in the list below) may give 1 to 2 %. The directions that could give much
more have odds from below even to very low, or unknown odds. Every estimate of a gain or of the odds here is
a guess and is marked so. None is a promise.

Three other notes belong with this one:

- [how-the-e8-unit-was-found.md](how-the-e8-unit-was-found.md) tells the path that worked.
- [what-did-not-work.md](what-did-not-work.md) lists the dead ends.
- [open-directions.md](open-directions.md) is the first note of open directions. It stays as it was written.
  Where the two overlap, section 9 says how.

I directed this work and chose its directions. The AI agents (Claude) that I ran wrote the programs, ran the
searches and wrote the reports I drew on. Those reports have not had a line-by-line human review.

## Words used here

- **unit**: the small fixed program that the whole method repeats at every scale; its size and its saving
  decide the final exponent.
- **data arrays / pairs**: the working copies of the input that the unit must transform; they come in pairs.
- **helper array**: an extra working copy that the unit borrows, uses and gives back unchanged; helpers cost
  time and save nothing themselves, so fewer is better.
- **saving** (written like 8.7625e-4): how far the exponent of log n is below 1; larger is better.
- **certificate**: the data file that lists every step of a unit, which the checking programs and Lean read.
- **price**: the saving that the Python checking programs compute for a certificate; it becomes a theorem only
  when Lean has built it.
- **label, source, target**: a pair is a source array and a target array, and it has a label, a string of
  bits. The number of bits is the width. Two pairs *know* each other when their labels overlap in an even
  number of places. E8 is the family of 120 labels of width 9 that the new unit uses. The block size is the
  number of directions along which every array must be transformed: 45 for E8, five times the width.
- **shared sum, hierarchy, chain**: helpers carry sums of sources to the targets that need them. A shared sum
  is built once and used more than once. The hierarchy is the list of shared sums and which target takes
  which. The chain is the sequence of programs (passes) that turns a hierarchy into a certificate. The pool
  is the list of candidate sums that the hierarchy is chosen from.
- **frame, move, stop, block**: a frame is the subspace of labels where an array stands; it only grows. A
  move is a climb of one dimension. A block is one climb between two stops. The full frame is the largest.
  Moves are what cost time.
- **scatter, totals, copy moves**: at the scatter a few sums, the totals, are added into all targets. The
  copy moves are the summed dimensions of the totals' frames. They are a cost: with v pairs and cst copy
  moves a five-stage unit saves D = 4v - 5 cst moves.

Three status words say how firm a statement is:

- A saving with a whole number in brackets, such as 8.7625e-4 (`8762479`), is a *price* unless it says *Lean*.
  The whole number is what the checking program prints: the exponent is 1 - 8762479/10^10.
- *What-if* is arithmetic with no design behind it.
- *On paper* means argued in writing by the agents, with computer checks.

A gain of 1 % means 1 % of the saving.

## The starting point

The unit of this revision, in numbers:

- 120 pairs and 783 helpers, that is 6.53 helpers per pair;
- labels 9 bits wide, block size 45;
- 72 copy moves and 120 moves saved;
- a saving of 8.7625e-4 (`8762479`).

The programs that rebuild it are in `tools/e8/`. `run.py` runs the whole chain. `solve.py` and `xmilp.py` are
the integer programme, the program that chooses the hierarchy.

## The eight directions in one table

Most promising first, by my own judgement. The last two are checks and give no saving.

| # | direction | what a success would be worth | odds (a guess) |
|---|---|---|---|
| 1 | tune the integer programme: a better objective, tie-breaks, new sums | about 1 % from the untried items; more only if new sums pay | fair |
| 2 | make the stops part of the design | 1 to 2 % from the first experiment; up to about 23 % only with a new kind of design | fair for the small gain, unknown for the large |
| 3 | a larger family of labels, on 13 points | a unit above E8, if at least one helper in eight can be removed | below even |
| 4 | fewer copy moves, between 55 and 72 | about 4 % of the saving per copy move | low |
| 5 | three arrays in eight moves at block size 3 | a new mechanism, over a hundred times the saving | very low; a proof that it is impossible is the likelier result |
| 6 | "fat" labels in the larger cases | a first fat family that saves | low |
| 7 | run the comparator, an outside checking program, in a real sandbox | a check, not a saving | high that it passes |
| 8 | a line-by-line human review | trust, not a saving | - |

## 1. The integer programme: objective, tie-breaks, new sums

The program that chose the design of this revision is not used up. I rank this direction first, and its
experiments are all small.

**What is known.** An integer programme chooses the hierarchy all at once, out of a pool of 1,470 candidate
sums, and its answer is the exact optimum for that pool. Two weights in its objective stand for the passes
that run afterwards.

The search ran in three rounds, and the programme came in round 3. That round gained 5.3 % in about an hour.
In the order of the steps:

| step of round 3 | gain |
|---|---|
| the optimal hierarchy | +2.2 % |
| the numbering of the sums | +1.0 % |
| one change in the compiler, the pass that turns the hierarchy into a first certificate | +1.5 % |
| the two weights | +0.5 % |

Three facts say that the programme is not used up:

- **Equal optima differ.** Optima with the same objective value differ by 0.3 to 0.6 % in the saving. And the
  re-scheduler, the pass that chooses the stops again, can change which of two close candidates is ahead: in
  one measured case a design priced at 8.680052e-4 came out at 8.747650e-4, and one priced at 8.676999e-4
  came out at 8.755827e-4. All four figures are prices. The certificate of this
  revision is one optimum, taken as the solver returned it.
- **The objective is a stand-in for the saving.** On three occasions a larger pool gave a better objective
  and a lower saving ([what-did-not-work.md](what-did-not-work.md), section 6).
- **The solver's own hints are unused.** The solver closes the gap after one to three nodes. So the linear
  relaxation, the same problem with fractions allowed, is nearly enough, and its dual values say which
  missing sum would pay. They have not been used.

**First experiments.** There are eight, each small. Three words first. An *in-place pair* is the two sums
a + b and a - b, formed on the two arrays that held a and b. A *delivery* is one sum handed to one target. A
*device pass* is a pass of the chain that removes helpers after the hierarchy is chosen. The ones meant here
let the readers of a helper read the source arrays directly. The first (`xsub.py`) does it for a helper that
holds a single source, the second (`xexp.py`) for a helper that holds a small sum.

1. **Sample the ties.** Solve the programme thirty times at the weights 0.55 and 0.85, each time with the
   costs disturbed by at most 0.001. Send each answer through the whole chain. This takes about 45 minutes
   single-threaded. Guess: a few tenths of a per cent.
2. **Fit the objective to the chain.** Take 20 to 30 optimal hierarchies and fit one coefficient for each of:
   deliveries; sums formed alone; single-source parts; and the frame dimension of a sum, as a stand-in for
   stops. A chain run takes about 90 seconds, the programme about 5.
3. **Price new sums by the dual values**, instead of enlarging the pool blindly. One closure step (adding
   every pair of pool sums that the earlier greedy rule would accept) makes the pool too large to solve:
   31,200 to 243,904 sums. Pricing is the way to those pools.
4. **Put the halves into the programme.** Some sums are held in halves. Give each sum one bit: whole or
   halves. An in-place pair needs at most one part in halves, and it gives halves. Now the compiler decides
   this after the fact. Guess: 5 to 10 helpers. One helper is worth about 0.1 % of the saving.
5. **Choose the direct reads inside the programme.** Add a constraint that the frames at which one source is
   read form a chain, each inside the next. 232 of 272 candidates of the second device pass fail on exactly
   that. Guess: 10 to 20 helpers.
6. **Extend the second weight** to parts of 3 to 6 sources. The second device pass takes parts up to 6.
7. **Two new terms in the objective.** One is a term for the six targets that still take 10 deliveries. The
   other is a credit for a sum that can be born on a helper that has just finished inside the same frame. On
   a finished file 99.9 % of such pairs fail on frames, so it has to be decided here.
8. **A faster faithful estimate.** 281 of 323 candidates of the second device pass fail on frames, each after
   a full rebuild. A frame test before the rebuild would bring the estimate from 13 seconds to about 3. That
   would allow a local search on the programme's solution.

**Odds (a guess).** Fair.

- For more of the same (weights, tie-break samples, the re-use pass), the agent that checked the certificate
  guessed under 0.5 %.
- For the two things nobody has done, the same agent guessed about 1 %. The first is item 4, the halves in
  the programme. The second is to make new sums for the eight targets that the pool serves worst: they take
  7 to 10 sums where the others take 4. It is not in the list as such. Item 7 only adds a term to the
  objective for the worst of those targets, and item 3 is the general way to new sums.
- Pricing new sums (item 3) is the one with a chance of more.
- The returns were falling: 10 % in round 2, 5.3 % in round 3.
- Each item alone is about the size of the noise. So every comparison needs repeats: 10 to 20 optima per
  setting.

## 2. Stops as part of the design

The helpers of this unit make more than four blocks each. A what-if with three blocks each would be worth
about 23 % more saving. Choosing the stops again cannot collect it, because the stops follow from the design.
So the design itself has to aim at fewer stops.

**What is known.**

- Helpers make 4.34 blocks each in the certificate of this revision.
- *What-if*: with three blocks per helper, and one per source and target, the same 783 helpers price at
  1.0752e-3 (`10752492`). That is 1.227 times the saving.
- Re-scheduling, which chooses the stops again for the same additions, cannot collect it. Its room was
  bounded on two earlier designs: about 1.5 % on the round-1 design, and at most 3.3 % on a 917-helper design
  of round 2 (a first estimate; the solver stopped at its time limit). Neither bound is for the certificate
  of this revision.
- Round 3 did not lower the blocks per helper at all. In this family of designs the stops follow from the
  hierarchy ([what-did-not-work.md](what-did-not-work.md), section 4).

**First experiment.** Let targets read low by design. A target's own frame has dimension 8, one below the
full frame. To read low is to read a helper at a frame below that one.

- The limit: a target has one chain of frames, each inside the next. So it can read at most one helper per
  dimension.
- The idea: hand a target S only sums that are the twins of sums handed to one partner S'. The twin of
  u + v is u - v. Then all these sums would stand at the one frame of dimension 7 inside both targets'
  frames. Both targets could read them there, and 480 helpers would skip the stop at dimension 8.
- In the integer programme this is a reward for deliveries of one target whose sums share a frame. It needs
  a variable per (target, frame).
- The agent that proposed it guessed a gain of 1.7 %.

**A further step, and a large one.** Charge the cost per array instead of per sum. The true cost is a sum
over each array's chain of rank x ln(45 / rank), where the rank of a block is the number of moves in it. That
is concave, and a cost per stop misses it. Column generation over the lives of arrays would be exact.

**Odds (a guess).** Fair for 1 to 2 %. The rest of the 23 % would need a new family of designs: one in which
a helper is not handed over one dimension below the full frame and then moved once more. Nobody has proposed
one.

## 3. Thirteen points as a design question

E8 has a larger relative on 13 points that promises more on paper. On a shortened run with the method of
round 2 it gave less than E8, because it needed more helpers per pair. The integer programme of round 3 has
not been tried on it.

**What is known.**

- The labels of E8 are all triples and all complements of pairs on 9 points. The same family exists on 11,
  13 and 15 points.
- *On paper*, the 13-point family promises 16 % more than E8 per bit of label width.
- The method of round 2 delivered less of that promise as the family grew
  ([what-did-not-work.md](what-did-not-work.md), section 8):

  | | 9 points (E8) | 13 points, on a shortened run |
  |---|---|---|
  | helpers per pair | 7.11 | 10.36 |
  | saving | 8.3202e-4 | 7.4774e-4 (`7477404`) |

- To pass the round-2 E8 saving of 8.3202e-4, a 13-point file needs 3,290 helpers or fewer. That run had
  3,770.
- The bar set by this revision's 8.7625e-4 is higher and was not computed.
- Fifteen points lands lower with the present rule and costs hours a run.

**What is untried.** The integer programme of round 3 at 13 points. Also a rule for the sums below the top
two levels, and for the 27 % of targets that the present rule leaves without any. (The present rule was read
off the best E8 design; see [what-did-not-work.md](what-did-not-work.md), section 8.)

**First experiment.** Two steps:

1. Run the two slowest steps of the search to the end at 13 points. (Steps of the search program are meant,
   not the stages of the unit.) They were cut short by a 600-second limit per job. A full run is 20 to 25
   minutes and about 200 MB.
2. Build a pool from that design's sums and solve the programme with no symmetry imposed. At 9 points that
   took 22,466 variables and seconds. The size at 13 points is unknown.

**Odds (a guess).** Below even that it passes E8. The gap to close is at least one helper in eight: 3,770
against 3,290, and more against this revision's saving. And E8 itself moved 5.3 % in round 3.

## 4. Copy moves between 55 and 72

Every copy move takes back five of the moves that the unit saves. The unit has 72 copy moves, and the floor
on paper is 55. Whether a unit with fewer than 72 exists is open.

**What is known.** With 120 pairs the unit saves D = 480 - 5 x cst moves, where cst is the number of copy
moves. Here cst is 72, so D is 120.

- *On paper*: cst is at least 55 for every circuit, so D is at most 205.
- *On paper*: with 8 totals and the scatter in use, 72 is forced.
- One copy move is worth about 4 % of the saving.
- The one alternative that was built has 9 totals one dimension lower. It also has 72 copy moves, and it
  needed more helpers.

What is left open depends on where the totals stand:

- **All totals at one dimension** (the published Lean generator needs this). Two cases are open. The first
  is 7-dimensional frames with 10 totals (cst = 70): 108 sets of frames meet the necessary conditions, and a
  numerical solver failed on all of them. The second is decompositions whose scatter is far from the known
  matrix.
- **Totals at mixed dimensions.** The only floor is 55.

**First experiment.** An exact attempt, not a numerical one, on the 108 sets of frames. For each set, solve
for a fitting matrix in rational arithmetic, or prove that none exists.

**Odds (a guess).** Low. The agent that worked on it advises treating 72 as fixed. A search with nine or
more totals may only prove that 72 is the floor. And a saved copy move has to survive the helpers it costs.

## 5. Three arrays in eight moves at block size 3

This is about a possible new mechanism at tiny size, not about labels. The question: can three arrays be
transformed at block size 3 in eight moves, where the plain schedule needs nine? A yes would be worth over a
hundred times the present saving. The agents expect a no, and a proof that it is impossible is the likelier
result.

**What is known.** There are two sets of rules, two "games":

- **The wider game** is what the Lean engine really applies. A free addition may combine any arrays with
  rational coefficients, with no condition on frames.
- **The narrower game** allows additions only between arrays at the same frame.

A pull request by eumemic to an outside repository states this difference between the two games (Sources,
item 1).

In the wider game, the two smallest cases at block size 3 are, *on paper*:

- **2 arrays in 5 moves: impossible.** They need 6. The argument gives every state a score. Free steps
  cannot raise the score, and a move raises it by a bounded amount. The amounts are read off a table that
  attaches a small space of vectors to each frame. Five moves supply at most 5.86, and two arrays need 6. The
  same outside pull request lists this case as open in its own model.
- **3 arrays in 8 moves: open.** The plain schedule needs 9. The same method gives at least 315/41 = 7.68
  moves. Excluding 8 needs a better table, with an error count of 118 or less. The error count adds up, over
  the 945 pairs of neighbouring frames, how far the table is above the ideal value; a perfect table would
  have 0. The best table found has 162. A table at 118 would settle the case: with 8 excluded, 9 moves are
  needed. What the method can never do is show that nothing saves for every number of arrays, because no
  perfect table exists once dividing by 2 is allowed.

An 8-move scheme would have to do all of this:

- make at least two moves on an array that holds a mixture of inputs;
- divide by 2 in some addition;
- use an address map or a phase;
- leave the equal-frame form;
- waste almost nothing.

Nobody has a mechanism for this. No finite search decides the case as it stands: the steps generate
371,589,120 matrices. A hit would be worth a saving of 0.107 to 0.127 (*what-if*, from the engine's formula).

**First experiment.** Search for a better table. It takes minutes of machine time, and it closes the case by
proof if it reaches 118. Before that, a human should read the half-page "averaging lemma" that the 2-array
argument rests on. Its weakest part is the reading of five Lean definitions as text.

**Odds (a guess).** Very low that a scheme exists; the agents expect none. Fair that a better table closes
it.

**In the narrower game** the smallest open case is 4 arrays in exactly 11 moves on one of 161 shapes. That
search has not been run or sized. A miss there says nothing about the engine's rules.

The outside pull request (Sources, item 1) lists the case of 3 arrays in 8 moves as open in its own model as
well, and has the earlier bound of 50/21 moves per array. The argument for 2 arrays, the table and its
checker are in `notes/network-coding/`. The same case is restated for readers from network coding, with what
that pull request has first and what is not claimed, in
[for-network-coding-readers.md](for-network-coding-readers.md).

## 6. Fat labels: the larger cases

The labels of E8 are "thin". A "fat" label has a dimension a of 2 or 3, and would save two or three moves
per data array. No fat family is known to save anything. The two smallest undecided cases were settled at
"no saving". Larger cases are open, and no construction is known.

**What is known.**

- No fat family (label dimension a of 2 or 3) is known to save anything.
- The two smallest undecided cases were settled exactly at "no saving"
  ([what-did-not-work.md](what-did-not-work.md), section 7).
- Two tools came out of that. The first is tables of rational numbers that prove a floor for all
  decompositions at once. The second is one extra step. It closed the 28-label case, where the tables alone
  stopped at 55 of 56.
- The earlier floors came from a narrow kind of table and may be far below the truth:

  | width | a | floor, as a part of the plain loss |
  |---|---|---|
  | 4 | 2 | 0.75 |
  | 5 | 3 | 0.744 |
  | 6 | 3 | 0.487 |

Where a gaining fat family would have to be, *on paper*:

- at width 6 with a = 3, away from the spread families (a simple kind of family named in
  [what-did-not-work.md](what-did-not-work.md), section 7);
- at widths 5 to 8 with a = 2;
- at widths 7 and 8 with a = 3;
- or wider.

What it would have to do, *on paper*:

- divide by 2;
- have labels that form a pattern that exists over bits and has no copy over ordinary fractions;
- avoid large groups of pairs that all know each other. (Two pairs know each other when their labels overlap
  in an even number of places.)

No construction is known.

**First experiment.** Not a search. Apply the two tools to the classes that are open by a floor only:

- all families at width 4 with a = 2;
- (width, a) = (4, 3) with 120 labels;
- (5, 4) with 496 labels;
- (5, 3).

Each proposed table takes seconds of exact computation. The work is finding the table.

**Odds (a guess).** Low that a fat family saves. Two things are unchecked, and they come first if one is
ever found: whether the engine admits fat labels at all, and the cost of the helpers, which the model leaves
out.

## 7. A comparator run in a real sandbox

This is a check, not a way to a larger saving. The comparator is an outside program that checks a Lean proof
against the statement it claims to prove. It is meant to run inside a sandbox, and the runs behind this
repository used a stand-in for the sandbox.

**What is known.**

- The comparator runs its builds inside the Linux sandbox `landrun`, which is built on Landlock, a feature of
  the Linux kernel.
- The runs behind this repository used a stand-in. On macOS the stand-in was `tools/landrun-macos-shim.sh`.
- On the Linux machine rented for this revision the Linux kernel has no Landlock. So `landrun` enforces
  nothing there, and a stand-in was used again.

For the certificate of this revision the comparator ran on the Mac (stand-in
`tools/landrun-macos-shim.sh`) and on the Linux machine (a stand-in built on namespaces). Both
configurations ended with `Your solution is okay!` on both machines. No run used the real `landrun`.

**First experiment.** Three things are needed:

1. a Linux machine whose kernel has Landlock;
2. `landrun` installed as VERIFY.md describes;
3. the comparator on the configurations of this revision, with `nanoda` switched on. `nanoda` is a second
   kernel, that is, a second program that checks Lean proofs.

**Odds (a guess).** High that it passes. The sandbox limits what a build may touch. It does not change what
Lean's kernel checks. It is still a gap in the checks, and I welcome the run.

## 8. A line-by-line human review

This is also a check, not a way to a larger saving. The new proofs of this repository have not had a
line-by-line human review, and I would welcome one.

**What is known, and where to start.** For this revision I would read in this order:

1. **What the theorem says.** The Lean statement gives an exponent for a program. It does not say "E8". The
   checkers accept, with the same saving, a copy of a circuit carried to another family of 120 labels.
2. **Fields that Lean does not read.** On an earlier E8 design, one record was made wrong on purpose: a frame
   number in the list of totals. Lean accepted it and the Python checkers refused it. By a text search, nine
   more fields are not read. The theorem is unchanged by this. But the certificate file holds data that the
   theorem does not depend on.
3. **The Python mirrors against the Lean check functions**, branch by branch. A mirror is a Python copy of a
   check that Lean makes. This has been read only in part. Every price in these notes rests on it.
4. **The hand-made part of the Fourier chain** for the new exponent. It is two files, both the text of the
   third revision's files with names and numbers replaced:
   - `Work/FourierE8/Seam.lean`: 96 lines changed against `Work/Fourier/Seam.lean`, 6 comment lines added;
   - `Work/FourierE8/Axioms.lean`: 10 lines changed.

   The other 13 Lean files of that folder are written by `tools/fourier/mkchain_e8.py`.
5. **The arguments on paper.** Each is a dozen lines to a page, and none is in Lean:
   - the floors for copy moves (55) and for helpers;
   - the argument that nothing saves in arithmetic modulo 2;
   - "2 arrays in 5 moves is impossible";
   - the two fat-label proofs.

   A formal version of the modulo-2 argument would be a small, self-contained Lean project.

## 9. Smaller items, and the first note

Eight smaller items.

- **Contents in quarters and twelfths.** The re-use pass lets a finished helper take a new value. A variant
  of it gains 0.06 to 0.12 % if the contents of arrays may be quarters (sources, helpers) and twelfths
  (targets). It passes the Python readiness tests. No Lean build of such a file was made.
- **E8 in three stages instead of five.** The unit of this revision has five stages. With three it would
  have block size 27 and save 24 moves. This is *on paper* only: there is no three-stage theorem for this
  certificate format. By a rough estimate, five stages wins while helpers per pair stay above about 4.5. So
  at 6.53 this is not the next step.
- **Targets read by targets.** This is legal and idle on this design, and it has never been through Lean. A
  design that used it would have a lower helper floor: 37 instead of 116, *on paper*.
- **The scratch copies of the first note.** Lead A1 of [open-directions.md](open-directions.md) is the same
  lever as section 4 here, at the size of the unit published before this one. It stands as written. For E8
  the numbers are 72 copy moves and a floor of 55.
- **Units smaller than E8.** E7 (width 7, block size 35, 7 moves saved) saves too little to matter. Width 8
  has no gain with three or five stages. Nothing was found at widths 3 to 6.
- **One search for stops and re-use.** The re-scheduler weighs every array alone. Counting a finished
  helper and the later helper that could take its place as one array would let it move frames toward a
  re-use. A guess of the agent that proposed it: a few more re-uses per file. *Untried.*
- **Sums that feed a total and are also read after the scatter.** In the files of the integer programme a
  third of the candidate pairs for re-use fall there. A compensating read is allowed before the scatter
  but costs the units; a variant that keeps halves has not been built. *Untried.*
- **Which sums of a larger pool do harm.** A larger pool gave a better objective and a lower saving three
  times ([what-did-not-work.md](what-did-not-work.md), section 6). Solving with each added group of sums
  alone would show which ones: five runs of about a minute. *Untried.*

## Sources

1. Pull request 288 of the outside repository CrocSwap/integer-mult-bounds, by eumemic:
   <https://github.com/CrocSwap/integer-mult-bounds/pull/288>. It states the difference between the two
   games of section 5 and lists "2 arrays in 5 moves" and "3 arrays in 8 moves" as open in its own model.
   The credit is spelled out in [what-did-not-work.md](what-did-not-work.md), in the section on tiny units.
2. The reports of the AI agents (Claude) that I ran. The results and the agents' guesses in this note are
   drawn from them. They have not had a line-by-line human review.
