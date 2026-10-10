# What did not work

Author: Jacob Sussman. Text of 2026-10-10, written with the fourth revision of this repository.

The subject of this repository: the Walsh-Hadamard and Fourier transforms of length n can be computed in
O(n (log n)^z) operations with an exponent z below 1, and the proof assistant Lean checks the proofs. The
saving is 1 - z, how far the exponent is below 1. Larger is better.

These are the dead ends of the work behind this revision. The work ended with a new small unit (the fixed
program that the method repeats at every scale) and a saving of 8.7625e-4; the third revision had 7.4745e-4.
Every dead end is told in the same four lines: what was tried, how large the attempt was, what came out, and
why it failed or was beside the point.

I publish them so that nobody, human or AI system, has to try these twice, and so that the ones that failed
for a weak reason can be tried again properly. The path that worked is in
[how-the-e8-unit-was-found.md](how-the-e8-unit-was-found.md) ("the first note" below). What is untried is in
[open-directions-2.md](open-directions-2.md) ("the next note").

**The most instructive failure.** A search for a much smaller unit found none. Most of the machine time
of that search, 11.6 hours by the clock, went to arithmetic modulo 2. A dozen-line argument shows that nothing
can be found there. (The argument is on paper: it is not in Lean and has had no human review.) Two lessons.
Before a search spends hours, spend an hour asking whether a short argument already settles it. And before
searching a game, check line by line that it is the game the checker plays. Section 9 tells this in full.

**The dead ends at a glance.** Sections 1 to 6 are about the new unit. Sections 7 and 8 are about other
families of labels (every pair of data arrays in a unit carries a label, a string of bits). Section 9 is
about tiny units.

| section | what was tried | what came out |
|---|---|---|
| 1 | five other designs for the helpers of the new unit | all stayed below the third revision's saving |
| 2 | designs with a large symmetry group | large groups gave the worst designs |
| 3 | cheaper scratch copies (the copies that hold the totals of the scatter) | none found; whether they exist is open |
| 4 | fewer stops per helper | not reached |
| 5 | four more devices (tricks that save helpers or moves) | all legal, none gains anything on this design |
| 6 | nine refinements of the search | little or nothing |
| 7 | "fat" labels, which would save more per array | no family that saves |
| 8 | other families of labels, smaller and larger | no unit with a larger saving |
| 9 | tiny units | nothing found, and most of the time went where nothing could be found |

## Words used here

- **unit**: the small fixed program that the whole method repeats at every scale; its size and its saving
  decide the final exponent.
- **data arrays / pairs**: the working copies of the input that the unit must transform; they come in pairs.
- **helper array**: an extra working copy that the unit borrows, uses and gives back unchanged; helpers cost
  time and save nothing themselves, so fewer is better.
- **saving** (written like 8.7625e-4): how far the exponent of log n is below 1; larger is better.
- **certificate**: the data file that lists every step of a unit, which the checking programs and Lean read.
- **price**: the saving that the Python checking programs compute for a certificate; it becomes a theorem
  only when Lean has built it.
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

## How to read the results

I directed this work and chose its directions. The AI agents (Claude) that I ran wrote the programs, ran the
searches and wrote the reports behind every result here. Those reports are one or two days old and have not
had a line-by-line human review.

- **A saving is a price unless it says *Lean*.** The whole number in brackets is the checker's exact number:
  8.3202e-4 (8320202) means the exponent 1 - 8320202/10^10.
- ***On paper*** means argued in writing by the agents with computer checks, not in Lean and without a human
  review.
- **Rounds.** The work on the new unit went in three rounds, told in the first note. "The winning design" is
  the design that won round 1. Round 2 ended at 8.3202e-4. Round 3 ended at the 8.7625e-4 of this revision.
- **"The third revision"** is the state of this repository before this revision. Its unit is the larger one
  (1,320 pairs, label width 22), and its saving is 7.4745e-4 (7474547), a Lean theorem of this repository.
  The "published" checkers, format and generator are the ones in `tools/` of this repository.
- **Outside work** is credited by name where it is used. The pull requests are listed under Sources at the
  end.

## 1. Designs for the E8 unit that stayed below the saving of the third revision

- **Tried.** Five ways to arrange the helpers of the E8 unit, other than the design that won.
- **Size.** Five designs, with one priced result each.
- **Result.** All five stayed below the 7.4745e-4 of the third revision. The best reached 6.8052e-4.
- **Why.** Each for its own reason, given under the table. Three agents built on one exact fact about
  triangles (the first note, section 3), and all three hit the same wall: the level below the top.

| design | helpers per pair | saving (price) |
|---|---|---|
| A. covering without shared sums | 19.72 | 3.9240e-4 (3923956) |
| B. the architecture of the third revision's circuit transplanted | 13.92 | 5.0732e-4 (5073193) |
| C. triangles on top, shared sums below | 12.29 | 5.5352e-4 (5535177) |
| D. triangle top with the machinery of the winning design | 9.66 | 6.8052e-4 (6805162) |
| E. top-down: seven 8-source sums per target, chosen target by target | 28.18 | 2.7028e-4 (2702762) |

**A. Covering without shared sums.** In designs where each helper meets its readers at one frame, at least
1,848 helpers are needed (15.4 per pair). That prices at about 4.87e-4. The best file is 28 % above that
floor. The reason: only coverings invariant under a group were searched, and no group tried admits a covering
without overlap.

**B. The architecture of the third revision's circuit transplanted.** The generator of that circuit cannot be
run at label width 9, the width of E8. Its accounting was transplanted into a small compiler. Every merge of
two partial sums at a frame new to both costs one array. The circuit of the third revision puts 3.6 additions
into each helper, this one 2.0.

**C. Triangles on top, shared sums below.** The 120 pairs split into 40 triangles, and a helper that holds a
sum of 27 sources serves two targets, which is cheap. Everything below those sums was a greedy search, one
that picks a sum at a time. That is where it stopped.

**D. Triangle top with the machinery of the winning design.** It took 9.3 helpers per 27-source sum, where
matching the saving of the third revision needs about 7. The cause is less sharing. Two sources meet in 24 or 28 of the
120 lists of what a target needs, but in only 4, 6 or 8 of the 120 lists of the triangle sums: about a
quarter of the sharing. Putting the triangle top on only a few triangles of the winning design made it worse:
+19 helpers for one triangle, +74 for four.

**E. Top-down.** There were 661 different sums for 840 uses, so almost nothing was shared. Sharing at the top
needs a global plan, not a choice target by target.

The triangle fact itself is true and exact (the first note, section 3). Three agents built on it and all three
hit the same wall, the level below the top.

## 2. Large symmetry groups

- **Tried.** Symmetric designs: every sum is taken together with all its images under a group of reorderings
  of the labels. The intuition was that a family this symmetric wants a symmetric design. The group of label
  maps that carry a design to a design has 348,364,800 elements.
- **Size.** One exact result for all permutations of the 9 points, exact optima for two large groups, greedy
  searches with groups of order 8 to above 3,000, a comparison of five more groups, and the integer programme
  with and without symmetry.
- **Result.** The opposite of the intuition. Large groups gave the worst designs.
- **Why.** As far as the agents' reports explain it: a large group forces images of a sum that overlap in the
  same targets, and the later ones are worth little. The full list of reasons ends this section.

**Exact results.**

- For all permutations of the 9 points, an invariant design built from shared sums does not exist at all,
  apart from the wasteful one, which gives every unwanted term its own helper (an exact result).
- Exact optima of the top level for two large groups:

  | order of the group | saving |
  |---|---|
  | 1,512 | 2.3773e-4 (2377262) |
  | 504 | 3.6989e-4 (3698864) |

  These are optima of one narrow class only: one level of sums of 6, 8, 10 or 16 sources, each taken with all
  its images under the group. The programme minimises a count of deliveries and sums, and the saving is the
  price of the certificate built from its answer. So they are no ceiling for designs with that symmetry. The greedy
  search below builds sums in several levels and came out a little higher with a group of order 504. The
  reports name the group only by its order in both places.

**The greedy search.** A greedy search picks one sum at a time. Groups whose order is a power of 2 between 8
and 128 gave the best designs, and groups of order above 500 the worst:

| group | order | helpers | saving |
|---|---|---|---|
| a group of order 504 | 504 | 2,356 | 3.8068e-4 (3806765) |
| the group of order 128 below, with one reflection added | 1,024 | 1,109 | 7.0813e-4 (7081283) |
| the group of order 128 below, with another reflection | above 3,000 | 2,160 | 4.1178e-4 (4117771) |
| no group at all | | 1,619 | 5.2369e-4 (5236901) |

How to read this table:

- The winning design uses a group of 64 reorderings of the 9 points. A group of 128 reorderings of 8 of the
  points gives exactly the same design. The two middle rows add one more map to that group of 128, which the
  agent's report calls a reflection. With it the group is much larger: order 1,024, as the report gives it.
  The new map is not a reordering of the 9 points: no group of such reorderings has order 1,024.
- The rows come from different versions of the greedy rule, so only the middle two compare directly with the
  winning design (1,028 helpers, 7.5590e-4). The order-504 run is an early one without in-place pairs. The
  run with no group has in-place pairs but not the last step, the removal of target stops that save no
  helper.

**Round 3, with all device passes applied.** A device is a trick that saves helpers or moves; a pass applies
it to a finished design. Five other groups of order 16 to 128 were compared again by a fast estimate. One
(order 128) gives the same design as the group of 64, the group of reorderings that the winning design uses.
The other four stay 11 to 40 % behind. The devices do not close the gap.

**The integer programme.** In round 3 an integer programme chose all the sums at once, from a fixed pool of
candidate sums. There, every bit of symmetry imposed costs (table in the first note, section 3). The optimum
that is exactly invariant under the group of 64 is worse than the greedy design: through the chain it gives
8.1515e-4 (8151461) against 8.3202e-4. The reason: the 8 targets fixed by many group elements must take their
sums by whole orbits, and they end with 10 to 16 deliveries each. (A delivery is one sum handed to one
target.)

**Why**, as far as the agents' reports explain it:

- Applying a sum with all its images under a large group forces images that overlap in the same targets, and
  the later ones are worth little.
- A few targets fixed by many group elements pay for exact invariance.
- A small group whose order is a power of 2 made the greedy search coherent. Once a pool of good sums
  existed, dropping the symmetry was better.
- The family also has two kinds of labels (triples and pairs), and a pair target needs choices that the full
  group of permutations cannot make.

## 3. Cheaper scratch copies for E8

- **Tried.** To bring the cost of the scratch copies below 72. A scratch copy is the array that holds one
  total of the scatter, and the cost is counted in copy moves.
- **Size.** Proofs on paper with exact checks, a numerical solver on 108 sets of frames, and one alternative
  built in full.
- **Result.** No set of totals costing fewer than 72 copy moves was found. A floor of 55 was proved on paper.
  The range 55 to 71 is open.
- **Why.** Some of the cheaper variants are excluded by proof. For the others nothing was found, which is
  not a proof.

The copy moves are the strongest lever on paper. The 72 of them take back 360 of the 480 moves the unit could
save. One copy move is worth about 4 % of the saving.

What is known now:

- **A floor** (*on paper*): at least 55 copy moves for every circuit the format allows. So the moves saved by
  the unit, called D, are at most 205. The earlier floor was 8.
- **72 is forced for 8 totals** (*on paper*, exact check). With exactly 8 totals and the scatter in use, all
  8 must stand at the full frame. The known matrix is locally the only fitting matrix of rank 8 (floating
  point).
- **Totals at frames of one dimension d**, which the published Lean generator needs. Values below 72 are
  possible only in the six cases of this table. The last column is for decompositions whose scatter delivers
  exactly the known matrix.

  | dimension d | totals | copy moves | status |
  |---|---|---|---|
  | 8 | 8 | 64 | excluded by proof |
  | 7 | 8 | 56 | excluded by proof |
  | 7 | 9 | 63 | exactly tight; no set of frames meets the necessary conditions. Not found, not a proof |
  | 7 | 10 | 70 | 108 sets of frames meet the necessary conditions; a numerical solver failed on all 108. Not found, not a proof |
  | 6 | 10 | 60 | excluded by proof |
  | 6 | 11 | 66 | excluded by proof |

- **An alternative of equal cost exists and was built**: 9 totals one dimension lower, also 72 copy moves.
  Both published checking programs accept it. It needs more helpers, 7,305. Its saving is 1.5295e-4
  (1529544), below the wasteful reference.

So the lever is not proved dead (55 to 71 is open) and nothing cheap opens it. For the size of the third
revision's unit, the note [scratch-copies.md](scratch-copies.md) stands as it was.

## 4. Stops

- **Tried.** To let each helper climb in fewer blocks, that is, with fewer stops.
- **Size.** Two rounds of work: a re-scheduler with an exact bound, then four attempts to make the stops part
  of the design.
- **Result.** Re-scheduling gave +0.74 %. The larger aim was not reached: blocks per helper went from 4.200
  to 4.320.
- **Why.** In this family of designs the stops follow from the hierarchy.

A helper makes about four blocks. A what-if is arithmetic with no design behind it. A what-if with three
blocks per helper priced the round-2 design at 1.214 times its saving: 1.0099e-3 against 8.3202e-4. Two
rounds tried to collect it.

**Re-scheduling the finished file** (the same additions, with frames and order chosen again).

- It gained +0.74 %.
- An integer programme over the schedule gave a bound. No schedule of those 3,244 additions has fewer than
  4,011 blocks (the file had 4,470), and the saving cannot rise more than about 1.5 %.
- Two further ideas brought nothing: a search guided by the relaxed optimum, and longer annealing.
- The re-scheduler is kept as the last pass of the chain.

**Round 3, stops as part of the design: the aim was not reached.** Blocks per helper were 4.200 at the
start of the round and 4.320 in that agent's best file. Four things were tried:

1. **Direct reads chosen together with the schedule.** In a direct read a gate reads a source array itself,
   not a helper copy of it. Of 229 candidates, 0 were repaired. A source array walks one chain of frames, and
   in every case the new read stands at a frame of the same dimension as an existing stop of that array.
2. **A better choice of copies in the compiler.** The greedy rule gives 2,108 chains of nested frames, and
   the minimum is also 2,108. Nothing to gain.
3. **Target stops chosen at compile time instead of by the pruning.** The pruning is the step that removes
   target stops that save no helper. There are 232 possible low reads. Both ends are worse than the pruning,
   and the range between them was not searched:

   | low reads kept | saving |
   |---|---|
   | none | 8.2324e-4 (8232351) |
   | all 232 | 8.2096e-4 (8209645) |
   | the pruning's 20 | 8.3202e-4 |

4. **A cost per stop in the integer programme.** Nine single runs gave between -0.20 % and +1.00 %, seven of
   them positive. At the settings closest to those of the final certificate it gave +0.02 %. Up to a moderate
   cost it only picks another optimum among ties. Four runs with identical counts span 0.34 %, so nine single
   runs prove little. Blocks per helper did not fall in any of them.

**Why**: in this family of designs the stops follow from the hierarchy.

- A pair (u + v, u - v) whose two sums go to two targets must be formed at a frame inside both targets'
  frames.
- A copy that two readers need at frames that do not contain each other must be born inside both.
- A helper handed to a target one dimension below the full frame has one move left.

The what-if with three blocks per helper is not a schedule or a compiler option of this family.

## 5. Devices that turned out to be idle on this design

- **Tried.** Four more devices. Each is a legal trick, as far as it was tested. Small test certificates that
  use the first and the fourth were accepted by both checkers. For the second, a file of the test harness was
  accepted when the harness was made, and it was not run again. The third is a way to make more direct reads
  possible, and the agent's report names no test certificate for it.
- **Size.** Counts on the design: 2,638 sets of deliveries for the first two devices, 297 candidates for the
  third. The fourth was tested on small files only and was not built into the E8 design.
- **Result.** All four gain nothing here.
- **Why.** The design has no place where the first three apply. The fourth is cost-neutral at best.

| device | places where it applies | why |
|---|---|---|
| targets read by targets | none | Of 2,638 sets of deliveries to single targets, no two targets receive the same sum, apart from 12 two-source sums that are already shared. This device has also never been through Lean. |
| one helper read by two targets at the meet of their frames | 0 pairs | The same reason. |
| moving a gate up to a larger frame so that a source array can reach it | 0 of 297 | Every such move puts a foreign label into a helper that some target reads later, and that target does not know the foreign source. |
| sources kept low until after the scatter | not built into the design | The totals of E8 stand at the full frame, so whatever feeds them cannot serve a target afterwards. Cost-neutral at best. |

The limit behind the first device pass, the one that reads sources directly, is the same everywhere: a
source array has one chain of frames. Of 303 further candidates for a direct read, 297 need the source at a
frame that is not nested with the one it already uses.

## 6. Search routes that gave little or nothing

- **Tried.** Nine ways to push the search further: in the greedy rule, in the integer programme, and in the
  passes that run after them.
- **Size.** From single runs to a full table of 301 changes. It is given with each route.
- **Result.** Little or nothing.
- **Why.** Given with each route.

**1. Random orders of the greedy choices.** They fail badly. With four random tie-breaks the winning group
gave 1,884 to 2,096 helpers instead of 1,028.

**2. Repairing the worst-served targets by hand.** A few targets read 7 to 18 helpers where the others read 4
to 6. The repair gained one helper. Every saved read became a helper that finishes early.

**3. The scan of single changed greedy choices, pushed to its end.**

- Size: a full table of 301 single changes from the round-2 design, and real scores for the top of it. Of the
  301, 122 repeat an earlier hierarchy. Of the other 179, 128 have exactly the cheap score of the base.
- Result: one more changed choice, +0.38 %. That is 8.3518e-4 (8351847) at 848 helpers.
- Pairs and a triple of the best single changes fell back onto a single result or below it. A change at one
  round of the greedy rule changes the lists of tied candidates at every later round, so the gains do not
  add.
- Scored by a better estimate, 26 changed choices in the first 13 rounds found nothing. 8 of the 26 turned
  out not to be real changes; 17 lost 6 to 32 %.
- The cheap score finds the large steps and is blind below 0.3 %.

**4. A device term inside the greedy rule.** +0.15 % on the plain greedy and -0.04 % on top of the changed
choices.

**5. A larger pool for the integer programme.** The pool is the set of candidate sums that the programme
chooses from. On three occasions a larger pool gave a better objective and a lower saving:

- 2,958 sums (the sums of 206 more greedy runs): objective 891.6 against 899.2, saving 8.6863e-4 (8686299)
  against 8.7168e-4 (8716750);
- 1,538 sums: objective 895.6, saving lower by about 1 %;
- with a stop cost, the larger pool was again below the smaller.

The objective is a stand-in for the saving. Beyond about 1 % it does not rank designs.

Pools made by one closure step (every pair of pool sums that the greedy rule would accept) have 31,200 to
243,904 new sums. They were not solved, for size: the programme used about 15 KB of memory per variable, so
about 80,000 variables fitted the memory allowed to one job.

**6. The plain helper count as the objective.** The worst weight tried (the first note, section 3).

**7. Re-use of finished helpers pushed further.** A finished helper is one that nothing reads any more;
re-use gives it a new value. Inside the compiler, two attempts brought 43 to 44 re-uses and no more. A new
pass on the finished file added less with each better file:

| file | gain from the new pass |
|---|---|
| the round-2 design | +1.17 % (18 more re-uses) |
| the first two files of the integer programme | +1.02 % and +0.40 % |
| its third file | +0.05 to +0.11 % (one or two re-uses) |
| the certificate of this revision | nothing |

Of 713,180 ordered pairs of a finished helper and a later one in the round-2 file, 712,556 fail at once on
frames. Each better hierarchy left the pass less to do. One variant gains 0.06 to 0.12 % by letting contents
be quarters and twelfths. It is not in the certificate, whose units are those of the files already checked
in Lean.

**8. Deleting helpers that are built and read once.** Of 827 helpers, no candidate.

**9. A second re-use pass after re-scheduling, or another re-scheduling seed before it.** Nothing.

**Two mistakes worth knowing about**, because they cost time and will recur:

- **A missing rule in the integer programme.** Its first two runs reported 922 where the design rebuilt from
  the solution had 1,048 to 1,054 helpers. The programme built the second sum of an in-place pair for
  nothing. (In an in-place pair, a + b and a - b are formed on the two arrays that held a and b.) The missing
  rule: every built sum is used at least once.
- **Fast estimates validated less than planned.** The estimate was compared against the whole chain on 10
  designs, not 50. Inside a family of optima of one programme the rank correlation was 0.60. The estimates
  were right to drop every candidate more than 0.5 % behind and are useless for the last 0.5 %.

## 7. Fat labels

- **Tried.** Families of fat labels. A fat label has dimension a of 2 or 3, so a data array would save two or
  three moves instead of one. The labels used everywhere else are thin.
- **Size.** 209 fat families: a = 2 at widths 4 to 8, a = 3 at widths 6 to 8.
- **Result.** No fat family is known to save anything at all, thin families in disguise aside.
- **Why.** For the smallest cases there are proofs on paper. For the rest the negative is weak: the exact
  step of the search hit its memory cap in 161 of the 209 families.

This was the search that flagged E8 by accident. For its own purpose it found nothing.

**The search.** None of the 209 families came below ratio 1, where a gain needs below 0.8. (The ratio is the
copy ratio: copy moves divided by pairs.) As a negative this
is weak: the exact step hit its memory cap in 161 of them, and with twice the memory on five it still did not
finish.

**Files that pass and gain nothing.** Eight files with a = 2 or 3 do pass the verifier at ratio 0.600 or
0.778. Each is the thin E8 or E7 written two or three times side by side. (E7 is the smaller family of
section 8.) The width grows as fast as the label, so nothing is gained.

**What is proved** (*on paper*, in the cost model of the search program, which prices shared sums only):

- without a division by 2 the least loss is exactly the plain one, for all families;
- families of two simple kinds (spreads and coordinate families) have exactly the plain loss;
- every family with width 4 and a = 2 loses at least 0.75 of the plain loss;
- every family with (width, a) = (3, 2), (4, 3) or (5, 4) loses at least 0.93 of it and so cannot gain.

**The two smallest undecided cases** were then settled exactly, each by a proof (*on paper*):

- A family of 6 labels at width 4 with a = 2 has least loss 12 of 12. The proof is a table of numbers checked
  in exact arithmetic on all 3,969 shapes of a shared sum, found independently by two agents.
- All 28 labels at width 3 with a = 2 have 56 of 56. This proof needed one extra step beyond the table and
  has one argument, not two independent ones.

**How far this goes.** No fat family is known to save anything at all, thin families in disguise aside. For
larger cells this is a weak negative. The two settled cases have almost no pairs that know each other (two
pairs know each other when their labels overlap in an even number of places), and the thin families that
save do it through many such pairs.

**Not checked at all:**

- whether the Lean engine, the checker that a unit must really pass, admits fat labels;
- the helper arrays that would build the shared sums, which are the main cost of the third revision's unit and are
  not in the model.

## 8. Other thin families: widths 3 to 8, and 11, 13 and 15 points

- **Tried.** Other families of thin labels: the small widths 3 to 8, fully symmetric families on 7 to 13
  points, and the larger relatives of E8 on 11, 13 and 15 points.
- **Size.** Proofs for widths 3 and 4, searches for widths 5 to 8, one search over the symmetric families,
  and runs of the round-2 method on 220, 364 and 560 pairs.
- **Result.** No unit with a larger saving than E8.
- **Why.** The small families gain nothing or too little. On the larger ones the method of round 2 does not
  keep its efficiency.

**Small widths.** The copy ratio is the number of copy moves divided by the number of pairs; a five-stage
unit gains when it is below 0.8.

| width | what is known |
|---|---|
| 3 and 4 | excluded by proof in the published format (*on paper*; width 4 in the checker's cost model, for three and five stages) |
| 5 and 6 | searched; best copy ratio 1.000 |
| 7 (E7) | a unit exists (63 pairs, block size 35, 7 moves saved), and both checkers accept a wasteful circuit for it, at about 3.8e-5: too small to matter |
| 8 | best copy ratio 0.893; no gain with three or five stages |

**Fully symmetric families on 7 to 13 points.** Nothing better than E7, E8 and "triples plus pairs" was found
in a search that is complete inside a stated class, with some very large orbits skipped.

**Triples plus pairs on 11, 13 and 15 points.** E8 is the family of this kind on 9 points. On paper the
larger families promise more per bit of label width than E8 (+2 %, +16 % and +19 %). I had the method of
round 2 carried to them as a test of whether it scales. It does not keep its efficiency:

| points | pairs | what was run | helpers per pair | saving (price) |
|---|---|---|---|---|
| 9 (E8), round 2 | 120 | the search of round 2 | 7.11 | 8.3202e-4 (8320202) |
| 11 | 220 | the same search, ported unchanged | 13.4 | 5.3608e-4 (5360784) |
| 11 | 220 | a rule read off the best E8 design for the top two levels, the old search below | 9.47 | 6.9548e-4 (6954826) |
| 13 | 364 | the same rule, two steps of the search cut short | 10.36 | 7.4774e-4 (7477404) |
| 15 | 560 | the same rule on a shortened run | 12.57 | 6.9348e-4 (6934825) |

How far these files were checked:

- The 11-point certificate (2,084 helpers) passes both checkers, the ten readiness tests for a Lean build and
  a stand-alone replay. Eight copies with one change each are all refused.
- No Lean was run on any of these files.
- The 13-point saving is level with the third revision's 7.4745e-4 by a margin of 0.04 %. It is not a claim.

**What carried over:** the format, the checkers, the compiler and the four device passes, and the stop
pattern (4.1 stops per helper against 4.2).

**What did not:** everything below the top two levels.

- At 9 points two thirds of the sums come in pairs that cost no helper. At 11 only half do, and the sums that
  cost a helper each go from 404 to 1,238.
- The rule read off the E8 design gives nothing at all for 27 % of the targets.
- Each step of two points costs about 3 to 5 times the time.

These runs used the method of round 2. The integer programme of round 3 has not been tried on them.

## 9. Tiny units, and eleven hours on arithmetic modulo 2

- **Tried.** A search for a tiny unit: a handful of arrays, each with only 2 to 4 steps to travel, served
  with fewer moves than the plain schedule. It would be a new mechanism. One hit would be worth far more than
  everything above.
- **Size.** 31.2 core-hours in arithmetic modulo 2, and 23 seconds in arithmetic modulo 3 and modulo 5.
- **Result.** Nothing found.
- **Why.** Most of the machine time went to a place where nothing could be found. And the game that the
  program searched was narrower than the game the real checker plays.

This is the most instructive failure of the two days. Here it is in order.

**The game that was searched.** W arrays each have to travel m steps on a graph of frames. The number m is
the block size. The rules:

- one step of one array is a move;
- combining arrays that stand at the same frame is free;
- extra arrays that start empty are free.

The plain schedule costs W m moves. A scheme with fewer would be a hit. The program searched schemes with
coefficients modulo 2, modulo 3 and modulo 5. It took modulo 2 first.

**What it cost.** The whole of an overnight run went to arithmetic modulo 2: 31.2 core-hours, 11.6 hours by
the clock. It covered 90 classes of starting positions and 29.8 billion states. Every one came back "no
scheme".

**The short argument that was missing.** In arithmetic modulo 2 no scheme uses fewer than W m moves, for
every m and W. The argument counts conditions:

1. Each array carries the condition "what I hold lies in the subspace of the frame where I sit".
2. At the start all conditions are free.
3. Free mixing adds none.
4. One move adds at most one, because neighbouring frames share all but one dimension.
5. The end needs W m of them, because a start frame and its target share nothing.

The argument needs 1 + 1 = 0. It is a dozen lines. It is *on paper*: six readings by agents found no gap, and
two independent programs replayed it on every state of the smallest cases. It has had no human review and it
is not in Lean.

**The consequence.** A saving scheme must divide by 2 somewhere. So a search modulo 2 could never have shown
one. The unit of the third revision says the same thing from the other side: all of its gain comes through
halves.

**Credit.** eumemic made these statements first, in a pull request to the outside repository
CrocSwap/integer-mult-bounds ("the hub") that was opened on the day the search was stopped (Sources, 1). It
states:

- that nothing saves at block size 2 for any frames, with a certificate script;
- as a reported row of a table, without a proof, that "gates without a (1+i) denominator never save, at any
  m";
- the difference between the two games described below.

The proofs for the game my program searched were made in this project on the day, three independent ones for
block size 2. They were made here because the model of eumemic's script is not that game, and the lemma that
script needs is not in its package.

**What was left of the search.** The real negatives are small: only 3 arrays at block size 3 on the subspace
frames, modulo 3 and modulo 5. That is 14 classes, 15.9 million states, 23 seconds. All 16 unfinished layers
of the queue are now decided without a search. The search is stopped.

**What else is excluded for that game** (*on paper*):

- block size 2: every number of arrays, every field;
- block size 3 or 4, on any frames: 3 arrays or fewer. This is a counting argument over a complete list of
  small graphs, 177,815 cases at block size 3 and 14,653,045 at block size 4;
- block size 3 on the subspace frames: 7 arrays or fewer;
- block size 3 on the engine's frames, with 4 arrays: 10 moves or fewer.

The smallest open case there is 4 arrays in exactly 11 moves, on one of 161 shapes.

**And the game was too narrow.** The search allowed additions only between arrays at the same frame. The
Lean engine, the checker that a unit must really pass, has no such rule. Its free addition takes any rational
combination of the arrays, with no condition on frames, and the proof of the third revision uses that freedom. So a miss
of the search says nothing about what the engine accepts.

For the engine's own rules (*on paper*):

- 2 arrays cannot be served in 5 moves at block size 3 (they need 6);
- a unit needs at least 105 W / 41 = 2.561 W moves at block size 3, so at most 14.6 % is saved;
- a unit in which no addition has an even denominator saves nothing, at every block size.

3 arrays in 8 moves is open (the next note).

**The lesson I take from it.**

- Before a search spends hours in the cheapest field, spend an hour asking whether a short argument empties
  that field.
- Before searching a game, check line by line that it is the game the checker plays.

## 10. What these negatives rest on

- Every negative here is a result of one or two days of work by AI agents. Most were checked by a second
  agent. None has had a human review.
- "Proved" means argued in writing and tested by a program, not checked in Lean.
- A "searched, not found" is only as strong as the search. Sections 6, 7 and 8 say where a memory or time
  limit cut a search short.
- None of these negatives changes any theorem of this repository.

## Sources

1. eumemic, pull request 288 to the outside repository CrocSwap/integer-mult-bounds ("the hub"):
   https://github.com/CrocSwap/integer-mult-bounds/pull/288. It was opened on the day my search of section 9
   was stopped. Section 9 says what it states and what was proved here.
2. The first note, [how-the-e8-unit-was-found.md](how-the-e8-unit-was-found.md): its section 3 has the
   triangle fact, the table of the cost of symmetry and the weights of the integer programme.
3. [scratch-copies.md](scratch-copies.md): the scratch copies at the size of the third revision's unit.
