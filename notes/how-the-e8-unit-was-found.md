# How the E8 unit was found

Author: Jacob Sussman. Text of 2026-10-10, written with the fourth revision of this repository.

The subject of this repository: the Walsh-Hadamard and Fourier transforms of length n can be computed in
O(n (log n)^z) operations with an exponent z below 1, and the proof assistant Lean checks the proofs. The
saving is 1 - z, how far the exponent is below 1. Larger is better.

This note tells, in order, how the new unit of this revision was found, with the saving reached at every
step. A unit is the small fixed program that the whole method repeats at every scale. The new one is called
E8 here. Its saving is **8.7625e-4**; the third revision had 7.4745e-4.

The path in five steps:

1. **The family.** A unit works on pairs of data arrays, and each pair carries a label, a string of bits. A
   search for something else turned up a family of 120 labels. Per bit of label width it promised about as
   much as the family of the third revision. So everything depended on the helper arrays: the extra working
   copies a unit borrows, which cost time.
2. **A wasteful first design**, with 56 helpers per pair, showed that a unit with these labels exists:
   1.6555e-4.
3. **Round 1.** Four designs side by side. The one that builds partial sums once on helpers and shares them
   passed the saving of the third revision: 7.5590e-4, with 8.57 helpers per pair.
4. **Round 2.** Four devices on the same sums, then a few changed choices: 8.3202e-4, with 7.11 helpers per
   pair. (A device is a trick that saves helpers or stops.)
5. **Round 3.** All shared sums chosen at once by an integer programme: 8.7625e-4, with 6.53 helpers per
   pair.

The main lesson is from round 3. The greedy rule that had chosen the sums one at a time was 5.5 % above the
optimum of its own pool of sums (the sums the greedy designs had produced, with their images under a symmetry
group). And the gain came from dropping the symmetry, not from using it. Waste remains in six places
(section 5), the number of helpers still among them.

The savings in this list are prices: the Python checking programs computed them. The designs of rounds 1 and
2 were also built and checked in Lean during the work. What Lean has checked of the last one is said under
"What is new here".

I wrote this for people and for AI systems that want to continue the work. Two notes beside it say what
failed ([what-did-not-work.md](what-did-not-work.md)) and what is untried
([open-directions-2.md](open-directions-2.md)). [e8-results.json](e8-results.json) is the same index for
machines.

## Words used here

The words are used as in the README.

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

## What is new here, and whose the parts are

New in this revision, all found on one day, 2026-10-10:

- **A much smaller unit, with a larger saving.**

  | | unit of the third revision | new unit |
  |---|---|---|
  | pairs of data arrays | 1,320 | 120 |
  | label width, in bits | 22 | 9 |
  | block size (five times the label width) | 110 | 45 |
  | helper arrays | 9,412 | **783** |
  | saving | 7.4745e-4 (7474547/10^10) | **8.7625e-4 (8762479/10^10)** |

  The new unit has 1,263 arrays in all. Its saving is 1.1723 times that of the third revision.
- **The family of 120 labels**, called E8 in this repository, used as a unit of this method, and a helper
  circuit for it (the arrangement of its helper arrays) designed in this project. The circuit inside the
  certificate of the third revision is the community's. The new one uses three devices of outside authors, named below.
- **The way the circuit was chosen.** An integer programme picks the shared sums all at once. Two weights
  make the programme aware of the passes that run after it. The greedy rule it replaced was 5.5 % above the
  optimum of its own pool of sums.
- **The search programs themselves**, in `tools/e8/`. `run.py` rebuilds the certificate from the recorded
  hierarchy and, on request, re-runs the pool and the solver stages. The integer programme is `xmilp.py` with
  `solve.py`. The README there says what each program does and what is not included.

Status of the figure of this revision:
The theorem `wht_main_block_B2Ge8x` at the exponent 1 - 8762479/10^10 was built by Lean's kernel in this
repository, on a Mac and on a Linux machine. It depends on the axioms `propext`, `Classical.choice` and
`Quot.sound` only. `tools/Compare.lean` printed `RESULT: PASS` on both machines, and the comparator printed
`Your solution is okay!` on both, each time with a stand-in for its sandbox.
`transform_mainE8` and `convolution_mainE8` state the exponent 1 - 8762478/10^10, one unit less than the
Walsh-Hadamard figure. They come from the same builds, depend on the same three axioms, and got
`RESULT: PASS` and `Your solution is okay!` on both machines.

Not mine, and used in the certificate of this revision. "The hub" below is the community repository
CrocSwap/integer-mult-bounds. The pull requests are listed under [Sources](#sources) at the end.

- **jamesyc and eumemic: the re-use of finished helpers with a compensating read.** A helper that nothing
  reads any more takes a new value. Whoever would receive the old content first reads it with the opposite
  sign.
- **Chafik Boukhalfa (account chafreaky): the late pairing rule.** The rule says which finished helper a new
  value takes: one that finished late, near the frame where the new value starts. The certificate uses his
  rule as the rule of choice in the re-use pass. His program, his pairs and his certificate are not used. On
  the first E8 design that had the pass, a control with the opposite rule kept 40 % less of the pass's gain.
- **icekylinx: the in-place pair.** The sums a + b and a - b are formed on the two arrays that held a and b.
  Two pull requests of ikeboy carry it.
- **The hub community: the frame networks and the circuit of the third revision.** This kind of network comes
  from the hub. The circuit inside the certificate of the third revision is ikeboy's, eumemic's and icekylinx's, as
  `RELATED-WORK.md` says.
- **The solver** is HiGHS, called through SciPy's `milp`.

I directed this work and chose its directions. The AI agents (Claude) that I ran wrote the programs, ran the
searches and wrote the reports this note is drawn from. Those reports are one day old and have not had a
line-by-line human review. I would welcome one.

## How to read the savings

A saving is 1 - z, where z is the exponent. It is written like 8.7625e-4, with the whole number of the
checker in brackets: 8762479 means the exponent 1 - 8762479/10^10. Each saving carries one of these status
words.

- *Lean*: a theorem of this repository, checked by Lean's kernel.
- *Lean during the work*: built and checked in Lean on the day, under a working name. The Lean files of that
  design are not included in this repository.
- *price*: computed by the published Python tools (`tools/gx/refcheck.py` and `tools/gx/gxdry.py`) for a
  certificate file that both accept. Not in Lean.
- *what-if*: arithmetic with no design behind it.
- *on paper*: argued in writing by the agents, with computer checks. Not in Lean, and without a human review.

## 1. The family came from a search for something else

**What was tried.** The night before, one of my searches was looking for "fat" labels: labels of dimension 2
or 3, which would save two or three moves per data array.

**What came out.** It found no fat family that saves. Beside the fat families it ran thin ones, and among
those it flagged three by their copy ratio. The copy ratio is cst/v: the copy moves (cst) divided by the
number of pairs (v). A unit of the five-stage kind used here saves D = 4v - 5 cst moves, so it gains when the
ratio is below 0.8. The three flagged families have ratios below 1: 0.600, 0.778 and 0.781.

An agent recognised the first two as cases of one construction. They are called E8 and E7 here, after the
root systems. E8 is this family:

- Take 9 points. The labels are all triples and all complements of pairs: 84 labels with three ones and 36
  with seven, 120 labels of width 9.
- Each pair knows 56 others.
- The scatter needs 8 totals, all at the full frame. So cst = 8 x 9 = 72, and D = 480 - 360 = 120 moves are
  saved per unit.
- The group of label maps that carry a design to a design has 348,364,800 = 960 x 9! elements (an exhaustive
  count). That is the order of the Weyl group of E8 divided by its centre.

**Why the next step was needed.** Per bit of label width the family promises about what the family of the
third revision does: 0.0278 against 0.0265 in the agents' measure (1 - 1.25 x copy ratio)/h, where h is the label width. So
everything depended on the helpers.

## 2. A deliberately wasteful circuit showed that a unit exists

**What was tried.** After the scatter each target holds its own source plus an unwanted half, with either
sign, of each of the 56 sources it knows. That makes 6,720 unwanted terms for helpers to take away. The
wasteful way is one helper per term.

**What came out.** That file has 6,728 helpers (56 per pair) and 7,208 arrays. Both published checkers accept
it: **1.6555e-4 (1655510), price**, 0.22 times the saving of the third revision.

**What it taught.** It proved nothing about efficiency. It showed that the published format, checkers and
generators take a unit at block size 45 as they are. And it set the question: two rough block models put
parity with the saving of the third revision at 6 to 10 helpers per pair, against 56.

## 3. The path, step by step

### Round 1: four designs side by side, and the floors

**What was tried.** I had several agents try different designs at once, on one scoreboard. One harness priced
every file with the published tools.

| design | helpers | per pair | saving (price) | why it stopped |
|---|---|---|---|---|
| covering without shared sums | 2,366 | 19.72 | 3.9240e-4 (3923956) | designs of this kind need at least 1,848 helpers (15.4 per pair), which prices at about 4.87e-4: below parity whatever the search does |
| the architecture of the third revision's circuit transplanted | 1,670 | 13.92 | 5.0732e-4 (5073193) | the top level is cheap, building the sums under it is not (below) |
| triangles on top, shared sums below | 1,475 | 12.29 | 5.5352e-4 (5535177) | the same wall, a little further |
| shared sums chosen by a small symmetry group, with in-place pairs | 1,028 | 8.57 | **7.5590e-4 (7558959)** | the first design above the saving of the third revision (1.011 times) |

**What came out.** One design passed the saving of the third revision: the last row, at 7.5590e-4. The other three
stopped below it.

**A fact worth keeping from the transplant.** The 120 pairs split into 40 triangles. Take a corner s of a
triangle (s, a, b). The 56 sources that s knows are exactly these: 27 shared with a, 27 shared with b, and a
and b themselves (checked on all 1,120 triangles). Given 120 "edge sums" of 27 sources, two helpers per pair
deliver everything. A circuit with only that top level would price at about 2.0e-3 (*what-if*). What breaks
is the level below: building the edge sums cost 11.9 helpers per edge sum, where parity needs about 7.

**The design that won the round.** It shares partial sums in a tree. A greedy rule picks, again and again,
the pair of partial sums that the most targets want together. It takes that pair with all its images under a
group of 64 reorderings of the 9 points (the translations of the 8-point space F_2^3 with two more
involutions). Where two groups of targets want u + v and u - v, both are formed in place on the two helpers
that held u and v. This happened 333 times. Its steps, all prices:

| step | helpers | saving |
|---|---|---|
| greedy sharing, no symmetry | 1,949 | 4.5955e-4 (4595534) |
| the same by orbits of a group of order 9 | 1,598 | 5.2484e-4 (5248392) |
| with in-place pairs, a group of order 16 | 1,299 | 6.1907e-4 (6190737) |
| ties broken toward the larger sums, the group of order 64 | 1,028 | 7.4679e-4 (7467861) |
| target stops that save no helper removed | 1,028 | 7.5590e-4 (7558959) |

Two things about it mattered later:

- It depends strongly on the tie-break rule. The same group with four random tie-breaks gave 1,884 to 2,096
  helpers.
- It is "guided" by the group, not invariant under it.

**The floors.** A floor is a lower bound: a count that no circuit can go below. In the same round two agents
worked on floors. All of them are *on paper*, with exact computer checks.

- Copy moves: at least 55 for every circuit the format allows, and exactly 72 for 8 totals with the scatter
  in use.
- Helpers: at least 8 for every circuit. With this scatter, at least 116 (0.97 per pair) when no target reads
  a target, and 37 when targets may.

The floors forbid little. At 116 helpers the saving would be about 2.8e-3 (*what-if*).

**Checked.** The 7.5590e-4 design was built and checked in Lean under a working name during the day (*Lean
during the work*).

### Round 2: four devices on the same sums, then better choices

**What was tried.** The winning design wasted helpers in ways that passes on the finished hierarchy could
repair. A pass is a program that reworks a finished design. The four devices are: sources read directly,
totals built on helpers that have already served, re-use of finished helpers, and re-scheduling. The table
and the list below count the second under the first.

**What came out.** 8.3202e-4 with 853 helpers, step by step:

| step | helpers | per pair | saving (price) |
|---|---|---|---|
| the round-1 design | 1,028 | 8.57 | 7.5590e-4 (7558959) |
| stops re-scheduled only (same additions, frames and order chosen again) | 1,028 | 8.57 | 7.6148e-4 (7614776) |
| re-use of finished helpers only | 985 | 8.21 | 7.6594e-4 (7659352) |
| sources read directly only | 910 | 7.58 | 8.0282e-4 (8028179) |
| all devices together | 868 | 7.23 | 8.2159e-4 (8215871) |
| the same with 8 of the 89 greedy choices changed | 853 | 7.11 | **8.3202e-4 (8320202)** |

What each step was:

- **Sources read directly.** The design made 414 helper copies of single sources. If the helper sums are
  formed before the scatter, a gate (one step of the circuit) may read the source array itself, at a frame
  above its line. The source walks on to the full frame as it must anyway. A second pass deletes a finished
  helper that holds a small sum (up to 6 sources) and lets its readers read the sources. A third builds each
  of the 8 totals on a helper that has already served. Together: 118 helpers.
- **Re-use of finished helpers.** 450 of the 1,028 helpers finished low and climbed idle to the full frame.
  All three certificates at label width 22 that were examined re-use such helpers; the one of the third
  revision does it 2,310 times. The E8 design did not yet. The pass took 43 or 44 of them, with the late pairing rule as
  its rule of choice. The credits are under "What is new here".
- **Re-scheduling.** The additions stay the same; their frames and order are chosen again. Alone it is worth
  +0.74 %. An integer programme over the schedule showed that no schedule of those additions can be more than
  about 1.5 % above the starting saving. The stops belong to the design.
- **Changed choices.** A scan changed single choices of the greedy rule and kept a change when a cheap score
  improved. With everything else equal it gave +1.36 % and 14 helpers: the same chain with the same options
  gives 8.2087e-4 (8208712) at 867 helpers on the unchanged choices, and 8.3202e-4 at 853 on the changed ones.
  The 867-helper file was priced in the work only and is not on the scoreboard. The row "all devices
  together" above is a neighbouring file made with other options. Against it the step is +1.27 % and 15
  helpers.

**Checked.** The 8.3202e-4 design was built and checked in Lean under a working name on two machines during
the day, with its Fourier and convolution twin at 1 - 8320201/10^10 (*Lean during the work*; files not
included).

### Round 3: the hierarchy chosen by an integer programme

**What was tried.** The greedy rule guesses one sum at a time. In round 3 one agent wrote the choice of the
whole hierarchy as an integer programme and solved it exactly.

- **The pool.** The programme chooses from a fixed pool of sums: the sums of the four greedy designs saved so
  far, with their images under the group of 64. That is 1,266 different sums, 1,470 with the images.
- **The programme.** It chooses which sums are built, from which two smaller sums each is formed, and which
  target takes which sum. The conditions: each of the 6,720 needed (target, source) terms is covered exactly
  once; every built sum is formed exactly once and used at least once. The costs: an in-place pair costs
  nothing; a sum formed alone costs 1 (one of its two parts finishes there); every delivery of a sum to a
  target costs 1. With no symmetry imposed it has 22,466 variables. The solver, HiGHS, reports gap 0 (its
  sign that the optimum is proved) after one to three nodes, in 2 to 27 seconds.
- **What it showed.** In the plain count the optimum over the pool is 938 where the greedy design has 990:
  **the greedy rule was 5.5 % above the optimum of its own pool of sums.**
- **Symmetry costs.** The programme was first proposed over orbits of the group of 64. That class is worse
  than the greedy design, which is not exactly invariant. The table gives the optimum for each symmetry
  imposed, with the weight of the next item:

  | group imposed | order | variables | optimum |
  |---|---|---|---|
  | the group of 64 | 64 | 1,034 | 944.4 |
  | a subgroup | 16 | 2,205 | 930.0 |
  | the translations | 8 | 3,287 | 925.2 |
  | a four-group | 4 | 6,318 | 915.6 |
  | one involution | 2 | 12,390 | 906.4 |
  | none | 1 | 22,466 | **899.2** |
  | the greedy design, for comparison | | | 932.4 |

  The gain came from dropping the symmetry, not from using it. The symmetry had been needed only to make the
  greedy search find anything. The pool it left behind is small enough to optimise without it.
- **The first weight.** The plain count ignores the device passes that run afterwards. Its optimum has only
  38 sums that the "sources read directly" pass can use, where the greedy design had 96. As an objective it
  was the worst weight tried: 0.6 to 1.3 % below the weight 0.4 after the device passes. The remedy: a sum
  formed alone whose finishing part is a single source is nearly free, because the source array can play that
  part. Counting such a sum as 0.4 instead of 1 fitted three runs of the chain within 4 helpers.
- **The second weight.** Another agent scanned that weight and added a second. A sum formed alone whose
  smaller part has two sources counts 0.85, since the second device pass can replace that helper by direct
  reads of the two sources. The certificate of this revision has the weights 0.55 and 0.85.
- **Two changes that are not the programme.** Both concern the compiler, the first pass of the chain, which
  turns the hierarchy into a file that can be priced.
  1. The sums that are parts of an in-place pair are numbered first, so that the compiler forms more pairs in
     place (338 instead of 319 of 367 in one measured case).
  2. The compiler now forms an in-place pair also when one of its two helpers holds halves. Both sums then
     come out in halves. Two helpers in halves stay excluded, because one sum would come out in quarters,
     which the checker's mirror refuses. This widens when the in-place pair of icekylinx applies; it adds no
     new device.

**What came out.** Where the gain of round 3 came from, measured (prices):

| step | helpers | per pair | saving |
|---|---|---|---|
| start of the round | 853 | 7.11 | 8.3202e-4 (8320202) |
| the optimal hierarchy (first weight 0.4) through the unchanged chain | 827 | 6.89 | 8.4996e-4 (8499603) |
| parts of in-place pairs numbered first | 815 | 6.79 | 8.5856e-4 (8585620) |
| in-place pairs also on a helper in halves | 794 | 6.62 | 8.7168e-4 (8716750) |
| weights 0.55 and 0.85 | 783 | 6.53 | **8.7625e-4 (8762479)** |

The compiler change alone, on the greedy hierarchy, gives 8.5775e-4 (809 helpers). It is worth +3.1 %, and
the optimal hierarchy adds +1.6 % on top of it. The two weights add +0.52 % and 11 helpers to 8.7168e-4.

**How firm the last step is.**

- Optima of the programme with the same objective value differ by 0.3 to 0.6 % in the saving, depending on
  which one the solver returns.
- Four agents ended the round within 0.1 % of each other (8762479, 8758761, 8757398, 8757259, all prices, all
  from this programme and chain).
- What supports the weights: all 10 weighted settings beat the unweighted one in a fast estimate (by 0.28 to
  0.53 %), and 4 of 4 beat it through the whole chain.
- Which of the top files is best is inside the noise. This revision takes the highest.

**The chain of the certificate of this revision**, pass by pass:

| pass | helpers | saving (price) |
|---|---|---|
| integer programme, weights 0.55 and 0.85, parts of in-place pairs numbered first | | |
| compiler with in-place pairs, target stops pruned | 927 | 8.1259e-4 (8125886) |
| re-use of finished helpers: 52 copies written on finished helpers | 875 | 8.2617e-4 (8261680) |
| single sources read directly | 819 | 8.5522e-4 (8552172) |
| small sums replaced by direct reads of their sources (up to 6) | 783 | 8.6861e-4 (8686103) |
| re-scheduling | 783 | 8.7625e-4 (8762479) |

In the programme's own terms the final hierarchy has 550 deliveries, 390 sums formed alone (68 of them with a
single-source part) and 367 in-place pairs. The whole chain runs in about 90 seconds.

## 4. Every approach, with helpers per pair and saving

The rows run from the first circuit to the certificate of this revision. Savings are prices unless the last
column says otherwise. "Times third" is the saving divided by that of the third revision (7474547).

| approach | helpers | per pair | saving | times third | outcome |
|---|---|---|---|---|---|
| one helper per unwanted term | 6,728 | 56.07 | 1.6555e-4 (1655510) | 0.22 | reference |
| nine totals one dimension lower | 7,305 | 60.88 | 1.5295e-4 (1529544) | 0.20 | still 72 copy moves |
| covering without shared sums | 2,366 | 19.72 | 3.9240e-4 (3923956) | 0.52 | cannot reach parity |
| architecture of the third revision's circuit transplanted | 1,670 | 13.92 | 5.0732e-4 (5073193) | 0.68 | failed |
| triangles on top, shared sums below | 1,475 | 12.29 | 5.5352e-4 (5535177) | 0.74 | failed |
| triangle top with the round-1 machinery | 1,159 | 9.66 | 6.8052e-4 (6805162) | 0.91 | clear negative |
| shared sums by a group of 64, in-place pairs | 1,028 | 8.57 | 7.5590e-4 (7558959) | 1.011 | *Lean during the work* |
| the same, re-scheduled only | 1,028 | 8.57 | 7.6148e-4 (7614776) | 1.019 | ceiling about 7.672e-4 |
| the same with re-use of finished helpers only | 985 | 8.21 | 7.6594e-4 (7659352) | 1.025 | |
| the same with sources read directly only | 910 | 7.58 | 8.0282e-4 (8028179) | 1.074 | |
| all devices, the greedy choices | 868 | 7.23 | 8.2159e-4 (8215871) | 1.099 | |
| all devices, 8 greedy choices changed | 853 | 7.11 | 8.3202e-4 (8320202) | 1.113 | best of round 2; *Lean during the work* |
| a ninth changed choice | 848 | 7.07 | 8.3518e-4 (8351847) | 1.117 | the greedy scan's end |
| integer programme, unchanged chain | 827 | 6.89 | 8.4996e-4 (8499603) | 1.137 | |
| parts of in-place pairs numbered first | 815 | 6.79 | 8.5856e-4 (8585620) | 1.149 | |
| in-place pairs also on a helper in halves | 794 | 6.62 | 8.7168e-4 (8716750) | 1.166 | |
| the same with a re-use pass on the finished file, units 1/4, 1/4, 1/12 | 792 | 6.60 | 8.7264e-4 (8726429) | 1.167 | units the Lean-checked files do not have |
| a stop cost in the programme | 785 | 6.54 | 8.7573e-4 (8757259) | 1.172 | |
| first weight 0.55, second re-scheduling pass | 782 | 6.52 | 8.7574e-4 (8757398) | 1.172 | |
| first weight 0.65 | 785 | 6.54 | 8.7588e-4 (8758761) | 1.172 | |
| **weights 0.55 and 0.85** | **783** | **6.53** | **8.7625e-4 (8762479)** | **1.172** | **the certificate of this revision** (status under "What is new here") |

For comparison, two outside certificates at label width 22, priced by the same Python tools:

| certificate | saving (price) |
|---|---|
| Chafik Boukhalfa's, from his pull request to this repository | 7.5474e-4 (7547361) |
| DreamingOfClouds's, from the hub | 7.6359e-4 (7635875) |

Both are above the 7.4745e-4 of the third revision. The two figures in this table are prices computed here
in Python. Chafik Boukhalfa's pull request has his figure as a Lean theorem. He reports the Lean kernel check,
and on 10 October 2026 the agents built his pull request on a Linux machine: all 54 Lean modules that it adds
compiled, and `tools/Compare.lean` printed `RESULT: PASS`. No Lean check of the other certificate was made
here. The pull requests are under [Sources](#sources).

## 5. Where the remaining waste is

The certificate of this revision still wastes in six places.

1. **Stops.** In the certificate of this revision the helpers make 3,396 blocks, 4.34 per helper. The round-2
   design had 4.20. Every gain of round 3 came from fewer helpers, none from fewer stops per helper.
   - A *what-if*: with three blocks per helper and one per source and target, the 783 helpers price at
     1.0752e-3 (10752492). That is 1.227 times the saving of the certificate. For the round-2 design the same
     what-if was 1.214 times.
   - Re-scheduling cannot reach it, because in this family of designs the stops are fixed by the hierarchy.
     Its room was bounded on two earlier designs: about 1.5 % on the round-1 design, and at most 3.3 % on a
     917-helper design of round 2 (a first estimate; the solver stopped at its time limit). Neither bound is
     for the certificate of this revision. Of the 853 helpers of the round-2 file, 480 are handed to a target one dimension
     below the full frame and then make one last move. And 1,989 of its 3,583 helper blocks have rank 1, the
     dearest per move.
2. **Copy moves.** There are 72. They take back 360 of the 480 moves the exchanges save. The floor is 55 (*on
   paper*), so D could be at most 205 instead of 120. Whether 55 to 71 can be had is open. One copy move is
   worth about 4 % of the saving.
3. **Helpers.** There are 783, that is 6.53 per pair, against floors of 0.97 or 0.31 per pair (*on paper*).
   The floors do not show that anything near them is reachable.
4. **Finished helpers that climb idle.** On a finished file the re-use is nearly used up. Frames stop 99.9 %
   of all candidate pairs, and on the certificate of this revision a re-use pass on the finished file found
   no candidate at all. What remains is decided when the sums are chosen.
5. **A few targets served badly.** In the round-2 design eight targets read 7 to 18 helpers where the others
   read 4 to 6. In the programme's hierarchies the worst are six targets with 10 and two or three with 7.
6. **The choice among equal optima.** 0.3 to 0.6 % of the saving hangs on which optimum the solver returns.
   Nothing yet picks among them by the true saving.

## 6. What was checked, and what was not

- The certificate of this revision is `tools/certificate/gcert1-e8-r783.json.gz`. The published checking
  programs accept it. Its labels are exactly the 120 E8 labels. Helper contents are whole numbers and halves,
  sources halves, targets sixths. For anyone who repeats the check:
  - sha256 of the file: `4b92f00fc7454b9b15e6d71a3b05062792eacca0a0ca4efee46f87827d57c8aa`
  - sha256 of the uncompressed JSON: `3594f19c4d01cf11fb930c5c61baeed399620acbc5b94096b16bcacb23997dd3`
  - the published `refcheck.py` prints `ACCEPTED by gx.check1: h=9 v=120 R=783 N=9039`
  - the mirror prints `MIRROR ACCEPTS`
  - `gxdry.py` prints `whole-block 8762479 (evaluated 8762482 holds / 8762483 fails)`
- A second agent checked the file before the Lean build, in three ways:
  1. A stand-alone replay and a check program written separately accept it with the same numbers.
  2. Seven copies with one fault each are each refused by both published checkers. The faults: a wrong
     in-place pair, a re-use without its compensating read, two helpers folded into one, a cheaper copy,
     exchanged labels, a wrong coefficient, cheaper stops.
  3. The making of the file was repeated from the pool of sums, with the same hierarchy and the same
     certificate byte for byte.
- Against the round-2 design that Lean checked during the work, the certificate has one new kind of gate
  that was introduced on purpose (the in-place pair on a helper that holds halves, 21 times), three kinds
  of step that came out of the same chain and that no earlier Lean build had (a source added into a helper
  with the coefficient +1/2, 4 times; a helper that climbs from a line to the full frame in one block, 2
  helpers; two additions there and back between two helpers with other signs, 5 times), and one parameter
  value that no earlier Lean build had (array tables of depth 10). The Lean text generated from it differs from that design's in numbers and
  in the data module only.
- What Lean checked is in the status lines under "What is new here" and in `VERIFY.md`.
- The Lean theorem states an exponent. It does not say "E8": the checkers accept, with the same saving, a
  copy of the circuit carried to another family of 120 labels. Only a separate test of the label set says E8.
- The optimum of the integer programme rests on the solver's report and on the constraints being the right
  ones. It is an optimum for a pool of 1,470 sums. It says nothing about hierarchies with other sums and is
  no lower bound on the helpers.
- The floors and every statement marked *on paper* are arguments by the agents with scripts beside them.
  They have had no human review and none is in Lean.
- The attributions to hub pull requests follow the descriptions of those pull requests and
  `RELATED-WORK.md`. The hub's commit history was not read for them.
- Whether a better file exists is open. Nothing here is claimed to be optimal.

## Sources

"Hub" is the community repository CrocSwap/integer-mult-bounds. Hub pull request n is at
`https://github.com/CrocSwap/integer-mult-bounds/pull/n`; for example, number 124 is at
<https://github.com/CrocSwap/integer-mult-bounds/pull/124>.

- Re-use of finished helpers with a compensating read: jamesyc, hub pull request 124; eumemic, hub pull
  request 143.
- Late pairing rule: Chafik Boukhalfa (account chafreaky), hub pull requests 200 and 233, and pull request 2
  to this repository.
- In-place pair: icekylinx, hub pull request 184; carried in ikeboy's hub pull requests 191 and 193.
- The circuit inside the certificate of the third revision: ikeboy (hub pull requests 191 and 193), eumemic (168) and
  icekylinx (144 and 184), as `RELATED-WORK.md` says.
- The two certificates priced for comparison at label width 22: pull request 2 to this repository (Chafik
  Boukhalfa) and hub pull request 304 (DreamingOfClouds).
- Dates: the work of this note was done on 2026-10-10; the search that flagged the family ran the night
  before.
