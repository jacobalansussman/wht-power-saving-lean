# The Reed-Muller family, and units that save a single move

This section reports two lines of research from 15 and 16 October 2026. The first asks whether one particular
family of small programs, which I call the Reed-Muller family, exists. If it did, its saving would tend to 1 as the
block size grows, far beyond anything known. The second line lowered the bar. It asked for any small program that
saves even one move, because a win there, however small, could teach how sharing works.

Neither line found a program. Both narrowed the question a great deal, and the programs that did the narrowing
are here, so the work can be repeated and pushed further.

I did this work with a team of AI agents (Claude) that I directed. I set the goals, chose which ideas to pursue and
decided at each stage what to do next. The agents wrote the programs, ran the searches and wrote the reports this
section is drawn from. Each main claim was then given to a separate AI agent (a "skeptic") whose only job was to
break it with its own code and its own reading of the proof. Most claims had one skeptic; the research rounds' own
rule asked for two. No person has read any of these proofs, and nothing in this section is checked in Lean.

## In one screen

- **No unit was found** in either line: no Reed-Muller program, and no program that saves even one move.
- **Where the Reed-Muller family is dead.**
  - At block 3, without phase letters, nothing saves at all, whatever the address maps (Theorem C3).
  - At block 4 the family cannot reach the bar in the engine game, the widened game, or any game inside them. The
    known lower bound allows at most 4 arrays on its 15 moves, and the bar needs 5.
  - At block 5, in its "natural form" (the form one would search first), 10 in 31 and 9 in 31 are impossible for
    every choice of different targets, and 8 in 31 for all but 127 of its 50,895 classes of eight different targets.
    Sets with repeated targets were only sampled.
- **Where it is still open.** In the full game at block 5 and above: address maps used in the middle of the program,
  or phase letters between moves. There, every test we have is blind, and one barrier shows that a whole kind of
  test can never close the question.
- **One move.** Where the searches can be trusted, the best two-array program one move short of plain still misses
  by 21 to 29 % of its target. Two arrays never save at block 4 when address maps act only at the ends (exhaustive).
  Symmetry sorts a search but does not create sharing.
- **What failed.** Several statements of the working notes did not survive their skeptics: side claims, coverage
  claims, and the reach of one set of symmetry exclusions. They are listed in section 4. Every exclusion in section 3
  that was given to a skeptic stood up to it; two short comparisons of results (RM4 and RM13) have no skeptic.
- **What is open**, each with a first experiment, is in section 6. The most useful single case is a 7-array problem
  at block 5 on which 24 of the 127 undecided classes hang.

## Words used here

- **Array, block, move.** The transform is computed by programs that work on many arrays of numbers at once. At
  *block m* a program works on m address bits of each array. A *paid move* (here simply a *move*) applies the
  transform along one *direction*, a nonzero m-bit vector, to a whole array. There are 2^m - 1 directions. The
  plain way costs m moves per array.
- **Unit "W in M", saving.** A unit serves W arrays with M moves. Its saving is ln(mW/M)/ln m. For example 8 in 31 at
  block 5 saves ln(40/31)/ln 5 = 1.58e-1.
- **The bar.** The long-standing aim of this research is a saving of 1e-1. At block m that allows m^0.9 moves per
  array: 2.69 at block 3, 3.48 at block 4, 4.26 at block 5.
- **Free steps and games.** Steps that cost nothing: shifts, *address maps* (relabelling the address bits by an
  invertible linear map) and *phase letters* (multiplying entries by powers of i). The *engine game* allows all
  three. The *widened game* also makes every address permutation (Toffoli gates included) and every pattern of
  fourth-root-of-unity phases free. The *no-phase game* is the engine game without phase letters. The *commutative
  slice* drops the address maps too. A unit in a smaller game is a genuine unit of the engine game.
- **Target.** What an array must end up with. Normally the transform itself; with address maps, an array can be owed
  the transform in a relabelled basis, written h h^T for an address map h.
- **Quadratic form, odd, tight, deletion.** The exact tests look at the program one binary quadratic form at a time.
  A form is *odd* for an array when that array's target needs a sign there. *Deletion* says that every form must
  act on at least as many moves as there are arrays odd at it. A form is *tight* when it acts on exactly that many.
- **The matrix E.** Every program has a matrix E built from its mixing coefficients, with one row and column per
  move. Its rank on the moves of a form equals that form's slack (moves minus odd arrays). Tight forms force blocks
  of zeros in E; the exact tests look for a contradiction between such zeros and the rank E must have elsewhere.
- **Certificate, cap, demand.** A certificate is our kind of lower-bound proof (Theorem R of the earlier rounds): a
  quantity that free steps never raise, which each move raises by at most its *cap*, while the targets need a total
  *demand*. Then moves >= demand / cap. *Hooks* and *lifts* are two families of certificates.
- **Singer cycle, placement.** Multiplication by a primitive element of the field with 2^m elements, used as an
  address map; it carries every direction to every other. Its *placement* is how it sits relative to the target;
  different placements give different problems.
- **Control.** A case with a known answer, run through the same test or solver. A known program must be found or
  passed (positive control); a case known to be impossible must be refused (negative control).

## 1. The question

**The Reed-Muller family.** Take one move in each of the 2^m - 1 directions and let W arrays share all of them. If
W can be as large as 2^(m-2) ("full sharing"), each array costs about 4 moves, whatever m is. The plain way costs m.
So at any fixed fraction of full sharing the moves per array stay constant while the bar allows m^0.9, and the
saving would tend to 1 as m grows. (The working notes call this family U_m.) Earlier rounds had already found that
the family cannot exist in the commutative slice at full sharing from block 4 on, and that every certificate we had
allowed 8 in 31 at block 5 in the engine game (those results, each with its status, are in
[../existence-routes/](../existence-routes/)). So I asked three things:

1. Does the family exist at blocks 4 and 5 once address maps are allowed? The members that matter are 4 in 15 at
   block 4 (saving 4.66e-2, the calibration case) and 8, 9 and 10 in 31 at block 5 (1.58e-1, 2.32e-1, 2.97e-1).
2. If not, can it be ruled out, and where?
3. Do phase letters change the answer?

**Units that save one move.** On 15 October I lowered the bar from 1e-1 to 1e-2 for one round. A unit that serves
W arrays with W·m - 1 moves saves one move. The targets were 2 in 7 at block 4 (9.63e-2), 2 in 9 at block 5
(6.55e-2), 2 in 11 at block 6 (4.86e-2), three to six arrays at block 4 (3 in 11, 6.28e-2, to 6 in 23, 3.07e-2), and
many arrays at block 3 under a symmetry group (5 in 14, 6.28e-2, to 21 in 62, 1.46e-2). For scale, the record for the
whole transform on the community's hub on 15 October was a saving of 1.6798e-3, so each of these targets was about 9
to 57 times larger.

## 2. What I tried, in order

1. **4 in 15, exactly, by address group** (15 October, the last numbered wave of the construction search). Decide
   4 in 15 at block 4 without phase letters, for each kind of address group: cyclic groups of order 3, 5 and 15, and
   the whole group GL(4,2).
2. **Round 1 of the Reed-Muller project** (15 October). Two tracks: build the Singer-invariant version of the family
   (the most symmetric placement), and work out the general-m theory of what address maps add.
3. **Round 2** (15 October). Three tracks: decide the family's natural form at block 5 exactly; search it with a
   calibrated solver; and study phase letters.
4. **Round 3** (15 and 16 October). Two tracks: the 1,266 target classes of 8 in 31 that survived round 2, and the
   smallest three-array shape next to plain at block 5 (3 in 14), by an exact algebra method that uses no theory.
5. **Tier 2: one move** (15 October). Three tracks: two arrays at blocks 4 to 6; three to six arrays at block 4;
   many arrays at block 3 under symmetry.

Each track wrote dated predictions before it computed, and each was followed by its skeptic.

## 3. What the work established

The state at the end, block by block. "Natural form" means: no phase letters, the Singer cycle in its canonical
placement, address maps only before an array's first move and after its last, one move per direction, and a
different target for each array.

| Block | Natural form | Address maps in the middle, no phase letters | Phase letters |
|---|---|---|---|
| 3 | no saving (Theorem C3) | no saving (Theorem C3) | at the ends: nothing new. Between moves: open (5 in 14, which is not a family shape) |
| 4 | 4 in 15 impossible for all 27 classes of four different targets | the family cannot reach the bar (engine and widened games) | at the ends: no unit above the bar with at most 4 arrays. Between moves: open |
| 5 | 9 and 10 in 31 impossible; 8 in 31 impossible except 127 classes; 3 in 14 next to plain impossible | open; every test is blind | at the ends: the same system as the natural form. Between moves: open |

Each finding below carries its status. "Checked by computer" means an exact computation (no floating point decides
anything) that an independent AI skeptic repeated with its own code, unless the line says otherwise. "Argued by
hand, not checked" means a written proof that an independent AI skeptic re-derived without finding an error; no
person and no proof checker has checked it. The identifiers (RM1 ...) match `results.json`.

### Address maps, and why block 3 is closed

**RM1. Address maps without phase letters are the same as giving each array its own target** (Theorem B,
Proposition R, Theorem MT). In the no-phase game with any group H of address maps, a commutative program in which
array l is owed h_l h_l^T instead of the plain target is already a unit with the maps used only at the two ends.
Conversely, a unit with group H becomes such a program with |H| times as many arrays and moves. Both reach the same
moves per array. So every no-phase question is a small system of scalar equations (a 3-in-11 trial takes about 2
seconds), and 8 in 31 becomes a system of 961 unknowns instead of 29,791.
*Status: argued by hand, not checked.* Three AI skeptics re-derived parts of it; one also tested the whole chain by
computer on a real unit with address maps (1 in 4 with the order-15 Singer group): the identity holds at all 1,023
forms, and with the targets put on the wrong copies it fails at 1,016, so the check can tell the difference.
Proofs: [proofs/address-maps-per-array-targets-and-block3.md](proofs/address-maps-per-array-targets-and-block3.md)
(section 8.1), [proofs/four-in-15-channels-and-rank-depth.md](proofs/four-in-15-channels-and-rank-depth.md).

**RM2. Address maps that keep the dot product never help, and one array always needs m moves** (Theorem O,
Proposition D). If every map h has h h^T = I, every copy has the plain target and the slice's theorems apply. So a
unit must use maps that change the dot product. And with any address maps and no phase letters, one array needs m
moves, because after k moves every register involves only matrices of rank at most k. So 1 in 3 at block 4 is
impossible, which makes it a clean negative control for any test.
*Status: argued by hand, not checked.* One skeptic read Proposition D line by line; Theorem O follows directly from
RM1 and was not checked on its own.

**RM3. Theorem C3: at block 3, no phase letters means no saving, whatever the address maps.** A block-3 certificate
found in an earlier round (the "face" certificate) has cap 8 for every direction and demand 24 in every one of the
28 bases of the address space. Through RM1 this gives M >= 3W for every address group. So any saving at block 3
must use address maps and phase letters together; either kind alone gives exactly the plain cost. It also closes
the Singer-invariant 5 in 14 at block 3 that round 1 had left open: its search floor belonged to an impossible shape.
*Status: argued by hand, not checked (one skeptic).* The certificate's numbers are checked by computer by two
separately written programs (ranks modulo two primes, and exact arithmetic over Q(i)); I re-ran both for this
revision and got the same output. The same averaging over the 840 bases at block 4 gives only 2.706 moves per array,
below the 41/13 = 3.154 that is already known there, so at block 4 it is a method, not a new bound.

### Block 4

**RM4. At block 4 the family cannot reach the bar in the engine game, the widened game, or any game inside them.**
One move per direction means 15 moves.
The bar needs 15/W <= 3.48, so at least 5 arrays. The block-4 lower bound of 41/13 = 3.154 moves per array (an
earlier result for the engine and widened games, checked by computer with certificates modulo a prime) allows at most
4 arrays on 15 moves. The bound is not known to hold in the everything-free game of the moonshot sections.
*Status: argued by hand, not checked.* It is two lines of arithmetic from the 41/13 bound; no skeptic read the
comparison itself.

**RM5. With address maps and phase letters only at the ends, block 4 has no unit above the bar with at most 4
arrays.** Exhaustive tables over every multiset of directions and every multiset of targets found no 2 in 7, 3 in 10,
3 in 11 or 4 in 13 (and no 3 in 8 at block 3). A second track repeated this with every invertible symmetric target,
the alternating ones that only phase letters can owe included, so it holds with phase letters at the ends too. With
the engine's bound, 3 in 10 and 4 in 13 are the only above-bar shapes with at most 4 arrays, so none exists in this
form. For pairs, a certified bound needs at least 56/15 = 3.733 moves per array for every one of the 20 kinds of
target pair, against the 3.5 that 2 in 7 would need.
*Status: checked by computer.* The skeptic of the pairs result enumerated all 116,280 multisets of 7 directions and
found none that passes deletion for any pair; I re-ran that program for this revision (20 seconds) and got the same
counts. I also re-ran the exhaustive tables for 3 in 8 at block 3 and 3 in 10 at block 4. The scope, set by a
skeptic: "at the ends" means letters on the inputs, a commutative program that holds all the mixing, and letters on
the outputs. Ends that mix arrays between their letters are not covered.
A side result, *argued by hand, not checked (one skeptic)*: phase letters used only at the ends add nothing at odd
blocks, because an invertible symmetric binary matrix is congruent to the identity unless its diagonal is zero, and
such alternating targets exist only at even blocks (28 at block 4, 13,888 at block 6).

**RM6. 4 in 15 without phase letters, group by group.** Excluded for 3 of the 5 kinds of fixed-point-free group of
order 3 (42 of their 56 groups) and for one kind of order 5 (24 groups); the linear-programming bound alone excludes
it for 76 of the 616 groups of order 3. Theorem E, half a page, excludes the evenly spread case of one order-3 group:
tight forms force zeros in E, and the form x1x2 + x1x3 then needs rank 8 from a block that holds at most 6. 2 in 7 and
3 in 11 are excluded in many of the same classes (2 in 7 for 4 of the 5 kinds of order 3 and 7 of the 10 kinds of
order 5; 3 in 11 for 4 and 3 of them). But with the order-15 Singer groups and with all of GL(4,2) the
tests are blind: at order 15 they even pass the impossible 1 in 3 (RM2) in 6 of the 10 kinds of placement, and with
GL(4,2) the linear program stops at 2.035 moves per array.
*Status: checked by computer.* The skeptic rebuilt all four exact tests with its own code; all 45 rows (15 kinds of
group, 3 shapes) agree, survivor lists included; it rebuilt the linear program and confirmed the 76 exclusions, and
wrote exact certificates in fractions for one group of each class that the linear program excludes. I re-ran its check of Theorem E and its count of placement classes for this revision; both match. A
small-group exclusion says nothing about a larger group that contains it.

### Block 5, the natural form

**RM7. 10 in 31 and 9 in 31 are impossible, and 8 in 31 is down to 1,266 classes** (round 2). Every set of
different targets was tested, one per class (the classes agree with Burnside's count):

| Shape | Saving | Classes | Excluded by | Left |
|---|---|---|---|---|
| 10 in 31 | 2.97e-1 | 286,155 | deletion 35,732; rank cover 248,363; star test 2,060 | 0 |
| 9 in 31 | 2.32e-1 | 130,065 | deletion 2,288; star test 127,501; zero-clique test 276 | 0 |
| 8 in 31 | 1.58e-1 | 50,895 | the star tool with pooling 49,629 | 1,266 |

The decisive new tool is the *star tool*: for each form that no array needs, it cuts that form's moves into pieces
along affine subspaces, uses a rank inequality (the star lemma) to show that some blocks of E must be zero, pools
the zeros of all forms, and finishes with two old facts (E can be zero on at most W moves at once; the rank on a set
of moves is at most the largest matching that avoids the zeros). Also: all 27 classes of 4 in 15 at block 4 with
four different targets are impossible, and all 300 random sets with repeated targets were excluded.
*Status: checked by computer.* The skeptic wrote its own field, forms, flats and recursion and verified every 9-
and 10-in-31 exclusion and 1,687 of the 8-in-31 ones; 47,942 of the 8-in-31 exclusions rest on one program. The
plain programs and the one known no-phase program with real sharing pass (U_3 = 2 in 7 at block 3: the family's own
block-3 member, two arrays on all seven directions; it costs one move more than plain, but its arrays really share
moves); the known impossible cases are refused. At block 5 there is no positive control with real sharing, so the zeros the tool
derives are tested against a genuine program only at block 3. For this revision I regenerated the tool's tables
(byte-identical), re-ran it on round 1's proposed first set (excluded, as recorded) and on three sampled classes and
one block-4 class (the same verdicts as recorded), and re-ran the skeptic's bookkeeping program on the packaged
outputs (identical counts).

**RM8. Round 3: 1,139 of the 1,266 are excluded, none is feasible, 127 are undecided.** The new idea is to bound
the rank of blocks of E from below as well as from above: if the derived zeros leave little room outside a small
block, that block must carry a lot of rank, and often it cannot. A worked example: in the first survivor, an even
form acts on 12 moves, so E has rank 12 there; six of those moves lie in a form of slack 3, so they carry rank at
most 3; outside that block the zeros leave a term rank of 8; and 12 <= 3 + 8 fails.
Of the 127 left, 24 are *hyperplane classes*: seven of the eight targets sit on a 3-dimensional subspace, no zero of
E can be derived, and every zero-based test is blind. All 24 hang on one 7-array problem (open question O1).
*Status: checked by computer.* The skeptic's program shares no code with the track's and excluded all 1,139 again,
so every round-3 exclusion has two programs behind it. The skeptic's program is stronger and also excludes 31 of the
127 (none of them a hyperplane class); only one program stands behind those 31, so they are not counted. If they
hold, 96 classes are left. I re-ran the first stage on two survivors and a control for this revision (same verdicts,
15 seconds).

**RM9. Sharing next to plain.** A shape "next to plain" is a plain program with one move removed. Two arrays next
to plain are excluded at blocks 3, 4 and 5. Three arrays at block 5 (3 in 14, 4.29e-2): all 44 candidate sets that
pass span and deletion are impossible, by two independent exact methods. Each fails on 5 moves: three of them pair up
into tight forms, so E has a 3 x 3 block of zeros, while a form even for all three arrays makes E invertible on the
5 moves. And, with no theory at all, for every set and every time order, 1 lies in the ideal of those moves' 63
equations in 40 unknowns over the rationals (Singular over Q and msolve over two primes agree). These equations are
part of the full system, so this rules out every time order of the full program, even over the complex numbers.
*Status: checked by computer.* The skeptic built its own model without a self-dual basis, found exactly the same 44,
ran all 120 time orders of each (5,280 in all) and multiplied out every cofactor certificate exactly. I re-ran the
zero-block test on all 44 sets in both coordinate systems for this revision: all excluded, and the certificate file
came out byte-identical. The Gröbner step was not re-run (it needs Singular or msolve; the systems are in
`data/three-in-14-groebner/systems.tar.gz`).

### Phase letters

**RM10. With phase letters between moves, only one kind of test sees anything, and it cannot refuse full sharing**
(Theorem T, Proposition H, "the angle barrier"). With phase letters free, the engine game maps onto a game on the
symplectic group Sp(2m, 2). Every certificate carries over with the same strength, and back, so the known bounds hold
there (3.154 moves per array at block 4, 2.992 at block 5), while deletion and the rank tests are blind. A
certificate's strength is a spectral factor, at most 2 - 2^(1-m), times an angle factor; the family's full-sharing
member costs exactly twice that spectral ceiling, at every block. So refusing it needs an angle factor above 2, and
no hook has one: no hook can ever refuse the full-sharing member, at any block, with or without phase letters. No lift
can refuse 4 in 15 either (it would need excess below 1/15 per edge; every lift has at least 0.075). At block 5, a
lift of the whole frame graph (frames are the coordinate systems an array can be written in) with excess below 0.290
per edge would refuse 8 in 31 in every game (below 0.452 for 9
in 31, below 0.613 for 10 in 31); the best such lifts known have 0.88, and nothing rules out a better one.
*Status: argued by hand, not checked (one skeptic).* The skeptic did not check the inputs taken from other rounds
(the hook bound, the lift floor 0.075, two facts of an earlier induction). Calibration: the plain schedule and the E8
unit pass; every block-2 saving and 2 in 5 at block 3, all impossible in the engine game, are refused.
Proofs: [proofs/phase-letters-and-angle-barrier.md](proofs/phase-letters-and-angle-barrier.md).

### Structure of the Singer placement

**RM11. The canonical Singer placement, and why symmetry does not split the problem.** Write the addresses in a
self-dual basis of the field with 2^m elements; then the Singer matrix is symmetric, each move reads a quadratic form
along the Singer orbit, and the target reads only its linear part. The copies' targets are the field's own
multiplication matrices, the configuration behind the Kerdock codes and complete sets of mutually unbiased bases.
The pieces of the program that carry the target are irreducible, of dimension 2^m - 1, so the hoped-for Fourier
split never reaches them. In every Singer placement a form is odd on at most 3 times as many copies as its weight, so
deletion and every test built on tight forms say nothing about a shape with more than 3 moves per array; the
canonical placement is where every test is weakest. Adding the Frobenius map changes nothing about savings (this last
point comes from comparing the tracks' notes and has no skeptic).
*Status: argued by hand, not checked*, with the structure checked exactly by computer at blocks 2 to 5 (some parts
also on samples of forms up to block 7). One skeptic per track.
Proofs: [proofs/singer-invariant-ansatz.md](proofs/singer-invariant-ansatz.md),
[proofs/address-maps-per-array-targets-and-block3.md](proofs/address-maps-per-array-targets-and-block3.md).

### One move

**RM12. No two-array unit saves at block 4 without phase letters when address maps act only at the ends** (part of
RM5), and three lemmas restrict where a one-move saver can share. Lemma S: if every chain of m moves passes
through one shared move, the program is not a unit; two m-move chains that share only their first move, or only
their last, never make one, so sharing must happen in the middle. Lemma B(k): in the natural one-move shape, both raw
inputs must be fed in early. Lemma Z: in the lifted block-5 game (the Singer cycle on four bits, the fifth left
alone), moves in the fifth direction can never be shared without moves that mix it with the others, so without them
a block-5 2 in 9 exists exactly when a block-4 2 in 7 does.
*Status: argued by hand, not checked.* The skeptic checked each lemma line by line; their one-line summaries go a
little further than the proofs.

**RM13. Block-4 triples are settled in the form tier 2 left open.** Tier 2 listed 58 kinds of target triple that
pass deletion at 10 moves and 452 at 11 as its cheapest next step. Comparing the notes the next day showed that the
exhaustive tables of RM5 already leave no 3 in 10 and no 3 in 11 at block 4 when address maps and phase letters act
only at the ends, which is exactly the setting of that step (by RM1).
*Status: argued by hand, not checked.* The comparison itself has no skeptic.

**RM14. Symmetry sorts a search but does not create sharing.** Split a symmetric program by the group's characters
(a Fourier transform over the group, which is free). If the group moves every move and every array around freely,
the symmetric unit is just a smaller ordinary one; for example a 7-fold symmetric 35 in 98 exists exactly when 5 in
14 does. Moves or arrays that part of the group leaves in place tie the pieces together, and they are expensive: with
the 7-fold group and no freely moving arrays, b fixed arrays need at least 7(b + 2) moves, against 3b plain. At
block 3 no free letter has order 5, so a 5-fold symmetry can only split a program into smaller ones.
*Status: argued by hand, not checked* (the core lemmas, one skeptic). The reach of the block-3 exclusions built on
this was refuted in part (section 4, item 4).
Proofs: [proofs/one-move-block3-symmetry.md](proofs/one-move-block3-symmetry.md).

## 4. What failed, or was withdrawn

Failures are as useful as successes, so here are the statements of the working notes that did not survive. Where
notes disagree, the later skeptic's verdict is the one reported in section 3.

1. **"Target sets with no tight form pass every exact test"** (round 1). *Refuted* by round 2: the star tool
   excludes them, round 1's proposed first search set and every 9- and 10-in-31 set with different targets included.
   Also from round 1: **"no tool we have sees phase letters"**, *refuted* by round 2 (certificates of Theorem R do;
   deletion and the rank tests are the blind ones), and two slips (9 in 31 saves 2.32e-1, not 2.16e-1; there are 840
   bases at block 4, not 420).
2. **The search track of round 2 was refuted as written; its search results stand.** Its two "negative controls that
   pass every exact test" (3 in 8 at block 3, 3 in 11 at block 4) fail the zero-block test. Its 44 three-array sets
   next to plain at block 5 were not open: all fail the zero-block test (RM9 then added a second proof). Its "716 or
   more" four-array sets came from a scan stopped early; the full scan has 15,139, and only 9 pass every test. The
   one candidate 8 in 39 next to round 1's first set is excluded exactly, by term rank. "No no-phase program with
   real sharing is known" is false: U_3 = 2 in 7 at block 3 is one, and the solver finds it. "Serving an array partly
   never lowers the residual" is too strong. Each was found by its skeptic.
3. **The two-array track of tier 2 overstated its coverage.** The Singer game holds every program for one placement
   only; there are 10 placements at block 4 and 92 at block 5, and the track searched 2 at each. The skeptic searched
   the other 8 at block 4 and found the same band (22 to 27 % of the target). The track's closest point was also not
   a resting point: a more careful method kept lowering it (230.76 to 230.24) while the coefficients grew. "Nowhere
   near a unit" still holds; "a local minimum with bounded coefficients" is *refuted*.
4. **The block-3 symmetry exclusions of tier 2 hold only under exact symmetry.** The track's rule "the number of
   moves must be a multiple of 7" is *refuted* for gauge symmetry (a different free letter for each stored array, the
   more natural notion): there one move can serve a single character on its own. Then 20 in 58, 20 in 59 (7-fold) and
   15 in 42 (21-fold) pass every screen the track had (open question O6). The same caveat applies to its block-4
   by-product (no 5-fold symmetric 5 in 19), which nobody has checked under gauge symmetry.
5. **Smaller corrections.** Only one kind of order-15 group (not two) is generated by a symmetric matrix. The
   order-3 solver control of wave 8 succeeded in 2 of 4 restarts, not 2 of 3. Tier 2's block-4 track counted 2,059
   kinds of triple that merge into 758 under the game's full symmetry, and about 124 recorded trials, not 150. An
   earlier claim that "address maps make the slice's exact tests blind" was *refuted* in part by wave 8: blindness
   depends on how the address group sits relative to the target, and only the order-15 Singer groups escape in every
   placement.

## 5. Searches that found nothing (evidence, not proof)

Every numerical search ended far from a unit. These are evidence only: a local solver cannot tell "impossible" from
"possible but isolated", and a unit at a saving would sit in special position. *Status of each: open.*

- **The natural form at block 5**, with a solver that does find the plain programs there (8 in 40 in 3 of 6 random
  starts) and the real-sharing control U_3 in 12 of 12: best 8-in-31 residual 3.08e4, where a unit has 0 and one
  unserved array costs about 1e4. Every stopping point serves some arrays the plain way and leaves the others nearly
  unserved. Every set it searched with different targets is now excluded by exact computation (RM7), and its floors
  on those sets look just like its floors on the open ones.
- **The Singer-invariant form at block 4**: 2 in 7 best residual 239.3, where the trivial split (one array served on
  4 moves, the other given the best of 3) costs 294.6. The arrays do share moves, but not enough.
- **4 in 15 by address group** (an exact-model solver): the order-3 targets never came near zero, while the solver
  reached its controls; at order 15 every 2-in-7 run stopped near the floor of 1 in 3 (one array served, the other a
  move short); at order 5 its control was found only once in 3 tries, so those runs say nothing.
- **One move, two arrays**: block 4, 2 in 7, best cost 230.7 of a target of 1,024 (22.5 %); lifted block 5, 2 in 9,
  0.21 to 0.26 of the target. At block 5 with the 31-step cycle and at block 6 the solver cannot find the plain
  programs, so nothing is known there.
- **One move, block 4, three and four arrays**: about 124 trials of 3 in 11, every floor 286 or more, that is at or
  just below the cost of "one array a direction short" (288; best 286.12); and 16 trials of 4 in 15, every floor 384
  or more. The skeptic's own 64 trials of 3 in 11 found nothing below 286.7.

The lesson the rounds drew (calibration decides what a floor means): a floor is evidence only when the same solver,
from the same kind of start, finds the known programs of that size; steer starts away from the trivial split, on
which every solver here sat for thousands of steps. Random-start searches on 8, 9 and 10 in 31 were stopped.

## 6. What is still open, with a first experiment each

**O1. The 127 undecided 8-in-31 classes, and the 7-array problem.** 24 of them hang on one problem: seven arrays at
block 5 whose targets sit on a 3-dimensional subspace, up to symmetry the single target set 0, 1, 2, 5, 11, 18, 19
(7 in 31, 7.54e-2). If it is impossible, all 24 fall; if it is possible, it is the first block-5 program with real
sharing, and the positive control the zero-based tests lack. Its arrays' characters are all the nonzero characters
of a 3-bit quotient, which no test uses yet. A later idea round (16 October, in the idea-rounds section) argued that
there is no block-4 shortcut: without scratch arrays and without interleaving, the moves that are not simple sign
flips of a single array must span all five directions (argued by one AI agent with a short computation, not checked
by another).
*First experiment:* run the local Gröbner test of RM9 (`three-in-14-groebner/gb2.py`) on the small move sets of the
7-array problem, every time order, starting with sets of 5 to 7 moves on which the most forms act. *Second:* confirm
or drop the skeptic's 31 extra exclusions with a second implementation of its three rules (`vr3.c` in
`eight-in-31-rank-bounds/independent-check/`); if they hold, 96 classes are left.

**O2. The rest of the natural form at block 5.** Not yet scanned: target sets with a repeated target (8 in 31 with at
most 7 different targets, about 235,000 classes, only sampled), the other 91 placements of the Singer cycle, larger
address groups up to GL(5,2), and general targets congruent to the identity. A second check of round 2's 47,942
one-program exclusions is also owed. *First experiment:* run `starrec2` over the repeated-target classes (a few hours
on an 8-processor machine), and run the round-3 skeptic's program `vr3.c` over all 50,895 classes as the second check.

**O3. Four arrays next to plain** (the smallest open sharing shapes). One set at block 4 (4 in 15, 4.66e-2) and nine
at block 5 (4 in 19, 3.19e-2) pass every exact test; a probe for small move sets where the Gröbner test bites found
none (block 4 fully, block 5 for at least 7,813 of 18,170 candidates). The later idea round found that the block-4
set cannot exist without scratch arrays, by a determinant count (argued by hand by one AI agent and re-derived by a
second; no person and no proof checker has checked it); with scratch it is open. Either answer helps: an impossible set is a negative control the old tests cannot see; a possible
one is the first control with real sharing. *First experiment:* write the full exact system of each of the nine block-5
sets for one time order and run it through Singular, or extend the star tool to repeated directions.

**O4. The general game: address maps in the middle, phase letters between moves.** This is the game that decides
savings. With address maps in the middle, RM1 turns a unit into a program 31 times larger in which every target is
used equally often, and the zero-based tests have nothing to feed on. With phase letters between moves, only
certificates see anything, and RM10 limits them. A certificate that works here must use what today's tests drop:
coefficients that commute with relabelling the copies, at most one move per register per physical move, and rank
depth. *First experiment:* add "a register is moved at most once per physical move" to the multi-target rank-cover
linear program of RM6 (`four-in-15-by-address-group/mtlp.py`) and check that it then refuses 1 in 3 in all ten
order-15 placement classes, which today's tests cannot. *For phase letters:* search for lifts at block 5 with excess
below 0.290 per edge among lifts unchanged by an odd-order subgroup of Sp(10, 2) of order 31 or 33, which shrinks the
search 31- to 165-fold; such a lift would refuse 8 in 31 in every game.

**O5. A proof route without phase letters at block 4.** A block-4 certificate whose average demand over all 840 bases
exceeds 3.482 times its average cap would rule out every no-phase unit above the bar at block 4, whatever the address
maps (above 3.75 it rules out 4 in 15; 4 closes block 4 without phase letters). The same at block 5 above 3.875
rules out 8 in 31 without phase letters. *First experiment:* rerun the earlier rounds' plane-subspace certificate
search at block 4 with its objective changed to (average demand over the 840 bases) / (average cap over the 15
directions); `address-maps-and-small-lift/t_cert4.py` computes the average for a given certificate.

**O6. Block 3 under gauge symmetry.** Four shapes pass every screen: 20 in 58 (3.09e-2) and 20 in 59 (1.53e-2) with a
7-fold symmetry, 15 in 42 (6.28e-2) with a 21-fold one, and the 5-in-14 shapes, whose list under gauge symmetry was
never enumerated. Each would need address maps and phase letters together (RM3). A 7-fold 15 in 42 comes down to two
6-move problems, each exactly tight on the hook bound; one lemma half-closes one of them. *First experiment:* decide
the two 6-move problems exactly, with one representative per orbit.

**O7. Pairs with phase letters, and pairs at block 5.** 2 in 7 at block 4 is excluded without phase letters when
address maps act only at the ends (RM5), but not with phase letters or with address maps inside. *First experiment:*
the exact formulation of 2 in 7 at block 4 with phase letters, through the Sp(8, 2) game of RM10, filtered first by
the two-array conditions of the frame-packing section. Before any block-5 pair search, a solver must find the plain
2 in 10 in the 31-step Singer game; none does yet.

**The question underneath.** Without phase letters, the family lives or dies on one question: do different targets
ever let arrays share better than one target does? Every case decided so far says no: block 3 in full, and block 4
with up to three arrays, and four arrays on 13 moves, when address maps act only at the ends. Nothing decides it from block 4 on in general, and the
cheapest above-bar test of it is 8 in 31.

## 7. How to run the programs

Everything here is Python 3.9 or later and C. Most programs need numpy; some need scipy (the linear programs use
HiGHS through scipy's `milp` and `linprog`); the two numerical solvers of `singer-invariant-search` and
`one-move-two-arrays` need torch. The Gröbner step needs Singular and msolve, which these programs only call. The C
programs build with `cc -O2 -o NAME NAME.c`.

The programs read and write their data files in the current folder, and a few read a sibling folder. So first build
a working copy:

```
sh programs/workdir.sh          # makes work/ with programs/<name> and data/<name> merged, archives unpacked
cd work/<name>
```

Each folder of `programs/` has a README with its commands, times and expected outputs. In brief:

| Folder | What it shows | Light command (re-run for this revision) | Heavy part (not re-run) |
|---|---|---|---|
| [address-maps-and-small-lift](programs/address-maps-and-small-lift/) | RM1, RM3, RM5 tables, RM11 | `python3 t_face3.py` (6 s), `independent-check/vk_face3.py` (14 s), `python3 t_small.py 4 3 10` (15 s) | 3 in 11 and 4 in 13 tables, the 2.706 transfer (hours) |
| [four-in-15-by-address-group](programs/four-in-15-by-address-group/) | RM2, RM6 | `independent-check/thmE.py` (<1 s), `independent-check/vchain.py classes` (40 s) | every split for every class; the solver |
| [natural-form-star-tool](programs/natural-form-star-tool/) | RM7 | `gen_tables.py 5`, `cc` `starrec2.c`, `./cal2.sh < input` (about a minute for an excluded class, several for one that passes); `independent-check/coverage.py` (10 s) | the full scans (many hours on several 8-processor machines) |
| [eight-in-31-rank-bounds](programs/eight-in-31-rank-bounds/) | RM8 | `starrec4` on a few classes (15 s) | all stages on all classes; the skeptic's cascade |
| [three-in-14-groebner](programs/three-in-14-groebner/) | RM9 | `python3 run_z.py` (80 s) | the Gröbner bases (Singular, msolve) |
| [phase-letters](programs/phase-letters/) | RM5 (alternating targets), RM10 | `python3 thresholds.py` (<1 s), `python3 t_setup.py` (4 min) | the 4-in-15 enumeration with alternating targets |
| [singer-invariant-search](programs/singer-invariant-search/) | RM11, section 5 | `python3 struct_check.py` (76 s) | the solver runs |
| [natural-form-solver](programs/natural-form-solver/) | section 4 item 2, section 5 | `python3 packing2.py` (about 1 min) | the solver runs |
| [one-move-block4](programs/one-move-block4/) | RM5 (pairs), RM13, section 5 | `independent-check/vpairs.py` (20 s) | the triple searches |
| [one-move-two-arrays](programs/one-move-two-arrays/) | RM12, section 5 | none | the solver runs |
| [one-move-block3-symmetry](programs/one-move-block3-symmetry/) | RM14, O6 | the 13 screens `s1_levi.py` ... `s14_fixedmove_check.py` (seconds to 90 s each) | none |
| [shared](programs/shared/) | genuine U_3 programs used as controls | | |

The heavy runs of the rounds were done on Linux machines of Google's Colab service (8 processors, about 50 GB each).

**Re-runs for this revision** (on an Apple M4 Mac mini with 16 GB, each under a 3 GB memory cap): the Theorem C3 certificate
(both programs, identical output); the star tool's tables (byte-identical), round 1's first set, three sampled
8-in-31 classes and one block-4 class (same verdicts); the skeptic's bookkeeping over the packaged round-2 outputs
(identical); round 3's first stage on two survivors and a control (same verdicts); the zero-block test on the 44
sets of 3 in 14 in both coordinate systems (all excluded, certificate file identical); the wave-8 skeptic's Theorem E
check and placement classes (same); the exhaustive pair enumeration (same); two exhaustive small-lift tables
(same, apart from timings); the phase-letter thresholds and the exact checks of the symplectic game (same); the exact
checks of the canonical Singer form (same); the frame count for three targets (6 frames each, as stated); and the 13
block-3 symmetry screens (same printed results; two recorded files also hold numerical warnings that the original
runs captured). Everything else was not re-run for this revision; the outputs in `data/` are the ones the
rounds produced.

## 8. Method warnings

- The exclusions of RM7 and RM8 that use derived zeros are tested on a genuine program only at block 3, because no
  block-5 program with real sharing is known (O1 and O3 could supply one).
- A pass of any exact test here means only that the test found no contradiction. At order 15 and with GL(4,2) the
  tests pass a case known to be impossible (RM6), so a pass there means nothing.
- msolve version 0.6.5 reads a product x*x as x. The msolve inputs of this section were written directly by their
  generator, with powers written as x^k, and in these systems every monomial has distinct variables anyway; no
  converter was used. (A converter of a different round, `sing2ms.py`, corrupts entries that are sums of two terms; it
  is not used here.)
- The outputs keep the working labels of the rounds in their file names and comments; the folder names here are
  descriptive. The map between them is in Sources.

## Sources

Notes and verdicts this section is drawn from (all in the private research archive, under `checks/wht43/`; not part
of this repository):
- The Reed-Muller project: `rmfam/NOTE.md` (rounds 1, 2 and 3, with the skeptics' corrections); the tracks
  `rm-singer`, `rm-theory` (round 1), `rm2-exact`, `rm2-solve`, `rm2-phase` (round 2), `r3/r3-classes`, `r3/r3-3in14`
  (round 3), and their skeptics `v-rm-singer`, `v-rm-theory`, `v-rm2-exact`, `v-rm2-solve`, `v-rm2-phase`,
  `r3/v-r3-classes`, `r3/v-r3-3in14`.
- Wave 8's exact 4 in 15: `moonshot/w8-rmexact/NOTE.md` and `PROOF.md`; skeptic `moonshot/v-w8-rmexact/NOTE.md`.
- Tier 2: `tier2/NOTE.md`; tracks `t2-pair`, `t2-block4`, `t2-block3`; skeptics `v-t2-pair`, `v-t2-block4`,
  `v-t2-block3`.
- The later verdicts on the 7-array problem and the block-4 four-array set: `checks/wht44/ideas/skeptic-2/NOTE.md`
  and `e7-rm/NOTE.md` (16 October).
- Background results used here: the commutative-slice theorems and Theorem 10 (`exist2/ex2-families/PROOF.md`), the
  star lemma (`exist3/e3-slice6`), the exact model of the engine game without phase letters (`exist3/e3-rmbuild`,
  Theorem A), the block bounds 2.8, 41/13, 2.99172 (`proofs/NOTE.md`, `wave16`), Theorem R (`checks/wht35/wider`),
  the face certificate (`checks/wht38/tiny/close`).
- The hub's record of 15 October quoted in section 1: `anatomy/NOTE.md`, as cited in `tier2/NOTE.md`.
- The private results catalogue's entries: L21, L22, L23, L24, L25, X36, X37, X38, X39, X42, N7, N14, N15, N16, N17,
  N19, T19, T20, M29, C27, C28, C29.

Folder map (descriptive folder here, working label of the round):
`address-maps-and-small-lift` = rm-theory and v-rm-theory; `singer-invariant-search` = rm-singer and v-rm-singer;
`four-in-15-by-address-group` = w8-rmexact (solver: its `clm/`) and v-w8-rmexact (linear program: its job folder jL);
`natural-form-star-tool` = rm2-exact (second programs: its `indep/`) and v-rm2-exact; `natural-form-solver` =
rm2-solve and v-rm2-solve; `phase-letters` = rm2-phase (block-4 enumeration: its `tile/`) and v-rm2-phase;
`eight-in-31-rank-bounds` = r3-classes and v-r3-classes; `three-in-14-groebner` = r3-3in14 and v-r3-3in14;
`one-move-two-arrays` = t2-pair (its `t2p-singer/`) and v-t2-pair; `one-move-block4` = t2-block4 and v-t2-block4;
`one-move-block3-symmetry` = t2-block3 and v-t2-block3; `shared` = two modules of ex2-families.

Tools: Singular, msolve, HiGHS (through SciPy), NumPy, PyTorch.

Not included, with reasons: the raw logs and remote job folders of the rounds; compiled binaries; three 5 MB linear
program row dumps of tier 2's block-4 track (regenerated by `certify_t.py`); the skeptic's 49 MB dumps of bounds on
genuine block-4 units (regenerated by `independent-check/u4mt.py` and `u3check.py` in `eight-in-31-rank-bounds`); the
duplicate copies of the round-2 outputs that the skeptic fetched. Large outputs are in the compressed files listed in
`data/SHA256SUMS`.
