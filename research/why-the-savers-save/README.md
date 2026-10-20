# Why the record savers save, and what caps them: anatomy, centre, copy wall, helpers

This section covers the rounds of 13 to 15 October 2026 that took apart the circuits behind the record figures,
OpenAI's, the community's and mine, to find out where their saving comes from, what limits it, and whether a
different centre, decoder or helper circuit could lift the record by a few per cent. I did this work with a team of
AI agents (Claude) that I directed: I set the goals, chose which ideas to pursue and decided at each stage what to
do next; the agents ran the searches, wrote the proofs and programs, and tried to refute one another's work. No
person has read any of the proofs in this section, and nothing in it is checked in Lean. (The record figures it
dissects, such as U8X0HQ, are Lean-checked; they are packaged elsewhere in this repository.)

**In one paragraph.** Every circuit that saves anything, OpenAI's unit, the community's record unit T74 and my
Lean-checked U8X0HQ, saves in the same way: each data array gets one of its directions free, because its data ends on a
different array whose frame is the complement of the line it started on, and a small shared *centre* that works like an
index code carries the data across and takes part of the gift back. Net of the centre, each data array in the
community's template saves about 0.36 of a move out of 73, and that number has not changed since the circuit that Ryan
Shea checked in Lean; every later gain, the community's and mine, came from grouping the paid work into larger pieces. I
then tried every lever on that mechanism that the rounds could think of: a second free direction per array, other port
designs (E7, Barnes-Wall, Leech, Coxeter-Todd), smaller and medium blocks, a cheaper centre, cheaper copies, fewer
helpers, symmetric decoders and rewriting the community's best circuit. None of them lifted the record: the one helper
that could be removed sits on the community's bit side, which does not set its headline figure, and the helper pairs
that could end one step narrower on the community's newest circuit had already been claimed by other contributors.
Several gave new lower bounds instead, the strongest being that no centre of the community's 112-port design that keeps
today's junk can cost less than 64 per call, for any number of terms. Several claims made on the way were refuted, and
they are listed as plainly as the results. A handful of levers are still open; each comes with a first experiment in
section 6.

## Contents

1. [The question](#1-the-question)
2. [Words used here](#2-words-used-here)
3. [What I tried](#3-what-i-tried)
4. [What came out](#4-what-came-out)
5. [What failed or was refuted](#5-what-failed-or-was-refuted)
6. [What is still open, with a first experiment for each](#6-what-is-still-open-with-a-first-experiment-for-each)
7. [The programs and how to run them](#7-the-programs-and-how-to-run-them)
8. [Sources](#sources)

## 1. The question

The Fourier and Walsh-Hadamard records in this repository, and on the community's leaderboard, all come from one
family of circuits. OpenAI's preprint on an explicit power saving for the exact discrete Fourier transform built
the first one; the community (the CrocSwap/integer-mult-bounds repository and its contributors) shrank it by many
orders of magnitude, took up the E8 unit that my agents found in October, and then improved it step by step; my Lean
records since the fourth revision are the community's circuits, checked in Lean on Ryan Shea's tree with overhead
removed. On the morning of 15 October the rounds compared against T74's open figure 1.6797550e-3; by that evening the
community's reviewed record was 1.6727e-3 and its newest open figure 1.7173e-3. Mine was U8X0HQ's 1.6256629e-3.

So I asked three questions:

1. **Why do these circuits save at all?** Is there one mechanism behind all three, and what exactly is it?
2. **What caps it?** How large could the saving of this mechanism become, inside the community's template and
   outside it? In particular, could it ever reach my standing bar of 0.1, or the 0.01 stepping-stone tier?
3. **Can any change to the centre, the decoder or the helpers lift the record?** Even a few per cent would beat the
   community's figure, which is my condition for publishing a record.

**Short answers.**

1. **One mechanism, "crossing slack" through an index-code centre.** It is the same in all three circuits; they
   differ in the port design (Alon's triples for OpenAI, E8 for the community and me), in how junk is cleaned and in
   how the paid work is priced (section 4.1).
2. **Net of the centre it is tiny, and nothing known scales it.** 0.36 of a move per data array at block 73. Inside
   the community's model the ceilings are +0.47 % for T74's own word, 2.88e-3 with every piece as large as the
   window allows, and 3.98e-2 under the most generous relaxation; these are ceilings of a model, not proofs. A claim
   that the whole class is capped near 3e-3 was refuted: under its own inputs the ceiling is about 5.6e-3, and the
   exclusion of 0.1 rests on an argued wall. At small blocks part of that wall now has a written argument that a
   skeptic checked step by step (sections 4.2 and 4.5).
3. **No lever moved the record.** Every lever priced in these rounds is closed or unknown. Open are centres that
   change the junk (costs down to 53 are not excluded), wider label tables, centres with several layers, a decoder for
   one point at block 22, a helper circuit that fits into five stages, and the 112 "M copies" with the cut values
   allowed to change (section 6). The record did move later, by a different route (the recursion lever in
   [../record-track/](../record-track/)), which leaves the centre and the free directions alone.

Earlier, shorter summaries of the copy-wall and helper rounds of 13 October are in
[notes/open-directions-3.md](../../notes/open-directions-3.md), sections 5 to 7, written that morning; this section
shows the work behind them and says where later rounds changed them.

## 2. Words used here

- **Array, block m, paid move, child.** A unit transforms W arrays, each along m binary directions (its *block*). A
  paid move is one recursive call on r of an array's directions, a *child* of rank r. Larger children cost less per
  direction, so the same free directions save more when the paid work comes in larger pieces (*grouping*).
- **Saving.** The exponent gain when the unit is used recursively, quoted in the 1.68e-3 style. To first order,
  saving ≈ (net free directions per array ÷ m) ÷ (the average of ln(m/r) over the paid children).
- **Frame.** The directions an array has already been transformed in: a subspace of the m binary coordinates. Paid
  moves enlarge it; a free addition may join only arrays whose frames are equal.
- **Template, replica, call, stage.** The community's template runs one fixed small program, the *word*, in 9
  *calls* (stages) per *replica*; 4 pairs of banks of 112 data arrays each take turns, and the idle ones keep
  climbing. Calls alternate between forward (y += x) and backward (x −= y). T74 uses 74 replicas; my words use 8.
- **Crossing slack.** A data array starts on a line of the label space, and its data ends on another array whose
  frame is the complement of that line. That pair of frames still delivers the whole transform, but no array pays
  for that one direction. The free directions per replica, net of what the centre takes back, are written *D*.
- **Port, label, centre.** Each data family is a *port* with a binary *label*. The *centre* is the few shared arrays
  that carry data from one bank to the other in each call; each must drop part of its frame and climb again, and
  that is its *charge*, in levels per call. Today's centre has 8 sums of 8 levels each: 64 per call.
- **Junk and helpers.** The centre also delivers small unwanted multiples (±1/2 or ±1/4) of neighbouring ports'
  data. *Helpers* (OpenAI's "carriers") cancel it along a second route. *Hosting* parks a helper's leftover climb
  inside an idle data array's climb, so the helper adds no array.
- **Decoder, total, copy.** The centre's recipe is its *decoder*: a list of *totals* (sums of sources). A total that
  must reach targets which cannot cancel its other terms is *copied* down to a lower frame; a copy of rank r costs r
  moves per call. The copy cost per call is *cst*; the *copy wall* is today's cst: 64 on the complex side, 63 on the
  bit side.
- **Star decoder, natural decoder.** In a *star* decoder each total is dropped by one copy all the way to frame 0,
  and costs the dimension spanned by the labels it reads. A *natural* decoder's totals are linear in the root
  vectors; the community's decoder and every change of basis of it are natural.
- **Keeps today's junk.** What a decoder delivers can be written as a 112 × 112 table. If it has rank 8, a theorem of
  the community says it is today's table up to rescaling, so today's helpers still clean it.
- **The frame game.** The rules every copy floor here assumes: additions only between registers at the same frame,
  frames only climb, content goes down only through a paid copy, frames are subspaces. The Lean engine of this
  repository allows additions across frames, so the frame game is narrower than the engine.
- **P112, P111, P120.** The 112 ports of the community's complex side (E8's 120 lines minus one octad of 8 mutually
  orthogonal lines), the 111 ports of its bit side, and all 120 lines of E8.
- **b, c, κ.** The community's complex-side saving b, its bit-side saving c, and its headline figure
  κ = min(b, c)/(1 + min(b, c)). On 13 October c = 1.7334e-3 sat above b, so only complex-side gains moved κ, and
  only up to about 1.7304e-3 unless the bit side improved too.
- **Price, what-if, ceiling.** A *price* is a saving computed by a program for a construction that exists; a
  *what-if* is a price that relies on a step nobody has built; a *ceiling* is the most a model allows. None of the
  three is a theorem.
- **Skeptic.** An independent AI agent that writes its own code and tries to refute a claim.

**Status words.** Every finding carries one: *proved in Lean*; *checked by computer* (an exact or exhaustive
computation, with whether anyone repeated it); *argued by hand, not checked* (a written proof, with whether a
skeptic re-derived it; no person and no proof checker has checked any of them); *refuted* (and by what); *open*.
Prices and ceilings are said to be computed by a program; they are not theorems. Most results here have one
skeptic, not the two the rounds' own rule asked for; I say so where it applies. When two notes disagree, the later
skeptic verdict wins, and I say where that happened.

## 3. What I tried

The work came in six rounds, each with builders and skeptics; the skeptics wrote their own code.

- **The copy wall (13 October).** The agents listed 17 assumptions behind the community's 64 and 63 copies per
  call, tested those they could, developed five ideas (relays, bit-side decoders off the community's pattern,
  isotropic cuts, small-span totals, reusing copies across stages), and sent every claimed proof to two skeptics,
  one re-running everything on a clean copy and one checking the model against the community's rules. A follow-up
  that day replayed stage sharing exactly (riders and a "permutation door"), each with a skeptic.
- **The helpers (13 October).** A census of the 504 complex-side helpers per call of the community's U8 word, an
  exact toy of its visit rules, seven ideas tested (narrower ends, half reads, bit exits and entrances, fewer M
  copies, one helper serving two planes, shared pair rows, retiring bit carriers that only pass one source along),
  each priced exactly and checked by two skeptics where a claim was made; then the same on the community's newer k8
  word.
- **Symmetry and other fields (13 October).** Six tracks: three borrowed tools from other fields (lattices invariant
  under the Clifford group, such as Barnes-Wall and Leech; designs with a prescribed symmetry, by the Kramer-Mesner
  method; superoptimisation of the community's best word), two confirmed earlier closures, and one asked whether a
  total built in one stage can be carried into the next.
- **The anatomy of savers, round 1 (15 October).** Three agents each took apart one source (OpenAI's unit, the
  community's T74, my Lean-checked words); one drew them together and proposed scaled versions; one built the best
  proposal exactly; two skeptics checked the synthesis and the build.
- **The anatomy, rounds 2 and 3 (15 October).** The two leads round 1 left (a word that crosses a second line only
  once; a small-block network), then networks at medium blocks 13 to 25 and a second line that differs from family
  to family. Each had a skeptic with its own code.
- **The centre (15 October).** One agent mapped what the community had already shown about centres, one searched
  the open shapes for a cheaper centre and then wrote an argument that the main one has none, and one skeptic
  checked both.

## 4. What came out

### 4.1 One mechanism behind every known saver

**Finding.** OpenAI's unit, T74 and U8X0HQ all save the same way: one free direction per data array, delivered
through a centre that works like an index code (it broadcasts a few sums; each target recovers its own value by
subtracting junk that helpers deliver along a second route, and that route is cheap only between ports whose
labels are orthogonal). My own words carry no saving object of their own: their saving is the community's network
and helper circuit, and what my Lean route added was removing overhead it had charged. *Status: argued by hand,
not checked. The synthesis skeptic confirmed the account of the mechanism and the design figures, and the skeptic
of a later round (the frame-packing section, [../frame-packing-two-arrays/](../frame-packing-two-arrays/)) agreed
that no known unit has a second way to save. The OpenAI and community dissections had no skeptic; mine was checked
in part.*

The three dissections, with the status of each figure:

- **OpenAI's unit.** About 1.87e21 arrays at block 10^6, saving 2.65e-13 (2.11e-13 after padding the array count to a
  power of two). Its paid move equals the Hadamard matrix up to free phases, so it is a Walsh-Hadamard unit in my
  own game. Two banks, 8.46e15 arrays in all, swap in three stages through 101 centre arrays built on Alon's 1998 design
  of triples on 100 points; the centre takes back 19 % of the gift; one private carrier per ordered pair cleans the
  junk, and carriers are 99.9995 % of all arrays, diluting the saving about 220,000-fold. The exact rule, paid moves
  = W·m − 2 × (starting-frame directions of data that changes arrays) + 2 × (directions the centre drops), was
  checked by brute force at small sizes. Taking 25 points instead of 100 would have given 4.18e-10, about 1,600 times
  more. *Status: checked by computer, one program each, not repeated by a skeptic.*
- **The community's T74.** One word of about 12,500 additions and halvings at block 73, 112 ports on E8 built from
  the Hamming code, 8 centres of rank 8, 480 helper roles, 896 data families per replica, every helper hosted. Per
  call 99.56 free directions minus the centre's 64 leaves 35.56 out of 7,267.6, that is 0.36 of a move per array
  (0.49 %). With every child at rank 1 the same D would give 1.143e-3; T74's average child rank of 3.94 lifts it to
  1.6797e-3, so a third of the saving is grouping. In the community's own checkpoint files the net free fraction
  (free directions net of the centre, as a share of all directions) is exactly 4.8924e-3 at every checkpoint from
  Jeremy Kothe's atm3 (the circuit Shea checked) to T74, so the whole +4.7 % climb since then is regrouping (average
  child rank 3.44 to 3.95). *Status: checked by
  computer from the community's checkpoint profiles, one program, not repeated by a skeptic. The round named this
  as the figure most worth checking. I re-ran it for this revision: `lineage.py` reproduces its recorded table.*
- **My Lean-checked words.** Every data array pays 72 of its 73 directions, so 896 are free per replica; the
  centres take back 576 (64 per call, 9 calls); helpers net nothing. So D = 320 per replica, 5/14 of a move per data
  array, the same on all 11 words of my ladder. Of my factor 1.318 over Shea's start, 1.032 is the community's better
  word and 1.277 is overhead that my Lean route had charged and my layer removed. *Status: checked by computer; the
  build skeptic reproduced D = 320 and the fold root of U8X0HQ with its own code. Re-run for this revision:
  `identity.py`, `anatomy_ours.py`, `knobs.py` and `stages.py` reproduce their recorded outputs.*

Why the small near misses of the moonshot rounds cannot borrow this mechanism (*argued by hand, not checked*): each
target receives two singular routes whose sum is the identity, so neither has to be invertible; every junk term
has a route that avoids the target; every gate joins arrays at one frame. The price is many arrays per saved
direction: with 2 arrays a centre alone costs at least 2m, and E7, with 63 ports, is the smallest root system whose
copy ratio is below 1.

### 4.2 What caps it

**Ceilings inside the community's model** (*prices computed by a program, ceilings of a model, not proofs. The first
three come from the dissection of T74, which had no skeptic; the build's skeptic reproduced the exact-fill price
1.6647327e-3 with its own code; the what-if with free centres was not repeated*):

| ceiling | saving | over T74 |
|---|---:|---:|
| T74's own word, best packing | 1.6876e-3 | +0.47 % |
| every in-call and hosted child as large as the 9-direction window | 2.88e-3 | ×1.7 |
| centre free and every child at the largest rank allowed | 3.98e-2 | ×24 |
| my U8 line, exact fill in the community's network | 1.6647327e-3 | ×0.991 |
| my U8 word with the centres free (what-if) | 4.67e-3 | ×2.8 |

With one free direction per family, 0.1 needs a block of 28 or less, even with the centre free and maximal
grouping; every saving design in the template needs a block of 25 or more, and at those blocks it nets far less
than one free direction. Only one knob scales in principle: more free directions per array. The anatomy says what
a new mechanism would have to deliver: about 0.5 to 1.5 free directions per array at block 8, or 2 to 4 at block
16. Today's savers give 0.36 at block 73.

**Refuted: "the whole crossing class is capped near 3e-3 and cannot reach 0.1".** The synthesis claimed it; its
skeptic showed that under the synthesis's own inputs (112 ports, 3 stages, block 25) the ceiling is about 5.6e-3,
that 3e-3 appears only if T74's own grouping is assumed, and that the exclusion of 0.1 rests on a "stage-window
wall" that was argued, not proved. Round 2 then gave a written argument for part of that wall, which a skeptic
checked step by step (section 4.5).

### 4.3 Which port design: E8 is the best centre among the root systems

**Finding.** Among the root systems, E8 has the cheapest centre per port (*checked by computer; the synthesis
skeptic rebuilt the table from scratch*):

| design | ports | centres | junk per port | window | copy ratio (centre cost per port) |
|---|---:|---:|---:|---:|---:|
| E8 | 120 | 8 | 56 | 9 | 0.600 |
| E8 minus one octad (the community's 112 ports) | 112 | 8 | 52 | 9 | 0.643 |
| E7 | 63 | 7 | 32 | 7 | 0.778 |
| E6 | 36 | 6 | 20 | 7 | 1.167 |
| D_n, A_n | | | | | 1.0 or more |

The complex (Gaussian) form of E8 halves the centres but needs a window of 12 or more, so it is worse; Alon's
triples tie E8 only with 13 to 29 times its junk. That E8 is the best among *all* two-angle line systems rests on a
classification theorem (Cameron, Goethals, Seidel and Shult) cited from memory: *argued, not checked*.

**Larger lattices fail their screens** (*checked by computer*):

- **Barnes-Wall (BW16).** 1,304 junk terms per port against a mark of 120; a clique of mutually junk lines forces a
  window of at least 30, and the build skeptic found one of 53. The symmetry round reached the same verdict
  independently on both sides: the BW16 complex ceiling at nine stages is 1.2863e-3, below the community's b; a real
  BW16 bit decoder built from 18 E8 blocks prices at 9.54e-4 with helpers and 1.38e-3 without, against 1.86e-3 for
  E8 in the same model (no skeptic report came back for that round). BW32 reaches at most 0.13 of E8's ceiling.
- **The full Leech set.** A 23-line clique excludes it; Leech reaches at most 0.06 of E8's ceiling. Leech
  sub-designs that leave those lines out were not screened (open, section 6).
- **Coxeter-Todd.** Real form: 378 lines, 242 junk per port, a clique of 45, copy ratio at least 1.43. Complex form:
  126 lines, window at least 15, copy ratio at least 0.71, and its phases are cube roots of unity, which are not
  free phases in the game. The round-2 skeptic rebuilt it from a different construction and every figure matched.

A side lesson of the symmetry round: the copy ratio alone misleads, because the window matters as much. The clique
screen (the largest set of pairwise non-orthogonal lines) rules out a proposed unit in seconds, before anything is
built.

### 4.4 A second free direction per array

This was the one knob that scales, so the rounds tried it three ways.

**Two-sided slack on my U8 word.** Give every Y array a second free direction s. On paper that is +896 moves per
replica. But the data must end in X arrays, so in each of the four backward calls per replica all Y information
leaves s: two "dirt gates" of exact rank 112 make 224 crossings per backward call, 896 in all, exactly the gain.
The skeptic reproduced the calibration, both ranks and the four backward calls, and found that the forward calls
pay too (at least 112 more frame changes each), so at this gate order two-sided slack loses outright. *Status:
checked by computer (exact crossing counts), one skeptic, who made the kill stronger. Re-run for this revision:
`calib.py`, `twosided_local.py`, `vcalib.py` reproduce their outputs.*

**Any word with one shared extra line: the crossing lemma.** No local word with one shared extra direction s can
beat the same word with s deleted, whatever its gate order. A free addition never joins arrays on both sides of s,
and each crossing carries one direction of data, so per replica 896 pieces must cross in and out: a tie at best.
The trivial lift reaches the tie (D = 320 at block 74 instead of 73), and a tie at a larger block prices lower:
1.6246236e-3 with the crossings priced as separate pieces, 1.6537530e-3 if they are absorbed into other work, both
below the one-sided ceiling 1.6647327e-3. *Status: argued by hand, not checked; checked exhaustively by computer on
small cases; one skeptic re-derived it and corrected the author's floor (see section 5). The deletion step is checked
by computer on 120,000 random steps;
re-run for this revision, `vproj.py` reproduces its output.* This also answers an open question the community had
listed (why only one unpaid direction per family): the endpoint contract does not forbid a second one; the way
back does.

**A second line per family (the mirrored line).** Give each family b its own line s + g_b, with g_b the port's
complement in the window. The lines exist: their binary Gram matrix equals the ports', so a second copy of the
community's code can serve as the second centre, and there are 69,616 "twisted" variants, all symmetric involutions
(*the lines and their Gram matrix are checked by computer, and the skeptic reproduced every fact; the count of twists
and their being involutions are argued by hand, with the key step reproduced by the skeptic*). The crossing counts stay
below the gain (9 or 564 per replica depending on how the gluing is read; *one program, not re-run by the skeptic,
because the verdict does not depend on it*). If only a second centre at 64 per call were paid, the price would be
3.2824847e-3, but to beat the community everything else must cost at most 34 levels per call. It fails on cleaning: each
x now ends on its own hyperplane, so a scratch array can visit at most 4 (natural mirror), 8 (any twist, checked
exhaustively) or 45 (any lines) of a total's 56 x's without a paid drop, and no single pass cleans a total. The best
cleaning scheme known, found and checked exactly by the skeptic, costs 120 per call, about 3.5 times the room. *Status:
dead at first order, argued by hand, not checked; no lower bound above 34 is proved for every scheme.*

### 4.5 Smaller blocks: a wall at small blocks, prices at medium blocks

**Small-block networks cannot reuse window directions** (*argued by hand, not checked; one skeptic checked each
step with its own code and found four limits of scope*). In the model the proposal was posed in (fixed addresses,
frames that only grow except at centres, one frame per gate, and each side-information copy erased against its own
source), each call uses up σ fresh directions in the arrays it touches, σ being one less than the window; in any
exchange some bank is a source twice (argued for any number of banks, and checked exhaustively on small cases), and
those two calls need orthogonal windows. So the block is at least 1 + 2σ: 13 for E7, 17 for E8 and the 112-port
design. A new floor on the centre's charge (for centres that gather, drop once and scatter) shows that D4 and A3 save
nothing. Generous ceilings with one free direction per array stay below 0.04 (0.034 for E7 at block 13 and for the
112 ports at block 17). The four limits: one-layer centres only; σ is computed for each design's natural code; the
floor needs the bank that sends twice to use shared centres both times; and the ceilings put each array's paid work
in pieces of at most σ directions, which is argued for about half the arrays. If the work came in one piece the same
model would give 0.18 to 0.33, so this is the assumption the "below 0.1" statement rests on most. The raw-array
checks (frame rules at blocks 1 to 5, a lone pair costs exactly 2m, two pairs sharing a centre end exactly right)
were *checked by computer by both agents with independent code*.

**Medium blocks 13 to 25: at most 7.39e-4 with any word anyone has.** The community's template shrunk to one pair of
banks (three calls) at block 25 gives 7.3920e-4 with T74's word and 7.3057e-4 with the word of my Lean-checked
line, both with helpers on ideally packed banks: ×0.44 of T74 and 13.5 times short of 0.01. With junk cleaning free
my word would give 2.3338e-3; the most generous ceiling at block 25 is 8.2607e-3. *Status: prices computed by a
program; the skeptic reproduced every price to every digit with its own fold. Re-run for this revision:
`price_p112.py`, `shrink.py`, `price_t74_p1.py`, `items.py`, `ceilings.py` and the skeptic's `vprice.py` reproduce
their outputs.* Shrinking the template only loses (1 pair 7.31e-4, 2 pairs 1.47e-3, 3 pairs 1.65e-3, 4 pairs
1.6647e-3, 5 pairs 1.65e-3); nine calls is the best layout for the U8 word.

**Centre floors at medium blocks** (for centres that gather, drop once and scatter; *argued by hand, not checked,
with the numbers computed by a program*): D5, D6, A4 and A5 cost exactly 0.8 per port per call, which equals the
two-pair break-even, so they never save with one or two pairs; E7 at least 35 per call (every decoder of the known
kind at least 42, the break-even; decoders on an orthonormal basis at least 48; the known one costs 49); the 112
ports at least 52.36; E8 at least 54.45. The skeptic rebuilt the floors of D5, D6, A4, A5 and E7 with its own code,
without building labels, and found them stronger than stated. For the 112 ports and E8 it upheld the argument but
did not recompute the numbers; the centre round later recomputed 52.36 exactly with two programs, and 54.45 has not
been recomputed. Free meetings, where side information rides on copies that never climb, can clean at most 58 % of
the 112-port junk and 44 % of E7's (*argued by hand, not checked; one skeptic upheld it in that scope*).

### 4.6 A cheaper centre for the 112 ports

**Finding: no centre that keeps today's junk costs less than 64 per call.** Every single-layer decoder of the 112
ports whose delivered table has rank 8 costs at least 64 per call under the charge that both my checker and the
community's shipped format use (each sum pays the dimension of the frame it is gathered from), for any number of
terms built anywhere. *Status: argued by hand, not checked; its counting steps and its last cases are checked by
computer. It is a written proof whose last cases are closed by exhaustive search; one skeptic re-derived the
rewriting step by hand, recomputed every computed input with its own code and checked the list of slices by a
separate argument, and found no gap (it had predicted, at 60 %, that it would).* For my own line the charge is
confirmed in my checker; for the community's own centre it rests on one sentence of the community's text.

How the proof goes. Using the community's theorem that a rank-8 table is today's up to rescaling (burntout76's and
Jose Manuel's, stated for real numbers and taken as given), any such decoder can be rewritten
at the same cost as a dictionary: a list of directions in 8 dimensions, each port's root a combination of the
directions it uses, each direction paying the dimension spanned by the labels of the ports that use it. Two
necessary conditions follow (the ports outside any k-dimensional slice of the roots must be reached by the other
directions; at most one port can rely on a single direction). Exact tables feed a count: the most ports a label
space of each dimension holds (1, 2, 4, 8, 12, 20, 36, 64, 112) and how few ports lie outside a 7-dimensional slice
(53, then 56, then 70). The count rules out every cost of 63 or less except a few cases with 9 to 12 directions, and
an exhaustive search shows that the covers those cases need do not exist. Calibration: the search found the
community's 64 decoder and the coordinate decoder at 72 from scratch, refused a cut-down decoder at 48 that fails
exactness in 2,333 places, and every one of the 13 exact decoders it built (costs 64 to 140) passes the counting
bounds at its own cost. Re-run for this revision: `calib.py`, `negctl.py`, `nd.py`, `coflat.py`, `coflats_all.py`,
`profiles.py`, `xcheck.py`, `m9.py`, `cover3.py`, `cover3b.py` and `soundness.py` reproduce their recorded outputs;
[programs/README.md](programs/README.md) says what happened to the last step, `h2k4sym.py`.

Before this round the community had 64 only for 8 terms or fewer, for terms inside the table's column space, for 9
terms with exactly one outside, and for terms built from the root code. Its results hold under a cheaper charge
(the pairing rank); mine holds under the shipped charge only. A centre at 60 would have put my line near 1.85e-3.
At +2.8 % per unit of centre cost, this was the only shape where a cheaper centre could be priced straight away,
because today's helpers would still clean its junk; it is closed.

**Two counting floors** (*checked by computer, each by two independent programs*): every single-layer decoder costs at
least 52.36 per call (a count over the 1,260 orthogonal frames); the window track's intermediate floor of 60 for
centres that keep today's junk (59.0733 from a weighted count over 9,300 bases, checked against all 8,283,457
subspaces) is superseded by the theorem above. Under the pairing charge (projected copies, which the community's
checker forbids) no counting floor above 16.9 exists (the skeptic's own, cruder search showed that the count can give
at most 20.1, the same conclusion). Re-run for this revision: `onesided.c` reproduces the 52.36 count.

### 4.7 The copy wall (13 October)

The round summary is in [notes/open-directions-3.md](../../notes/open-directions-3.md), section 5; here is the
evidence and what later rounds changed.

- **Relays never beat the static floor.** A relay is a target reading another target's register. A target can read
  only totals that were already copied and sources it already knows, so the same totals fit statically at the same
  cost, in any order, any number of rounds, scalar or vector, over any field. *Status: argued by hand, not checked,
  in the frame game; upheld by two skeptics; and checked by computer with a SAT solver (rank 8 is impossible with
  relays on P120 and P111; a skeptic re-ran it on the community's true bit puncture and got the same answers).* The
  proof uses the community's frame-sandwich rule, which the community itself marks as argued.
- **Bit side, Theorem A.** Every F_2 star decoder whose read span contains the label space and has dimension 11 or
  less costs at least 63 (at least 64 at dimension 10), and every decoder with 9 or fewer totals has such a span.
  *Status: checked by computer (exact, on both punctures); upheld by both skeptics within that scope (one skeptic's
  own encoding confirmed the true puncture at 9 or fewer totals; the other's did not finish in 20 minutes).* Not
  covered: pairing prices, 10-total decoders whose reads miss the label space (they exist; the best of 42 sampled
  costs 65), 12 or more dimensions, payloads over F_3.
- **Complex side.** Any decoder on P112 costs at least 53 (a one-line rank argument, that every octad needs 8
  different totals, followed by an exact count of the octads each frame meets over all subspaces), and natural
  decoders cost at least 64, only the community's reaching it (exhaustive over subspaces). On the full E8 unit
  natural decoders need 72 (exhaustive). *Status: checked by computer, one program each, not repeated by a skeptic.*
  An older bound of at least 55 copies for every decoder on P120 (10 October), by the same kind of octad count, is
  *argued by hand, not checked, with its counts computed by a program; no skeptic has read it*. So the open complex
  region is 53 to 63: decoders with at least 9 totals that are not natural.
- **The octad Hall test**, a cheap necessary condition (each octad's 8 ports must be matched to 8 different totals).
  In the copy-wall round, shapes at cost 64 and 70 passed it easily and no cost-63 shape was found to pass it (the
  smallest shortfalls were 10 and 4; the SAT runs were undecided), so the round planned to decide it for the three
  cost-63 shapes as its next step. The symmetry round then showed that the test is too weak to lean on: a cost-63
  family of 9 frames of dimension 7 passes it on all 1,260 octads, yet fails on 21 other maximal cliques and is not
  a decoder; and counting over all 1,996 maximal cliques gives only 52. *Checked by computer, no skeptic report.* So
  any future screen of this region should use all maximal cliques or a matroid-partition test, not the octad test
  alone.
- **Stage sharing ("can stage i's copies serve stage i + 1?").** The count of frame-sharing pairs was corrected (ports
  S and T share a frame iff their labels are orthogonal: 6,470 of 12,210 pairs), and one single-carrier escape was
  closed by linear algebra (*upheld by a skeptic*). An exact two-stage replay of "riders" (registers that carry
  several ports' content on climbs paid anyway) in a toy at block 3 gives an exact lower bound of 18 paid moves for
  every word; riders with their cleanup cost 18 or 21, and without cleanup 15 with six wrong registers (*checked by
  computer; the skeptic agreed with independent code that shares nothing with the author's*). The toy is too small to
  show a gain of any kind (an earlier lemma rules out every saving at block 3 with these moves), so this closes riders
  inside the toy only; at the community's scale they are not decided. The "permutation door" fails at the stage
  boundary on all 32 tiny pairs and all 6,470 pairs of the community's real bit geometry. The author's cheapest repair
  cost 8 moves more than the standard copies on the tiny unit; the skeptic found a cheaper repair that costs exactly
  as much as the standard copies (168 moves), still above the plain word's 160, so it saves nothing (*checked by
  computer, upheld*). The symmetry round's stage-equivariant track then argued in writing that in the community's
  construction no total can cross any of the five stage boundaries: a total's 56 ports continue into 56 different
  invocations, and paying to re-route costs at least 220 moves per total per boundary against 8 for the star copy
  (*argued by hand and checked against the community's own boundary files, all 120 ports and 1,080 stage translations;
  the small case h = 3, m = 5 checked exhaustively; no skeptic report came back*). Re-run for this revision, the same
  program also completed the whole cover at h = 3, m = 7 (1,451,520 invocations) and the case h = 4, m = 7, with no
  counter-case; the recorded run had covered a sample of 2,001 invocations of the first.
- **Monotonicity.** Adding one port raises the copy cost by at most 1, so 72 for the full E8 unit would give at
  least v − 48 for every punctured unit, which is exactly where every known decoder sits. *Argued by hand, not
  checked; no skeptic has read it.* This is the 72-copy question of section 6.

### 4.8 The helpers (13 October)

The round summary is in [notes/open-directions-3.md](../../notes/open-directions-3.md), section 6.

- **Census.** On the community's U8 word (in Jonathan Capalbo's packing) the 504 complex-side helpers per call are 392
  "plane visitors" and 112 "M copies". The plane visitors meet Dominik Scholz's packet-plane bound for this packet
  design exactly, and the agents' integer programme reproduces the community's 36 helpers per block, so 504 is
  optimal inside the community's visit model. The M copies are the only registers that reach each target's final
  frame; nobody has a lower bound for them. On the k8 word the census is 480 = 368 + 112. *Status: checked by
  computer, reproduced by a skeptic.*
- **Nothing removed on the complex side.** Narrower ends: an exact no for all 378 full-width helpers (*upheld, and
  stronger than stated: the floor holds even with the entrance frames emptied*). Two pairs can instead end one rank
  lower, which makes them cheaper without removing them. On U8 one pair was Alexander Schwartzman's and the other was
  new; on the k8 word, by then the community's newest, both had already been claimed by others (section 5). Half
  reads: 0 of 168 second copies (U8 already uses them at every target; *no skeptic ran*). One helper serving two
  planes: 0 pairs, exact obstruction at the first step (*no skeptic ran*). Bit exits are impossible (*exact,
  upheld*).
- **One bit-side carrier removed.** Carrier 470 retires (439 to 438 helpers) in a word that passes the community's own
  bit checker, worth +0.012 % to +0.119 % on c; κ does not move, because the complex side sets it. *Checked by
  computer with the community's checker; no skeptic ran, and the hosting step was not re-run, so the price is a
  range.*
- **The M copies.** A claimed floor of 8 per group was *refuted as stated*: every target already holds all 64 packet
  sources after the cut, so the M copies cancel content rather than deliver it, and the floor holds only with the
  centre values and the cut held fixed. The Walsh-packet idea for them breaks the frame rule in all 14 groups.
  Removing all 112 would be worth +4.9 % to +12.7 % of b (what-if). This is the largest open helper lever.

### 4.9 Symmetry, superoptimisation and the community's best word (13 October)

Nothing moved: the community's figure stayed at κ = 1.6722e-3 with b = 1.6750e-3 and c = 1.7334e-3. Five tracks
ended in a written argument or an exact search that rules out their idea within a stated scope. *No skeptic report
came back on any track of this round; every result is negative, so the risk is a missed gain, not a false claim.*

- **Symmetric decoders.** About 8.7 billion frame families invariant under a symmetry of order 7 or 8 were tested
  at cost 63 or less with exact necessary tests; none survives, while 21 cost-64 controls do. Orders 4 and 3 were
  not finished. In each Hamming block the 12 face planes are forced, and the exact best per block is 28, the
  community's. *Checked by computer.*
- **Rewriting the k8 word.** The word cannot be improved by local rewrites: its frames at its own gate order are
  exactly optimal (an integer programme returns the word in 0.4 s), and with commuting gates reordered, integer
  programmes show it optimal for dependency gaps up to 4 (98 % of the gates). At gap 5 the solver stalled after about
  an hour with at most 1.66 cost units (about 1.3e-7 of κ) not excluded. Four rewritten words passed the community's
  audit and priced the same or worse. *Checked by computer.*
- **Confirmations.** Re-graining the idle remainder is not allowed in the community's model (its one legal variant,
  reversing one entrance promotion, is seslly's, already open in the community); Walsh packets for the M copies break
  the frame rule (the pipeline re-run on the k8 word gave 35 byte-identical files).

## 5. What failed or was refuted

Failures are as useful as successes, so here they are in one place. In each case the later verdict is the one used
in this section.

- **The anatomy synthesis's cap.** "The whole crossing class is capped near 3e-3 and no construction in it reaches
  0.1": refuted by its skeptic; under its own inputs the ceiling is about 5.6e-3, and the exclusion of 0.1 rested on
  an argued wall.
- **Two-sided slack (proposal 1).** Predicted ×2.03 (3.4e-3); killed by the exact build (the backward calls pay the
  gain back) and the kill made stronger by its skeptic (the forward calls pay too). The skeptic also faulted the
  proposal's test: its kill mark was set at an extra charge of 128 per call while its own model breaks even at 99.9,
  and the build's measured 99.6 fell into that gap.
- **The single-crossing word.** The round-1 ceiling of 3.2e-3 for a word that crosses s once was withdrawn when the
  crossing lemma showed no such word beats its one-sided version. The round-2 author's floor "D ≤ 208, price at most
  1.05e-3" assumed that no X array holds s between calls; the skeptic corrected it to "a tie at best".
- **The mirrored line's cleaning floor.** "Cleaning costs at least 272 per call" was the cost of one greedy scheme,
  not a floor; the skeptic built a scheme at 120. The lead stays closed at first order, but only by argument.
- **Medium blocks.** "The ×0.44 comes from the missing idle levels" was refuted: it is mostly the centre's share of
  each call (×0.52); hosting accounts for ×0.84. "Five template points, one saves; medium blocks are closed" was
  refuted: the skeptic found points at blocks 16 and 22, and the block-22 point is open (section 6).
- **Copy-wall "proofs".** The isotropic-cut proof and the stage hand-over proof were both refuted as proofs by both
  of their skeptics (the first left out star totals of width 7 or less, which is the open region, and forbade
  off-chain steps instead of pricing them; the second had an untested escape, needed a new rule of its own and
  hard-coded several "computed" steps). What survives of the first is an exact counting lemma for isotropic line
  cuts with 9 or more totals. The escape the second left open was later closed by the stage-equivariant argument
  (section 4.7).
  The stage-sharing toy proof of 13 October ("no saving on a restricted toy") was likewise refuted as a proof by its
  skeptic; the exact rider replay then settled the toy.
- **Shared pair rows: "56 fewer helpers".** The visit model's optimum of 448 helpers was real, but both skeptics
  refuted the gain: the model charged U8 8 helpers per block for work U8 does for free, left out the M copies, and
  its recipe breaks the frame rule (about 336 extra face copies). The same idea on the k8 word ("32 to 40 fewer")
  was refuted by its skeptic: like for like, the model is worse (448 against the word's 368 in-block helpers).
- **The rank-7 helper pair on the k8 word.** The narrower-ends track found that helpers 580 and 581 can end one rank
  lower and said nobody had claimed them; a build stacked them with Kothe's pair 574 and 575, on Kothe's owner map for
  that pair repacked by my agents, and passed the community's whole admission pipeline at κ = 1.673020532e-3 (its
  skeptic reproduced that figure). The claim of novelty was refuted by the skeptic: Thomas DiFiore had already claimed
  580 and 581 in a pull request opened earlier that day, with six more, as a strict superset with a higher admitted
  figure. Nothing was posted.
- **The M-copy floor of 8 per group.** Refuted as stated (see section 4.8).
- **A method lesson.** The octad Hall test, on which the copy-wall round relied, turned out too weak (section 4.7).

## 6. What is still open, with a first experiment for each

None of these is a unit. Most would lift the community's figure and mine, by a few per cent up to about ×2 at the
very most, if they worked out; items 8 to 11 are closures that would tell people where not to look. Only two come
anywhere near the 0.01 tier, and only under assumptions nobody has met: the block-22 point (item 4, and then only
with junk cleaning free) and leaving the frame game (item 13, and then only with the helpers gone as well).

1. **Centres that change the junk.** For the 112 ports, decoders whose delivered table has rank 9 or more: the
   lowest cost not excluded is 53. A hit would also need new helpers, priced again; at 53 my line would be near
   2.18e-3 if the new junk cleaned as cheaply as today's, which nobody knows. *First experiment:* for the community's
   fixed nine-term completion, compute how few ports lie outside a 7-dimensional slice of its own vector
   configuration, then rerun `centre-floor-search/profiles.py`, `m9.py` and `cover3*.py` with those numbers. Seconds
   to minutes on a laptop.
2. **Wider label tables (rank 10 and up)** outside the one case the community tried (tie prices 59 at rank 10 and 55
   at rank 11). *First experiment:* extend `centre-floor-window/onesided.c` (the 52.36 count) to label tables of
   rank 10 and compare with the tie prices.
3. **Centres with more than one layer.** No floor exists at all, but the shipped verifier restricts the extra
   descents this needs, so a hit means a format change. *First experiment:* price the smallest two-layer centre on
   E7 (63 ports) in the medium-block model, where exhaustive search is feasible, before touching E8.
4. **The block-22 point.** E8's 56 lines at 60 degrees to one root: centre floor 28 per call against a break-even of
   37.3; the best decoder known costs 60. A decoder below 37.3 would make block 22 save before junk cleaning (the
   skeptic's generous estimates: 1.6e-2 with free junk cleaning, 2.95e-3 with U8X0's helper cost). A 52-line section
   of the 112 ports at block 22 is also flagged (formula charge 28 against 34.7). *First experiment:* run the
   centre search's dictionary model (`centre-floor-search/dsearch.py`, `sas.py`) on the 56-line section with target
   cost 37, calibrated first on its known decoder at 60.
5. **A helper circuit that fits five stages (block 41).** About 2.9e-3 at first order, but the calibrated
   medium-block model gives at most 3.44e-3 even with free junk cleaning, and the U8 word run in five calls prices at
   1.47e-3, because only 46 % of its helper leftovers fit into idle climbs. *First experiment:* re-price the
   five-stage layout in the calibrated model (`medium-block-networks/shrink.py` with the helper leftovers halved) and
   search for the circuit only if it still clears the community's figure.
6. **The M copies with the cut re-valued.** Removing all 112 is worth +4.9 % to +12.7 % of b (κ takes only about the
   first 3.5 % unless the bit side improves too). *First experiment:* one exact small case set up as an index-coding
   problem: one group of 4 targets, 8 centres and 8 super-packets, chains fixed, centre values and cut coefficients
   free, least number of relays. If one group needs fewer than 8, scale it up and price it exactly on the k8 word.
   Index coding gives the best linear answer (minrank) and has cases where vector schemes beat scalar ones; check
   first that the frame rule fits that framing.
7. **The complex copy region 53 to 63.** Decoders with 9 or more totals that are not natural, reading 6- or
   7-dimensional subspaces. *First experiment:* screen the three cost-63 shapes (9 totals of width 7; widths 6, 7^7,
   8; widths 6, 7^8) with the test on all 1,996 maximal cliques, or a matroid-partition test, instead of the octad
   Hall test, with the first frame fixed up to the octad's symmetry group (order 344,064); then finish the symmetric
   searches of orders 4 and 3 (`symmetry-invariant-decoders/`).
8. **Does the full E8 unit need 72 copies?** Known: 55 ≤ cst(P120) ≤ 72, natural decoders exactly 72. A proof would
   close the community's problem of beating the star floor for the whole puncture family (with monotonicity, also
   unproved). It is a closure, not a gain. *First experiment:* a SAT or integer-programming run for "a decoder of
   cost at most 71 exists on P120", with the Weyl group of E8 (order 696,729,600) broken, overnight on an 8-core
   machine.
9. **Leech sub-designs without the 23-line clique.** *First experiment:* run the clique screen
   (`anatomy-synthesis-centres/clique.py`) on them; seconds. Expect a fail.
10. **Closing the crossing class formally.** Carry out the third-window count (a swap needs 1 + 3σ, outlined, not
    done), prove the σ limit on piece sizes for every array, prove a centre floor near 0.8 for D5, D6, A4 and A5 in
    more than one layer, and redo the window argument with per-array address maps between calls (the community's words
    use them; the proofs exclude them). *First experiment:* the third-window count by exhaustive search on E7 at
    block 19, with `small-block-window-wall/net.py`.
11. **A lower bound on cleaning for the mirrored line** (room 34 per call, best scheme 120). Low priority given the
    margin. *First experiment:* an exact search over cleaning schedules for one total of 56 x's in the visiting-rule
    model (`second-line-mirrored-skeptic/v5perx.py` is the starting point).
12. **Bit decoders with 10 totals whose reads miss the label space** (the best of 42 sampled costs 65). They move only
    c, so they matter only after b has gained more than about 3.5 %. *First experiment:* the remaining light-pair SAT
    probes in `copy-wall-bit-side/`.
13. **Outside the frame game.** Every copy floor here assumes it. Removing the copies outright would need additions
    between registers at different frames (which the Lean engine of this repository allows), frames that are not
    subspaces, or another unit. *First experiment:* the smallest unit where the wall bites is E7; test in the engine's
    wider rules whether one E7 call with one addition across frames can deliver a target with fewer copies than the
    star floor.

## 7. The programs and how to run them

Everything is Python 3.9 or later unless marked; most programs need only the standard library and run in seconds on
a laptop. Where a program needs numpy or scipy (install them with pip) or python-sat (`pip install python-sat`), its
folder's entry in [programs/README.md](programs/README.md) says so. Each program folder also holds the outputs as
recorded during the rounds; running a program overwrites them. Run each program from its own folder.

- [programs/README.md](programs/README.md) lists every folder, what it shows, its command, its time and memory, and
  whether it was re-run for this revision.
- [data/](data/) holds the shared inputs (my lowering of the community's U8 word, the profiles of my words, the
  community's checkpoint profiles, as numbers) and copies of the key results and certificates;
  [data/hub-inputs/README.md](data/hub-inputs/README.md) lists the community's files that some programs read and that
  are not copied here, with where to get them and their sha256.
- [notes/](notes/) holds the agents' working notes behind each track, as they wrote them, with private paths removed.
  Where a note disagrees with this README, this README wins: it follows the later skeptic verdicts.
- [results.json](results.json) lists every finding with its status, figures and evidence files.

Some programs price constructions with `pricer.py`, my agents' exact re-implementation of the community's moment
engine (it imports none of the community's code). It is shipped once, with the tools of the lever package that the
record-track section ([../record-track/](../record-track/)) describes, and these programs find it through the
environment variable `PRICER_DIR` (the folder that holds `pricer.py`).

**Re-runs for this revision.** On 19 October I re-ran the light programs of the anatomy, centre and copy-wall rounds,
and a few of the helper and symmetry rounds, in a scratch copy of this folder on an Apple M4 Mac mini; the table in
[programs/README.md](programs/README.md) gives each result. Programs that need the community's files, a SAT solver or
more than a few minutes were not re-run for this revision; their recorded outputs are shipped as they were.

## Sources

The rounds behind this section, in my private research archive (the working folders of 13 to 15 October 2026): the
copy-wall round (`checks/wht43/copywall/`, with the parked note OPEN-72.md), the stage-sharing replays in the same
folder, the helper round (`checks/wht43/helpers/`), the symmetry and cross-field round (`checks/wht43/symmetry/`), the
anatomy of savers, rounds 1 to 3 (`checks/wht43/anatomy/`, with the round note NOTE.md), and the centre round
(`checks/wht43/centre/`, with NOTE.md and WINDOW.md). The catalogue entries are R18, R31 to R38, N18, T21, C23 and C24
of the archive's results catalogue. The private provenance file maps every staged file to its source.

Track labels and the folders they became (programs/ unless noted):

| round | track (builder; skeptics) | folder |
|---|---|---|
| anatomy 1 | an-openai | anatomy-openai-unit |
| anatomy 1 | an-hub | anatomy-hub-record-unit |
| anatomy 1 | an-ours | anatomy-our-lean-words |
| anatomy 1 | an-synth; v-an-synth | anatomy-synthesis-centres; anatomy-synthesis-skeptic |
| anatomy 1 | an-build; v-an-build | second-line-two-sided-build; second-line-two-sided-build-skeptic |
| anatomy 2 | an2-word; v-an2-word | second-line-single-crossing-word; second-line-single-crossing-word-skeptic |
| anatomy 2 | an2-smallnet; v-an2-smallnet | small-block-window-wall; small-block-window-wall-skeptic |
| anatomy 3 | cx-medium; v-cx-medium | medium-block-networks; medium-block-networks-skeptic |
| anatomy 3 | cx-mirror; v-cx-mirror | second-line-mirrored; second-line-mirrored-skeptic |
| centre | ctr-search; ctr-window; v-centre | centre-floor-search; centre-floor-window; centre-floor-skeptic |
| copy wall | u-floor; u-ours | copy-wall-floor-checker; copy-wall-our-unit |
| copy wall | d-relay (d-R1); v-R1-recheck, v-R1-model | copy-wall-relays; copy-wall-relays-skeptic-rerun, -skeptic-model |
| copy wall | d-bitspan (d-g-9hyp); v-g-9hyp-recheck, v-g-9hyp-model | copy-wall-bit-side; -skeptic-rerun, -skeptic-model |
| copy wall | d-isocut (d-g-break-1); v-g-break-1-recheck, v-g-break-1-model | copy-wall-isotropic-cuts; -skeptic-rerun, -skeptic-model |
| copy wall | d-R7 | copy-wall-small-span-totals |
| copy wall | d-g-crossstage; v-g-crossstage-recheck, v-g-crossstage-model | copy-wall-stage-hand-over; -skeptic-rerun, -skeptic-model |
| copy wall | g-break, g-codes, g-small, g-recombine, p-pick | copy-wall-idea-probes/ |
| stage sharing | s-fix, s-wide, r-riders, r-swaps; v-s-fix, v-s-wide, v-r-riders, v-r-swaps | stage-sharing-* |
| helpers | u-helpers-hub; u-helpers-ours | helpers-census; helpers-census-ours |
| helpers | d-hpe-multigate; v-hpe-multigate-reprice, -rules | helpers-narrower-ends; -skeptic-reprice, -skeptic-rules |
| helpers | d-halfread; d-bit-exits; v-bit-exits-* | helpers-half-reads; helpers-bit-exits; helpers-bit-exits-skeptic-* |
| helpers | d-mcopy-code; v-mcopy-code-* | helpers-m-copies; helpers-m-copies-skeptic-* |
| helpers | d-aboveplane | helpers-one-helper-two-planes |
| helpers | d-pair-market, pm-build; v-pair-market-* | helpers-shared-pair-rows, -build; helpers-shared-pair-rows-skeptic-* |
| helpers | g-break, g-fields, g-small, g-search, g-recombine | helpers-idea-probes/ |
| helpers | k8-pairs; v-k8-0 | helpers-k8-census; helpers-k8-census-skeptic |
| helpers | k8-twin; v-k8-1 | helpers-k8-rank7-pair; helpers-k8-rank7-pair-skeptic |
| symmetry | t-bw16; t-km; t-egraph | symmetry-lattice-units; symmetry-invariant-decoders; symmetry-word-rewrites |
| symmetry | d-I1-idle; d-I2-walsh; d-I3-stage-equivariant; t-open | symmetry-idle-regrain; symmetry-walsh-packets; symmetry-stage-equivariant; symmetry-open-tests |
| (shared) | the E8 label tables of 10 October | e8-labels |

Whose work this rests on, from the community (CrocSwap/integer-mult-bounds; numbers are its pull requests):

- The template, the words and the circuits: Jeremy Kothe (#588, the circuit Shea checked in Lean, "atm3"; #799; #869 and
  #870, the U8 word; #1018 and #1042, the seven-slot k8 word and its owner packing; #1121; #1839, the T74 line), on
  Chafik Boukhalfa's frame changes (#741) and T21 owner packing (#890) and icekylinx's gauge44 construction (#645);
  Jonathan Capalbo's T37 owner packing of U8 (#942); Rohan Arun's #1082 and #1094 (the reviewed record on 13 October).
  T74 itself: Rohan Arun, Marius-Paul Dumitrean, Jeremy Kothe and LeverRank, on Kothe's k8 line, with the contributors
  the lever package's NOTICE names. My E8 unit entered the community in DaysSky's #352.
- The bit side: DreamingOfClouds' #589 on eumemic's #555.
- Copy and decoder floors: DreamingOfClouds (#717, its decoder-floor note D1 being the community's 64; #2228, #2338,
  #2639, #2724, the decoder-floor pull requests), burntout76 (#1212, rank-8 tables are today's up to rescaling),
  Jose Manuel (jmbarrancoml, #1815), joaqcoo (#1224, the normal form of a decoder), woklerx (#1180), mihir-s-05
  (#845), Josie Elliston (#1787), commdata2338 (#2021, #2385, #2398), Hanhong Zhao (Hope7Happiness, #2321), Sricharan
  Vagvala (theon199, #2382), ywz65 (#2515), Sagar Thakrar (thakrarsagar, #2581), themilkman216 (#968); the co-stars
  discussion (Discussion 1267, in a fork); the community's corrected format-limits notes (star floor, sandwich lemma,
  uniqueness, descents).
- Helpers: Dominik Scholz's packet-plane bound; Alexander Schwartzman (AlexS-scale, #1030, the first rank-7 helper
  pair on U8); Cinna Davis (1337isnot1337, #1090); Thomas DiFiore (#1148, the eight rank-7 ends on k8, which include
  580 and 581); Sagar Thakrar (#1023); ToikY (#988, the U8 terminal-cone census); Hao Yan (#858, the T16 checker);
  seslly (#1111).
- The newest open figures cited: Kothe's #2745 and #2767.

Literature: OpenAI, "An explicit power saving for the exact discrete Fourier transform", "Finite tensor savings and
exact Fourier circuits" and "Integer multiplication below n log n" (preprints, September 2026; the Lean development
is openai/math, family 130); N. Alon (1998), the triple design behind OpenAI's centre; C. H. Bennett (1973) and
H. Buhrman, R. Cleve, M. Koucký, B. Loff and F. Speelman (2014), restoring dirty scratch space; P. J. Cameron,
J.-M. Goethals, J. J. Seidel and E. E. Shult (1976), line systems with angles 60 and 90 degrees (cited from memory,
not re-read); E. S. Barnes and G. E. Wall (1959), and G. Nebe, E. M. Rains and N. J. A. Sloane, the Barnes-Wall
lattices as Clifford-group invariants; E. S. Kramer and D. M. Mesner (1976), designs with a prescribed automorphism
group; Z. Bar-Yossef, Y. Birk, T. S. Jayram and T. Kol, index coding and minrank; M. A. Maddah-Ali and U. Niesen,
coded caching; N. Alon, A. Hassidim, E. Lubetzky, U. Stav and A. Weinstein, and E. Lubetzky and U. Stav, vector and
non-linear index codes that beat scalar linear ones; J. Edmonds, matroid partition.

The community's figures quoted here (T74's 1.6797550e-3, #1094's 1.6727e-3, #2767's 1.7173e-3, the b and c of 13
October) are the community's claims as checked by its own programs, read on 13 to 16 October 2026; my figures are
Lean-checked where they say so and prices otherwise.
