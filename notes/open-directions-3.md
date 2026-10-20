<!-- swap-certificate: keep figures. This note records the figures of 13 October 2026; tools/swap-certificate.sh leaves it unchanged. -->

# Open directions after the fifth result

This note hands on what may come after the fifth result: the Fourier transform of every length with the saving
1.5699e-3, proved in Lean on Ryan Shea's tree, from Jeremy Kothe's circuit with mutable sources and same-run
hosting. Each lead says what is known, what it might give, and a first experiment. It is written for people and
for AI agents who want to go further. [fifth-result-story.md](fifth-result-story.md) tells how the day went; the
two earlier notes, [open-directions.md](open-directions.md) and [open-directions-2.md](open-directions-2.md),
still hold for my own circuits, except that section 8 below settles the block-size-3 question they leave open.

*Added on 19 October 2026, for the fifth revision.* Much of what sections 5 to 8 summarise from the agents'
working notes is now shown in [research/](../research/README.md), and later results there supersede some
statements of section 8: the block-3 bound of 2.8 moves per array is proved in Lean with no hypothesis left, and
block 2 never saves in any game whose free steps touch boundedly many entries per row, chain maps included
(Theorem 2W, argued by hand, not checked), so 2 in 3 at block 2 is no longer open there. Where this note and
research/ disagree, research/ has the later verdict. The leads of the lever are in
[open-directions-4.md](open-directions-4.md).

*Added on 14 October 2026.* Two things changed after this note was written. First, the full stack of section 2
was built and checked in Lean that night; it is now the headline of this repository, 1.6197705e-3
(`community/u8xmhh`), and the README calls it the fifth result. Second, section 8 settles a question that
[open-directions-2.md](open-directions-2.md), section 5, [what-did-not-work.md](what-did-not-work.md) and
[for-network-coding-readers.md](for-network-coding-readers.md), section 3, present as open: every unit of block
size 3 needs at least 14/5 moves per array, so 3 arrays in 8 moves is impossible (the words are defined just
below). Section 8 is brought up to the
state of that morning. The rest keeps the figures and words of 13 October: there, "the fifth result" means the
result of that afternoon, 1.5699305e-3 on Kothe's circuit, and the percentages are over it. Against the new
headline, the next rungs of section 1 are worth +0.36 % (no source copies) and +0.59 % (hosting across runs).

**Words.** As in the story note: a *saving* is 1 - z; a *price* is a saving computed by a program for a
construction that exists, without a Lean proof; a *what-if* is a price that relies on a step nobody has built
or proved; the *community* is CrocSwap/integer-mult-bounds; a *run* is one use of the small circuit inside the
program; *sources*, *targets* and *helpers* are its input arrays, output arrays and spare arrays; a *move* is
one costly step on one array, which takes it to a larger *frame* (the coordinate system it is written in),
and a *climb* is a run of moves; a *bank* is an array of its own that keeps helper pieces until the end of a
run; *hosting* keeps them in idle data arrays instead; an *exit frame* lets a helper stop its climb early. A
*unit* is the finite network that the whole program repeats; its *block size* is the number of directions
along which it transforms each array. Every figure below says what it is. Percentages are over the fifth
result, 1.5699305e-3, unless the line says otherwise. Sections 5 to 7 are for specialists: they use the
community's own terms and those of coding theory and quantum information, and a reader new to these can skip
them.

## 1. What is left on the ladder

The gap between Shea's Lean figure for Kothe's circuit (1.2331685e-3) and the community's own price of the same
circuit (about 1.6049e-3) is made of steps that a Lean proof pays for and the community assumes. The fifth result
closed about nine tenths of it with two proofs, mutable sources and same-run hosting. Exit frames, the third
proof of the day, work on the community's newer U8 circuit. These are the next rungs.

| Step | Saving | Status when this was written |
|---|---|---|
| the fifth result: mutable sources and same-run hosting, on Kothe's circuit | 1.5699305e-3 | Lean-checked |
| hosting across runs, on Kothe's circuit | 1.5785076e-3 (+0.55 %) | price; needs a change in the order of the network |
| all three steps on the U8 circuit: exit frames, mutable sources (5 copies kept), same-run hosting | 1.6197705e-3 (+3.17 %) | price when this was written; Lean-checked that night, now the headline (`community/u8xmhh`) |
| the same with no source copies at all | 1.6256506e-3 (+3.55 %) | price; needs one generalisation in Lean |
| the same three steps with hosting across runs | 1.6293239e-3 (+3.78 %) | price; no Lean format for such tables yet |

Each price reproduces, as a control, Shea's figure and every Lean-checked figure of the day (1.3552606e-3,
1.3868856e-3 and 1.5699305e-3, the last with byte-identical tables). An earlier what-if for all three steps,
1.6361e-3, was made before real host tables existed and is superseded by the lines above.

## 2. All three steps on the U8 circuit (done that night)

**What it is.** The U8 circuit of the community (Jeremy Kothe's, on earlier work by Chafik Boukhalfa, icekylinx,
DreamingOfClouds and others that README.md credits) already has two Lean-checked results here: exit frames
alone, 1.2704023e-3 (`community/u8x`), and exit frames with mutable sources, 1.3868856e-3 (`community/u8xm`).
Hosting is proved for Kothe's "atm3" circuit. Putting hosting on the U8 circuit stacks all three.

**What was done.** The hosting proof in `Work/XHost` is written for Mamano's gauged interface, which the U8
packages also use, so its network modules went over unchanged. The agents' hosting programs chose hosts within
each run for the U8 certificate: 310 banks, against 3,771 without hosting. The price was 1.6197705e-3, +31.35 %
over Shea's figure. That night the agents built the instance on a Colab machine over a copy of the tree of
`community/u8xm`, and it passed every check: Lean's kernel, the comparison with OpenAI's statement and its
failure one unit stronger, the official comparator, a host table with one entry off by one rejected inside
Lean, an offline re-creation of the leaderboard's verifier, and two agents whose only job was to refute it. I
moved the headline to it by hand. The package, with every log, is `community/u8xmhh`; the README, "What is new
in the fifth revision", item 1, lists the checks.

**Credit.** The circuit is the community's, as above. Ryan Shea's notes name exit frames as
the next step for the Lean route, and the charge for them follows the community's written lemmas about where a
helper may end. Hosting is DaysSky's idea, and the hosted network is built from Nil Mamano's files in Shea's
tree.

## 3. The last five source copies

**What it is.** On the U8 circuit, five of the eight central totals sit in source arrays, and the step that
reads the totals reads them only from helpers. So `community/u8xm` keeps the copies of those five sources.

**What is known.** Generalising that reading step so that a total may sit in a source is an estimated 200 to 400
lines of Lean. With it, every copy goes: 1.3928618e-3 without hosting, 1.6256506e-3 with same-run hosting
(prices).

**First experiment.** State the generalised reading step for one central total in a source, prove it on the
smallest case, and check that the existing control still reproduces 1.3868856e-3.

## 4. Hosting across runs

**What it is.** The fifth result hosts every piece in an array that is idle within the same run. Allowing a
piece to be hosted by an array that is idle in a neighbouring run as well gives more hosts and fewer banks.

**What is known.** On Kothe's circuit the best such assignment prices at 1.5785076e-3, +0.55 % over the fifth
result. The Lean network fixes the order in which runs and stages happen, and its host table names hosts of the
same run only. Hosting across runs needs that order changed so that a piece and its host meet, and a table
format that says which run each host belongs to.

**First experiment.** Extend the table format and the checker `hostSchedK`, keep its soundness proof, and build
a control in which every host is in the same run; it must reproduce 1.5699305e-3 exactly. Only then add hosts
in other runs.

## 5. The copy wall, mapped

The proofs, searches and skeptics' checks behind sections 5 to 7 are in the agents' working notes, which are
not part of this repository, so none of them is shown here.

**What it is.** Every run of the community's circuit copies 8 central totals into lower frames before it can
hand them to the targets. The copies cost 64 moves per run on the complex side and 63 on the bit side. These
are the *copy wall*. One move fewer per run would be worth +2.8 % of the community's complex side; a whole
total fewer, +22 %.

**What the rounds of 13 October 2026 established.** The agents listed 17 assumptions behind the wall and tested
those they could. Inside the community's rules (additions only between arrays in the same frame, frames that
only climb, frames that are subspaces):

| Region | State |
|---|---|
| decoders whose totals are linear in the root vectors | at least 64; only the community's decoder reaches it (exhaustive) |
| any decoder on the complex unit | at least 53 (an exact count of the octads each frame meets) |
| relays: a target reads another target's register | never beats the static floor (proved, any order, any field; upheld by two skeptics) |
| bit side, at most 9 totals, or reads spanning at most 11 dimensions | at least 63 (computer-assisted, star model; upheld) |
| copies shared between stages | none: in the community's construction as built, no total can cross a stage boundary (proved) |
| alignment on frames that are not subspaces | no gain under the community's rules; a pair saves nothing under the engine's wider rules |
| decoders invariant under a symmetry of order 7 or 8 | none below 64 among about 8.7 billion; orders 3 and 4 not finished |
| complex totals that read 6- or 7-dimensional subspaces, cost 53 to 63 | open |

**What might give.** A decoder in the open complex region would be worth +2.8 % to +31 % of the community's
complex side (what-if). Outside the community's rules, my engine allows additions between arrays in different
frames; none of the proofs above says anything there. Removing the copies outright needs either that, frames
that are not subspaces, or another unit.

**First experiment.** Decide the octad Hall test for the three cost-63 shapes on the complex side (9 totals of
width 7; one of width 6, seven of width 7 and one of width 8; one of width 6 and eight of width 7), with the
first frame fixed to one representative per orbit of the octad's symmetry group (order 344,064). It is one
overnight job on a Colab machine. "No" closes those shapes for every decoder over every field in the star model;
"yes" gives the first family of frames worth trying to decode.

**Credit.** The star floor, the sandwich lemma and the decoders are the community's (their corrected
format-limits notes and the decoder notes of several contributors); the closures above are the agents'.

## 6. The helpers, mapped

**What it is.** Helpers are 93 % of the distance between the community's figure and the ceiling of its format.
On the U8 circuit, the 504 complex-side helpers per run split into 392 *plane visitors*, which carry three-
dimensional bundles of inputs (packet planes) between registers, and 112 *M copies*, the only registers that
reach the targets' final frames.

**What is known.**

- The 392 plane visitors meet a published lower bound (Dominik Scholz's packet-plane bound) for this packet
  design exactly. The agents' integer programme reproduces the community's 36 helpers per block (504 per run),
  so that count is optimal inside the community's model of visits.
- None of the 504 could be removed by any of the developed ideas (narrower ends, half reads, one helper for two
  planes, shared rows, Walsh packets). One pair can end a step earlier (about +6e-8 on kappa), and that pair was
  already claimed by another contributor.
- The 112 M copies have no published lower bound. A claimed floor of 8 per group was refuted as stated: they
  cancel content the cut has put into every target, rather than deliver it, so the floor holds only with the
  central values and the cut held fixed.
- Removing all 112 would be worth +4.9 % to +12.7 % of the community's complex side (what-if).

**First experiment.** Take one exact small case: one group of 4 targets, 8 centres and 8 super-packets. Let the
central values and the cut coefficients vary while the chains stay fixed, and ask for the least number of
relays, set up as an index-coding problem. If one group needs fewer than 8, scale it up, on the community's
newest circuit (k8), where a first census counts 480 helpers per run.

**Credit.** The circuit and its packet design are the community's (Kothe, Boukhalfa, Capalbo, Scholz and the
contributors their notices name).

## 7. Does the full E8 unit need 72 copy moves?

**The question.** For the E8 complex unit with all 120 ports, does every decoder cost at least 72 copy moves?
(A *port* is one source-target pair with its label; the E8 unit has one port for each of the 120 lines of the
E8 lattice.)

**Why it matters.** Adding one port raises the copy cost by at most 1 (an argument by one of the agents, not yet
reviewed). So 72 for the full unit would give at least v - 48 for every unit made by removing ports from it.
Every known decoder sits exactly there: 56 for 104 ports, 63 for 111, 64 for 112 (the community's complex
side), 72 for 120. A proof would close the community's open problem of beating the star floor for this whole
family. It would be a closure, not a gain.

**What is known.** 72 is reached. Any decoder needs at least 55 (on 120 ports) and 53 (on 112). Decoders linear
in the root vectors need exactly 72 and 64 (exhaustive). With the community's helpers fixed, every star decoder
with 9 or 10 totals on 112 ports costs at least 64. The hard part is decoders that are not linear: lower
bounds of this kind are minrank problems, which are hard in general.

**First experiment.** A computer proof: SAT or integer programming for "a decoder of cost at most 71 exists
on the 120 ports", with the symmetry of the Weyl group of E8 (order 696,729,600) broken. One overnight run on a
Colab machine. The other route is a hand proof that lifts the octad count from 55 to 72, which needs a new idea.
The quantum reading may help: the copy floor is a minrank problem on the commutation graph of the 120
four-qubit Pauli operators with an odd number of Y's, and work on Kochen-Specker sets of the 120 E8 rays may
give bounds (a guess).

## 8. Small units: the state on the morning of 14 October 2026

In the evening of 13 October 2026 I asked the agents for new ideas again, allowed this time to change the model
itself rather than tune the community's circuit. Most of what came back concerns small units. A unit of block
size m with W arrays transforms each array along m directions, so the plain schedule needs mW moves; a unit
that needs only k < mW moves gives the saving log(mW/k) / log m. "2 in 7" below means 2 arrays in 7 moves.
Most bounds below come from a *certificate*: a family of linear tests on the arrays, with a count of how many
of them one move can settle and of how many each array must have settled by the end; the ratio is a lower
bound on the moves per array. The *engine's game* allows, for free, what the engine of this repository counts
as free: linear relabellings of the addresses and *phase patterns* (multiplications by powers of i that depend
on the address). Apart from the Lean theorem of the first item, none of what follows is in Lean. The proofs,
the skeptics' checks and the searches are in the agents' working notes, which are not part of this repository
(not shown). Each item says who checked it. Later: the fifth revision publishes them, cleaned and each with its
status, in [research/](../research/README.md); the bound of block size 3 below was proved in Lean with no hypothesis
left on 14 October ([research/block-bounds-and-new-rules/](../research/block-bounds-and-new-rules/README.md)).

- **Block size 2: no saving, proved in Lean.** Every unit of block size 2 needs at least 2 moves per array,
  which is the plain count. This is a Lean theorem stated on the engine's own definitions (its words, moves and
  costs), with no `sorry` and only the three standard axioms. The agents built it on a Colab machine, and two
  skeptics rebuilt it on another, one checking its logic and one its numbers; neither found a fault. The same
  Lean files also prove the bound of block size 3 below from one numeric fact, which they take as a hypothesis.
  That fact is computed exactly elsewhere, but it is not yet proved in Lean, so the bound of block size 3 is
  not claimed as a Lean theorem.
- **Block size 3: at most 0.063.** Every unit of block size 3 needs at least 14/5 moves per array, 14W/5 in
  all, in the engine's game and also when every relabelling of addresses and every phase pattern is free.
  Three skeptics upheld it independently: one recomputed its numbers from scratch, one re-derived the theorem
  it rests on against the Lean source of the engine, and one checked its remaining hypothesis. So 2 in 5,
  3 in 7, 3 in 8 and 4 in 11 are impossible, and block size 3 can save at most 0.063. The bound leaves open
  every shape with at least 14/5 moves per array. The best, 5 in 14 (and its multiples), would be worth 0.063.
  A round on it narrowed it without settling it: a 5-in-14 unit would have to meet the bound with no slack at
  all, in a precise sense, and no certificate of the families tried goes above 14/5 (three skeptics upheld
  these findings). 6 in 17 (0.052), 7 in 20 (0.044) and the other shapes above 14/5 are open as well. The
  search for 3 in 8 that was to run on a Colab machine for up to three days was stopped, since the bound
  excludes its target.
- **Block size 4: at most 0.171.** Every unit of block size 4 needs at least 41/13, about 3.154 moves per
  array, in the engine's game and also with every relabelling and phase pattern free. Three programs written by
  different agents found the same count behind it (one in floating point), and a skeptic rebuilt it with its
  own code in exact arithmetic. So 2 in 6, 3 in 9 and every "W arrays in 3W moves" are impossible, and block
  size 4 can save at most about 0.171. Open: 2 in 7 (worth 0.096), 3 in 10 (0.132) and 4 in 13 (0.150). The
  family of certificates behind 41/13 can never exceed 10/3 at block size 4, and no certificate of the other
  families tried reaches 3.5, so none of them can exclude 3 in 10 or 2 in 7 (a skeptic upheld this). A search
  for 2 in 7 found none, and a proof upheld by a skeptic shows that a 2-in-7 unit would have to hold, in some
  array at some point, a sum of two different terms of the same input. 3 in 10 has not been searched.
- **Every block size from 4 on.** In the engine's game, every unit needs at least 2.9447 moves per array, and
  at least W + m - 1 moves in all, so a unit for one array needs all m of its moves. Two skeptics upheld each
  bound. Both are far below the plain m moves per array when m is large. The two main families of certificates
  the agents use can never certify 4 moves per array or more: one at any block size, the other at block sizes
  3 to 10 (both proved; a skeptic upheld them). So from block size 5 on, those families cannot rule out a saving
  of 0.1; other kinds of certificate might.
- **Free relabellings that are not linear.** In OpenAI's RAM model a non-linear relabelling of the addresses,
  such as a Toffoli gate on the address bits (it flips one bit when two others are both 1), costs one pass over
  the data, like a linear one, which is of lower order than the transform. The engine could count it as free
  in that model only if it gets a table for AND next to the table for XOR that it already uses (an estimated
  400 to 450 lines of Lean, not written). They do not help at block sizes 2 and 3: the bounds above hold with
  them free (three skeptics upheld this). The bound 41/13 at block size 4 also allows them; searches with them
  at block size 4 have not finished.
- **Free maps that mix neighbouring digit columns.** Inside the recursion, an address splits into columns of m
  bits. A *chain map* XORs bits of each column into the next column. It costs a fixed amount of work per entry,
  so it is free in OpenAI's model, and a saving that uses it carries over to the Fourier figure through
  OpenAI's reduction (a skeptic upheld both). Chain maps break the certificates above, so the agents proved the
  bounds again. The bounds of block sizes 2 and 3 (2 and 14/5 moves per array) still hold when every free map
  acts on all columns alike, Toffoli gates and phase patterns included, and when any maps are added on the
  first and last few columns. For chain maps that another chain map undoes, block size 2 keeps its bound of 2,
  and block size 3 keeps 105/41, about 2.56 per array. Two skeptics upheld these; one found a step false as
  written and repaired it. Open is the plain chain map that XORs each bit into the same bit of the next column,
  whose inverse is a running XOR over all earlier columns, together with quarter-turn phases (multiplications
  by i). For it, a certificate built on the 27 lines of a cubic surface argues the bound of block size 2, but
  two lemmas it needs are not proved, and a skeptic is still checking it. At block size 3, with chain maps free,
  only 105/41 carries over, so 3 in 8 (worth 0.107) is open in that game. No unit has been found: least-squares
  searches at block sizes 2 and 3 found none, which is evidence, not proof (a skeptic re-ran them).
- **Closed the same evening:** folding small blocks into the idle stretches of the network (no block
  qualifies); a "border unit" at block size 2 (an artefact of a model that forgets part of each term); units
  built from Kerdock sets and from mutually unbiased bases, two classic families of codes and bases (being
  unbiased removes the side information a unit needs); and larger E8-like units built from lines whose inner
  products are all odd or zero (killed by a bound on the number of such lines).

If a unit at block size 4 turns up, the write-up will also say where else free non-linear relabellings might
help: transforms built on the Walsh-Hadamard transform, other transforms with large symmetry groups, and
lower-bound proofs that rely on a group of cheap operations.

## 9. Small things, nearly free

- **The Walsh-Hadamard corollary** of any new headline: one term, built and checked in under an hour on a Colab
  machine with the agents' script (`community/wht-c588mhh/scripts-as-run/`). It has not been built for the
  full stack yet.
- **The fourth-order rate inequality** (`term_le4`) in the certificates that do not use it yet: +1.7e-9 on the
  mutable-sources certificate, not measured on the others. Worth adding only at the next rebuild.
- **The quantum dictionary as a search tool.** A normal form from quantum computing (Bravyi and Maslov's, for the
  operations of the story note's dictionary) says that between two additions that touch an array, one run of moves
  of the right size suffices (checked by computer at small sizes only; not shown). That shrinks any search over
  small units.

## Sources

- The packages of this repository named above, each with its `result.json`: `community/c588mhh`,
  `community/c588m`, `community/u8xm`, `community/u8x`, and the hosting programs in
  `community/c588mhh/generators/`.
- The agents' working notes of 13 October 2026, which are not part of this repository: the ladder and the price
  table of the full stack, the copy-wall map with its 17 assumptions, the helper map, the alignment and symmetry
  rounds, the parked question on 72 copy moves, the evening's new-ideas round, the bound of 14/5 moves per
  array at block size 3 with its three skeptics' checks, the certificates and searches at block size 4 with
  their skeptics' checks, and the study of non-linear relabellings in OpenAI's RAM model. Their figures are
  quoted here with their status. Section 8 also rests on the agents' notes of the night of 13 to 14 October
  2026: the Lean files of the block-size-2 theorem and the two skeptics' rebuilds of them, the round on 5 in 14,
  the bounds of block sizes 4 to 6, the bounds for every block size, and the rounds on maps that mix columns.
  None of them is part of this repository.
- Ryan Shea's notes on exit frames: `research/complex-transfer-588-596/HANDOFF.md` of his repository at the
  commit below, section "Next mathematical test".
- The community's repository, CrocSwap/integer-mult-bounds: Kothe's circuits #588, #799 and #869 (the U8
  circuit), with Boukhalfa's frame changes #741 and icekylinx's gauge44 construction #645; the community's
  owner maps for the U8 circuit, which the Lean results here do not use (#870 by Kothe with Hao Yan's #858,
  Boukhalfa's #890 and #685, Capalbo's #942), and DreamingOfClouds' bit supplier #589; the k8 word #1018 and
  its packing #1042; DaysSky's hosting #442; the endpoint lemmas #593 (Kothe) and #599 (DreamingOfClouds);
  Zhihao Chen's #425; the corrected format-limits notes (`research/format-limits-corrected/`),
  the decoder-floor notes and Scholz's packet-plane bound (`research/structural-gain/`); the claimed helper pair
  in #1090 and #1148; the open problem of the star floor and the decoder work in #717, #939, #954 and #968.
- Ryan Shea's repository, shea256/fourier-transform-below-nlogn, at `52394849c7ecf3412d39b925d75fc7fdb98994d8`,
  folder `lean/complex588`.
- OpenAI's cost model: the RAM namespace of `UniformFourier.lean` in openai/math.
- Quantum information: Aaronson and Gottesman (quant-ph/0406196); Bravyi and Maslov (arXiv:2003.09412);
  Maslov and Roetteler (arXiv:1705.09176); Vandaele et al. (arXiv:2302.07040); Dehaene and De Moor
  (quant-ph/0304125); Garcia, Markov and Cross (arXiv:1210.6646, arXiv:1711.07848); Lisoněk, Raussendorf and
  Singh (arXiv:1401.3035) for the 120 E8 rays; Koetter and Kschischang (arXiv:cs/0703061) for the subspace
  distance.
