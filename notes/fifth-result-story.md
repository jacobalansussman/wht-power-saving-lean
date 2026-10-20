<!-- swap-certificate: keep figures. This note records the figures of 13 October 2026; tools/swap-certificate.sh leaves it unchanged. -->

# How the fifth result was reached

*Added on 14 October 2026.* That night the agents built the last Lean line of the ladder in section 3: exit
frames, mutable sources and same-run hosting together on the community's U8 circuit. It passed every check,
and two agents whose only job was to refute it found no fault (section 6). Its saving, 1.6197705e-3, is now the
headline of this repository, and the README calls it the fifth result. This note keeps the figures and the
words of 13 October: below, "the fifth result" means the result of that afternoon, 1.5699305e-3 on Kothe's
circuit.

This note tells how I reached the fifth result on 13 October 2026: a Lean proof that the discrete Fourier
transform of every length, and convolution, can be computed with the saving 1.5699e-3. It came in two steps on
the same day: first 1.3553e-3, with a proof that the circuit may write into its sources, then 1.5699e-3, with a
proof that leftover helper storage may live in arrays that are idle. It is written for people and AI agents who
want to understand the route and go further. It tells what failed as plainly as what worked, and it lists the
claims of the day that turned out wrong.

I worked as on the earlier results: I set the goals and chose the routes, and a team of AI agents (Claude)
did the searches, wrote the proofs, built them and tried to refute one another. Agents whose only job was to
refute a result checked both steps; for the fifth result one of them compiled its new modules again from a
fresh copy of the sources. Some of the other results had no such check, or its log is not in this repository;
their packages mark those checks "not run" or "not shown".

## Words used here

- **Saving.** The program needs O(n (log n)^z) operations; the saving is 1 - z. Larger is better.
- **The community.** The contributors of the public repository
  [CrocSwap/integer-mult-bounds](https://github.com/CrocSwap/integer-mult-bounds). They work on bounds for
  integer multiplication. One part of their figure, the *complex side*, is the same kind of quantity as my
  Fourier saving. Their headline figure, *kappa*, also needs a second part, the *bit side*.
- **Circuit.** The small network of additions that runs over and over inside the program. Its *sources* hold
  what one run reads, its *targets* receive what it writes, and its *helpers* are spare arrays that hold
  partial sums. A *register* is one array that the circuit writes into. A *gate* is one step of the circuit.
- **Frame, move, climb.** Each array is written in a coordinate system, its *frame*, which grows as the work
  on it proceeds. A *move* is one costly step on one array: it takes the array to a larger frame. A *climb* is
  a run of moves. By the end of each run every array must reach the *full frame*. Additions are free, but only
  between arrays at the same frame (the *frame rule*). The *dimension* of a frame is its size.
- **Network.** The arrangement of stages in which the circuit runs over the data. Shea's proof uses a network
  of nine stages (the *nine-stage network*); the edge between two stages is a *stage boundary*.
- **Bank.** An array of its own in which the conversion keeps the pieces of the helpers that must last to the
  end of a run (the *leftover pieces*: the part of a helper's climb still ahead of it). Every bank widens every
  run, so fewer banks give a larger saving. A bank's unused room is climbed at the end, and that *padding* is
  paid.
- **Engine, kernel program, Fourier chain.** The *engine* is the part of my Lean proof that turns a checked
  network into one program for the lengths 2^k with a time bound, the *kernel program*. Its last theorem is the
  *end theorem*. The *Fourier chain* is the set of files that carries the kernel program through OpenAI's
  reduction to every length; only its first file, the *seam*, names the end theorem. The *rate* step is the
  inequality that turns the cost of one unit into the exponent.
- **Unit.** The finite network that the whole program repeats. My own units are built on the *E8* lattice, a
  highly symmetric set of 120 *lines* (directions) in 8 dimensions; each pair of arrays in a unit carries a
  *label*, a word of bits.
- **Central totals.** Sums formed in the middle of a run of the community's circuit that every target reads.
  Copying them to where the targets can read them costs 64 moves per run, the *copy wall*. A *decoder* is a way
  of making those copies; a *relay* passes a total from one target to another.
- **Ceiling of the format.** The saving that the community's circuit format would give if its helpers cost
  nothing.
- **Lean-checked.** A Lean proof that passed Lean's kernel, the comparison against the statement, its failure
  one unit stronger, and the official comparator.
- **Price.** A saving computed by a program for a construction that exists, without a Lean proof.
- **What-if.** A price that relies on a step nobody has built or proved yet.
- **Skeptic.** An agent whose only job is to refute a claimed result.
- **Control.** A run whose answer is known in advance, made to show that a check can tell right from wrong.
  The *baseline* is the figure a new step is measured against.

## 1. Where things stood that morning

On 10 October my fourth result, 8.7625e-4, was the largest saving with a Lean proof that I knew of. Three days
later it was third. On the leaderboard [bigospeedrun.com/exact-dft](https://bigospeedrun.com/exact-dft),
which lists Lean proofs of OpenAI's Fourier statement, Ryan Shea was first with 1.2331685e-3 and Nil Mamano
second with 1.0221e-3. Shea's proof takes Jeremy Kothe's circuit from the community and carries it through my
Lean engine, with an extension of the engine by Mamano.

The community itself had moved further. Its selected figure that morning had a complex side of 1.6644e-3, and
during the day it moved twice more, to about 1.675e-3 by the afternoon. Those figures are conditional: the
community checks a finite construction exactly, but the figure also assumes steps that nobody has proved.

So I asked the agents a plain question: can the community's best work and mine, together, give a larger
number? I wanted the largest number, not polish on my own circuit. Building on others is the point.

## 2. The survey: where the gap is

The agents read the community's records, Shea's repository and the leaderboard, and priced everything with an
exact copy of the community's own pricing program. The copy reproduced the community's figures to the last
digit before it was trusted.

The finding that shaped the day: the community prices Kothe's circuit, the one Shea proved, at a complex side
of about 1.6049e-3. Shea's Lean figure for the same circuit is 1.2332e-3, 23 % lower. The circuit is the
same. The difference is made of costs that a Lean proof has to pay and the community's price assumes away. The
community's repository lists those assumptions by name, among them "mutable-source realization" and
"whole-call hosting".

The agents split the gap into four costs:

| Cost a Lean proof pays | Worth about | What would remove it |
|---|---|---|
| copies of the 112 sources, made because the engine required that a run never write its sources | 10 % | a Lean proof that the engine allows *mutable sources* |
| helpers kept in banks of their own, instead of in data arrays that sit idle | 16 % | a Lean proof of *hosting* |
| helpers that must climb to the end of a run, where newer community circuits let them stop early | 3 % | a Lean proof of *exit frames* |
| the layout of the helper banks | 0.4 % | a better table; no new proof |

Each of these, once proved, is also one of the community's assumptions turned into a theorem.

## 3. The ladder

The agents priced each step and each combination, and every price was first checked against Shea's own figure
as a control. This is the ladder at the end of the day; the figures are Fourier savings.

| Step | Saving | Over Shea | Status at the end of the day |
|---|---|---|---|
| Shea's proof, as published | 1.2331685e-3 | | Lean-checked (his) |
| a better bank table, with a sharper rate inequality | 1.2383291e-3 | +0.42 % | Lean-checked (`community/c588-repack`) |
| exit frames, on the community's U8 circuit | 1.2704023e-3 | +3.02 % | Lean-checked (`community/u8x`) |
| mutable sources, on Kothe's circuit | 1.3552606e-3 | +9.9 % | Lean-checked (`community/c588m`) |
| the same with a packed bank table | 1.3566302e-3 | +10.01 % | price; data ready, not built in Lean |
| exit frames and mutable sources, on U8 (5 copies kept) | 1.3868856e-3 | +12.47 % | Lean-checked (`community/u8xm`) |
| **mutable sources and same-run hosting, on Kothe's circuit** | **1.5699305e-3** | **+27.3 %** | **Lean-checked: the fifth result (`community/c588mhh`)** |
| mutable sources and hosting across runs, on Kothe's circuit | 1.5785076e-3 | +28.0 % | price; needs a change in the order of the network |
| exit frames, mutable sources (5 copies kept) and same-run hosting, on U8 | 1.6197705e-3 | +31.35 % | price; Lean data generated, build started that evening; built and checked in Lean that night (section 6), now the headline (`community/u8xmhh`) |
| the same with no source copies at all | 1.6256506e-3 | +31.83 % | price; needs one generalisation in Lean |
| the same three steps with hosting across runs, on U8 | 1.6293239e-3 | +32.13 % | price; no Lean format yet |

(The percentages are over Shea's 1.2331685e-3. Only the lines marked Lean-checked are theorems.)

I took the steps in order of gain for the work. Mutable sources came first: its Lean proof was one idea, and
the data pipeline already existed. Hosting came next: it was the largest step, and the day's pace made it look
like an afternoon rather than the weeks first guessed (section 8). Exit frames ran alongside on the U8 circuit.

## 4. Mutable sources

**The idea.** My engine proves that a run of the circuit does the right thing, and then undoes the run's side
effects with an undo step derived from the run. Until now it proved this under one rule: no step of the run
writes a source. Kothe's circuit, as the community designed it, does write its sources. To fit my rule, Shea's
conversion copied each source into a helper first and let the circuit work on the copy. That cost 112 extra
helpers and 1,008 extra moves in every run.

The new proof drops the rule. The undo step then restores the sources together with the helpers. What the
proof needs is that no output of the run feeds back into a source or a helper, that the run can be inverted,
and that at the end of the run the sources and helpers that the undo step touches stand at the same frame, so
that the step is legal there. Zhihao Chen had written the community's argument for this step; the
agents wrote the Lean proof for my engine, in 7 modules and 1,429 lines, in about two hours.

**The certificate.** With the copies gone, Kothe's circuit has 504 helpers per run instead of 616, and 5,916
moves instead of 6,924. The agents generated the certificate and all of the Lean data with Shea's pipeline,
adapted where the copies had been built in. An annealing search chose the table of helper banks: 4,032 banks
instead of 5,040. Exact Python mirrors of every check that Lean would make passed before anything was built.

**The build and the checks.** The 70 new modules compiled on a Colab machine over Shea's built tree in 29
minutes. Then came the axioms, the comparison with OpenAI's statement (`RESULT: PASS`), the same comparison one
unit stronger (it fails, as it must), the official comparator with both of its kernels, negative controls
inside Lean, and a scan of the sources. Two skeptics then worked independently. One compared the statement
files with OpenAI's file, found that only the exponent and seven names differ, and read the logs; it did not
run Lean, and its written verdict is its only record. The other rebuilt all 70 modules on a second Colab
machine, and on its build the comparison passed and failed one unit stronger. Both upheld the result, and both
reports are in the package (`community/c588m/logs/`). Last, the agents re-created the leaderboard's own
verifier offline: on a packaged submission they ran its three phases with scripts of their own, inside the same
stand-in sandbox, because its controller needs tools the Colab machines lack, and all three passed. That made
1.3553e-3, +9.9 % over Shea, by mid-afternoon.

**One incident.** To run the comparator's stand-in sandbox, the agents gave the build folder to an
unprivileged user. Lake, the build tool, then no longer recognised its package clones as its own, deleted
them and fetched them again without their compiled files. The agents restored them from the built tree.
VERIFY.md, section 11.3, warns about it.

## 5. Hosting

**The idea.** In every run of the nine-stage network, the data arrays sit idle during some stages. DaysSky
proposed, for the community's circuits, to keep the helpers' leftover pieces in those idle arrays rather than in
banks of their own ("idle-axis hosting", with a written proof). The community's figures assume it. For a Lean
proof it is the largest single step: on Kothe's circuit it takes the banks from 4,032 to 111.

**The lemmas first.** Before building anything large, the agents proved two lemmas that their plan for the
proof needed:

- L1 says how a run acts on the place where a helper's piece is kept, whatever that place held before. While
  writing it, the agents found that its content was already proved in Shea's tree (`B401.gauged_step`,
  `B401.gauged_Wn`, `GX.entW_spec`, `GX.entW_form`, `GX.gauge_step`), on Dugongue's entrance gauges (which let
  a helper enter a run part of the way up its climb), Mamano's files and my lemma `core_nested`. The new file
  only restates it in the form the plan needed, and the credit for it belongs there. Before formalising it, the
  agents checked it on 98 exact instances.
- L2 says that a relabelling of addresses, which the engine counts as free because it is one pass over the
  data, can carry a helper's piece into the part of its host's frame that the host has not yet climbed, and
  that the host's climb then finishes with the piece inside it.

Two skeptics rebuilt L1 and L2 on the Mac from their own copies and upheld them; their logs and written
verdicts are in `community/c588mhh/supporting/`, with the two files. In the end the checked result uses
neither of them (next paragraph).

**The hosted network.** Then the agents built the network itself, `Work/XHost`, 17 modules, and they built it
from Mamano's files in Shea's tree rather than from scratch. The hosted stage and its certificate are copies of
Mamano's network files with hosting added. The step that moves a piece into its host is adapted from his bank
lemmas: a host is treated as a bank whose climb starts and ends at the host's own frames, so the step L2 proves
was not needed in the form L2 gives it. The checker `hostSchedK`, which reads the table of hosts and banks
inside Lean's kernel, is modelled on his bank-table checker, and comes with a proof that any table it accepts
is sound. The network ends in a theorem of the same form as the one Shea's proof ends in, adapted from Mamano's
file, so that the Fourier chain of the 1.3553e-3 result runs on it with only its seam changed. What is new is
hosting added to that network, the checker's soundness, the host table and the instance. The hosts were chosen
within each run, by an integer programme laid out after DaysSky's.

**A control first.** The same network with an empty host table, in which every piece stays in a bank, had to
give exactly 1.3552606e-3, the figure without hosting. It did, and it passed the comparison and the comparator.
That showed, before anything depended on the host table, that the hosted network gives back the result without
hosting when nothing is hosted.

**The build.** The modules were compiled on the first Colab machine with plain `lean`, one at a time, in
rounds while the agents fixed errors; the logs of every round, failures included, are in the package.
The machine's disk was saturated by other jobs (one module that compiles in about 12 seconds took 19 minutes), so the agents
compiled from a copy of the compiled files in memory. Once, a build run as the machine's administrator made
Lake delete its copy of Mathlib, for the same ownership reason as in the incident of section 4; the agents
restored it.

**The checks.** Standard axioms; `RESULT: PASS` against OpenAI's statement at 1 - 15699305/10^10; `RESULT:
FAIL` one unit stronger; the official comparator accepted it with both kernels after compiling the hosted
modules again from source (14 minutes), and rejected the statement one unit stronger. A first comparator
attempt had stopped after one second with "Child exited with 1"; its cause is not in the log, and the accepted
run started a minute later. Then a third skeptic copied the sources afresh and compiled the 32 modules of the
result into a build folder of its own: all compiled, the comparison passed, the stronger statement failed, and
the axioms were the standard three. It did not run the comparator.

That made 1.5699305e-3, +27.3 % over Shea's figure and +15.8 % over the mutable-sources result. The
Walsh-Hadamard corollary on the same kernel program, 1.5699306e-3, was built and checked on the second Colab
machine the same evening.

What was not done: the leaderboard's own verifier has not been run on the hosting result, and the comparator
ran only inside a stand-in for its sandbox, because the Colab machines lack the kernel feature (Landlock) it
needs.

## 6. The other things that were built and checked

- **The better bank table on Shea's own certificate, 1.2383291e-3.** An integer programme found a table of
  5,009 banks whose padding costs the least possible, by the agents' exact count, and a sharper (fourth-order)
  inequality in the rate step took the figure to the exact limit for that table. The same programme showed
  that no table of 4,999 banks, the fewest that could hold all the pieces, exists; whether a table between
  4,999 and 5,009 banks exists was not decided. A skeptic re-checked the built tree: the build was up to date,
  the axioms standard and the comparison passed; its comparison one unit stronger did not finish
  (`community/c588-repack/logs/colab/skeptic-v-c-pack/`).
- **Walsh-Hadamard corollaries.** My theorem `BlockWHT.wht_main_of_engine` turns any kernel program of the
  engine into a Walsh-Hadamard theorem in one term. Applied to Shea's kernel program, to the repacked one, to
  that of the mutable-sources result and to that of the hosting result, it gives 1.2331686e-3,
  1.2383292e-3, 1.3552607e-3 and 1.5699306e-3, each compared with my published Walsh-Hadamard statement.
- **Exit frames on the community's U8 circuit, 1.2704023e-3.** The agents wrote the exit-frame theory in Lean
  (about 1,400 lines), translated the community's circuit into the engine's form and built it on the Colab;
  the comparison with OpenAI's statement passed and the comparator accepted it.
- **Exit frames and mutable sources together, 1.3868856e-3.** Five Lean modules join the two theories. The
  certificate keeps 5 source copies: on U8, five of the eight central totals sit in source arrays, and the step
  that reads the totals reads them only from helpers. It passed the same checks. No independent rebuild of
  either U8 package is shown in the repository.
- **All three steps on the U8 circuit, that night, 1.6197705e-3.** The hosted network needed no change: its
  twelve modules went over unchanged, byte for byte. The agents' hosting programs chose 2,016 hosts for 4,092
  pieces in each run, and 310 banks remained, against 3,771 without hosting. The agents built the 32 new
  modules on the second Colab machine in 15 minutes, over a second build of the exit-frames-and-mutable-sources
  package. The comparison with OpenAI's statement passed, and failed one unit stronger; the official comparator
  accepted it and rejected the stronger statement; a host table with one entry off by one was rejected inside
  Lean; and the three phases of the leaderboard's verifier, re-created offline with the agents' scripts,
  passed. Then one agent re-checked the statement and replayed the comparator, and another compiled the 32
  modules again on the first machine over its own copy of an earlier build. Neither found a fault. That is
  +31.35 % over Shea's figure and +3.17 % over the fifth result of the afternoon. Its logs are in
  `community/u8xmhh/`.

## 7. What failed, and why

Most of the day's other rounds were aimed at the community's own figure, not at a Lean proof. None of them beat
it. That is worth saying plainly, because the community's construction turned out to be tightly optimised
inside its own model, and the routes below are now mapped for anyone who tries them next.

The proofs, searches and skeptics' checks of these rounds are in the agents' working notes, which are not part
of this repository, so none of them is shown here.

**Five grafts onto the community's circuit.** The agents tried ideas from my earlier idea waves on the
community's best circuit: targets that read other targets to retire helpers, merging helper slots, reordering
gates to remove steps where nothing can proceed, a cheaper decoder for the copies, and a cheaper copy on the
bit side. Every one gave
exactly nothing or a small loss. The best legal reordering came out 0.035 % lower. The community's own
attempts on the same levers had also failed.

**The copy wall.** Every run of the community's circuit copies 8 central totals at a cost of 64 moves. One
move less per run would be worth +2.8 % on the complex side. No decoder, relay or schedule went below 64.
Two closures survived both skeptics. First, a relay never beats the *static floor*, the least cost when every
target gets its totals with no relay (proved in the community's rules, for any order of the steps and over any
number system). Second, on the bit side every decoder with at most 9 totals, or with reads that span at most 11
dimensions, costs at least 63 (computer-assisted, in the community's "star" model, in which every copy is made
from the totals directly). The region between the proved floor of 53 and the known 64 stays open.
Open directions has the map.

**Sharing copies between stages.** The output of one stage of the network is the input of the next, so a
copy made in one stage might serve the next one too. It would have been worth twice the complex side as a
what-if, and capped by the bit side in practice. A first proof that it cannot work had a counting bug and was
refuted. A second round then replayed the two most promising routes exactly, with the cost of cleaning up
counted, and both failed. The skeptics upheld that. Later a proof showed that in the community's construction,
as it is built, no total can cross a stage boundary at all.

**Alignment on other frames.** The agents asked whether frames outside the ones every checker uses could get
past the copy wall, by letting a target cancel content it does not know, as interference alignment does in
wireless networks. Under the community's rules this is impossible (proved). Under my engine's wider rules a new
kind of cancellation exists, but for a pair it costs exactly what two separate totals cost, and at full scale
it leaves the 64 where it was. The same round worked out an exact dictionary with quantum information: a frame
is a stabilizer state (a standard kind of quantum state), the free steps of the engine are the Hadamard-free
Clifford group (the quantum operations built without the Hadamard gate), and a paid move is one Hadamard
gate. The dictionary is exact, but almost all of it was already known: in the quantum literature, in
the community and in my own notes of 9 October.

**The helper round.** Helpers are 93 % of the distance between the community's figure and the ceiling of its
format, so the agents went after them. None of the 504 complex-side helpers per run could be removed. One
helper on the bit side could (+0.012 % to +0.119 % on that side, which does not move kappa; not checked by a
skeptic), and a pair of helpers could end one step earlier, worth about +6e-8 on kappa. That pair turned out to
be already claimed inside another contributor's open submission, so there was nothing to send.

**Symmetry and other fields.** In the afternoon I asked for ideas that change the model rather than tune it:
lattices with more symmetry than E8, designs with a prescribed symmetry group, and superoptimisation of the
community's circuit. All came back negative. The 16-dimensional Barnes-Wall lattice does worse than E8 as a
unit: 72 of its lines are pairwise non-orthogonal, which forces wide labels, and its ceiling on the complex side
is 1.2863e-3. Grouping the *delivery helpers* (those that carry the central totals to the targets) into
"Walsh packets", worth 56 helpers per run as a what-if, breaks the frame rule in every group. A legal
regrouping of idle steps gains +0.029 %, but another contributor had already proposed it. About 8.7 billion
decoders that a symmetry repeating after 7 or 8 steps leaves unchanged were searched, and none costs less than
64. A superoptimiser (a program that searches exhaustively for shorter equivalent code) found no rewrite that
improves the community's circuit.

**New ideas, in the evening.** After those rounds I asked for ideas again, this time allowed to change the
model itself. They are leads, not results, and open directions lists them with their state.

## 8. Claims that turned out wrong, and the corrections

These are the day's wrong claims, the agents' and mine. Each was caught by a skeptic or by a later check.
Where a correction concerns the rounds of section 7, its check is in the agents' working notes (not shown).

| Claim | What was wrong | What stands |
|---|---|---|
| "Natural exits" give 1.486e-3 on the Lean route | The engine's undo step needs every helper restored at its full frame; helpers that simply stop where their last gate leaves them break it | Exit frames priced as the community charges them, 1.2704e-3, which is Lean-checked |
| Two proofs that the copy wall holds (stage sharing, and cuts through *isotropic* lines, those orthogonal to themselves) | One had a counting bug; both left cases open or forbade steps instead of pricing them | Both count as no progress. The counting lemma for isotropic cuts survives as a lemma |
| The quantum link "appears new" | The Clifford reading of frames was already in the community, in a contributor's pull request and in my own notes | The dictionary is exact and checked by computer; what may be new is narrow |
| Interference alignment "exactly doubles" the pairs that cancel for free | True only at odd sizes; at size 4 it is 56 against 32 | The doubling at size 9, the case that matters, holds |
| Hosting would take "1.5 to 3 weeks" of Lean | An agent's human-scale guess, which I repeated unchecked | Hosting was built and checked in Lean the same afternoon |
| The first hosting lemma, L1, was new work | Its content was already proved in Shea's tree, on Dugongue's and Mamano's work and my engine | The file restates it; the credit is theirs, and the text says so |
| The second lemma, L2, is what makes hosting legal in the checked proof, and the hosted network is new | The checked network imports neither lemma. It is Mamano's network with hosting added, and its step into the host is adapted from his bank lemmas | The text says what is copied or adapted from whose files; L2 stays as a supporting lemma |
| The end theorem of the hosting follows a file of mine, `Work/Bridge377/End.lean` | No such file of mine exists. Shea's `End.lean`, which the hosting's end theorem adapts, says only that it follows a file of that name, and that the engine it uses is mine | The text credits the end theorem to Mamano's files, built on my engine. The header comment of one compiled file still carries the wrong attribution, and its package says so |
| The leaderboard's verifier was "re-run" on the mutable-sources result | Its own controller never ran: the agents re-created its three phases with scripts of their own, inside the stand-in sandbox | The text says "re-created" and gives the limits |
| Two skeptics' checks of the mutable-sources result had no log in the repository, and one skeptic's check of the repacked table had no verdict | The logs were on the agents' machines, not in the packages | The logs are now in the packages; what has no verdict is marked "not shown" |
| All three steps together on U8 give 1.6361e-3 | An early what-if, before real host tables existed | Re-priced with real tables: 1.6198e-3 within one run, 1.6293e-3 across runs (prices) |
| Sharing rows between targets saves 56 helpers per run | The model charged the circuit for work it already does for free, and left out the delivery helpers | No gain; the recipe would also break the frame rule |
| The 112 delivery helpers are already minimal (8 per group) | They cancel content rather than deliver it, so the claimed floor holds only with more held fixed | Fewer than 8 per group is open |
| A re-run of the helper pairing on the newest circuit | It used the wrong baseline, a second time | Its count of the circuit's helpers is still useful |
| A summary of the isotropic-cut counts | Two numbers were swapped (135 exact classes, 120 numerical) | Corrected in the working notes |
| The list of ports removed on the bit side (a *puncture*) that one agent's search used | Numbered in the wrong file's order | A skeptic re-ran the search on the true list and got the same answers |

Two corrections concern the process. A watcher script that the agents ran marked the second Colab machine ready
when its setup had in fact failed silently, because a timing tool was missing. The mark was removed within a
minute, and from then on the machine was marked ready only when a build and a comparison had passed. And the
release audits found that the text named the wrong commit of the leaderboard's verifier, that the script that
swaps the headline certificate could change the text on checks it had not run, and several of the wrong claims
in the table above; all were fixed before this note.

## 9. What I take from the day

- The largest gains on the Lean side come from proving what the community assumes. Each such proof is also a
  contribution back: one assumption fewer behind the community's own figure. Two of them, mutable sources and
  same-run hosting, took the Lean figure for Kothe's circuit from 1.2332e-3 to 1.5699e-3, about nine tenths of
  the way to the community's own price of that circuit. With exit frames, a third such proof, the same steps
  on the community's U8 circuit gave 1.6198e-3 that night.
- Credit has to be checked as carefully as proofs. One of the lemmas I thought new was already in Shea's tree,
  and the step I had credited to the other was, in the checked proof, adapted from Mamano's files.
- Tuning the community's construction inside its own model gave nothing lasting. The community's many
  contributors move faster at that than my agents did, and their circuit is tight.
- A result is not finished until it has been compared with the exact statement, rejected one unit stronger,
  accepted by the official comparator and attacked by agents whose job is to refute it. Every wrong claim of
  the day was caught by those checks or by the release audits that followed (section 8).

## Sources

- The fifth result of the afternoon, its checks and logs: `community/c588mhh/` (`result.json` lists every check
  with its status and, for each check that passed, its log), and VERIFY.md, section 11. The full stack of that
  night: `community/u8xmhh/`, with its runs of 13 and 14 October 2026 (UTC) in its `logs/`. The other packages
  of the day: `community/c588m/`, `community/u8xm/`, `community/u8x/`, `community/c588-repack/`,
  `community/wht-c588/` and `community/wht-c588mhh/`, each with its `result.json`.
- The README of this repository, section 9, "For the fifth result": the commits, pull requests and dates
  behind the statements above.
- Ryan Shea's repository: shea256/fourier-transform-below-nlogn at commit
  `52394849c7ecf3412d39b925d75fc7fdb98994d8`, folder `lean/complex588`, with its `NOTICE.complex588`. Its
  `Work/Bridge401/` and `Work/Bank401/` are Nil Mamano's files, taken from his submission to the leaderboard
  (bigospeedrun/verifier commit `025396033eede56e3af6e02ed2679ba828381bd5`, as its `SOURCE_MANIFEST.json`
  records); the files of `Work/XHost` that copy or adapt them say so in their first lines.
- The community's repository, CrocSwap/integer-mult-bounds, read on 13 October 2026: its README (the list of
  assumptions), `certificates/selected-result.json`, `certificates/checkpoints/family588.json` (kappa
  0.001602278921983212 for Kothe's circuit #588, complex side 1.6049e-3), and the selected results of that day,
  in turn: #799 (complex side 1.6644e-3, in the morning), #942 (1.6647273e-3) and #1042 (kappa 1.672188e-3,
  complex side 1.6750e-3, in the afternoon).
- People and pull requests behind the statements above: Jeremy Kothe (#588; #799; #869, the U8 circuit), Chafik
  Boukhalfa (#741, the frame changes the U8 circuit carries), icekylinx (#645, the gauge44 construction it is
  built on), the endpoint lemmas behind exit frames (#593 by Kothe, #599 by DreamingOfClouds) and helpers that
  end early (#601 by Kothe, #604 by DreamingOfClouds), Zhihao Chen / jacklightChen (#425, the argument for
  mutable sources; #377, the compressed nine-stage layout), Dugongue (#401, the entrance gauges and banks),
  lydakis (#275, the bank patterns), DaysSky (#363, the bridged layout; #442, idle-axis hosting, with
  `research/idle-axis-hosting/PROOF.md`), the claimed helper pair (#1090 by 1337isnot1337 and #1148 by
  tomdif), the idle re-grain (#1111 by seslly), the community's own failed attempts on the same levers (#909,
  #910, #961), and the earlier Clifford reading of frames (#130 by icekylinx, #288 by eumemic).
- Not used by the Lean results, though the community's figures for the U8 circuit use them: the owner maps of
  #870 (Kothe, with Hao Yan's #858), #890 (Chafik Boukhalfa) and #942 (Jonathan Capalbo), and the bit supplier
  #589 (DreamingOfClouds). README.md, section 9, says more.
- The leaderboard and its verifier: bigospeedrun.com/exact-dft and bigospeedrun/verifier (commit `c92b9b6`,
  the one the offline re-creation used; `community/scripts/bos/` holds the agents' scripts that ran its
  phases).
- The quantum literature named in section 7: Aaronson and Gottesman (quant-ph/0406196), Bravyi and Maslov
  (arXiv:2003.09412), Maslov and Roetteler (arXiv:1705.09176), Koetter and Kschischang (arXiv:cs/0703061),
  and the others listed in `notes/open-directions-3.md`.
- The figures of the ladder that are not Lean-checked are prices from the pipeline's own generators, the
  agents' hosting programs (`community/c588mhh/generators/design-tools/`) and the agents' copy of the
  community's pricing program; each reproduces Shea's 1.2331685e-3, and the Lean-checked figures above it, as
  controls. The prices of the steps not yet built are in the agents' working notes, not in this repository.
- The rounds that failed (section 7) are recorded in the agents' working notes of 13 October 2026, which are
  not part of this repository; their figures are quoted with their status.
