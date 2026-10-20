# How the lever was found, and how the proof was built

This note tells how I reached this result: where the idea came from, how it was priced and checked, why it waited
five days, and how the Lean proof was built. It is written for people and AI agents who want to understand the
route and go further. It tells what failed as plainly as what worked, and it lists the claims that did not stand.

I worked as on my earlier results: I set the goals and chose the routes, and a team of AI agents (Claude) did the
searches, wrote the proofs, built them and tried to refute one another. The words *array*, *block*, *rank*,
*head*, *job* and *balance* are those of the README, in "The lever, with the numbers".

## 1. The idea: a question about matrix multiplication (14 October 2026, morning)

On the morning of 14 October my best Lean-checked saving was 1.6197705e-3 (by midday, 1.6256629e-3), built on
the community's U8 circuit, and the community's own figure was higher, 1.6797550e-3 on their circuit T74.
Everything my stack proved in Lean was already assumed in the community's prices, so a faithful Lean port of their
best circuit could at most tie with it. To beat it I needed something outside their model.

That morning I had the agents look at abstract routes. One round asked whether two classical ideas from fast
matrix multiplication carry over to this problem: the *laser method*, which builds a good object out of weaker
pieces, and Schönhage's *τ-theorem*, which says what a mixture of different-sized problems is worth. The answer
came in two halves:

- **The laser half is impossible here.** The agents proved that, in a simple model of jobs, combining, powering
  or mixing circuits can never create a saving that no single circuit has. So the hope that a cheap object could
  be powered up to a big saving died that morning, with a proof.
- **The τ-theorem half holds, exactly.** A circuit that asks for jobs of several sizes at once is worth the root of
  a balance with all those sizes on the job side. And there was a corollary, a *truncation lemma*: if a circuit's
  first (or last) operation on some array is a block on that array alone, the circuit is worth strictly more as a
  mixed circuit with that block removed and the array's job made smaller.

On paper, applied to an earlier community circuit, the truncation lemma promised +8.7 %. Two things were open:
whether the community's blocks really were at the start of their arrays, and whether the theorem, proved for
the simplest kind of block, extends to the community's blocks of rank 48.

## 2. Counting, proving and pricing (14 October 2026, late morning)

Three agents worked in parallel, and three more tried to refute them.

- **Counting.** On T74 the agents followed every idle block of rank 48 through the community's schedule. Of 22,748,
  only 11,010 (48.4 %) are heads; the other 11,738 sit on arrays that other steps work on first. None is a
  *tail* (the last operation on its array), and running the circuit backwards only swaps the two kinds. Two
  independent programs gave the same count.
- **Proving.** The agents extended the τ-theorem to blocks of any rank. The result, *Theorem 1G*, needs one
  condition: every block that is paid for must be smaller than every job (the rule is *separated*). The
  construction behind it is a tree of calls run in lockstep, in which every call gets exactly the jobs of each size
  it needs, and the left-over jobs are done by fresh calls of the whole program. A naive order of the calls can
  deadlock; the agents' order cannot. Two skeptics read it and did not refute it.
- **Pricing.** Cutting every head down to a job of 25 axes priced at +3.79 %. A skeptic then pointed out that this
  rule is not separated: the 11,738 other blocks of rank 48 are larger than the new job of 25, so Theorem 1G does
  not cover it. The rescue was to split those blocks into two of rank 24 each, which keeps the rule separated but
  costs back most of the gain: +0.90 %.

Then the planning agent found the variant that is now this result: cut each head by only 24 of its 48 axes. The
job on a head array becomes 49 axes, every paid block stays at rank 48 or below, and the rule is separated with
nothing to split: +1.02 % on T74. Another skeptic checked it, including the question whether a block may be cut
along half of its axes in both the community's semantics and mine.

The plan for the Lean proof was large: a new recursion in my engine, the head lemma, a new seam for the Fourier
chain, and a port of T74 with a network that hosts helper storage across copies of the circuit. The agents
estimated 7,000 to 13,000 lines of Lean. Three cheap gates came first, and all three passed by half past twelve:
the head count in my own conversion of T74 (11,010 again), the skeptic on the half cut, and a design note for the
new recursion with its own skeptic. The build started.

At about 1:20 pm I paused it: I wanted every machine on the abstract routes for a while, and the smaller gains over
the community could wait. The folders were kept, so that the build could resume where it stopped.

## 3. Five days later: the community had moved, and the lever moved with it (19 October 2026, morning)

On the morning of 19 October the agents read the community's repository again. The community had gone past the
figure the build had aimed at: their best circuit without a known refutation was now Dugongue's #3945, at
1.7306632e-3, and the lever on T74 (1.6968621e-3) was 1.95 % below it.

But the newest circuits use the same schedule as T74, and the agents re-priced the lever on them. On #3945, 57,312
of the 137,952 idle blocks of rank 48 are heads (41.5 %), checked by the head test, by its negative control and by
a recount with no chronology. The half cut gives 1.7501109e-3, +1.12 % over #3945 and above every complex-side
construction claimed by others on the community's repository, the disputed ones included (some research notes there print higher ceilings or targets, which are not constructions). The full cut, by contrast, is worthless there:
its break-even is 57,262 heads, only 50 below the count.

I said go, with the instruction to run the cheap gates first and then build everything at once, in parallel, on
five Colab machines.

## 4. The gates (19 October 2026, 08:30 to 10:15)

- **Is #3945 sound where nobody had checked it?** Two of its parts had not been re-run by anyone: the native check
  of its 26 widened terminal frames and the phase check of its new width-5 helpers. An agent ran the community's
  own code for both on a Linux machine, then the exact operator replay over all 360 copies, then #3945's complete
  stated verification from scratch, with nothing reused (1 hour 8 minutes). Everything passed, and so did the
  fallback, #3815. Nothing on the community's repository refuted either one, and nobody there used heads or a
  recursion with two job sizes.
- **Are the heads heads in my own conversion?** An agent converted both words of #3945 into my certificate format
  and followed every array through my own Lean network, not the community's schedule: the same 57,312 heads, with
  every one of the 1,321,920 hosted roles inside a window of my network.
- **How big is the port?** The sizing agent found the one thing nobody expected: #3945's two words use *borrowed*
  helpers with *warps*, free changes of coordinates by a fixed linear map, a device of Cinna Davis's that my Lean
  route had never handled. It also found that the T74 words of the paused build had never actually been built in
  Lean, so the word pipeline had to be tested on them first.

## 5. The build (19 October 2026, from 08:35)

The work ran as 23 agents in parallel, each on its own part, with interfaces written down in files and
the heavy compiles spread over the five machines. A few things are worth telling:

- **The paused code was resumed, not rewritten.** The flow and the cost bound of 14 October (838 lines) compiled
  unchanged; the head lemma of 14 October (718 lines) needed four rounds of fixes and then compiled.
- **The figure did not lose a grid step.** The plans expected the Lean figure to sit one unit of 10^-10 below the
  price, 1.7501108e-3, as in every earlier result. The agent who proved the rate fact found that the rate point and
  the seam's exponent no longer need to lie on the grid, so the figure is the price itself, rounded down:
  1.7501109e-3. Lean also proves that the rate fact fails one unit higher.
- **The seam needed no new envelope.** OpenAI's statement already allows a constant factor in the cost, which the
  first plan had missed and the design note had noticed; the agents proved the bridge and tested the whole chain on a control instance, through Compare
  and the official comparator, before the real engine existed.
- **The table of hosts had to be cut finely.** In a first test on T74's table, Lean's kernel checked the cross-copy
  table quickly (77 seconds on one machine) only when each piece held one kind of host; pieces with several kinds
  blew up.
- **The pieces did not fit at first.** When every part had compiled, the certificate that the heads come first had been written for the plain network, but #3945's two words borrow helpers through warps and need the warp version of the network. An agent of a parallel round, on the older circuit #2979, had just ported that certificate to the warp network in a form that works for any circuit; the integration reused those files unchanged and wrote only the part specific to #3945, about 480 lines. The whole proof then compiled at the first attempt, in about half an hour.
- **One scare of the agents' own making.** While the official comparator held the build folder, a Lean command run there as the wrong user made Lake delete its library packages and fetch them again. The comparator had already read what it needed; the packages were restored from an identical copy, and the run that followed rebuilt nothing. VERIFY.md warns about it.
- **None of the stopping rules fired, and the fallback was not needed.** #3945's words passed the community's complete checks on the agents' machines, so the fallback circuit #3815 was never built; just before the final comparison, and again at the end, no complex-side construction claimed by others on the community's repository was at or above the result (some research notes print higher ceilings or targets, which are not constructions).
- **The checks took longer than the build.** The comparator's two kernels needed an hour and a half to re-check every kernel decision of the table of hosts and of both words; the re-created leaderboard verifier, which compiles the whole proof from source and then runs both kernels, took 3 h 21 min.
- **One limit for the leaderboard.** Its verifier takes at most 500 files in a submission, and the proof's closure has 700 modules, so a submission would first need the generated data files merged into fewer files, and the merged proof checked again.
- The last check passed on the evening of 19 October.

### The smaller circuit, and every head kept (19 October 2026, afternoon)

- **#2979.** While #3945 was being built, the agents priced a second circuit: Marius-Paul Dumitrean's #2979, about
  one fifth of #3945's size (666 runs and 271,728 hosted helper roles, against 3,240 and 1,321,920). In my own
  conversion it has 16,086 heads, and the half cut priced at 1.7556473e-3, above the lever on #3945. Smaller and
  higher, it was the better target, so I had it built at once by the same route, in parallel with #3945; most of
  the Lean files are shared.
- **Every head kept.** An agent then reopened the community's allocation models. In both circuits an empty block of
  rank 48 can sit in a position where it is a head or in one where it never is, at the same price, and the
  community's solvers had split them as they happened to. Re-running the community's own allocation programs with
  every head position kept empty changed only the table of where helper storage goes: 80,640 heads on #3945
  (a price of 1.7581520e-3) and 16,576 on #2979 (1.7564708e-3). The community's complete checks passed on both
  new tables, and both were built in Lean as well, each as its own package.

### The full cut (19 October 2026, night)

- **The break-even moved.** On #3945's published table the full cut had been worthless: its break-even is 57,262
  heads and the table had 57,312. With every head kept the table has 80,640, far above it, and the agent that built
  the new table noticed that the full cut there prices at 1.7806935e-3, 1.28 % above the half cut.
- **The gate.** Another agent checked, before any Lean, that the full cut is a real program in my lowering: the
  57,312 other blocks of rank 48 are whole idle blocks, splitting a block in two is an identity that needs no
  condition on the block, and with every paid block at rank 24 or less below the smaller job of 25 axes, Theorem 1G
  applies. It priced the cut by four independent routes, which agreed to the last digit.
- **The build.** The pieces were small: a general bridge for the full cut, its instance on the same table, the rate
  module, and the Fourier chain at the new exponent; the rest of the half-cut package was reused unchanged. They
  compiled the first time, and the same checks as for the half cut followed (`package/leverF`).

## 6. What did not work, and the claims that did not stand

- **+8.17 % on T74**, cutting every idle block: only 48.4 % of them are heads, so the structure it assumed is false.
- **+3.79 % on T74**, cutting the heads to a job of 25 and keeping the other blocks of rank 48: not separated; no
  theorem covers it.
- **+8.7 % on the earlier circuit #1042**: 46 % of its idle blocks fail the head test, and #1042 is worse than T74
  under every separated variant.
- **Tails.** The other half of the idle blocks cannot be used without changing the community's controller: the
  bridges that read their arrays would have to be reordered, which would be a new construction.
- **Packing the mixed rule into a single circuit for the old recursion.** It loses a fixed fraction of the whole
  saving per level; only a new recursion gets the full value.
- **The full cut on #3945's published table** is at break-even and was not built; with every head kept it pays, and it was built (section 5).
- **The T74 target** was overtaken by the community within five days. The lever itself was not: it carried over to
  the newer circuit unchanged, because the newer circuits keep T74's schedule.
- **An apparent contradiction.** An earlier census had found no idle block that could be merged into a larger one.
  That was a different question (merging, inside the old recursion) and does not conflict with the head test.

## 7. What has not been tried

[open-directions-4.md](open-directions-4.md) has each lead with a first experiment.
