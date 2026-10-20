# The record track: the Lean ladder, the levers that failed, the recursion lever and the frontier

This section covers the race for the record from 13 to 20 October 2026: how I climbed from Ryan Shea's Lean-checked
figure to my own, why every attempt to beat the community's best circuit inside its own model failed, the
mathematics of the *recursion lever* that finally gave a price above every complex-side construction claimed by others on the community's repository,
the checks that came before each build, and where the frontier stood on 19 October. The records themselves are in
their packages, and the story of how the lever was found and built is told in the lever's notes
([notes/how-the-lever-was-found.md](../../notes/how-the-lever-was-found.md)); this section holds the theory, the pricing programs, the gates and the failures.

I did this work with a team of AI agents (Claude) that I directed: I set the goals, chose which ideas to pursue and
decided at each stage what to do next; the agents ran the searches, wrote the proofs and programs, and tried to refute
one another's work. No person has read any of the proofs in this section. What is proved in Lean is said finding by
finding; everything else here is not checked in Lean.

**In one paragraph.** On 13 and 14 October I proved in Lean a ladder of records on the community's circuits, from a
reproduction of Shea's 1.2332e-3 to 1.6257e-3, by proving steps that the community's prices only assume (no source
copies, helpers hosted in idle registers, exit frames). That ladder could never beat the community: its figure for
the same circuits already assumes all of those steps, so a faithful Lean port of its best circuit only ties, and
every lever I found for a margin inside its model died. What worked was a change of recursion. Many of the
community's idle blocks of rank 48 are the first operation on their array; cutting such a *head* by half of its axes
gives a recursion with jobs of two sizes, and Theorem 1G, an extension of Schönhage's τ-theorem to grouped blocks
that I worked out for this work (argued by hand, not checked, in a simple model of jobs; its form in OpenAI's model
is proved in Lean, sections 4.7 and 4.8), prices that recursion by a balance equation. On
19 October this priced at 1.7556e-3 on the community's circuit #2979 and 1.7501e-3 on #3945, above every complex-side
construction claimed by others on the community's repository (some research notes there print higher ceilings or targets, which are not constructions); with the community's allocation re-solved so that every head bin stays empty, the same
circuits give 1.7565e-3 and 1.7582e-3.
Which of these are proved in Lean is stated in sections 4.7 and 4.8. The *full cut*, which removes heads entirely
and splits the other idle rank-48 blocks 24 + 24, prices at 1.7797e-3 and 1.7807e-3 on those allocations; on #3945's
it was built in Lean on the night of 19 October, at 1.7806935e-3, the headline of this revision (section 4.8). On #2979's words a re-solved table (A-fc, section 4.10) prices at 1.7887046e-3 under the full cut, *checked by computer*, not proved in Lean.

## Contents

1. [The question](#1-the-question)
2. [Words used here](#2-words-used-here)
3. [What I tried](#3-what-i-tried)
4. [What came out](#4-what-came-out)
   - [4.1 The Lean ladder](#41-the-lean-ladder)
   - [4.2 A faithful port only ties](#42-a-faithful-port-only-ties)
   - [4.3 Every lever for a margin inside the community's model is dead](#43-every-lever-for-a-margin-inside-the-communitys-model-is-dead)
   - [4.4 The recursion lever: the mathematics](#44-the-recursion-lever-the-mathematics)
   - [4.5 Heads, and the head test](#45-heads-and-the-head-test)
   - [4.6 The half cut and the full cut](#46-the-half-cut-and-the-full-cut)
   - [4.7 The gates and the prices on #2979 and #3945](#47-the-gates-and-the-prices-on-2979-and-3945)
   - [4.8 Every head bin empty: higher designs on the same circuits](#48-every-head-bin-empty-higher-designs-on-the-same-circuits)
   - [4.9 The frontier on 19 October](#49-the-frontier-on-19-october)
5. [What failed or was refuted](#5-what-failed-or-was-refuted)
6. [What is still open, with a first experiment for each](#6-what-is-still-open-with-a-first-experiment-for-each)
7. [The programs and how to run them](#7-the-programs-and-how-to-run-them)
8. [Credits](#8-credits)
9. [Sources](#sources)

## 1. The question

The theorem this repository proves in Lean is OpenAI's statement for the discrete Fourier transform of every length
(and convolution): the transform can be computed in O(n (log n)^z) operations, and the *saving* is 1 - z. A
community leaderboard (the repository CrocSwap/integer-mult-bounds, "the community" below) builds circuits of the
same kind and prices them with its own engine. Its *complex side* b is the same kind of number as my saving, but it
is a price under assumptions that nobody there has proved: that helpers can be hosted, that sources can be written in
place and restored, and other steps of the same kind. Its headline figure κ is for integer multiplication and is
derived from the smaller of b and a *bit side* c. My rule for publishing was: only a Lean-checked figure that beats
the community's figure, and the fair comparison for a DFT figure is b.

On the morning of 13 October the Lean-checked frontier was Ryan Shea's 1.2332e-3 on Jeremy Kothe's circuit #588 (my
own published figure was 8.7625e-4), and the community's selected b was 1.6644e-3 (#799) and rising: by the
afternoon #942 and #1042 had taken it to 1.6750e-3. So the questions were:

1. How far does proving in Lean the steps the community only assumes take a Lean-checked figure?
2. Can anything beat the community's best circuit by a margin, inside its own model?
3. If not, can a different recursion price the same circuits higher?
4. Where was the frontier when the result was built?

## 2. Words used here

- **Saving.** The b in O(n (log n)^(1-b)); larger is better. Written 1.7501e-3, never as a checker integer.
- **Circuit (word), call, replica.** One full pass of a community circuit is a *unit*: a fixed program that transforms
  W *arrays* together. It is run as many *calls*, grouped in *replicas* (T74 means 74 replicas; #3945 has 360).
- **Axis, block, rank.** Every array has 73 coordinates (*axes*). A costly step is a *block*: a transform along r of
  one array's axes, done by a recursive call; r is its *rank*. Additions between arrays are free.
- **Job.** What the recursion asks of one array in one unit: a transform along k of its axes. In the community's
  recursion every job has 73 axes.
- **The balance.** For a recursion of this kind the saving is fixed by one equation: τ is the root of
  Σ (rank of each paid block)^τ = Σ (job size of each array)^τ, and the saving is 1 - τ. With every job of size 73 this
  is the community's *moment formula*. A *mixed* rule has jobs of several sizes.
- **Idle block, head, tail.** The community's allocation leaves some arrays with axes that no helper uses; those axes
  are transformed by one *idle block* (of rank 48 on the circuits here). An idle block is a *head* if it is the first
  operation on its array, acts on that array alone, and nothing is stored in that array before it; a *tail* is the same
  at the end.
- **Separated.** A mixed rule is separated if every paid block is smaller than every job.
- **Half cut, full cut.** Two ways to use heads (section 4.6): cut each head by 24 of its 48 axes (its array's job
  becomes 49), or remove it entirely (job 25) and split every other rank-48 block in two blocks of rank 24.
- **Fill.** The community's allocation of helper storage to idle cells, as a table; changing the fill changes which
  rank-48 bins stay empty, and so which blocks are heads.
- **Price.** A figure computed by a program from a profile; it is not a theorem. A *what-if* is a price that rests on
  a step not yet built or proved.
- **Lean-checked record.** A Lean build with only the three standard axioms (propext, Classical.choice, Quot.sound),
  Compare.lean printing PASS on the positive statement and FAIL one unit higher, and the official comparator accepting
  it; every comparator run here used a stand-in sandbox, not the real Landlock one.
- **Skeptic.** An independent AI agent asked to refute a result, with its own code where possible.

## 3. What I tried

- **13 October: the ladder.** The agents surveyed where Shea's Lean figure lost against the community's price of the
  same circuit and found the gap made of costs a Lean proof must pay and the community's price assumes away: source
  copies (about 10 %), helpers in their own banks instead of hosted in idle registers (about 16 %), exit frames
  (about 3 %) and packing (0.4 %). I had each step proved in Lean in order of gain for the work (section 4.1; the
  story of that day is in `notes/fifth-result-story.md` of this repository). In the same round five grafts of my
  earlier ideas onto the community's best circuit were tried (section 4.3).
- **14 October, morning: beat the community.** Could the ladder's methods give a margin on the community's best
  circuit, T74? A faithful port only ties (section 4.2), and every lever for a margin died (section 4.3).
- **14 October, same morning: an abstract route.** A theory round asked whether the laser method and Schönhage's
  τ-theorem from fast matrix multiplication carry over. The laser half is impossible here; the τ-theorem half holds,
  with a truncation lemma that made idle blocks at the start of their arrays worth more
  ([../existence-routes/](../existence-routes/), section 3). Three tracks then counted the heads of T74, extended the
  τ-theorem to grouped blocks (Theorem 1G) and priced the lever, and three skeptics tried to refute them; a planning
  agent found the half cut and sized the Lean route at 7,000 to 13,000 lines. Three cheap gates passed by early
  afternoon and the build started.
- **14 to 16 October: paused.** I paused the build that afternoon to put every machine on the abstract routes,
  approved it again on 15 October and cancelled it on 16 October before it started. The folders were kept.
- **19 October: the frontier, and the build.** A frontier check found the community past the T74 target, and the
  lever re-priced on its newest circuit #3945. I said go; the gates were re-run on #3945, a sizing agent found the
  community's *warps* (free linear changes of coordinates in borrowed helpers) that my route had never handled, and
  the build ran on five Linux machines of Google's Colab service (8 processors and about 50 GB each). During the day
  a second circuit, #2979, priced higher at one fifth of the size, and I gave it priority. I then found, with a
  re-optimisation round, that the community's allocations can leave every head bin empty at no cost, which gives
  higher designs on both circuits (sections 4.7 and 4.8).

## 4. What came out

### 4.1 The Lean ladder

Every rung below is **proved in Lean**: only the three standard axioms, Compare.lean PASS on the claim and FAIL one
unit higher, and the official comparator accepted it (in a stand-in sandbox). The skeptics named rebuilt or re-read
each one and did not refute it; one rung has no skeptic of its own, as marked. The packages are elsewhere in this
repository; this table only places them.

| rung | construction (on whose circuit) | DFT saving | over Shea | skeptics | package |
|---|---|---|---|---|---|
| Shea's result, reproduced through my engine | Kothe's #588 with 112 source copies (Shea's conversion) | 1.2331685e-3 | - | every price below reproduces it first | Shea's own repository |
| optimal re-packing of Shea's bank table, sharper rate inequality | #588 | 1.2383291e-3 | +0.42 % | five (the packing, its re-pricing, its truth, the statement, the WHT corollary) | [`community/c588-repack`](../../community/c588-repack/) |
| exit frames | the community's U8 circuit (#869, #870, #890) | 1.2704023e-3 | +3.02 % | two (statement, rebuild) | [`community/u8x`](../../community/u8x/) |
| mutable sources (written in place, then restored) | #588 | 1.3552606e-3 | +9.90 % | two, one a rebuild of all 70 modules | [`community/c588m`](../../community/c588m/) |
| exit frames and mutable sources (5 centre copies kept) | U8 | 1.3868856e-3 | +12.47 % | two | [`community/u8xm`](../../community/u8xm/) |
| mutable sources and helpers hosted in idle data registers within one run | #588 | 1.5699305e-3 | +27.3 % | three, one a clean rebuild | [`community/c588mhh`](../../community/c588mhh/) |
| all three on U8 (5 centre copies kept) | U8 | 1.6197705e-3 | +31.35 % | two, one an independent rebuild | [`community/u8xmhh`](../../community/u8xmhh/) |
| the same with no centre copies | U8 | 1.6256506e-3 | +31.83 % | none of its own (replaced within hours) | not shipped |
| the same with a one-pad host table and the sharper lemma | U8 | **1.6256629e-3** | +31.83 % | a fresh rebuild that regenerated every file byte for byte, and a statement skeptic | [`package/u8x0hq`](../../package/u8x0hq/) |

Of the factor 1.318 between Shea's start and the top rung, about 1.032 is the community's better circuit (U8 against
#588) and about 1.277 is overhead that the Lean route had charged and the ladder removed (*argued by hand, not
checked*: from a dissection of my words that a skeptic checked in part). My words carry no saving object of their
own; the saving comes from the community's network and helper circuit ([../why-the-savers-save/](../why-the-savers-save/)).

Higher rungs were priced and never built: on U8 with no copies, hosting across runs at 1.6438659e-3 and deferred
finals at 1.6459856e-3 (*a price computed by a program, not proved*; both need a hosted network across replicas that
was not written at the time). An early what-if for the whole U8 stack, 1.6360855e-3, came before real host tables
existed and was re-priced lower with them. All of these were below the community's figure of the day.

### 4.2 A faithful port only ties

On 14 October the community's best circuit was T74, with b = 1.6797550e-3. Every step my ladder proves in Lean is
already inside the community's price, and my Lean rate theorem certifies the same moment inequality the community's
b solves. So a faithful Lean port of T74 lands at 1.6797548e-3, about 2e-10 *below* b, and on T74 every extra I had
built or priced (hosting across runs, no copies, deferred finals, packing, the sharper lemma, mutable sources, exit
frames) adds exactly 0 (*checked by computer*: prices of each extra on T74's own profile, and the port through my Lean grid; one
AI skeptic re-read the plan and did not refute it). A port would have beaten only the community's *selected* record of
that day (b 1.6754627e-3), not its newest open figure.

### 4.3 Every lever for a margin inside the community's model is dead

Each of the first five was a candidate for a margin of 0.4 % or more on T74, and each failed; the grafts of
13 October were tried on the community's best circuits of that day. One skeptic each upheld the verdicts.

- **Walsh packets** (grouping the delivery helpers so that additions need no common frame): wrong at the matrix level
  as built; every correct repair costs at least 24 in rank per group, against a break-even of 50 moves in all
  (*checked by computer*).
- **Copies and plane ports in my engine:** the shortest climbs never leave the community's frames, so my wider engine
  adds nothing there (*argued by hand, not checked*, in part as a written proof). The open "bundle" question is worth at most
  +1.48 % as a what-if, with no construction.
- **Idle blocks merged into larger blocks** (my engine has no rank cap): the idle blocks of each family never overlap
  in time, and bridges wall them in (*checked by computer*: 0 qualifying blocks).
- **Better fills:** no other fill of T74's idle cells, and no other split of its 666 calls between its two helper
  circuits (the community uses 417 and 249), gives a higher community price, by a nine-way balance count and an exact
  LP-duality certificate (*argued by hand and checked by computer*). This is about the community's own price; a fill
  does change the lever's price, through the number of heads (section 4.8).
- **Thirteen generated levers:** each fails against the same wall: one paid rank-1 move added per call costs -2.82 %
  (*checked by computer*).
- **Grafts on the circuits #799 and #942 (13 October):** letting targets read targets retires exactly 0 helpers;
  helper-slot reuse passes exactly 0 merges; the best legal reordering of gates loses 0.035 % (248 attempts); no
  9-line puncture of the bit side beats copy cost 63 (by proof, by SAT and by enumeration; a 10-line lead reaches 62
  only in a finer model, and its skeptic weakened it). κ can rise at most 4.13 % before the bit side binds
  (*checked by computer*; the skeptics' texts are not shown).

The day's verdict was plain: no route beat the community that day except one lever outside its model.

### 4.4 The recursion lever: the mathematics

**The job model.** Count cost per entry: a *job of size s* is the Walsh-Hadamard transform on s chosen bits of a
register, and a rule turns jobs of the target sizes into jobs of the paid sizes. The community's recursion is the case
with one job size: a paid block is one recursive call, and the arrays of the next level are slices of that call's
batch. The full definitions are in [proofs/theorem-1g.md](proofs/theorem-1g.md), section 0.

**The τ-theorem (rank-1 paid moves).** A mixed unit with W_k arrays owed k axes and M paid moves of rank 1 is worth
exactly the exponent τ that solves Σ W_k k^τ = M, with no loss (*argued by hand, not checked*; two AI skeptics found
no error). The proof and its programs are in [../existence-routes/](../existence-routes/), section 3, with the
companion result that combining or powering pieces never creates a saving no piece has (the laser half fails), and
the *truncation lemma*: if a unit's first (or last) operation on some array is a block on that array alone, removing
it and shrinking that array's job always raises the unit's value (*argued by hand, not checked*; no skeptic of its
own).

**Theorem 1G (grouped paid blocks).** The community's paid blocks have ranks up to 48, not 1. Theorem 1G extends the
τ-theorem to them under one condition, *separation*: every paid rank below every job size. Then the cost is at most
C n^τ (ln n)^κ, so every saving below 1 - τ is reached with a fixed constant (*argued by hand, not checked*; two AI
skeptics, one re-deriving the logic with its own simulation, one building the construction explicitly with a generic
deadlock test, found no error). The construction is a *dilated lockstep tree*:

- units are indexed by round sequences; the units of a class run round j together, and round j of a class feeds one
  child class and nothing else;
- the Pascal weights u(a) make every child's demand equal its parent's supply exactly (Lemma B); this is where
  separation is used;
- new input enters only at the top, as frontier jobs of sizes between n/m and n (Lemma C);
- what leaves the tree at the bottom (*exits*, jobs of size at most r_max (n/m)^θ with θ = ln r_max / ln k_min < 1)
  carries exactly the Φ of the sources, and is done by fresh calls of the whole program (Lemma D);
- the number of units has no log factor (Lemma E);
- the order "issue round j of every unit, then finish the child class, then round j + 1" cannot deadlock (Lemma F);
  a careless order does, and a skeptic's random matching across classes deadlocked as it should.

Two side results bound what the construction can do: the induction closes with leftovers exactly under a condition
on their Φ (Lemma A), and a flow whose leftovers carry all the Φ can never give exact O(n^τ) through it (Lemma A4).
So *exact* O(n^τ) with no log factor is **open** for grouped blocks; it costs nothing in the exponent. Every
numerical check of the construction (exact balance, sources in range, Φ conserved, exits within the bound, unit count,
the cost of the whole schedule up to n = 12,800, and hypothesis (E) of the induction at every n up to 12,800 on five
test rules and, for the T74 what-if rule, at every n up to 20,000 and at every 997th n up to 3,000,000) was re-run for
this revision and reproduced the saved outputs exactly (*checked by computer*).

**The Lean form.** On 19 October the recursion was written in Lean inside OpenAI's RAM model as part of the lever
build: the flow, the cost bound, the program of the lockstep tree, the head lemma and the rate facts. The Lean cost
bound needs no log factor: for every exponent above a rate point the exits' share of Φ shrinks, so the induction
closes directly; the constants exist but are astronomically large, which OpenAI's big-O goal allows. Status:
*proved in Lean*, as part of each of the five lever packages (in `package/leverF` and `package/leverB`, the modules under `Work/Mixed/` and `Work/MixedCost/`), with the official comparator accepting each result and three skeptics upholding each ([proofs/theorem-1g.md](proofs/theorem-1g.md), section 6).

**Why a head can be cut even though its array is read later.** Transforms along different axes commute. If the
first operation on array w is a block H_S on w alone, then the rest of the circuit, applied to y = H_S x, outputs the
transform along the other axes on w and the correct outputs everywhere else. So the head's axes can be taken out of
this unit's job on w and done elsewhere in the recursion (*argued by hand, not checked*; two AI skeptics read it on
14 October and found no error, one of them including the version for the complex kernel, where the reduced job is a
standard job up to free maps on the raw input).

### 4.5 Heads, and the head test

The head test follows every array through the community's schedule (local sweeps, bridges between arrays, idle
sweeps, hosted helper roles, terminal events) and marks a rank-48 idle block a *literal head* if nothing touches its
array before it.

- **T74.** 11,010 of the 22,748 idle rank-48 blocks are heads (48.4 %): children 8 and 10 of the circuit's pair 3
  (6,150 + 4,860). The other 11,738, in children 12 and 13 (pair 0), are neither heads nor tails, because a bridge
  that adds pair 0 into pair 1 reads their arrays after them. No idle block is a tail, and running the circuit
  backwards only swaps the two kinds (*checked by computer*: two independent programs on the community's controller,
  then a third in my own conversion of the words and my own Lean network, which also checks that every one of the
  319,680 hosted roles lies inside a window of my network).
- **Controls.** Moving the bridge that reads pair 3 ahead of pair 3's idle sweep drops the heads to 0; deleting the
  gate that reads pair 0 afterwards turns all 11,738 into tails; the derived idle and paid profiles equal the
  community's own, digit for digit.
- **#1042**, the circuit of the first what-if: 2,329 literal heads, 1,735 "net" heads (only hosted roles come before
  them) and 3,486 neither, of 7,550 idle blocks of several ranks. 46 % fail the test (*checked by computer*, two
  programs).
- **#3945.** 57,312 of 137,952 (41.5 %): child 8 28,672 and child 10 28,640; pair 0 (children 12 and 13, 80,640 blocks)
  is never hosted and never a head. A recount with no chronology gives the same numbers, the template file is
  byte-identical to T74's, and a skeptic recomputed everything with its own method, inferring every helper width from
  the fill itself (*checked by computer*, three programs). In my own lowering and network: the same 57,312, with every
  one of 1,321,920 hosted roles inside a window of my network.
- **#2979.** 16,086 of 28,190 (57.1 %): child 8 8,064 and child 10 8,022. The circuit ships no native word, so the
  first count used a stub for the one field the test needs (the list of helper origins), checked against #3945, where
  the stub equals the real field; my own transport of #2979's words then produced the real field, which equals the
  stub, and my lowering and network give the same 16,086 (*checked by computer*; the first count was re-run for this
  revision and reproduced exactly).

Why the counts differ between circuits: the community's allocation places hosted roles on pair 0 or on pair 3, and
every pair-3 bin it leaves fully empty is a head. #2979 leaves more of them empty than #3945 does. That observation is
the start of section 4.8.

### 4.6 The half cut and the full cut

Write I for the idle rank-48 blocks, h for the heads among them and W for the arrays. The two separated ways to use
heads change the balance as follows (everything else in the circuit is unchanged):

    half cut:  (other blocks) + (I - h)·48^τ + h·24^τ   =  (W - h)·73^τ + h·49^τ
    full cut:  (other blocks) + 2(I - h)·24^τ           =  (W - h)·73^τ + h·25^τ

- **Half cut.** Each head keeps 24 paid axes and its array's job shrinks from 73 to 49. Both sides fall by about the
  same amount, but x^τ grows more slowly for large x, so the step from 48 to 24 removes slightly more than the step
  from 73 to 49, and τ falls. The largest paid rank is 48 and the smallest job 49: separated, with nothing to split.
  Every single head helps, so there is no break-even.
- **Full cut.** Each head disappears and its array's job becomes 25. Then the other rank-48 blocks are larger than
  the smallest job, so they must be split into two blocks of rank 24, and each split costs something
  (2·24^τ > 48^τ). The full cut beats the community's own figure only when at least 41.5 % of the idle rank-48 blocks
  are heads, on every circuit here.
- **Why not more.** Removing the heads and keeping the other rank-48 blocks whole gives +3.79 % on T74 and +4.20 % on
  #3945, but that rule is not separated (48 is not below 25) and no theorem covers it (section 6).
- **Is a block of rank 48 the same as two of rank 24?** In the community's semantics the 48 axes of an idle block
  are orthonormal lines of the port's basis, so any 24 of them and the other 24 are again orthonormal and span the
  right spaces; in mine a block acts by its flat list of axis transforms, which a 24 + 24 split leaves unchanged. A
  gate skeptic checked both on all 112 ports and 2,016 random splits (*checked by computer*, and *argued by hand* for
  the general case); the Lean form proves the split as a lemma of the head-lemma layer (its status is that of the
  Lean form, section 4.4).

Every price in this section and the next is the root of the balance, bracketed with exact rational interval
arithmetic to a width of 1e-16 by [programs/lever-prices/lever_table.py](programs/lever-prices/lever_table.py), which
also reproduces every figure stated by the earlier rounds and finds the break-evens (*checked by computer*). The T74
and #3945 figures and design B's half cut were also reproduced by AI skeptics with their own code; the #2979
figures, #3815's half cut, design A's half cut and design B's full cut were computed by at least one earlier program
of the build rounds, and this program reproduces them; the full cut on design A is new here and has no second check.

| circuit | heads (of idle rank-48) | its own b | half cut | full cut | full-cut break-even |
|---|---|---|---|---|---|
| T74 | 11,010 (of 22,748) | 1.6797550e-3 | 1.6968621e-3 (+1.02 %) | 1.6948358e-3 (+0.90 %) | 9,443 |
| #2979 | 16,086 (of 28,190) | 1.7290340e-3 | **1.7556473e-3** (+1.54 %) | 1.7744822e-3 (+2.63 %) | 11,702 |
| #3945 | 57,312 (of 137,952) | 1.7306632e-3 | **1.7501109e-3** (+1.12 %) | 1.7307680e-3 (+0.006 %) | 57,262 |
| #3815 (fallback for #3945) | 57,312 (of 137,952) | 1.7305784e-3 | 1.7500242e-3 (+1.12 %) | 1.7306832e-3 | 57,262 |

These are prices, not theorems; which of them are proved in Lean is section 4.7.

### 4.7 The gates and the prices on #2979 and #3945

Before any heavy Lean, each platform had to pass three gates: is the community's circuit sound where nobody had
checked it, are the heads heads in my own conversion and network, and does the price stand.

**#2979** (Marius-Paul Dumitrean, merged into the community's repository on 16 October as research, built on Jeremy
Kothe's #2957 and Mason Hale's #2917). It is about one fifth of #3945's size: 666 calls and 271,728 hosted roles,
against 3,240 and 1,321,920.

- *Soundness.* The community's registered check of #2979 and its complete replay from scratch (fresh bootstrap,
  fresh native chains for both helper circuits, complete admission and release, with the author's exact library
  versions) both pass on a Linux machine; the complete replay took 55 minutes. Twelve merged independent audits on the
  community's repository found no defect and nothing there refutes #2979 or its lineage (*checked by computer*, the
  community's own programs run unchanged; *read on the community's repository*).
- *My conversion.* My transport of #2979's two words into native form replays exactly to the community's endpoint on
  all 632 columns, its block histogram and helper widths equal the community's stated ones, 7 of 7 deliberately wrong
  versions are rejected, and the canonical words equal the community's own fresh native words in all 25 fields
  (*checked by computer*).
- *Heads and price.* 16,086 heads in my network, the global profile equal to the community's digit for digit, and the
  half cut 1.7556473e-3, +5.5e-6 above the half cut on #3945 (*checked by computer*).
- *Lean.* *Proved in Lean* (`package/lever2979`): OpenAI's statement at 1 - 17556473/10^10, a saving of 1.7556473e-3; the official comparator accepts it with both kernels and rejects it one unit higher, the leaderboard's verifier re-created offline accepts the merged form, and three skeptics (statement, rebuild, price) upheld it. It is a fallback of the revision's headline, design B.

**#3945** (Dugongue, open on the community's repository until it was archived on its main branch on 20 October, on the staged line of #3815 by ChaloupeMusic).

- *Soundness.* Two parts of #3945 had not been re-run by anyone: the native check of its 26 widened terminal frames
  and the phase check of its new width-5 helpers. Both pass when the community's own code is run on its words, and
  #3945's complete stated verification passes from scratch in one uninterrupted run (1 hour 8 minutes), as does the
  fallback #3815's with fresh natives (*checked by computer*, the community's own programs).
- *Heads and price.* 57,312 heads in my own lowering and network, the profile equal to the community's, and the half
  cut 1.7501109e-3 (*checked by computer*). The full cut is at break-even there (50 heads clear) and was not built.
- *Warps.* Every circuit on which the lever was priced near #3945's (eight of them, #2979 among them) uses the same
  72 warps from Cinna Davis's borrowed helpers, so changing platform avoids none of that work and none of it needs
  redoing (*checked by computer*: a comparison of the word files of these circuits); the best circuit without
  warps is T74, whose half cut is below #3945's b (a keyword search of the community's pull requests found no circuit
  without warps above T74's b; a search, not a proof).
- *Lean.* *Proved in Lean* (`package/lever`): OpenAI's statement at 1 - 17501109/10^10, a saving of 1.7501109e-3; the official comparator accepts it with both kernels and rejects it one unit higher, the leaderboard's verifier re-created offline accepts the unmerged form, and three skeptics (statement, rebuild, price) upheld it.

One thing the plans got wrong, usefully: every earlier record lost one unit of 10^-10 between its price and its Lean
figure, and the plans expected the same here. The rate point and the seam's exponent no longer need to lie on that
grid, so the Lean figure is the price rounded down (the rate fact is proved in Lean at that figure, and proved false
in Lean one unit higher). Design A in section 4.8 is the one exception.

Other circuits priced the same day, with heads measured on their published fills: #3653 1.7543510e-3 (95,387 heads)
and #3682 1.7529511e-3 (66,304 heads), both below #2979 (*checked by computer*, one program). With every head bin
empty both rise to 1.7578999e-3 (section 4.8).

### 4.8 Every head bin empty: higher designs on the same circuits

The head counts in section 4.5 depend on where the community's allocation puts its hosted roles. I found that in both
of the community's allocation models an empty rank-48 bin can sit in children 8 and 10 (where it is a head) or in
children 12 and 13 (where it never is) at exactly the same base price, and the community's solvers split them
arbitrarily. Forcing every non-empty rank-48 bin into children 12 and 13 keeps the circuit's profile and price exactly
and makes every bin of children 8 and 10 a head: 224 per replica, the most possible. In designs A and B only the fill
changes; the words, the order of calls and the dispatch are the community's, unchanged. Design C keeps #2979's words
but runs them at 60 replicas instead of 74, so its order of calls and its fill are new. (The scheduling programs are
the community's; the change is one constraint or one column cost.)

| design | circuit | heads | its own b | half cut | full cut |
|---|---|---|---|---|---|
| A | #2979 (T74 scale) | 16,576 (published fill: 16,086) | 1.7290340e-3 (unchanged) | 1.7564708e-3 | 1.7797095e-3 |
| B | #3945 (T360 scale) | 80,640 (published fill: 57,312) | 1.7306632e-3 (unchanged) | 1.7581520e-3 | 1.7806935e-3 (built in Lean) |
| C | #2979's words at T60 | 13,440 | 1.7290358e-3 (the cap of #2979's words) | 1.7564726e-3 | - |

The same re-allocation, found by a local search from the published fills, gives 80,640 heads on #3682 (T360) and
109,312 on #3653 (T488), both at 1.7578999e-3: above design A, below design B (*checked by computer*, one program:
the fills pass the community's owner engine with its five rejection controls, but not those circuits' own complete
admissions; no skeptic).

- *Validity and price of design B* (*checked by computer*; one AI skeptic upheld it): it passes #3945's own owner
  verification unchanged (plain and optimised, with its four rejection controls), the earlier staged audit, and
  #3815's whole verifier at every fill-level step (it stops only where the recorded admission binds the published
  fill, as any new fill must); the head test finds 80,640 literal heads and its bridge-first control 0. The skeptic
  re-ran the community's verifiers with its own harness, recounted the heads with code of its own (the published fill
  gives 57,312, a fill with one role moved gives 80,639), derived the profile itself and priced design B by two routes,
  one of them certified on a grid of 1e-24.
- *Validity of design A* (*checked by computer*, no skeptic yet): it passes #2979's registered allocation check with
  its five rejection controls, and the head test finds 16,576 heads (bridge-first control 0). A gate agent independent
  of the one that made the fills re-ran the community's three checks of #2979 on a fresh copy, and counted the same
  16,576 heads with its own code and in my own lowering and network.
- *Validity and price of design C* (*checked by computer*; one AI skeptic upheld it): the owner engine that #2979
  vendors passes it with seven rejection controls, and the head test finds 13,440 heads. The skeptic also ran #2979's
  allocation check with only its three constants for 74 replicas set for 60 (it passes, five controls rejected),
  counted the heads with its own code, derived the profile itself and priced the design by four routes, which agree
  to 24 digits. A new number of replicas needs new admission inputs and a new Lean port; neither exists.
- *The operator replay.* The exact operator replay with which the community's audit refuted the idle-axis chain
  (section 4.9), run with the community's modules on every owner of both designs, finds 0 differences (2,029,945 moved
  characters on B, 594,132 on A), and its wrong-time and wrong-stream controls change the operator as they should
  (*checked by computer*, by the same independent gate agent).
- *Admission.* The community's complete admission replay with the new fills (fresh natives, owners, discovery, outer
  checks, release) passes for both designs, from scratch, with the community's own figure for the circuit and its
  release receipt: #3945's full replay with design B's fill (119 of 119 controls rejected, 1 hour 29 minutes) and
  #2979's complete replay with design A's fill (fresh bootstrap too; 139 of 139 controls rejected, 1 hour 18 minutes)
  (*checked by computer*, the community's own programs, run by the same independent gate agent).
- *Optimality* (*checked by computer*, two independent programs, no skeptic): an exact LP dual, checked in 90-digit
  arithmetic over all 841 patterns, shows that, with the half cut, no fill of #3945's words, with any chronology, exceeds design B by more than 4e-34; a second certificate, built separately from exact basis duals in 80-digit arithmetic, gives the same
  ceilings. On #2979's words, with the half cut, nothing exceeds 1.7564726e-3, which design C reaches at 60 replicas; at #2979's own 74
  replicas the best whole-number fill is at most 1.7564710e-3, 2.4e-10 above design A. These bounds hold for the half cut only: with the full cut the same LP bounds are
  1.7823049e-3 on #3945's words and 1.7893658e-3 on #2979's words (*checked by computer*), above the full cut of design B. Fills come close to them: 1.7822982e-3 on #3945's words and A-fc, 1.7887046e-3, on
  #2979's (section 4.10), both *checked by computer* and not built in Lean.
- *Lean.* Design A's rate enclosure reaches 1.7564707e-3 and not 1.7564708e-3 (its root, 1.75647080068e-3, sits 0.007
  grid steps above that point), so a Lean proof of design A states 1.7564707e-3, one unit below its rounded price;
  there the rate fact still holds one unit higher, and Lean proves it false from two units higher.
  Design A: *proved in Lean* at 1.7564707e-3 (`package/leverA`); the comparator accepts it, the re-created verifier accepts the merged form, and three skeptics upheld it. Design B, the full cut: *proved in Lean* at 1.7806935e-3 (`package/leverF`), the headline of this revision; the comparator accepts it, the re-created verifier accepts the merged and the unmerged forms, a separate check found the merged form identical to the unmerged one theorem by theorem, and three skeptics upheld it. Design B, the half cut: *proved in Lean* at 1.7581520e-3 (`package/leverB`); the comparator accepts it, the re-created verifier accepts the merged and the unmerged forms, a separate check found the merged form identical to the unmerged one theorem by theorem, and three skeptics upheld it.
- *History.* The gate note of 14 October had already proposed a head-maximising fill on T74 as an optional upside
  (all 16,576 head bins empty, priced then on T74's own profile at 1.7506306e-3 under the full cut), but it was never
  run, because the build was paused that afternoon. It was found again on 19 October.

On the published fills #2979 is both the smaller circuit and the higher lever. With every head bin empty that is no
longer so: #3945's design B is the highest of all, 1.7e-6 above design A, at five times the size, and the smallest
design, C, is 1.7e-6 below B too. Which of them leads this revision is said in the top README.

### 4.9 The frontier on 19 October

Read on the morning of 19 October and again at midday and in the afternoon (all read-only):

- **The community's selected record** was unchanged: κ 1.6726602e-3 (#1094; b 1.6754627e-3, with the bit side
  1.7334270e-3 of #589). Selection had last run on 13 October; hundreds of later pull requests were merged for the
  archive only.
- **The highest open claims** were the idle-axis chain #3934 to #3959 (up to κ 1.7367387e-3, b 1.7411726e-3). Its base
  #3251 had been retracted by its author on 18 October after exact operator replays, an open audit (#4046) showed by
  exact replay that every one of #3934's 2,374 new moves changes the operator, and the chain was retracted and closed
  on the morning of 19 October (*read on the community's repository*; nobody on my side re-ran that replay on the
  chain).
- **The best claims without a known refutation** were #3945 (1.7306632e-3, reproduced here to all 24 digits by my own
  pricer and a skeptic's) and, later in the day, claims of about 1.7310e-3 built on #3945's words (#4076, #4085,
  #4088, and #4091, a relaxation optimum rather than a literal circuit). No complex-side construction claimed there by
  others, refuted or not, was at or above 1.7501e-3 at the last read; some research notes print higher ceilings or
  targets, which are not constructions.
- **The public Lean board** (bigospeedrun.com, exact DFT) was topped by 1.2654e-3 (Hrishik Sai Bojnal, sharpening the
  constants of canghaimeng's C588AT, 1.26e-3, which is built on my Lean engine and credits it). C588AT does a partial
  version of my mutable-sources step (84 of 112 copies) and none of hosting, exit frames or the lever. Every rung of
  section 4.1 from exit frames (1.2704023e-3) up is above it.
- **κ.** The lever raises the complex side above the community's selected bit side, 1.7334270e-3 (#589). κ is capped
  both by the complex side and by a bit side, so with that bit side it would rise only to 1.7334270e-3, about +0.16 %
  over #3945's κ (under the selected record's own transfer the #589 bit side caps κ at 1.7304275e-3). Higher bit sides
  have been claimed since (#3647, ai-nanachi, 1.7366460e-3, audited by #3700, archived on 20 October; and #2423,
  thakrarsagar, 1.7348463e-3, merged on 16 October; neither is selected); they too are below the lever's complex
  side. Most of the lever's gain would not reach κ until the bit side improves; this repository does not touch the
  bit side.

Every community figure above is the author's claim checked by the community's programs, and by mine where I say so;
only my own figures are Lean-checked, and only those marked so.

### 4.10 A higher table on #2979's words with the full cut (A-fc), checked by computer, not proved in Lean

On #2979's words, unchanged, I re-solved the table of where helper storage goes for the full cut (a class solution
of a mixed-integer program, then a literal placement with the construction of #3815 by Erwan Secousses (ChaloupeMusic) in
a tolerant mode, which needed 9 spare roles). The full cut on that table prices at **1.7887046e-3**, 8.0e-6 above the Lean result on
#3945 (*checked by computer*, three routes that agree to 3e-60; the grid figure holds and one step higher fails). It
has 16,576 heads and no other rank-48 block, every paid rank is at most 23 < 25, and Theorem 1G applies (Theorem 1G and the general theorem of the full cut are
*proved in Lean*, as used on #3945 with every head kept). In this table's own instance, its heads-first certificate and
the rate fact 17,887,045/10^10 are *checked by computer*, *not proved in Lean*; the cut fact on #2979's words is
*argued by hand, not checked*, from the template it shares with #3945 (where it is *proved in Lean*). The community's
own figure for it, with no lever, is 1.6897030e-3, below #2979's 1.7290340e-3, so their admission rejects it by
design at its strict-gain steps; it passes every other check of #2979's line, run unmodified on a fresh clone of
their repository, including #2979's owner engine with its corruption controls. Files, checks and the commands to
re-create it: [data/a-fc/](data/a-fc/README.md); programs: `programs/a-fc/`. A first experiment: build it in Lean
the way `package/leverF` was built (the bridge of the full cut is general in the table, but its cut fact is proved in
Lean only on #3945's words, so it must first be ported to #2979's).

## 5. What failed or was refuted

- **"Natural exits", 1.4862904e-3 and 1.4836786e-3** (13 October): letting helpers stop wherever their last gate
  leaves them, with no restore. **Refuted** by a skeptic: the engine's undo step needs helpers restored at full frames.
  Both figures were withdrawn before anything used them; the exit frames the community charges (1.2704e-3) stand.
- **A 4,999-bank table for #588** (1.2401918e-3): **refuted** as unbuildable (the packing programme is infeasible).
- **The lever at +8.17 % on T74** (1.8169468e-3): assumed every idle block a head or tail. **Refuted** by the head
  test: only 48.4 % are heads, and none is a tail.
- **The lever at +8.7 % on #1042** (1.821393e-3): **refuted** the same way: 46 % of its idle blocks fail the test, and
  #1042 is worse than T74 under every separated variant (best 1.6918825e-3).
- **The lever at +3.79 % on T74** (1.7434773e-3), and +4.20 % on #3945 (1.8033130e-3): heads cut to 25 with the other
  rank-48 blocks kept. Not covered: the rule is not separated, a skeptic's caveat that held. It is open, not refuted
  (section 6).
- **Tails.** The pair-0 blocks (51.6 % of the idle mass on T74) cannot be used without reordering the community's
  bridges, and moving one bridge changes what the next one reads: a new construction, not a property of the current
  circuit (*open*: a guess, not tested; section 6, item 4).
- **Packing the mixed rule into one word for the old recursion** loses about 0.5/ln B to 0.75/ln B of the whole
  saving per level on the test rules of [../existence-routes/](../existence-routes/) (there the loss in the exponent
  is about 0.09/ln B); with a lever worth about 1 % of the saving that would need ln B > 50. Only a new recursion gets
  the value (*argued by hand, not checked*, from the measured drift).
- **The full cut on #3945's published fill** is at break-even (+0.006 %) and was not built.
- **The T74 target** (1.6968621e-3) was overtaken by the community by 15 October (its newest open figure that day was
  1.7173e-3); the lever was not, because the newer circuits keep T74's schedule.
- **Every margin lever inside the community's model** (section 4.3), and the five grafts of 13 October.
- **A Walsh-packet what-if of up to 1.7107451e-3** (+1.84 % on T74): **refuted** as built (wrong at the matrix
  level).
- **An apparent contradiction**: an earlier census found 0 idle blocks that could be merged into larger blocks. It
  asked a different question (merging inside the old recursion) and does not conflict with the head test.
- **The 5-copy hosting re-solve** of 13 October was cancelled while waiting for a machine; its floor 1.6277775e-3, a
  what-if that needs hosting across runs, was never re-solved or built.
- **A misreading caught by a skeptic:** one reader labelled a κ headline bound by the bit side (1.7349004e-3) as a
  complex-side figure. Corrected before use; no lever figure depended on it.
- **A "safe lower end" of heads on the T740 circuits** (idle48 - 224 T) rests on two unchecked assumptions (the same
  template, no early sibling role); it is *argued by hand, not checked*, not a safe price. It does not touch #2979 or #3945, whose heads were
  measured.

## 6. What is still open, with a first experiment for each

1. **The full cut on design A** (the largest lead on #2979). The full cut is built on design B (`package/leverF`,
   1.7806935e-3, proved in Lean): its bridge from the network certificate to the rule's numbers
   (`Work/Mixed/Bridge/CoreF.lean`) is general, and the split of a block is an identity with no condition on the
   block. Prices on #2979: 1.7797095e-3 on design A (+1.32 % above its half cut) and 1.7744822e-3 on the published
   fill, whose rate fact is already enclosed (refused one unit higher) by the rate generator. **First experiment:**
   generate design A's full-cut instance on that bridge (the instance file, the rate module at the largest grid
   numerator, the Fourier chain at its exponent) on design A's built tree, and run the release checks.
2. **Rules that are not separated** (+4.20 % on #3945, heads cut to 25 with the other rank-48 blocks whole). Theorem 1G
   uses separation to bound the exits; a different construction might not need it, and the existence-routes section
   shows that pieces with larger paid blocks can boost each other. **First experiment:** simulate the lockstep tree for
   small non-separated rules in the job model (the rule X of `cost_fit.py` is one, and runs in seconds) and look for
   left-over jobs that grow faster than the theorem allows; a counterexample closes the lead, none suggests a proof.
3. **Heads with hosted helpers in front of them** (27,648 partly hosted blocks of rank 3 and 5 on #3945's published
   fill; the half cut plus these blocks cut fully would price at 1.7600092e-3, +1.70 % over #3945's b against the half
   cut's +1.12 %, *argued by hand, not checked*; on design B only 4,272 remain, the rest became literal heads). The head lemma needs nothing before the
   head; a version that allows helper roles disjoint from the head's axes would cover them, and the Lean proof already
   has a commuting route for a block on one array. **First experiment:** for each such block, list the hosted roles in
   front of it and check that their cells are disjoint from the block's axes (seconds of Python on the head test's
   output and the owner fill).
4. **Tails.** **First experiment:** search over orders of the bridges that read pair 0 (three on T74), with the
   community's exact operator replay as the test of correctness.
5. **Exact O(n^τ) for grouped blocks**, with no log factor. Lemma A4 shows this flow cannot give it. **First
   experiment:** look for a flow whose exits carry a Φ fraction below η - β (Lemma A2), for example by feeding exits
   back into a second tree, and test it with `check_E.py`'s hypothesis check on rules A to J.
6. **The next circuit of the community.** The lever is not tied to one circuit: every newer circuit that keeps T74's
   schedule has heads in the same places, about +1.1 % to +1.6 % under the half cut. A circuit that pays idle blocks
   above rank 48 (one community note says 57 or 72 is possible under another reading of the halving guard) would break
   the half cut's separation. **First experiment:** run the head test (`heads_fill.py`) and `lever_table.py` with the
   new circuit's profile and fill; then check whether its allocation leaves head bins empty that could be filled
   elsewhere at the same price (section 4.8).
7. **κ, the bit side.** Out of this repository's scope; the lever's gain reaches κ only when the bit side rises.

## 7. The programs and how to run them

Everything here is Python 3.9 or newer with the standard library only, and runs on a laptop in under a minute and
under 300 MB. Run each with `python3 -I -B` from its own folder. The head test and the conversion programs are shipped
once, in the lever's tools ([`tools/lever/`](../../tools/lever/): `heads_fill.py`, `recount.py`, `g0_staged.py`,
`lower_b.py`, `lever_price.py`, `pricer.py`, the community's exact moment engine re-implemented, and `mixed.py`);
they are not copied here.

| folder | what it shows | command | re-run for this revision |
|---|---|---|---|
| [programs/theorem-1g/](programs/theorem-1g/) | the construction of Theorem 1G and both skeptics' checks (section 4.4) | see its README | yes, all outputs identical |
| [programs/lever-prices/](programs/lever-prices/) | the certified table of sections 4.6 and 4.8; the T74 figures and break-evens of 14 October; the gate's rigorous re-pricing | `python3 -I -B lever_table.py out.json` | yes, every stated figure reproduced |
| [programs/head-test-2979/](programs/head-test-2979/) | how the #2979 head test was assembled from the community's files (stub of helper origins, fetch, warp scan) | see its README | the head test, yes, identical |

Outputs are in [data/](data/): `theorem-1g/` (the construction's checks), `lever-prices/` (the certified table, the
14 October prices, the gate's brackets, the #3945 and #2979 price reports), `heads/` (every head test: T74 in the
community's controller and in my network, the skeptic's recount, #3945 and #2979 with their controls). The proof is
[proofs/theorem-1g.md](proofs/theorem-1g.md). A machine-readable summary is [results.json](results.json).

The community's files (words, fills, profiles) are not copied into this section, except #2979's order and dispatch
files in data/a-fc/design/ (NOTICE, section 10.i); the programs' READMEs name each one by pull request, path and
sha256.

## 8. Credits

- **The circuits.** T74: Rohan Arun (#1297, #1319), Marius-Paul Dumitrean (#1625), Jeremy Kothe and LeverRank (#1672),
  on Jeremy Kothe's k8 leaver word and seven-slot line (#1018, #1042) with Thomas DiFiore (#1148), Jeremy Diviney (#1193)
  and Chafik Boukhalfa (chafreaky, #1196); outer and assembly by Vlado Vojdanovski (#1673), Chafik Boukhalfa (#1699,
  #1786) and Jeremy Kothe (#1839). #2979: Marius-Paul Dumitrean (prepared, as his notice says, with substantial AI
  assistance from Anthropic's Claude), on Jeremy Kothe's #2957, #2948 and #2966 and Mason Hale's #2917, with Jeremy
  Diviney (#2908), Dugongue (#2969, #2863, #2839), Yiyang Lu (#2857, #2848), Ross Richardson (#2918), Tekalign,
  icekylinx, themilkman216 and Gareth Davidson (#1995, #2013); its notice keeps every original notice of its
  predecessors in force. #3945: Dugongue and the lineage its notice names (dataisfire #3734, zxv #3682, Ross Richardson
  #3605, ChaloupeMusic #3815, and the others quoted in the lever package's NOTICE). The borrowed helpers and their warps:
  Cinna Davis (#2226) and Jeremy Kothe (#2312). #1042 and #588: Jeremy Kothe.
- **The ladder.** Ryan Shea (the Lean instantiation of #588, its generators and bank and seam pipeline), Jeremy Kothe
  (#588, #799, and the U8 word of #869 and #870), Chafik Boukhalfa (#890, #741), Nil Mamano (the gauged interface, the
  nine-stage network and the record files), Dugongue (#401, the entrance gauges and banks), Zhihao Chen /
  jacklightChen (#377, the compressed nine-stage layout; #425, the argument for mutable sources), DaysSky (#442,
  hosting), DreamingOfClouds (#509, #589), eumemic (#555), the community's endpoint lemmas (#593 by Jeremy Kothe, #599
  by DreamingOfClouds), and OpenAI for the model and the Fourier reduction. Each package's NOTICE gives the exact files.
- **The community's moment engine**, re-implemented (not imported) by the pricing programs: Rohan Arun's native complex
  moment package on Dugongue's #552, eumemic's #543 and #558 and icekylinx's notes, with eumemic's outer assembly.
- **The frontier.** The audit #4046 by ChaloupeMusic; the idle-axis chain by Rohan Arun; the public Lean board entries
  of Hrishik Sai Bojnal and canghaimeng (Qiuyu Chen), built on my engine.
- **The ideas.** The τ-theorem is Schönhage's; the laser method is Coppersmith and Winograd's. Theorem 1G, the head
  test, the half cut and the every-head-bin-empty fills are mine, made with the agents.

## Sources

Everything below is on my own machines unless it names the community's repository; the agents' working folders are
not part of this repository. Paths are relative to my verification folder.

- The ladder: `wht43/bound-check/LADDER.md` (13 Oct; rows R0-R4 priced, later built), `NOTE-round2.md` (the grafts,
  the natural exits, the hosting what-ifs), the build folders `m-gen`, `k-combine`, `h-L3-build`, `fs-build`,
  `fsp-cc`, `fsp-build` and their skeptics `v-c-pack`, `v-m-statement`, `v-m-rebuild`, `v-k-statement`, `v-k-rebuild`,
  `v-host-0..2`, `v-fs-statement`, `v-fs-rebuild`, `v-fsp-rebuild`; the anatomy dissection `wht43/anatomy/` (the factor
  1.318). Private catalogue entries R1-R7, R22, R31, C6.
- A faithful port only ties: `wht43/bound-check/bh-plan/PLAN.md`, `bh-gap`, `bh-extras`, skeptic `bh-v-plan`
  (catalogue R21). The dead levers: `hm-gate/GATE.md`, `hm-walsh`, `hm-copies`, `hm-fills`, `hm-generate` and their
  skeptics `v-hm-walsh`, `hm-copies-skeptic`, `v-hm-fills`, `v-hm-generate` (R24); the grafts `g-p112-decoder`,
  `g-x4-targets`, `g-x4-slotreuse`, `g-x4-order`, `g-bit-puncture` (R25).
- The mathematics: `wht43/meta/ex-laser/PROOF.md` (Theorem 1, Theorem 3, Lemma 6; shipped in existence-routes),
  `wht43/bound-check/rl-prove/PROOF.md` (Theorem 1G; here as `proofs/theorem-1g.md`), skeptics `rl-v-logic` and
  `rl-vprove-compute` (verdicts recorded in my handoff log, 14 Oct about 15:20 EDT), the plan `rl-plan/PLAN.md`, the gate
  folders `lg-g0`, `lg-g1` (the half cut and the 24 + 24 split), `lg-g2` (the design note whose exit-decay argument
  removed the log factor) (catalogue M18, R23). The Lean form: round `wht45/lever` tracks `e-r12` (flow and cost
  bound), `e-r3prog` and `e-r3cost` (the program and its cost), `e-r45m1` (seam, rate facts, head lemma).
- Heads: `wht43/bound-check/rl-check` (controller head test, T74 and #1042), `rl-vprice` (independent recount),
  `lg-g0` (T74 in my network); `wht45/frontier/lever-reprice` (#3945 head test, recount, prices) and its skeptic
  `wht45/frontier/v-lever`; `wht45/lever/g-heads` (#3945 in my network), `wht45/lever/p-alt` (#2979 and other circuits,
  the warp scan), `wht45/lever2979/x-gate` (#2979 gate: transport, heads in my network, the community's complete replay,
  credits).
- Soundness of #3945: `wht45/lever/g-sound` (native checks of the widened terminals, complete stated verification,
  #3815 with fresh natives).
- The full cut on #2979: `wht45/lever2979/x-gate/GATE.md` section 4 and `x-inst/NOTE.md` section 4 (the rate enclosure
  N = 17744822); on design B: `wht45/allheads/bB-data/NOTE.md` section 3; on design A: `lever_table.py`, computed for
  this section.
- Every head bin empty: `wht45/reopt/FOUND.md` (designs A, B, C, written 19 Oct 14:02-14:42 EDT, with addenda to
  15:40 EDT), `o-model/MODEL.md`, `o-mip/NOTE.md` (the first optimality certificate), `o-local/NOTE.md` (the local
  search, the #3653 and #3682 fills and the second certificate), `o-param/NOTE.md` (design C), `o-eval` (the
  evaluator), skeptics `v-o-mip` (design B, 15:08 EDT), `v-o-local` (the local-search fill of #3945, 15:40 EDT) and
  `v-o-param` (design C, 15:50 EDT); the builds `wht45/allheads` (`ah-gate`: `GATE.md`, final verdict 15:58 EDT, with
  the fill checks, the operator replay and the complete admission replays; `aA-data`, `bB-data`); the 14 October proposal `wht43/bound-check/hm-gate/GATE.md` section 3 (gate G3).
- The frontier: `wht45/frontier/NOTE.md` and its readers `claims`, `selection`, `lean-frontier`, `overlap`,
  `lever-reprice` and skeptics `v-claims`, `v-lean`, `v-lever` (read 19 Oct 07:40-08:10 EDT); the re-reads in
  `wht45/lever2979/x-gate` (13:20 and 14:17 EDT) and `wht45/lever/i-final` (14:41 EDT); my own read-only listing of the
  community's newest pull requests at about 15:05 EDT (newest #4096). An early read gave #4088 as 1.7309783e-3; the
  later reads and its title give 1.7309791e-3. The idle-axis chain was retracted and closed at 13:55 UTC (09:55 EDT),
  as recorded by `wht45/lever2979/x-gate` and `x-final`.
- The community's repository: CrocSwap/integer-mult-bounds, main at 20ecf289 (19 Oct morning) and bb4bccbf (19 Oct
  afternoon); #2979 merged at be191346 (`research/rank7-joint-endframes`); #3945 at head 91f81d13
  (`research/staged544-terminal-gap-fill`); #3815 at 77ffbf28; TEMPLATE sha256 21297b20…, PLAN sha256 d566cfb8….
  The public Lean board: bigospeedrun.com/exact-dft, read 19 Oct 07:40 EDT.
- Lean status of the lever records: the lever package's VERIFY file and the skeptic notes of rounds `wht45/lever`,
  `wht45/lever2979` and `wht45/allheads` (the package folders `package/lever`, `package/lever2979`, `package/leverA`, `package/leverB` (`result.json`, `logs/`); skeptics `v-statement`, `v-rebuild`, `v-price`, `v2979-statement`, `v2979-rebuild`, `v2979-price`, `vA-statement`, `vA-rebuild`, `vA-price`, `vB-statement`, `vB-rebuild`, `vB-price`, 19 Oct 2026, 18:05 to 21:54 EDT; round `wht45/fullcut` for `package/leverF`: `fc-gate`, `fc-engine`, `fc-rate`, `fc-final`, skeptics `vF-statement`, `vF-rebuild`, `vF-price`).
- Folder names of this section: `programs/theorem-1g` holds the programs of `rl-prove`, `rl-v-logic` and
  `rl-vprove-compute`; `programs/lever-prices` those of `rl-plan`, `lg-g1` and `hm-gate`, plus `lever_table.py`
  (new); `programs/head-test-2979` those of `p-alt`.
