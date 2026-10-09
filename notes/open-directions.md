# Open directions

Author: Jacob Sussman. Text of 2026-10-09, written with the second revision of this repository. With the
third revision (2026-10-09) lead A3 was carried out: A3, the paragraph "Where things stand" and item 8 of
section C were rewritten for it. Everything else is the text of the second revision.

I expect to leave this project alone for a while: my weekly Claude limits ran out with the list below still
open. So this note hands the project to anyone who wants to continue it. It has every lead I know of that is
worth trying, large and small, with what is known, what it might give, how to start, and what was already
tried and failed, so that nobody has to try it twice.

The designs, searches and derivations below were made by the AI agents (Claude) that I directed, as the
proofs of this repository were. The words (array, frame, move, block, helper, pair) are those of the README,
sections 5 and 6: v is the number of pairs that one run of the helper circuit serves, R its number of helper
arrays, h the dimension of its frames. Savings are whole-block figures in the bridged five-stage word unless
said otherwise; figures given to eight decimals (such as 5.4143e-4) are truncated. Pull requests (#n) are
those to CrocSwap/integer-mult-bounds, whose "kappa" and "complex side" the README explains in section 2;
the full credits are in RELATED-WORK.md. Every figure carries a label:

- *proved in Lean*: checked by the kernel; "comparator" is added where the comparator passed.
- *tested design*: computed from a complete design that exists as a certificate file accepted by a Python
  reference checker; tests on arrays are named where they exist. Not in Lean.
- *what-if*: arithmetic with no design behind it.
- *outside claim*: claimed by outside authors; how they checked it is said. Not checked here unless said.
- *on paper*: derived on paper by the agents, with computer checks at small sizes. Not in Lean.

Where things stand (all four proved in Lean, comparator): for OpenAI's Fourier statement of every length
7.474546e-4, from the kernel program of the next figure; for the Walsh-Hadamard transform 7.474547e-4 with
the community's circuit of #193 at h = 22; 5.399225e-4 with my first circuit at h = 16; in per-rank
accounting 3.155781e-4, first circuit.

## A. Major leads

### A1. The scratch copies

Every run of the helper circuit copies h totals into scratch arrays, and the copies make moves of their own
that buy nothing: a unit saves D = 4v - 5 cst moves, where cst is the number of copy moves of one run. This is
the largest lead. The README, section 7, and `notes/scratch-copies.md` have the full account. In outline:

| circuit | cst | share of 4v eaten | proved in Lean | what-if: no copies at all |
|---|---|---|---|---|
| second result (community circuit, h = 22) | 440 = h(h-2) | 41.7 % | 7.474547e-4 | 1.28147e-3 (x1.71) |
| first result (my circuit, h = 16) | 241 = (h-1)^2 + h | 53.8 % | 5.399225e-4 | 1.16867e-3 (x2.16) |

- **Outside lower bound** (outside claim: #203, huxint, a draft, read at 16:25 UTC on 2026-10-09; a written
  proof in a research package, with 59,048 finite controls). For the paired ports with p >= 5, every non-zero
  combination of coordinate stars has a support that spans at least h - 2 dimensions, so a scatter through
  exactly h copied centres costs at least h(h-2). The community circuit stands exactly at that floor. In the
  author's words: "Overcomplete center factors and different center implementations remain open."
- **Floor from this project** (on paper). Copies that move along coordinate directions and are read only by
  the outputs need cst >= h(h-3) for every identity over triples. At h = 16, in the three-stage word, even
  that floor gives at most x1.25.
- **Open.** More than h copies of lower rank each (the overcomplete case); totals that are not combinations of
  stars; a hand-over of the totals without copies; one set of copies for all five stages. None has a design.
- **How to start.** `python3 tools/whatif_copies.py <certificate>` prices a block list without its copy
  blocks. Copies at other frames are data for the checker, which reads their ranks from the frames of the
  retained totals. A change of the scatter step itself probably needs a change in `Work/GCert/Chain`.

### A2. A per-rank theorem for the generalised engine

Per-rank accounting prices every single move as one recursive call on 1/m of the coordinates, as OpenAI's own
recursion does. The generalised engine has a whole-block theorem only (`BlockWHT.wht_main_of_ggroup_list`,
`Work/GFrame/Engine/Group.lean`), so the second result has no per-rank companion.

- **What exists.** The rate inequality of the community circuit at 1 - 4058344/10^10 is a Lean lemma already:
  `foldRank_B2Gp193x` in `Work/GCert/Data/RateB2Gp193x.lean` (proved in Lean; it is arithmetic about the block
  list, not yet a theorem about a program). It has the form of the numeric hypothesis of the first engine's
  per-rank theorem, `RAM.engine_program_group_perRank` in `Work/FoldRate/Group.lean`.
- **What is missing, and what it would give.** That theorem for certificates of the generalised engine (`GCert`
  in place of `BCert`): a per-rank theorem at 4.058344e-4, 1.286 times the present 3.155781e-4.
- **How to start.** By a reading of the Lean source made for this note: the first engine's proof cuts every
  block of the word into unit moves (`w.flat`) and hands the result to the engine for unit moves
  (`engine_program_scratch`). The generalised word needs the same cut for blocks that stand between two
  adapter steps (`Work/GFrame/Engine/Letters.lean`, `Adapt.lean`). The work has not been sized.

### A3. The Fourier transform of every length: done, and what remains

This lead is no longer open. The third result of this repository is OpenAI's headline statement, the discrete
Fourier transform of every length (`DFTProgram`, `TimeBounds`), with saving 7.474546e-4, and convolution
with it (proved in Lean, comparator: `transform_mainZ` and `convolution_mainZ` in `Work/Fourier/Main.lean`).
OpenAI's own statement has 10^-13.

- **How it went.** The first reading for this note had the exponent enter OpenAI's chain through one
  definition. The scoping made that exact: the chain uses the network through one theorem, `hills_program`
  in `TensorProgram.lean`, with one call site, and the generalised engine proves the same sentence at the
  exponent of the second result (`hillsZ_program`, `Work/Fourier/Seam.lean`). OpenAI's 11 files between that
  theorem and `transform_main` are repeated under `Work/Fourier/` with that one constant changed (README,
  "Whose reduction this is"). The other route named here before, through danadran01's development, was left
  aside: it gives the same exponent and would have moved the Lean of this repository into theirs.
- **Outside** (outside claims; none was built or checked here). danadran01/exact-dft-power-saving: Lean,
  saving 3.2e-6 on the same statement, the first Lean improvement of it that my scans found. Written
  transfers, not formalised: shea256/fourier-transform-below-nlogn, the first of them (6.7e-4 when
  last read, "proposed conditional transfer"), and eumemic/exact-dft-bounds (every saving below 4.856e-4 with
  an extra log log factor).

What remains of it:

| | what | size | label | what it needs |
|---|---|---|---|---|
| a | The last digit: the Fourier statement at the kernel's own 7.474547e-4 | one unit of 10^-10 | on paper | OpenAI's reduction costs a factor (log log n)^2, so the Fourier exponent has to lie strictly above the kernel's. A kernel program one unit better gives it, and that is B3: 7474548 and 7474549 hold by exact arithmetic and are not yet proved |
| b | A transfer theorem with the kernel program as a hypothesis: a kernel program at exponent z gives the Fourier statement at every exponent above z | no new figure; after it, a better network is one line for both transforms | on paper: an estimate of the scoping agent, not checked | the 11 repeated files with the kernel program passed through as a hypothesis (167 declarations; Lean universe levels are the risk). `envelope_program_ge` in `Work/Fourier/Seam.lean` is the kernel side at every exponent above the present one. The twin for the Walsh-Hadamard transform exists: `wht_main_of_engine`, `Work/Block/WHT.lean` |
| c | A per-rank companion for the Fourier statement | just below 4.058344e-4 (what-if) | what-if | A2 first; then the same seam |
| d | The transfer without padding (eumemic's batched route) with this network | below 3e-10 (the padding of the kernel program has fill parameter s = 40) | on paper: an estimate of the scoping agent | much more Lean for the ninth decimal; not recommended |
| e | Why OpenAI's Lean has the factor (log log n)^2 where their manuscript has (log log n)^(4 - theta) | no figure | not determined | a reading of `Asymptotics.lean` beside the manuscript |
| f | An independent build and comparator run of danadran01's theorem | no figure of mine; their record reports a clean build and no comparator run | not attempted | their tree builds on the same Lean and Mathlib versions as this one |

### A4. The room under the outside ceilings

Three outside notes bound what this family of designs can reach. All three are outside claims (written proofs
of their authors, evaluated by their own exact scripts). Two do not cover the design of the second result,
and the third is more than ten times above it.

- **#192** (DaysSky): kappa below 7.010e-4 for every frame layout of #168's word. It fixes the additions, and
  #193 changes them. One agent reproduced its three figures on the rebuilt data of #168 (not rechecked).
- **#201** (DaysSky; read at 16:25 UTC on 2026-10-09): complex side below 2.5657e-3 for every word on #144's
  paired-cube design with its three stages. Its list of what is not covered names "Source-assisted words
  (#184, #191, #193, #194)", so the circuit used here, and "A different decoder, port family or stage
  structure", so the five-stage word. Its author calls the ceiling "far from tight".
- **#203** (huxint, a draft): with the data calls of #193's three-stage complex word kept, the saving stays
  below 8.867e-3 even if every helper and every copy were free.

Computed here with the first lemma of #192 on the unit of the five-stage word (what-if): one helper array per
source would allow 2.54377e-3 at h = 22, and no helper at all 3.24122e-3. So a factor 3.4 lies between what is
proved and what one helper per source would allow, and no known argument closes it. The number of helper
arrays per pair decides: 12.7 in my first circuit, 8.3 in #168, 7.1 in #193. Where to look:

- **Outside the ceilings.** Another decoder than the coordinate star, another family of ports, another number
  of stages (three and seven were priced, and both are below five with the same circuit).
- **Fewer helpers.** In my first circuit 3,519 of the 7,122 helper arrays hold a copy of one value (6.3 per
  pair, against 3 in #168 and 2 in #193); a redesign of its local sums in the manner of #193 has not been
  searched. Untried first steps: the exact minimum of R at h = 6, 7, 8 by a SAT solver (CaDiCaL ships with the
  Lean toolchain; no encoding exists; known bounds at h = 8: 112 and 568), and sums in which every value is
  used at nested frames only (one change of that kind took an early circuit from 11,402 helpers to 10,666).
- **What fewer helpers would buy** (what-if; first circuit family, h = 16, block shapes of an earlier circuit):
  6.9141e-4 at 10 helper arrays per pair, 8.2949e-4 at 8, 1.03649e-3 at 6. On paper, under the frame rules of
  the first engine, at least 8/3 per pair are needed for even h >= 10; that floor does not cover general frames.

### A5. The Lean frame lemma and the checker, for the multiplication project

The community's figures belong to a theorem about integer multiplication, and their authors say what is still
assumed or unchecked: the frame theorem is "stated here but not machine-checked"
(Swapnil-jain/integer-mult-kappa), "All of #144's interfaces remain assumptions" (#168), and for #193 "There
is no globally renumbered scalar transcript of the new complex word and no full Clifford/router replay". Parts
of this repository may answer parts of that:

- the frame lemma as a statement about matrices, `GF.frame_lemma` in `Work/GFrame/Clifford/Base.lean` (proved
  in Lean), with `rep_change` and `core_eq_iff` beside it. That statement uses nothing of the RAM model;
- the complex word of #193 as one explicit numbered program, `tools/certificate/gcert1-p11-pr193.json.gz`,
  with its frames and its scalar identity, for arbitrary old contents of the helper arrays, checked by the
  kernel (proved in Lean for the rebuilt circuit; that it is their circuit rests on a Python conversion);
- the same check for their other circuits, as data. #168 (eumemic) is half done (frames and rates accepted, the
  scalar check run on 165 of 1,320 sources); #181 (chafreaky) and #186 (Dugongue) passed a Python dry run only.

What is not known: whether their obligations can be stated so that these Lean facts discharge them. Their
unit has three stages, and their theorem has a second ingredient, the bit side (in #168 with 2,200 helpers
that start off the zero frame, a rule the checker here does not have), and analytic and routing interfaces.
None of that has been looked at. A first step: put the statement of their note on general frames
(`notes/general-clifford-frames.tex` in their repository) beside `frame_lemma`.

## B. Minor leads

Everything on my first circuit (B4 to B11) is below the 7.474547e-4 of the second result. It matters for the
per-rank figure, which still comes from that circuit, and as a second, independent circuit.

| | lead | size | label | what it needs |
|---|---|---|---|---|
| B1 | Re-tune the frames of the community circuit for five stages; the published frames were chosen for three | at most about 2.7 percent (to about 7.68e-4); at most about 7 percent if the chains of the input and output arrays are re-cut too | what-if: an upper bound by one agent with the method of #192 (DaysSky), not rechecked | a frame optimiser for the five-stage price, then data only. Outside frame work to draw on: #179 (chafreaky), #186 (Dugongue), #195 (huxint), by their figures 0.005 to 0.12 percent each in their layout |
| B2 | Neighbouring sizes of the community circuit (p = 10, 12, 13; used here: p = 11) | not priced in the five-stage word. For my first circuit the step from three to five stages moved the best size from h = 18 or 20 to h = 16, where an earlier version gave 5.1284e-4 against 4.7362e-4 at h = 14 and 5.0211e-4 at h = 18 | tested design (the three first-circuit figures) | the modules at another p. For my first circuit the kernel checks at h = 18 were about twice as heavy as at h = 16 |
| B3 | The last two units of the rate lemma | 7.474549e-4: by the auditing agent's exact arithmetic 7474548 and 7474549 hold, and 7474550 fails (proved in Lean, `fold_B2Gp193x_sharp`). For the first result 5399226 holds as well | tested design: exact arithmetic on the block lists of the two published certificates; not proved, not claimed | tighter rational bounds in the generated rate lemma (`tools/gx/gxrate.py`): the Lean proof leaves a slack of 2e-12 where the exact slack is 6e-10 |
| B4 | First circuit: retain the total N_q in place of the grand total, so that all 16 copies have rank 15 (cst 241 to 240) | 5.4143e-4 (+0.28 percent); per-rank 3.1635e-4 | tested design (reference checker and a Python mirror of the three Lean checks; arrays at h = 8 and 10) | data only on the first pipeline (`tools/gen/`). Claim 5414343/10^10: the slack at 5414344 is thin |
| B5 | First circuit: slot reuse (the technique of #124, jamesyc, and #143, eumemic), 11 pairs | 5.404984e-4 (+0.11 percent); per-rank 3.159496e-4 | proved in Lean, comparator, audited by a second agent, in my working tree; the files are not among the published ones | nothing. It yields little there because the output arrays stop at frames that lie inside almost no helper's first frame. With those chains redesigned: 337 to 391 pairs, about 5.53e-4 (what-if) |
| B6 | B4 and B5 together | 5.4200e-4; per-rank 3.1672e-4 | tested design | data only on the first pipeline |
| B7 | First circuit: dead-copy recycling (after round eleven of Swapnil-jain/integer-mult-kappa). 269 later helpers are hosted on arrays that end with the same content as another array | about 5.5000e-4 (+1.9 percent); per-rank 3.2488e-4 | what-if: the hosts were counted on the certificate, no certificate was written | a generator; then data only on the first pipeline, by a reading of the checker's source |
| B8 | First circuit: degenerate frames. The sources and outputs of a cube meet at a degenerate frame of dimension 3, and the helper additions move to cheaper frames | 5.6714e-4 (+5.0 percent); on top of B6: 5.6930e-4 | tested design (reference checker for general frames; arrays at h = 8 and at h = 12, there 3,061 arrays of 4,096 entries, both directions) | data for the new chain (`Work/GCert`); the kernel run has not been made. Its per-rank figure, 3.2313e-4, needs A2 |
| B9 | First circuit: a divisor rule. Scatter coefficients 1/3 and -1/6 (or -11/26 and 1/13), no grand total, R unchanged | 5.4252e-4 | what-if | the lever of B4 without its 21 extra helpers: 150 to 300 lines in the first checker; the new checker takes 1/3 and -1/6 already |
| B10 | First circuit: one-target sinks. 512 helpers that only collect for one output write into it directly | up to 5.7072e-4 | what-if, optimistic | the redesign of the output chains named in B5: 0 of 1,001 candidates fit today |
| B11 | First circuit: helper blocks in better shapes, R unchanged | about 6.6e-4 if every helper had the cheapest shape (computed for an earlier version with 13.5 helpers per pair) | what-if | a question about the circuit: frame choice alone is used up on the first engine (section D) |
| B12 | Other community circuits in the five-stage word, as further machine checks | #168 (eumemic) 6.9620e-4; #186 (Dugongue) 7.0397e-4; #196 (chafreaky) 7.0801e-4; round eleven of Swapnil-jain 7.1508e-4. All are below #193's | what-if: prices computed from their block lists | data; #168 is half checked (A5) |
| B13 | The extension rule of the certificate format: helpers that start off the zero frame | 6.1895e-4 to 6.5715e-4 on an early form of #168's circuit; nothing in per-rank accounting; unused by #193 | what-if | a new network theorem, estimated at 2,000 to 2,500 lines (section D) |
| B14 | The constants, not the exponent: a top level that first transforms a of the k coordinates in the ordinary way would replace the factor 2^a of the README, section 2, by an additive a n | changes what the theorem says about real lengths, not the saving | on paper: one agent's remark, not checked | a new top-level theorem. A smaller table than the whole orthogonal group is not known (section D) |

## C. Checking that would strengthen the result

1. **A run on Linux**: the comparator with its real sandbox `landrun`, with the second kernel
   (`"enable_nanoda": true`) and with `leanchecker`, best directly after `lake exe cache get`, so that the
   comparator compiles every module itself (that order may need more than 16 GB). VERIFY.md, sections 5 and 8.
2. **A build of the second result from nothing** (VERIFY.md, section 3, and section 8, point 12).
3. **A human reading of the statement**: lines 271 to 331 of `Work/GCert/Data/ChallengeB2Gp193x.lean`, the
   about 60 lines that define `wht`, `WHTProgram` and `WHTTimeBoundsAt`; and a check that this challenge is
   not vacuous (one was run for an earlier package with the same text).
4. **A human reading of the new proofs**: `Work/GFrame` (54 modules, 7,628 lines) and `Work/GCert/Chain`,
   `Labels` and `Scalar` (38 modules, 6,336 lines), all written on 2026-10-09; the auditing agents read
   statements, not proofs. A possible order: `GX.wht_of_halves` (`Work/GCert/Chain/End.lean`), the engine
   (`Work/GFrame/Engine/Group.lean`), the frame lemma (`Work/GFrame/Clifford/Base.lean`), and the cost of
   the two adapter steps (`Work/GFrame/Engine/Adapt.lean`), which was argued on paper only before.
5. **Runs on arrays.** No file of the rebuilt community circuit has been run on arrays at any size (smallest
   test size: p = 6, 1,896 arrays of 4,096 entries), and by their own statement the community's frame replay
   covers p = 6 and 7 only. The in-place step of #193 stands at a degenerate frame in all 1,980 places and has
   not been run on arrays at such a frame. At full size one run is 12,052 arrays of 2^22 entries.
6. **An independent review of the conversion of the outside circuit.** The programs that rebuilt it are not in
   this repository; `tools/gx/ORIGIN.md` names the outside files and commits that they read. What agrees with
   #193: 9,412 helper arrays, 12,052 arrays per vertex of their word, rank mass 794,112, the block histogram.
7. **Smaller**: a negative control and a walk of the proof term for the per-rank companion of the first
   result; an audit of the Python generators (`tools/gx/gxgen.py`, `gxrate.py`, `gxconv.py`).
8. **For the third result.** The comparator on Linux as in item 1; a build from nothing; a reading of
   `Work/Fourier/Seam.lean`, the one new proof file; a run of `diff` on the challenge (VERIFY.md, section 1) and
   of the script that makes the copies. And for the neighbouring work: an independent build and comparator run
   of danadran01's theorem, which their record does not report.

## D. What was tried and did not work

Each with its reason, in a line or two. Every negative statement holds inside the rules it was derived under
(free additions only between arrays at the same frame, one frame per array, frames only grow). None is in Lean.

1. **Two pairs exchanged by four runs of the helper circuit, not five.** Excluded on paper, derived twice:
   k twins on one line need 2k + 1 runs. The proof treats a run as a black box.
2. **Helpers that earn a saving themselves.** Moving helper duty onto data arrays only raises m; all priced variants lose.
3. **General steps on the unchanged first engine** (such as moves between frames that are not nested). Each
   costs extra rank, and one unit of it per run costs 0.48 percent. Of 1,144 candidate sites none gains.
4. **Another dot product in place of new frames.** No single form exists for 330 of #168's helper frames.
5. **Odd sizes.** Three-stage word: 3.4717e-4 at h = 17 and 3.5976e-4 at h = 19, against 3.6317e-4 at h = 18
   (tested designs). Untried: odd h with frames chosen by parity.
6. **Halving the scratch loss inside the present identity.** At most x1.25 in the class of A1; built: x1.008.
7. **Helpers starting off the zero frame, with a merged tail** (B13). Not built: a new theorem for nothing per-rank.
8. **A seven-stage word.** 5.1160e-4 (what-if) against 5.1284e-4 (tested design) for five stages, same circuit.
9. **Merging helper blocks across stages; arrays that change their role.** No network has the array the merge
   needs; no array can change kind (exhaustive search up to five arrays).
10. **Frame choice alone, fixed circuit, first engine.** Used up: at most 0.9 percent more (general frames: B1).
11. **Slot reuse at the outside rate on my first circuit.** 11 pairs fit; the output frames block the rest (B5).
12. **The recount of #193 on my first circuit as it is.** 13 helper arrays (what-if); a redesign is open (A4).
13. **Merging duplicate helpers once frames are free.** 3 to 5 of 512: the copies rest at frames not nested.
14. **A smaller table.** The isometries of the pairs generate the whole orthogonal group (exact for h = 4 to 8).
15. **Combining the two best searched circuits by parts.** 5.1273e-4 against 5.1284e-4 (tested designs).
16. **Erasing donors, reusing emptied arrays** (#184, #191). #193 dropped both: an erased donor costs a helper.
17. **Outside.** Boolean-lattice networks (#172): refuted in its comments. "Lockstep" pairing (#132): withdrawn by its author.

## E. Practical notes

**How a better circuit becomes a theorem here.** No new proof is needed if the circuit fits the certificate
format gcert/1 (README, section 5; conditions in `tools/gx/gx.py`; example: the file of #193's circuit in
`tools/certificate/`). Try a small size first: at h = 14 the 23 Lean modules of a test circuit built in 146 s.

1. Check the file in Python. `python3 -B tools/gx/refcheck.py <file>` is the reference checker (17 to 19 s at
   h = 22). `python3 -B tools/gx/gxdry.py <file>` mirrors both Lean checks and prints the sizes, the plan of
   the kernel evaluations and the savings to expect, whole-block and per-rank.
2. Generate the Lean modules. `python3 tools/gx/regen_check.py <file> <Name> <Inst> --keep <dir>` runs
   `gxgen.py` (data, one module per kernel evaluation) and `gxrate.py` (rate lemmas, final theorem, challenge,
   solution, comparator configuration) into `<dir>`. For a new circuit its last step has nothing to compare.
3. Build the certificate modules one at a time (VERIFY.md, section 3). Measured at h = 22: the data 72 s; the
   frame table 4 evaluations of 39 to 54 s and 0.74 GB; the replay of the additions 6 evaluations of 6 to
   186 s and up to 3.8 GB; the scalar replay 2 evaluations of 122 to 180 s and 7.0 GB each; the counts 54 s.
4. Run `tools/Compare.lean` (148 to 227 s, up to 7.2 GB) and the comparator (1058 to 1814 s, 9.1 to 9.6 GB;
   for the first result, at h = 16: 171 s and 6.3 GB). State one unit less than the arithmetic allows when the
   slack of the rate inequality is thin. The first result's format goes through `tools/gen/` in the same way.

**Limits that were met.**

- The machine: an Apple M4 with 16 GB. One run of the comparator needs about 9 GB at this size, so nothing
  else can run beside it, and two kernel checks of a certificate must never run together.
- The certificate at h = 22 does not fit into one kernel evaluation. It is cut into 13 evaluations for the
  label half, 3 for the scalar half and 63 counts. The scalar replay is linear, so it is cut by sources and no
  state passes between the batches. A larger circuit needs more batches, not a new proof.
- What the checker covers: certificates without the extension rule, every scratch copy at a frame of dimension
  at least 1, whole-block accounting, and the scatter step of the present invocation theorem.
- Not in this repository: the programs that rebuilt the community circuit, the search programs of the first
  circuit, and the certificate files of B4 to B8. `tools/certificate/` has the two that the theorems rest on.

If one of these leads works out, I would be glad to see it.
