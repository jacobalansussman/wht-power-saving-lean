# Walsh-Hadamard transform in O(n (log n)^z), z = 1 - 7.474547e-4, checked in Lean

Author: Jacob Sussman. Repository: https://github.com/jacobalansussman/wht-power-saving-lean. Licence: Apache-2.0 ([LICENSE](LICENSE), [NOTICE](NOTICE)).
Text of 2026-10-09, second revision. The first text, of 2026-10-09, is commit `024f763`.

A Lean 4 proof that one fixed program computes the Walsh-Hadamard transform of every length n = 2^k in
O(n (log n)^z) operations with z < 1, in the exact-arithmetic RAM cost model of OpenAI's "Exact Fourier
transforms below n log n" ([openai/math](https://github.com/openai/math), family 130, commit `fd4aeeb`).
OpenAI's result is about the discrete Fourier transform. This repository builds on OpenAI's Lean files, applies
their method to the Walsh-Hadamard transform and proves a larger saving 1 - z for it: 7.474547e-4, where
OpenAI's own tensor engine gives 2/10^11 for the same statement.

**This is the second result of this repository.** The first, 5.399225e-4, was published here on the morning of
2026-10-09. It is still in the repository and still stands (section 1, at the end). The second is 1.384
times the first.

**This is not an OpenAI project.** OpenAI did not write, review or endorse it. The 89 files under `OAI/` are
OpenAI's, unmodified. The Lake package name `OAI` and the namespace `OAI.PowerSaving` are theirs and are kept
only because the new modules extend their development.

I built this on 8 and 9 October 2026 with a team of AI agents (Claude) that I directed (section 4): the first
result in about a day, the second in the hours after it. Every proof here is checked by Lean's kernel. The new
proofs have not yet had a line-by-line human review, and I would welcome one. Please read section 2 before
quoting the number. Credits and the relation to other work: [RELATED-WORK.md](RELATED-WORK.md).

## What is new here

- **The largest saving for this problem that I could find anywhere at the time of writing, and it is
  machine-checked:** 7.474547e-4. The largest outside figure for the same quantity (the "complex side" of the
  community's networks) was 7.009184e-4 (pull request #193; read again at 16:06 UTC), certified by scripts.
- **A five-stage "bridged" layout that, as far as my scans found at the time of writing, is in no other
  work.** The same helper circuit gives 7.009e-4 in the three-stage layout its authors use and 7.474547e-4 in
  this one.
- **As far as my scans found at the time of writing, the first Lean proof of the general frame lemma for this
  machine model**, including the cost of its "adapter" steps, which until now was argued on paper only.
- **A certificate checker with soundness proofs.** A better circuit in its format is now a data file for it,
  not a new proof.
- **As far as my scans found at the time of writing, the first complete machine check of the community's
  circuit of pull request #193**, the one with the largest complex-side saving I know of, as rebuilt here
  from their published data (section 3 says what that rests on). Its own validation was local; its author's
  words are quoted in the next section.

The helper circuit itself is not mine: it is the work of the people credited in the next section. Section 5
says in plain words what each of these items is. "My scans" are the read-only looks at public sources by the
AI agents I directed; RELATED-WORK.md says what they covered and when, the last one at 16:06 UTC on
2026-10-09.

## Whose circuit this is

The helper circuit (section 6) is the part of the network that decides the size of the saving. In the second
result it is the community's in every part:

- the paired-cube construction: icekylinx, [#144](https://github.com/CrocSwap/integer-mult-bounds/pull/144)
  (in their main branch since 2026-10-09 05:06 UTC). Its notice credits an664 (#128) for the workspace-sharing
  principle and eumemic (#117) for the producer that its modules restrict;
- the modules, carrier links, frames and physical layer at p = 11: eumemic,
  [#168](https://github.com/CrocSwap/integer-mult-bounds/pull/168) ("PR168 v4");
- the local circuit inside each cube that takes a + b and a - b from the same two arrays, and the frame flow
  that goes with it: icekylinx, [#184](https://github.com/CrocSwap/integer-mult-bounds/pull/184);
- that construction on #168's newer modules: ikeboy, [#191](https://github.com/CrocSwap/integer-mult-bounds/pull/191);
- the choice of hand-over pairs with which that step runs in place and erases nothing: ikeboy,
  [#193](https://github.com/CrocSwap/integer-mult-bounds/pull/193). In their accounting this last step is the
  whole jump from 6.6307e-4 to 7.009184e-4.

#193 states its own credits so: "icekylinx (PR184, with GPT-6 Astra and Codex assistance) for the method, the
tools and the bit supplier. eumemic (PR168 v4, with Claude assistance) for the query modules and the physical
layer. Package by Avi Eisenberg with Claude assistance." Techniques of other authors that reach the circuit
through these pull requests (slot reuse at birth, frame descent) are credited in RELATED-WORK.md, section 1.

The agents rebuilt the circuit from the files that these pull requests publish (the module files, carrier
links and frames of #168 and the pair list of #193; [NOTICE](NOTICE) names each file), with a generator of their
own, as one explicit list of additions. No outside program was run. That repository is under Apache-2.0.

What is mine: the Lean proof of the general frame lemma for this RAM model, the generalised engine, the
extended certificate checker (section 5), the bridged five-stage layout in which their circuit is placed
(section 6), and the machine check of the whole.

As far as my scans found at the time of writing, this is the first full machine check of that circuit. The
author of #193 states the
scope of that pull request's own validation in these words: "The exact lift and the contract checks establish
the local maps and the flow ledger. There is no globally renumbered scalar transcript of the new complex word
and no full Clifford/router replay." What is checked here is the circuit inside my statement. Their theorem is
about integer multiplication and has further parts that this repository does not touch (section 2, point 1).

7.474547e-4 is above the 7.009184e-4 that #193 states for the same circuit because the layout differs: their word
has three stages, mine has five. By the agents' computation the same circuit in a three-stage word gives 7.0091e-4.

## 1. The result

One fixed program, in the RAM model of OpenAI's proof, computes the Walsh-Hadamard transform of every length
n = 2^k in O(n (log2 n + 1)^z) operations with z = 1 - 7474547/10^10 = 0.9992525453, and in o(n log n).
The saving 1 - z is 7.474547e-4.

    theorem OAI.PowerSaving.WHT.wht_main_block_B2Gp193x :
        ∃ solve W, WHTProgram solve W ∧ WHTTimeBoundsAt (1 - 7474547/(10:ℝ)^10) W

Proof: module `Work.GCert.Data.SolutionB2Gp193x`. The statement with every definition it depends on: module `Work.GCert.Data.ChallengeB2Gp193x`
(`Work/GCert/Data/ChallengeB2Gp193x.lean`, 331 lines, imports only Mathlib). It is the statement of the first result
with the exponent replaced and nothing else: the theorem has a new name, and the comments of the file quote
the new names and the new exponent.

- `WHTProgram solve W`: `solve` is one program of the RAM, in the mode without general products. For every k
  and every complex vector x of length N = 2^k, its run on (k, i, x) is valid (no division by zero), returns
  exactly `wht k x` (Sylvester order: entry j is the sum over l of (-1)^(number of common 1-bits of j and l) x_l),
  costs at most W(k) units of work, and every integer it produces, and W(k) + 10(N+2), is at most (N+2)^c for
  one constant c. The imaginary unit i is an input because the RAM has no instruction that produces it.
- `WHTTimeBoundsAt β W`: W(k) = O(2^k (k+1)^β) and W(k) = o(2^k k) as k grows.

**The cost model is OpenAI's, with its caveats.** Numbers are exact complex numbers. One addition, one
subtraction or one multiplication by a prepared scalar costs one unit of work, whatever the size of the numbers.
Index arithmetic and array access are charged too. The coefficients are unrestricted. Only the integers
(indices, lengths) are bounded, by a polynomial in n. OpenAI describes its own result in these words: "the model
permits exact complex arithmetic and unrestricted coefficients, while integer values and indices remain
polynomially bounded in n" (their `lean/docs/130.md`). So the theorem says nothing about bit complexity, about
computations with bounded coefficients, or about numerical stability (section 2, point 6).

**Whose statement this is.** OpenAI's repository states and proves the discrete Fourier transform of every
length (`DFTProgram`, `TimeBounds`, stated saving 10^-13). Its family 130 has no Walsh-Hadamard statement. The RAM
cost model in my statement (challenge lines 22-270) is OpenAI's, byte for byte (their
`lean/ComparatorChallenges/UniformFourier.lean`, lines 9-257). The about 60 lines after it (271-331), which
define `wht`, `WHTProgram` and `WHTTimeBoundsAt`, are my project's, written on the pattern of OpenAI's
`DFTProgram`, and should be read. The agents first proved the statement with the exponent 1 - 2/10^11 of OpenAI's
tensor engine (`alpha` in their `TensorSaving.lean`) as a corollary of that engine (`WHTCheck/Solution.lean`).
Here only the exponent is replaced.

**The first result.** `wht_main_block_B2Ke16x`, saving 5.399225e-4 (z = 1 - 5399225/10^10), with my project's
own helper circuit and the first engine. Proof: `Work.CarrierCheck.SolutionB2Ke16x`; statement:
`Work/CarrierCheck/ChallengeB2Ke16x.lean` (331 lines). It was published on 2026-10-09 (commit `024f763`) and its
files are in this repository as they were published. So is its companion with OpenAI's own recursion,
`wht_main_rank_B2Ke16x`, the same statement with 1 - 3155781/10^10 (`Work.CarrierCheck.SolutionB2Ke16xR`). The two
figures of the first result are two accountings of the same network. Per-rank: every single move is one
recursive call on 1/m of the coordinates, as in OpenAI's proof. Whole-block: r consecutive moves of one array
are one call on r/m of the coordinates. That idea is the community's "whole-residual batching"
(RELATED-WORK.md). Its Lean proof for this RAM model is in this repository. The second result exists in
whole-block accounting only (section 2, point 3).

## 2. What is not claimed

1. **It is not a better circuit, and it does not verify the community's theorem.** The second result places
   their circuit in my layout and checks it for my statement. Nothing here is about integer multiplication.
   Their headline figure, kappa, also needs a second ingredient (the "bit side") and analytic and routing
   interfaces that their authors call assumptions. None of that is checked or used here.
2. **The outside figures are a different quantity, and they keep moving.** The titles of the pull requests to
   CrocSwap/integer-mult-bounds give kappa, the saving for integer multiplication: roughly the smaller of the
   "complex side" saving, which is the quantity of my theorem, and the bit side. For #193: complex side
   7.009184e-4, kappa 6.647872e-4.
   - At 16:06 UTC on 2026-10-09 pull requests up to #207 existed. The largest kappa in a title was 6.831905e-4
     (#207, Dugongue, open, 16:00). The largest complex-side saving stated in any of the bodies of #193 to #207
     was still 7.009184e-4: the circuit of #193, which #194, #197, #202, #204, #205, #206 and #207 keep while
     they change the bit side. The front page of the hub's main branch (`3b6b668`) named #186 as its reviewed
     result (kappa 6.61885549e-4).
   - Their authors call each figure conditional. The figures are certified by Python scripts, are not peer
     reviewed and are not Lean theorems.
   - [#192](https://github.com/CrocSwap/integer-mult-bounds/pull/192) (DaysSky) states a ceiling of 7.010e-4
     "for every frame layout of #168's word". #193 says of it: "The PR192 frame ceiling covers only PR168's
     fixed word and pairs, so it does not apply to this word." My figure is above it for one more reason, the
     five-stage layout.
3. **Whole-block accounting only.** The generalised engine has a whole-block theorem and no per-rank theorem.
   So the second result has no companion with OpenAI's own recursion. In per-rank accounting the result of this
   repository is still 3.155781e-4, from the first network.
4. **It is the Walsh-Hadamard statement, not the Fourier transform of every length.** The exponent has not
   been carried over to OpenAI's `DFTProgram` statement. For that statement an outside Lean development claims
   3.2e-6 ([danadran01/exact-dft-power-saving](https://github.com/danadran01/exact-dft-power-saving)), and an
   outside written proof, not formalised, claims every saving below 4.856e-4 with an extra log log factor
   ([eumemic/exact-dft-bounds](https://github.com/eumemic/exact-dft-bounds)). Neither was built or checked for
   this repository.
5. **It is a statement about growth, not a usable algorithm.** (c) marks numbers computed from a formula, not in
   Lean. The figures are those of the second result; the first result's are in brackets.
   - The proof works with a table of 2^a arrays indexed by all orthogonal matrices of a 110-dimensional space
     over F_2 [80-dimensional]. The table is never written down. a = 6049 (c) [3213], and 2^a has
     1821 digits (c) [968].
   - Below log2 n = m (a + 1) = 665,500 (c) [257,120] the program is the ordinary n log n algorithm.
     Nothing changes for any input whose length has fewer than 200,336 digits (c) [77,401].
   - The constant is huge. As one of the agents read in the Lean proofs, the program fills all 2^a arrays with
     the input, which puts a factor 2^a into the constant. With it the bound is below n log2 n only when
     log2 log2 n exceeds a/(1 - z), about 8e6 (c) [6e6]. A second explicit factor is 1 + 1/slack of the
     rate inequality, with slack 6e-10 (the auditing agent's arithmetic; the rational bounds of the Lean proof leave 2e-12)
     here [2e-10 to 4e-10]. Machine constants were never computed, here or by OpenAI.
   - With every constant set to 1, a gain of 1% over the ordinary algorithm needs log2 n near 2^39
     (c, heuristic) [2^45].

   In Lean: `tableExp n s = Nat.clog 2 n + s`, `threshold u a = u * (a + 1)`, m = 110, s = 40, 14,692
   arrays per unit [m = 80, s = 40, 9362], one unit per orthogonal matrix. Not in Lean: the order of that group,
   which gives a and so every (c).
6. **It is not a statement about bit complexity or about numerical computation.** One unit of work is one exact
   operation on complex numbers of any size, and the coefficients are unrestricted (section 1). OpenAI writes of
   its circuit result that it makes "no all-length, bounded-coefficient, conditioning, or bit-complexity claim".
   The last three limits hold here as well. The classical lower bound of order n log n for linear computations
   with bounded constants (Morgenstern, 1973) is not contradicted: this model does not bound the constants.

What is real is the exponent: a theorem about growth for all n in this model, machine-checked.

## 3. How it was checked, and the limits

All times are UTC. The AI agents ran every check, on one machine.

**Second result.**

- **Comparator.** The [comparator](https://github.com/leanprover/comparator), which OpenAI's repository uses
  for its own results, was run with the configuration `comparator/B2Gp193x.json`. It checks that the solution module proves
  exactly the challenge statement with identical definitions, replays the proof in the Lean kernel, and admits
  only the axioms `propext`, `Classical.choice`, `Quot.sound`, so a `sorry` or a `native_decide` is rejected.
  Its verdict "Lean default kernel accepts the solution" / "Your solution is okay!" was obtained three times:
  twice in the project's working tree, on Lean files that are byte for byte those of this repository, by the
  agent that built the proof at 15:14 (1814 s) and by an auditing agent with its own copy of the configuration
  at 15:43 (1116 s); and in the build copy of the repository at 16:13 (1058 s; next item).
- **Build in a copy of the repository.** A copy with the Lean files, `lakefile.lean`, `lean-toolchain`,
  `lake-manifest.json` and `comparator/` of this revision was built on the compiled files of the first
  publication's rebuild from scratch (Mathlib from its cache, the 194 published modules compiled from nothing).
  The 102 new modules that do not load the new certificate were compiled from source in the copy. For the 31
  modules that load it, the compiled files of the project's working build were cloned: Lake accepted the 25 data
  and check modules by their traces and compiled the 6 modules above them from source. `lake build --no-build`
  of all eight roots then reported every target up to date (9250 jobs, 15:48), with the three permitted axioms
  for the theorem and "declaration uses `sorry`" for the three challenge modules only. So the kernel evaluations
  of the certificate in those 25 modules were made by the compiler in the working tree, not in the copy, and the
  second result has not yet had a build from nothing. The comparator was then run in the copy and passed, with
  the same two verdict lines (16:13, 1058 s): it replays the whole proof in the kernel, those evaluations
  included. `tools/Compare.lean` then said `RESULT: PASS` there as well (16:15, 148 s).
  [VERIFY.md](VERIFY.md) has every command, time and memory figure, and its own list of limits.
- **Certificate.** The certificate is literal data in `Work/GCert/Data/Gen/`: the 43,036 gates of one
  invocation, each a group of single additions at one frame, and the table of its 18,248 frames. It is too large
  for one evaluation, so the Lean kernel checks it in segments (each `by decide +kernel`, axiom `propext` only):
  13 evaluations for the label half (4 for the frame table, 5 for the replay of the gates and 1 for its final
  segment, 3 for the end conditions and the frames of the copies), 3 for the scalar half (2 for the replay, 660
  sources each, and 1 for the outputs), and 63 block counts.
  Lean theorems, written by the agents and not generated from the certificate, join the segments and say that
  every certificate accepted by them gives the network at the counted price. Five files of `Work/GCert/Chain/`,
  which carry the bridged network over to the generalised engine, were made by a script from the corresponding
  files of the first result; ORIGIN.md names them.
- **Audits by separate agents** (other Claude sessions with their own scripts; not humans).
  I had one auditing agent read the statements of the new theorems before the run (the statements, not the
  proofs): it found every hypothesis of the end theorem discharged inside the proof, by a kernel evaluation on
  the certificate, by the rate inequality or by a size fact. I had a second one re-run the comparator with its
  own copy of the configuration and the second comparison script (`tools/Compare.lean`); both passed, and both
  reject the same proof against a statement whose saving is one unit larger. It walked the proof term: 38,678
  constants, among them every kernel check of this certificate, and no `sorryAx`, no compiler axiom and no
  constant of the first certificate. It scanned the 124 new modules that the proof imports for `sorry`, `axiom`,
  `native_decide`, `unsafe`, macros and similar: none, apart from one file-local notation in two proof files;
  the statement file contains one `sorry`, on purpose. It recomputed the rate arithmetic with independent
  integer code: 7474547 holds; 7474548 and 7474549 would also hold but are neither proved nor claimed; 7474550
  fails (also a Lean lemma). The packaging agent ran the generators again on the certificate and got all 36
  generated files of this repository byte-identical. Not examined by the auditing agents: the proof texts, the
  Python generators, the conversion from the community's files, and whether the statement is vacuous (it is the
  text of the earlier challenges; see the limits).

**First result** (the record of the first publication, 2026-10-09; nothing in it has changed).

- **Comparator**, with the configurations `comparator/B2Ke16x.json` and `comparator/B2Ke16xR.json`. Its verdict
  "Lean default kernel accepts the solution" / "Your solution is okay!" was obtained
  - by the agent that built the proof: whole-block theorem 09:59, per-rank companion 10:03;
  - by an auditing agent with its own configuration: whole-block theorem 10:09;
  - in a rebuild from scratch: whole-block theorem 11:32, per-rank companion 11:35.
- **Rebuild from scratch.** A copy of exactly the files of the first publication, with no build directory, was
  built from nothing: the compiled Mathlib files were downloaded from Mathlib's cache and the 194 modules were
  compiled (11 minutes). Then `tools/Compare.lean` and the comparator were run in that copy for both theorems.
  All four runs passed.
- **Certificate.** The certificate is 3.9 MB of literal data in `Work/CarrierCheck/Gen/Cr2h16.lean`: the
  addition circuit of one invocation (section 6). Three Boolean checks of it and 43 block counts are evaluated
  by the Lean kernel (`by decide +kernel`, axiom `propext` only). Lean theorems, written by the agents and not
  generated, say that every certificate accepted by the three checks gives the network.
- **Audits by separate agents.** For the whole-block theorem I had an auditing agent re-run the comparator with
  its own configuration and a second comparison script (`tools/Compare.lean`); both passed, and both reject the
  same proof against a statement whose saving is one unit larger. It walked the proof term: 36,438 constants,
  resting on the kernel checks of this certificate and on no `sorryAx` or compiler axiom. It scanned the modules
  written here that the two proofs import (102 for the whole-block theorem, and the per-rank solution module)
  for `sorry`, `axiom`, `native_decide`, `unsafe`, macros and similar (none); the two statement files contain
  one `sorry` each, on purpose. It ran the generators again and got all 21 files that they write,
  byte-identical (16 of them are in this repository). It recomputed the rate arithmetic with independent
  integer code: 5399225 holds, 5399226 would also hold but is neither proved nor claimed, and 5399227 fails
  (also a Lean lemma). The 180 older files of the proof, with the checker and its soundness proofs, belong to
  proofs examined by two earlier audits of the same kind.

**Limits. Please weigh them.**
- The new proofs are checked by the Lean kernel but have not yet had a line-by-line human review; the
  auditing agents checked statements, not proof texts. So the claim rests on the kernel and on the statement
  being the right one. I would be glad to have both read.
- The comparator ran on macOS with a stand-in for its Linux sandbox `landrun`. The second kernel (`nanoda`) was
  off, as in OpenAI's configuration of its own Fourier challenge. `leanchecker` was not run on these modules.
- In every comparator run the modules of the proof had been compiled beforehand, outside the sandbox, so the
  comparator replayed existing build products. Mathlib was not rebuilt from source, and the Lean toolchain and
  the comparator programs were the ones already installed. See VERIFY.md, section 8.
- New in this revision: the second result has not yet had the rebuild from nothing that the first result had.
  The copy of the repository in which the comparator passed took the compiled files of its 25 certificate
  modules from the project's working tree ("Build in a copy of the repository", above), and the other two
  comparator passes were in that working tree.
- New in this revision: the certificate of the second result is checked by the kernel in segments. The theorems
  that join the segments are Lean proofs and are replayed by the comparator like everything else. The frame
  lemma, the generalised engine, the checker and its segmentation were all written on 2026-10-09; they are the
  larger part of what has not yet been read by a person.
- New in this revision: the circuit is the outside authors' design, and it was rebuilt from their published
  data by Python programs of the agents. That the rebuilt list of additions is the circuit #193 describes
  rests on that Python conversion and on equal counts (9,412 helper arrays, the same histogram of block
  ranks), not on Lean. The theorem does not depend on this: the kernel checks the certificate as it is.
  The programs that rebuilt the circuit from the community's files are not in this repository;
  `tools/gx/ORIGIN.md` names the outside files and commits they read. The last step, the conversion to the
  certificate format, is `tools/gx/gxconv.py`.
- New in this revision: the second result exists in whole-block accounting only (section 2, point 3), and the
  Lean checker covers the certificate format without its extension rule (unused by this certificate) and with
  every scratch copy at a frame of dimension at least 1.
- For the per-rank companion of the first result the comparator and `tools/Compare.lean` passed (the two runs
  above), but no negative control and no walk of the proof term were run. Its files were covered by the text
  checks and the arithmetic.
- A check that the statement is not vacuous was run for an earlier package (same text, other exponent), not
  for the published ones. It showed that `wht` at lengths 2 and 4 is the Sylvester matrix and that
  W(k) = 2^k k fails.
- The Python programs that searched for the first certificate and that rebuilt the second are not part of the
  trusted base. The search programs were not audited. The kernel checks the certificates, not the programs.

**I welcome re-runs**, most of all on Linux with `landrun` and with `nanoda` switched on. See VERIFY.md.

## 4. How it was made

I made this on 8 and 9 October 2026, in about a day, by directing a team of AI agents (Claude). I began by
having them verify OpenAI's proof independently, and then kept pushing on one question: how much larger can
the saving be made while every step stays machine-checked? I set the goals, chose which ideas to pursue,
and decided at each stage what to do next. The agents wrote the Lean proofs, designed the network, found
the certificate by computer search, and audited one another's work. Several of the design ideas were found
and published independently by other people during the same hours; RELATED-WORK.md gives the public record
and the credit.

The second result came in the hours after the first was published. The first text named the next thing to
try: a Lean proof of the community's general frame lemma in this framework. I had the agents attempt it, and
it went through the same day, together with an engine and a checker that use it. From then on the best
circuit the community had published could be checked as it stands, so I chose to check theirs. The agents
rebuilt it from the published data, wrote the Lean proofs and ran the checks.

## 5. What is new in this revision, in plain words

Words. An *array* is one working copy of the data, of length n. Its *frame* is the coordinate system it is
written in at a given moment. A *move* takes an array to a larger frame, and a *block* is a run of moves of one
array; blocks are the only thing that costs, a block of rank r being one recursive call on r/m of the
coordinates. An addition between two arrays is free and is allowed only when both stand at the same frame.

- **The general frame lemma.** In the first result a frame is given by a set of mutually orthogonal unit
  directions of the label space F_2^m, and a block adds such directions. The lemma says that every subspace U of
  F_2^m can serve as a frame, and that an array goes from the frame of U to the frame of any larger subspace V
  in one block of rank dim V - dim U. "Every" includes the *degenerate* subspaces, those that contain a non-zero
  vector orthogonal to the whole subspace; they have no basis of orthogonal unit directions. The community
  states the lemma in a written note (`notes/general-clifford-frames.tex` in the CrocSwap repository), and
  every one of their figures above 7.4e-5 rests on it. `Work/GFrame` proves it in Lean for the family-130 RAM
  model (`frame_lemma`; 54 modules; by the auditing agent's reading the final proof uses it in the form
  `SMv.of_sub`, `Work/GFrame/Labels/StageBSem.lean`). Also proved there: two frames of the same subspace
  differ by a free step, and a block between two subspaces of which neither contains the other has rank
  dim U + dim V - 2 dim(U ∩ V).
- **The generalised engine.** The engine is the part of the proof that turns a schedule of blocks and free
  additions into one RAM program with a time bound. OpenAI's engine, and the first engine here, know frames of
  orthogonal unit directions. The generalised engine accepts every frame of the lemma. For that it has two more
  free steps on a single array (the "adapter" steps), a permutation of addresses and a multiplication by a
  power of i that depends on the address, and it proves that each is one linear pass, like the free steps that
  existed before. It has a whole-block theorem only.
- **The extended checker.** A certificate is now one ordered list of single additions "target += c * source",
  each with the frame at which it happens, and a table of the frames, each a subspace given by a basis. The
  checker is a Lean function with two halves. The label half follows every array through the list: its frames
  must be nested, the arrays of an addition must stand at the same frame, and every climb is counted as one
  block of its rank. The scalar half replays the additions in exact arithmetic and tests that they compose to
  the exchange the network needs, whatever old content the helper arrays hold. Theorems in `Work/GCert` say
  that a certificate accepted by both halves gives one invocation at the counted price, and they carry the
  bridged five-stage network of the first result over to the generalised engine.
- **The circuit.** The community's paired-cube circuit at p = 11 ("Whose circuit this is", above; sizes in
  section 6), written as a certificate in that format. By the agents' count 91 percent of its moves pass
  through frames that the first engine cannot use. That is why the three items above had to exist first.

## 6. The design

OpenAI's proof reduces the saving to a finite network. In the form used here: W arrays, each to be transformed
along all m directions of a label space F_2^m. One *move* transforms one array along one direction and is what
costs. Additions between arrays with equal labels are free. The plain schedule needs W m moves. A network that
needs only W m - D moves gives a positive saving, larger when D / (W m) is larger. The steps to the network of
the first result:

- **Bank pairs and the helper circuit.** A bank pair is two data arrays that exchange their contents through
  additions. The exchanges are what saves moves. One *invocation* serves v pairs and needs R *helper* arrays
  that hold the partial sums. Helpers make all m moves and save nothing, so R is the main cost. In the first
  result v = C(h,3), one pair per 3-subset of h points, and the circuit started from the outside circuit NStar3
  with carrier links.
- **Folded layout.** Earlier layouts index the pairs by a grid, which forces m = h^2. Indexing pairs and
  invocations by a group of isometries needs only m of about s h for s stages of invocations, linear in h.
- **One helper set for all stages.** The same R helpers serve one invocation in every stage. A re-framing lemma
  lets an invocation start with helpers that still hold old values.
- **Bridged five-stage word.** Every pair gets a twin, a second pair whose labels equal its own at three moments.
  Free additions between twins then exchange two pairs with five invocations in place of six. m = 5h, W = 4v + R.
- **Carrier certificates.** A wider certificate format in which the data arrays themselves take part in the
  additions, carry intermediate values and are read in place. This removes helpers. It needed a new invocation
  theorem, a checker and its soundness proof.
- **Searched addition order.** The order and grouping of the additions in the helper circuit was optimised by
  computer search. Together with the carrier format this brought R from 11,500 to 7,122.

Who published which of these ideas, and when, is in RELATED-WORK.md, sections 1 and 2.

**The second result** keeps the bridged five-stage word and replaces the helper circuit by the community's
paired-cube circuit ("Whose circuit this is", above). There the h = 2p coordinates come in p = 11 pairs. A
*port* takes one coordinate from each of three pairs, so v = 8 C(p,3). The 8 ports on the same three pairs form
a *cube*, and the circuit is 165 cubes, each with the same small local circuit on its 8 inputs, and the sums
that connect the cubes. Helper arrays that are no longer needed are handed to new users, and inside a cube two
arrays that hold a and b are turned into a + b and a - b in place.

| | h | v | m | R | W = 4v + R | D | moves per unit, in place of W m | one invocation |
|---|---|---|---|---|---|---|---|---|
| second result | 22 | 1320 | 110 | 9412 | 14,692 | 3080 | 1,613,040 of 1,616,120 | 70,169 blocks, 262,944 moves |
| first result | 16 | 560 | 80 | 7122 | 9362 | 1035 | 747,925 of 748,960 | 37,286 blocks, 130,993 moves |

The blocks of one invocation are counted by the kernel, rank by rank. For the second result the two totals
are the Lean theorems `inv_count` and `inv_moves` (`Work/GCert/Data/P193Price.lean`).

## 7. An open direction: the scratch copies

My weekly Claude limits ran out before I could chase this one, so here it is for anyone who wants it.

**What the copies are.** Every run of the helper circuit hands h totals to all the outputs at one moment (the
scatter). For that, each total is copied into a scratch array, and the copy has to make moves of its own: 22
copies of 20 moves in the circuit of the second result (cst = 440 moves per invocation), 15 copies of 15 moves
and one of 16 in my circuit of the first (cst = 241). These moves buy nothing. In the five-stage word the number of moves
saved per unit is D = 4v - 5 cst: the exchanges could save 4v, and the copies of the five invocations take
5 cst of it back.

**What they cost.** Between two fifths and a half of the saving this design could reach is eaten by the
scratch copies that every run of the helper circuit needs. A what-if computation says that removing them
entirely would raise the second result by a factor 1.71 and would roughly double the first:

| circuit | 4v | 5 cst | share eaten | saving proved | WHAT-IF: no copies |
|---|---|---|---|---|---|
| second result (community circuit, h = 22) | 5280 | 2200 | 41.7 % | 7.474547e-4 | 1.28147e-3 (x1.71) |
| first result (my circuit, h = 16) | 2240 | 1205 | 53.8 % | 5.399225e-4 | 1.16867e-3 (x2.16) |

The last column is a what-if and not a result. It is the block list of the same network with the copy blocks
and their loss deleted and everything else unchanged, priced by the rate calculator (whole-block accounting,
eight decimals, truncated; the same calculator gives 7.4745e-4 and 5.3992e-4 for the two networks as they
are). No circuit is known that reaches it. `python3 tools/whatif_copies.py tools/certificate/gcert1-p11-pr193.json.gz tools/certificate/c2-combine-best-h16.json` recomputes both rows.

**What is known.**

- One unit of copy rank costs 5 units of D. In the first network that is 0.48 percent of the saving per unit,
  as much as 41 helper arrays.
- Small steps exist. Retaining a different total in my circuit brings cst from 241 to 240 (worth 151 units of
  10^-8; found by the agents, not in Lean). The community's star rule, h copies of h - 2 moves, is already
  below the (h-1)^2 + h of my circuit.
- An earlier round of this project settled one class negatively. If the copies move along coordinate
  directions and are read only by the outputs, every identity over triples (the kind my circuit uses) needs cst
  of at least h(h-3), and the identity of my circuit needs at least h(h-1). Even the floor h(h-3) would give at
  most a factor 1.25 at h = 16 (computed for the three-stage word). Halving the loss is excluded in that class.
  This is a derivation on paper with computer checks at small sizes, and one step quotes an outside theorem.
  It is not a Lean theorem.
- That finding predates the general frame lemma in this repository. The copies may now stand at any subspace,
  degenerate ones included. What that allows for the copies has not been examined in this project.

No design for that is known yet. If you find one, the checker and pipeline in this repository will take it as
data: a circuit in the certificate format of section 5 is a data file. The note `notes/scratch-copies.md` has the
formulas and says where each of the four points above comes from.

### The other open directions

I expect to leave this project alone for a while, so [notes/open-directions.md](notes/open-directions.md) hands
it on: every lead I know of, with what is known, what it might give, how to start, and what already failed.

- **Major leads.** The scratch copies above (what-if 1.28147e-3, no design known); a per-rank theorem for the
  generalised engine (its rate inequality at 4.058344e-4 is a Lean lemma already); the Fourier transform of
  every length, which OpenAI headlines (largest outside Lean claim that my scans found: 3.2e-6); the room under
  the outside ceilings (#201, DaysSky, claims 2.5657e-3 and covers neither this circuit nor this layout); the
  Lean frame lemma and the checker, for what the multiplication project still assumes.
- **Minor leads** (fourteen, each with its computed size and its label), **the checks I would welcome most**,
  **seventeen things that were tried and did not work**, and **how a better circuit becomes a theorem here**.

## 8. Layout

    OAI/          89 files of openai/math (lean/OAI/Computability), unmodified
    WHTCheck/     the first Walsh-Hadamard corollary; defines `wht` and `WHTProgram`
    Work/         237 modules written here.
                  First result: recursion engines and rate arithmetic (Scratch, Block*, FoldRate), invocation
                  and network theorems (SharedSumStructured, Combine, Reframe, Bridge, BridgeGeom, Carrier),
                  certificate checkers, data and final theorems (SharedSumChecker, CarrierCheck).
                  Second result: GFrame (frame lemma, generalised engine, 54 modules, 48 of them imported
                  by the final proof), GCert (extended checker in the parts Labels, Scalar and Chain; Data
                  with the certificate and the final theorem)
    comparator/   the comparator configurations
    tools/        `Compare.lean`, the certificates as JSON (`tools/certificate/`), their reference checkers,
                  the generators (`tools/gen/` first result, `tools/gx/` second), `whatif_copies.py` (section 7)
    lakefile.lean, lean-toolchain, lake-manifest.json    Lean v4.34.1 and the Mathlib commit pinned by openai/math
    VERIFY.md     how to rebuild everything and re-run every check
    ORIGIN.md     where every file comes from (with MANIFEST.sha256); `tools/gx/ORIGIN.md` for the tools and
                  the certificate of the second result
    RELATED-WORK.md, NOTICE, LICENSE
    notes/, third-party/    the notes of section 7 (scratch copies; all open directions); the NOTICE file of CrocSwap/integer-mult-bounds (NOTICE, section 6)

The comments inside the Lean files were not edited for this release, so that the files are byte for byte the
ones that were checked. They still use path names and working titles of the tree they were written in ("Scratch
file", `checks/...`, "(key: ...)"). ORIGIN.md translates them.
