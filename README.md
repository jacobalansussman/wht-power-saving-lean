# Walsh-Hadamard transform in O(n (log n)^z), z = 1 - 5.399225e-4, checked in Lean

Author: Jacob Sussman. Repository: https://github.com/jacobalansussman/wht-power-saving-lean. Licence: Apache-2.0 ([LICENSE](LICENSE), [NOTICE](NOTICE)).
Text of 2026-10-09.

A Lean 4 proof that one fixed program computes the Walsh-Hadamard transform of every length n = 2^k in
O(n (log n)^z) operations with z < 1, in the exact-arithmetic RAM cost model of OpenAI's "Exact Fourier
transforms below n log n" ([openai/math](https://github.com/openai/math), family 130, commit `fd4aeeb`).
OpenAI's result is about the discrete Fourier transform. This repository builds on OpenAI's Lean files, applies
their method to the Walsh-Hadamard transform and proves a larger saving 1 - z for it: 5.399225e-4, where
OpenAI's own tensor engine gives 2/10^11 for the same statement.

**This is not an OpenAI project.** OpenAI did not write, review or endorse it. The 89 files under `OAI/` are
OpenAI's, unmodified. The Lake package name `OAI` and the namespace `OAI.PowerSaving` are theirs and are kept
only because the new modules extend their development.

I built this in about a day, on 8 and 9 October 2026, with a team of AI agents (Claude) that I directed
(section 4). Every proof here is checked by Lean's kernel. The new proofs have not yet had a line-by-line
human review, and I would welcome one. Please read section 2 before quoting the number. Credits and the
relation to other work: [RELATED-WORK.md](RELATED-WORK.md).

## 1. The result

One fixed program, in the RAM model of OpenAI's proof, computes the Walsh-Hadamard transform of every length
n = 2^k in O(n (log2 n + 1)^z) operations with z = 1 - 5399225/10^10 = 0.9994600775, and in o(n log n).
The saving 1 - z is 5.399225e-4.

    theorem OAI.PowerSaving.WHT.wht_main_block_B2Ke16x :
        ∃ solve W, WHTProgram solve W ∧ WHTTimeBoundsAt (1 - 5399225/(10:ℝ)^10) W

Proof: module `Work.CarrierCheck.SolutionB2Ke16x`. The statement with every definition it depends on:
`Work/CarrierCheck/ChallengeB2Ke16x.lean` (331 lines, imports only Mathlib).

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
computations with bounded coefficients, or about numerical stability (section 2, point 4).

**Whose statement this is.** OpenAI's repository states and proves the discrete Fourier transform of every
length (`DFTProgram`, `TimeBounds`, stated saving 10^-13). Its family 130 has no Walsh-Hadamard statement. The RAM
cost model in my statement (challenge lines 22-270) is OpenAI's, byte for byte (their
`lean/ComparatorChallenges/UniformFourier.lean`, lines 9-257). The about 60 lines after it (271-331), which
define `wht`, `WHTProgram` and `WHTTimeBoundsAt`, are my project's, written on the pattern of OpenAI's
`DFTProgram`, and should be read. The agents first proved the statement with the exponent 1 - 2/10^11 of OpenAI's
tensor engine (`alpha` in their `TensorSaving.lean`) as a corollary of that engine (`WHTCheck/Solution.lean`).
Here only the exponent is replaced.

**Companion with OpenAI's own recursion:** `wht_main_rank_B2Ke16x`, the same statement with 1 - 3155781/10^10
(`Work.CarrierCheck.SolutionB2Ke16xR`). The two figures are two accountings of the same network. Per-rank:
every single move is one recursive call on 1/m of the coordinates, as in OpenAI's proof. Whole-block (the
headline): r consecutive moves of one array are one call on r/m of the coordinates. That idea is the community's
"whole-residual batching" (RELATED-WORK.md). Its Lean proof for this RAM model is in this repository.

## 2. What is not claimed

1. **It is not the largest saving claimed anywhere.** Larger figures for the same kind of network have been
   posted as pull requests to [CrocSwap/integer-mult-bounds](https://github.com/CrocSwap/integer-mult-bounds), a
   community project on integer multiplication. Their authors call each of them conditional. They are certified
   by Python scripts, are not peer reviewed and are not Lean theorems.
   - The comparable figure (their "complex side" saving) was 6.5634843e-4 in the newest pull request that the
     agents read in full ([#178](https://github.com/CrocSwap/integer-mult-bounds/pull/178), a draft, 2026-10-09
     10:29 UTC). Pull-request titles show another number, kappa (RELATED-WORK.md, rules).
   - At 12:57 UTC the project's main branch named one of them,
     [#186](https://github.com/CrocSwap/integer-mult-bounds/pull/186), as its current reviewed result, with
     kappa = 6.61885549e-4 (its README, which also says that finite replay and written review "do not formally
     verify the multiplication theorem").
   - By 13:12 UTC pull requests up to #195 existed, with kappa up to 6.674e-4 in a title (#194). Of
     #179 to #195 the agents read only the titles. The figures were still rising when this was written.
   - Such figures passed mine at 06:40 UTC, three hours before the theorem here passed its check.
   - They rest on a "general frame lemma" that the Lean framework here does not have, and on further rules that
     its certificate checker does not have. A first look by the agents on 2026-10-09 (derivations on paper and
     tests on small arrays; no Lean proof; not part of this repository) supports the lemma itself and found
     that it does not fit this framework as it is. Whether it can be proved here is the next thing I plan to
     try. RELATED-WORK.md, section 4.
2. **It is the Walsh-Hadamard statement, not the Fourier transform of every length.** The exponent has not
   been carried over to OpenAI's `DFTProgram` statement. For that statement an outside Lean development claims
   3.2e-6 ([danadran01/exact-dft-power-saving](https://github.com/danadran01/exact-dft-power-saving)), and an
   outside written proof, not formalised, claims every saving below 4.856e-4 with an extra log log factor
   ([eumemic/exact-dft-bounds](https://github.com/eumemic/exact-dft-bounds)). Neither was built or checked for
   this repository.
3. **It is a statement about growth, not a usable algorithm.** (c) marks numbers computed from a formula, not in Lean.
   - The proof works with a table of 2^a arrays indexed by all orthogonal matrices of an 80-dimensional space
     over F_2. The table is never written down. a = 3213 (c), and 2^a has 968 digits (c).
   - Below log2 n = m (a + 1) = 257,120 (c) the program is the ordinary n log n algorithm. Nothing changes for
     any input whose length has fewer than 77,401 digits (c).
   - The constant is huge. As one of the agents read in the Lean proofs, the program fills all 2^a arrays with
     the input, which puts a factor 2^a into the constant. With it the bound is below n log2 n only when
     log2 log2 n exceeds a/(1 - z), about 6e6 (c). A second explicit factor is 1 + 1/slack of the rate
     inequality, with slack 2e-10 to 4e-10 here. Machine constants were never computed, here or by OpenAI.
   - With every constant set to 1, a gain of 1% over the ordinary algorithm needs log2 n near 2^45 (c, heuristic).

   In Lean: `tableExp n s = Nat.clog 2 n + s`, `threshold u a = u * (a + 1)`, m = 80, s = 40, 9362 arrays per
   unit, one unit per orthogonal matrix. Not in Lean: the order of that group, which gives a and so every (c).
4. **It is not a statement about bit complexity or about numerical computation.** One unit of work is one exact
   operation on complex numbers of any size, and the coefficients are unrestricted (section 1). OpenAI writes of
   its circuit result that it makes "no all-length, bounded-coefficient, conditioning, or bit-complexity claim".
   The last three limits hold here as well. The classical lower bound of order n log n for linear computations
   with bounded constants (Morgenstern, 1973) is not contradicted: this model does not bound the constants.

What is real is the exponent: a theorem about growth for all n in this model, machine-checked.

## 3. How it was checked, and the limits

All times are UTC on 2026-10-09. The AI agents ran every check, on one machine.

- **Comparator.** The [comparator](https://github.com/leanprover/comparator), which OpenAI's repository uses
  for its own results, was run with the configurations in `comparator/`. It checks that the solution module
  proves exactly the challenge statement with identical definitions, replays the proof in the Lean kernel, and
  admits only the axioms `propext`, `Classical.choice`, `Quot.sound`, so a `sorry` or a `native_decide` is
  rejected. Its verdict "Lean default kernel accepts the solution" / "Your solution is okay!" was obtained
  - by the agent that built the proof: headline theorem 09:59, per-rank companion 10:03;
  - by an auditing agent with its own configuration: headline theorem 10:09;
  - in the rebuild from scratch described next: headline theorem 11:32, per-rank companion 11:35.
- **Rebuild from scratch.** A copy of exactly the files of this repository, with no build directory, was built
  from nothing: the compiled Mathlib files were downloaded from Mathlib's cache and the 194 modules were
  compiled (11 minutes). Then `tools/Compare.lean` and the comparator were run in that copy for both theorems.
  All four runs passed. [VERIFY.md](VERIFY.md) has every command, time and memory figure, and its own list of
  limits.
- **Certificate.** The certificate is 3.9 MB of literal data in `Work/CarrierCheck/Gen/Cr2h16.lean`: the
  addition circuit of one invocation (section 5). Three Boolean checks of it and 43 block counts are evaluated
  by the Lean kernel (`by decide +kernel`, axiom `propext` only). Lean theorems, written by the agents and not
  generated, say that every certificate accepted by the three checks gives the network.
- **Audits by separate agents** (other Claude sessions with their own scripts; not humans). For the headline
  theorem I had an auditing agent re-run the comparator with its own configuration and a second comparison
  script (`tools/Compare.lean`); both passed, and both reject the same proof against a statement whose saving is
  one unit larger. It walked the proof term: 36,438 constants, resting on the kernel checks of this certificate
  and on no `sorryAx` or compiler axiom. It scanned the modules written here that the two proofs import (102
  for the headline theorem, and the per-rank solution module) for `sorry`, `axiom`, `native_decide`, `unsafe`,
  macros and similar (none); the two statement files contain one `sorry` each, on purpose. It ran the
  generators again and got all 21 files that they write, byte-identical (16 of them are in this repository).
  It recomputed the rate arithmetic with independent integer code: 5399225 holds, 5399226 would also hold but
  is neither proved nor claimed, and 5399227 fails (also a Lean lemma). The 180 older files of the proof, with
  the checker and its soundness proofs, belong to proofs examined by two earlier audits of the same kind.

**Limits. Please weigh them.**
- The new proofs are checked by the Lean kernel but have not yet had a line-by-line human review; the
  auditing agents checked statements, not proof texts. So the claim rests on the kernel and on the statement
  being the right one. I would be glad to have both read.
- The comparator ran on macOS with a stand-in for its Linux sandbox `landrun`. The second kernel (`nanoda`) was
  off, as in OpenAI's configuration of its own Fourier challenge. `leanchecker` was not run on these modules.
- In every comparator run the modules had been compiled beforehand, outside the sandbox, so the comparator
  replayed existing build products. Mathlib was not rebuilt from source, and the Lean toolchain and the
  comparator programs were the ones already installed. See VERIFY.md, section 8.
- For the per-rank companion the comparator and `tools/Compare.lean` passed (the two runs above), but no
  negative control and no walk of the proof term were run. Its files were covered by the text checks and the
  arithmetic.
- A check that the statement is not vacuous was run for the previous package (same text, other exponent), not
  for this one. It showed that `wht` at lengths 2 and 4 is the Sylvester matrix and that W(k) = 2^k k fails.
- The Python programs that searched for the certificate are not part of the trusted base and were not audited.
  The kernel checks the certificate, not the search.

**I welcome re-runs**, most of all on Linux with `landrun` and with `nanoda` switched on. See VERIFY.md.

## 4. How it was made

I made this on 8 and 9 October 2026, in about a day, by directing a team of AI agents (Claude). I began by
having them verify OpenAI's proof independently, and then kept pushing on one question: how much larger can
the saving be made while every step stays machine-checked? I set the goals, chose which ideas to pursue,
and decided at each stage what to do next. The agents wrote the Lean proofs, designed the network, found
the certificate by computer search, and audited one another's work. Several of the design ideas were found
and published independently by other people during the same hours; RELATED-WORK.md gives the public record
and the credit.

## 5. The design

OpenAI's proof reduces the saving to a finite network. In the form used here: W arrays, each to be transformed
along all m directions of a label space F_2^m. One *move* transforms one array along one direction and is what
costs. Additions between arrays with equal labels are free. The plain schedule needs W m moves. A network that
needs only W m - D moves gives a positive saving, larger when D / (W m) is larger. The steps to this network:

- **Bank pairs and the helper circuit.** A bank pair is two data arrays that exchange their contents through
  additions. The exchanges are what saves moves. One *invocation* serves v = C(h,3) pairs, one per 3-subset of
  h points, and needs R *helper* arrays that hold the partial sums. Helpers make all m moves and save nothing,
  so R is the main cost. The circuit started from the outside circuit NStar3 with carrier links.
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

Final sizes: h = 16, v = 560, m = 80, R = 7122, W = 4v + R = 9362 arrays per unit, D = 1035, so 747,925 moves
per unit in place of 748,960. One invocation is 37,286 blocks and 130,993 moves, counted by the kernel.

## 6. Layout

    OAI/          89 files of openai/math (lean/OAI/Computability), unmodified
    WHTCheck/     the first Walsh-Hadamard corollary; defines `wht` and `WHTProgram`
    Work/         104 modules written here: recursion engines and rate arithmetic (Scratch, Block*, FoldRate),
                  invocation and network theorems (SharedSumStructured, Combine, Reframe, Bridge, BridgeGeom,
                  Carrier), certificate checkers, data and final theorems (SharedSumChecker, CarrierCheck)
    comparator/   the two comparator configurations
    tools/        `Compare.lean`, the certificate as JSON, its reference checker, the generators
    lakefile.lean, lean-toolchain, lake-manifest.json    Lean v4.34.1 and the Mathlib commit pinned by openai/math
    VERIFY.md     how to rebuild everything and re-run every check
    ORIGIN.md     where every file comes from (with MANIFEST.sha256)
    RELATED-WORK.md, NOTICE, LICENSE

The comments inside the Lean files were not edited for this release, so that the files are byte for byte the
ones that were checked. They still use path names and working titles of the tree they were written in ("Scratch
file", `checks/...`). ORIGIN.md translates them.
