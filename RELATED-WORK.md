# Related work and credits

Author of this repository: Jacob Sussman. Repository: https://github.com/jacobalansussman/wht-power-saving-lean.

This file says what this work takes from others, what it shares with work that others published first, and
where the other results stand. It rests on three read-only looks at public sources by the AI agents I directed:
a scan on 2026-10-08 at about 21:10 UTC, a scan on 2026-10-09 ending 10:39 UTC, and a re-read of 20 pull requests
and 4 repositories on 2026-10-09 from 11:46 to 11:55 UTC. **All times are UTC.**

Rules of this file:
- Every outside figure is **a claim of its authors**. How they checked it is stated next to it. The authors of
  the pull requests quoted here label each of their figures "conditional".
- **Which number is quoted.** The titles of the pull requests to CrocSwap/integer-mult-bounds give kappa, the
  saving for integer multiplication. The figure that compares with mine is the "complex side" saving, which is
  in the body or the notes of a pull request. This file quotes the complex-side saving, and writes "kappa"
  where only a title was read. Example: the title of #178 says kappa = 6.5592e-4, and its body gives the
  complex-side saving 13126968687/(2 * 10^13) = 6.5634843e-4. In the five pull requests where the agents compared
  the two (#130, #144, #155, #168, #178) the complex-side saving is a little above the kappa of the title.
- No outside checker was run and no outside Lean development was built. The agents read sources and recomputed
  some arithmetic.
- People are named by their GitHub handles, as in the sources.
- "Not found" means: searched for and not found. It does not mean "does not exist".
- The field moved every 10 to 30 minutes while this was written. Section 5 is out of date by the time you read it.

"Saving" is 1 - z in a bound n (log n)^z. My figure is 5.399225e-4 (whole-block accounting), a Lean theorem
for the Walsh-Hadamard statement described in the README.

## 1. What this work takes from others

**OpenAI, [openai/math](https://github.com/openai/math), family 130, "Exact Fourier transforms below n log n"**
(commit `fd4aeeb`, Apache-2.0). Everything here is built on it: the theorem and the idea of the proof; the RAM
cost model, which is copied byte for byte into my statement; the recursion engine and program framework, which
this repository imports unchanged (89 files) and generalises. Their result is the discrete Fourier transform of
every length with saving 10^-13, proved in Lean. Their two comparator challenges passed when I had the agents
re-run them. The tensor engine in their Lean has saving 2/10^11. The two manuscripts of family 130 are "Finite
tensor savings and exact Fourier circuits" and "An explicit power saving for the exact discrete Fourier
transform" (OpenAI, 25 September 2026, in `preprints/` of openai/math). Section 7 reproduces their BibTeX
entries. OpenAI did not write, review or endorse this repository.

**Whole-residual batching** (the accounting behind my headline figure). Community follow-ups to OpenAI's
sibling result on integer multiplication apply the whole residual of a network edge, of rank r, as one
recursive call. The earliest that the agents found is pull request #10 to
[CrocSwap/integer-mult-bounds](https://github.com/CrocSwap/integer-mult-bounds) (icekylinx, 2026-10-08 08:26):
"Apply whole complex residuals as recursive calls". shea256's manuscript (below) credits further complex
batching to eumemic (#15). The "whole-block" recursion of this repository is that idea, proved in Lean for the
family-130 RAM model and for runs of unit moves only. The Lean development of danadran01 (section 5) formalises
the same accounting under the name "grouped recursion".

**The helper circuit NStar3 and the two-stage word.** This work takes them from round six of
[Swapnil-jain/integer-mult-kappa](https://github.com/Swapnil-jain/integer-mult-kappa) (NStar3, two stages,
copied centres), the "complex side" network of that integer-multiplication follow-up. That repository names its
own sources: it builds on OpenAI's "Integer multiplication below n log n" and on the framework of the CrocSwap
repository (Douglas Colkitt); it has the copied centres "after PR #36's copied retained-centre schedule"
(CrocSwap #36, icekylinx, 2026-10-08 12:27); and it attributes its two-stage construction to "Paureel's
two-stage motif". #36 itself says that it "adopts the attributed two-stage and corner contributions from
Paureel". shea256 transferred the round-six
network to the family-130 RAM model in
[shea256/fourier-transform-below-nlogn](https://github.com/shea256/fourier-transform-below-nlogn) (commit
`bf38c00`) and proposed a saving of 7.3e-5. Checked by its author with a written argument and finite Python
checks, not in Lean. Its README describes the network as "Swapnil Jain's round-six complex network, building on
work by Douglas Colkitt, icekylinx, eumemic, and Aurel Prosz / Paureel". The agents put that network into
Lean before they changed the layout. The helper circuit in this repository still has NStar3's sums as its
starting point.

**Carrier links.** At an addition n = a + b one argument slot becomes the sum and the other rests. A link hands
the resting slot to a later use of the same value, which saves one helper array. The idea is from CrocSwap pull
requests #24, #32 and #36 (icekylinx), with a weighted version in #44 (rohanarun) and link-friendly circuits in
#53 and #62 (ikeboy). Checked by their authors in Python and C++. The agents applied it to NStar3 with a maximum
matching.

## 2. Ideas developed here that others published first

Four of the design steps in the README were worked out inside my project. Its first scan (2026-10-08, about
21:10) found none of them in public, and the project did not look outside again until the second scan. In that
interval others published the work below. Nothing of mine was published. The outside work is therefore
independent of mine, and **priority in the public record is theirs.** I claim none.

| idea (README, section 5) | public record | how close |
|---|---|---|
| Layout on a group of isometries with label dimension linear in h; a third stage | CrocSwap [#130](https://github.com/CrocSwap/integer-mult-bounds/pull/130) (icekylinx), created 2026-10-09 02:29:02, merged 05:06:02: "three invocation banks on a regular Cayley cover in dimension `3h-2`" | same idea, three stages |
| One helper set for all stages | [#137](https://github.com/CrocSwap/integer-mult-bounds/pull/137) (eumemic), 03:18:33, open; [#144](https://github.com/CrocSwap/integer-mult-bounds/pull/144) (icekylinx), 04:19:35, merged 05:06. Inside one stage: #128 (an664), 01:35:20 | same idea. In theirs each core finishes its clean-up before the next stage starts |
| Data arrays as carriers: sources used in place, pair labels, the cube of eight ports | #144 (icekylinx): "1,760 ports grouped into eight-port cubes"; "This adds no K auxiliary bank". The construction is in its note `notes/paired-cube-construction.tex` (CrocSwap main `d1d6c07`): "Ports with the same three pair labels form an eight-port cube". Cube producers appear earlier, in #115 (icekylinx, 00:26:56) | same idea on the input side. Output side related: terminal sinks, #166 (jamesyc), used in #168 |
| Searching over the order of the additions | #111 (rohanarun), 00:04; #167 (GamingPuzzled), 08:15; [#168](https://github.com/CrocSwap/integer-mult-bounds/pull/168) (eumemic), 08:22 | related: the same lever on other circuits, none of them NStar3. An earlier outside test of other tree shapes on NStar3 itself reported no gain (sobakadog8, 2026-10-08 20:13) |

The timestamps that I can show for my project are local ones. The private note that first describes the folded
layout, with one helper set serving three stages, has the file time 2026-10-08 23:32:10. Pull request #130 was
created on 2026-10-09 at 02:29:02. My project's first Lean theorem with a three-stage layout passed the
comparator at 02:38 on 2026-10-09, nine minutes after #130 appeared; the one with a shared helper set passed at
02:47. A local file time or log time is not a public record, and nobody outside can verify it.

## 3. What the scans did not find elsewhere

- **The bridged word.** A twin for every pair on the same line, so that two pairs are exchanged with five
  invocations in place of six, and the five-stage layout with m = 5h. Not found by any of the ten scanning agents.
  The nearest outside relative is "lockstep" pairing (#132, ikeboy), which its author withdrew as invalid.
- **A Lean proof at this exponent.** The largest saving in an outside Lean theorem that the agents found is 3.2e-6
  (section 5, row D).

Both statements are as of 10:39 on 2026-10-09 and are limited by what the scans covered (section 6). The
re-read at 11:46 looked at the titles of the newer pull requests only and found no new commit in the
repositories of Swapnil-jain, eumemic, danadran01 and shea256.

## 4. Why the larger outside figures are not like for like

Python-certified outside figures are larger than mine (section 5, row A). The comparison needs care.

- Theirs are per-vertex histograms certified by Python programs inside the integer-multiplication project. Nobody
  has turned the figures above mine into a Fourier theorem, and none is in Lean. (Row C does it on paper for
  the smaller figure of row B.) Mine is a Lean theorem of a transform statement.
- **Their counts rest on a frame lemma for arbitrary binary subspaces, "including degenerate subspaces"**
  (CrocSwap `notes/general-clifford-frames.tex`), with one recursive child per nested step. Every outside figure
  above 7.4e-5 depends on it. The Lean framework of this repository does not have that lemma.
- **A first look at that lemma.** On 2026-10-09, from about 11:10 to 12:10, the agents examined it: three
  reports, three counter-checks and a summary, by derivations on paper and tests on small arrays. The lemma has
  not been proved in Lean, and none of this work is in this repository. What they report:
  - The lemma itself held up, on paper and in tests on small arrays.
  - The recursion engine of this repository does not need unit directions. The layers above it do: its label
    calculus and its certificate checker work with unit directions and orthogonal projectors.
  - The lemma does not fit those layers as they are. For degenerate subspaces the outside frames are not in the
    family of frames that this framework uses. Using them needs new free operations (a permutation of addresses
    and a phase that depends on the address), whose cost in the RAM model is argued on paper only.
  - A substitute that stays inside that family of frames (another dot product) was tried on the frames
    published with #168. It works for one part of that circuit. With one dot product for the whole circuit
    it cannot work for the helper frames.
  - The lemma is necessary for the outside figures and not sufficient. The 6.56e-4 design also uses four further
    rules that the engine and checker here do not have, among them entrance gauges and slot reuse at birth. By
    the agents' count, 13,416 of the 15,149 operation frames published with #168 are degenerate subspaces.

  Whether the lemma and those rules can be proved in this framework is open. It is the next thing I plan to try.
- Their stock of arrays is indexed by the whole orthogonal group O(66,2), described in writing and never built.
  Mine is indexed by O(80,2) in the same way: defined in Lean, never enumerated.
- Outside figures use whole-residual accounting. They compare with my 5.399225e-4, not with my per-rank
  3.155781e-4. The agents found no outside per-rank figure.
- In the other direction: rows C and D below are about the Fourier transform of every length. Mine is about the
  Walsh-Hadamard transform.

Outside techniques that this work does not use: slot reuse at birth (#124, #143), operation-frame descent in
bundles (#131, #168, #178), partial and pair-aware source gauges (#115, #144), extended carrier matching (#162),
merged output reads (#161), the paired-cube producer at p = 11, and the transfer to every length without
padding to a power of two (eumemic, danadran01).

## 5. Where the other results stood on 2026-10-09

**Read in full, state at 10:39.** The column "saving claimed" is the complex-side saving (see the rules).

| | saving claimed | what and by whom | how its authors checked it |
|---|---|---|---|
| A | 6.5634843e-4 | CrocSwap [#178](https://github.com/CrocSwap/integer-mult-bounds/pull/178) (rohanarun), a draft, created 10:28:58. Same network family ("complex side"), h = 22, m = 66. Title: kappa = 6.5592e-4 | Python checkers. Author: "Full repository verification and GitHub CI pending." CI was green at 10:39. No comments. Its interfaces "remain conditional assumptions" |
| A' | 6.5632e-4 | CrocSwap [#168](https://github.com/CrocSwap/integer-mult-bounds/pull/168) (eumemic), open, head `4a3c769`, 10:09:48. Title: kappa = 6.558894e-4 | Python gates, CI green, no review. "All of #144's interfaces remain assumptions." Two of the agents recomputed the figure from the published histogram (arithmetic only) |
| A'' | 6.1609676e-4 | [Swapnil-jain/integer-mult-kappa](https://github.com/Swapnil-jain/integer-mult-kappa), round ten, `d2f6146`, 09:38 | Python replays. Lean for the moment arithmetic only. Parts "stated here but not machine-checked" |
| B | 4.856569e-4 | CrocSwap main `d1d6c07` ([#144](https://github.com/CrocSwap/integer-mult-bounds/pull/144), merged 05:06), h = 24, m = 72. Title: kappa = 4.609169e-4 | Python, written notes, a maintainer review assisted by Codex: "no independent human peer review" |
| C | every saving below 4.856e-4, with an extra (log log n) factor | [eumemic/exact-dft-bounds](https://github.com/eumemic/exact-dft-bounds), `6f87d1a`, first commit 08:33. Fourier transform of every length, and convolution | "a paper proof with an exact finite certificate, and it is not formally verified" |
| D | 3.2e-6 | [danadran01/exact-dft-power-saving](https://github.com/danadran01/exact-dft-power-saving), `1a7b25e`, 03:48. Fourier transform of every length, on OpenAI's own `DFTProgram` and `TimeBounds` | **Lean.** Their record: a clean build, standard axioms only, `leanchecker`. A comparator challenge file is shipped; no comparator run is reported. **It was not built here** |
| E | 7.3e-5; 7.4026e-5; 5.1e-10; conditional 1e-9 | shea256; sobakadog8; teal-sea; a fork by whyihaveyou. Family 130 | written arguments or Python; the last one Lean arithmetic with open premises |
| | 5.399225e-4 | this repository, Walsh-Hadamard statement, 09:59 | Lean, comparator, audits by AI agents; limits in the README |

**Titles only, state at 11:46.** Pull requests #179 to #185 had appeared. The agents read their titles, not
their bodies, so the figures are kappa and not complex-side savings: 6.5592e-4 (#179, chafreaky, open, 10:41),
6.5592e-4 (#180, M6LI, a draft, 11:10), 6.5918e-4 (#181, chafreaky, open, 11:10), 6.6022e-4 (#182, chafreaky,
open, 11:27), 6.5664e-4 (#184, icekylinx, closed, 11:32), 6.6063e-4 (#185, rohanarun, a draft, 11:41). #183
claims no new figure. #178 was still a draft. The repositories of rows A'', C, D and of shea256 had no new commit.

**Titles only, state at 12:37.** Pull requests #186 to #193 had appeared; again the agents read titles only.
Largest kappa in a title: 6.6479e-4 (#193, ikeboy, open, 12:36); before it 6.6263e-4 (#191, ikeboy, open, 12:16)
and 6.6189e-4 (#186, Dugongue, open, 11:54). #192 (DaysSky, open, 12:28) states a ceiling, "7.010e-4 for every
frame layout of #168's word", and no new figure. The hub's main branch and the repositories of rows A'', C, D had
no new commit.

**State at 13:12, front page and titles only.** The hub's main branch moved to `27d82eb` at 12:57 and now names
#186 (Dugongue) as its current reviewed result, with kappa = 6.61885549e-4; row B above is the state before
that. Pull requests #194 and #195 had appeared: 6.6743e-4 (#194, ikeboy, open, 12:43) and 6.6265e-4 (#195,
huxint, open, 13:00), both kappa from titles. Swapnil-jain/integer-mult-kappa published round eleven at 13:03
(`1a580dc`, "conditional witness kappa = 2^-10.547"). The repositories of rows C and D had no new commit. None
of this was read beyond the front page, the titles and the commit subjects.

How the largest outside complex-side claim for this network family moved during the night (all
Python-certified claims): 7.4e-5 (first scan) -> 3.148e-4 (#130, 02:29) -> 4.139e-4 (#137, 03:18) ->
4.856569e-4 (#144, 04:19) -> 5.1135e-4 (#152, 06:05) -> 5.514156e-4 (#155, 06:40, the first above my final
figure; closed at 06:54 as superseded by #157) -> 5.6228e-4 (#157) -> 5.8857e-4 (#161, 07:23) -> 6.5632e-4
(#168, 10:09) -> 6.5635e-4 (#178, 10:29).

The agents found no lower bound or impossibility result for this problem in public sources. Only limits for
particular routes exist.

## 6. How this record was made, and its gaps

- Two scans by AI agents (the second: five scanners, five sceptics, one synthesis), read-only: public GitHub
  pages and API reads, arXiv, MathOverflow, Hacker News. Nothing of mine was sent anywhere. Nothing downloaded
  was executed. Each figure in rows A to D was read by at least two agents and fetched again at 10:35 to 10:39.
  Row E rests on the scanners' readings alone.
- A third look by one reviewing agent (11:46 to 11:55): 36 unauthenticated reads of the public GitHub API. It
  confirmed the number, author and time of 20 of the pull requests named in this file and read the titles of
  the newest ones.
- Covered: openai/math and its forks (every branch of the roughly 65 forks pushed since 2026-10-07); the
  CrocSwap pull requests up to #178 (#169 to #177 by title and summary only); the repositories named above.
- Not reached: X / Twitter, the openai.com page, a paywalled newspaper article. GitHub code search does not
  index day-old repositories. Private repositories cannot be seen.
- The medium-confidence part: whether the larger outside claims are correct. The agents checked their arithmetic
  only.
- If you know of work that should be credited here and is not, please open an issue at https://github.com/jacobalansussman/wht-power-saving-lean.

## 7. How to cite OpenAI's manuscripts

openai/math asks that its manuscripts be cited with the BibTeX entries of their directories. The two entries
of family 130, copied from `preprints/` of openai/math at commit `fd4aeeb`:

    @misc{OAI:Finite-tensor-savings-and-exact-Fourier-circuits-September-25-2026,
      author = {{OpenAI}},
      title = {{Finite tensor savings and exact Fourier circuits}},
      howpublished = {OpenAI Math Release preprint
                      \href{https://github.com/openai/math/blob/main/preprints/Finite-tensor-savings-and-exact-Fourier-circuits-September-25-2026/main.pdf}{OAI:Finite-tensor-savings-and-exact-Fourier-circuits-September-25-2026}},
      year = {2026}
    }

    @misc{OAI:An-explicit-power-saving-for-the-exact-discrete-Fourier-transform-September-25-2026,
      author = {{OpenAI}},
      title = {{An explicit power saving for the exact discrete Fourier transform}},
      howpublished = {OpenAI Math Release preprint
                      \href{https://github.com/openai/math/blob/main/preprints/An-explicit-power-saving-for-the-exact-discrete-Fourier-transform-September-25-2026/main.pdf}{OAI:An-explicit-power-saving-for-the-exact-discrete-Fourier-transform-September-25-2026}},
      year = {2026}
    }
