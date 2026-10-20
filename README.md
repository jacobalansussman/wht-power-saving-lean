# Discrete Fourier transform of every length in O(n (log n)^z), z = 1 - 1.7887045e-3, checked in Lean

Author: Jacob Sussman. Repository: https://github.com/jacobalansussman/wht-power-saving-lean. Licence: Apache-2.0 ([LICENSE](LICENSE), [NOTICE](NOTICE)).

**What this is.** A proof, checked by the proof assistant Lean 4, that one fixed program computes the discrete
Fourier transform of every length n in O(n (log n)^z) operations with z < 1. This is OpenAI's own headline
statement, from their "Exact Fourier transforms below n log n" ([openai/math](https://github.com/openai/math),
family 130). The cost model is theirs too: a random-access machine (RAM) with exact arithmetic. I call 1 - z
*the saving*. OpenAI's statement has a saving of 10^-13.

**The new result.** On 20 October 2026 I proved the saving 1.7887e-3, exactly
z = 1 - 17887045/10^10, for the Fourier transform of every length. OpenAI's statement file also states
convolution, and that is proved at the same exponent. It is built on the community's circuit
[#2979](https://github.com/CrocSwap/integer-mult-bounds/pull/2979), Marius-Paul Dumitrean's, with two steps of
mine: a new table of where helper storage goes, which I solved for this purpose with the allocation model and the
placement program of the community's #3815 (Erwan Secousses / ChaloupeMusic), and the *recursion lever* in its
full form (the *full cut*). I call the table A-fc. The saving is above every complex-side construction claimed by
others on the community's repository, proved or not (the complex side is the part of their figure that measures
the same thing as my saving; some research notes there print higher ceilings or targets, which are not
constructions), and 0.45 % above my own previous Lean-checked saving, 1.7806935e-3, the full cut on the
community's #3945, which was the headline of the fifth revision. This is the eighth result of this repository; the
table "All the results" below lists all of them.

| | saving | how it is checked |
|---|---|---|
| **this result** (the full cut on #2979's words with the table A-fc) | **1.7887045e-3** | Lean proof of OpenAI's statement |
| the same, its price rounded to the nearest unit of 10^-10 | 1.7887046e-3 | computed; not proved in Lean, and Lean proves the rate fact behind it false there |
| the full cut on #3945 with every head kept, also proved in Lean here (the seventh result) | 1.7806935e-3 | Lean proof |
| the half cut on #3945 with every head kept, also proved in Lean here | 1.7581520e-3 | Lean proof |
| the lever on #2979 with every head kept, also proved in Lean here | 1.7564707e-3 | Lean proof |
| the lever on #2979 as published, also proved in Lean here | 1.7556473e-3 | Lean proof |
| the lever on #3945 as published, also proved in Lean here | 1.7501109e-3 | Lean proof |
| the circuit it is built on, Marius-Paul Dumitrean's [#2979](https://github.com/CrocSwap/integer-mult-bounds/pull/2979), as the community prices it | 1.7290340e-3 | the community's own programs; merged there as research |
| the same words with the table A-fc and no lever, as the community prices them | 1.6897030e-3 | the community's own programs |
| Dugongue's [#3945](https://github.com/CrocSwap/integer-mult-bounds/pull/3945), as the community prices it | 1.7306632e-3 | the community's own programs |
| the highest complex-side construction claimed on the community's repository on 19 October ([#4088](https://github.com/CrocSwap/integer-mult-bounds/pull/4088), Cinna Davis; the same figure in [#4100](https://github.com/CrocSwap/integer-mult-bounds/pull/4100)) | 1.7309791e-3 | the community's own programs; open then, archived on 20 October; not selected |
| the highest complex-side construction ever claimed there (the chain of pull requests #3934 to #3959) | 1.7411726e-3 | withdrawn: it carried moves that their author had retracted, and another contributor's audit (#4046) disputed it |
| the community's selected record ([#1094](https://github.com/CrocSwap/integer-mult-bounds/pull/1094)) | kappa 1.6726602e-3 (complex side 1.6754627e-3) | reviewed and selected by the community |
| the largest saving on the public Lean leaderboard, [bigospeedrun.com/exact-dft](https://bigospeedrun.com/exact-dft) | 1.2654e-3 | Lean proof (Bojnal's, on canghaimeng's C588AT, which runs on my engine) |
| my previous Lean-checked saving before the lever, the sixth result (first published in the fifth revision) | 1.6256629e-3 | Lean proof |

The result is 3.45 % above the figure of the circuit it is built on and 3.35 % above #3945's. Its price is
1.78870456e-3: rounded to the nearest unit of 10^-10 that is 1.7887046e-3, and rounded down it is 1.7887045e-3,
the figure proved in Lean; Lean also proves that the rate fact behind it is false one unit higher. The community's
figures in the table are those of 19 and 20 October 2026; the community moves fast, and section 9 says where each
one was read.
Ross Richardson's [#4091](https://github.com/CrocSwap/integer-mult-bounds/pull/4091) certifies 1.7309803e-3 as the
optimum of a pooled (relaxed) model of #4088's words: a ceiling, not a construction. On 20 October
[#4181](https://github.com/CrocSwap/integer-mult-bounds/pull/4181) (Rohan Arun, open) claimed a fill that attains
that ceiling; my figures are above both.
The community's figures compare with mine only on what they call the *complex side*: their headline figure for
integer multiplication also depends on a second ingredient, the *bit side*, which this work does not touch.

**The idea, in one paragraph.** The program is recursive. It cuts a long transform into a fixed pattern of
smaller transforms, the community's circuit, and does each of those by calling itself. How fast such a program
is depends on one balance between the sizes of the smaller transforms and the size of the job. In the
community's recursion every array of the circuit always takes a full-size job. But in each round of the circuit
used here, 16,576 of the smaller transforms are *heads*: each one is the very first thing done to its
array, and it touches no other array. On those arrays the circuit can take a smaller job instead, and the head
itself leaves the circuit: its work is done elsewhere in the recursion. The recursion then has to handle jobs of
two sizes at once, in a fixed ratio, and the proof shows that it still runs at the speed the balance predicts,
provided every piece of paid work is smaller than every job. With the whole head removed, the smaller job has 25
axes, and in the table used here every other block has rank 23 or less, so nothing else has to change. This step
is what I call *the lever*, and removing the whole head is its *full cut*. The table matters: the community's
tables for #2979 leave thousands of other idle blocks of rank 48, which the full cut has to split into two blocks
of rank 24 each, and the table A-fc was solved so that none is left. Without the lever it is worse than #2979's
own table; with the full cut it is better. "What is new in the sixth revision", below, gives the numbers, and
"What is new in the fifth revision" explains the lever itself.

**What is new in this revision.** This is the sixth revision of the repository. It adds one result, the one above
(`package/leverAFc`), and keeps everything of the fifth revision, whose headline was the full cut on #3945,
1.7806935e-3. The research finding for the table A-fc changes from "checked by computer" to "proved in Lean" (at
1.7887045e-3; the rounded price 1.7887046e-3 is not proved in Lean). The regeneration scripts of three packages
run again from the published files, and their instructions now say that they need Python 3.11 or earlier ("What
is new in the sixth revision"). The fifth revision came after the fourth, published on 10 October 2026 with a
saving of 8.7625e-4, and added:

- **The lever, proved in Lean** (19 October 2026): a recursion with jobs of two sizes (Theorem 1G and its Lean
  proof inside OpenAI's machine model), the *head lemma* (removing half of each head, or all of it, leaves a
  circuit that does the smaller job on the head's array and the full job everywhere else), the split of a block
  of rank 48 into two blocks of rank 24 (an identity, proved in Lean with no condition on the block), a Lean check that the heads are heads
  (a *heads-first certificate*), and the community's circuits #3945 and #2979 ported to Lean, with a device the Lean route
  had never handled: *borrowed* helpers whose coordinates are changed for free by a fixed linear map (the
  community calls them *warps*). Proved in Lean on four pairs of circuit and table of hosts, five results in all: #3945 with every head kept, with the full cut (1.7806935e-3, the fifth revision's headline) and with the half cut (1.7581520e-3), and, with the half cut, #2979 with every head kept (1.7564707e-3), #2979 as published (1.7556473e-3) and #3945 as published (1.7501109e-3).
- **Every head kept.** The community's programs that decide where helper storage goes put each empty block of
  rank 48 either where it is a head or where it never is, at the same price, as their solvers happened to. I re-ran the community's own
  allocation programs with one more condition, that every block that can be a head stays empty, and the same
  circuits then have more heads: 80,640 instead of 57,312 on #3945, and 16,576 instead of 16,086 on #2979. Only
  the table of where helper storage goes changes; the words, the order of the runs and the price of the circuit
  stay the community's, and the community's own checks pass on the new tables.
- **Four results proved in Lean on 13 and 14 October and not published until now**, each a rung of a ladder on
  the community's circuits: 1.3552606e-3 (mutable sources on Jeremy Kothe's circuit), 1.5699305e-3 (with same-run
  hosting), 1.6197705e-3 (the fifth result: exit frames, mutable sources and hosting on the community's U8
  circuit) and 1.6256629e-3 (the sixth result: the same with no source copies left and a tighter bank table),
  with the smaller rungs of the same days. They stayed unpublished because they did not beat the community's own
  figures. The lever is built on top of them.
- **The research of 10 to 20 October**, in [research/](research/README.md): twelve sections, 585 findings, each with its status (24 proved in Lean, 211 checked by computer, 216 argued by hand and not checked, 92 refuted, 42 open),
  and 122 open questions, each with a first experiment. It holds the moonshot waves of 14 and 15 October,
  which looked for a small unit with a large saving at blocks 2 to 5 in every game I could define and found
  none; the lower bounds found along the way (three proved in Lean); why the record circuits save; the
  routes to proving that a unit exists without finding it; the links to other fields; the speed waves on real
  processors and FPGAs; and the record track that led to the lever. Failures are told as plainly as successes.

**Whose work it rests on.** Most of what this result stands on is other people's, and the credit for it is
theirs. [NOTICE](NOTICE), section 10, lists every contribution with its pull request; in short:

- **The circuits** of the lever are the community's. #3945 is Dugongue's (made with OpenAI Codex), the
  newest in a line of circuits by many contributors: dataisfire made the staged circuit it starts from, zxv the
  staged pair, Ross Richardson the staged entrance, Erwan Secousses (ChaloupeMusic) the clock and the table of
  hosts (the programs that made the every-head table on #3945 are his), hexapode the endpoint combination,
  roger-sato the retimed releases, Jeremy Diviney and dcharb78 the coupled carriers, and Marius-Paul Dumitrean
  the joint end frames. #2979 is Marius-Paul Dumitrean's (made with Anthropic's Claude), on Jeremy Kothe's gate
  rewrites (#2957) and Mason Hale's words, chronology and allocation checker (#2917), with Jeremy Diviney's
  paired extraction, Yiyang Lu's shared kernels and allocation model, Dugongue's owner engine, Ross Richardson,
  Tekalign, icekylinx, themilkman216 and Gareth Davidson; the programs that made the every-head table on #2979
  are #2979's own. The table A-fc of the result above was made with the allocation model and the placement program
  of #3815, Erwan Secousses's (ChaloupeMusic). Before them came Jeremy Kothe, Lyy-iiis, philippeVV, Mason Hale, Sagar Thakrar (with Arham
  Tahir's audit), Tekalign, themilkman216, Vlado Vojdanovski and the others NOTICE names.
- **The borrowed helpers and their warps**, with the exact Clifford machinery that checks them, are Cinna Davis's
  (1337isnot1337, with OpenAI Codex), on icekylinx's native line and the wrappers of DreamingOfClouds and
  chafreaky. Tekalign added the lifted targets, zxv the source-partner aliases, and Jeremy Kothe and philippeVV
  the hosting of a whole helper inside another one before it starts.
- **The layout of copies underneath** grows from the community's circuit T74: Rohan Arun's, with Marius-Paul
  Dumitrean, Jeremy Kothe and LeverRank, on Jeremy Kothe's earlier words with Thomas DiFiore, Jeremy Diviney and
  Chafik Boukhalfa (chafreaky); Vlado Vojdanovski, Chafik Boukhalfa and Jeremy Kothe made its assembly.
- **The price of a circuit** is the community's moment formula. The agents re-implemented the community's exact
  moment engine (the package "native complex moment" by Rohan Arun with OpenAI Codex, on work of Dugongue,
  eumemic and icekylinx) and reproduced the community's figures to 24 digits before using it.
- **Checks of the circuits by others:** Chris R (Clear1nsight) audited #3945, and ChaloupeMusic replayed its exact
  operator, on the operator model of Nick Petalas and Matarener.
- **Ryan Shea** turned Jeremy Kothe's earlier circuit into a Lean proof of OpenAI's statement
  ([shea256/fourier-transform-below-nlogn](https://github.com/shea256/fourier-transform-below-nlogn), folder
  `lean/complex588`). His tree is the Lean project every result since the fifth is laid over, and his conversion
  pipeline and generators are the base of mine.
- **Nil Mamano** wrote the Lean files of the nine-stage network and the templates of the Fourier chain that
  Shea's tree takes over; **Dugongue** made its entrance gauges and banks, **Zhihao Chen** (jacklightChen) its
  compressed layout, **lydakis** the bank patterns and **DaysSky** the layout it starts from. DaysSky proposed
  hosting helpers in idle arrays and Zhihao Chen wrote the argument for mutable sources, both in the community's
  terms; Jeremy Kothe and DreamingOfClouds wrote the community's lemmas on where a helper may end.
- **The ladder of 13 and 14 October** is built on the community's U8 circuit (Jeremy Kothe's, with Chafik
  Boukhalfa's frame changes, on icekylinx's construction "gauge44") and on Jeremy Kothe's circuit of Shea's proof;
  "The results of 13 October 2026", below, gives its credits in full.
- **Schönhage's τ-theorem** for matrix multiplication is the idea that Theorem 1G extends.
- **OpenAI** made the cost model, the statement, and the step from the lengths 2^k to every length.

What is mine:

- the mixed-target recursion: Theorem 1G and its Lean proof in OpenAI's model;
- the lever: finding the heads, the half cut and the full cut with its split, which keep the recursion's
  condition, the head lemma, and the tables with every head kept (made with the community's own allocation programs and one added condition);
- the table A-fc on #2979's words (its class solution, solved for the full cut, and its tolerant placement);
- the Lean ports of the community's newest circuits: the conversion of their borrowed-helper words into my
  certificate format, the Lean checker for warps with its soundness proof, the networks of many copies with
  hosting across copies, and the generators of the Lean data;
- the Lean engine underneath, and its extensions for exit frames, mutable sources, centre reading and hosting
  (the fifth and sixth results);
- the research in [research/](research/README.md).

The agents (Claude) wrote the proofs and programs under my direction ("Who made it, and how", below).

**How it was checked, and the limits.** Every check below passed for the result above, on Linux machines of
Google's Colab service (8 processors, about 50 GB of memory). The package's `README.txt` and `result.json` list
each check with its log:

- Lean's kernel accepted the 280 new modules, compiled from source over Shea's tree with `package/u8x0hq`, on the
  machine where they were built and again, from the package's files alone, on a second machine. The two theorems
  depend only on the three standard axioms (`propext`, `Classical.choice`, `Quot.sound`), and the source files
  contain no `sorry` other than the two lines that OpenAI's statement file has too.
- `tools/Compare.lean` printed `RESULT: PASS`: the theorems proved are those of OpenAI's statement file,
  `UniformFourier.lean`, with only the exponent and seven names changed. The same comparison against a statement
  one unit of 10^-10 stronger fails, as it must.
- The official comparator, the program that OpenAI's repository uses, accepted the proof with both of its kernels
  (nanoda and Lean's) and rejected the statement one unit stronger.
- The agents re-created the leaderboard's verifier offline, on the proof packaged as a submission, and its phases
  passed: the Lean check, which compiled the packaged proof from source in 63 minutes, the audit of the axioms, and
  the nanoda kernel; the same run one unit stronger was rejected. For the leaderboard's limit of 500 files the
  packaged proof merges small modules into larger ones (449 files); the unmerged form (577 modules) passed the same
  run. A separate check compared the merged form, constant by constant, with the proof built from the package's
  files: the 21,665 constants that are not generated by Lean itself are the same in both, and all 9,601 theorems
  have the same statement and the same proof.
- The rate fact behind the figure is proved in Lean at 1.7887045e-3, and Lean also proves that it fails at every
  saving one unit higher. Lean's kernel counts the heads in the table of hosts (16,576) and refuses one more, and
  Lean proves that the full cut's profile is wrong with one head fewer, and wrong for the table with every head kept
  of the fifth revision, whose heads removed still leave 11,614 blocks of rank 48.
- The one fact of the full cut that depends on the circuit's words, that the 48 directions of each head are the 48
  lines that the cut removes, is proved in Lean on #2979's words. In the fifth revision it was proved only on
  #3945's words, and on #2979's it was argued by hand.
- The agents ran the community's own checks of #2979's line on the new table, unmodified, on a fresh clone of the
  community's repository: #2979's owner engine with its corruption controls, its verify program, the complete
  replay's fresh native runs and the other steps of its admission pass. The table is rejected only by the steps
  that ask for a higher figure without the lever (the strict-gain steps and a near-optimality bound), as it must
  be: without the lever it is worth 1.6897030e-3, against #2979's 1.7290340e-3.
- Three more agents, each with the sole job of finding a fault, checked the result again on other machines: one
  fetched OpenAI's statement afresh, derived the statement files from it once more, rebuilt the proof from the
  package's files, ran the official comparator again, and walked the constants the proof uses, to see that it
  reaches the new table's data and the cut fact on #2979's words and no data of any other table; one rebuilt every
  new module from the package alone on a fresh copy of the base and ran `tools/Compare.lean` and the official
  comparator again; one decoded the Lean data files with programs of its own, rebuilt the table from them byte for
  byte, re-derived the heads, the profile, the separation and the rate fact (true at the figure, false one unit
  higher), checked the hypotheses of Theorem 1G, and re-read the community's figures. None found a fault in the
  proof. They found one fault in the package: the program that regenerates its Lean data gives different bytes for
  36 data files on Python 3.12 or later, because Python 3.12 changed how it adds floating-point numbers; the files
  Lean checked are not affected, and the package's instructions now say to use Python 3.11 or earlier, with a
  workaround. Known cosmetic issues, listed with the correct text in `package/leverAFc/KNOWN-ISSUES.txt` and left as
  compiled, because the files are checked byte for byte: comments that name the earlier builds a file was made
  from; the comment above the Fourier theorem calls it a sibling of a Walsh-Hadamard theorem
  `wht_main_block_L2979Fcx` that does not exist (it should say only that the theorem follows from the package's
  kernel program); and the header of `Goal.lean` cites OpenAI's lines 75-78, 80 and 82 where the definitions are at
  75-77, 79 and 81.
- This repository's own scripts (`scripts/`), run as VERIFY.md section 12 says from a fresh clone of this revision
  on a fresh copy of Shea's tree, laid `package/leverAFc` over the tree, built it and passed every check of
  `scripts/build-and-check.sh` (`ALL CHECKS PASSED`;
  [scripts/logs-2026-10-19-fresh-clone/](scripts/logs-2026-10-19-fresh-clone/README.md)). That script runs five of
  the package's twelve Lean controls; the other seven need files of other packages or side files, and VERIFY.md
  section 12.3 says how to run them.

The limits are properties of how it was built and checked:

- Like every figure in this race, the result is asymptotic. The constant hidden in the O is astronomically
  large, so it says nothing about running times at sizes anyone can run.
- The proof lives in Shea's Lean tree, laid over it as two packages (`package/u8x0hq`, then
  `package/leverAFc`), not in this repository's own Lean project.
- Every Lean run was on Colab machines. The comparator's sandbox needs a feature of the Linux kernel (Landlock)
  that those machines lack, so it ran inside a stand-in; the re-creation of the leaderboard's verifier used the
  agents' own scripts in the same stand-in, and the verifier's own controller never ran.
- The community's figures that the result is compared with are their claims, checked by their programs. My
  result does not rest on those checks or on the community's premises: the circuit enters the Lean proof only as
  data, which Lean checks.
- Theorem 1G was first argued by hand in a simple model of jobs; what is proved in Lean is its form inside
  OpenAI's model for rules with two job sizes that are separated, as used here. Rules that are not separated, where some paid block is larger than the
  smaller job, are not covered by any theorem ([notes/open-directions-4.md](notes/open-directions-4.md)).
- No person outside this project has reviewed the new proofs, and I would welcome a line-by-line review.

**All the results.**

| Result | What is computed | Saving | Where the proof is | Date |
|---|---|---|---|---|
| eighth | Fourier transform of every length, and convolution: the community's #2979 words with the table A-fc, solved for the full cut, and the recursion lever (16,576 heads, full cut, nothing to split) | 1.7887045e-3 | `package/leverAFc`, on Shea's tree | 20 October 2026 |
| seventh | Fourier transform of every length, and convolution: the community's #3945 with every head kept, and the recursion lever (80,640 heads, full cut, the other 57,312 blocks of rank 48 split) | 1.7806935e-3 | `package/leverF`, on Shea's tree | 19 October 2026 |
| with the seventh | Fourier transform of every length, and convolution: the community's #3945 with every head kept, and the recursion lever (80,640 heads, half cut) | 1.7581520e-3 | `package/leverB`, on Shea's tree | 19 October 2026 |
| with the seventh | Fourier transform of every length, and convolution: the community's #2979 with every head kept, and the recursion lever (16,576 heads, half cut) | 1.7564707e-3 | `package/leverA`, on Shea's tree | 19 October 2026 |
| with the seventh | Fourier transform of every length, and convolution: the community's #2979 as published, with the recursion lever (16,086 heads, half cut) | 1.7556473e-3 | `package/lever2979`, on Shea's tree | 19 October 2026 |
| with the seventh | Fourier transform of every length, and convolution: the community's #3945 as published, with the recursion lever (57,312 heads, half cut) | 1.7501109e-3 | `package/lever`, on Shea's tree | 19 October 2026 |
| sixth | Fourier transform of every length, and convolution: the community's U8 circuit with exit frames, mutable sources without centre copies and same-run hosting (one-pad bank table) | 1.6256629e-3 | `package/u8x0hq`, on Shea's tree | 14 October 2026 |
| fifth | Fourier transform of every length, and convolution: the community's U8 circuit with exit frames, mutable sources and same-run hosting | 1.6197705e-3 | `community/u8xmhh`, on Shea's tree | 13 October 2026 |
| with the fifth | Fourier transform of every length: Kothe's circuit with mutable sources and same-run hosting | 1.5699305e-3 | `community/c588mhh`, on Shea's tree | 13 October 2026 |
| with the fifth | Walsh-Hadamard transform, lengths 2^k, on that kernel program | 1.5699306e-3 | `community/wht-c588mhh`, on Shea's tree | 13 October 2026 |
| with the fifth | Fourier transform of every length: the community's U8 circuit with exit frames and mutable sources | 1.3868856e-3 | `community/u8xm`, on Shea's tree | 13 October 2026 |
| with the fifth | Fourier transform of every length: Kothe's circuit with mutable sources only | 1.3552606e-3 | `community/c588m`, on Shea's tree | 13 October 2026 |
| with the fifth | Walsh-Hadamard transform, on that kernel program | 1.3552607e-3 | `community/wht-c588`, on Shea's tree | 13 October 2026 |
| with the fifth | Fourier transform of every length: the U8 circuit with exit frames | 1.2704023e-3 | `community/u8x`, on Shea's tree | 13 October 2026 |
| with the fifth | Fourier transform of every length: Shea's certificate with a repacked bank table | 1.2383291e-3 | `community/c588-repack`, on Shea's tree | 13 October 2026 |
| with the fifth | Walsh-Hadamard transform, on Shea's kernel program and on the repacked one | 1.2331686e-3 and 1.2383292e-3 | `community/wht-c588`, on Shea's tree | 13 October 2026 |
| fourth | Fourier transform of every length, and convolution | 8.762478e-4 | this repository's Lean project | 10 October 2026 |
| fourth | Walsh-Hadamard transform, lengths 2^k | 8.762479e-4 | this repository's Lean project | 10 October 2026 |
| third | Fourier transform of every length, and convolution | 7.474546e-4 | this repository's Lean project | 9 October 2026 |
| second | Walsh-Hadamard transform | 7.474547e-4 | this repository's Lean project | 9 October 2026 |
| first | Walsh-Hadamard transform | 5.399225e-4 | this repository's Lean project | 9 October 2026 |

Every line is checked by Lean's kernel, by `tools/Compare.lean` and by the official comparator. The lines on
Shea's tree differ in what else was done; "What is new in the fifth revision" and "The results of 13 October in
detail", items 1 and 4, give each one's checks, logs and gaps. The sixth result's package also holds the Fourier
chains of four instances built on the way to it (U8XM, U8XMHH, U8X0 and U8X0HH; its `README.txt`). There is no Walsh-Hadamard line for
the kernel programs of the lever, of the sixth result or of the fifth: that corollary is one term, as for the
others, but it has not been built. The Walsh-Hadamard lines on Shea's tree are compared with this repository's own
published Walsh-Hadamard statement, with only the exponent and the theorem name changed. Each Walsh-Hadamard
figure sits one unit of 10^-10 above the Fourier figure on the same kernel program: the step to every length costs
a factor (log log n)^2, so the Fourier statement is taken one unit lower (section 2, point 4).

**Who made it, and how.** I built this repository between 8 and 20 October 2026 with a team of AI agents
(Claude) that I directed. I set the goals, chose which ideas to pursue, and decided at each stage what to do
next. The agents wrote the Lean proofs, designed the networks, ran the searches and the community's checks, and
tried to refute one another's work. For the results since the fifth I chose to build on the best work the
community had published rather than on my own circuit, and to prove the steps its figures assume. The lever's
idea came out of a round on abstract routes on 14 October 2026; the agents counted the heads, proved Theorem 1G on
paper and planned the Lean route that day, and then I paused it for other directions. On 19 October the agents
priced the lever on the community's newest circuits, and I had the proofs built on five machines at once. On
20 October the agents solved the table A-fc for the full cut, and I had it built in Lean the same day, on the same
machines, with three more agents trying to refute it.
[notes/how-the-lever-was-found.md](notes/how-the-lever-was-found.md) tells that story and
[notes/fifth-result-story.md](notes/fifth-result-story.md) the story of 13 October, each with what failed and the
claims that did not stand.

**This is not an OpenAI project, and not the community's.** OpenAI did not write, review or endorse it. Their
material is here in these forms:

- The 89 files under `OAI/` are OpenAI's, unmodified.
- Under `Work/Fourier/` and `Work/FourierE8/` there are copies of OpenAI's files with changes. Each says so in
  its first lines, and [NOTICE](NOTICE), section 7, lists them.
- The packages under `community/` and `package/` hold more such copies. Each copy of OpenAI's statement file among
  the files of `community/` that are laid over Shea's tree (the folders `overlay/` and `variants/`) says so in a
  notice on its lines 3 to 12 (3 to 14 in `community/c588-repack`); the notice was added after the build and the
  checks, and it gives the checksum of the file as it was compiled, which removing the notice gives back. The
  statement copies in `package/` are kept byte for byte as they were compiled, because the scripts derive them
  again from OpenAI's file and compare the two; [NOTICE](NOTICE), section 10, lists them and what was changed. A
  few more changed copies are kept among the logs of `community/` exactly as they were compiled or fetched, so
  they carry no such notice; a file `NOTICE.txt` in each of their folders says what was changed. The copies of
  OpenAI's proof files say so in their first lines. [NOTICE](NOTICE), sections 9 and 10, lists them all.
- The name of the Lean (Lake) package, `OAI`, and the namespace `OAI.PowerSaving` are theirs. They are kept
  only because the new modules extend their development.

Ryan Shea, Nil Mamano, DaysSky and the community's authors have not reviewed it either.

**Where to look next.**

- "What is new in the sixth revision", just below: the table A-fc with the numbers, the theorem and the other
  changes of this revision.
- "What is new in the fifth revision", after it: the lever with the numbers, the theorem, what the proof is made
  of and what is not claimed; the sixth result; the research.
- [research/README.md](research/README.md): the research of 10 to 20 October, section by section, and
  [research/open-problems.md](research/open-problems.md), every open question with a first experiment.
- [notes/how-the-lever-was-found.md](notes/how-the-lever-was-found.md) and
  [notes/open-directions-4.md](notes/open-directions-4.md): the lever's story, and what was left before the full cut
  (written for #3945 as published), each lead with a first experiment.
- [VERIFY.md](VERIFY.md), section 12: how to rebuild and check the lever and the sixth result yourself, and how to
  re-run the searches and prices on a laptop.
- "The results of 13 October 2026", further down: the fifth result and its rungs, with every check and its log.
- The fourth result and the earlier ones follow, from "What was new in the fourth revision" on. Section 9 has
  the dates, commits and outside sources behind every statement, kept out of the prose.

## What is new in the sixth revision (2026-10-20)

### The full cut on #2979's words with a table solved for it (A-fc)

The words used here (array, block, rank, job, unit, the balance, heads, the half cut and the full cut) are
explained in "What is new in the fifth revision", below.

**The circuit and the table.** #2979, Marius-Paul Dumitrean's circuit, runs 74 copies of a pattern of 9 calls in
each unit: 666 calls of two provider words, over 66,304 arrays per copy of its gauge, with 271,728 helper roles that
a table of hosts assigns to calls. The words, the clock, the order of the calls, the dispatch and the providers stay
#2979's here; only the table is new. The fifth revision's table with every head kept on #2979 (`package/leverA`)
has 16,576 heads, but it also leaves 11,614 other idle blocks of rank 48, which the full cut has to split into two
blocks of rank 24; its full cut prices at 1.7797095e-3, below the full cut on #3945. The table A-fc was solved for
the price after the full cut instead. With the allocation model of the community's #3815 (Erwan Secousses /
ChaloupeMusic), a linear program for each tick with the full cut's objective, and then an integer program on the
columns it generated, gave a class solution; #3815's own placement program, which works in reverse time, then
placed the roles. I made that program tolerant: a role with no planned slot left at its call goes to a free place
that obeys the same rules of windows, calls and capacity. Nine roles needed that (eight of width 6 at call 553 and
one of width 4 at call 0).

**The numbers.** In A-fc every bin of the two head positions of each copy is empty, so its 16,576 blocks of rank 48
(8,288 in each position) are all heads, and no other idle block has rank 48: the other idle blocks have ranks 1, 6,
7, 17, 19 and 23. One copy of the gauge does 1,630,639 blocks of rank 1, 470,196 of rank 2, 154,512 of rank 3,
71,928 of rank 4, 11,988 of rank 5, 5 of rank 6, 3 of rank 7, 5,328 of rank 8, 8 of rank 17, 1 of rank 19,
25,905 of rank 23 and the 16,576 heads of rank 48, on W = 66,304 arrays. With the full cut the heads leave, their
arrays take jobs of 25 axes, and the balance becomes

    (the blocks of ranks 1 to 23 above, each counted as r^τ)  =  49,728 · 73^τ + 16,576 · 25^τ .

Every paid block has rank at most 23 and the smaller job has 25 axes, so the rule is separated, Theorem 1G
applies, and nothing has to be split. The root gives a saving of 1.78870456430e-3. Rounded down to the grid of
10^-10 that is 17,887,045/10^10, the figure proved in Lean; rounded to the nearest it is 1.7887046e-3, which is not
proved in Lean, and Lean proves that the rate fact fails there. Without the lever the table is worse than #2979's
own: the community's formula gives it 1.6897030e-3, against 1.7290340e-3. An exact bound from a linear program
(checked by computer; [research/record-track/](research/record-track/README.md), section 4.8) puts the full cut on
#2979's words, in the allocation model used there, at no more than 1.7893658e-3, so A-fc is within 6.7e-7 of that
bound.

| table on #2979's words | heads | other idle blocks of rank 48 | without the lever | price with the full cut | in Lean |
|---|---|---|---|---|---|
| A-fc (this revision) | 16,576 | 0 | 1.6897030e-3 | 1.78870456e-3 | 1.7887045e-3, proved in Lean (`package/leverAFc`) |
| every head kept (fifth revision) | 16,576 | 11,614 | 1.7290340e-3 | 1.7797095e-3 | the half cut, 1.7564707e-3, proved in Lean (`package/leverA`); the full cut not built |
| as published | 16,086 | 12,104 | 1.7290340e-3 | 1.7744822e-3 | the half cut, 1.7556473e-3, proved in Lean (`package/lever2979`); the full cut not built |

**What is new in the Lean proof.** Almost nothing. The proof reuses the fifth revision's parts: Theorem 1G in
OpenAI's model, the head lemma with its full cut, the heads-first certificate, the port of #2979's words with their
warps and the network of copies. New are the table's data (the hosts of all 271,728 roles, checked by Lean's
kernel), its instance of the full cut, its rate fact and the Fourier chain at its exponent, and one fact: the full
cut needs the 48 directions of each head to be the 48 lines that the cut removes. On #3945's words that followed
from a fact proved there; on #2979's words it is now proved in Lean as well (the gate of the build,
`package/leverAFc/logs/afc-gate/`). The full cut's bridge always splits the blocks of rank 48 that are not heads;
on A-fc there are none, so the split changes nothing (`KNOWN-ISSUES.txt`, item 6).

**The theorem.**

| | |
|---|---|
| statement | OpenAI's `UniformFourier.lean` (family 130), with the exponent `1 - 17887045/10^10` and seven names changed |
| theorems | `OAI.PowerSaving.transform_mainL2979Fc` (the transform of every length) and `OAI.PowerSaving.convolution_mainL2979Fc` (convolution) |
| file | `Work/FourierL2979Fc/Main.lean` in `package/leverAFc/overlay/` |
| statement files compared | `Work/FourierL2979Fc/UniformFourierChallenge.lean` (this exponent) and `Work/FourierL2979FcNeg/UniformFourierChallenge.lean` (one unit stronger) |
| axioms | `propext`, `Classical.choice`, `Quot.sound` |
| base | Ryan Shea's `lean/complex588` at commit `5239484`, with `package/u8x0hq` and then `package/leverAFc` laid over it |

Its checks are listed at the top ("How it was checked, and the limits"), and each with its log in the package's
`README.txt` and `result.json`.

**Other changes in this revision.**

- The regeneration scripts (`generators/REGENERATE.sh`) of `package/lever2979`, `package/leverA` and
  `package/leverAFc` run again from the published files. Two reports in their generators name the machine paths of
  the run that made them; those paths were shortened for publication, which changed the two files' checksums, so
  in the fifth revision the scripts of `package/lever2979` and `package/leverA` stopped at their first check. The
  checksum lists now give the published bytes, and the scripts put the reports' original checksums back where the
  data records them. On a Mac with Python 3.9, all three end with `REGENERATE: ALL AS EXPECTED`, every generated file
  byte for byte as compiled. They need Python 3.11 or earlier: with Python 3.12 or later, 36 data files of the
  provider words come out different (the packages' `generators/README.txt` gives a workaround).
- `package/leverAFc/result.json` lists under `lean_controls` only the five Lean controls that this repository's
  `scripts/build-and-check.sh` can run in a tree with this package alone; the other seven need files of
  `package/leverA` or side files of the package, and are listed apart with what each needs. All twelve were
  compiled on the build machine.
- In [research/](research/README.md), the record track's finding for A-fc is now "proved in Lean at 1.7887045e-3";
  its price, 1.7887046e-3 rounded, stays a computed figure.

**What is not claimed.** The list "What is not claimed for the lever" in the fifth revision, below, holds for this
result too: nothing about integer multiplication (the community's kappa is capped by a bit side that this work does
not touch), no speed-up in practice, nothing about the community's premises, and no claim to be the largest saving
anyone will have. In addition, the figure 1.7887046e-3 is a price, not a Lean result.

## What is new in the fifth revision (2026-10-20)

### The lever, with the numbers

**Words used here.**

- **Saving.** The program needs O(n (log n)^z) operations; the saving is 1 - z. Larger is better.
- **Array, axis.** The program keeps its data in many *arrays*. In the circuits used here every array carries 73
  coordinates, its *axes*, and the work on an array is a transform along some of its axes.
- **Block, rank.** A costly step of the circuit is a *block*: a transform along some number r of one array's
  axes, which the program does by calling itself on a smaller problem. r is the block's *rank*. Additions between
  arrays are free.
- **Unit.** One full pass of the community's circuit: on #3945, 360 copies of a smaller circuit, 3,240 runs in all,
  over W = 322,560 arrays; on #2979, 74 copies and 666 runs.
- **Job.** What the recursion asks of one array in one unit: a transform along k of its axes. In the community's
  recursion k is always 73, the *full job*.

**The balance.** For a recursive program of this kind the saving is fixed by one equation. If one unit does
blocks of ranks r_1, r_2, ... on W arrays that each take a job of 73 axes, the operation count grows like
n (log n)^τ, where τ is the number that balances

    r_1^τ + r_2^τ + ...  =  W · 73^τ ,

and the saving is 1 - τ. This is the community's moment formula. For #3945 the unit's blocks are 8,041,680 of
rank 1, 2,306,880 of rank 2, 729,648 of rank 3, 353,160 of rank 4, 69,120 of rank 5, 25,920 of rank 8 and
137,952 of rank 48, and the balance gives the community's figure, 1.7306632e-3.

**Heads, and the half cut.** The blocks of rank 48 finish off parts of arrays in which no helper storage was
placed. The agents followed every array through the unit and found that some of them are *heads*: the first
operation on their array, acting on that array alone, with nothing kept in that array before it. In #3945 as the
community published it, 57,312 of the 137,952 are heads; the other 80,640 are not, because other steps work on
their arrays before them.

A transform along 73 axes is the same as one along 24 axes followed by one along the other 49, because transforms
along different axes commute. So an array whose first step is a head along 48 axes can be given a smaller job:
the 24 axes that the head would start with are taken out of this unit's job on that array and done elsewhere in
the recursion. On that array the unit then has a job of 49 axes, and the head needs only its other 24 axes.

The balance changes on both sides. With h heads among the rank-48 blocks, each head turns a 48 on the left into a
24, and each head array turns a 73 on the right into a 49. On #3945 as published:

    (other blocks) + 80,640 · 48^τ + 57,312 · 24^τ  =  265,248 · 73^τ + 57,312 · 49^τ .

Both sides lose about the same amount, but not exactly: x^τ grows a little more slowly for large x when τ is
below 1, so the step from 48 to 24 removes slightly more than the step from 73 to 49. The left side falls further,
the balance moves to a smaller τ, and the saving grows from 1.7306632e-3 to 1.7501109e-3. More heads give more:
that is why the tables with every head kept are worth building.

**Why only half of each head.** A recursion with jobs of two sizes is not the community's recursion, so it needs
its own theorem. **Theorem 1G** (the agents' proof, first in a simple model of jobs) says: if every block that is
paid for is smaller than every job (the rule is *separated*), the program runs in O(n (log n)^τ') operations for
every τ' above the balance root τ, with a constant that depends on τ'. The construction is a tree of calls run in
lockstep, so that every call receives exactly the jobs of both sizes it needs, and the few jobs that are left over
are done by fresh calls of the whole program. Here every paid block has rank at most 48 and the smaller job has 49
axes, so the rule is separated. Cutting the whole head instead (a job of 25 axes) would leave the other blocks of
rank 48 larger than the smaller job. They can be split into two blocks of 24: a block along 48 axes is the same as
one along 24 of them followed by one along the other 24. On #3945 as published that split costs back almost all
of the gain: the full cut gives 1.7307680e-3 against the half cut's 1.7501109e-3. With every head kept (below)
there are fewer blocks to split and more heads to remove, and the full cut is the best: 1.7806935e-3.

**Every head kept.** In the community's allocation model an empty block of rank 48 can sit in two places in each
copy at the same price: in one it is a head, in the other it never is. The community's programs split them as
their solvers happened to. Re-running those programs with every head position kept empty changes only the table
of where helper storage goes. On #3945 the heads become 80,640 (the most there can be, 224 per copy) and the other
57,312 are not heads; on #2979 the heads become 16,576 of 28,190 (8,288 in each of the two head positions,
against 16,086 as published). The paid blocks, the community's figure for the circuit and its checks stay the
same.

**The figures.** Every figure here except those of the last column is a price by the balance, computed to at
least ten digits and rounded down to the grid of 10^-10 ([tools/lever/README.md](tools/lever/README.md), item 2).
The last column gives the figure proved in Lean, a point of the same grid at or below the price.

| circuit and table of hosts | heads | the circuit's own figure | with the half cut | with the full cut | in Lean |
|---|---|---|---|---|---|
| #3945 (Dugongue) with every head kept | 80,640 | 1.7306632e-3 | 1.7581520e-3 (+1.59 %) | 1.7806935e-3 (+2.89 %) | both proved in Lean: 1.7806935e-3 (full cut, `package/leverF`) and 1.7581520e-3 (half cut, `package/leverB`) |
| #2979 (Marius-Paul Dumitrean) with every head kept | 16,576 | 1.7290340e-3 | 1.7564708e-3 (+1.59 %) | 1.7797095e-3 | 1.7564707e-3, proved in Lean (`package/leverA`) |
| #2979, as published | 16,086 | 1.7290340e-3 | 1.7556473e-3 (+1.54 %) | 1.7744822e-3 | 1.7556473e-3, proved in Lean (`package/lever2979`) |
| #3945, as published | 57,312 | 1.7306632e-3 | 1.7501109e-3 (+1.12 %) | 1.7307680e-3 | 1.7501109e-3, proved in Lean (`package/lever`) |
| #3815 (ChaloupeMusic), the fallback for #3945 | 57,312 | 1.7305784e-3 | 1.7500242e-3 (+1.12 %) | - | not built |
| #4076 (Rohan Arun, a draft built on #3945's words) | - | 1.7308090e-3 | 1.7504268e-3 (+1.13 %) | 1.7315054e-3 | not built |
| T74 (#1839, Jeremy Kothe), where the lever was first found | 11,010 | 1.6797550e-3 | 1.6968621e-3 (+1.02 %) | 1.6948358e-3 | not built |

The figures of the full cut are prices, not built, except on #3945 with every head kept, where it is proved in
Lean (`package/leverF`). The full cut on #2979 with every head kept (1.7797095e-3) would need the data of that
table and the cut fact proved on #2979's words: the Lean bridge of the full cut (`Work/Mixed/Bridge/CoreF`, listed in
`package/leverF/README.txt`) is general in the table, but in the fifth revision its cut fact was proved in Lean only on #3945's words. On #2979's words a re-solved table, A-fc, prices at 1.7887046e-3 under the full cut (the price 1.7887045643e-3
rounded to the nearest; on the grid of 10^-10, rounded down as above, 17,887,045/10^10), above the Lean result on #3945.
In the fifth revision it was checked by computer only; the sixth revision proves it in Lean at 1.7887045e-3
(`package/leverAFc`, "What is new in the sixth revision", above; its table and the programs that made it are also in
[research/record-track/data/a-fc/](research/record-track/data/a-fc/README.md)).

**The theorem.**

| | |
|---|---|
| statement | OpenAI's `UniformFourier.lean` (family 130), with the exponent `1 - 17806935/10^10` and seven names changed |
| theorems | `OAI.PowerSaving.transform_mainL3945F` (the transform of every length) and `OAI.PowerSaving.convolution_mainL3945F` (convolution) |
| file | `Work/FourierL3945F/Main.lean` in `package/leverF/overlay/` |
| statement files compared | `Work/FourierL3945F/UniformFourierChallenge.lean` (this exponent) and `Work/FourierL3945FNeg/UniformFourierChallenge.lean` (one unit stronger) |
| axioms | `propext`, `Classical.choice`, `Quot.sound` |
| base | Ryan Shea's `lean/complex588` at commit `5239484`, with `package/u8x0hq` and then `package/leverF` laid over it |

The other lever packages have the same form with their own suffix (the table "All the results").

**What the proof is made of.** Each lever package's `README.txt` lists its folders and `MODULES.txt` every module
in build order. In short:

| part | what it proves | where |
|---|---|---|
| the flow and the cost bound | Theorem 1G's bookkeeping: the weights of the tree of classes, exact balance, sources and exits, and the cost bound for every exponent above the rate point | `Work/Mixed/` (Flow, Phi, LemmaA, Bound, TwoTarget) |
| the program | the lockstep tree of calls as a program in OpenAI's machine model, its correctness and its cost | `Work/Mixed/Prog/`, `Work/MixedCost/` |
| the seam and the rate fact | OpenAI's statement accepts an engine whose cost has a constant factor; the balance holds at the figure and fails one unit above | `Work/Mixed/` (SeamConst, Rate*, RateMono) |
| the head lemma | removing heads, or half of each, leaves a circuit that does the smaller job on the head arrays and the full job elsewhere; any block of rank a + b is two blocks | `Work/Mixed/HeadM1` and, in some packages, `HeadComm` |
| the heads-first certificate | the circuit equals the same circuit with its heads moved to the very front, as a first step on arrays that nothing has touched yet; Lean's kernel counts the heads in the table of hosts | `Work/Mixed/Bridge/`, `Work/XHostT/HF*/` |
| warps | the Lean checker for the community's borrowed helpers, and what a warp means as a matrix | `Work/KLowW/`, `Work/Warp/`, `Work/PWGen/` |
| the network of copies | the hosted network with any number of copies, any order of runs, and helper storage hosted in idle arrays of any copy | `Work/XHostD/`, `Work/XHostR/`, `Work/XHostT/` |
| the circuit as data | both words of the circuit, the dispatch of its runs and the table of hosts, each checked by Lean's kernel | `Work/GCert/Data/`, `Work/XHostT/Data/` (generated) |
| the Fourier chain | the step from the kernel program to every length, at this exponent | `Work/Fourier<suffix>/` (generated) |

Every generated file can be made again from the community's files: each package's `generators/REGENERATE.sh`
does it and compares the result with the package's `SHA256SUMS`. Each generator that writes a figure refuses a
figure one unit higher.

**What is not claimed for the lever.**

- Nothing about integer multiplication. The community's headline figure, kappa, is capped by the complex side b (through b/(1+b): the selected record's b of 1.6754627e-3 gives its kappa of 1.6726602e-3) and by a bit side; the lever raises only the complex side, so the community's kappa would rise only as far as its
  bit side.
- No speed-up in practice (the constant is astronomically large).
- Nothing about the community's premises beyond what the Lean proof itself establishes. The proof is of OpenAI's
  statement in OpenAI's model; the community's own figures rest on premises that they state, which this
  repository neither uses nor checks.
- No claim to be the largest saving anyone will have: the community's figures that are not Lean-checked are listed
  in the table at the top, and higher claims may appear at any time.

**How the full cut on #3945 was checked.** Every check below passed for the full cut on #3945 with every head
kept (`package/leverF`), the headline of the fifth revision, on Linux machines of
Google's Colab service (8 processors, about 50 GB of memory). Each lever package's `README.txt` and `result.json`
list its own checks with their logs:

- Lean's kernel accepted the new modules (410 for the full cut), compiled from source over Shea's
  tree with `package/u8x0hq`. The two theorems depend only on the three standard axioms (`propext`,
  `Classical.choice`, `Quot.sound`), and the source files contain no `sorry` other than the two lines that
  OpenAI's statement file has too.
- `tools/Compare.lean` printed `RESULT: PASS`: the theorems proved are those of OpenAI's statement file,
  `UniformFourier.lean`, with only the exponent and seven names changed. The same comparison against a statement
  one unit of 10^-10 stronger fails, as it must.
- The official comparator, the program that OpenAI's repository uses, accepted the proof with both of its kernels
  (nanoda and Lean's) and rejected the statement one unit stronger.
- The agents re-created the leaderboard's verifier offline, on the proof packaged as a submission, and its phases
  passed: the Lean check, which compiled the packaged proof from source, the audit of the axioms, and the nanoda
  kernel, and the same run one unit stronger was rejected. The leaderboard admits at most 500 files, so the
  packaged proof merges small modules (generated data, and chains of modules each importing the one before) into
  larger ones (450 files). A separate check compared the merged form, constant by constant, with the proof built
  from the unmerged files and found the same 29,540 constants in both, each with the same statement and definition,
  and every one of the 12,330 theorems with the same proof.
- The rate fact behind the figure is proved in Lean at 1.7806935e-3, and Lean also proves that it fails at
  every saving one unit higher. Lean's kernel counts the heads in the table of hosts (80,640) and
  refuses one more, and Lean proves that the full cut's profile is wrong with #3945's own head count, with one head
  more, and without the split.
- The agents also re-ran the community's complete checks of each circuit from scratch, with the new tables where
  the tables changed, and they passed.
- Three more agents, each with the sole job of finding a fault, checked the full cut again on other machines: one
  fetched OpenAI's statement afresh, derived the statement files from it once more and rebuilt the proof from the
  package alone; one rebuilt every new module from the package on a fresh copy of the base and ran the official
  comparator again; one re-derived the head count, the split blocks, the separation and the rate fact from the Lean
  data files with programs of its own, checked the hypotheses of Theorem 1G, and re-read the community's figures.
  None found a fault in the proof. They found cosmetic points only, listed with the correct text in the package's `KNOWN-ISSUES.txt`; they are left as
  compiled, because the files are checked byte for byte. Known cosmetic issues: names left from the half-cut package
  that the full cut reuses; in every package's `Main.lean` the comment above the Fourier theorem calls it a sibling of
  a Walsh-Hadamard theorem `wht_main_block_...` that does not exist (it should say only that the theorem follows from
  the package's kernel program); line 15 of design B's `Seam.lean` says 57,312 heads where it should say 80,640; the
  header of `Rate3945F.lean` gives a relative slack of 8.926e-11, which is measured against W (against the right side
  it is 1.07e-10); and the header of each `Goal.lean` cites OpenAI's lines 75-78, 80 and 82 where the definitions are
  at 75-77, 79 and 81. For the half cut on the same table, three other agents did the same; none found a fault in the proof, and one found two comments in Lean files that are out of date (one gives #3945's own head count); they are left as compiled, because the files are checked byte for byte.
- This repository's own scripts (`scripts/`), run as VERIFY.md section 12 says from a fresh clone on a fresh copy
  of Shea's tree, rebuilt `package/leverF` (the full cut), `package/leverB` (the half cut on the same table) and
  `package/leverA` (the lever on #2979 with every head kept), passed every check on all three, and the comparator
  accepted all three again ([scripts/logs-2026-10-19-fresh-clone/](scripts/logs-2026-10-19-fresh-clone/README.md)).
  Preparing that run found three faults in `scripts/build-and-check.sh` with the lever packages; they are fixed.

The limits listed with the result at the top hold for it too, with `package/leverF` in place of `package/leverAFc`.

### The sixth result: 1.6256629e-3 (14 October 2026)

On 14 October 2026 I proved the saving 1.6256629e-3, z = 1 - 16256629/10^10, for the Fourier transform of every
length and for convolution, on the community's U8 circuit. It is the fifth result (below) taken further: the
central totals are read from the sources themselves, so no source copies remain (`CSX/`), and the hosting uses a
bank table with one pad and rate facts with a quartic bound (`Work/BlockAccounting/Moment4`). It stayed unpublished
because it did not beat the community's own figure for that circuit, and it is the engine every lever package is
laid over. The package is `package/u8x0hq`: 257 Lean files, with the Fourier chains of four instances built on
the way (U8XM, U8XMHH, U8X0 and U8X0HH).

Its checks, on Colab machines (the package's `result.json` lists each with its log): the build (299 s, 8.9 GB at
most); the three standard axioms; the host-table checker proved by Lean's kernel to reject a table with one
remainder changed; both statement files OpenAI's file with only the exponent changed (16256629 and 16256630) and a suffix added to seven names; `tools/Compare.lean`
`RESULT: PASS`, and one unit stronger `RESULT: FAIL`; the official comparator in the stand-in sandbox, both kernels
accepting and the stronger statement rejected; the leaderboard's verifier re-created offline (the Lean check
compiled the packaged proof from source in 78 minutes); and an agent whose only job was to refute it rebuilt every
module from a fresh copy of Shea's tree on another machine and re-ran the generators byte for byte, and did not
refute it. On 19 October it was rebuilt on all five machines before the lever build, and this repository's scripts
were run on it from a fresh copy of Shea's tree (`ALL CHECKS PASSED`; the comparator accepted).

One limit of its own: comments in seven of its files describe the U8 circuit with an earlier, inexact list of
makers that names the community's owner maps; those owner maps are not used, and the package's `README.txt` names
the files. The files are left exactly as they were compiled and checked.

### The research since the fourth revision, and the moonshots

[research/](research/README.md) holds everything else I did between 10 and 19 October 2026, in twelve sections,
each with its programs, data and a machine-readable `results.json`. In short:

- **No new unit was found**, at any open shape and in any game, not even at the stepping-stone tier of a saving of
  0.01. The moonshot waves of 14 and 15 October searched blocks 2 to 5 in the engine game, the widened game (every
  address permutation free, Toffoli gates included) and the everything-free game (column-mixing maps and the
  running XOR free as well), and mixed units that finish some arrays only partly. The nearest miss is a unit of 2
  arrays in 7 moves at block 4, allowed by every bound by a margin of 640 out of 137,280; it is still open.
- **What came out instead are lower bounds.** Proved in Lean: the counting rule behind every lower bound, that
  block 2 never saves, and that block 3 needs at least 2.8 moves per array (a block-3 unit saves at most 6.28e-2).
  Checked by computer: block 4 needs at least 41/13 moves per array (at most 0.171). Argued by hand, not checked:
  block 2 never saves in any game whose free steps touch boundedly many entries per row, and the methods behind
  these bounds can never rule out a saving of 0.1 from block 5 on.
- **How the savers save.** Every known saving is frame packing (argued by hand, not checked), and the record
  circuits, OpenAI's, the community's and mine, all save the same way: one free direction per array, delivered
  through a small shared centre. Every lever on that mechanism that the rounds could think of failed to lift the
  record; that is what sent the record track to the recursion instead.
- **Routes, fields and speed.** No route proved that a unit exists without finding it, although such a unit would
  count in full (proved in Lean). The problem is not a famous problem in disguise. On real machines, one result
  moved an application: a number-theoretic transform inside Microsoft SEAL made a loop of encryption operations
  21.2 % faster.
- Failures are told as plainly as successes (92 findings are marked refuted, each with what refuted it), and every
  open question comes with a first experiment ([research/open-problems.md](research/open-problems.md)).

## The results of 13 October 2026 (the fifth result and its rungs)

This part is the text I wrote on 13 and 14 October 2026 for the fifth result, which was not published until
the fifth revision; its figures, checks and credits are unchanged. The sixth result and the lever, above, are built
on it.

**The fifth result.** On 13 October 2026 I proved the saving 1.6198e-3, exactly z = 1 - 16197705/10^10, for the
Fourier transform of every length. OpenAI's statement file also states convolution, and that is proved as well.
Before that day, the largest saving with a Lean proof that I knew of was Ryan Shea's 1.2331685e-3, the top
entry that morning of the leaderboard [bigospeedrun.com/exact-dft](https://bigospeedrun.com/exact-dft). The new
saving is 31.35 % above his (1.3135 times). It is built in his Lean tree, on the community's newer "U8" circuit.
On the way I proved smaller savings the same day, each a rung of the same ladder: 1.3553e-3 and then 1.5699e-3
on the circuit of Shea's proof, Jeremy Kothe's, and 1.3869e-3 on the U8 circuit. The fifth result is 3.17 %
above the best of them. The table "All the results" at the top lists all of them.

**What is new in it.** Three proofs in Lean, each of a step that the community behind the circuit assumes in its
own figures, and the certificate that stacks all three on the community's U8 circuit. A few words first:

- The program works on many *arrays*, working copies of the data. A small circuit of additions runs over and
  over inside it, in a network of nine stages. In each run some arrays are *sources*: they hold what that run
  reads.
- Each array is written in a coordinate system, its *frame*, which grows as the work on the array proceeds. A
  *move* is one costly step on one array: it takes the array to a larger frame. A *climb* is a run of such
  moves. By the end of each run every array must have climbed to the *full frame*, the largest one. Additions
  between arrays are free, but only between arrays that stand at the same frame. Section 5 says more.
- A *helper* is a spare array that holds partial sums. Between stages a helper still has part of its climb
  ahead of it, its *leftover piece*. Shea's pipeline keeps those pieces in *banks*: arrays of their own that
  the program carries through every run. Every bank widens every run, so the fewer banks, the larger the
  saving.
- A *certificate* is a data file that lists the additions of a circuit one by one. A Lean checker reads it, and
  the checker's soundness is proved.
- The *engine* is the part of the Lean proof that turns a checked network into one program for the lengths
  2^k with a time bound. That program is the *kernel program*, and the theorem that delivers it is the *end
  theorem*. The *Fourier chain* is the set of files that carries the kernel program through OpenAI's reduction
  to every length. Only its first file, the *seam*, names the end theorem.
- *The community* is the contributors of the public repository
  [CrocSwap/integer-mult-bounds](https://github.com/CrocSwap/integer-mult-bounds), where the circuit comes from.

The result is four layers, one on top of the other. The first is the community's circuit; each of the other
three is a step that I proved in Lean:

- **The circuit: the community's U8 circuit.** Jeremy Kothe built it in the community's repository, from
  another circuit of his there, published earlier the same day, which carries frame changes by Chafik
  Boukhalfa and is built on icekylinx's construction that the community calls "gauge44". The agents rebuilt it
  with the community's own program, and its checksum matched the one the community records
  (`community/base-logs/u8-word/`).
- **Exit frames, proved in Lean.** The engine checks that a run does the right thing and then cleans up the
  helpers with an *undo step* derived from the run. So every helper climbs to the full frame by the end of a
  run, where the undo step can restore it. The community's newer circuits let some helpers stop
  earlier, at an *exit frame*, a smaller frame, and charge only for the climb that the undo step then needs. I
  proved in Lean that the engine stays correct with exit frames (`XExit`, `Work/KLow`). On the U8 circuit, 126
  helpers end at exit frames. Ryan Shea's notes name exit frames as the next step for the Lean route, and the
  charge for them follows the community's written lemmas about where a helper may end.
- **Mutable sources, proved in Lean.** The engine asked that a run never write into its sources. So Shea's
  conversion first copies each of the circuit's 112 sources into a helper, and every run pays for those
  copies. I proved that the engine stays correct when a run writes into its sources, provided the undo step at
  the end of the run restores them too. This is 7 modules, 1,429
  lines, in `Work/XMut`, and `Work/KMX` joins it with exit frames. Zhihao Chen had written the argument for
  this step in the community's terms; the Lean proof is new. On the U8 circuit five copies stay: five of the
  eight *central totals* (sums formed in the middle of a run that every target reads) sit in source arrays,
  and the engine's step that reads those totals reads them only from helpers.
- **Same-run hosting, proved in Lean.** In each run, some data arrays sit idle for part of the time. DaysSky
  proposed, for the community's circuits, to keep the helpers' leftover pieces in those idle arrays instead of
  in banks ("idle-axis hosting"), and wrote a proof of it in the community's terms. I proved in Lean that the
  engine stays correct when it does so, with every piece hosted by an array that is idle within the same run.
  This is `Work/XHost`, about 5,400 lines. It takes the nine-stage network that Shea's proof runs on, from Nil
  Mamano's files, and adds hosting to it. It adds a checker, `hostSchedK`, that reads the table of hosts and
  banks inside Lean, with a proof that a table it accepts is sound (`HostSched.ofTab_fitsR`). And it ends in an
  end theorem, `XHost.engine_of_ginv_host`, of the same form as the one Shea's proof ends in, so his Fourier
  chain runs on it with only its seam changed. Which parts are adapted from whose files is set out below.

**The certificate that stacks them.** With exit frames and mutable sources, and every leftover piece kept in a
bank, the U8 circuit needs 3,771 banks. With same-run hosting, 2,016 idle arrays host 4,092 pieces in each run,
and 310 banks remain. The agents generated the host table and the Lean data with programs adapted from those
of the hosting result on Kothe's circuit, and Lean's kernel checks the table with `hostSchedK`. The network modules of
`Work/XHost` are the same files, byte for byte, in both hosting results.

**Whose work it rests on.** Most of what this result stands on is other people's, and the credit for it is
theirs:

- **Jeremy Kothe** made the U8 circuit, from another circuit of his and the line of circuits before it. The
  community's notice credits its last step, the "triangle-hybrid" construction, to NewBArch, and the driver
  and composition of that step to ColdWitness/icekylinx.
- **Chafik Boukhalfa** made the frame changes that the U8 circuit carries from that other circuit of Kothe's.
- **icekylinx** made the "gauge44" construction that the U8 circuit is built on.
- **Kothe, DreamingOfClouds and icekylinx**, among others, made the community's helpers that end early, at
  low frames, in the line of circuits the U8 circuit grows from. **Kothe and DreamingOfClouds** wrote the
  community's lemmas about where a helper may end, on which the charge for exit frames rests.
- All of them built on the earlier work that the community's notices name; copies of the notices are in
  `third-party/CrocSwap-integer-mult-bounds/`.
- **Jeremy Kothe** also made the circuit of Shea's proof, "atm3" in the community's repository, on which the
  smaller hosting result is built.
- **Ryan Shea** turned Kothe's circuit into a Lean proof of OpenAI's statement: the conversion pipeline, the
  generators and the Lean instantiation, in
  [shea256/fourier-transform-below-nlogn](https://github.com/shea256/fourier-transform-below-nlogn), folder
  `lean/complex588`. His proof is the record this result improves, and the new files, those of the U8 circuit
  included, are laid over his tree.
- **Nil Mamano** published, in his submission to the leaderboard, the Lean files that Shea's tree takes over
  for the nine-stage network: the extension of the engine for *entrance gauges* (the "gauged" interface
  `GInvE`), the network itself, written by his generator, its bank lemmas and bank-table checker, its end
  theorem `GX.engine_of_ginv_401`, and the templates of the Fourier chain. An entrance gauge lets a helper
  enter a run already part of the way up its climb, so that its leftover piece is smaller. The hosting proof is
  built from these files (next paragraph).
- **Dugongue** made the entrance gauges and the banks of that network; **Zhihao Chen** (jacklightChen) its
  compressed nine-stage layout with its banks and padding; **lydakis** the bank patterns; and **DaysSky** the
  bridged layout it starts from.
- **DaysSky** also proposed hosting helpers in idle arrays, with a written proof that the Lean proof follows
  for hosts within one run, and laid out the integer programme that chooses the hosts.
- **Zhihao Chen** also wrote the argument for mutable sources in the community's terms.
- **The authors of CrocSwap/integer-mult-bounds** built the line of circuits that both circuits belong to.
- **OpenAI** made the cost model and the statement, and the step from the lengths 2^k to every length.

Two parts of the community's own figures for the U8 circuit are not part of this result. One is its
*owner maps*, which say which array keeps each helper's leftover piece: Hao Yan's controller with sixteen
copies of the data, which Kothe applied to the U8 circuit, and the newer owner packings of Chafik Boukhalfa and
Jonathan Capalbo. The other is DreamingOfClouds' bit supplier, the second ingredient of the community's figure
for integer multiplication. My certificates take from the U8 circuit only its additions and their frames; they
keep the leftover pieces in the banks of Shea's pipeline and, for the fifth result, in the agents' host table.

**What the hosting proof takes from whose files.** The hosted network is Mamano's nine-stage network with
hosting added. `HDef`, `HStage`, `HGeom` and `HCert` are copies of its files `GNetDef`, `GNetStage`, `GGeom`
and `GNetCert` with hosting added, and `HEnd` adapts its end theorem (`Work/Bridge401/End.lean` of Shea's
tree), which runs on my engine. The step that moves a helper's leftover piece into its host is `host_adapt`,
`host_pad`, `slot_frames` and `host_frames` in `Work/XHost/Frames.lean`. They are adapted from the bank lemmas
`bank_adapt`, `bank_pad` and `bank_frames` of the same files, with the host treated as a bank whose climb starts
and ends at the host's own frames. The checker `hostSchedK` is modelled on the bank-table checker of those
files and copies its lemmas about lists. Each file says in its first lines what it copies, and
[NOTICE](NOTICE), section 9, lists them.

What is mine:

- the original Lean engine of this repository, which Shea's and Mamano's proofs run on: the *frame lemma*,
  which widens the frames an array may use, and the *generalised engine*, which accepts those frames (section 5
  has both); the certificate checker of the second result; and the lemma `core_nested` that the entrance gauges
  build on;
- the Lean proofs of exit frames (`XExit`, `Work/KLow`), of mutable sources (`Work/XMut`) and of the two
  together (`Work/KMX`);
- for same-run hosting: adding hosting to Mamano's network, the checker `hostSchedK` and its soundness proof,
  the host tables and the instances that use them (`U8XMHH`, and `C588MHH` on Kothe's circuit);
- the generators that wrote the new Lean files, the smaller results of the table below, and the checks.

The agents wrote the proofs under my direction ("Who made it, and how", below).

**Two lemmas that the checked results do not use.** Before the network was built, the agents proved two
lemmas, L1 and L2, as a plan for the proof. They are compiled separately, on a Mac, in
`community/c588mhh/supporting/`, and no module of either hosting result imports them. **L1 is not new**: it restates, in
the form the plan needed, what Shea's tree already proves (`B401.gauged_step`, `B401.gauged_Wn`,
`GX.entW_spec`, `GX.entW_form`, `GX.gauge_step`), on Dugongue's entrance gauges, Mamano's files and my
`core_nested`. L2 shows that a relabelling of addresses, a step the engine counts as free because it is one
pass over the data, can carry a helper's leftover piece into the part of its host's frame that the host has
not yet climbed. The checked result does this with `host_adapt` instead, adapted from Mamano's files as above.

**How it was checked, and the limits.** Every check below passed, on Linux machines of Google's Colab service
(8 processors, 51 GiB):

- Lean's kernel accepted the 32 new modules, compiled from source on the second Colab machine over Shea's tree
  with the files of the 1.3869e-3 result on the U8 circuit. The two theorems depend only on the three standard
  axioms (`propext`, `Classical.choice`, `Quot.sound`).
- `tools/Compare.lean` printed `RESULT: PASS`: the theorem proved is OpenAI's own statement file,
  `UniformFourier.lean`, with only the exponent and seven names changed. The same comparison against a statement
  one unit stronger fails, as it must.
- The official comparator, the program that OpenAI's repository uses, accepted the proof with both of its
  kernels (nanoda and Lean's). It rejected the statement one unit stronger.
- The agents re-created the leaderboard's verifier offline, on the proof packaged as a submission like the
  entries the leaderboard holds, and its three phases passed: the Lean check, which compiled the packaged
  proof from source, the audit of the axioms, and the nanoda kernel.
- A control inside Lean: a copy of the host table with one entry off by one is rejected by the checker.
- Two agents whose only job was to refute the result checked it on their own. One compared the statement with
  OpenAI's file, checked the theorems in Lean in its own copy of the tree and ran the comparator again. The
  other compiled the 32 modules again from the uploaded sources on the first Colab machine, over its own copy
  of an earlier build of the 1.3869e-3 result; all compiled, the comparison passed, the stronger statement
  failed and the axioms were the standard three. Neither found a fault.

The limits are properties of how it was built and checked:

- The proof lives in Shea's Lean tree, not in this repository's own Lean project. This repository carries the
  new files as a package, `community/u8xmhh`, which is laid over his tree after `community/u8xm`; the two need
  one line added to the tree's `lakefile.lean`.
- Every Lean run of the result was on Colab machines. Only a first run of the host-table checker, on the
  generated table and on a copy with one entry off by one, was on a Mac.
- The comparator's sandbox needs a feature of the Linux kernel (Landlock) that those machines lack, so it ran
  inside a stand-in. Its build found the compiled files up to date, so it checked them with both kernels
  without compiling them again; the verifier's re-creation and the second agent's rebuild compiled them from
  source (Mathlib came already compiled).
- The re-creation of the verifier used the agents' own scripts inside the same stand-in. The verifier's own
  controller was never run, and nothing was submitted to the leaderboard.
- The build ran over a second build of the 1.3869e-3 result that the agents made on the same machine; the logs
  of that second build are not in the package, apart from the lines that show its comparator verdicts.
- Some comments inside the new seam file are wrong: they name lemmas of the 1.3869e-3 result that the code
  does not call. The package's `README.txt` lists them, and the rebuilding agent's log shows which lemmas the
  proof does reach.
- Two comments in its new files, and one in the exit-frame files of the U8 results, credit the U8 circuit
  with an earlier, inexact list that names the community's owner maps among its makers. The credits above are
  the right ones; each package's `README.txt` names those comments.
- No person outside this project has reviewed the new proofs, and I would welcome a line-by-line review.

"The results of 13 October in detail", just below, has every check with its log, and [VERIFY.md](VERIFY.md),
section 11, has the commands.

### The results of 13 October in detail

There are eight new things. The first is the result; the others are what it needs and what came with it.

**1. The theorems of the fifth result, and their checks.**

| Theorem | What it is about | Exponent z |
|---|---|---|
| `OAI.PowerSaving.transform_mainU8XMHH` | Fourier transform of every length | 1 - 16197705/10^10 |
| `OAI.PowerSaving.convolution_mainU8XMHH` | convolution | 1 - 16197705/10^10 |

The statement file is `Work/FourierU8XMHH/UniformFourierChallenge.lean`. It is OpenAI's `UniformFourier.lean`
with the exponent `1 - 1/10^13` replaced by `1 - 16197705/10^10` and the suffix `U8XMHH` put on the seven names
that carry the exponent (`decimalExponent`, `decimalTime`, `TimeBounds`, `DFTGoal`, `ConvGoal`,
`transform_main`, `convolution_main`). The rest is OpenAI's.

The checks, on the second Colab machine unless the line says otherwise:

| Check | Result | Log, under `community/u8xmhh/` |
|---|---|---|
| build: the 32 new modules with Lake, in six steps, over Shea's tree with `community/u8xm` built | every step ended without error (15 minutes, 8.6 GB at most); then every target up to date | `logs/colab/build.out` |
| axioms of the two theorems and of the program behind them | `propext`, `Classical.choice`, `Quot.sound` only | `logs/colab/axioms.log`, and `logs/colab/axioms-glue.log` for the host-table check and the end theorem |
| `tools/Compare.lean` against OpenAI's statement with the exponent of the table above | `RESULT: PASS` | `logs/colab/compare-pos.log` |
| the same against the statement one unit of 10^-10 stronger | `RESULT: FAIL`, as it must | `logs/colab/compare-neg.log` |
| the two statement files against the agents' own derivation from OpenAI's file | the same, byte for byte | `logs/colab/build.out` |
| a copy of the host table with one entry off by one | Lean's kernel proves that `hostSchedK` rejects it | `logs/colab/hsneg.time`, `logs/colab/ax/HsU8XMHHNeg.lean` |
| the sources scanned for `sorry`, `admit` and other forbidden words | only the placeholders of the statement files | `logs/colab/sorry-grep.txt`, `admit-grep.txt`, `forbidden-grep.txt` |
| official comparator, with a stand-in for its sandbox; its build found the compiled files up to date | nanoda and Lean's kernel accept: `Your solution is okay!` (15 minutes, 7.0 GB) | `logs/colab/cmp/fsb-U8XMHH-pos.log`, `logs/colab/cmp.out` |
| official comparator, one unit stronger | rejected | `logs/colab/cmp/fsb-U8XMHH-neg.log` |
| the leaderboard's verifier, re-created offline with the agents' scripts in the same stand-in, on the proof packaged as a submission | its three phases pass: the Lean check, compiling the packaged proof from source (68 minutes, 8.4 GB); the axiom audit; nanoda and Lean's kernel (14 minutes) | `logs/colab/bos-H.out`, `logs/colab/bos-H/` |
| Python checks on the Mac before the build: an independent audit of the generated tables, the placement and the rate | pass | `logs/mac-preflight/` |
| static checks on the Mac (`tools/community/package_check.py`), as for every package | pass | `logs/static-check.txt` |
| an agent checks the statement on its own, in its own copy of the tree on the same machine: OpenAI's files, the differences of the statement files, the theorems in plain Lean, the comparison and the comparator again | all as above; the comparator accepted, and rejected one unit stronger | `logs/skeptic-statement/` (`report.txt` is its written verdict; `logs/out/check.log`, `logs/cmp/`) |
| an agent compiles the 32 modules again from the uploaded sources on the first Colab machine, one per step, over its own copy of an earlier build of `community/u8xm` | all compiled; `RESULT: PASS`, one unit stronger `RESULT: FAIL`; standard axioms; the proof reaches the end theorem of the hosted network and the host-table checker; it did not run the comparator | `logs/skeptic-rebuild/` (`report.txt`; `logs/build.out`, `compare-pos.log`, `compare-neg.log`, `uses.log`) |
| `community/u8xm`, laid under it, built a second time on the second machine | its comparator accepted it, and rejected one unit stronger | `logs/colab/setup.out` (first lines); the logs of that build are not shown |

**The same two proofs on Kothe's circuit** (`community/c588mhh`, 1.5699305e-3). The first hosting result, built
on the first Colab machine the same afternoon, has these theorems:

| Theorem | What it is about | Exponent z |
|---|---|---|
| `OAI.PowerSaving.transform_mainC588MHH` | Fourier transform of every length | 1 - 15699305/10^10 |
| `OAI.PowerSaving.convolution_mainC588MHH` | convolution | 1 - 15699305/10^10 |
| `OAI.PowerSaving.WHT.wht_main_block_C588MHHx` | Walsh-Hadamard transform, lengths 2^k | 1 - 15699306/10^10 |

The statement file of the first two is `Work/FourierC588MHH/UniformFourierChallenge.lean`, made the same way
with the suffix `C588MHH`. The third theorem is stated in the form of this repository's published
Walsh-Hadamard statement, with its exponent and its name changed. Its proof is one term: my theorem
`BlockWHT.wht_main_of_engine` applied to the kernel program of the first two.

The checks of the first two, on the first Colab machine unless the line says otherwise:

| Check | Result | Log, under `community/c588mhh/` |
|---|---|---|
| build: the 49 new modules, one at a time with plain `lean`, over Shea's tree with `community/c588m` built | the last compile of every module ended without error (earlier rounds, while the agents fixed errors, are kept) | `logs/colab/r2.summary` to `logs/colab/r14.summary` |
| axioms of the two theorems and of the program behind them | `propext`, `Classical.choice`, `Quot.sound` only | `logs/colab/r14-Work.FourierC588MHH.Axioms.log` |
| `tools/Compare.lean` against OpenAI's statement with the exponent of the table above | `RESULT: PASS` | `logs/colab/final-C588MHH-compare-pos.log` |
| the same against the statement one unit of 10^-10 stronger | `RESULT: FAIL`, as it must | `logs/colab/final-C588MHH-compare-neg.log` |
| official comparator, with a stand-in for its sandbox; its own Lake build compiled the hosted modules again | nanoda and Lean's kernel accept: `Your solution is okay!` (14 minutes, 8.6 GB) | `logs/colab/final-C588MHH-cmp-pos.out` |
| official comparator, one unit stronger | rejected | `logs/colab/final-C588MHH-cmp-neg.out` |
| control: the hosted network with an empty host table | `RESULT: PASS` at 1 - 13552606/10^10, the figure without hosting; one unit stronger fails; comparator accepts, and rejects one unit stronger | `logs/colab/final-C588MH.out`, `logs/colab/cmp-C588MH.out` |
| static checks on the Mac (`tools/community/package_check.py`): checksums; the statement against OpenAI's, with the notice of lines 3 to 12 removed; the notice present; the negative statement; the comparator configurations; the module list; no forbidden word in the Lean code; every check of `result.json` either passed with a log inside the package or marked "not shown" or "not run" | pass | `logs/static-check.txt` |
| an agent copies the sources afresh and compiles the 32 modules of the result into its own build folder | all compiled; `RESULT: PASS`, one unit stronger `RESULT: FAIL`; standard axioms, also for `XHost.engine_of_ginv_host`; it did not run the comparator | `logs/rebuild-skeptic/a.summary`, `cmp-pos.log`, `cmp-neg.log`, `ax.log` |
| the lemmas L1 and L2, which the result does not import, compiled on the Mac | no `sorry`, standard axioms | `supporting/L1/logs/`, `supporting/L2/logs/` |
| two agents, one per lemma, rebuild L1 and L2 on the Mac from their own copies, with the files of Shea's tree they need taken from his tree | all compiled, standard axioms; for L2 also an audit module of its own, for L1 a link to the network's windows and an exact check of the seam | `supporting/L1/skeptic/`, `supporting/L2/skeptic/` (`report.txt` in each is its written verdict) |
| negative control of the host-table checker inside Lean | ran on the Mac on an earlier version of the data files | not shown |
| the leaderboard's verifier (bigospeedrun/verifier) | not run on this result | not run |

The Walsh-Hadamard theorem was built and checked on the second Colab machine: 21 modules compiled, standard
axioms, `tools/Compare.lean` (and its failure one unit stronger), and the official comparator with both kernels
(and its rejection one unit stronger). Its logs are in `community/wht-c588mhh/logs/colab/` (`run.out` prints
`ALL CHECKS PASSED`; `cmp.out` has the comparator runs).

**2. Same-run hosting: what the proof says.** A run of the circuit is a product of simple matrices that act on
the arrays, in nine stages. Between stages, each helper still has a leftover piece of its climb to the full
frame. Shea's pipeline gives every such piece a place in a bank. But in each run, the data arrays sit idle
during some of the stages: the nine-stage network in Shea's tree has 18 such idle stretches, and at the end of
each the idle array still has a climb of its own to make. Hosting keeps a helper's leftover piece inside such
an idle array instead of in a bank. A relabelling of addresses, which the engine counts as free, moves the
piece into the part of the host's frame that the host has not yet climbed. The host's own climb then finishes
with the piece inside it, and only the rest of that climb is paid. In the proof this is `host_adapt` and
`host_pad` (`Work/XHost/Frames.lean`), adapted from Mamano's bank lemmas `bank_adapt` and `bank_pad`: a host is
handled as a bank whose climb starts and ends at the host's own frames.

The Lean proof in `Work/XHost` has these parts:

- the roles of the arrays in a hosted stage, and a proof that no two pieces are given the same place
  (`HRoles`); the hosted stage and its certificate (`HDef`, `HStage`, `HCert`), which are copies of Mamano's
  `GNetDef`, `GNetStage` and `GNetCert` with hosting added; the frames of a host (`Frames`), adapted from his
  bank lemmas;
- the geometry of the 18 idle stretches (`HGeom`, a copy of the instance part of Mamano's `GGeom` with hosting
  added);
- host schedules (`HSched`) and the checker `hostSchedK` (`HSchedK`), which reads the host table and the bank
  table inside Lean's kernel. It is modelled on the bank-table checker of Mamano's files and copies their
  lemmas about lists. Its soundness, `HostSched.ofTab_fitsR`, turns any table it accepts into a proved
  schedule, so a wrong table cannot pass;
- the end theorem `XHost.engine_of_ginv_host` (`HEnd`). It has the same conclusion as Mamano's
  `GX.engine_of_ginv_401`, and its file is adapted from his `Work/Bridge401/End.lean` in Shea's tree, which
  is built on my engine. So the Fourier chain of `community/c588m` runs on it with only its seam changed
  (`Work/FourierC588MHH`);
- the instances. `U8XMHH`, the fifth result, is the certificate of `community/u8xm` with its host table; its
  data, glue and Fourier chain are in `community/u8xmhh`, whose twelve network modules are the same files as
  here. `C588MHH` is the certificate of `community/c588m` with its host table, and the control `C588MH` is the
  same network with an empty host table.

*Same-run* means that every piece is hosted by an array that is idle within the same run. Letting a piece wait
in an array that is idle during the next run would give a little more and needs a change in the order of the
network (open directions).

What the two proofs do to Kothe's circuit, per run:

| | Shea's conversion | mutable sources (`community/c588m`) | and same-run hosting (`community/c588mhh`) |
|---|---|---|---|
| sources | 112 | 112 | 112 |
| helpers | 616 (112 of them copies of sources) | 504 | 504 |
| banks | 5,040 | 4,032 | 111 |
| saving, Fourier transform of every length | 1.2331685e-3 | 1.3552606e-3 | 1.5699305e-3 |

And on the U8 circuit, per run:

| | exit frames and mutable sources (`community/u8xm`) | and same-run hosting (`community/u8xmhh`) |
|---|---|---|
| helpers | 509 (504, and the 5 copies of sources that stay) | 509 |
| hosted pieces | none | 4,092, in 2,016 idle arrays |
| banks | 3,771 | 310 |
| saving, Fourier transform of every length | 1.3868856e-3 | 1.6197705e-3 |

The community's repository lists "whole-call hosting" among the steps its headline figure assumes. The Lean
proof here is of hosting within one run, for this engine. It does not by itself check the community's own
construction.

**3. Mutable sources: what the proof says** (the package `community/c588m`, 1.3552606e-3). The engine checks
that a run does the right thing and then undoes the run's side effects with an undo step derived from it. Before,
it proved this under the rule that no step of the run writes a source, so that the sources come out of the run
unchanged. The new modules drop that rule. The undo step then restores the sources together with the helpers.
The proof needs only that no output feeds back into a source or a helper and that the run is invertible. It also
checks that the undo step is legal where it stands: at the end of a run, the sources and the helpers it touches
all stand at the same frame. The theorems are `XMut.dirt_fwdZ` and `XMut.dirt_bwdZ`
(`Work/XMut/Dirt.lean`); the other six modules carry them through the checker and into Mamano's gauged
interface.

Without the copies, Kothe's circuit has 5,916 moves per run instead of 6,924. It breaks the old rule in 535
places, by the exact Python check, and the old check inside Lean rejects it too (`community/c588m`, the negative
controls), so the gain comes from the new proof and from nothing else. The package passed the same kinds of
check as the hosting results, with its logs in `community/c588m/logs/colab/` (`mb-build.out`, `axioms.log`,
`compare-fourier.log`, `compare-negative.log`, `cmp-A1.out`, `cmp-A2neg.out`,
`68-Work.GCert.Data.C588MControls.log`, `source-scan.txt`).

It is also, with the fifth result, one of the two packages on which the agents re-created the leaderboard's
verifier. On a submission packaged like the entries the leaderboard already holds, they ran the verifier's three
phases offline (the comparator's Lean check, the audit of the axioms and the nanoda kernel) with scripts of their
own, inside the same stand-in sandbox as the comparator, and all three passed (`logs/colab/bos/`). The verifier's
own controller, which needs sandboxing tools those machines lack, was not run. Two agents also tried to refute the
result. One rebuilt all 70 modules on the second Colab machine against that machine's own build of Shea's
tree; every module compiled, the comparison passed on its build and failed one unit stronger, and the axioms
were the standard three (`logs/skeptic-rebuild/`); it did not run the comparator. The other compared the
statement files with OpenAI's file and read the logs; it did not run Lean, and its written verdict is the only
record of its check (`logs/skeptic-statement/report.txt`).

**4. The supporting results, with their checks and their gaps.** Every package has a `result.json` that lists
each check with its status. A check marked "pass" names its log inside the package; any other check is marked
"not shown" (it ran, but its log is not in the package or has no verdict) or "not run". The static checks of
`tools/community/package_check.py` refuse a package that breaks this rule.

| Package | What is proved | Saving | Passed, with logs | Not done or not shown |
|---|---|---|---|---|
| `community/c588mhh` | Fourier transform of every length: Kothe's circuit with mutable sources and same-run hosting | 1.5699305e-3 | all of the second check table of item 1, with its logs | not shown: a negative control of the host-table checker on the Mac; not run: the leaderboard's verifier |
| `community/wht-c588mhh` | Walsh-Hadamard transform on the kernel program of `c588mhh` | 1.5699306e-3 | build, axioms, comparison and its failure one unit stronger, comparator accepted and rejected one unit stronger (`logs/colab/`) | not run: an independent check; the leaderboard's verifier |
| `community/u8xm` | Fourier transform of every length: the U8 circuit with exit frames and mutable sources (5 of the 112 copies kept) | 1.3868856e-3 | build, axioms, comparison and its failure, comparator accepted and rejected, negative controls in Lean (`logs/colab/lake.out`, `checks2.out`, `cmp/`) | not shown: a second build by the agents and two independent agents' checks of it, which found no fault (their logs are not in the package); not run: the leaderboard's verifier |
| `community/c588m` | Fourier transform of every length: Kothe's circuit with mutable sources | 1.3552606e-3 | all of item 3: the checks of the build, the re-created verifier (`logs/colab/bos/`), an agent's rebuild of the 70 modules with the comparison and its failure (`logs/skeptic-rebuild/`), and an agent's written comparison of the statement (`logs/skeptic-statement/report.txt`) | not run: the verifier's own controller; the comparator on the agent's rebuild |
| `community/wht-c588` | Walsh-Hadamard transform on the kernel program of `c588m` | 1.3552607e-3 | build, axioms, comparison and its failure, comparator accepted and rejected (`logs/colab/c588m/`) | not run: an independent check; the leaderboard's verifier |
| `community/u8x` | Fourier transform of every length: the U8 circuit with exit frames | 1.2704023e-3 | build, axioms, comparison and its failure, comparator accepted and rejected (`logs/colab/direct.out`, `checks.out`, `cmp/`); two controls in Lean (`controls.log`; their file was compiled on the Colab and is not in the package) | not shown: a second build by the agents and two independent agents' checks of it, which found no fault (their logs are not in the package); not run: the leaderboard's verifier |
| `community/c588-repack` | Fourier transform of every length: Shea's certificate with a repacked bank table | 1.2383291e-3 | build, axioms, comparison and its failure, comparator accepted and rejected (`logs/colab/`); an agent's re-check of the built tree, mounted read-only: up to date, standard axioms, the exponent, the comparison (`logs/colab/skeptic-v-c-pack/vcp2.log`, `compare-pos.log`) | not shown: that agent's comparison one unit stronger, which did not finish; not run: the leaderboard's verifier |
| `community/wht-c588` | Walsh-Hadamard transform on Shea's kernel program and on the repacked one | 1.2331686e-3 and 1.2383292e-3 | build, axioms, comparison and its failure, comparator accepted and rejected (`logs/colab/`); on Shea's kernel program, an agent's rebuild with the comparison and its failure (`logs/colab/skeptic-v-c-wht/base.log`) | not shown: that agent's check on the repacked kernel program, stopped at its time limit; not run: the leaderboard's verifier |

- **Exit frames** (`u8x`, `u8xm`, and under the fifth result). In the engine every helper climbs to the full frame at the end of a run, so
  that the undo step can restore it. The community's newer circuits let some helpers stop earlier, at an *exit
  frame* of lower dimension, and charge for the climb that the undo step then needs. The agents wrote that theory
  in Lean (`XExit`, `Work/KLow`) and applied it to the community's U8 circuit: 126 helpers end at exit frames of
  dimension 4 or 8. `Work/KMX` joins it with mutable sources. Five source copies stay. In the U8 circuit, five
  of the eight *central totals* (sums formed in the middle of a run that every target reads) sit in source
  arrays, and the engine's step that reads those totals reads them only from helpers. Ryan Shea's notes name
  exit frames as the next step for the Lean route, and the charge for them follows the community's written
  lemmas about where a helper may end. These two packages, and the fifth result laid over `u8xm`, need one line
  added to the lakefile of Shea's tree (VERIFY.md, 11.2c).
- **The repacked bank table** (`c588-repack`). The circuit keeps its 112 copies. At the end of every unit each
  bank climbs whatever room its pieces did not fill, and that *padding* is paid. An integer programme found a
  table with 5,009 banks whose padding costs the least possible, by the agents' exact count, and a sharper
  inequality in the step that turns the cost of one unit into the exponent (`term_le4`, in a new module
  `Moment4.lean`) takes the figure to its exact limit for that table. It replaces 38 of Shea's generated files,
  so it needs its own copy of his tree.
- **The Walsh-Hadamard corollaries** (`wht-c588`, `wht-c588mhh`). Each is the one-term corollary of item 1 on
  another kernel program. None has been built on the kernel program of the fifth result.

**5. The programs, and a script to swap the certificate.** `tools/community/` holds the programs that generated
the packages of the mutable-sources result and their inputs: `c588m-gen` regenerates all 147 Lean and
comparator files of `community/c588m` byte for byte (`REPRODUCE.sh`, a few seconds on a Mac); `mutable` re-derives
the circuit without the copies from Shea's generated data, with its checks and controls; `repack` remakes the
bank table of item 4. The programs that wrote the hosting files and chose the host tables are inside the
packages, in `community/c588mhh/generators/` and `community/u8xmhh/generators/`; regenerating the packages from
them was not re-run. A better certificate in
the same form can replace the headline one with `tools/swap-certificate.sh`, which checks the new package,
installs it next to the old one and changes the text only when every check has passed with a log (VERIFY.md,
11.5). I moved the headline from `community/c588mhh` to `community/u8xmhh` by hand: the script would also have
renamed the many places that name `community/c588mhh` itself, which stays.

**6. The account of the day, and the open directions.**
[notes/fifth-result-story.md](notes/fifth-result-story.md) tells how the result was found: a survey of the
community's work, a ladder of steps that a Lean proof must pay for and the community's figures assume, the
three proofs that climbed it, and what failed and why, including the claims of the day that turned out wrong.
It keeps the figures of the afternoon, when the fifth result was the hosting result on Kothe's circuit; a note at
its top says how the full stack was built that night.
[notes/open-directions-3.md](notes/open-directions-3.md) lists what may come next, each with a first experiment.

**7. What is not claimed for the fifth result.**

- **The circuit is not mine, and neither are the ideas of the three steps.** The U8 circuit is the
  community's, put into Shea's Lean tree; exit frames are the community's, with their charge from
  its written lemmas; mutable sources had Zhihao Chen's written argument; hosting is DaysSky's idea. What is
  new is the Lean proofs of the three steps, and the certificate that uses them.
- **It is not the community's figure, and it is below it.** The community's figure has two parts, and the one
  of the same kind as my Fourier saving is its *complex side*. Its selected figures on 13 October 2026 had
  complex sides of 1.6644e-3, for the circuit of Kothe's that the U8 circuit is built from, and then
  1.6647e-3, for the U8 circuit with Capalbo's owner packing; its newest circuits were priced higher still,
  around 1.675e-3. Its own price of Kothe's earlier circuit has a complex side of 1.6049e-3. Those figures
  assume several steps that have no Lean proof; exit frames, mutable sources and hosting within one run were
  three of them. The fifth result is compared only with figures that have a Lean proof.
- **Nothing here is about integer multiplication.** Section 2, point 1, holds for this result as well.
- **It has the limits of every result here:** it is a statement about growth (section 2, point 5), and its
  Fourier saving sits 1/10^10 below that of its kernel program (point 4).
- **The builds ran over trees that already held compiled modules**: Shea's, and for the fifth result those of
  `community/u8xm`, from a second build of it whose own logs are not in the package. The offline re-creations
  of the leaderboard's verifier, on `community/c588m` and on the fifth result, and the independent rebuild of
  the fifth result compiled the new modules from source; Mathlib came already compiled.
- **The figures of later steps are prices or what-ifs.** The figures in the open-directions note (1.6256506e-3
  without the last five source copies and 1.6293239e-3 with hosting across runs, on the U8 circuit; 1.5785076e-3
  with hosting across runs on Kothe's circuit) are computed by Python programs for constructions whose Lean
  proofs are not written.

**8. Whose work it is, once more, in one table.** Section 9, "For the fifth result", gives the pull request or
file behind each line.

| Who | What |
|---|---|
| Jeremy Kothe | the U8 circuit and the line of circuits it grows from; the circuit of Shea's proof ("atm3"), with the predecessors its notices name; with DreamingOfClouds, the lemmas about where a helper may end, on which the charge for exit frames rests |
| Chafik Boukhalfa | the frame changes that the U8 circuit carries |
| icekylinx | the "gauge44" construction that the U8 circuit is built on |
| DreamingOfClouds | with Kothe, icekylinx and others, the helpers that end early in the line of the U8 circuit; with Kothe, the lemmas about where a helper may end |
| NewBArch, ColdWitness/icekylinx | the triangle-hybrid step of the U8 circuit, as the community's notice credits it |
| Ryan Shea | the conversion pipeline, the generators, the bank and seam pipeline and the Lean instantiation of Kothe's circuit, whose tree the U8 results are laid over too; the lemmas that L1 restates; the notes that name exit frames as the next step for the Lean route |
| Nil Mamano | the Lean files of the nine-stage network that Shea's tree takes over from his submission to the leaderboard: the gauged extension of the engine (`GInvE`), the network (`GNetDef`, `GNetStage`, `GGeom`, `GNetCert`), the bank lemmas (`bank_adapt`, `bank_pad`, `bank_frames`), the bank-table checker, the end theorem `GX.engine_of_ginv_401`, and the Fourier chain templates. The hosted network copies or adapts these files |
| Dugongue | the entrance gauges and the banks of that network |
| Zhihao Chen (jacklightChen) | the compressed nine-stage layout with its banks and padding; the written argument for mutable sources in the community's terms |
| lydakis | the bank patterns |
| DaysSky | the bridged layout the network starts from; hosting helpers in idle arrays, with the written proof that the Lean proof follows for hosts within one run; the layout of the integer programme that chooses the hosts |
| the authors of CrocSwap/integer-mult-bounds | the line of circuits both circuits belong to, and the earlier work their notices name |
| OpenAI | the cost model, the statements and the step to every length |
| me, with the agents | the original Lean engine these proofs run on; the Lean proofs of exit frames, of mutable sources and of the two together; for same-run hosting, adding hosting to Mamano's network, the checker `hostSchedK` and its soundness proof, the host tables and the instances (`U8XMHH`, `C588MHH`); the generators, the certificates, the repacked bank table, the Walsh-Hadamard corollaries and the checks |

Not in this result, though the community's own figures for the U8 circuit use them: the owner maps of Hao Yan,
Kothe, Chafik Boukhalfa and Jonathan Capalbo, and DreamingOfClouds' bit supplier ("Whose work it rests on",
above).

## What was new in the fourth revision (2026-10-10)

**The fourth result.** The saving proved on 10 October 2026 was 8.7625e-4, for two statements:

- the Fourier transform of every length, with z = 1 - 8762478/10^10; the file that states it also states
  convolution, and that is proved as well;
- the Walsh-Hadamard transform of the lengths n = 2^k, with z = 1 - 8762479/10^10. That is 1.1723 times the
  7.4745e-4 of my second result (+17.2 %).

Both are proved in this repository's own Lean project. The file that states the Fourier theorem is OpenAI's own
statement file, which they call a challenge file, with the exponent replaced (section 1 names the two other
differences, a comment and a suffix on seven names).

**Where the improvement came from.** OpenAI's proof reduces the saving to a finite network: a fixed number of
working copies of the data (*arrays*), each to be transformed along a fixed number of directions, with a
schedule that needs fewer steps than the plain one. I call one such network a *unit*, and the program for the
lengths 2^k that is built on it the *kernel program*. The fourth result came from a new and much smaller unit,
found by computer search on 10 October 2026. It has 1,263 arrays with 45 directions each (in the words used
below, a label space of dimension 45). The unit of the second and third results has 14,692 arrays with 110
directions each. The Fourier theorem is the same kernel program carried to every length; its saving is
smaller by 1/10^10.

**Whose work it rests on.** The step from the lengths 2^k to every length is OpenAI's reduction, repeated here
with one constant changed. The helper circuit of the second and third results, that is the scheme of additions
between their arrays, was made by the contributors of the public repository CrocSwap/integer-mult-bounds, whom
I call *the community*. Three devices in the circuit of the fourth result are other people's too. The sections
"Whose ideas the new unit uses", "Whose reduction this is" and "Whose circuit this is" give the credit, and
[RELATED-WORK.md](RELATED-WORK.md) has the credits and the relation to other work.

**Where to look.** Section 1 has the exact statements. Section 2 has what is not claimed. Section 3 and
[VERIFY.md](VERIFY.md) have the checks and how to repeat them. Section 9 has the dates, commits and outside
sources behind the statements, which are kept out of the prose.

There are five new things.

**1. A machine-checked saving of 8.7625e-4 for the Walsh-Hadamard transform, and for the Fourier transform of
every length and for convolution.** The three theorems are in this repository:

| Theorem | What it is about | Exponent z |
|---|---|---|
| `wht_main_block_B2Ge8x` | Walsh-Hadamard transform | 1 - 8762479/10^10 |
| `transform_mainE8` | Fourier transform of every length | 1 - 8762478/10^10 |
| `convolution_mainE8` | convolution | 1 - 8762478/10^10 |

All three are checked in three ways:

- by Lean's kernel;
- by the official comparator, the program that OpenAI's repository uses for its own results. It tests that a
  proof proves the statement of a given challenge file;
- by `tools/Compare.lean`, a second comparison script, which is in this repository.

All three checks were run on a Mac and on a Linux machine, and all passed on both. On both machines the 39
new Lean modules were compiled from source, in a tree that already held the compiled modules of the earlier
results. On both, the comparator ran with a stand-in for its sandbox, not with the real one. Section 3 has
the runs and their limits.

As of 10 October 2026 this was the largest saving with a Lean proof that I knew of, for either statement. Here
is what that rests on.

- On that day the agents read the titles and descriptions of the newest pull requests of the community's
  repository, CrocSwap/integer-mult-bounds. Its contributors work on bounds for integer multiplication. One
  part of their figure, the "complex side", is the quantity of my Walsh-Hadamard theorem.
- The largest figure there for that quantity with a reported Lean check was 7.795e-4, by DreamingOfClouds
  ([pull request 327](https://github.com/CrocSwap/integer-mult-bounds/pull/327)). It is for the five-stage
  layout of this repository, the arrangement of the network in five stages that this repository uses
  (section 6). Its author reports a Lean kernel check of it, made with the
  generators of this repository (the Python programs that write its Lean files). It was not built here.
- The largest kappa in a title there was 7.77596e-4, by eumemic
  ([pull request 334](https://github.com/CrocSwap/integer-mult-bounds/pull/334)). Kappa is the saving for
  integer multiplication, which is another quantity (section 2, point 2).
- The claim covers figures with a Lean proof only. Prices without one, in that repository or elsewhere, are
  outside it. (A *price* is the saving that a program computes for a network. Without a Lean proof it is a
  computed number, not a theorem.)
- Those figures move within hours.
- For the Fourier statement the agents read the outside repositories again on that day. The largest outside
  figure with a Lean proof was 7.547360e-4. It is in pull request #2 of this repository, by Chafik Boukhalfa
  (the account chafreaky).
- That figure is his result. It is the circuit of his pull request 233 in the community's repository, in the
  five-stage layout of this repository. He reports the Lean build and the comparator runs.
- On 10 October 2026 the agents built his pull request on a Linux machine. All 54 Lean modules that it adds
  compiled. His three theorems depend only on the three standard axioms. `tools/Compare.lean` printed
  `RESULT: PASS` for his Walsh-Hadamard statement and for his Fourier statements. After the build, the
  agents ran the official comparator on the two configurations for which he reports runs, with a stand-in
  for its sandbox. It accepted both. He reports his own comparator runs, made with the real sandbox. On my
  Mac the build could not finish, because two of his modules need about 13 GB of memory each. Section 9 has
  the details.
- shea256's repository now gives that same figure, where it had 6.7e-4 on 9 October. It reports its own
  rebuild of that proof and says that it "does not claim a new exponent of our own". That rebuild is
  reported by its author. It was not repeated here.

Section 9 lists what was read, and when.

**2. A unit of the kernel program that is an order of magnitude smaller.** Six words first.

- A *bank pair* is two data arrays that exchange their contents through additions. The exchanges are what
  saves steps.
- A *helper array* is an extra array that holds partial sums. Helpers save nothing, so their number is the
  main cost.
- A *move* transforms one array along one direction, and is what costs. A *block* is a run of moves of one
  array.
- The *label* of a pair is a word of bits, 9 bits wide in the new unit. The *label space* is the space of the
  m directions along which every array has to be transformed.

| | New unit | Unit of the second and third results |
|---|---|---|
| bank pairs | 120 | 1,320 |
| width of their labels | 9 | 22 |
| helper arrays | 783 | 9,412 |
| arrays in all | 1,263 | 14,692 |
| dimension m of the label space | 45 | 110 |
| one pass of the circuit | 9,039 moves in 3,974 blocks | 262,944 moves in 70,169 blocks |

Section 6 has the full table.

**3. The new unit needed no new Lean proof for the Walsh-Hadamard theorem.** The second result came with a
*certificate* format and a checker. A certificate is a data file that lists the additions of a circuit one by
one. The checker is a Lean function, with proofs that what it accepts is sound. The circuit of the new unit is
one more data file for that checker: `tools/certificate/gcert1-e8-r783.json.gz`. The Lean modules that hold
and check it were written from that file by the two generators of this repository, the Python programs
`tools/gx/gxgen.py` and `tools/gx/gxrate.py`, which are unchanged. The new circuit does use parts of the
certificate format that the
certificate of the second result never used (section 5).

**4. The search programs, and the account of the search.** `tools/e8/` holds the programs that found the
certificate. Four notes say how it went, and a fifth is for readers from another field. They are written for
people and for AI agents who want to take this further.

- [notes/how-the-e8-unit-was-found.md](notes/how-the-e8-unit-was-found.md): the design, and each step with its
  figure.
- [notes/what-did-not-work.md](notes/what-did-not-work.md): the routes that failed, told as plainly as the one
  that worked.
- [notes/e8-results.json](notes/e8-results.json): every figure in machine-readable form, each marked as
  Lean-checked here or as a price computed in Python.
- [notes/open-directions-2.md](notes/open-directions-2.md): what has not been tried and may be promising.
- [notes/for-network-coding-readers.md](notes/for-network-coding-readers.md) describes, for readers from
  network coding, the small coding problems that sit behind schedules of this kind, with what is checked,
  what is open and what I do not claim; the reading it starts from is that of pull request #288 of
  CrocSwap/integer-mult-bounds (eumemic), not mine, and nothing in that note is checked by a proof assistant.

**5. The scripts that rebuild the certificate of the second result from the community's data, and a
correction.** The scripts are in `tools/rebuild/`. Issue #1 of this repository asked for them. The correction:
six of those scripts are adaptations of the community's own scripts, and the earlier text called the generator
the agents' own (section "Whose circuit this is"; [NOTICE](NOTICE), section 6).

**Whose ideas the new unit uses.** The circuit of the fourth result is not the community's circuit of the
second result: the agents I directed designed it and found it by search. Three of the devices in it are other
people's, and the saving would be smaller without them:

- **The in-place pair.** Two arrays hold a and b, and a + b and a - b are formed on those same two arrays. It
  comes from icekylinx ([pull request 184](https://github.com/CrocSwap/integer-mult-bounds/pull/184) of the
  community's repository), and two pull requests by ikeboy carry it. It is in the circuit of the second result
  too.
- **The re-use of finished helpers.** A helper that nothing reads any more takes a new value, and each target
  concerned first reads the old content with the opposite sign. It comes from jamesyc
  ([pull request 124](https://github.com/CrocSwap/integer-mult-bounds/pull/124)) and eumemic
  ([pull request 143](https://github.com/CrocSwap/integer-mult-bounds/pull/143)).
- **The late pairing rule.** This is the choice of which finished helper a new value takes. It comes from
  Chafik Boukhalfa, the account chafreaky (pull requests
  [200](https://github.com/CrocSwap/integer-mult-bounds/pull/200) and
  [233](https://github.com/CrocSwap/integer-mult-bounds/pull/233) there, and pull request #2 of this
  repository). The search programs here apply his idea with a greedy rule of their own. His program, his pairs
  and his certificate were not used.

Which partial sums the circuit shares was chosen by an integer programme that the agents wrote and solved with
the open-source solver HiGHS. The Fourier and convolution theorems rest, as in the third result, on OpenAI's
reduction (section "Whose reduction this is").

What is mine in the fourth result:

- the new unit: the family of 120 labels of width 9 as a unit of this mechanism, the circuit on it, and its
  certificate;
- the search that found it, whose direction I chose, and the search programs in `tools/e8/`, which the
  agents wrote;
- from the second result: the certificate format, the checker with its soundness proofs, the generators, the
  five-stage layout and the generalised engine (the engine is the part of the proof that turns a schedule
  into one program with a time bound; section 5);
- the machine check of the whole.

[RELATED-WORK.md](RELATED-WORK.md), section 8, has the fuller credit, with what was and was not read at the
community's repository. The pull requests named here are listed once more in section 9.

## What was new in the first three revisions

All three revisions were published here on 9 October 2026. The searches of outside work behind this list were
made on that day.

- **Third result: a machine-checked proof of OpenAI's own headline statement, with the saving 7.474546e-4.**
  The statement is the discrete Fourier transform of every length. OpenAI's statement has 10^-13, the
  fourth result has 8.7625e-4 and the fifth 1.6198e-3. The statement file is OpenAI's challenge file with
  the exponent replaced (section 1). As far as I could find, the largest outside figures were these:
  - in an outside Lean proof of that statement: 3.2e-6
    ([danadran01/exact-dft-power-saving](https://github.com/danadran01/exact-dft-power-saving); not built
    here);
  - in an outside written argument: 6.7e-4 (shea256; conditional, not formalised). That repository gives
    another figure since 10 October 2026 (item 1 above).
- **Second result: a machine-checked saving of 7.474547e-4 for the kernel program and for the Walsh-Hadamard
  transform.** When I published it on 9 October 2026 it was the largest figure for this quantity that I could
  find anywhere. The largest outside figure was 7.009184e-4, certified by scripts: the "complex side" of the
  community's networks, in their pull request 193. It is not the largest any more. Larger figures for the
  five-stage layout of this repository (next item) are public. They are not theorems of this repository. Two
  of them are prices that its own Python tools compute:
  - 7.547361e-4 for the same circuit with Chafik Boukhalfa's hand-over pairs (the community's pull request
    233). The rebuild in `tools/rebuild/` produces that circuit and `gxdry.py` prices it. His own pull
    request #2 to this repository has that figure as a Lean theorem, which the agents built on a Linux
    machine (section "Whose circuit this is").
  - 7.635875e-4, the price that `gxdry.py` gives to the certificate of
    [pull request 304](https://github.com/CrocSwap/integer-mult-bounds/pull/304) by DreamingOfClouds. The pull
    request states 7635/10^7 for it.

  The fourth result of this repository, 8.7625e-4, is above both and is a Lean theorem.
- **A five-stage "bridged" layout.** As far as I could find, it was in no other work. The same helper circuit
  gives 7.009e-4 in the three-stage layout its authors use and 7.474547e-4 in this one. Since then others have
  taken the layout up: pull request 304, for one, applies the five-stage formula and credits this repository
  for the layout.
- **A Lean proof of the general frame lemma for this machine model.** As far as I could find, it was the
  first. It includes the cost of the lemma's "adapter" steps, which until now was argued on paper only.
  Section 5 says what the lemma is.
- **A certificate checker with soundness proofs.** A better circuit in its format is now a data file for it,
  not a new proof. The fourth result is the first result of this repository that uses it that way. Others
  used it first: by their own reports, DreamingOfClouds had the generators of this repository write Lean
  files for circuits of their own, and Lean's kernel accepted them (pull requests 315 and 327 of the
  community's repository; not built here).
- **A complete machine check of the community's circuit of pull request 193**, as rebuilt here from their
  published data (section 3 says what that rests on). As far as I could find, it was the first. Its own
  validation was local; its author's words are quoted in the section "Whose circuit this is".

Two things behind this list are not mine. The step from the kernel program to every length is OpenAI's
reduction, and the helper circuit of the second and third results is the community's. The next two sections
give the credit. Section 5 says in plain words what each of these items is.

The searches behind "as far as I could find" were read-only looks at public sources by the AI agents I
directed. They were not repeated in full for the fourth revision. On 10 October 2026 the agents read the
outside repositories named below again, and the newest pull requests of the community's repository.
RELATED-WORK.md says what the looks covered and when, and section 9 has the dates.

## Whose reduction this is

The step from the kernel program to the Fourier transform of every length is OpenAI's, in every part. The
kernel program is the core program, for the lengths 2^k; the Walsh-Hadamard theorem of the second result
rests on it too. The step from those lengths to every length is called the reduction here.

OpenAI's Lean proof uses the network through one theorem (`hills_program` in their `TensorProgram.lean`),
and that theorem is used in one place. The generalised engine of this repository proves the same sentence
with my exponent in place of theirs. So the third result needed no new mathematics after the kernel program.

**What is repeated from OpenAI.** OpenAI has 11 files between that theorem and their final statement. They
are repeated here with one constant changed: the exponent. The proofs are OpenAI's.

- **The files.** `SectorAlgorithm`, `SynchronizedAlgorithm`, `WorkingTransform`, `WorkingCompiler`,
  `WorkingPreparation`, `ArbitraryLength`, `TransformProgram`, `ConvolutionProgram`, `UniformBounds`,
  `Asymptotics` and `Main` (1,825 lines in OpenAI's repository). The copies are under `Work/Fourier/`. A script
  made them (`tools/fourier/mkchain.py`).
- **The constant.** In the copies:

  | in OpenAI's files | in the copies |
  |---|---|
  | the exponent `alpha` = 1 - 2/10^11 | `alphaZ` = 1 - 7474547/10^10 |
  | the work bound `hills`, defined from it | `hillsZ` |
  | the one use of `hills_program` | `hillsZ_program` |
  | the exponent of the final statement, 1 - 1/10^13 | 1 - 7474546/10^10 |

- **The names.** Every declaration of the copies has the letter `Z` in its name, so that Lean cannot take
  OpenAI's original for the copy, and the import lines point to the copies. Of the 167 declarations, 115 have
  the letter at the end of their own name. The other 52 have it on the name of the declaration they belong
  to, as in `WeftZ.transfer`.
- **Two further edits.** The first comment of each file names its own.
  - In the copy of `Main` the same constant stands once as a literal number, 2/10^11, and is replaced there
    too. The two final theorems of that copy also received a comment of their own.
  - In the copy of `SynchronizedAlgorithm` the name `canopy` is written out as `Grove.canopy` in 5 places,
    because in the larger environment of the copy the short name would mean another declaration.
- **The comments.** The script also put the letter on three names inside OpenAI's comments. The comments are
  otherwise OpenAI's and describe OpenAI's network.
- **A twelfth small file.** `Work/Fourier/Goal.lean` repeats the five definitions of OpenAI's `Goal.lean` that
  contain the exponent of the statement.
- **The copies are marked as modified.** Each begins with a comment that names the OpenAI file it is derived
  from and says what was changed, as the Apache-2.0 licence asks, and [NOTICE](NOTICE), section 7, lists
  them. OpenAI's own files under `OAI/` are unmodified, as before, and the proof still imports them.

**What is new Lean for the third result.** One file, `Work/Fourier/Seam.lean` (245 lines), written by the
agents. It states the kernel program of the generalised engine in the form that OpenAI's reduction asks for,
and proves the few facts about the work bound that the reduction uses.

**What the reduction costs.** A factor (log log n)^2 at the kernel's own exponent. That is why the Fourier
saving is 7.474546e-4, just below the kernel's 7.474547e-4 (section 2, point 4).

**The fourth result repeats the same reduction once more, for the new kernel program.** `Work/FourierE8/` is a
twin of `Work/Fourier/`:

- It has the same 11 files of OpenAI, copied by the same kind of script (`tools/fourier/mkchain_e8.py`, a
  copy of `mkchain.py` with the constants changed).
- It has the suffix `E8` where the third result has `Z`.
- It has `alphaE8` = 1 - 8762479/10^10, and the exponent 1 - 8762478/10^10 in the final statement.
- Its seam file `Work/FourierE8/Seam.lean` (251 lines) is the seam of the third result with names and numbers
  replaced. No proof step had to be changed.

Everything this section says about the copies of the third result holds for these copies too, and
[NOTICE](NOTICE), section 7, lists them. The two chains share no name, so both are in the repository side by
side.

**Others reached this statement before me, and the credit for that is theirs.**

- **danadran01** wrote the first Lean proof that improves OpenAI's Fourier statement, as far as I could
  find: [danadran01/exact-dft-power-saving](https://github.com/danadran01/exact-dft-power-saving). Its saving is
  3.2e-6, on OpenAI's own `DFTProgram` and `TimeBounds`, with whole-block recursion in Lean under the name
  "grouped recursion". Their development changes three of OpenAI's files, all three at the theorem named
  above, and by their record the rest of OpenAI's chain then goes through as it is. Their record reports a
  clean build, the standard axioms and `leanchecker`, and no comparator run. It was not built here.
- **shea256** wrote the first written transfer of a community network to this statement:
  [shea256/fourier-transform-below-nlogn](https://github.com/shea256/fourier-transform-below-nlogn). Its first
  commit is of 8 October 2026. eumemic writes of it: "That work has priority for the idea of the transfer".
  On 9 October 2026 its figure was 6.7e-4. It was described there as "a proposed conditional transfer,
  supported by a written argument and finite checks", and it was not formalised. On 10 October 2026 it gives
  7.547360e-4 instead. That is the figure of Chafik Boukhalfa's pull request #2 of this repository, a Lean
  proof. The agents built that pull request on a Linux machine on the same day ("Whose circuit this is").
  shea256 reports a rebuild of that proof in its own repository, which was not repeated here, and writes
  that the repository "does not claim a new exponent of our own".
- **eumemic** wrote a written transfer with batched recursion:
  [eumemic/exact-dft-bounds](https://github.com/eumemic/exact-dft-bounds). It gives every saving below
  4.856e-4, with an extra log log factor. In their words it is "a paper proof with an exact finite
  certificate, and it is not formally verified". Their README also states in writing the observation that
  the route here rests on: OpenAI's later sections use the kernel theorem "only through three things: the
  bound O(2^k (k+1)^θ), the fact that only rational constants and i are needed (their §5.3), and its
  word-size accounting".
- **The community** made the helper circuit behind the kernel program (next section). The accounting is
  their whole-residual batching (RELATED-WORK.md, section 1).

Section 9 has the commits at which these three repositories were read, and what "as far as I could find"
rests on.

**What is mine in the third result:** the kernel program at this exponent with its machine check (the second
result), the seam file, and the machine check of the whole against OpenAI's statement.

**How the figures compare.** The table has the three outside figures as they stood on 9 October 2026.

| | outside figure | my 7.474546e-4 is |
|---|---|---|
| danadran01, in Lean | 3.2e-6 | 233 times that |
| eumemic, written | every saving below 4.856e-4 | 1.54 times that |
| shea256, written | 6.7e-4 | 1.12 times that |

The two written figures rest on other networks, and they were moving within hours while this was written.
The agents read the three repositories again on 10 October 2026. Two had not changed. shea256's now gives
7.547360e-4, the figure of Chafik Boukhalfa's pull request #2 of this repository ("What is new in the fourth
revision", item 1).

## Whose circuit this is

This section is about the circuit of the second and third results. The circuit of the fourth result is another
one; its credits are in "What was new in the fourth revision", under "Whose ideas the new unit uses". The
circuit of the fifth result is the community's U8 circuit; its credits are in "What is new in the fifth
revision".

The helper circuit is the part of the network that decides the size of the saving. It is the additions by
which pairs of data arrays exchange their contents, with the helper arrays (extra arrays that hold the
partial sums) that this needs; section 6 describes it. In the second result, and so in the third, it is the
community's in every part. The community is the group of people who publish their work as pull requests to
the public repository [CrocSwap/integer-mult-bounds](https://github.com/CrocSwap/integer-mult-bounds). A
number with # below is a pull request there, unless the sentence says otherwise.

- **icekylinx** made the paired-cube construction
  ([#144](https://github.com/CrocSwap/integer-mult-bounds/pull/144)). Its notice credits an664 (#128) for the
  workspace-sharing principle and eumemic (#117) for the producer that its modules restrict.
- **eumemic** made the modules, carrier links, frames and physical layer at p = 11
  ([#168](https://github.com/CrocSwap/integer-mult-bounds/pull/168), "PR168 v4").
- **icekylinx** made the local circuit inside each cube that takes a + b and a - b from the same two arrays,
  and the frame flow that goes with it ([#184](https://github.com/CrocSwap/integer-mult-bounds/pull/184)).
- **ikeboy** put that construction on #168's newer modules
  ([#191](https://github.com/CrocSwap/integer-mult-bounds/pull/191)).
- **ikeboy** made the choice of hand-over pairs with which that step runs in place and erases nothing
  ([#193](https://github.com/CrocSwap/integer-mult-bounds/pull/193)). A hand-over pair says which helper
  array that is no longer needed is handed to which new user. In their accounting this last step is the whole
  jump from 6.6307e-4 to 7.009184e-4.

**How #193 states its own credits.** "icekylinx (PR184, with GPT-6 Astra and Codex assistance) for the
method, the tools and the bit supplier. eumemic (PR168 v4, with Claude assistance) for the query modules and
the physical layer. Package by Avi Eisenberg with Claude assistance."

**Techniques of other authors** reach the circuit through these pull requests: slot reuse at birth, frame
descent, the closure-compiled carrier links of DaysSky (#162), and the fusion of the single-target outputs
after chafreaky (#163). They are credited in RELATED-WORK.md, section 1.

**How the circuit came into this repository.** The agents rebuilt the circuit as one explicit list of
additions, from the files that these pull requests publish: the module files, carrier links and frames of
#168 and the pair list of #193 ([NOTICE](NOTICE) names each file). The programs that did this are in
`tools/rebuild/`.

- The agents made them with the community's scripts open as text.
- Six of the thirteen files are adaptations of scripts by icekylinx, eumemic, DaysSky and ikeboy, in parts
  statement by statement. [NOTICE](NOTICE), section 6, says which follows what.
- The other seven were written for this repository. One of them, `rewrite.py`, is the step that turns the
  word of #193 into one explicit list of gates. (The word is the network written as a sequence of stages.)
- No outside program was run. That repository is under Apache-2.0.

**Another choice of hand-over pairs, by Chafik Boukhalfa.** Another choice of the 2,310 hand-over pairs for
the same word is by chafreaky (Chafik Boukhalfa),
[#233](https://github.com/CrocSwap/integer-mult-bounds/pull/233), in his words "PR #200's maximum-weight late
compensated pairing" ([#200](https://github.com/CrocSwap/integer-mult-bounds/pull/200)). He proposed the swap
in issue #1 of this repository. With his pairs in place of #193's, the rebuild in `tools/rebuild/` gives a
certificate that my checker accepts. `gxdry.py` prices that certificate at 7.547361e-4; the published
certificate has 7.474547e-4. That is a price computed in Python: this repository has no Lean modules for that
certificate, and its theorems do not use it (`tools/gx/ORIGIN.md`, "Another file of hand-over pairs"). His
pull request #2 to this repository goes further. It holds a certificate of his own making for the same
circuit, with the same block histogram and the same price, the Lean modules that the unchanged generators
of this repository write from it, and a Fourier chain. He reports a Lean kernel check of the Walsh-Hadamard
statement at 1 - 7547361/10^10 and of the Fourier statement at 1 - 7547360/10^10, with `tools/Compare.lean`
and a comparator run on Linux. With the unchanged generators I reproduced his 41 generated files byte for
byte, and my checker and `gxdry.py` accept his certificate at 7547361. On 10 October 2026 the agents built
his pull request on a Linux machine. All 54 Lean modules that it adds compiled. His three theorems depend
only on the three standard axioms. `tools/Compare.lean` printed `RESULT: PASS` for his Walsh-Hadamard
statement and for his Fourier statements. After the build, the agents ran the official comparator on the
two configurations for which he reports runs, with a stand-in for its sandbox. It accepted both, each with
`Your solution is okay!` (74 and 73 minutes; 14.4 and 14.5 GB). He reports his own comparator runs, made
with the real sandbox. On my Mac the build could not finish, because two of his modules need about 13 GB of
memory each. So in his pull request 7.547361e-4 is a Lean theorem, and it is his result. Section 9 has the
details of the build and of those runs. The pull request is not merged: I would rather this repository make
that circuit with its own scripts, which `tools/rebuild/` now does, than carry generated files.

**What is mine:** the Lean proof of the general frame lemma for this RAM model, the generalised engine, the
extended certificate checker (section 5), the bridged five-stage layout in which their circuit is placed
(section 6), and the machine check of the whole.

**What the machine check here adds.** This is the first full machine check of that circuit that I know of, as
of 9 October 2026 (section 9 says what that rests on). The author of #193 states the scope of that pull
request's own validation in these words: "The exact lift and the contract checks establish the local maps and
the flow ledger. There is no globally renumbered scalar transcript of the new complex word and no full
Clifford/router replay." What is checked here is the circuit inside my statement. Their theorem is about
integer multiplication and has further parts that this repository does not touch (section 2, point 1).

**Why my figure for the same circuit is above theirs.** #193 states 7.009184e-4 for this circuit, and the
second result has 7.474547e-4. The layout differs: their word has three stages, mine has five. By the agents'
computation the same circuit in a three-stage word gives 7.0091e-4.

## 1. The results

### The eighth result: 1.7887e-3, on Shea's Lean tree (new on 2026-10-20)

The full cut on #2979's words with the table A-fc gives the largest saving of this repository: 1.7887e-3 for the
Fourier transform of every length and for convolution. Its proof is the package `package/leverAFc`, laid over
Ryan Shea's tree `lean/complex588` after the package `package/u8x0hq` (VERIFY.md, section 12). The definitions that
carry the exponent, in `Work/FourierL2979Fc/Goal.lean` and in the statement file
`Work/FourierL2979Fc/UniformFourierChallenge.lean`:

    def decimalExponentL2979Fc : ℝ := 1 - 17887045/(10^(10:ℕ))
    def decimalTimeL2979Fc (n : ℕ) : ℝ := (n:ℝ) * (Real.log (n:ℝ)) ^ decimalExponentL2979Fc
    def TimeBoundsL2979Fc (W : ℕ → ℕ) : Prop :=
      let W' := fun n : ℕ => (W n : ℝ)
      W' =O[atTop] paperTime ∧ W' =O[atTop] decimalTimeL2979Fc ∧ W' =o[atTop] nlogn
    def DFTGoalL2979Fc : Prop := ∃ order solve W, DFTProgram order solve W ∧ TimeBoundsL2979Fc W
    def ConvGoalL2979Fc : Prop := ∃ order solve W, ConvProgram order solve W ∧ TimeBoundsL2979Fc W

    theorem OAI.PowerSaving.transform_mainL2979Fc : DFTGoalL2979Fc
    theorem OAI.PowerSaving.convolution_mainL2979Fc : ConvGoalL2979Fc

"What is new in the sixth revision", at the top, has the numbers, and the opening has the checks and the credits.

### The lever on #3945: 1.7807e-3, on Shea's Lean tree (new on 2026-10-19)

The full cut on #3945 with every head kept gave the largest saving of the fifth revision: 1.7807e-3 for the
Fourier transform of every length and for convolution, with the full cut. Its proof is the package `package/leverF`, laid over Ryan Shea's
tree `lean/complex588` after the package `package/u8x0hq` (VERIFY.md, section 12). The definitions that carry the
exponent, in `Work/FourierL3945F/Goal.lean` and in the statement file
`Work/FourierL3945F/UniformFourierChallenge.lean`:

    def decimalExponentL3945F : ℝ := 1 - 17806935/(10^(10:ℕ))
    def decimalTimeL3945F (n : ℕ) : ℝ := (n:ℝ) * (Real.log (n:ℝ)) ^ decimalExponentL3945F
    def TimeBoundsL3945F (W : ℕ → ℕ) : Prop :=
      let W' := fun n : ℕ => (W n : ℝ)
      W' =O[atTop] paperTime ∧ W' =O[atTop] decimalTimeL3945F ∧ W' =o[atTop] nlogn
    def DFTGoalL3945F : Prop := ∃ order solve W, DFTProgram order solve W ∧ TimeBoundsL3945F W
    def ConvGoalL3945F : Prop := ∃ order solve W, ConvProgram order solve W ∧ TimeBoundsL3945F W

    theorem OAI.PowerSaving.transform_mainL3945F : DFTGoalL3945F
    theorem OAI.PowerSaving.convolution_mainL3945F : ConvGoalL3945F

`DFTProgram`, `ConvProgram`, `paperTime` and `nlogn` are OpenAI's, unchanged. The other lever packages of the
table "All the results" have the same form with their own suffix (the half cut on the same table: `package/leverB`,
suffix `L3945H`, at 1 - 17581520/10^10). "What is new in the fifth revision" has the checks and the credits.

### The sixth result: 1.6257e-3, on Shea's Lean tree (14 October 2026)

The sixth result is `package/u8x0hq`, laid over Ryan Shea's tree (VERIFY.md, section 12), with the theorems
`OAI.PowerSaving.transform_mainU8X0HQ` and `OAI.PowerSaving.convolution_mainU8X0HQ` at the exponent
1 - 16256629/10^10, in the same form as above with the suffix `U8X0HQ`. "What is new in the fifth revision" has
its checks.

### The fifth result: 1.6198e-3, on Shea's Lean tree (new on 2026-10-13)

The fifth result was the largest saving of this repository until the sixth: 1.6198e-3 for the Fourier transform of every length
and for convolution. Its proof is not in this repository's own Lean project. It is the package
`community/u8xmhh`, laid over Ryan Shea's tree `lean/complex588` after the package `community/u8xm`
(VERIFY.md, section 11). The definitions that carry the exponent, in `Work/FourierU8XMHH/Goal.lean` and in the
statement file `Work/FourierU8XMHH/UniformFourierChallenge.lean`:

    def decimalExponentU8XMHH : ℝ := 1 - 16197705/(10^(10:ℕ))
    def decimalTimeU8XMHH (n : ℕ) : ℝ := (n:ℝ) * (Real.log (n:ℝ)) ^ decimalExponentU8XMHH
    def TimeBoundsU8XMHH (W : ℕ → ℕ) : Prop :=
      let W' := fun n : ℕ => (W n : ℝ)
      W' =O[atTop] paperTime ∧ W' =O[atTop] decimalTimeU8XMHH ∧ W' =o[atTop] nlogn
    def DFTGoalU8XMHH : Prop := ∃ order solve W, DFTProgram order solve W ∧ TimeBoundsU8XMHH W
    def ConvGoalU8XMHH : Prop := ∃ order solve W, ConvProgram order solve W ∧ TimeBoundsU8XMHH W

    theorem OAI.PowerSaving.transform_mainU8XMHH : DFTGoalU8XMHH
    theorem OAI.PowerSaving.convolution_mainU8XMHH : ConvGoalU8XMHH

`DFTProgram`, `ConvProgram`, `paperTime` and `nlogn` are OpenAI's, unchanged. The packages of the same day with
smaller savings (`community/c588mhh` at 1 - 15699305/10^10, `community/u8xm`, `community/c588m`,
`community/u8x`, `community/c588-repack`) have the same form with their own suffix. The largest Walsh-Hadamard
saving of the day is on the kernel program of `community/c588mhh`:

    theorem OAI.PowerSaving.WHT.wht_main_block_C588MHHx :
      ∃ solve W, WHTProgram solve W ∧ WHTTimeBoundsAt (1 - 15699306/(10:ℝ)^10) W

with `WHTProgram` and `WHTTimeBoundsAt` as in the second result below; it is in `community/wht-c588mhh`.
"The results of 13 October 2026", above, has the checks and the credits.

### The fourth result: both statements at 8.7625e-4 (new on 2026-10-10)

The fourth result is the largest saving proved in this repository's own Lean project: 8.7625e-4, for the Walsh-Hadamard transform, for
the Fourier transform of every length and for convolution. The saving is 1 - z, where the program needs
O(n (log n)^z) operations. The second and third results have 7.4745e-4.

**What is new is the unit.** The unit is the finite network that the saving is computed from: a fixed set of
arrays (working copies of the data, each of length n) and a fixed list of steps on them. The figure is larger
for one reason only: the unit is another one, and it is much smaller.

| | fourth result | second and third results |
|---|---|---|
| bank pairs (a bank pair is two data arrays that exchange their contents) | 120 | 1,320 |
| helper arrays (extra arrays that hold partial sums) | 783 | 9,412 |
| arrays in all | 1,263 | 14,692 |
| dimension m of the label space (each array has to be transformed along m directions) | 45 | 110 |
| moves in one pass of the circuit (a move transforms one array along one direction, and is what costs) | 9,039 | 262,944 |
| blocks in one pass of the circuit (a block is a run of moves of one array) | 3,974 | 70,169 |
| saving proved for the Walsh-Hadamard transform | 8.7625e-4 | 7.4745e-4 |

Section 6 describes the new unit and has the full table.
[notes/how-the-e8-unit-was-found.md](notes/how-the-e8-unit-was-found.md) says how it was found.

**What the theorems say.** They are the statements of the third and second results with another number, and
with other names so that all of them can stand in one repository.

    theorem OAI.PowerSaving.transform_mainE8 : DFTGoalE8
    theorem OAI.PowerSaving.convolution_mainE8 : ConvGoalE8
    def decimalExponentE8 : ℝ := 1 - 8762478/(10^(10:ℕ))

    theorem OAI.PowerSaving.WHT.wht_main_block_B2Ge8x :
        ∃ solve W, WHTProgram solve W ∧ WHTTimeBoundsAt (1 - 8762479/(10:ℝ)^10) W

- **Fourier and convolution.** One fixed program computes the discrete Fourier transform of every length
  n >= 1 in O(n (log n)^z) operations with z = 1 - 8762478/10^10, and in o(n log n). A second fixed program
  does the same for convolution.
  - Proof: module `Work.FourierE8.Main`.
  - Statement: module `Work.FourierE8.UniformFourierChallenge` (348 lines, imports only Mathlib). It is
    OpenAI's challenge file with the same three changes as in the third result: a comment, the exponent, and
    a suffix on the same seven names (here `E8`).
  - `DFTGoalE8`, `ConvGoalE8` and `TimeBoundsE8` are the definitions shown below for the third result, with
    that suffix and that exponent.
- **Walsh-Hadamard.** One fixed program computes the Walsh-Hadamard transform of every length n = 2^k in
  O(n (log2 n + 1)^z) operations with z = 1 - 8762479/10^10 = 0.9991237521, and in o(n log n).
  - Proof: module `Work.GCert.Data.SolutionB2Ge8x`.
  - Statement: module `Work.GCert.Data.ChallengeB2Ge8x` (`Work/GCert/Data/ChallengeB2Ge8x.lean`, 331 lines,
    imports only Mathlib). It is the statement of the second result with the exponent replaced, a new theorem
    name, and comments that quote the new names.

**What the theorems rest on.**

- The Walsh-Hadamard theorem needed no new Lean proof. The circuit of the new unit is written down as a
  certificate: a data file that lists every addition of the circuit. The checker of the second result, a
  Lean function with proofs that what it accepts is sound, reads that file
  (`tools/certificate/gcert1-e8-r783.json.gz`). Section 3 says how the result was checked, and section 5 says
  which parts of the certificate format this circuit is the first to use.
- The Fourier and convolution theorems are the same kernel program carried through OpenAI's reduction. The
  kernel program is the core program, for the lengths 2^k. The reduction is the step from those lengths to
  every length ("Whose reduction this is", above).

**Where the number comes from.** `tools/gx/gxdry.py` prints the figure 8762479 for this unit, in whole-block
accounting at ten decimals, and the generated rate lemma proves that figure. (Whole-block accounting counts a
run of r moves of one array as one recursive call; section 2, point 3.) The Python tool evaluates the rate
inequality as true up to 8762482 and as false at 8762483. The three figures in between are neither proved nor
claimed. The Fourier statement is one unit lower, for the reason given under the third result (section 2,
point 4).

**What the theorem does not say.** The statement gives the exponent. It says nothing about the unit. That
the certificate behind it has the 120 labels described in section 6 is tested by the search programs
(`tools/e8/`). It is not tested in Lean, and not by the two checkers of `tools/gx/`: those two accept, with
the same figure, a copy of the circuit carried to another family of 120 labels. The theorem does not depend
on which family it is.

Everything said below about the cost model, about whose statements these are and about what `DFTProgram`,
`ConvProgram` and `WHTProgram` mean holds for the fourth result word for word.

### The Fourier transform of every length (third result)

One fixed program, in the RAM model of OpenAI's proof, computes the discrete Fourier transform of every length
n >= 1 in O(n (log n)^z) operations with z = 1 - 7474546/10^10, and in o(n log n). The saving 1 - z is
7.474546e-4. The logarithm is the natural one, as in OpenAI's statement. A second fixed program does the same
for convolution.

    theorem OAI.PowerSaving.transform_mainZ : DFTGoalZ
    theorem OAI.PowerSaving.convolution_mainZ : ConvGoalZ

    def DFTGoalZ : Prop := ∃ order solve W, DFTProgram order solve W ∧ TimeBoundsZ W
    def ConvGoalZ : Prop := ∃ order solve W, ConvProgram order solve W ∧ TimeBoundsZ W
    def TimeBoundsZ (W : ℕ → ℕ) : Prop :=
      let W' := fun n : ℕ => (W n : ℝ)
      W' =O[atTop] paperTime ∧ W' =O[atTop] decimalTimeZ ∧ W' =o[atTop] nlogn
    def decimalTimeZ (n : ℕ) : ℝ := (n:ℝ) * (Real.log (n:ℝ)) ^ decimalExponentZ
    def decimalExponentZ : ℝ := 1 - 7474546/(10^(10:ℕ))

Proof: module `Work.Fourier.Main`. The statement with every definition it depends on: module
`Work.Fourier.UniformFourierChallenge` (`Work/Fourier/UniformFourierChallenge.lean`, 348 lines, imports only
Mathlib). It is
OpenAI's challenge file `lean/ComparatorChallenges/UniformFourier.lean` (336 lines) with three changes, which
its first comment lists: that comment; the exponent, `1 - 1/(10^(13:ℕ))` in OpenAI's file; and the letter `Z`
at the end of seven names (`decimalExponent`, `decimalTime`, `TimeBounds`, `DFTGoal`, `ConvGoal`,
`transform_main`, `convolution_main`). The letter is forced: the proof imports OpenAI's unmodified files, in
which these names carry OpenAI's values. Every other byte of the file is OpenAI's.

- `DFTProgram order solve W`: `order` and `solve` are two programs of the RAM, `solve` in the mode without
  general products. For every n >= 1, `order` computes from n a number d with 0 < d < 1024 n^3. For every
  complex vector x of length n, the run of `solve` on (n, exp(2 pi i / d), x) is valid (no division by zero)
  and returns exactly `dft n x` (entry j is the sum over k of exp(2 pi i j k / n) x_k). The work of `order`,
  one unit for the root, and the work of `solve` are together at most W(n). Every integer the two programs
  produce, and W(n) + 10(n+2), is at most (n+2)^c for one constant c. The root of unity is an input because
  the RAM has no instruction that produces it.
- `ConvProgram order solve W`: the same for the convolution of two vectors of length n (`conv n x y`, of
  length 2n - 1), with `solve` in the mode with general products and a root of order below 1024 (2n-1)^3.
- `TimeBoundsZ W`, three bounds as n grows: W(n) = O(n (log n)^(1 - d)) with d = 7.474546e-4
  (`decimalTimeZ`; the one line whose value differs from OpenAI's, where d = 10^-13); W(n) =
  O(n (log n)^theta (log log n)^(4 - theta)), the bound of OpenAI's manuscript with its own theta, about
  1 - 2.1e-13 (`paperTime`, unchanged); and W(n) = o(n log n).

The cost model and its caveats are the same for all four results. They are stated below, under "The cost
model is OpenAI's, with its caveats".

**Whose statement this is.** OpenAI's, with one number replaced. No definition in it was written for this
repository.

**Where the exponent comes from.** The kernel program runs at saving 7.474547e-4; the theorem of the second
result rests on it. OpenAI's reduction turns it into a Fourier program whose work is
O(n (log n)^z (log log n)^2) at that z. A pure power of log n therefore needs an exponent slightly above z,
and the statement has the saving 7.474546e-4 (section 2, point 4).

### The Walsh-Hadamard transform (second and first results)

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

**Whose statement the Walsh-Hadamard one is.** Part OpenAI's, part my project's.

- OpenAI's repository states and proves the discrete Fourier transform of every length (`DFTProgram`,
  `TimeBounds`, stated saving 10^-13). Its family 130 has no Walsh-Hadamard statement.
- The RAM cost model in my statement is OpenAI's, byte for byte (from their
  `lean/ComparatorChallenges/UniformFourier.lean`).
- The about 60 lines after it, which define `wht`, `WHTProgram` and `WHTTimeBoundsAt`, are my project's,
  written on the pattern of OpenAI's `DFTProgram`. They should be read. Section 9 has the line ranges.
- The agents first proved the statement with the exponent 1 - 2/10^11 of OpenAI's tensor engine (`alpha` in
  their `TensorSaving.lean`) as a corollary of that engine (`WHTCheck/Solution.lean`). Here only the exponent
  is replaced.

**The first result.** `wht_main_block_B2Ke16x`, saving 5.399225e-4 (z = 1 - 5399225/10^10), with my project's
own helper circuit and the first engine.

- Proof: `Work.CarrierCheck.SolutionB2Ke16x`. Statement: `Work/CarrierCheck/ChallengeB2Ke16x.lean` (331 lines).
- It was published on 9 October 2026, and its files are in this repository as they were published.
- So is its companion with OpenAI's own recursion, `wht_main_rank_B2Ke16x`: the same statement with
  1 - 3155781/10^10 (`Work.CarrierCheck.SolutionB2Ke16xR`).

The two figures of the first result are two accountings of the same network.

- Per-rank: every single move is one recursive call on 1/m of the coordinates, as in OpenAI's proof.
- Whole-block: r consecutive moves of one array are one call on r/m of the coordinates. That idea is the
  community's "whole-residual batching" (RELATED-WORK.md). Its Lean proof for this RAM model is in this
  repository.

The second result exists in whole-block accounting only (section 2, point 3).

## 2. What is not claimed

**The lever has a list of its own** ("What is new in the fifth revision", "What is not claimed for the lever"),
and so does **the fifth result** ("The results of 13 October in detail", item 7). In short: the circuit
is the community's U8 circuit, and the ideas of the three steps are the community's, Zhihao Chen's and
DaysSky's; the proof lives in Ryan Shea's tree; the community's own
figures are higher and rest on steps without a Lean proof; the leaderboard's own verifier was not run on it
(the agents re-created its three phases offline); and points 4 and 5 below (the step to every length, and growth) hold for it too.

0. **The fourth result has the limits of the second and third. Its figure is larger for one reason only: the
   unit.** Every point below holds for it with 8.7625e-4 in place of 7.4745e-4, except where a point says
   otherwise. In particular:
   - It is in whole-block accounting only (point 3).
   - Its Fourier saving is strictly below its kernel's: 8762478 against 8762479 (point 4).
   - It is a statement about growth (point 5). The sizes given there are those of the earlier units. They
     were not recomputed for the new unit, the one with m = 45.
   - It says nothing about integer multiplication. No bit side is offered here (point 1).
1. **The second result is not a better circuit, and it does not verify the community's theorem.** The second
   result places their circuit in my layout and checks it for my statement. Nothing here is about integer
   multiplication. Their headline figure is kappa, their saving for integer multiplication. Kappa also needs
   a second ingredient (the "bit side") and analytic and routing interfaces that their authors call
   assumptions. None of that is checked or used here. The third result uses the same circuit and has no
   circuit of its own. The fourth result has a circuit of its own, for my statement only.
2. **The outside figures are a different quantity, and they keep moving.** The titles of the pull requests to
   CrocSwap/integer-mult-bounds give kappa, the saving for integer multiplication. Kappa is roughly the
   smaller of two numbers: the "complex side" saving, which is the quantity of my theorem, and the bit side.
   For #193 the complex side is 7.009184e-4 and kappa is 6.647872e-4.
   - **On 9 October 2026** the pull requests were read twice. Both times the largest complex-side saving
     stated as a result was still 7.009184e-4, the circuit of #193, which several later pull requests keep.
     The largest kappa in a title was 6.831905e-4 the first time (Dugongue) and 6.839217e-4 the second
     (eumemic). One pull request prices further steps, up to kappa 1.226488e-3, and its title calls them
     "target, not built".
   - **On 10 October 2026** the 100 newest pull requests were read, and the list was read once more later
     that evening. The largest complex-side saving with a reported Lean check was 7795/10^7
     (DreamingOfClouds, [#327](https://github.com/CrocSwap/integer-mult-bounds/pull/327)), for a word at
     h = 20 in the five-stage layout of this repository. Its author reports a Lean kernel check of it with
     this repository's generators. It was not built here. A closed pull request there, by maxime-fleury,
     prices 8.3465e-4 for the complex side and itself calls that "a conditional priced target, not a
     witness". The largest kappa in a title was 7.77596e-4 (eumemic).
   - Section 9 has the full record of these readings: the hours, the pull-request numbers and what was
     read.
   - What kind of figure each of these is. Most are prices: their authors' Python scripts compute them, they
     are not peer reviewed, and they have no Lean proof. That holds for every kappa figure above, for the
     7.009184e-4 of #193 and for the 8.3465e-4 of the closed pull request. Their authors call them
     conditional. The one exception in this list is the 7795/10^7 of #327: its author reports a Lean kernel
     check of it, and that check was not built here.
   - [#192](https://github.com/CrocSwap/integer-mult-bounds/pull/192) (DaysSky) states a ceiling of 7.010e-4
     "for every frame layout of #168's word". #193 says of it: "The PR192 frame ceiling covers only PR168's
     fixed word and pairs, so it does not apply to this word." My figure is above it for one more reason, the
     five-stage layout.
3. **Whole-block accounting only.** There are two ways to count the cost of a network. Per-rank accounting
   counts every single move as one recursive call, as in OpenAI's proof. Whole-block accounting counts a run
   of r moves of one array as one call. The generalised engine has a whole-block theorem and no per-rank
   theorem. So the second result has no companion with OpenAI's own recursion. In per-rank accounting the
   result of this repository is still 3.155781e-4, from the first network. The third result rests on the same
   whole-block kernel program, so it has no such companion either: in OpenAI's own recursion the Fourier
   statement stands where OpenAI proved it.
4. **The Fourier saving is strictly below the kernel's.** The final statement of the third result carries
   no extra factor: it is O(n (log n)^(1 - d)) with d = 7.474546e-4. Behind it, OpenAI's reduction costs a
   factor (log log n)^2 on top of the kernel program. At the kernel's own exponent, z = 1 - 7.474547e-4, the
   chain therefore gives O(n (log n)^z (log log n)^2). A pure power of log n needs an exponent a little
   above z. The same argument would give every saving below 7.474547e-4. The theorem here states
   7.474546e-4, and it does not state 7.474547e-4 itself.
   - Outside figures for the same statement, as the agents found them on 10 October 2026:
     - In Lean: 7.547360e-4, in pull request #2 of this repository, by Chafik Boukhalfa. It is the circuit
       of his pull request 233 in the community's repository, in the layout of this repository.
       [shea256/fourier-transform-below-nlogn](https://github.com/shea256/fourier-transform-below-nlogn) now
       gives the same figure, with its own rebuild of that proof.
     - In Lean: 3.2e-6
       ([danadran01/exact-dft-power-saving](https://github.com/danadran01/exact-dft-power-saving)).
     - In a written argument that is not formalised: every saving below 4.856e-4 with an extra log log
       factor ([eumemic/exact-dft-bounds](https://github.com/eumemic/exact-dft-bounds)).
     - The agents built the proof of pull request #2 on a Linux machine on 10 October 2026, and every
       module of it compiled ("Whose circuit this is"). They did not build danadran01's proof or shea256's
       rebuild. Each is reported by its author. The last two figures rest on other networks.
     - On 9 October 2026 shea256's repository had 6.7e-4 in place of today's figure. That was a conditional
       written argument, not formalised.
   - The Walsh-Hadamard results are for the lengths 2^k only. The first of them, with my project's own
     circuit, has not been carried over to the Fourier statement.
   - The statement is OpenAI's `DFTProgram` and `ConvProgram`: transforms over the complex numbers, with a
     root of unity supplied to the program. Nothing is claimed for other transforms or other models.
5. **It is a statement about growth, not a usable algorithm.** The items below give the sizes involved. A
   number marked (c) was computed from a formula, not in Lean. The figures are those of the kernel program of
   the second and third results. The first result's are in square brackets.
   - The proof works with a table of 2^a arrays indexed by all orthogonal matrices of a 110-dimensional space
     over F_2 [80-dimensional]. The table is never written down. a = 6049 (c) [3213], and 2^a has
     1821 digits (c) [968].
   - Below log2 n = m (a + 1) = 665,500 (c) [257,120] the program is the ordinary n log n algorithm.
     Nothing changes for any input whose length has fewer than 200,336 digits (c) [77,401].
   - The constant is huge. As one of the agents read in the Lean proofs, the program fills all 2^a arrays with
     the input, which puts a factor 2^a into the constant. With it the bound is below n log2 n only when
     log2 log2 n exceeds a/(1 - z), about 8e6 (c) [6e6].
   - A second explicit factor is 1 + 1/slack of the rate inequality. The slack is 6e-10 here
     [2e-10 to 4e-10], by the auditing agent's arithmetic; the rational bounds of the Lean proof leave 2e-12.
     Machine constants were never computed, here or by OpenAI.
   - With every constant set to 1, a gain of 1% over the ordinary algorithm needs log2 n near 2^39
     (c, heuristic) [2^45].
   - For the third result these are figures of the kernel program, not of the Fourier program. The
     Fourier program calls the kernel program on arrays whose lengths it chooses itself: by the agents'
     reading of OpenAI's files, about (log log n)^2 sweeps on arrays of total length at most 8n. The
     lengths n at which the saving starts, and the constants of the reduction, were not computed, here
     or by OpenAI.

   In Lean: `tableExp n s = Nat.clog 2 n + s`, `threshold u a = u * (a + 1)`, m = 110, s = 40, 14,692
   arrays per unit [m = 80, s = 40, 9362], one unit per orthogonal matrix. Not in Lean: the order of that group,
   which gives a and so every (c).
6. **It is not a statement about bit complexity or about numerical computation.** One unit of work is one exact
   operation on complex numbers of any size, and the coefficients are unrestricted (section 1). OpenAI writes of
   its circuit result that it makes "no all-length, bounded-coefficient, conditioning, or bit-complexity claim".
   The last three limits hold here as well. The classical lower bound of order n log n for linear computations
   with bounded constants (Morgenstern, 1973) is not contradicted: this model does not bound the constants.
   All of this holds for the third result too. Its root of unity is an exact input.

What is real is the exponent: a theorem about growth for all n in this model, machine-checked.

## 3. How it was checked, and the limits

The AI agents ran every check. The checks of the first three results were made on one machine, a Mac, on
9 October 2026. One more was added on 10 October 2026: a build of the published third revision from nothing
on a rented Linux machine (it is reported under "Fourth result", below). The checks of the fourth result were
made on that Mac and on that Linux machine, on 10 October 2026. The limits are at the end of this section.
Please weigh them.

Three checking tools are named throughout:

- **Lean's kernel** is the part of Lean that accepts or refuses a proof.
- **The comparator** ([leanprover/comparator](https://github.com/leanprover/comparator)) is the tool that
  OpenAI's repository uses for its own results. It checks that the proof proves exactly the given statement,
  replays the proof in Lean's kernel, and admits three standard axioms only ("Second result", below, has the
  details). It is meant to run inside a sandbox, a confinement of the build; the limits say why no run so far
  had the real sandbox.
- **`tools/Compare.lean`** is a second comparison script.

The clock times and durations of the runs of the first three results are in section 9.

**Fourth result.** [VERIFY.md](VERIFY.md), section 10, has every command, the expected lines and the measured
times and memory.

- **Build.** The new unit has 24 Lean modules: 14 with the certificate data and its kernel evaluations, under
  `Work/GCert/Data/Gen/`, and 10 above them. The two unchanged generators wrote them from the certificate.
  They were compiled from source, on the compiled modules of the earlier results. Then came the 15 modules of
  `Work/FourierE8`. On the Mac: all 39 builds ended with exit code 0, one module per `lake build`: 183 s for the 24
  modules and 168 s for the 15, with 6.15 GB as the largest resident size of one Lean process. The same 39
  modules were built once more on the Mac from the Lean files of this revision alone, in a copy of the
  complete build of the published tree: 146 s for the 24 and 47 s for the one Fourier command of VERIFY.md
  10.1, exit code 0 throughout. The three theorems depend on the axioms `propext`,
  `Classical.choice`, `Quot.sound` and on no other: each of the three lines that Lean printed on the Mac
  ends with `depends on axioms: [propext, Classical.choice, Quot.sound]`.
- **`tools/Compare.lean`**, on the Mac: `RESULT: PASS` for `wht_main_block_B2Ge8x` (107 s; largest resident
  size 6.61 GB) and `RESULT: PASS` for the pair `transform_mainE8` / `convolution_mainE8` (103 s; 7.32 GB).
  The process's own memory (footprint, sampled every four seconds) peaked at about 9.8 GB in both runs.
- **Comparator**, with the configurations `comparator/B2Ge8x.json` and `comparator/UniformFourierE8.json`, on
  the Mac, with the stand-in for its sandbox: `Lean default kernel accepts the solution` and
  `Your solution is okay!` for both, exit code 0: 85 s and 6.05 GB for the first, 89 s and 6.09 GB for the
  second. Every module was replayed from its compiled file; none was compiled inside the run.
- **Linux.** The same 39 modules were compiled from source on a rented Linux machine (x86-64, Ubuntu
  24.04, 51 GiB of memory), in a copy of that machine's build of the published third revision (next point).
  `tools/Compare.lean` and the comparator were run there too, the comparator with a stand-in for its
  sandbox: all 39 builds ended with exit code 0 (1185 s in all; largest process 6.69 GB); the three axiom
  lines are the same as on the Mac; `tools/Compare.lean` printed `RESULT: PASS` for both commands (99 s and
  150 s; 15.7 and 16.0 GB); the comparator printed `Your solution is okay!` for both configurations (406 s
  and 432 s). The stand-in there is built on Linux namespaces, and the comparator replayed modules that had
  been compiled before the run. The machine was shared with other jobs, so its times are upper bounds.
  Memory on that machine, here and in the next point, is the kbytes figure that GNU `time` printed, divided
  by a million. A kbyte there is 1024 bytes, so in units of 10^9 bytes each figure is 2.4 % higher.
- **The third revision, built from nothing on Linux.** On 10 October 2026 the whole published third revision
  was also built from nothing on that machine. It was cloned from GitHub, and all 342 Lean modules of the
  repository were compiled from source with the commands of VERIFY.md (94 minutes; largest process 10.5 GB).
  Mathlib was not compiled from source: its compiled files came from Mathlib's cache server.
  `tools/Compare.lean` passed there for all five theorems of the first three results. That is the build
  from nothing that the limits below, written on 9 October 2026, say had not been made. The agents then ran
  the official comparator on that machine for the first three results, each time with a stand-in for its
  sandbox (the same one, built on namespaces), in a tree whose modules had been compiled before. It accepted
  all three, each with `Your solution is okay!` and exit code 0:
  - first result, `comparator/B2Ke16x.json` and `comparator/B2Ke16xR.json`: 37 and 25 minutes, 6.3 and
    6.4 GB;
  - second result, `comparator/B2Gp193x.json`: 124 minutes, 11.3 GB;
  - third result, `comparator/UniformFourier.json`: 108 minutes, 11.4 GB.

  The runs shared the machine, so the times are upper bounds. Two negative controls, one on the statement of
  the second result and one on the Fourier statement of the third, were refused with exit code 1.

  The agents also ran the comparator for the third result in a tree where none of the repository's modules
  had been compiled before. Mathlib's compiled files came from its cache. The comparator accepted the
  result, with `Your solution is okay!` and exit code 0. That was the third attempt in that tree, and each
  attempt continued the build of the one before. In the three attempts together, 314 modules of the
  repository were compiled inside the stand-in. The last attempt took 127 minutes and 11.4 GB. This verdict
  came after the fourth revision was published (section 9).
- **No real sandbox on that Linux machine.** The comparator's sandbox needs Landlock, a feature of the Linux
  kernel. The Linux kernel of that machine has no Landlock, so the comparator's real sandbox could not run
  there either (limits, below).
- **The generated files.** `python3 tools/gx/regen_check_e8.py` and `python3 tools/fourier/regen_check_e8.py`
  run the generators again and compare the result with the files of this repository. Both end with `RESULT:
  REPRODUCED`: 25 files for the unit (the 24 Lean modules and the comparator configuration) and 14 for
  `Work/FourierE8`, byte for byte.
- **The certificate, in Python.** Both published checkers accept it:
  - `refcheck.py` prints `ACCEPTED by gx.check1: h=9 v=120 R=783 N=9039`;
  - `gxdry.py` prints `MIRROR ACCEPTS` and `whole-block 8762479`.

  A stand-alone replay written by another agent, which reads only the certificate, gives the same numbers.
  Both checkers refuse seven deliberately damaged copies: a wrong coefficient in an in-place pair, a removed
  compensating read of a re-used helper, two helpers folded into one, a total declared one dimension lower,
  two labels exchanged, a source added with 1 in place of 1/2, and three blocks booked as one cheaper block.
- **Negative controls in Lean were run on the two earlier designs of the same day and then on this
  certificate.**
  A negative control is a check that has to fail: a damaged input, or a statement with a larger saving, that
  the tools must refuse.
  - Two earlier units of the same kind were built in Lean on 10 October 2026 during the work, at 7.5590e-4
    (7558959) and at 8.3202e-4 (8320202; its Fourier twin at 8320201).
  - Both were Lean-checked during the work: on the Mac (kernel build, `tools/Compare.lean`, comparator) and
    on the Linux machine (kernel build, `tools/Compare.lean`). Their Lean files are not included, since the
    certificate of this revision supersedes them.
  - On the first, on the Mac, an auditing agent gave Lean six damaged inputs and two raised figures. Lean
    refused five of the six and both raised figures, and the comparator refused both challenges with a
    larger saving. The sixth is the record of the next point. `tools/Compare.lean` could not be run to its
    end in that audit, because of memory.
  - On the second, on Linux, another agent made 13 builds that had to fail and 8 that had to pass, each
    outcome written down before the run. All 21 came out as written, and `tools/Compare.lean` failed both
    challenges with a larger saving.

  On the certificate of this revision an agent ran negative controls in Lean on the Mac on 10 October 2026,
  in a separate copy of the build tree. A second agent then audited that record. It compared the damaged
  files with the true ones and read the logs. It ran no Lean itself. This is what the audit supports:

  - **Twelve wrong Lean modules were refused**, and none of them left a compiled file. The unchanged text
    under other names built to the same theorem, on the three standard axioms.
  - **Six of the twelve had one wrong coefficient each.** Three were in half pairs (two additions in a row
    between two helpers, one of them with the factor 1/2). One was in an addition of a source into a helper
    with 1/2. One was in an addition of a helper into a helper at a frame of dimension 8 after the scatter
    (the step that hands the totals to the targets). One was in the scatter itself. Each differs from the
    accepted text in one token, and Lean's check of the coefficients refused each.
  - **Two had a wrong label.** In one, a helper is sent to a line (a frame of dimension 1) that does not lie
    in its next frame. In the other, a helper that climbs from a line to the full frame in one block has a
    wrong rank. Lean's check of the labels refused both, and it built a legal twin of the first.
  - **Three claimed a larger saving.** The rate lemma did not build at 8.762483e-4 or at 8.8e-4. A solution
    file that states 8.762483e-4 and cites the real theorem failed with a type mismatch.
  - **The twelfth** tried to profit from the gap of the next point, and Lean refused it.
  - `tools/Compare.lean` printed `RESULT: FAIL` for a challenge at 8.762483e-4 against the real solution.
  - The published Python checkers `refcheck.py` and `gxdry.py` refused every damaged certificate.

  Limits of these controls:

  - The published generator writes no Lean text for a damaged certificate. So the damaged Lean text came
    from an unpublished driver, or from an edit of generated text.
  - They touch 3 of the 21 half pairs, 1 of the 4 additions of a source into a helper with 1/2, 1 of the 5
    additions of a helper into a helper at dimension 8 after the scatter, and 1 of the 2 climbs from a line
    to the full frame. None is aimed at the table depth 10, the digit width 7 or the largest digit 36,
    three sizes that the check of the coefficients works with for this certificate.
  - A refusal by Lean says only that a check is false. It does not say which test inside the check failed.
  - There is no control on the Fourier modules and none on Linux, and the official comparator was not run
    on a wrong pair.
  - That each expected outcome was written down before its run rests on the time stamps in the agent's own
    file.
  - The audit is a reading of files and logs. No Lean run and no run of the Python checkers was repeated.
- **One damaged record that Lean accepts.** In those controls one wrong input went through: a wrong frame
  number in the certificate's list of totals. Lean does not use that number: the frame of a total is taken
  from the replayed state. For this certificate that was seen in a run for the first entry of the list. For
  the other seven entries it is a reading of the Lean functions that receive the list, not a run. With that
  one number wrong and nothing else changed, Lean built the same theorem, and the Python checkers
  `refcheck.py` and `gxdry.py` refuse the file. A version that tries to profit from the wrong frame, with
  all 8 totals priced one dimension lower, is the twelfth module above: Lean refused it. The entries of the
  published certificate are correct. So this is a field of the JSON file that only the Python checkers
  check, and about which "Lean accepted it" says nothing.

**Third result.**

- **Comparator.** The comparator (described under the second result, below) was run with the configuration
  `comparator/UniformFourier.json`. Its challenge is `Work/Fourier/UniformFourierChallenge.lean`, OpenAI's
  challenge file with the exponent replaced (section 1). It checks both theorems, `transform_mainZ` and
  `convolution_mainZ`. Its verdict "Lean default kernel accepts the solution" / "Your solution is okay!" was
  obtained twice, in the project's working tree, on Lean files that are byte for byte those of this
  repository: once by the agent that built the proof, and once by an auditing agent with its own
  configuration.
- **Audit by a separate agent** (another Claude session with its own scripts; not a human).
  - Besides its own comparator run it ran the second comparison script (`tools/Compare.lean`): pass, on the
    three standard axioms.
  - Both tools reject the same proof against a statement whose saving is one unit larger.
  - It compared the eleven repeated files with OpenAI's originals token by token. Apart from the renames, the
    differences are the seam, one constant, and one name written out in full in five places, which keeps
    OpenAI's meaning.
  - It walked both proofs, that is, followed everything they depend on. They reach every kernel check of the
    certificate and no `sorryAx` (the axiom that a `sorry`, an unproved step, leaves behind).
  - It did not read the repeated proofs line by line and did not rebuild from source.
- **What the kernel replays.** The Fourier theorem rests on the certificate of the second result. The
  comparator replays the whole proof in the kernel, the segmented evaluations of that certificate included.
- **Build.** The modules of `Work/Fourier` were compiled from source in the project's working tree, on the
  compiled modules of the second result (under two minutes; VERIFY.md, section 3). What is said below about the
  build of the second result holds here too: on 9 October 2026 no build from nothing had been made. One was
  made on 10 October 2026, on Linux ("Fourth result", above). In a copy of the
  repository the 15 new modules were compiled from source on the compiled modules of the second revision. The
  comparator and `tools/Compare.lean` were then run in that copy, and both passed.
- **The repeated files.** `tools/fourier/mkchain.py` (242 lines of Python) makes the 11 copies, `Goal.lean`, the
  challenge and the comparator configuration from OpenAI's unmodified files. `python3
  tools/fourier/regen_check.py` runs it again and compares the result with the files of this repository: all 14
  were reproduced byte for byte (VERIFY.md, section 6). So every difference between a copy and its OpenAI
  original is one that the script makes, and the script can be read. A second comparison, by another agent's own
  script, of each copy with its original line by line made the table in ORIGIN.md. Beyond the renamed names, the
  import lines and the added comments it found the two edits named in "Whose reduction this is", the exponent
  line of the challenge, and no other change to the code.
- **Audits by separate agents** (other Claude sessions; not humans). I had an auditing agent compare the new
  statements with OpenAI's as kernel terms, not as text. In the seam file, `hillsZ_program` has the type of
  OpenAI's `hills_program` with `hills` replaced by `hillsZ`; `hillsZ` is OpenAI's `hills` with `alpha` replaced
  by `alphaZ`; `alphaZ` is OpenAI's `alpha` with 2/10^11 replaced by 7474547/10^10; and the six side facts have
  OpenAI's statements. A walk of the proof term of `hillsZ_program` reached 38,544 constants: all 106 names of
  the certificate of the second result and of its soundness chain that the audit of the second result had listed
  (every kernel check of the certificate among them), no other certificate, none of OpenAI's `hills_program`,
  `hills` or `alpha`, no `sorryAx`, and the three permitted axioms only. The agent compared statements and
  walked proof terms; it did not read the proofs.

**Second result.**

- **Comparator.** The [comparator](https://github.com/leanprover/comparator), which OpenAI's repository uses
  for its own results, was run with the configuration `comparator/B2Gp193x.json`. It checks that the solution module proves
  exactly the challenge statement with identical definitions, replays the proof in the Lean kernel, and admits
  only the axioms `propext`, `Classical.choice`, `Quot.sound`, so a `sorry` or a `native_decide` is rejected.
  Its verdict "Lean default kernel accepts the solution" / "Your solution is okay!" was obtained three times:
  - twice in the project's working tree, on Lean files that are byte for byte those of this repository: once
    by the agent that built the proof, and once by an auditing agent with its own copy of the configuration;
  - once in the build copy of the repository (next item).
- **Build in a copy of the repository.** A copy with the Lean files, `lakefile.lean`, `lean-toolchain`,
  `lake-manifest.json` and `comparator/` of this revision was built on the compiled files of the first
  publication's rebuild from scratch (Mathlib from its cache, the 194 published modules compiled from nothing).
  The 102 new modules that do not load the new certificate were compiled from source in the copy. For the 31
  modules that load it, the compiled files of the project's working build were cloned: Lake accepted the 25 data
  and check modules by their traces and compiled the 6 modules above them from source. `lake build --no-build`
  of all eight roots then reported every target up to date (9250 jobs), with the three permitted axioms
  for the theorem and "declaration uses `sorry`" for the three challenge modules only. So the kernel evaluations
  of the certificate in those 25 modules were made by the compiler in the working tree, not in the copy, and the
  second result has not yet had a build from nothing. The comparator was then run in the copy and passed, with
  the same two verdict lines: it replays the whole proof in the kernel, those evaluations
  included. `tools/Compare.lean` then said `RESULT: PASS` there as well.
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

**First result** (the record of the first publication, of 9 October 2026; the clock times of its runs are now
in section 9, and nothing else in it has changed).

- **Comparator**, with the configurations `comparator/B2Ke16x.json` and `comparator/B2Ke16xR.json`. Its verdict
  "Lean default kernel accepts the solution" / "Your solution is okay!" was obtained
  - by the agent that built the proof: for the whole-block theorem and for the per-rank companion;
  - by an auditing agent with its own configuration: for the whole-block theorem;
  - in a rebuild from scratch: for the whole-block theorem and for the per-rank companion.
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
- **No line-by-line human review yet.** The new proofs are checked by the Lean kernel but have not yet had a
  line-by-line human review. The auditing agents checked statements, not proof texts. So the claim rests on
  the kernel and on the statement being the right one. I would be glad to have both read.
- **New in the fourth revision: one day, and read by AI agents only.** The fourth result was found, built
  and checked on one day, 10 October 2026. The certificate, the search programs and the notes have been read
  by AI agents only. The theorem does not depend on the search programs or on the notes: the kernel checks
  the certificate as it is.
- **New in the fourth revision: no run with the comparator's real sandbox yet.** No run with the comparator's
  real sandbox exists yet, for any result of this repository.
  - The real sandbox is the program `landrun`. It needs Landlock, a feature of the Linux kernel.
  - On macOS the comparator ran with a stand-in for `landrun`.
  - The Linux machine has a kernel without Landlock, on which `landrun` enforces nothing. So a run there is a
    run with a stand-in or with no sandbox.
  - The sandbox protects the checker against a solution file that attacks the build. It is not part of the
    logical check.
  - A run on a Linux machine with Landlock is the re-run I would welcome most.
- **New in the fourth revision: no new Lean proof in the Walsh-Hadamard chain.** The fourth result rests on
  the generalised engine, the checker and the segmented kernel check of the second result, and every limit
  stated here for those holds for it. Its circuit uses parts of the certificate format that the certificate
  of the second result never used (section 5). That the checker's Lean functions cover those parts was first
  a reading of the agents. It is now a matter of the kernel having evaluated them on this certificate.
- **New in the fourth revision: the Fourier chain is a second set of copies.** The Fourier chain of the
  fourth result is a second set of copies of OpenAI's 11 files, made by a copy of the script of the third
  result with other constants. Its seam file is the seam of the third result with names and numbers replaced.
  What is said below about the copies of the third result holds for these too.
- **New in the fourth revision: the negative controls.** The negative controls in Lean on the certificate of
  this revision are one agent's run on the Mac. A second agent audited its record by reading the files and
  the logs, and ran no Lean itself. The controls try one to three instances of each new kind of step, not
  all of them. There is none on the Fourier modules, none on Linux and none with the official comparator.
  The other sets were run on two earlier units of the same day (the list above). And one field of the
  certificate is not read by Lean ("One damaged record that Lean accepts", above).
- **New in the fourth revision: the earlier units of that day.** The Lean files of the earlier units of
  10 October 2026 are not in this repository. Where the notes give their figures, the figures are prices
  computed by the Python tools. Two of them (7558959 and 8320202) were also Lean-checked during the work, as
  said above.
- **New in the third revision: no independent check of the kernel program.** The Fourier theorem adds no
  independent check of the kernel program. It rests on the certificate, the generalised engine and the
  checker of the second result, and every limit stated here for the second result holds for the third.
- **New in the third revision: the repeated files.** The 11 repeated files are OpenAI's proofs with every
  declaration renamed by a script. The kernel checks the result as it is. That the copies differ from
  OpenAI's files only as their first comments and NOTICE say rests on the file comparison named above, not on
  Lean.
- **New in the third revision: the statement.** The statement of the third result is OpenAI's challenge file
  with three changes (section 1). `diff` shows them in a few lines (VERIFY.md, section 1). The one that
  carries the claim is the exponent.
- **New in the third revision: the runs.** The third result has had three comparator runs, two in the
  project's working tree and one in a copy of the repository. `tools/Compare.lean` passed in both places. The
  negative control (the checkers must reject this proof against a statement whose saving is one unit larger)
  was run in the working tree. On 9 October 2026 a build of this revision from nothing had not been made.
  One was made on 10 October 2026, on Linux ("Fourth result", above).
- The comparator ran on macOS with a stand-in for its Linux sandbox `landrun` (for the Linux machine of
  10 October 2026 see the limit "no run with the comparator's real sandbox yet", above). The second kernel
  (`nanoda`) was off, as in OpenAI's configuration of its own Fourier challenge. `leanchecker` was not run on
  these modules.
- In every comparator run the modules of the proof had been compiled beforehand, outside the sandbox, so the
  comparator replayed existing build products. Mathlib was not rebuilt from source, and the Lean toolchain and
  the comparator programs were the ones already installed. See VERIFY.md, section 8. One later run is the
  exception: on 10 October 2026, on the Linux machine, the comparator for the third result ran in a tree
  where none of the repository's modules had been compiled before, and it accepted the result at the third
  attempt (section 3).
- New in the second revision: on 2026-10-09 the second result had not yet had the rebuild from nothing that
  the first result had.
  The copy of the repository in which the comparator passed took the compiled files of its 25 certificate
  modules from the project's working tree ("Build in a copy of the repository", above), and the other two
  comparator passes were in that working tree. The build from nothing was made on 2026-10-10, on Linux
  ("Fourth result", above).
- New in the second revision: the certificate of the second result is checked by the kernel in segments. The theorems
  that join the segments are Lean proofs and are replayed by the comparator like everything else. The frame
  lemma, the generalised engine, the checker and its segmentation were all written on 2026-10-09; they are the
  larger part of what has not yet been read by a person.
- New in the second revision: the circuit of the second result is the outside authors' design, and it was
  rebuilt from their published data by Python programs that the agents made, in part by adapting the outside
  authors' own scripts ([NOTICE](NOTICE), section 6). That the rebuilt list of additions is the circuit #193
  describes
  rests on that Python conversion and on equal counts (9,412 helper arrays, the same histogram of block
  ranks), not on Lean. The theorem does not depend on this: the kernel checks the certificate as it is.
  The programs that rebuilt the circuit from the community's files are in `tools/rebuild/`. I added them
  after issue #1, whose author could not run the rebuild without them. With the community's files fetched
  at the two commits, they give the published certificate again, byte for byte once it is unpacked.
  `tools/gx/ORIGIN.md` says how to run them, and which two checking programs are still not here. The last
  step, the conversion to the certificate format, is `tools/gx/gxconv.py`.
- New in the second revision: the second result exists in whole-block accounting only (section 2, point 3), and the
  Lean checker covers the certificate format without its extension rule (unused by this certificate) and with
  every scratch copy at a frame of dimension at least 1.
- For the per-rank companion of the first result the comparator and `tools/Compare.lean` passed (the two runs
  above), but no negative control and no walk of the proof term were run. Its files were covered by the text
  checks and the arithmetic.
- A check that the statement is not vacuous was run for an earlier package (same text, other exponent), not
  for the published ones. It showed that `wht` at lengths 2 and 4 is the Sylvester matrix and that
  W(k) = 2^k k fails.
- The Python programs that searched for the first certificate, that rebuilt the second and that searched for
  the certificate of the fourth result are not part of the
  trusted base. The search programs were not audited. The kernel checks the certificates, not the programs.

**I welcome re-runs**, most of all on a Linux machine with Landlock, so that `landrun` really confines the
build, and with `nanoda` switched on. See VERIFY.md.

## 4. How it was made

I made the first three results on 8 and 9 October 2026, in about a day, and the fourth on 10 October, by
directing a team of AI agents (Claude). I began by
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

The third result came later the same day. The note on open directions had named it as a lead: carry the exponent
over
to the statement OpenAI headlines. I had the agents scope it before anything was built. They found that
OpenAI's proof uses the network through a single theorem and that the generalised engine already proves that
theorem's sentence at my exponent, and they recommended the route that leaves OpenAI's files unmodified and
repeats OpenAI's reduction in this repository. I took that route and asked for speed. The agents wrote the
seam file, repeated the 11 files by script, and ran the checks.

The fourth result came on 10 October 2026. With the checker of the second result a better circuit is a data
file, so I had the agents look for a better unit instead of a better proof.

A first round found that a family of 120 labels of width 9 can serve as a unit of this mechanism with 8
totals. That was in a wasteful form, with one helper for every unwanted term: 6,728 helpers, a price of
1.6555e-4. Three rounds of search on that family followed, all on the same day.

| round | what it did | helpers | figure |
|---|---|---|---|
| 1 | tried several architectures; the one that shares partial sums in a tree, chosen by a greedy rule under a group of 64 reorderings of the label positions, reached these figures | 1,028 | 7.5590e-4, just above the published 7.4745e-4 |
| 2 | added four devices (sources read directly, the re-use of finished helpers, totals built on used helpers, a re-scheduler) and a search over the greedy choices | 853 | 8.3202e-4 |
| 3 | replaced the greedy rule by an integer programme over a pool of 1,470 sums, with two weights chosen by the figure after the device passes | 783 | 8.7625e-4 |

The designs of rounds 1 and 2 were built in Lean as they came, and I locked the design of round 3 as the
certificate of this revision. The figures above, other than the last one, 8.7625e-4, are prices computed by
the Python tools; the notes say for each whether it was also Lean-checked.

Several routes failed, among them the transplant of the published architecture (5.0732e-4, also a price
computed by the Python tools) and a device-aware greedy rule. [notes/what-did-not-work.md](notes/what-did-not-work.md) tells them. I asked for
the search programs and that account to be published with the result, so that people and AI agents can see
how it was reached and go on from it.

The fourth revision also answers issue #1 of this repository, whose author could not run the rebuild of
the community's circuit because its scripts were not here. I had the agents tidy those scripts for
`tools/rebuild/` and compare them statement by statement with the community's own. The comparison showed
that six of the thirteen files are adaptations of the community's scripts. The earlier text called the
generator the agents' own. This revision corrects that and names the authors (NOTICE, section 6).

The fifth result came on 13 October 2026. By then the community had moved far past my fourth result, and Ryan
Shea had carried Jeremy Kothe's circuit into a Lean proof of OpenAI's statement through my engine, with Nil
Mamano's extension. I chose to build on that work rather than on my own circuit. I had the agents survey what
the community had published and find which steps a Lean proof must pay for while the community's figures assume
them away. Mutable sources was the cheapest of those steps to prove, so I put it first: the agents wrote the
proof, generated the certificate without the copies, and built and checked it on Colab machines (1.3553e-3).
Hosting was the largest step left, and the agents built it the same afternoon: two lemmas first, then the hosted
network, its checker and the instance, with a control that had to reproduce the figure without hosting exactly
(1.5699e-3). Exit frames, on the community's U8 circuit, came alongside: alone (1.2704e-3), then with mutable
sources (1.3869e-3). That evening the agents put the hosted network on the U8 circuit too, and the three steps
together gave the fifth result, 1.6198e-3. Other rounds of agents tried to beat the community's own figure
meanwhile and found nothing lasting.
[notes/fifth-result-story.md](notes/fifth-result-story.md) tells the whole day, with its failures and its
corrections, and `tools/community/` has the programs.

The sixth result came on 14 October 2026: the agents took the stack two steps further on the same circuit (no
source copies left, and a tighter bank table with a quartic rate bound), and it gave 1.6257e-3. It was still
below the community's own figure, so I kept it unpublished and turned to research. For five days I ran rounds of
agents on the question behind all of this: is there a small unit with a much larger saving, or a proof that there
is none? The moonshot waves of 14 and 15 October searched for one at blocks 2 to 5 in every game I could define;
other rounds asked how the known savings work, which lower bounds can be proved and checked by machine, whether a
unit can be shown to exist without being found, what other fields know, and, for two days, how to make real
transforms faster on processors and FPGAs. No unit was found. [research/](research/README.md) holds all of it,
with its programs, the claims that were refuted and the questions still open.

The lever came out of one of those rounds. On the morning of 14 October I had the agents ask whether two ideas
from fast matrix multiplication carry over: the laser method does not, but Schönhage's τ-theorem does, and with it
came the observation that a block at the very start of an array can be traded for a smaller job. The agents found
such heads in the community's circuit, proved Theorem 1G on paper and planned the Lean route, and I paused it for
the research. On 19 October the agents priced it on the community's newest circuits, checked those circuits with
the community's own programs, and found that the community's allocation programs could keep every head with one
added condition. I had the four proofs built in parallel on five machines.
[notes/how-the-lever-was-found.md](notes/how-the-lever-was-found.md) tells it, with what failed.

## 5. What is new in the second and third revisions, in plain words

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
- **The seam and the repeated reduction (third result).** OpenAI's route from the kernel program to every
  length asks one thing of the network: a kernel program within a work bound. In their files the bound is the
  fixed function `hills`, the ceiling of (k+1)^alpha. `Work/Fourier/Seam.lean` states the kernel program of
  the generalised engine in exactly that form, with the bound `hillsZ` at my exponent (`hillsZ_program`), and
  proves what the reduction uses of the bound: it is at least 1, it grows with k, it is at most k + 1 and at
  most 2 (k+1)^z, and z lies between 0 and 1 and below the exponents of OpenAI's statement. The 11 files
  after it are OpenAI's, repeated with that bound ("Whose reduction this is", above).
- **The new unit as a data file (fourth result).** The certificate of the fourth result is a file in the
  format of the extended checker, and the checker needed no change. It uses the format more freely than the
  certificate of the second result did. The scatter, the one step that hands the totals to the targets, is a
  table (in the second result: one fixed rule). There are 8 totals, all at the full frame of dimension 9 (22
  totals, two dimensions below the full frame). Gates are written in a form that adds arrays into one helper
  or target (in the second result that form occurs only for source into source). The sources are read by
  helpers before the scatter. Helpers are added into helpers before the scatter, with halves among the
  coefficients, and a helper that has served takes a new value. The agents first
  read the checker's Lean functions and found them general: one function for every addition, one for every
  move. The builds then had the kernel evaluate them on these certificates. The certificate of this revision
  has a few kinds of step that neither of the two earlier units built in Lean on 2026-10-10 had (a source
  added into a helper with the coefficient 1/2, in-place pairs with a half as the first coefficient, a helper
  that climbs from a line to the full frame in one block); the Lean build of this revision is the first
  kernel evaluation of those.

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
| fourth result | 9 | 120 | 45 | 783 | 1,263 | 120 | 56,715 of 56,835 | 3,974 blocks, 9,039 moves |
| second and third results | 22 | 1320 | 110 | 9412 | 14,692 | 3080 | 1,613,040 of 1,616,120 | 70,169 blocks, 262,944 moves |
| first result | 16 | 560 | 80 | 7122 | 9362 | 1035 | 747,925 of 748,960 | 37,286 blocks, 130,993 moves |

The blocks of one invocation are counted by the kernel, rank by rank. For the second result the two totals
are the Lean theorems `inv_count` and `inv_moves` (`Work/GCert/Data/P193Price.lean`); for the fourth they are
in `Work/GCert/Data/E8Price.lean`.

**The fourth result** keeps the bridged five-stage word and replaces the helper circuit once more, this time
by a circuit for another family of pairs. There are h = 9 coordinates and v = 120 pairs. The label of a pair
is a word of 9 bits with three ones (84 of them) or with seven ones (36). The agents named the family E8. Two
pairs know each other when their labels overlap in an even number of places; each pair knows 56 others. After
the scatter each target holds its own source and, with either sign, an unwanted half of each of the 56 sources
it knows: 6,720 unwanted terms, which the helpers have to take away. The wasteful way is one helper per term.
The circuit of the fourth result shares partial sums among the targets: a helper holds a sum of several
sources with signs, larger sums are made from two smaller ones, and where two targets want u + v and u - v
both are formed in place on the two helpers that held u and v. Which sums to form is chosen by an integer
programme over a pool of 1,470 candidate sums. Four passes then save helpers and blocks: a source is read
directly where a helper would only hold a copy of it; a helper that nothing reads any more takes a new value;
totals are built on helpers that have already served; and the frames and the order of the gates are chosen
again. The result is 783 helpers, 6.5 per pair, where the second result has 7.1 per pair, on labels less than
half as wide: the label space has dimension m = 5h = 45 in place of 110. Of the 480 moves that the 120
exchanges could save per unit, the scratch copies of the 8 totals take back 360 (section 7).
[notes/how-the-e8-unit-was-found.md](notes/how-the-e8-unit-was-found.md) has the whole account, step by
step with figures.

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
| fourth result (h = 9) | 480 | 360 | 75.0 % | 8.762479e-4 | not computed for this text |
| second result (community circuit, h = 22) | 5280 | 2200 | 41.7 % | 7.474547e-4 | 1.28147e-3 (x1.71) |
| first result (my circuit, h = 16) | 2240 | 1205 | 53.8 % | 5.399225e-4 | 1.16867e-3 (x2.16) |

In the unit of the fourth result the copies weigh more than in either earlier one: 8 totals at dimension 9,
cst = 72, so the copies take back 360 of the 480 moves the exchanges could save, and D = 120. One copy move
is about 4 percent of that saving. A paper argument with scripts, which nobody has read yet, says that 72 is
forced for 8 totals with the scatter used here and that no circuit for this family can go below 55; the range
55 to 71 is open. [notes/open-directions-2.md](notes/open-directions-2.md) has this and the other leads of the
fourth result. The rest of this section was written on 2026-10-09 about the first two circuits.

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

[notes/open-directions.md](notes/open-directions.md), written on 2026-10-09, hands the project on: every lead I
knew of then, with what is known, what it might give, how to start, and what already failed.
[notes/open-directions-2.md](notes/open-directions-2.md) adds the leads of the fourth result: what the search
of 2026-10-10 did not try and what may be promising.
[notes/open-directions-3.md](notes/open-directions-3.md) adds the leads of the fifth result: the steps the
community assumes that a Lean proof could still pay for (hosting across runs, the last source copies), with
their prices, what the rounds of 13 October 2026 mapped about the copies and the helpers, and what the agents
found about small units that evening and night, each finding with its status (the working notes behind it
are now in [research/](research/README.md), where later results supersede some of its statements).
[notes/open-directions-4.md](notes/open-directions-4.md) adds the leads of the lever: rules that are not
separated, heads with hosted helpers in front of them, the other blocks of rank 48, the full cut, and the next
circuits of the community. [research/open-problems.md](research/open-problems.md) lists the 122 open questions
of the research since 10 October, each with a first experiment.

- **Major leads.** The scratch copies above (what-if 1.28147e-3, no design known); a per-rank theorem for the
  generalised engine (its rate inequality at 4.058344e-4 is a Lean lemma already); the room under the
  outside ceilings (#201, DaysSky, claims 2.5657e-3 and covers neither this circuit nor this layout); the
  Lean frame lemma and the checker, for what the multiplication project still assumes. The Fourier transform
  of every length was on this list and is now the third result; the note says what remains of it (A3).
- **Minor leads** (fourteen, each with its computed size and its label), **the checks I would welcome most**,
  **seventeen things that were tried and did not work**, and **how a better circuit becomes a theorem here**.

## 8. Layout

    OAI/          89 files of openai/math (lean/OAI/Computability), unmodified
    WHTCheck/     the first Walsh-Hadamard corollary; defines `wht` and `WHTProgram`
    Work/         291 modules: 252 of the first three results (237 written here for the first two, and
                  Work/Fourier) and 39 of the fourth.
                  First result: recursion engines and rate arithmetic (Scratch, Block*, FoldRate), invocation
                  and network theorems (SharedSumStructured, Combine, Reframe, Bridge, BridgeGeom, Carrier),
                  certificate checkers, data and final theorems (SharedSumChecker, CarrierCheck).
                  Second result: GFrame (frame lemma, generalised engine, 54 modules, 48 of them imported
                  by the final proof), GCert (extended checker in the parts Labels, Scalar and Chain; Data
                  with the certificate and the final theorem)
                  Third result: Fourier (15 modules: the seam file written here, the modified
                  copies of 11 OpenAI files and of five definitions of a twelfth, the challenge, and a file that prints the axioms of the final theorems;
                  the final theorems are in the copy of `Main`)
                  Fourth result (39 modules, in addition to the 252): 24 under GCert/Data, generated from the
                  new certificate (namespace `GXD.E8`, instance `B2Ge8`; the final theorem is in
                  `SolutionB2Ge8x`), and FourierE8 (15 modules, the twin of Fourier with the suffix `E8`)
    comparator/   the comparator configurations
    community/    the fifth result and the results that came with it, as packages laid over Ryan Shea's
                  Lean tree (VERIFY.md, section 11): u8xmhh (the fifth result: the U8 circuit with exit frames,
                  mutable sources and hosting, with its host table, generators, logs and the two independent
                  checks), c588mhh (hosting on Kothe's circuit, with its host table, generators, logs, the
                  independent rebuild and the lemmas L1 and L2), c588m (mutable sources,
                  with its certificate, data, logs and two bank tables not built in Lean), u8xm and u8x (exit
                  frames), c588-repack, wht-c588 and wht-c588mhh (the Walsh-Hadamard corollaries), scripts
                  (apply, build and check, comparator, the leaderboard's packaging),
                  base-logs (the build of Shea's tree on the Colab, the agents' rebuild of the community's U8
                  circuit, and the scripts used there);
                  current.json names the headline package of 13 October
    package/      the sixth result and the lever, as packages laid over Ryan Shea's Lean tree with `scripts/`
                  (VERIFY.md, section 12): u8x0hq (the sixth result, laid first; every lever package requires it),
                  leverAFc (the full cut on #2979's words with the table A-fc, 1.7887045e-3; the result of the
                  title, added in the sixth revision),
                  leverF (the full cut on #3945 with every head kept, 1.7806935e-3; the fifth revision's headline),
                  leverB (the half cut on #3945 with every head kept, 1.7581520e-3),
                  leverA (the lever on #2979 with every head kept, 1.7564707e-3),
                  lever2979 (the lever on #2979 as published, 1.7556473e-3),
                  lever (the lever on #3945 as published, 1.7501109e-3),
                  each with its overlay, `MODULES.txt`, `result.json`, `SHA256SUMS`, logs, controls and generators
    scripts/      apply-overlay.sh, build-and-check.sh, run-comparator.sh (VERIFY.md, section 12) and the scripts
                  as the agents ran them (base/, cc-jail/, bos/, lever-run/)
    research/     the research of 10 to 20 October 2026 in twelve sections, with its index (README.md,
                  index.json, open-problems.md, make_index.py)
    data/, results/   the lever's inputs and recorded outputs: the community's #3945 words in canonical form and
                  their conversion into this project's certificate format; the head tests, prices and the
                  community's checks as re-run (`results/reproduced/` compares a fresh run with them)
    tools/        `Compare.lean`, the certificates as JSON (`tools/certificate/`), their reference checkers,
                  the generators (`tools/gen/` first result, `tools/gx/` second and fourth), `whatif_copies.py` (section 7),
                  `tools/fourier/` (the scripts that made the copies of the third and of the fourth result, and their checks),
                  `tools/rebuild/` (the programs that rebuild the certificate of the second result from the
                  community's data files; the data files are fetched, not included),
                  `tools/e8/` (the search programs that found the certificate of the fourth result),
                  `tools/community/` (the generators of the packages under `community/` and
                  `package_check.py`), `tools/swap-certificate.sh` (replaces the headline package of `community/`),
                  `tools/lever/` (the lever's search and pricing programs, with their README),
                  `tools/lever-hub/` (the programs shared with the research package for the community), and
                  `tools/reproduce-search.sh` and `tools/compare_results.py` (re-run them and compare with `results/`)
    lakefile.lean, lean-toolchain, lake-manifest.json    Lean v4.34.1 and the Mathlib commit pinned by openai/math
    VERIFY.md     how to rebuild everything and re-run every check
    ORIGIN.md     where every file comes from (with MANIFEST.sha256); `tools/gx/ORIGIN.md` for the tools and
                  the certificates of the second and fourth results
    RELATED-WORK.md, NOTICE, LICENSE
    notes/, third-party/    the notes of section 7 (scratch copies; all open directions) and the four notes of the
                  fourth result (how the unit was found, what did not work, open directions 2, the results as JSON);
                  the note for readers from network coding (`for-network-coding-readers.md`, with the folder `network-coding/`);
                  the two notes of the fifth result (`fifth-result-story.md`, `open-directions-3.md`);
                  the two notes of the lever (`how-the-lever-was-found.md`, `open-directions-4.md`);
                  the NOTICE file of CrocSwap/integer-mult-bounds (NOTICE, section 6), and its notices at the
                  commit the U8 circuit was rebuilt from (NOTICE, section 9, item g);
                  OpenAI's challenge file `UniformFourier.lean`, unmodified, to compare the challenges of the third and fourth results with;
                  unmodified files of Ryan Shea's repository that the generators of the fifth result use, with his
                  licence and notices (`third-party/shea256-fourier-transform-below-nlogn/SOURCE.md`);
                  unmodified copies of 32 files of the community's #3945 that the lever's programs read, with
                  their licence, notices and git blob hashes
                  (`third-party/CrocSwap-integer-mult-bounds/staged544-terminal-gap-fill/`, NOTICE, section 10)

The comments inside the Lean files were not edited for this release, so that the files are byte for byte the
ones that were checked. They still use path names and working titles of the tree they were written in ("Scratch
file", `checks/...`, "(key: ...)"). ORIGIN.md translates them.

## 9. Sources, dates and what was looked at

This section holds the dates, the commits and the outside pull requests behind the statements above. They are
kept out of the prose so that the prose can be read.

### For the opening and what is new

The six texts of this README, and one follow-up:

- First text (first result): commit `024f763`, 2026-10-09.
- Second text (second result): commit `fdfb781`, with its note on open directions in commit `e0bbe1c`,
  2026-10-09.
- Third text (third result): commits `6a1d04f` and `f010392`, 2026-10-09.
- Fourth text (fourth result): commit `5ba3046`, 2026-10-10.
- Follow-up to the fourth text: commit `9c94857`, 2026-10-10. It added the comparator verdicts on the Linux
  machine that arrived after the fourth text was published, and it corrected wording about the runs on that
  machine. It changed no Lean file and no figure of a result.
- Fifth text (fifth revision): commit `563b4e6`, 19 and 20 October 2026. Its part "The results of 13 October 2026" was written
  on 13 October 2026 for a fifth text that was not published; before the fifth revision, its headline had been moved
  on 2026-10-14 from `community/c588mhh` to the full stack, `community/u8xmhh`. Its sources are in "For the fifth
  result", below; those of the lever, the sixth result and the research are in "For the lever, the sixth result
  and the research".
- Sixth text (sixth revision): this one, 20 October 2026. It keeps the fifth text below its new opening and its
  section "What is new in the sixth revision"; its sources are in "For the sixth revision", below.

OpenAI's work that the statement and the cost model come from: [openai/math](https://github.com/openai/math),
family 130, commit `fd4aeeb`.

The searches of public sources:

- They were read-only looks by the AI agents I directed. RELATED-WORK.md says what they covered and when.
- For the first three results they were made on 2026-10-09, the last at 18:11 UTC. Every "as far as I could
  find" and every "as of 9 October 2026" about those results speaks of that day.
- They were not repeated in full for the fourth revision. What was read again on 2026-10-10 is listed next.

What was read on 2026-10-10 for the sentence "the largest saving with a Lean proof that I know of" (all times
UTC):

- At 20:49 UTC the agents read the titles and descriptions of the 100 newest pull requests of
  CrocSwap/integer-mult-bounds, #233 to #332, read-only.
- The largest figure there with a reported Lean check, for the quantity of my Walsh-Hadamard theorem (the
  "complex side", five-stage layout), was 7795/10^7, in
  [#327](https://github.com/CrocSwap/integer-mult-bounds/pull/327) (DreamingOfClouds, 18:51 UTC). Its author
  reports a Lean kernel check of it at 1 - 7795973/10^10, made with the generators of this repository. It
  was not built here.
- The largest kappa in a title was then 7.77268e-4
  ([#332](https://github.com/CrocSwap/integer-mult-bounds/pull/332), eumemic, 20:17 UTC). Kappa is the saving
  for integer multiplication, another quantity (section 2, point 2).
- At 21:33 the agents read the list of that repository's pull requests again. Two more existed, #333 and
  [#334](https://github.com/CrocSwap/integer-mult-bounds/pull/334) (both eumemic). #334 has the largest kappa
  in a title, 7.77596e-4 (21:15), and keeps #327 as its complex side. Neither states a Lean build.
- A closed pull request there, [#246](https://github.com/CrocSwap/integer-mult-bounds/pull/246)
  (maxime-fleury), prices 8.3465e-4 for the complex side and calls it "a conditional priced target, not a
  witness". It is a price, so the sentence does not cover it.
- Those figures move within hours.
- The sentence is about figures with a Lean proof only. Prices without one, in that repository or elsewhere,
  are outside it.
- Between 21:25 and 21:35 the agents read the outside repositories for the Fourier statement again,
  read-only. These readings built nothing.
- Pull request #2 of this repository (Chafik Boukhalfa, the account chafreaky; open, head `837dc31`) states
  1 - 7547360/10^10 for the Fourier transform of every length and 1 - 7547361/10^10 for the Walsh-Hadamard
  transform. Its author reports the Lean build and the comparator runs for both. It was the largest outside
  figure with a Lean proof for the Fourier statement. The agents built it on a Linux machine that evening,
  and every module of it compiled (the next list).
- shea256/fourier-transform-below-nlogn was read at commit `62554f9` (04:53). Since commit `a101b36` (03:16)
  its selection is 7.547360e-4, marked "formalized", in place of 6.7e-4. Its words: "an incorporated upstream
  result with a reproduced Lean proof" and "it does not claim a new exponent of our own". The rebuild is
  reported by its author and was not built here.
- danadran01/exact-dft-power-saving was still at `1a7b25e`, with 3.2e-6.
- eumemic/exact-dft-bounds was still at `6f87d1a`, with 4.856e-4, written. Its open pull request #1
  (chafreaky) has written figures up to 7.099e-4, "not formally verified".
- openai/math was still at `fd4aeeb`.
- xangma/exact-fourier-circuits, which this README and RELATED-WORK.md had not named before: a Lean proof for
  every length with theta < 1 - 2/10^13 and an extra log log factor, in a machine model of its own. That is
  the size of OpenAI's saving. It is reported by its author and was not built here.
- Also looked at: the forks of openai/math pushed since 9 October (not every branch), GitHub's repository
  search, arXiv by title and abstract, and Hacker News.
- Not looked at: X, Zulip, Discord, forums, other code hosts, and the branches of forks of the community's
  repository that are not pull requests.
- In all of this the agents found no Lean proof, reported or built, of a saving of 8.7625e-4 or more for
  either statement.

The build of pull request #2 of this repository, on 2026-10-10 (all times UTC). It is Chafik Boukhalfa's
result on his own hand-over pairs:

- The agents built it at its head `837dc31`, between 20:38 and 21:51, on the rented Linux machine of
  section 3 (x86-64, Ubuntu 24.04, 51 GiB of memory). They built it on top of that machine's from-source
  build of the published tree of this repository (commit `f010392`).
- All 54 Lean modules that the pull request adds compiled, with no error.
- The two heaviest modules, `Work.GCert.Data.Gen.P233Scal0` and `Work.GCert.Data.Gen.P233Scal1`, built in
  574 s and 642 s. The largest process of each held about 13 GB (12.7 GiB and 12.8 GiB). The 54 modules
  took 66 minutes in all, on a machine that other jobs shared.
- His three theorems, `wht_main_block_B2Gp233x`, `transform_mainY` and `convolution_mainY`, depend on
  `propext`, `Classical.choice` and `Quot.sound`, and on no other axiom.
- `tools/Compare.lean` printed `RESULT: PASS` for his Walsh-Hadamard challenge and solution and for his
  Fourier challenge and solution, with the constant counts that he states (38871; 43381 and 43404).
- Not run in that build: the official comparator. The agents ran it afterwards on
  `comparator/B2Gp233x.json` and `comparator/UniformFourier233.json`, the two configurations for which he
  reports runs. They used a second copy of that tree and the stand-in for the sandbox (section 3). The
  comparator accepted both, each with `Lean default kernel accepts the solution`, `Your solution is okay!`
  and exit code 0. Each log has 176 `Replayed` lines and no `Built` line, so no module was compiled inside
  the runs.
  - `comparator/B2Gp233x.json`: started at 22:58 on 2026-10-10, ended at 00:12 on 2026-10-11; 74 minutes,
    14.4 GB.
  - `comparator/UniformFourier233.json`: started at 23:00 on 2026-10-10, ended at 00:13 on 2026-10-11;
    73 minutes, 14.5 GB.
  - The runs shared the machine, so the times are upper bounds.
  - The memory is counted as in section 3: the kbytes figure that GNU `time` printed (14443212 and
    14501748), divided by a million.
  - These are not runs with the comparator's real sandbox. He reports his own runs, with the real `landrun`
    sandbox.
- Mathlib was not compiled there. Its compiled files came from Mathlib's cache.
- The Python checks of his certificate were made on the Mac earlier that day, not on the Linux machine
  ("Whose circuit this is").
- On the Mac (16 GB of memory) 36 of the 54 modules built. The build could not finish there, because the
  two heaviest modules need about 13 GB each.
- His files are not merged, so these three theorems are in his pull request and not in this repository.

The pull requests behind "Whose ideas the new unit uses" (all in CrocSwap/integer-mult-bounds unless marked):

- The in-place pair: icekylinx, [#184](https://github.com/CrocSwap/integer-mult-bounds/pull/184); carried in
  ikeboy's #191 and #193.
- The re-use of finished helpers: jamesyc, [#124](https://github.com/CrocSwap/integer-mult-bounds/pull/124),
  and eumemic, [#143](https://github.com/CrocSwap/integer-mult-bounds/pull/143).
- The late pairing rule: Chafik Boukhalfa (the account chafreaky),
  [#200](https://github.com/CrocSwap/integer-mult-bounds/pull/200) and
  [#233](https://github.com/CrocSwap/integer-mult-bounds/pull/233), and pull request #2 of this repository.

The outside figures in "What was new in the first three revisions":

- 7.009184e-4, the largest outside figure for the quantity of the second result when it was published:
  CrocSwap/integer-mult-bounds, pull request #193; read again at 18:11 UTC on 2026-10-09.
- 7.547361e-4: the price that `gxdry.py` computes for the same circuit with chafreaky's hand-over pairs, which
  are in pull request #233. His pull request #2 of this repository has that figure as a Lean theorem, which
  the agents built on a Linux machine (the list above).
- 7.635875e-4: the price that `gxdry.py` gives to the certificate of pull request
  [#304](https://github.com/CrocSwap/integer-mult-bounds/pull/304) (DreamingOfClouds), which states 7635/10^7
  for it. #304 is also the pull request that applies the five-stage formula and credits this repository for
  the layout.
- 3.2e-6: [danadran01/exact-dft-power-saving](https://github.com/danadran01/exact-dft-power-saving), an
  outside Lean proof of the Fourier statement; not built here. 6.7e-4: shea256, an outside written argument;
  conditional, not formalised. Both as the searches of 2026-10-09 found them. On 2026-10-10 shea256's
  repository gave 7.547360e-4 instead (list above).

### For the fifth result

**The text.** The part of this README on the fifth result was written on 13 October 2026. The commit of this
revision adds, with it, the packages under `community/`, `tools/community/`, `tools/swap-certificate.sh`, the
third-party folder of Ryan Shea's files and VERIFY.md, section 11.

**The base.**

- Ryan Shea's repository [shea256/fourier-transform-below-nlogn](https://github.com/shea256/fourier-transform-below-nlogn)
  at commit `52394849c7ecf3412d39b925d75fc7fdb98994d8`, folder `lean/complex588`. Its `lean/complex588` is the
  same at `9749945` (the commit Nil Mamano's leaderboard entry names) and at `0e61af6`. Lean v4.34.1, Mathlib
  `d13f23b723b8a846827a245b89c10fc7d3f11612`.
- Jeremy Kothe's circuit is pull request #588 of CrocSwap/integer-mult-bounds ("atm3"). The community's
  record of it (`certificates/checkpoints/family588.json`) gives kappa 0.001602278921983212; its complex side,
  kappa/(1 - kappa), is 1.6049e-3.
- Nil Mamano's gauged extension and Fourier templates: Shea's `NOTICE.complex588` says that the Lean files of
  the #401 layout come from Nil Mamano's published submission to the leaderboard (bigospeedrun/verifier commit
  `025396033eede56e3af6e02ed2679ba828381bd5`), unchanged apart from an import prefix, and that those files name
  "DFT Bridge3" and "DFT New Sums", and the record instance "DFT Lean". Shea's `lean/complex588/SOURCE_MANIFEST.json`
  lists `Work/Bank401/Adapt.lean`, `Work/Bank401/Frames.lean` and every file of `Work/Bridge401/` (among them
  `GNetDef`, `GNetStage`, `GGeom`, `GNetCert`, `SchedK` and `End`) as taken from
  `submissions/dft/nmamano-dft-record/Submission/` at that commit. The headers of those files credit the
  entrance gauges and the residual banks to Dugongue (#401), the bank patterns to lydakis (#275), the
  compressed nine-stage layout with its banks and padding to jacklightChen (#377) and the bridged word to
  DaysSky (#363). The entry of Mamano and Shea, `submissions/dft/nmamano-shea-c588` of bigospeedrun/verifier,
  was merged on 12 October 2026 (pull request #7 there).
- Shea's notice also describes what his conversion added: the "source-copy lowering of the mutable-source
  provider", the handling of scales, the cyclic bank packing, the generated checks and the Fourier instance. The
  fifth result removes that lowering.
- Zhihao Chen (jacklightChen): pull request #425 of CrocSwap/integer-mult-bounds, the written scalar-interface
  argument for mutable sources.
- Dugongue: pull request #401 there, the entrance gauges and the banks of the nine-stage network; Zhihao Chen
  (jacklightChen): pull request #377 there, its compressed nine-stage layout with its banks and padding.
- DaysSky: pull request #442 there, "idle-axis hosting", with its written proof
  (`research/idle-axis-hosting/PROOF.md`, sections 1 and 2, the hosting lemma) and the integer programme that
  chooses the hosts. The lemmas L1 and L2 (`community/c588mhh/supporting/`) and the hosting programs
  (`generators/design-tools/host_price.py` of `community/c588mhh` and of `community/u8xmhh`) cite its sections;
  the checked modules of `Work/XHost` cite the pull request. The files of `Work/XHost` that
  copy or adapt files of Shea's tree say so in their first lines: `HDef`, `HStage`, `HGeom` and `HCert` (from
  `Work/Bridge401/GNetDef`, `GNetStage`, `GGeom`, `GNetCert`), `HEnd` (from `Work/Bridge401/End.lean`),
  `Frames` (from `Work/Bank401/Adapt.lean` and `Frames.lean`) and `HSchedK` (lemmas about lists copied from
  `Work/Bridge401/SchedK.lean`).
- The U8 circuit of `community/u8x`, `community/u8xm` and `community/u8xmhh`: pull request #869 of
  CrocSwap/integer-mult-bounds (Jeremy Kothe, "eight triangle-hybrid B deliveries in block 0 on PR799's word"),
  on Kothe's #799, which carries the eight frame changes of Chafik Boukhalfa's #741, and on icekylinx's #645
  (the gauge44 construction). #869 credits chafreaky (Boukhalfa's account) for #741, #685 and #660 ("the
  frames, the owner map and the exact-flow packing"), icekylinx for #645, and "the PR601/604/617/622/627/635
  lineage and every upstream contributor credited in PR645"; the community's
  notice at the commit below (`research/complex-t21-owners/NOTICE.md`) says "PR869/PR870 supplies the
  unchanged U8 local word" and credits its triangle-hybrid construction to NewBArch, with its driver and
  composition to ColdWitness/icekylinx. The helpers that end early come into that line with, among others,
  #591, #592 and #601 (Kothe), #604 (DreamingOfClouds) and #645 (icekylinx). The charge for exit frames follows
  the community's endpoint lemmas in Lean, #593 (Kothe) and #599 (DreamingOfClouds).
- Not used by those certificates, though the community's figures for the U8 circuit use them: the owner maps
  of #870 (Kothe, applying Hao Yan's sixteen-replica controller of #858), #890 (Boukhalfa, "T21 whole-owner
  packing") and #942 (Jonathan Capalbo, the account jon314159, "T37 whole-owner packing"), the last two of
  which state that the U8 local word is unchanged, and the bit supplier #589 (DreamingOfClouds); #869 states
  "The complex side binds; PR589's bit supplier gives 0.0017334". I read the authors, titles and texts of
  these pull requests on 14 October 2026.
- The rebuild of the U8 circuit: CrocSwap/integer-mult-bounds, main branch at commit
  `407d52ff15b2217b6e22bce44728b4c3907b6515`, program
  `research/complex-t21-owners/vendor/pr870/construction/barch/reproduce_barch.py`, recipe U8 (#869's builder,
  in the copy of #870's package that #890 keeps), on Kothe's #799 X4 word; canonical sha256 of the word
  `24db880b7dab31e6af4eec8a1aba1bbbae29cd50d9a7b77cf0e9b571b3e5ae7c`, pinned by the community in
  `LOCAL-PINS-barch-U8.json` in the same folder at that commit. Logs: `community/base-logs/u8-word/`. The
  community's notices at that commit are copied in `third-party/CrocSwap-integer-mult-bounds/at-407d52ff/`.
- Ryan Shea's notes on exit frames: `research/complex-transfer-588-596/HANDOFF.md` of his repository at the
  commit above, section "Next mathematical test", which names as "the decisive follow-up" to generalise
  `GInvE` and its bank-frame construction from the full frame to the actual exit frames. The helpers that end
  early are the community's (its pull requests #596 and #601 by Kothe and #604 and #610 by DreamingOfClouds,
  which that note names).
- The community's selected results on 13 October 2026 were, in turn, #799 (complex side 1.6644e-3, in the
  morning), #942 (1.6647273e-3) and #1042 (in the afternoon).
- The community's list of assumptions: its README, read on 13 October 2026 at commit `f567884a`, names
  "common-label/mutable-source realization" and "whole-call hosting" among them. The newest figures quoted in
  "The results of 13 October in detail", item 7, are those of its selected result of that afternoon, pull request
  #1042 (kappa 0.001672188119205421; complex side 0.001674989015937672, reproduced by the agents' copy of its
  pricing program).

**The leaderboard.** bigospeedrun.com/exact-dft, read on 13 October 2026 (the site's data file
`data/exact-dft.json` at commit `d48436fa2c35` of bigospeedrun/site): first Ryan Shea, 1.2331685e-3, Lean
Verified; second Nil Mamano, 1.0221e-3; third this repository's fourth result, 8.7625e-4.

**The machines and the runs.** Every Lean run of the fifth revision's results was on two Google Colab
machines: Linux x86-64, Intel Xeon 2.20 GHz, 8 logical processors, 51 GiB, kernel 6.6.122+ without Landlock,
no systemd. On a Mac (Apple silicon, 16 GB) ran only the lemmas L1 and L2 that the hosting results do not
import, the two agents' rebuilds of them, the negative control of `hostSchedK` for `community/c588mhh` on an
earlier version of the data files, whose log is not in the package, and the first run of `hostSchedK` on the
tables of `community/u8xmhh` (`logs/generation/mac-FsPrice_HsU8XMHH*.log` there).

- The 32 modules of `community/u8xmhh`, its checks, the comparator and the offline re-creation of the
  leaderboard's verifier ran on the second machine, on 13 October 2026: the build and the comparator between
  22:05 and 22:42 UTC, the re-creation between 23:31 and 00:54 UTC (14 October), with bigospeedrun/verifier
  at `c92b9b6` and `community/u8xmhh/scripts-as-run/build/fsb_submission.py`, which packages the proof like
  Nil Mamano's entry. The comparator configurations are `comparator/fourier-u8xmhh.json` and
  `comparator/fourier-u8xmhh-neg.json`. The tree was a copy of a second build of `community/u8xm` that the
  agents made on that machine the same evening. The independent check of the statement ran on the same
  machine in its own copy of the tree; the independent rebuild ran on the first machine, over its own copy of
  the agents' first build of `community/u8xm` there, between about 01:00 and 01:50 UTC on 14 October.
- The agents' rebuild of the U8 circuit with the community's program ran on the first machine.
- The 49 modules of `community/c588mhh`, its checks and the comparator ran on the first machine. The
  comparator configurations are `comparator/fourier-c588mhh.json` and `comparator/fourier-c588mhh-neg.json`
  (and `fourier-c588mh*.json` for the control). The independent rebuild of its 32 modules ran on the same
  machine, in a build folder of its own. The lemmas L1 and L2 were compiled on a Mac.
- The packages `community/u8x` and `community/u8xm` were built and checked on the first machine.
- The 70 modules of `community/c588m`, the checks and the comparator ran on the first machine. The comparator
  configurations are `comparator/fourier-c588m.json` and `comparator/fourier-c588m-neg.json`; the sandbox
  stand-in is `community/base-logs/scripts/jail2.sh`.
- The offline re-creation of the leaderboard's verifier used bigospeedrun/verifier at commit `c92b9b6` (the
  preparation script, `community/scripts/bos/mb-bos-prep.sh`, stops on any other commit), on the submission
  folder packaged by `community/scripts/bos/mb_submission.py`, the script that reproduces the format of Nil
  Mamano's entry byte for byte. `community/scripts/bos/mb-bos-run.sh` ran the verifier's three phases (lean,
  audit, nanoda) one after the other, each inside the stand-in sandbox `jail2.sh` with `LEAN_NUM_THREADS=1`,
  because the verifier's own controller uses bwrap and systemd-run, which those machines lack. Its Lean phase
  took 108 minutes at 8.2 GB, the axiom audit 18 minutes and the nanoda phase 16 minutes. Logs:
  `community/c588m/logs/colab/bos/`.
- The Walsh-Hadamard corollaries at the kernels of `c588m` and `c588mhh`, and the independent rebuild of
  `community/c588m`, ran on the second machine.
- Tools at the pins of the verifier (`profiles/tools.json`): comparator `d03acab`, lean4export `076e8e5`,
  nanoda `3a24072`, landrun `811cfff`.

**The rounds of 13 October 2026** that the story note tells were run by the agents in the working folders of
this project; their figures and verdicts are in the story note and in `notes/open-directions-3.md`, each with
its status.

### For the lever, the sixth result and the research

**The text.** The opening, "What is new in the fifth revision" and the lever's entries elsewhere were written on
19 October 2026 and revised on 20 October 2026 for the fifth revision. The sixth result was built and checked on 14 October 2026
and rebuilt on 19 October 2026 (`package/u8x0hq/logs/`).

**The community's figures** (the table at the top), read on 19 October 2026 from the pull requests of
CrocSwap/integer-mult-bounds, read-only: #3945 at head `91f81d13` (complex root
0.001730663236534549309146...); #2979 on main at `20ecf289`, merged there as research (complex root
0.001729033982646513517717...); the highest complex-side construction among the open pull requests, #4088 (and #4100), 1.7309791e-3; #4091's
1.7309803e-3 is the optimum of a pooled (relaxed) model, not a construction; the chain
#3934 to #3959 (b 1.7411726e-3 in #3955, kappa 1.7367387e-3 in #3959), which carried the idle-axis moves of #3251
that their author retracted on 18 October after exact operator replays, was disputed by the open audit #4046 and
closed; the selected record #1094 (kappa 1.6726602e-3, complex side 1.6754627e-3); main at `bb4bccbf`.
A re-read on the evening of 19 October (the newest pull request was #4121) found no complex-side construction at or above the lever's; some research notes there print higher ceilings or targets, which are not constructions.
On the morning of 20 October (05:00 to 05:30 UTC) the community archived a large batch of pull requests on its main branch, among them #3945, #4091, #4085, #4076 and #4088, with an explicit review scope and none selected; #3945's files on main are the same as at its head. A re-read after that (main `12fa8c5c`, newest pull request #4167) found no complex-side construction at or above the lever's. #4181 (Rohan Arun, opened on 20 October at 11:22 UTC, open at the last read, newest pull request #4195) claims a fill at 1.7309803e-3 that attains #4091's pooled ceiling.

**The leaderboard.** bigospeedrun.com/exact-dft (the data file `data/exact-dft.json` of bigospeedrun/site), read on
19 October 2026: top entry 1.2654e-3, Bojnal's, on canghaimeng's C588AT, which runs on my engine.

**The community's checks re-run by the agents**, from scratch, on Colab machines (`results/community-checks/` and
the logs of each lever package):

- #3945: its complete stated verification `admission/full_replay.py` with fresh native replays
  (`PASS_TERMINAL_WIDTH5_COMPLETE_CONDITIONAL_RELEASE`, 119 of 119 controls rejected), its eight registered
  commands, the fresh native replay of both words (`PASS_COMPLETE_T74_NATIVE_REPLAY`) and the exact eight-stream
  operator replay over all 360 owners (0 differences); #3815's complete admission with fresh natives.
- #2979: the registered runner `scripts/research_package.py` (`PASS_REGISTERED_RESEARCH_PACKAGE`) and the complete
  replay `reproduce.py` with the author's recorded pins (`PASS_RANK7_JOINT_ENDFRAMES_COMPLETE_REPLAY`).
- The tables with every head kept: on #3945, its `owner/verify.py` (`PASS_LITERAL_ENLARGED_OWNER_FILL`, plain and
  `-O`) and `full_replay.py` (`PASS_TERMINAL_WIDTH5_COMPLETE_CONDITIONAL_RELEASE`, 119 of 119 controls), fill
  sha256 `9c06677d...`; on #2979, its `check_allocation.py` (`PASS_RANK7_LITERAL_OWNERS_AND_EXACT_MOMENT`) and
  `reproduce.py` (`PASS_RANK7_JOINT_ENDFRAMES_COMPLETE_REPLAY`), fill sha256 `764395fc...`. The exact operator
  replay over every owner gave 0 differences on both.

**The circuits' files.** #3945: `research/staged544-terminal-gap-fill/` at head
`91f81d131252db704f945891da1ea54e29d46b8b`; the 32 files the programs read are in `third-party/`, each checked
against its git blob hash (`BLOBS.txt`). #2979: `research/rank7-joint-endframes/` at main `20ecf289`; 12 of its
files are copied unmodified into the generators of `package/lever2979` and `package/leverA` (NOTICE, section
10.c), each checked by its hash.

**Credits.** The pull request bodies, `CREDITS.md`, `NOTICE.md` and native lineage notices of #3945, #2979 and
their predecessors (copies of #3945's in `third-party/`), the community's ledger of pull requests, and, for the
T74 line, the bodies of #1297 to #1839.

**The Lean builds and their skeptics.** The five lever packages were built and checked on 19 and 20 October 2026 (UTC) on five Colab machines by the integration agents `i-final` (`package/lever`), `x-final` (`package/lever2979`), `aA-final` (`package/leverA`), `bB-final` (`package/leverB`) and `fc-final` (`package/leverF`, round `wht45/fullcut`, with `fc-gate`, `fc-engine` and `fc-rate`), on Lean data generated by `aA-data`, `bB-data` and the earlier tracks of rounds `wht45/lever` and `wht45/lever2979`; the merged leaderboard form and its identity check are `h-bos`'s; the run of `scripts/` from a fresh clone of this repository is `r-finalize`'s. Each package's `logs/` and `result.json` give the
machine, date, times and peak memory of each step. The skeptics, three per package (the statement, an independent rebuild, the price): `v-statement`, `v-rebuild`, `v-price` (`package/lever`); `v2979-statement`, `v2979-rebuild`, `v2979-price` (`package/lever2979`); `vA-statement`, `vA-rebuild`, `vA-price` (`package/leverA`); `vB-statement`, `vB-rebuild`, `vB-price` (`package/leverB`); `vF-statement`, `vF-rebuild`, `vF-price` (`package/leverF`). All fifteen upheld their result; no fix round was needed.

**The bases.** Ryan Shea's tree as for the fifth result (below). OpenAI's statement:
`lean/ComparatorChallenges/UniformFourier.lean` of openai/math at commit
`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`, sha256 `4c571f27...`, copied in `third-party/openai-math/`. The
leaderboard's verifier: bigospeedrun/verifier at `c92b9b6`, re-created offline; its file limit of 500 read at main
`1e009d7`.

**The research.** [research/README.md](research/README.md) and each section's own Sources section name the notes,
rounds and skeptic verdicts behind it.

### For the sixth revision (the table A-fc in Lean)

**The text.** The opening and "What is new in the sixth revision" were written on 20 October 2026 for this
revision; the entries of the fifth revision were changed only where they said that A-fc was not proved in Lean.

**The table A-fc.** Made by round `wht45/fullcut` (agent `fc-reopt`, design `u1-2979sep`) on 20 October 2026;
`FILL.json.gz` sha256 `7fc766c6c3e1561cee60b258138a918020b96d85ea782a2b01ad0ee035fd082e`, `PROFILE.json` sha256
`5df79411...`. Its admission by the community's checks, its prices and its first skeptic: round `wht45/hubx`
(`hx-admit`, `hx-price`, `hx-skeptic`), on a fresh clone of CrocSwap/integer-mult-bounds at main `12fa8c5c`; #2979's
own files at main `20ecf289` (read at `bb4bccbf`).

**The Lean build and its skeptics.** Round `wht45/afc`, 20 October 2026 (about 10:20 to 16:05 EDT), on four Colab
machines: the data layer `afc-data`, the gate of the cut fact on #2979's words `afc-gate`, the full-cut instance and
the rate fact `afc-cut`, and the integration `afc-final` (the final run, the official comparator and the release
controls on C4; the merged re-creation of the leaderboard's verifier, a fresh build from the package's files and
the identity check on C5; the unmerged re-creation on C2), with `afc-prep`'s machine setup. The skeptics:
`vAFc-statement` (C2), `vAFc-rebuild` (C3) and `vAFc-price` (C4); all three upheld the result, and the statement
skeptic found the Python 3.12 property of the generators. The package as checked is the staging tarball of the
package, sha256 `fcad92871cc93b83fa0ad71bf91ef5b3a47bec364eecb4f361ebcafdaeb1f336` (its overlay
`061188433b0d...`, the merged leaderboard stage `1cb208e5d920...`). The copy in this repository differs from it in
text files only, none of them under `overlay/` or `controls/`: machine paths shortened in eight files (two price
records, `generators/PATCHES.txt`, the two lowering reports, `logs/afc-gate/GATE.md` and the two logs of
`logs/regenerate/`), `README.txt`, `KNOWN-ISSUES.txt` and `result.json` (the controls listed apart), and the
generators' `README.txt`, `PATCHES.txt`, `REGENERATE.sh` and `tree.sha256`; `SHA256SUMS` lists the files as
published. The run of `scripts/` from a fresh clone of this revision is agent `fu-prep`'s (C3).

**The regeneration scripts.** The two lowering reports
`generators/tree/wht45/lever2979/x-gate/out/lower-base.json` and `lower-extra.json` had sha256 `f219f479...` and
`5e5d7945...` as the run wrote them; `generators/PATCHES.txt` of `package/lever2979`, `package/leverA` and
`package/leverAFc` gives the published values. The Python property: `x-words/tools/pack_ceil2.py` and the change to
`sum()` over floats in Python 3.12; the 36 files are `Work/GCert/Data/Gen/S2979BRows0..17.lean` and
`S2979ERows0..17.lean`.

**The community's figures**, read on 20 October 2026, read-only: by the skeptics between 15:00 and 15:14 EDT (all
pull requests and issues updated since 17 October, with their bodies and comments; newest #4221 and #4223; main
`12fa8c5c`), and again by me at 16:36 EDT (the pull requests opened since, up to #4243; main `12fa8c5c`). No
complex-side construction claimed by others was at or above 1.7887045e-3. The highest construction claims were
#4181's kappa 1.7309803e-3 (Rohan Arun, open), #4100's 1.7309791e-3 and #3945's b 1.7306632e-3; the chain #3934 to
#3959 (b 1.7411726e-3) was withdrawn. Higher figures there are ceilings, caps, targets or hypothetical ledgers (for
example #3729's 4.7377e-3, a ceiling), not constructions.

**The leaderboard.** bigospeedrun.com/exact-dft, read on 20 October 2026: top entry 1.2654e-3, as on 19 October.
The leaderboard's verifier was re-created offline at `c92b9b6`, as for the fifth revision; its main was `7f4f6bd` on
20 October.

### For the credits

What "as far as I could find", "that I know of" and "as of 9 October 2026" rest on, in the two sections "Whose reduction this is" and
"Whose circuit this is":

- The AI agents I directed looked at public sources, read-only. RELATED-WORK.md says what they covered and
  when.
- The last of those looks for the first three results was on 2026-10-09 at 18:11 UTC.
- They were not repeated in full for the fourth revision. The repositories below and the newest pull
  requests of the community's repository were read again on 2026-10-10 ("For the opening and what is new",
  above).

The three outside repositories for the Fourier statement:

- danadran01/exact-dft-power-saving: read at commit `1a7b25e` (2026-10-09 03:48 UTC). It was not built here.
- shea256/fourier-transform-below-nlogn: first commit 8 October 2026. Read on 2026-10-09 at 18:11 UTC, at
  commit `a2840b1`; its figure was 6.7e-4 then.
- eumemic/exact-dft-bounds: read at commit `6f87d1a`.
- All three were read again on 2026-10-10, between 21:25 and 21:35 UTC. The first and the third had no new
  commit. shea256's was at commit `62554f9` and gave 7.547360e-4 ("For the opening and what is new", above).

The community's repository, CrocSwap/integer-mult-bounds:

- #144 (icekylinx, the paired-cube construction) has been in its main branch since 2026-10-09 05:06 UTC.
- The pull requests behind the circuit of the second and third results: #144, #168, #184, #191 and #193.
  Credited through them: #117, #128, #162 and #163.
- The other choice of hand-over pairs: #233 and #200 there, and issue #1 and pull request #2 of this
  repository.

### For the results, the limits and the checks

**The readings of the community's pull requests** (section 2, point 2). The repository is
CrocSwap/integer-mult-bounds. All times are UTC. Kappa is the saving for integer multiplication; the
complex-side saving is the quantity of my theorem.

- 2026-10-09, 16:06. Pull requests up to #207 existed.
  - Largest kappa in a title: 6.831905e-4 (#207, Dugongue, open, 16:00).
  - Largest complex-side saving stated in any of the bodies of #193 to #207: still 7.009184e-4, the circuit
    of #193. #194, #197, #202, #204, #205, #206 and #207 keep that circuit while they change the bit side.
  - The front page of that repository's main branch (`3b6b668`) named #186 as its reviewed result (kappa
    6.61885549e-4).
- 2026-10-09, 18:11. Pull requests up to #215 existed.
  - Largest kappa in a title: 6.839217e-4 (#210, eumemic, open, 16:58).
  - Largest complex-side saving stated as a result in the bodies of #208 to #215: still 7.009184e-4, the
    circuit of #193, which #210, #211 and #213 keep.
  - #208 prices further steps, up to kappa 1.226488e-3. Its title calls them "target, not built".
- 2026-10-10, 20:49. Pull requests up to #332 existed. The 100 newest, #233 to #332, were read.
  - Largest kappa in a title: 7.77268e-4 (#332, eumemic, open, 20:17).
  - Largest complex-side saving with a reported Lean check: 7795/10^7 (#327, DreamingOfClouds, open, 18:51),
    for a word at h = 20 in the five-stage layout of this repository. Its author reports a Lean kernel check
    of it with this repository's generators. It was not built here.
  - A larger complex-side figure is priced in #246 (maxime-fleury, closed): 8.3465e-4, which its own text
    calls "a conditional priced target, not a witness". It states no Lean check.
- 2026-10-10, 21:33. Pull requests up to #334 existed. The list was read again, with the descriptions of #246,
  #327 and #334.
  - Largest kappa in a title: 7.77596e-4 (#334, eumemic, open, 21:15). It keeps #327 as its complex side and
    states no Lean build.
  - Largest complex-side saving with a reported Lean check: still 7795/10^7 (#327).

**The runs of the first three results** (section 3). All were made on 2026-10-09, on the Mac. Times are UTC;
the seconds are the duration of the run.

- First result, comparator:
  - by the agent that built the proof: whole-block theorem 09:59, per-rank companion 10:03;
  - by an auditing agent with its own configuration: whole-block theorem 10:09;
  - in a rebuild from scratch: whole-block theorem 11:32, per-rank companion 11:35.
- Second result:
  - comparator in the project's working tree, by the agent that built the proof: 15:14 (1814 s);
  - comparator in the project's working tree, by an auditing agent with its own copy of the configuration:
    15:43 (1116 s);
  - `lake build --no-build` of all eight roots in the build copy of the repository, every target up to date:
    15:48 (9250 jobs);
  - comparator in the build copy: 16:13 (1058 s);
  - `tools/Compare.lean` in the build copy, `RESULT: PASS`: 16:15 (148 s).
- Third result:
  - comparator in the project's working tree, by the agent that built the proof: 18:27 (1528 s);
  - comparator in the project's working tree, by an auditing agent with its own configuration: 18:55
    (1236 s);
  - comparator in a copy of the repository: 19:01 to 19:20 (1137 s);
  - `tools/Compare.lean` in that copy: 19:23.

**The first result** was published on 2026-10-09 as commit `024f763` (section 1).

**The Walsh-Hadamard statement, line by line** (section 1, "Whose statement the Walsh-Hadamard one is").

- Lines 22 to 270 of my challenge file are the RAM cost model. They are OpenAI's, byte for byte: lines 9 to
  257 of their `lean/ComparatorChallenges/UniformFourier.lean`.
- Lines 271 to 331, about 60 lines, define `wht`, `WHTProgram` and `WHTTimeBoundsAt`. They are my project's.

**The fourth result.** It was found, built and checked on 2026-10-10. Its runs, with commands, times and
memory, are in [VERIFY.md](VERIFY.md), section 10.

**The comparator for the third result in a tree compiled from nothing, on the Linux machine** (section 3).
All times are UTC, on 2026-10-10.

- The first attempt started at 20:40 and the second at 21:17. Neither gave a verdict: the agents stopped
  both during the build, to keep the shared machine within its memory.
- The third attempt started at 21:23 and ended at 23:30, with `Your solution is okay!` and exit code 0.
- That was after the fourth text was published, so the fourth text said that the run had no verdict. The
  follow-up text has the verdict.
- The exact time and memory are in [VERIFY.md](VERIFY.md), section 8, point 14.
