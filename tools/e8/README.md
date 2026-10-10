# tools/e8: the programs that found the E8 certificate

These are the programs that found `tools/certificate/gcert1-e8-r783.json.gz`. That file is a *certificate*: a
data file that lists every step of a *unit*, the small fixed program that the whole method repeats at every
scale. This unit is on the E8 label family (label width 9, 120 pairs of data arrays, block size 45) and has 783
helper arrays, which are extra working arrays. Its figure is 8762479, that is a saving of 8.7625e-4 for the
Walsh-Hadamard transform. The section "Words used here", below, explains these words.

One command, run from the top folder of the repository, rebuilds that file byte for byte from the data in this
folder and checks it with the repository's own checkers:

    python3 -B tools/e8/run.py

- **Time:** 98 s on the machine used, an Apple-silicon Mac.
- **What you need:** Python 3.9 or later and its standard library (tested on 3.9.6 only). Only the longer run with `--pool --solve` (section 1)
  also needs numpy and scipy.
- **What "reproduced" means:** the run ends with `RESULT: REPRODUCED` and exit code 0 only if two things hold.
  The JSON it wrote is byte for byte the JSON of the certificate named above. And all three checkers accept it:
  `tools/gx/refcheck.py` with 783 helpers, and `tools/gx/gxdry.py` and the stand-alone `replay.py` of this
  folder at the figure 8762479. Anything else ends with a line starting `FAILED` and exit code 1. Section 3 says
  what was tested and where it can break.

I found the certificate on 10 October 2026 with a team of AI agents (Claude) that I directed. I chose the
directions of the search; the agents wrote the programs here and ran them.

**Credit.** Three devices in these programs come from circuits that other people published in the public
repository CrocSwap/integer-mult-bounds, whose contributors work on bounds for integer multiplication. The
solver is not mine either. "There" below means that repository.

- **jamesyc (pull request 124 there) and eumemic (pull request 143 there)**: slot reuse with a compensating
  read. A copy of a sum is written on top of a helper that is already dead, that is, a helper that nothing reads
  any more, instead of on a new helper. The dead helper's old content is still in it and travels with the copy.
  Each target that ends up with that old content reads that helper once, with the opposite sign, which cancels
  it. Here: `comp2.py`, `reuse2.py`.
- **Chafik Boukhalfa (account chafreaky), pull requests 200 and 233 there and pull request #2 of this
  repository**: the late pairing, that is, which dead helper a new value takes, chosen late. Here it is the rule
  in `reuse2.py` that picks the dead helper. It is a greedy rule written for this search; his program, his pairs
  and his certificate are not used.
- **icekylinx (pull request 184 there), carried in ikeboy's pull requests 191 and 193 there**: the in-place
  pair. Two arrays hold a and b, and a + b and a - b are formed on those same two arrays. Here: `comp.py`,
  `comp2.py`.
- **HiGHS**, the open-source solver that `scipy.optimize.milp` calls, solves the integer programme of `xmilp.py`.

**What kind of figure this is.** Everything this folder prints is a price computed by the Python tools
(`tools/gx/refcheck.py`, `tools/gx/gxdry.py`). Nothing here is a proof.
The Lean check of this same file is the theorem `wht_main_block_B2Ge8x` of this repository (VERIFY.md).

## Words used here

- **unit**: the small fixed program that the whole method repeats at every scale; its size and its saving
  decide the final exponent. The unit here is called E8, after its label family.
- **pairs, source, target**: the working copies of the input that the unit must transform come in pairs. In
  each pair one array is the *source* and the other the *target*. E8 has 120 pairs.
- **label, label family**: a label is the string of bits that names a pair, and the family is the set of all
  120 labels. A label's length is the *width* (9 for E8). The *block size* is the number of directions along
  which every array must be transformed (45 for E8, five times the width).
- **helper array**, or helper: an extra working copy that the unit borrows, uses and gives back unchanged. Here
  each one holds one partial sum of sources for a while. Helpers cost time and save nothing themselves, so fewer
  is better. A *dead* helper is one that nothing reads any more.
- **saving** (written like 8.7625e-4): how far the exponent of log n is below 1; larger is better.
- **figure**: the same saving as the whole number that the checking programs print. The figure 8762479 means
  the exponent 1 - 8762479/10^10, that is the saving 8.7625e-4.
- **price**: the saving that the Python checking programs compute for a certificate; it becomes a theorem only
  when the proof checker Lean has built it.
- **frame, stop**: the frame of an array is the subspace of labels where it stands; it only grows. A *move*
  makes it one dimension larger, and moves are what cost time. An array makes its moves in runs called *blocks*,
  with a *stop* between two blocks.
- **scatter, totals, copy moves**: at one point of the unit, the *scatter*, a few arrays called *totals* (8
  here) are added into all targets. The copy moves are the summed dimensions of the totals' frames. They take
  back part of what the unit saves.
- **shared sum, hierarchy, pool**: a shared sum is a partial sum that a helper holds. The *hierarchy* is the
  list of shared sums and which target takes which. The *pool* is the list of candidate sums it is chosen from.
- **pass**: one of the five programs that, one after the other, turn a hierarchy into a certificate (steps 1
  to 5 of section 2).
- **round**: the search went in three rounds. The steps of section 2 are the third; the first two produced the
  four hierarchies that the pool is collected from (`data/pool-sources/`).

## 1. How to run it

| command | what it does | needs | time on the machine used |
|---|---|---|---|
| `python3 -B tools/e8/run.py` | starts from the recorded hierarchy, runs the five passes, compares checksums, runs the three checkers | Python 3.9 or later, standard library | 98 s |
| `python3 -B tools/e8/run.py --pool --solve` | first rebuilds the pool of 1,470 sums and re-runs the solver, then the same | also numpy and scipy | 107 s, peak memory 258 MB |
| `python3 -I -B tools/e8/replay.py <certificate> e8` | the stand-alone replay alone | standard library | under 1 s for this file |

`--work DIR` chooses the working directory (default `./e8-run`); every pass writes its whole output to
`DIR/logs/`. The run ends with `RESULT: REPRODUCED` and exit code 0, or with a line starting `FAILED` and exit
code 1. The machine used is an Apple-silicon Mac with Python 3.9.6, numpy 2.0.2 and scipy 1.13.1.

## 2. The path, step by step

A *hierarchy* is the list of shared sums and which target takes which. Each of the 120 targets has to lose 56
unwanted terms; the hierarchy says which partial sums of sources are formed, from which two smaller sums, and
which target reads which sum. A *helper* is an array that holds one such sum for a while. The figure grows as
helpers and stops go away.

| step | program | what it does | helpers | figure after it |
|---|---|---|---|---|
| 0a | `pool1.py` | collects the sums of four hierarchies of the two earlier rounds (`data/pool-sources/`) and their images under a group of 64 symmetries: 1,266 sums, 1,470 with their images (`data/pool-a.json`) | | |
| 0b | `solve.py`, `xmilp.py`, `xlib.py` | integer programme: choose a whole hierarchy from the pool. It minimises uses + ordinary sums (a use is one target reading one sum; an ordinary sum is one formed alone, not as half of an in-place pair); an ordinary sum with a single-source part counts 0.55 and one whose smaller part has two sources counts 0.85, because steps 3 and 4 can remove the helper such a sum costs. No symmetry is imposed. HiGHS reports the optimum 886.3 with gap 0. Output: `data/hierarchy-l0.55-m0.85.json.gz` | | |
| 1 | `prune.py`, `ilib.py`, `comp.py` | compiles the hierarchy into a circuit (which helper holds what, at which frame; in-place pairs where both sums of a pair are wanted) and drops the target stops that do not pay | 927 | 8125886 |
| 2 | `reuse2.py`, `comp2.py` | slot reuse: 52 copies are written on dead helpers | 875 | 8261680 |
| 3 | `xsub.py` | a helper that only ever holds one source is replaced by the source array itself | 819 | 8552172 |
| 4 | `xexp.py 6` | a helper that lives only before the scatter is replaced by direct reads of the source arrays it held (parts of up to 6 sources); the 8 totals first serve as such helpers | 783 | 8686103 |
| 5 | `sched.py 1 60 1.5` | the same additions at other frames: local search with annealing over where each array stops | 783 | **8762479** |

The figures of steps 1 to 4 are prices computed in process (`e8.quick`, which calls the two published checkers
as functions). The last one is also the verdict of `refcheck.py` and `gxdry.py` run as programs, and of
`replay.py`.

Shared modules: `e8.py` (the E8 family, the builder that turns a design into a certificate file, the in-process
check), `lib.py`, `gram.py` and `groups.py` (the symmetry groups; `gram.py` also holds the greedy of the first
round), `hio.py` (the JSON formats of hierarchies and pools, described in its first lines).

## 3. What "reproduces" means here

What was run, and what it gave:

- From a clean copy of `tools/` alone, `run.py` wrote a file whose uncompressed JSON has the sha256
  `3594f19c4d01cf11fb930c5c61baeed399620acbc5b94096b16bcacb23997dd3`, the checksum of the E8 certificate;
  the `.gz` bytes were the same too. `refcheck.py` printed `ACCEPTED`, `gxdry.py` and `replay.py` gave 8762479.
- With `--pool --solve` the pool came out equal to `data/pool-a.json`, the solver returned the recorded
  hierarchy, and the certificate was again the same bytes. This was run once when the folder was made, and once
  more by a checking agent on the committed files.
- Three damaged copies of the recorded hierarchy (one coefficient with the wrong sign, one sum formed from the
  wrong part, the file cut in half) each ended with `FAILED` and exit code 1.

What makes a run repeatable, and where that can break:

| stage | seeds and limits | what could change the result |
|---|---|---|
| solver | time limit 30 s (the solve takes 9 s); costs not disturbed (tie-break seed 0); gap 0 | many optima have the same value. Another HiGHS version may return another one, and then the later steps give another certificate. Other optima of equal value differed by up to 0.26 % in an estimate of the final figure during the search |
| `prune.py` | none | comparisons of sums of logarithms with a tolerance of 1e-9 |
| `reuse2.py` | 8 fixed plans and 20 seeded ones (seeds 0 to 9), Python's `random.Random(seed)`; stops trying plans after 400 s of wall-clock time (needs 20 s) | a machine about 20 times slower would stop early |
| `xsub.py` | none | |
| `xexp.py` | no randomness; stops after 420 s (needs 13 s) | a much slower machine |
| `sched.py` | seed 1, 60 annealing sweeps, start temperature 1.5; stops after 400 s (needs 48 s) | a machine about 8 times slower; another floating-point library could round a tie differently |

All of this was tested on one machine and one Python version. On another platform the checksum may differ.
The three checkers then still say whether what you got is a valid certificate and what its figure is; only the
claim "the same bytes" would fail, and `run.py` says so.

What `run.py` does not re-run:

- **How the two weights were chosen.** Eleven settings of the two weights were solved and valued by a fast
  estimate of the final figure, and four of them went through steps 1 to 5 (figures 8740087, 8746256, 8753727,
  8762479; with the first weight alone, at 0.4, the figure was 8716750). The scan programs are not in this folder; `solve.py` takes
  the two weights as arguments, so any setting can be tried.
- **The two earlier rounds.** The four hierarchies in `data/pool-sources/` are recorded outputs of the greedy
  searches of rounds 1 and 2 (shared sums chosen by whole orbits of a symmetry group; then single greedy choices
  changed). The greedy itself is in `gram.py`; the round-2 program that changes choices is not included, and
  re-running those searches from this folder is untested.

## 4. The stand-alone replay, `replay.py`

A second opinion on a certificate file. It imports nothing from `tools/gx` and reads only the certificate. It
replays the circuit, checks every rule of the two published checkers as their code states
them (each rule cites the line it was read from), and computes the figure with its own interval arithmetic.
It does not check Lean, the generators, or whether the rules themselves are the right ones: it shows that the
file obeys the rules, by a program that shares no code with the checkers. It refuses files with a non-empty
`ext` field.

## 5. Data that documents the search

`data/scoreboard.tsv` (11,492 bytes): all 61 certificates registered during the search, in order, with helpers, copy moves,
figure and checksum. Every figure in it is a price by `gxdry.py` for a file that both published checkers
accept. Two of them (7558959 and 8320202) were also Lean-checked during the work; their Lean files are not
included. Only five of the 61 files are in this repository: the certificate and the four below.

`data/others/`: one certificate for each approach that did not win. All four are accepted by `refcheck.py` and
`gxdry.py` at the figure in their name.

| file | bytes | approach | helpers | figure | times 7474547, the figure of the certificate published before this one |
|---|---|---|---|---|---|
| `covering-R2366-3923956.json.gz` | 25,176 | every target served by pieces of complete bicliques, no shared sums | 2,366 | 3923956 (3.9240e-4) | 0.52 |
| `transplant-R1670-5073193.json.gz` | 24,089 | the architecture of the certificate published before this one, moved to E8 | 1,670 | 5073193 (5.0732e-4) | 0.68 |
| `triangles-R1475-5535177.json.gz` | 24,779 | triangles on top (one 27-source sum serves two targets), shared sums below | 1,475 | 5535177 (5.5352e-4) | 0.74 |
| `triangle-top-R1159-6805162.json.gz` | 20,226 | the triangle top with the symmetric greedy and in-place pairs below, then `xsub.py` | 1,159 | 6805162 (6.8052e-4) | 0.91 |

In these four files the field `derived_from` was rewritten for publication (it held working names); the
circuits are unchanged, and the scoreboard gives both checksums.

One working name is left in the certificate itself: its field `derived_from` reads
`i-stops: re-scheduled c-oh-l0.55-m0.85-s0-xexp.json.gz (seed 1)`, the name the last pass had during the search
and the name of its input file. `sched.py` and `run.py` keep both, because any other text gives the same circuit
with another checksum.
