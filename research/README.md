# Research since the fourth revision (10 to 20 October 2026)

This folder holds all the research I did between the fourth revision of this repository (10 October 2026) and
this one, apart from what is already public elsewhere. It is written so that people and AI systems can use it: for
each line of work it says what I tried, what worked, what failed and what is still open, with the programs that
produced every result, the commands to run them, and a machine-readable summary.

I did this work with a team of AI agents (Claude) that I directed: I set the goals, chose which ideas to pursue and
decided at each stage what to do next; the agents ran the searches, wrote the proofs and programs, and tried to
refute one another's work. **No person has read any of the proofs in this folder.** Where a section says "a
skeptic", it means one AI agent whose only job was to break a result with its own reasoning and its own code; it is
not a human referee. Most results had one skeptic, some two or three, some none, and each section says which.

Already public, and not repeated here: the 5/6 circuits for the exact Walsh-Hadamard transform and the ten idea
waves of the night of 10 to 11 October 2026 are in
[jacobalansussman/wht-fewer-operations](https://github.com/jacobalansussman/wht-fewer-operations); the notes of the
fourth revision and earlier are in [`../notes/`](../notes/).

## What came out, in one screen

- **The record moved, by a change of recursion.** The search for a better small unit found nothing, but the work
  on why the community's circuits save led to the *recursion lever*: many idle blocks of the community's circuits
  are the very first operation on their array, so the recursion can give those arrays a smaller job. The repository's
  [README](../README.md) has the result; [record-track](record-track/README.md) has the theory (Theorem 1G), the
  pricing programs, the checks before each build, and every lever that failed before it.
- **No new unit was found**, at any open shape, in any game, not even at the stepping-stone tier of a saving of
  0.01. The moonshot waves of 14 and 15 October ([small blocks](moonshot-small-blocks/README.md),
  [blocks 4 and 5](moonshot-block4-block5/README.md)), the Reed-Muller family and the one-move units
  ([reed-muller-and-one-move](reed-muller-and-one-move/README.md)) and the two-array designs
  ([frame-packing-two-arrays](frame-packing-two-arrays/README.md)) all came back empty.
- **No existence proof reached an open shape** ([existence-routes](existence-routes/README.md)), although a unit
  known only to exist would count in full (proved in Lean).
- **The main gains are lower bounds and an understanding of the mechanism.** Proved in Lean: the counting rule
  behind every lower bound (Theorem R), that block 2 never saves, and that block 3 needs at least 2.8 paid moves per
  array, so a block-3 unit saves at most 6.28e-2 ([block-bounds-and-new-rules](block-bounds-and-new-rules/README.md)).
  Checked by computer: block 4 needs at least 41/13 moves per array (a saving of at most 0.171). Argued by hand, not
  checked: block 2 never saves in any game whose free steps touch boundedly many entries per row (Theorem 2W), and
  every unit with two or more arrays needs at least W + m moves. Also argued by hand, not checked: every saving of
  every unit is *frame packing* ([frame-packing-two-arrays](frame-packing-two-arrays/README.md)), and the record
  circuits all save in one way ([why-the-savers-save](why-the-savers-save/README.md)).
- **The methods have limits, now stated** (each argued by hand, not checked). No hook, cut or coordinate
  certificate can rule out a saving of 0.1 from block 5 on; logic gives no shortcut to an existence proof
  ([theory-of-theories-and-logic](theory-of-theories-and-logic/README.md)). The problem is not a famous problem in
  disguise, and the dictionaries to quantum circuits and network coding are textbook
  ([other-fields](other-fields/README.md)).
- **Speed on real machines** ([speed-waves](speed-waves/README.md)): one application-level gain, a number-theoretic
  transform inside Microsoft SEAL that made a loop of real encryption operations 21.2 % faster on an Apple M4, and a
  list of kernel gains and failed leads.
- **Failures are told as plainly as successes.** Across the twelve sections, 92 findings are marked refuted, each
  with what refuted it. When two sources disagreed, the later skeptic verdict is the one used, and the section says
  so.
- **What is open**: 122 questions, each with a concrete first experiment, collected in
  [open-problems.md](open-problems.md). The nearest are the 2-in-7 near miss at block 4 (open by 640 out of 137,280,
  a saving of 9.63e-2 if it exists), Proposition 6 for two arrays, the full cut of the lever on other tables (on #3945 with every head kept it is the
  headline, proved in Lean), and rules of the lever that are not separated.

## The sections

Each folder has a `README.md` (plain English, statuses in words, a Sources section at its end), a `results.json`
(every finding with its status, every open question with its first experiment, every program with its command),
`programs/` (grouped by finding, each folder with a short README: what it shows, the command, time and memory,
expected output) and `data/` (small certificates and outputs). Some also have `proofs/` or `notes/`.

<!-- the table below is written by make_index.py -->
| section | dates (2026) | what it is about | findings: proved in Lean / checked by computer / argued by hand / refuted / open | open questions |
|---|---|---|---|---|
| [moonshot-small-blocks](moonshot-small-blocks/README.md) | 14-15 Oct | the moonshot waves at blocks 2 and 3: the everything-free game and mixed units | 0 / 4 / 20 / 4 / 2 | 7 |
| [moonshot-block4-block5](moonshot-block4-block5/README.md) | 14-16 Oct | the moonshot waves at blocks 4 and 5: the 2-in-7 near miss and the shape searches | 0 / 21 / 18 / 14 / 11 | 10 |
| [reed-muller-and-one-move](reed-muller-and-one-move/README.md) | 15-16 Oct | the Reed-Muller family, and units that save a single move | 0 / 5 / 10 / 5 / 4 | 7 |
| [frame-packing-two-arrays](frame-packing-two-arrays/README.md) | 15 Oct | how savings happen: crossing slack, frame packing, the frame game, Proposition 6 | 0 / 12 / 24 / 13 / 1 | 13 |
| [why-the-savers-save](why-the-savers-save/README.md) | 13-15 Oct | why the record circuits save, and what caps them: anatomy, centre, copy wall, helpers | 0 / 19 / 10 / 9 / 0 | 13 |
| [block-bounds-and-new-rules](block-bounds-and-new-rules/README.md) | 13-16 Oct | lower bounds at each block size, in Lean; Toffoli and column-mixing maps | 5 / 18 / 31 / 6 / 10 | 10 |
| [theory-of-theories-and-logic](theory-of-theories-and-logic/README.md) | 14-15 Oct | what our methods can and cannot do: growth exponents, certificates, logic | 0 / 22 / 21 / 9 / 2 | 12 |
| [existence-routes](existence-routes/README.md) | 14-15 Oct | proving that a unit exists without finding it | 2 / 19 / 23 / 9 / 4 | 9 |
| [other-fields](other-fields/README.md) | 10-15 Oct | famous problems, network coding and quantum computing | 1 / 16 / 23 / 8 / 1 | 13 |
| [idea-rounds](idea-rounds/README.md) | 10 and 16 Oct | the two thinking rounds: closed doors, structural facts, open leads | 0 / 9 / 29 / 6 / 4 | 10 |
| [speed-waves](speed-waves/README.md) | 11-13 Oct | speed-first waves 11 to 15: real processors and FPGAs | 0 / 35 / 0 / 3 / 1 | 12 |
| [record-track](record-track/README.md) | 13-20 Oct | the Lean ladder, the levers that failed, the recursion lever and the frontier | 15 / 32 / 7 / 6 / 2 | 6 |
<!-- end of the table written by make_index.py -->

The counts are of the findings listed in each `results.json`; one finding can be a single number or a whole
theorem, so they say where the weight of a section lies, not how much it achieved. A reading order for a newcomer:
start with [why-the-savers-save](why-the-savers-save/README.md) (how the known savings work), then
[block-bounds-and-new-rules](block-bounds-and-new-rules/README.md) (what is impossible), then the two moonshot
sections (what was tried at the open shapes), then the others in any order.

[index.json](index.json) is the machine-readable index of this table (titles, dates, counts, paths), and
[open-problems.md](open-problems.md) lists every open question with its first experiment. Both are written by
[make_index.py](make_index.py) from the sections' `results.json` (`python3 -I -B research/make_index.py`, from the
root of the repository; `--check` compares them with the sections and changes nothing).

## The status words

Every finding carries one of five statuses, in these words, everywhere in this folder:

| status | what it means here |
|---|---|
| **proved in Lean** | a Lean 4 build accepted it using only the three standard axioms (`propext`, `Classical.choice`, `Quot.sound`) |
| **checked by computer** | an exact or exhaustive computation; the section says whether an AI skeptic repeated it with its own code, or whether one program, run once, is all there is |
| **argued by hand, not checked** | a written argument; the section says whether an AI skeptic re-derived it and found no error, or whether no skeptic read it. No person and no proof checker has checked it |
| **refuted** | shown false, with what refuted it: a skeptic's counterexample, a later theorem or a bug |
| **open** | nobody settled it either way; it comes with a first experiment |

A search that found nothing is evidence, not proof, and is marked open. A price of a construction that has not
been built is a price computed by a program, not a proved figure. No status is ever upgraded: when a later verdict
is weaker than an earlier one, the later one is used.

## Words used

- **Saving.** An algorithm that uses about N (log N)^(1 - s) operations has saving s. Larger is better. Savings are
  written in the 7.47e-4 style.
- **Array, block m, paid move, free step.** The programs work on many *arrays* of data. A unit works on m address
  bits ("directions") at a time; m is the *block* size. A *paid move* applies one direction's two-point transform to
  one array and costs one recursive call. A *free step* (relabelling addresses, a phase, adding arrays together) costs
  a constant amount per entry.
- **Unit, "W in M".** A program that transforms W arrays with M paid moves. Done one at a time they would cost W m
  moves; a unit saves ln(mW/M)/ln m. A bound of c moves per array at block m rules out every saving above
  1 - ln c / ln m.
- **Frame.** The coordinate system an array is written in, which grows as moves are made; additions between arrays
  are free only between arrays at compatible frames.
- **Certificate, cap, demand, hook.** A *certificate* is a quantity that every free step keeps. The counting rule,
  *Theorem R*, turns it into a lower bound: one paid move raises it by at most a *cap*, and each array needs a
  *demand*. The *hook* is the best family of certificates found.
- **Games**, which say what is free:
  - the *engine game* has the free steps of this repository's Lean engine: rational mixing of arrays, shifts,
    Clifford address maps and phases;
  - the *widened game* also has every address permutation (Toffoli gates included) and every phase pattern of
    powers of i;
  - the *everything-free game* also has *chain maps* (a row's bits mixed with its neighbour's, `row ^= row << 1`),
    the running XOR across digit columns, controlled-S phases and Toffoli gates;
  - the *chi game* is the engine game with the running XOR;
  - the *commutative slice* has only shifts, paid moves and mixing;
  - the *frame game* is the engine with its Pauli part dropped.
- **Digit columns f.** A unit must work for every number f of digit columns of the data.
- **The bar.** A saving of 0.1 is the standing target for any new route; since 15 October a saving of 0.01 counts as
  a stepping-stone tier.
- **Near miss.** A shape that every known bound allows but no search has filled.
- **The community.** The contributors of the public repository
  [CrocSwap/integer-mult-bounds](https://github.com/CrocSwap/integer-mult-bounds), whose circuits the records are
  built on.

## Running the programs

Python 3.9 or later; numpy for most programs and scipy for some; sympy, python-flint, Singular, msolve or a C
compiler where a program's README says so. Every program folder says what it needs, the command, its time and its
memory, and what it should print. The heavy jobs ran on Linux machines of Google's Colab service (8 processors,
about 50 GB of memory) and the light ones on an Apple M4 Mac mini; each README says which.

Two warnings apply to every section:

- **msolve 0.6.5 reads a product written `x*x` as `x`.** Write powers as `x^k` and cross-check every msolve verdict
  with Singular. Each section that used msolve says how its verdicts were cross-checked.
- **`ledger.py` certifies a word only when the number of digit columns is 1 more than a multiple of 4.** Every search
  whose last check is `ledger.py` must also check the other three residues of f (or check raw arrays at f = 1 to 6,
  or certify in the frame game) ([frame-packing-two-arrays](frame-packing-two-arrays/README.md)).

## Sources

- The sections were written from the agents' working notes of the rounds held between 10 and 19 October 2026, the
  verdicts of their skeptics, and a private catalogue of results that I keep; each section's own Sources section
  names the notes and rounds behind it. The working notes themselves, the agents' instructions and logs, and the
  catalogue are not published: they carry machine details and bookkeeping. Programs and data were copied from those
  notes with paths made relative and machine names removed, and each section's programs were re-run for this
  revision where its README says so.
- Paper text and preprint copies are not copied into this folder, and files of the community's repository only in
  the few cases that NOTICE (section 10.i) names; everything else is cited by title, arXiv identifier, or pull
  request, commit and file hash.
- The counts in the table of sections are those of the `results.json` files as they stand in this revision
  (`index.json`, written by `make_index.py`).
