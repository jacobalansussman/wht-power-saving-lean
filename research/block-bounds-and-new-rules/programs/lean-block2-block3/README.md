# The block-2 and block-3 lower bounds in Lean

**What it shows.** In the engine game of this repository, every unit at block 2 needs at least 2W paid moves
(`block2_bound`), and every unit at block 3 needs at least 14W/5 (`block3_bound`), with no hypothesis left. Along
the way it proves Theorem R for the engine's own words (`engine_lower_bound` in `TheoremR.lean`), that every free
letter keeps the hook subspace, and the hook's cap at every block size (`Hook.lean`).

**Files.**
- `PfLean/TheoremR.lean` (396 lines), `Hook.lean` (449), `Dem2.lean` (156), `Bounds.lean` (73): Theorem R, the hook,
  the block-2 demand and the two bounds; block 3 is proved here given the demand fact `224 ≤ raise ...`.
- `PfLean/Dem3.lean` (432 lines): the demand fact, from eight 28 × 28 character blocks (the split is proved).
- `PfLean/Dem3Defs.lean` (50), `Dem3Data.lean` (17, generated certificate data: rows, columns and inverses mod
  65537; untrusted), `Dem3Check.lean` (24: the eight checks, each `decide +kernel`), `Dem3Final.lean` (30: `dem3` and
  `block3_bound`, with `#print axioms`).
- `PfLean/Statements.lean`: prints the exact statements and their axioms; `logs/statements.txt` is its output.
- `PfLean/NegControl.lean`: one inverse entry changed by 1; Lean must prove that the check fails.
- `certificate-generator/blocks.py` (numpy): ranks of the 8 character blocks modulo 65537, the chosen rows and
  columns, and the inverses; `emit_lean.py` writes `Dem3Data.lean` from its JSON output
  (`../../data/lean-block3-certificate/cert65537.json`). Neither is trusted by the proof.
- `skeptic-restatement/Check.lean`, `NonVac3.lean`: a skeptic's own statement of the theorem, which Lean accepted with
  the same proof, and a one-array program that meets the hypothesis (so the theorem is not vacuous).

**Build.** The files import `Work.GFrame.Engine.Letters` (the engine of this repository) and Mathlib. In a checkout of
this repository (Lean 4.34.1, the Mathlib revision of `lakefile.lean`):
```
mkdir -p Work/PfLean && cp PfLean/*.lean Work/PfLean/
lake build Work.PfLean.Dem3Final Work.PfLean.Statements Work.PfLean.NegControl
```
(`skeptic-restatement/*.lean` import `Work.PfLean.Dem3Final`; copy them next to it to build them.) The eight kernel
checks take 10 to 45 minutes of one core, depending on the machine's load, and up to 10 GB of memory; everything else
builds in minutes once Mathlib is compiled. Expected: `dem3` and `block3_bound` depend on axioms `[propext,
Classical.choice, Quot.sound]`, the block checks on `[propext]`, and `neg_control` is proved (the corrupted
certificate is rejected).

**Regenerate the certificate:** `python3 -I certificate-generator/blocks.py 65537 cert.json` (numpy, 20 seconds), then
`python3 -I certificate-generator/emit_lean.py cert.json Dem3Data.lean`. Re-run on 19 October 2026: 8 blocks of rank
28 modulo 65537 (total 224), and the output is identical to the published certificate
(`../../data/lean-block3-certificate/rerun_2026-10-19.txt`).

**Where it was built.** On 13 and 14 October 2026 the files were built on two Linux machines of Google's Colab service
in copies of Ryan Shea's Lean tree (shea256/fourier-transform-below-nlogn, `lean/complex588`, commit 5239484), which
pins the same Lean and Mathlib versions as this repository and carries byte-identical copies of every engine file of
this repository that it contains; a skeptic rebuilt everything in a fresh copy and replayed it with Lean's independent
checker.

**Re-run for this revision:** yes. On 19 October 2026 these exact files (same SHA-256) were built again in a fresh
copy of the same tree on a Linux machine of Google's Colab service: exit code 0 in 11 minutes (the eight kernel checks
562 seconds), peak memory 10.06 GB, `dem3`, `block3_bound` and `block2_bound` on the three standard axioms only,
`neg_control` proved, no `sorry` and no `native_decide` (`logs/rebuild_2026-10-19.txt`). They have not yet been built
inside a checkout of this repository itself; they should build there unchanged, because the tree used pins the same
Lean and Mathlib and its copies of the engine files are byte-identical to this repository's.

**What it does not cover.** The widened game (the Lean engine's relabellings are affine only), the refined bound for
grouped moves, and block 4.
