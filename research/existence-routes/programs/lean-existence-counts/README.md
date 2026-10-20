# Lean: a unit known only to exist counts

**What it shows.** `ExLegit.lean` proves five theorems: the engine's kernel program from a unit given only by "there exists"
(`exlegit_engine_of_exists`), the same with the unit named by `Classical.choose`, the price of a mixing letter does not
depend on its matrix (`exlegit_costR_localM`, by `rfl`), a fixed shape with unknown coefficients whose price is proved
once for the shape (`exlegit_engine_of_coeff_exists`), and the same from the engine's existential certificate `GCert`.
`ExLegitClosure.lean` is OpenAI's `ComparatorChallenges/UniformFourier.lean`, unchanged, followed by
`ExLegitClosure.tail` (252 lines): every value of paint `scalar` in OpenAI's RAM lies in any subfield containing the
scalar inputs (`run_good`), so a Fourier program computes in Q(ζ_d), and `zeta_mem_of_DFTProgram`: every DFT program
has ζ_n in Q(ζ_d). OpenAI's file is not copied here; `check.sh` takes it from `third-party/openai-math/` and checks
its sha256.

**Commands.** From the repository root, after `lake build Work.GFrame.Engine.Main Work.GFrame.Engine.Cert`:

```
sh research/existence-routes/programs/lean-existence-counts/check.sh
```

or `sh check.sh <lean project> <path of UniformFourier.lean>` for another tree. It writes only into a temporary folder.

**Needs.** The repository's Lean toolchain and Mathlib build. Never on a small machine: each file needs about 6.6 GB.

**Time and memory.** 13 s and 2 minutes on an idle machine (14 October); 8 minutes on a loaded one (19 October).

**Expected output.** `logs/check-2026-10-19.log`: no errors; two `declaration uses sorry` warnings at lines 326 and 330 of
`ExLegitClosure.lean`, which are OpenAI's open challenge theorems and are used by nothing here; every `#print axioms`
line shows only `[propext, Classical.choice, Quot.sound]`. `logs/original-2026-10-14-*.log` are the first builds.

**Re-run for this revision.** On 19 October 2026 `check.sh` ran on a Linux machine of Google's Colab service (8 processors, about 50 GB), in Ryan Shea's Lean tree, whose `Work/GFrame/Engine/{Main,Letters,Cert}.lean`, `Work/Block/Word.lean` and `third-party/openai-math/lean/ComparatorChallenges/UniformFourier.lean` have the same sha256 as this repository's: PASS (`logs/check-2026-10-19.log`).

Every program here is the program as it ran in October 2026, with private paths removed; working names inside the files are mapped in `../../proofs/README.md`.

**Sources.** Shea's tree was at commit 5239484 for the re-check of 19 October 2026; the sha256 of OpenAI's `UniformFourier.lean` is 4c571f275bd7506c3ccd1e6d78f2fca77e07112c969c1c62ad8359d61544ec74 (also in `check.sh`).
