# How to rebuild this repository and re-run the checks

The AI agents I directed ran everything below on one machine on 2026-10-09, from a fresh copy of exactly these
files with no build directory: Apple M4 (10 cores, 16 GB memory), macOS 26.6.2, the copy on an external SSD. What
was run and what was not is said in each section; the section "Limits" at the end collects the caveats. Linux
commands are given where they differ; **none of these runs was on Linux.**

## 1. What is checked

| theorem (namespace `OAI.PowerSaving.WHT`) | statement | challenge module | solution module | comparator config |
|---|---|---|---|---|
| `wht_main_block_B2Ke16x` | `∃ solve W, WHTProgram solve W ∧ WHTTimeBoundsAt (1 - 5399225/(10:ℝ)^10) W` | `Work.CarrierCheck.ChallengeB2Ke16x` | `Work.CarrierCheck.SolutionB2Ke16x` | `comparator/B2Ke16x.json` |
| `wht_main_rank_B2Ke16x` | the same with `1 - 3155781/(10:ℝ)^10` | `Work.CarrierCheck.ChallengeB2Ke16xR` | `Work.CarrierCheck.SolutionB2Ke16xR` | `comparator/B2Ke16xR.json` |

A challenge module states the theorem with `sorry` and imports only Mathlib: it is the text to read to see what
is claimed (331 lines; the cost model, `wht`, `WHTProgram`, `WHTTimeBoundsAt`). A solution module proves a
theorem of the same name. Two tools compare the two:

* **comparator** (github.com/leanprover/comparator): the solution proves the same statement as the challenge, uses
  only the permitted axioms `propext`, `Quot.sound`, `Classical.choice`, and is accepted by the Lean kernel.
* **`tools/Compare.lean`** (this repository, about 220 lines): an independent second check of "same statement
  as a kernel term, every definition the statement depends on identical, permitted axioms only".

## 2. What you need

* Linux or macOS, `git`, `curl`, and `elan` (the Lean version manager, github.com/leanprover/elan). The file
  `lean-toolchain` selects Lean v4.34.1; `lakefile.lean` pins Mathlib at `d13f23b723b8a846827a245b89c10fc7d3f11612`.
* Disk: about 9.5 GB. Measured: `.lake/packages` 7.6 GB (Mathlib and its dependencies with compiled files),
  `.lake/build` 1.3 GB (this repository), Mathlib's download cache 0.45 GB (in `~/.cache/mathlib` unless
  `MATHLIB_CACHE_DIR` is set; the agents set it to an empty directory).
* Memory: 16 GB was enough here when the ten modules that load the certificate are built one at a time (as
  below). Largest resident size of one process measured here: 5.8 GB for a Lean process during the build,
  6.2 GB for `tools/Compare.lean`, 6.3 GB during a comparator run. These figures include the mapped Mathlib
  files. Less than 16 GB was not tried.
* Python 3 (standard library only) for section 6; the agents used Python 3.9.6.
* For section 5: the comparator and its helper programs, see there.

## 3. Rebuild

    lake exe cache get

downloads the compiled Mathlib files for the pinned revision (Lake first clones Mathlib and its 8 dependencies).
Measured: 59 s, 8908 files, exit code 0. Then, on a 16 GB machine:

    LEAN_NUM_THREADS=4 lake build Work.CarrierCheck.B2 Work.CarrierCheck.RateB2Ke16 Work.CarrierCheck.RateB2Ke16x \
        Work.FoldRate.Group Work.FoldRate.SharpTools \
        Work.CarrierCheck.ChallengeB2Ke16x Work.CarrierCheck.ChallengeB2Ke16xR
    for m in Gen.Cr2h16 Gen.Cr2h16Hist Gen.Cr2h16Lab Gen.Cr2h16Scal Gen.Cr2h16Shape Cr2h16Cert \
             InstB2Ke16 B2Ke16 SolutionB2Ke16x SolutionB2Ke16xR; do
      lake build Work.CarrierCheck.$m || break
    done

The first command builds the 184 modules that do not load the certificate (89 from openai/math, 95 of this
project, both challenge modules among them); the loop builds the other 10 one at a time, so that two kernel
checks of the certificate never run together. Every `lake build` must end with
`Build completed successfully`, and no line may start with `error`. The build prints many linter warnings
(deprecated lemma names, unused `simp` arguments); they do not touch correctness. The warning "declaration uses
`sorry`" must appear for the two challenge modules and for no other module (a challenge states its theorem with
`sorry` on purpose). The last two builds print

    'OAI.PowerSaving.WHT.wht_main_block_B2Ke16x' depends on axioms: [propext, Classical.choice, Quot.sound]
    'OAI.PowerSaving.WHT.wht_main_rank_B2Ke16x' depends on axioms: [propext, Classical.choice, Quot.sound]

Measured here (all exit code 0, 194 modules, no error; "declaration uses `sorry`" for the two challenge modules only):

| part | time | largest resident size of one Lean process |
|---|---|---|
| the 184 modules (run as 16 `lake build` calls of 12 modules each, `LEAN_NUM_THREADS=4`) | 403 s | 5.8 GB |
| `Gen.Cr2h16` (the data, 3.9 MB of Lean text) | 30 s | 2.3 GB |
| `Gen.Cr2h16Hist` (43 block counts, `decide +kernel`) | 50 s | 0.9 GB |
| `Gen.Cr2h16Lab` (`labelCheck_ok`) | 17 s | 4.5 GB |
| `Gen.Cr2h16Scal` (`scalarCheck_ok`) | 38 s | 5.4 GB |
| `Gen.Cr2h16Shape` (`shapeCheck_ok`) | 9 s | 3.1 GB |
| `Cr2h16Cert` | 63 s | 5.7 GB |
| `InstB2Ke16`, `B2Ke16`, `SolutionB2Ke16x`, `SolutionB2Ke16xR` | 8 to 9 s each | 5.3 to 5.5 GB |
| **whole build** | **644 s (about 11 minutes)** | |

The agents ran the first command in 16 pieces (to share the machine); the single command above names the same 184
modules and was afterwards run as written, with the loop, on the finished build (all eleven calls ended with
`Build completed successfully` and built nothing). With more memory one command does
everything, `lake build Work.CarrierCheck.SolutionB2Ke16x Work.CarrierCheck.SolutionB2Ke16xR
Work.CarrierCheck.ChallengeB2Ke16x Work.CarrierCheck.ChallengeB2Ke16xR`, but it runs up to four kernel checks
of the certificate at the same time; it was not run in that form.

`lake build` without a module name builds nothing: `lakefile.lean` declares no default target. (With one, a
plain `lake build` would start several kernel checks of the certificate at the same time.) Run here on the
finished build, it printed `warning: no targets specified and no default targets configured` and
`Nothing to build.` (exit code 0). Always name the modules, as above.

## 4. Check with `tools/Compare.lean`

    lake env lean --run tools/Compare.lean Work.CarrierCheck.ChallengeB2Ke16x \
        Work.CarrierCheck.SolutionB2Ke16x OAI.PowerSaving.WHT.wht_main_block_B2Ke16x
    lake env lean --run tools/Compare.lean Work.CarrierCheck.ChallengeB2Ke16xR \
        Work.CarrierCheck.SolutionB2Ke16xR OAI.PowerSaving.WHT.wht_main_rank_B2Ke16x

Expected last lines (exit code 0); the second run says `33452 constants`:

    OK: every statement dependency is identical in the solution
    solution proof reaches 36438 constants; axioms: #[Classical.choice, Quot.sound, propext]
    OK: only permitted axioms

    RESULT: PASS

Measured here: both `RESULT: PASS`, exit code 0, 140 s and 117 s. Memory was not recorded in these two runs; a
repeat of the first one with a memory sampler gave `RESULT: PASS` again, 144 s, largest resident size 6.2 GB.

## 5. Check with the official comparator

The comparator needs three programs on `PATH`: `comparator` itself (github.com/leanprover/comparator),
`lean4export` (github.com/leanprover/lean4export) built for the Lean version of this repository, and the Linux
sandbox `landrun` (github.com/Zouuup/landrun). Its README says how to install them and under which assumptions
its verdict can be trusted. Used here: comparator at commit `d03acab154d269c06e60e4de7e4cc85deebff94b`,
lean4export at `076e8e57707e813375e8f9da8bf989799ace9680`, in both checkouts `lean-toolchain` changed from
v4.34.0 to v4.34.1 (the Lean version of this repository) and nothing else.

    lake env comparator comparator/B2Ke16x.json
    lake env comparator comparator/B2Ke16xR.json

On Linux the comparator's README recommends to run it as

    systemd-run --property=RestrictAddressFamilies=~AF_UNIX --user --pty -E PATH="$PATH" \
        --working-directory $(pwd) -- bash -c 'lake env comparator comparator/B2Ke16x.json'

Expected last lines of each run, exit code 0:

    Running Lean default kernel on solution.
    Lean default kernel accepts the solution
    Your solution is okay!

The comparator starts with `lake build` of the challenge and of the solution inside the sandbox. After section 3
nothing is left to build, so it only replays (0 modules built in the runs here) and then exports and re-checks the
proof in the kernel.

**macOS.** `landrun` does not exist for macOS. The agents ran the comparator with `tools/landrun-macos-shim.sh`
installed under the name `landrun` in a directory on `PATH` (it translates the sandbox request into a
`sandbox-exec` profile: no network, writes only where the comparator allows them, in the system's temporary
directories and in `/dev`). For example:

    mkdir -p "$HOME/landrun-shim" && cp tools/landrun-macos-shim.sh "$HOME/landrun-shim/landrun"
    chmod +x "$HOME/landrun-shim/landrun" && export PATH="$HOME/landrun-shim:$PATH"

Measured here, in the fresh copy after section 3 (macOS, the stand-in above, `enable_nanoda` false):

| config | verdict lines | exit | time | largest resident size |
|---|---|---|---|---|
| `comparator/B2Ke16x.json` | `Lean default kernel accepts the solution` / `Your solution is okay!` | 0 | 171 s | 6.3 GB |
| `comparator/B2Ke16xR.json` | `Lean default kernel accepts the solution` / `Your solution is okay!` | 0 | 178 s | 6.2 GB |

## 6. Regenerate the generated files (no Lean needed)

Fourteen Lean files and the two comparator configurations are generated from the certificate
`tools/certificate/c2-combine-best-h16.json` (ORIGIN.md lists them).

    python3 tools/regen_check.py

runs the generators into a temporary directory and compares the result with the files of this repository.
Expected last line, exit code 0, about 2 seconds:

    RESULT: REPRODUCED (16 of 16 repository files reproduced byte for byte)

The generators also write five files that are not in this repository (a comparator package with eight decimals,
`1 - 53992/10^8`, that the ten-decimal theorem makes superfluous, and two logs); the script lists them.

    python3 -I -B tools/ccheck.py tools/certificate/c2-combine-best-h16.json general
    python3 -I -B tools/ccheck.py tools/certificate/c2-combine-best-h16.json lpc

is a reference checker of the certificate written without Lean; expected first line `ACCEPTED (general): ...`
and `ACCEPTED (lpc): ...`, about 1 second each. Both were run here on this tree with these results.

None of the Python is part of the proof. The theorems rest on the Lean data module, which the Lean kernel checks
(`labelCheck_ok`, `shapeCheck_ok`, `scalarCheck_ok`, the block counts); the scripts only show where that module
comes from.

## 7. File integrity

    shasum -a 256 -c MANIFEST.sha256        # Linux: sha256sum -c MANIFEST.sha256

`MANIFEST.sha256` lists every file of the repository except itself.

## 8. Limits

Read these before relying on the result.

1. **The sandbox on macOS was a stand-in.** The comparator runs the build inside `landrun` (Linux Landlock).
   Here it ran inside `tools/landrun-macos-shim.sh`, which forbids network access and writes outside the allowed
   directories but, unlike `landrun`, does not restrict what can be read. The sandbox protects the checker
   against a solution file that attacks the build; it is not part of the logical check. None of the comparator
   runs here was on Linux.
2. **One kernel only.** Both configurations have `"enable_nanoda": false`: the proof was checked by Lean's own
   kernel (when the modules were compiled, and again by the comparator on the exported proof), not by a second,
   independent kernel. `leanchecker` was not run either. A re-run with `"enable_nanoda": true` is welcome.
3. **The comparator replayed files that were already compiled.** Section 3 compiles the solution modules before
   the comparator runs, outside its sandbox; in the comparator runs here no module was built. The comparator's
   README lists as an assumption that the solution was not compiled before. To meet it, run section 5 directly
   after `lake exe cache get`: the comparator then builds all 194 modules itself. It was not run in that order
   here; it compiles several kernel checks of the certificate at the same time and may need more than 16 GB.
4. **Mathlib was not rebuilt from source.** Its compiled files came from Mathlib's cache server for the pinned
   revision (`lake exe cache get`), as is usual.
5. **Not fetched fresh for this test:** the Lean toolchain v4.34.1 and the `comparator` and `lean4export`
   programs were already installed on the machine. Lake did clone Mathlib and its dependencies and download the
   compiled Mathlib files into an empty directory.
6. **One machine, one run.** Each step was run once, on the machine named at the top. No negative control (a
   challenge with a larger saving, which must be rejected) was repeated in this test. Times and memory are
   single measurements; "resident size" is what `ps` reports once a second and includes mapped files.
7. **Files added after the build.** The copy that was built held the Lean files, `lakefile.lean`,
   `lean-toolchain`, `lake-manifest.json`, `comparator/` and `tools/` of this repository. Documentation files
   (`README.md`, `RELATED-WORK.md`, `NOTICE`, `LICENSE`, `ORIGIN.md`, this file, `.gitignore`) were added or
   edited afterwards; Lean does not read them. After the last edit every `.lean` file, `lean-toolchain`,
   `lake-manifest.json` and both files in `comparator/` were compared by sha256 with the files of the copy that
   was built and checked: all identical.
8. **What the checks do not say.** They say that the solution proves the statement written in the challenge
   module. Whether that statement means what a reader expects has to be judged by reading the challenge module:
   the cost model in it is OpenAI's, the definitions `wht`, `WHTProgram`, `WHTTimeBoundsAt` are my project's.
   The checks say nothing about the speed of any program on inputs of realistic size.
9. **`tools/Compare.lean` is my project's own tool.** It is a second opinion written independently of the
   comparator, not an independent party. The comparator is the outside tool.
