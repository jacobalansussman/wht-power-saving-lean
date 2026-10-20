# How to rebuild this repository and re-run the checks

This file says how to rebuild the proofs and re-run the checks yourself. The lever and the sixth result,
built on Ryan Shea's Lean tree with the packages under `package/`, are checked in section 12, and section 12.6
re-runs the lever's searches and prices on a laptop, with no Lean. The results of 13 October 2026, also built on
Shea's tree, are checked in section 11. Section 13 says how to run the programs of `research/`. To check the
fourth result, start at section 10. It sends you to sections 2 to 5 for what you need and for the first
commands.

The AI agents I directed ran everything in sections 1 to 9 on one machine, on 2026-10-09: Apple M4 (10 cores,
16 GB memory), macOS 26.6.2, the copy on an external SSD. What was run and what was not is said in each
section; the section "Limits" collects the caveats. Linux commands are given where they differ. **None of the
runs of 2026-10-09 was on Linux. Two things were done on 2026-10-10 (section 10). The published third
revision was cloned on a Linux machine, and all its 342 Lean modules were compiled there from source; Mathlib's
compiled files came from its cache server. The 39 new modules of the fourth result were compiled from source
on that machine and on the Mac, each time in a tree that already held the compiled modules of the earlier
results. A build of the whole fourth revision from an empty directory, with nothing compiled before, has not
been made (10.7, point 4).** Section 9 keeps the measurements of the first publication (2026-10-09, 194
modules).

**State of the checks of the fourth result (2026-10-10).** The three theorems of the fourth result are
`wht_main_block_B2Ge8x`, `transform_mainE8` and `convolution_mainE8`. On the Mac and on the Linux machine:

* the 39 new modules were compiled from source, in a tree that already held the compiled modules of the
  earlier results;
* `tools/Compare.lean` passed for the Walsh-Hadamard theorem and for the Fourier pair;
* the official comparator passed for both configurations. Each time it ran with a stand-in for its sandbox,
  on modules that had been compiled before.

Three facts about how far these builds go:

* A build of the whole fourth revision from an empty directory has not been made. Such a build starts with
  nothing compiled and compiles every module of the revision in one run.
* On the Mac the 39 modules were compiled a second time. For that run the files were exported from git and
  laid over a copy of a complete build of the third revision. That tree held none of the working files of the
  day. All 39 modules compiled there from the exported files.
* The revision has 391 files that Lean reads: the Lean sources, the comparator configurations,
  `lakefile.lean`, `lake-manifest.json` and `lean-toolchain`. By sha256, all 391 are identical to the files
  in the trees where the builds and checks were run, on the Mac and on Linux.

10.7, point 4, has the details.

The comparator's real sandbox ran on neither machine. Section 10 has the commands, the expected lines, what
was measured, and the limits that are particular to the fourth result (10.7).

**State of the checks of the second result when this text was written.** The theorem passed the official
comparator three times (section 5) and `tools/Compare.lean` twice (section 4): in the project's working tree, on
Lean files that are byte for byte those of this repository, and then in a copy of the repository built with the
`lakefile.lean` of this repository, partly from compiled files of the working build (section 3). On that day,
2026-10-09, no build of this revision from nothing had been made (section 8, point 12). The third revision,
which contains this result, was built from nothing on Linux on 2026-10-10 (section 8, point 14, and 10.6).

**State of the checks of the third result.** Both theorems passed the official comparator three times
(section 5): twice in the project's working tree and once in a copy of the repository, on Lean files that are
byte for byte those of this repository. `tools/Compare.lean` passed in the working tree and in the copy. The
times of these runs are in the README, section 9. On that day, 2026-10-09, no build of this revision from
nothing had been made (section 8, point 13). One was made on Linux on 2026-10-10 (section 8, point 14, and 10.6).

## 1. What is checked

| theorem (namespace `OAI.PowerSaving.WHT` unless written out) | statement | challenge module | solution module | comparator config |
|---|---|---|---|---|
| `OAI.PowerSaving.transform_mainE8` and `OAI.PowerSaving.convolution_mainE8` (fourth result) | `DFTGoalE8` and `ConvGoalE8`: OpenAI's `DFTGoal` and `ConvGoal` with `decimalExponentE8 := 1 - 8762478/(10^(10:ℕ))` | `Work.FourierE8.UniformFourierChallenge` | `Work.FourierE8.Main` | `comparator/UniformFourierE8.json` |
| `wht_main_block_B2Ge8x` (fourth result) | `∃ solve W, WHTProgram solve W ∧ WHTTimeBoundsAt (1 - 8762479/(10:ℝ)^10) W` | `Work.GCert.Data.ChallengeB2Ge8x` | `Work.GCert.Data.SolutionB2Ge8x` | `comparator/B2Ge8x.json` |
| `OAI.PowerSaving.transform_mainZ` and `OAI.PowerSaving.convolution_mainZ` (third result) | `DFTGoalZ` and `ConvGoalZ`: OpenAI's `DFTGoal` and `ConvGoal` with `decimalExponentZ := 1 - 7474546/(10^(10:ℕ))` | `Work.Fourier.UniformFourierChallenge` | `Work.Fourier.Main` | `comparator/UniformFourier.json` |
| `wht_main_block_B2Gp193x` (second result) | `∃ solve W, WHTProgram solve W ∧ WHTTimeBoundsAt (1 - 7474547/(10:ℝ)^10) W` | `Work.GCert.Data.ChallengeB2Gp193x` | `Work.GCert.Data.SolutionB2Gp193x` | `comparator/B2Gp193x.json` |
| `wht_main_block_B2Ke16x` (first result) | the same with `1 - 5399225/(10:ℝ)^10` | `Work.CarrierCheck.ChallengeB2Ke16x` | `Work.CarrierCheck.SolutionB2Ke16x` | `comparator/B2Ke16x.json` |
| `wht_main_rank_B2Ke16x` (first result, per-rank) | the same with `1 - 3155781/(10:ℝ)^10` | `Work.CarrierCheck.ChallengeB2Ke16xR` | `Work.CarrierCheck.SolutionB2Ke16xR` | `comparator/B2Ke16xR.json` |

A challenge module states the theorem with `sorry` and imports only Mathlib: it is the text to read to see what
is claimed (331 lines; the cost model, `wht`, `WHTProgram`, `WHTTimeBoundsAt`). The four
Walsh-Hadamard challenge modules differ in the exponent, in the name of the theorem and in the comments that quote the two
and the module names, and in nothing else. A solution module proves a theorem of the same name. Two tools compare the two:

* **comparator** (github.com/leanprover/comparator): the solution proves the same statement as the challenge, uses
  only the permitted axioms `propext`, `Quot.sound`, `Classical.choice`, and is accepted by the Lean kernel.
* **`tools/Compare.lean`** (this repository, about 220 lines): an independent second check of "same statement
  as a kernel term, every definition the statement depends on identical, permitted axioms only".

The challenge module of the third result (348 lines) is OpenAI's challenge file with three
changes, which its first comment lists: that comment, the exponent, and the suffix `Z` on seven names. To see
them (the first file is OpenAI's, unchanged; ORIGIN.md, last section, says how to compare it with a checkout
of openai/math):

    diff third-party/openai-math/lean/ComparatorChallenges/UniformFourier.lean Work/Fourier/UniformFourierChallenge.lean

Expected: 9 hunks. One block of 12 added lines after line 2 (an empty line and the comment), and 8 changed
lines, OpenAI's lines 269, 274, 317, 319, 321, 323, 326 and 330: the exponent, and the seven names where they
are defined or used. So 20 lines are marked `>` and 8 are marked `<`, and nothing else differs.

## 2. What you need

* Linux or macOS, `git`, `curl`, and `elan` (the Lean version manager, github.com/leanprover/elan). The file
  `lean-toolchain` selects Lean v4.34.1; `lakefile.lean` pins Mathlib at `d13f23b723b8a846827a245b89c10fc7d3f11612`.
* Disk: about 10.2 GB. Measured: `.lake/packages` 7.6 GB (Mathlib and its dependencies with compiled files),
  `.lake/build` 2.0 GB (this repository), Mathlib's download cache 0.45 GB (in `~/.cache/mathlib` unless
  `MATHLIB_CACHE_DIR` is set; the agents set it to an empty directory).
* Memory: **plan for 16 GB free for one run of `tools/Compare.lean` on Linux, and run nothing beside it on a
  16 GB Mac.** Corrected on 2026-10-10: do not read the Mac record of 2026-10-09 below as "16 GB is enough".
  What was measured:

  - **Mac, 2026-10-09.** 16 GB was enough when the modules that load a certificate were built one at a time
    (as below) and nothing else ran. Largest resident size of one process: 7.0 GB (in the project's working
    tree) for a Lean process during the build, 7.2 GB for `tools/Compare.lean`, 9.6 GB during a comparator
    run. These figures include the mapped Mathlib files. Less than 16 GB was not tried.
  - **Linux machine (51 GiB, no swap), 2026-10-10.** The largest process of the build was 10.5 GB
    (`Gen.P193Scal0`, `Gen.P193Scal1`), against 7.0 GB on the Mac. One run of `tools/Compare.lean` was
    reported at 15.7 to 16.0 GB by `/usr/bin/time -v`. That is more than twice the Mac's figure. The likely
    reason is the way Linux counts the two mappings of the Mathlib files, challenge side and solution side.
    That reason is a guess and was not tested.
  - **Mac, 2026-10-10.** Two audit runs of `tools/Compare.lean` were stopped by the auditing agents' own
    memory rule, at 9.3 and 9.7 GB of footprint. Footprint is a process's own memory, without shared library
    files. Both runs were on an earlier unit of the fourth result's kind (a unit is the finite network that
    the saving is computed from).
* Python 3 (standard library only) for section 6; the agents used Python 3.9.6. The rebuild of the second
  certificate needs Python 3.9 at least and network access for nine small downloads.
* For section 5: the comparator and its helper programs, see there.

## 3. Rebuild

    lake exe cache get

downloads the compiled Mathlib files for the pinned revision (Lake first clones Mathlib and its 8 dependencies).
Measured at the first publication: 59 s, 8908 files, exit code 0. Then, on a 16 GB machine, first the modules
that load no certificate:

    LEAN_NUM_THREADS=4 lake build Work.CarrierCheck.B2 Work.CarrierCheck.RateB2Ke16 Work.CarrierCheck.RateB2Ke16x \
        Work.FoldRate.SharpTools Work.CarrierCheck.ChallengeB2Ke16x Work.CarrierCheck.ChallengeB2Ke16xR \
        Work.GFrame.Top.All Work.GCert.Chain.End Work.GCert.Scalar.Main Work.GCert.Data.Glue \
        Work.GCert.Data.Price Work.GCert.Data.InstB2Gp193 Work.GCert.Data.ChallengeB2Gp193x

and then the modules that load a certificate, one at a time, so that two kernel checks of a certificate never
run together:

    for m in CarrierCheck.Gen.Cr2h16 CarrierCheck.Gen.Cr2h16Hist CarrierCheck.Gen.Cr2h16Lab \
             CarrierCheck.Gen.Cr2h16Scal CarrierCheck.Gen.Cr2h16Shape CarrierCheck.Cr2h16Cert \
             CarrierCheck.InstB2Ke16 CarrierCheck.B2Ke16 CarrierCheck.SolutionB2Ke16x \
             CarrierCheck.SolutionB2Ke16xR GFrame.Top.Regress \
             GCert.Data.Gen.P193.S0 GCert.Data.Gen.P193.S1 GCert.Data.Gen.P193.S2 GCert.Data.Gen.P193.S3 \
             GCert.Data.Gen.P193.S4 GCert.Data.Gen.P193.S5 GCert.Data.Gen.P193.S6 GCert.Data.Gen.P193.S7 \
             GCert.Data.Gen.P193 GCert.Data.Gen.P193Par \
             GCert.Data.Gen.P193Tab0 GCert.Data.Gen.P193Tab1 GCert.Data.Gen.P193Tab2 GCert.Data.Gen.P193Tab3 \
             GCert.Data.Gen.P193Seg0 GCert.Data.Gen.P193Seg1 GCert.Data.Gen.P193Seg2 GCert.Data.Gen.P193Seg3 \
             GCert.Data.Gen.P193Seg4 GCert.Data.Gen.P193SegF GCert.Data.Gen.P193End GCert.Data.Gen.P193Hist \
             GCert.Data.Gen.P193YChk GCert.Data.Gen.P193Scal0 GCert.Data.Gen.P193Scal1 \
             GCert.Data.P193Price GCert.Data.P193Cert GCert.Data.P193ScalOK \
             GCert.Data.B2Gp193 GCert.Data.B2Gp193Main GCert.Data.SolutionB2Gp193x; do
      lake build Work.$m || break
    done

The first command builds 285 modules (89 from openai/math, the others of this project, the three
challenge modules among them); the loop builds the other 42. Every `lake build` must end with
`Build completed successfully`, and no line may start with `error`. The build prints many linter warnings
(deprecated lemma names, unused `simp` arguments); they do not touch correctness. The warning "declaration uses
`sorry`" must appear for the three challenge modules and for no other module (a challenge states its theorem
with `sorry` on purpose). The builds of the three solution modules print

    'OAI.PowerSaving.WHT.wht_main_block_B2Gp193x' depends on axioms: [propext, Classical.choice, Quot.sound]
    'OAI.PowerSaving.WHT.wht_main_block_B2Ke16x' depends on axioms: [propext, Classical.choice, Quot.sound]
    'OAI.PowerSaving.WHT.wht_main_rank_B2Ke16x' depends on axioms: [propext, Classical.choice, Quot.sound]

Measured (327 modules; "declaration uses `sorry`" for the three challenge modules only). The commands
above were not run in this form for this revision: each row says where and how its part was built.

| part | time | largest resident size of one Lean process |
|---|---|---|
| the 285 modules that load no certificate | 184 of them at the first publication: 403 s (section 9). The 101 new ones, compiled from source in the copy: 470 s for 90 of them (nine `lake build` calls of ten modules), 8 to 169 s per call for the other 11 in separate calls; a later change of one source file recompiled it and what imports it (71 s) | 5.8 GB (section 9); 5.7 GB for the new ones |
| the 10 certificate modules of the first result | section 9 (first publication; not rebuilt for this revision) | section 9 |
| `Work.GFrame.Top.Regress` (the first result once more, through the generalised engine; loads the first certificate), compiled from source in the copy | 220 s | not sampled |
| `Gen.P193.S0` to `S7`, `Gen.P193`, `Gen.P193Par` (the data of the second certificate, 7.5 MB of Lean text). This row and the next five: measured in the project's working tree, one module per `lake build`; in the copy the compiled files were cloned (below) | 72 s for the eight data modules | 1.9 GB |
| `Gen.P193Tab0` to `Tab3` (frame table, 4 evaluations) | 39 to 54 s each | 0.74 GB |
| `Gen.P193Seg0` to `Seg4` (replay of the gates, 5 evaluations) | 6, 17, 55, 156 and 186 s | up to 3.8 GB |
| `Gen.P193SegF`, `Gen.P193End`, `Gen.P193YChk` | 10 s, 9 s, 2 s | up to 1.7 GB |
| `Gen.P193Hist` (63 block counts) | 54 s | 1.2 GB |
| `Gen.P193Scal0`, `Gen.P193Scal1` (scalar replay, 660 sources each) | 122 to 180 s each | 7.0 GB each |
| `P193Price`, `P193Cert`, `P193ScalOK`, `B2Gp193`, `B2Gp193Main`, `SolutionB2Gp193x` | working tree: 5, 5, 24, 8, 42 and 10 s. In the copy, from source, in one `lake build` with `LEAN_NUM_THREADS=3`: 276 s | not sampled in the copy |
| **whole build** | **not run as one build for this revision** | |

**How the copy of this revision was built.** It was not built from nothing. It started from the compiled files
of the first publication's rebuild from scratch (section 9: Mathlib from its cache, the 194 published modules
compiled from source), cloned into a new copy that holds the files of this revision. Then, in the copy, the 101
new modules that load no certificate and `Work.GFrame.Top.Regress` were compiled from source (table). For the 31
modules that load the certificate of the second result the fast path was used: their compiled files were cloned
from the project's working build (same source bytes), and a `lake build` of the eight roots accepted the 25
modules under `Work/GCert/Data/Gen/` by their traces ("Replayed") and compiled the 6 modules above them from
source. Finally `lake build --no-build` of the eight roots (the four published ones, `Work.GFrame.Top.All`,
`Work.GFrame.Top.Regress`, the new challenge and the new solution) reported "All targets up-to-date (9250 jobs)"
at 15:48, exit code 0, no line starting with `error`, with the axiom lines above. The times of the certificate
modules in the table are therefore those of the working tree, where the agent that built the proof compiled each
of them from source. On that day, 2026-10-09, no build of the whole revision from nothing, with the commands
above, had been made; the one of 2026-10-10, on Linux, is in 10.6. The
runs of sections 4 and 5 in this copy came after these steps.

The certificate of the second result is checked by the kernel in segments (`Work/GCert/Data/Gen/`, one
module for each large evaluation):
13 evaluations for the label half (4 for the frame table, 5 for the replay of the gates and 1 for its final
segment, 3 for the end conditions and the frames of the copies; two of the last three are in `P193Cert`), 3 for
the scalar half (2 for the replay, 660 sources each, and 1 for the outputs), and 63 block counts in
`Gen.P193Hist`.
One of these modules needs up to 7.0 GB and about 3 minutes (table). Build them one at a time.

`lake build` without a module name builds nothing: `lakefile.lean` declares no default target. (With one, a
plain `lake build` would start several kernel checks of a certificate at the same time.) Always name the
modules, as above.

**Third result.** The modules of `Work/Fourier` hold no certificate data and make no kernel evaluation of a
certificate. They import the final module of the second result, so they are built after the loop above:

    lake build Work.Fourier.UniformFourierChallenge Work.Fourier.Main Work.Fourier.Axioms

The warning "declaration uses `sorry`" then appears for a fourth module, the challenge
`Work.Fourier.UniformFourierChallenge`. The build of `Work.Fourier.Axioms` prints, among other lines:

    'OAI.PowerSaving.transform_mainZ' depends on axioms: [propext, Classical.choice, Quot.sound]
    'OAI.PowerSaving.convolution_mainZ' depends on axioms: [propext, Classical.choice, Quot.sound]
    'OAI.PowerSaving.RAM.hillsZ_program' depends on axioms: [propext, Classical.choice, Quot.sound]

Measured in the project's working tree, on the compiled modules of the second result: `Work.Fourier.Seam` 35 s;
the other 14 modules in one `lake build`, 58 s (9238 jobs, exit code 0). Memory was not sampled.

## 4. Check with `tools/Compare.lean`

    lake env lean --run tools/Compare.lean Work.GCert.Data.ChallengeB2Gp193x \
        Work.GCert.Data.SolutionB2Gp193x OAI.PowerSaving.WHT.wht_main_block_B2Gp193x
    lake env lean --run tools/Compare.lean Work.CarrierCheck.ChallengeB2Ke16x \
        Work.CarrierCheck.SolutionB2Ke16x OAI.PowerSaving.WHT.wht_main_block_B2Ke16x
    lake env lean --run tools/Compare.lean Work.CarrierCheck.ChallengeB2Ke16xR \
        Work.CarrierCheck.SolutionB2Ke16xR OAI.PowerSaving.WHT.wht_main_rank_B2Ke16x

Expected last lines (exit code 0); the first run says `38678 constants`, the second `36438 constants`,
the third `33452 constants`:

    OK: every statement dependency is identical in the solution
    solution proof reaches 38678 constants; axioms: #[Classical.choice, Quot.sound, propext]
    OK: only permitted axioms

    RESULT: PASS

Measured: for the second result `RESULT: PASS`, exit code 0, twice: 227 s, 4.3 GB by an auditing agent in the
project's working tree (15:43 to 15:47), and 148 s, 7.2 GB in the build copy (16:13 to 16:15). The two runs for
the first result were not repeated for this revision (section 9).

**Third result.** `tools/Compare.lean` takes both theorems in one run:

    lake env lean --run tools/Compare.lean Work.Fourier.UniformFourierChallenge \
        Work.Fourier.Main OAI.PowerSaving.transform_mainZ OAI.PowerSaving.convolution_mainZ

It passed for the third result in the project's working tree and, at 19:23 UTC, in a copy of the repository (163 s).

## 5. Check with the official comparator

The comparator needs three programs on `PATH`: `comparator` itself (github.com/leanprover/comparator),
`lean4export` (github.com/leanprover/lean4export) built for the Lean version of this repository, and the Linux
sandbox `landrun` (github.com/Zouuup/landrun). Its README says how to install them and under which assumptions
its verdict can be trusted. Used here: comparator at commit `d03acab154d269c06e60e4de7e4cc85deebff94b`,
lean4export at `076e8e57707e813375e8f9da8bf989799ace9680`, in both checkouts `lean-toolchain` changed from
v4.34.0 to v4.34.1 (the Lean version of this repository) and nothing else.

    lake env comparator comparator/UniformFourier.json
    lake env comparator comparator/B2Gp193x.json
    lake env comparator comparator/B2Ke16x.json
    lake env comparator comparator/B2Ke16xR.json

On Linux the comparator's README recommends to run it as

    systemd-run --property=RestrictAddressFamilies=~AF_UNIX --user --pty -E PATH="$PATH" \
        --working-directory $(pwd) -- bash -c 'lake env comparator comparator/B2Gp193x.json'

Expected last lines of each run, exit code 0:

    Running Lean default kernel on solution.
    Lean default kernel accepts the solution
    Your solution is okay!

The comparator starts with `lake build` of the challenge and of the solution inside the sandbox. After section 3
nothing is left to build, so it only replays and then exports and re-checks the proof in the kernel.

**macOS.** `landrun` does not exist for macOS. The agents ran the comparator with `tools/landrun-macos-shim.sh`
installed under the name `landrun` in a directory on `PATH` (it translates the sandbox request into a
`sandbox-exec` profile: no network, writes only where the comparator allows them, in the system's temporary
directories and in `/dev`). For example:

    mkdir -p "$HOME/landrun-shim" && cp tools/landrun-macos-shim.sh "$HOME/landrun-shim/landrun"
    chmod +x "$HOME/landrun-shim/landrun" && export PATH="$HOME/landrun-shim:$PATH"

Measured (macOS, the stand-in above, `enable_nanoda` false; each row says where the run was made):

| config | verdict lines | exit | time | largest resident size |
|---|---|---|---|---|
| `comparator/UniformFourier.json` (third result, both theorems), in the working tree, by the agent that built the proof, ended 18:27 | "Lean default kernel accepts the solution" / "Your solution is okay!" | 0 | 1528 s | 9.0 GB (9,660,301,312 bytes, the largest process as `time -l` reports it) |
| `comparator/B2Gp193x.json`, in the working tree, by the agent that built the proof, ended 15:14 | "Lean default kernel accepts the solution" / "Your solution is okay!" | 0 | 1814 s | 9.1 GB |
| the same, in the working tree, by an auditing agent with its own copy of the configuration, 15:25 to 15:43 | the same two lines | 0 | 1116 s | 9.4 GB (all Lean processes together) |
| the same, in the build copy of this repository, 15:55 to 16:13 (every module replayed, none built) | the same two lines | 0 | 1058 s | 9.6 GB |
| `comparator/B2Ke16x.json`, `comparator/B2Ke16xR.json` | not re-run for this revision; section 9 | | | |

## 6. Regenerate the generated files (no Lean needed)

    python3 tools/regen_check.py
    python3 tools/gx/regen_check.py

run the generators into a temporary directory and compare the result with the files of this repository. The
first covers the first result (certificate `tools/certificate/c2-combine-best-h16.json`), the second the second
result (certificate `tools/certificate/gcert1-p11-pr193.json.gz`). Expected last lines, exit code 0:

    RESULT: REPRODUCED (16 of 16 repository files reproduced byte for byte)
    RESULT: REPRODUCED (36 repository files reproduced byte for byte, 0 different, 3 generated files not in the repository)

Both lines were obtained on the files of this revision (15:22, and again at 16:11). The three generated files
that are not in the repository are the eight-decimal variant of the second result (`ChallengeB2Gp193`,
`SolutionB2Gp193` and its comparator configuration). It compiles in the working tree and has not been put
through the comparator, so it is not published.

For the third result, the copies of OpenAI's files under `Work/Fourier/`:

    python3 tools/fourier/regen_check.py

makes the 11 copies, `Goal.lean`, the challenge and the comparator configuration again, in a temporary
directory, from the unmodified files under `OAI/` and `third-party/openai-math/`, and compares them with the
files of this repository. It needs Python 3 and no Lean. Expected last line, exit code 0:

    RESULT: REPRODUCED (14 repository files reproduced byte for byte, 0 different, 0 generated files not in the repository)

Reference checkers of the two certificates, written without Lean:

    python3 -I -B tools/ccheck.py tools/certificate/c2-combine-best-h16.json general
    python3 -I -B tools/ccheck.py tools/certificate/c2-combine-best-h16.json lpc
    python3 -B tools/gx/refcheck.py tools/certificate/gcert1-p11-pr193.json.gz

Expected first lines `ACCEPTED (general): ...`, `ACCEPTED (lpc): ...` and `ACCEPTED by gx.check1: h=22 v=1320 R=9412 N=262944 gates A=12158 B=30878; ...`
(17 to 19 s here; this script needs its own directory on the module path, so it is run without `-I`, and
it creates the empty directory `tools/gx/out/`, which git ignores).

    python3 tools/whatif_copies.py tools/certificate/gcert1-p11-pr193.json.gz tools/certificate/c2-combine-best-h16.json

recomputes the what-if of the README, section 7. It is a what-if and no part of any proof.

None of the Python is part of the proof. The theorems rest on the Lean data modules, which the Lean kernel
checks; the scripts only show where those modules come from. **Where the chain starts:** for the second result
the Lean files start from the rebuilt certificate `tools/certificate/gcert1-p11-pr193.json.gz`.
The programs that rebuilt it from the community's files are in `tools/rebuild/`. The community's files are
those of the outside repository github.com/CrocSwap/integer-mult-bounds. They are not in this repository;
they are fetched at two fixed commits and checked by sha256:

    python3 -I -B tools/rebuild/fetch_inputs.py <inputs>
    python3 -I -B tools/rebuild/rebuild.py <inputs>

`<inputs>` is a directory of your choice. Expected last lines, exit code 0:

    RESULT: all 9 input files have the expected sha256
    RESULT: REPRODUCED (the uncompressed content is that of the published certificate)

The second command compares the uncompressed JSON. The `.gz` file itself differs from the published one at
least in its gzip time stamp. On the Mac named at the top of this file (Python 3.9.6), the second command
took 70 to 95 s and about 0.6 GB of memory. The rebuild has not been run on another system. It needs no Lean.
The first command needs the network. The last step of the second command, the conversion to the certificate
format, is `tools/gx/gxconv.py`. `tools/gx/ORIGIN.md` names the
outside files and commits, and says what is still not in the repository.

The second command, `tools/rebuild/rebuild.py`, can also be given another file of hand-over pairs; its step
`tools/rebuild/aligned.py` is the script that reads the file (`tools/gx/ORIGIN.md`, "Another file of hand-over
pairs"). A hand-over pair says which helper array that is no longer needed is handed to which new user. What
comes out of that is priced in Python only. It is no part of any proof here.

For the fourth result the Lean files start from the certificate `tools/certificate/gcert1-e8-r783.json.gz`,
which the search programs in `tools/e8/` found (section 10; `tools/gx/ORIGIN.md`).

## 7. File integrity

    shasum -a 256 -c MANIFEST.sha256        # Linux: sha256sum -c MANIFEST.sha256

`MANIFEST.sha256` lists every file of the repository except itself. It is written by

    python3 -I -B tools/community/package_check.py manifest .

which lists every file that git would commit (tracked files and untracked files that `.gitignore` does not
exclude), sorted by name, with its sha256. Each package under `community/` and `package/` also has its own
`SHA256SUMS`, which `scripts/apply-overlay.sh` and `tools/swap-certificate.sh` check before they lay it.

## 8. Limits

Read these before relying on the result.

1. **The sandbox on macOS was a stand-in.** The comparator runs the build inside `landrun` (Linux Landlock).
   Here it ran inside `tools/landrun-macos-shim.sh`, which forbids network access and writes outside the allowed
   directories but, unlike `landrun`, does not restrict what can be read. The sandbox protects the checker
   against a solution file that attacks the build; it is not part of the logical check. None of the comparator
   runs of 2026-10-09 was on Linux. The Linux machine of 2026-10-10 does not close this point. Its kernel was
   built without Landlock. `landrun` refuses to start there without `--best-effort`. With `--best-effort`,
   which is how the comparator calls it, `landrun` starts the command and enforces nothing: in a probe, a
   forbidden write succeeded and the network was open. So no run of the comparator inside its real sandbox
   exists yet for any result of this repository (section 10).
2. **One kernel only.** Every configuration has `"enable_nanoda": false`: the proof was checked by Lean's own
   kernel (when the modules were compiled, and again by the comparator on the exported proof), not by a second,
   independent kernel. `leanchecker` was not run either. A re-run with `"enable_nanoda": true` is welcome.
3. **The comparator replayed files that were already compiled.** Section 3 compiles the solution modules before
   the comparator runs, outside its sandbox. The comparator's README lists as an assumption that the solution
   was not compiled before. To meet it, run section 5 directly after `lake exe cache get`: the comparator then
   builds every module itself. It was not run in that order here; it compiles several kernel checks of a
   certificate at the same time and may need more than 16 GB. On 2026-10-10 it was run in that order on the
   Linux machine, for `comparator/UniformFourier.json` only, and it accepted the third result at the third
   attempt, each attempt continuing the build of the one before (point 14).
4. **Mathlib was not rebuilt from source.** Its compiled files came from Mathlib's cache server for the pinned
   revision (`lake exe cache get`), as is usual.
5. **Not fetched fresh for this test:** the Lean toolchain v4.34.1 and the `comparator` and `lean4export`
   programs were already installed on the machine. For this revision nothing was fetched: the copy started from
   the packages and compiled files of the first publication's rebuild, for which Lake had cloned Mathlib and its
   dependencies and downloaded the compiled Mathlib files into an empty directory.
6. **One machine on 2026-10-09, two on 2026-10-10.** Every run of sections 3 to 5 was on the machine named at
   the top. Section 10 adds one Linux machine.
   Each step was run once, except the comparator on the second result (twice in the working tree, by the agent
   that built the proof and by an auditing agent, and once in the copy) and `tools/Compare.lean` on it (once in
   the working tree, once in the copy). The negative control (a challenge whose saving is one unit larger, which
   must be rejected) was run once, by the auditing agent in the working tree: the comparator and
   `tools/Compare.lean` both rejected it. It was not repeated in the copy.
   Times and memory are single measurements; "resident size" is what `ps` reports once a second and includes
   mapped files.
7. **Files added after the build.** The documentation files (`README.md`, `RELATED-WORK.md`, `NOTICE`,
   `ORIGIN.md`, this file), `tools/whatif_copies.py`, `notes/` and `third-party/` were added or edited after the
   builds and the checks; Lean does not read them. After the last edit all 327 `.lean` files, `lean-toolchain`,
   `lake-manifest.json` and `comparator/B2Gp193x.json` of this repository were compared byte for byte with the
   files of the working tree in which the comparator and `tools/Compare.lean` passed: all identical.
   `lakefile.lean` is the one of the first publication; the working tree has its own, and the copy was built
   with the published one. The build copy was compared with this repository after its comparator run: the same
   `.lean` files, `lakefile.lean`, `lean-toolchain`, `lake-manifest.json`, `comparator/` and `tools/` (without
   `tools/whatif_copies.py`).
8. **What the checks do not say.** They say that the solution proves the statement written in the challenge
   module. Whether that statement means what a reader expects has to be judged by reading the challenge module:
   the cost model in it is OpenAI's, the definitions `wht`, `WHTProgram`, `WHTTimeBoundsAt` are my project's.
   For the third result the whole challenge module is OpenAI's, with the three changes named in section 1.
   The checks say nothing about the speed of any program on inputs of realistic size.
9. **`tools/Compare.lean` is my project's own tool.** It is a second opinion written independently of the
   comparator, not an independent party. The comparator is the outside tool.
10. **The segmented kernel check is new.** The certificate of the second result does not fit into one kernel
    evaluation on this machine. It is cut into segments by a Python generator, each segment is one kernel
    evaluation, and Lean theorems join them. The joining is proved in Lean. The cut is data and is not trusted.
11. **The circuit of the second result is outside data, rebuilt.** The kernel checks the rebuilt certificate as
    it is. That it is the circuit of pull request #193 rests on a Python conversion and on equal counts, not on
    Lean (README, section 3).
12. **The copy of this revision was not built from nothing.** (This point and the next describe the runs of
    2026-10-09. Point 14 says what was done afterwards.)
    It started from the compiled files of the first publication's rebuild, and for the 25 certificate modules of
    the second result it took compiled files of the project's working build, which Lake accepted by their traces
    (section 3). So for the second result no build from nothing has been made; the comparator pass in this copy
    (section 5) replayed compiled files of that origin. A rebuild from nothing with the commands of section 3,
    followed by sections 4 and 5, closes this gap; re-runs are welcome.

13. **The third result.** The comparator passed three times: twice in the project's working tree and once in a
    copy of the repository, where the 15 new modules had been compiled from source. Every run replayed compiled
    files, those of the certificate of the second result among them (point 12). `tools/Compare.lean` passed in
    both places; the negative control was run in the working tree. No build from nothing has been made. A
    rebuild from nothing with the commands of section 3, followed by sections 4 and 5, closes this gap;
    re-runs are welcome.

14. **Points 12 and 13 after 2026-10-10.** The rebuild from nothing with the commands of section 3 was made on
    2026-10-10, on Linux. Section 4 was run on that build. Both are in 10.6, "The published third revision on
    Linux". Section
    5, the comparator, on that machine: the agents ran it with a stand-in for its sandbox, built on Linux
    namespaces (that kernel has no Landlock), in a tree whose modules had been compiled before. State when this
    text was written: all four configurations of the first three results were accepted, each with
    `Lean default kernel accepts the solution`, `Your solution is okay!` and exit code 0. No module was compiled
    inside these runs. The time is the wall clock of the run, and the memory is the maximum resident set size
    that GNU `time` printed:

    * `comparator/B2Ke16x.json` (first result): 2196 s, 6330636 kB;
    * `comparator/B2Ke16xR.json` (first result): 1500 s, 6405516 kB;
    * `comparator/B2Gp193x.json` (second result): 7413 s (2:03:32), 11339956 kB;
    * `comparator/UniformFourier.json` (third result): 6496 s (1:48:16), 11396748 kB.

    The runs shared the machine with each other, so the times are upper bounds. The run for the third result
    and the two runs for the first result used the same copy of the tree, partly at the same time. The run for
    the second result used a second copy. The two negative controls were
    refused with exit code 1 (a challenge with 7474548 for the second result, and one for the Fourier statement
    of the third).

    The agents also ran the comparator with `comparator/UniformFourier.json` in a tree where none of the
    repository's modules had been compiled before. Mathlib's compiled files came from its cache. The comparator
    accepted the third result there: `Build completed successfully (9236 jobs)`,
    `Lean default kernel accepts the solution`, `Your solution is okay!` and exit code 0. That was the third
    attempt in that tree, and each attempt continued the build of the one before. The first two attempts gave
    no verdict: the agents stopped both during the build, to keep the shared machine within its memory. In the
    third, a helper that the agents ran outside the stand-in paused and resumed Lean processes for the same
    reason. The logs of the three attempts have 132, 1 and 181 `Built` lines, so 314 modules of the
    repository were compiled inside the stand-in. The last attempt: 7639 s (2:07:18), 11402808 kB. This verdict
    came after the fourth revision was published; the text published then said that the run had no verdict.

    None of these runs used the comparator's real sandbox (point 1).

## 9. Measurements of the first publication (2026-10-09, 194 modules)

Kept as a record; the commands were those of the first text of this file (commit `024f763`).

| part | time | largest resident size of one Lean process |
|---|---|---|
| the 184 modules that load no certificate (16 `lake build` calls of 12 modules each, `LEAN_NUM_THREADS=4`) | 403 s | 5.8 GB |
| `Gen.Cr2h16` (the data, 3.9 MB of Lean text) | 30 s | 2.3 GB |
| `Gen.Cr2h16Hist` (43 block counts, `decide +kernel`) | 50 s | 0.9 GB |
| `Gen.Cr2h16Lab` (`labelCheck_ok`) | 17 s | 4.5 GB |
| `Gen.Cr2h16Scal` (`scalarCheck_ok`) | 38 s | 5.4 GB |
| `Gen.Cr2h16Shape` (`shapeCheck_ok`) | 9 s | 3.1 GB |
| `Cr2h16Cert` | 63 s | 5.7 GB |
| `InstB2Ke16`, `B2Ke16`, `SolutionB2Ke16x`, `SolutionB2Ke16xR` | 8 to 9 s each | 5.3 to 5.5 GB |
| **whole build** | **644 s (about 11 minutes)** | |

`tools/Compare.lean`: both `RESULT: PASS`, exit code 0, 140 s and 117 s; a repeat of the first with a memory
sampler: 144 s, 6.2 GB. Comparator (macOS, the stand-in, `enable_nanoda` false): `comparator/B2Ke16x.json`
171 s, 6.3 GB; `comparator/B2Ke16xR.json` 178 s, 6.2 GB; both "Lean default kernel accepts the solution" /
"Your solution is okay!", exit 0. Disk then: `.lake/build` 1.3 GB, about 9.5 GB in all.

## 10. The fourth result (2026-10-10)

This section is how to check the fourth result yourself: the three theorems `wht_main_block_B2Ge8x`,
`transform_mainE8` and `convolution_mainE8`.

* **What you check.** You build the new Lean modules from source (10.1). You confirm with two tools that each
  theorem proved is the theorem stated (10.2, 10.3). You regenerate the generated files and check the
  certificate in Python, without Lean (10.4). The certificate is the data file that lists the additions of the
  circuit; the Lean kernel checks it.
* **Operating system.** macOS or Linux. "Where it was run", below, names the two machines.
* **Memory.** Plan for 16 GB free for one run of `tools/Compare.lean` on Linux. On a 16 GB Mac, run nothing
  beside it (section 2).
* **Time.** The measured times of the Lean steps are in 10.1 to 10.3. The two Python checkers of 10.4 answer
  in under a second on the Mac.
* **Before you start.** Section 2 lists what you need. Section 3 has the first build command, which 10.1
  starts from.

**What the fourth result adds.** It is 39 Lean modules under `Work/`, on top of the 252 modules that the
first three results have there. (10.6 counts 342 modules for the published third revision. That figure is the
whole repository: the 252 under `Work/`, the 89 unmodified files of OpenAI under `OAI/`, and one module under
`WHTCheck/`.) The 39 are:

* 24 are generated from the certificate `tools/certificate/gcert1-e8-r783.json.gz`. The two generators of the
  second result make them, and the generators are unchanged (`tools/gx/gxgen.py`, `tools/gx/gxrate.py`).
  14 of the 24 are under `Work/GCert/Data/Gen/` and hold the data and its kernel evaluations. The other 10 are
  above them, up to `Work.GCert.Data.SolutionB2Ge8x`.
* 15 are `Work/FourierE8/`, the twin of `Work/Fourier/`.

No published file was changed for them, and `lakefile.lean` is the published one.

**Where it was run.** The AI agents I directed made the runs of this section on 2026-10-10, on two machines.
The build (10.1), `tools/Compare.lean` (10.2), the comparator (10.3) and the two regeneration scripts of 10.4
were run on both. The comparator had a stand-in for its sandbox on both. The two Python checkers of 10.4 were
run on the Mac. 10.5 says where each negative control was run. 10.6 is a run on the Linux machine only. The
two machines:

* the Mac named at the top;
* a rented Linux machine: x86-64, Ubuntu 24.04.5, 8 logical processors on 4 cores, 51 GiB of memory
  (`MemTotal: 53467192 kB`), no swap.
  It is a Docker container. Its kernel, 6.6.122+, was built without Landlock, which the comparator's sandbox
  program `landrun` needs (section 8, point 1). It runs no systemd.

sha256 of the certificate file: `4b92f00fc7454b9b15e6d71a3b05062792eacca0a0ca4efee46f87827d57c8aa`; of the
uncompressed JSON: `3594f19c4d01cf11fb930c5c61baeed399620acbc5b94096b16bcacb23997dd3`.

### 10.1 Build

This step compiles the 39 modules from source and prints the axioms that each theorem depends on.

Run the first command of section 3 first; it builds the modules that load no certificate. Then build the new
modules one at a time:

    for m in Gen.E8.S0 Gen.E8 Gen.E8Par Gen.E8Tab0 Gen.E8Tab1 Gen.E8Tab2 Gen.E8Tab3 Gen.E8Seg0 Gen.E8Seg1 \
             Gen.E8SegF Gen.E8End Gen.E8Hist Gen.E8YChk Gen.E8Scal0 RateB2Ge8 RateB2Ge8x InstB2Ge8 E8Price \
             E8Cert E8ScalOK B2Ge8 B2Ge8Main ChallengeB2Ge8x SolutionB2Ge8x; do
      lake build Work.GCert.Data.$m || break
    done
    lake build Work.FourierE8.Axioms Work.FourierE8.UniformFourierChallenge Work.FourierE8.Main

The certificate modules of the first and second results need not be built first. The reason: by their import
lines, the 24 modules import the checker and the engine of the second result, and no data module of another
certificate.

What to expect:

* Every `lake build` must end with `Build completed successfully`, and no line may start with `error`.
* The warning "declaration uses `sorry`" must appear for the two challenge modules, which state the theorems
  without proof: one line for `Work.GCert.Data.ChallengeB2Ge8x`, two lines for
  `Work.FourierE8.UniformFourierChallenge`. It must appear for no other of the 39 modules.

The builds of `Work.GCert.Data.SolutionB2Ge8x` and `Work.FourierE8.Axioms` print, among other lines:

    'OAI.PowerSaving.WHT.wht_main_block_B2Ge8x' depends on axioms: [propext, Classical.choice, Quot.sound]
    'OAI.PowerSaving.transform_mainE8' depends on axioms: [propext, Classical.choice, Quot.sound]
    'OAI.PowerSaving.convolution_mainE8' depends on axioms: [propext, Classical.choice, Quot.sound]
    'OAI.PowerSaving.RAM.hillsE8_program' depends on axioms: [propext, Classical.choice, Quot.sound]

The kernel checks the certificate in segments, in the same way as the certificate of the second result. Each
large evaluation has its own module:

* `Gen.E8Tab0` to `Tab3`: the frame table;
* `Gen.E8Seg0`, `Seg1` and `SegF`: the replay of the gates;
* `Gen.E8End`: end conditions;
* `Gen.E8Scal0` and `Gen.E8YChk`: the scalar half;
* `Gen.E8Hist`: the block counts.

The data is 0.24 MB of Lean text; the second certificate has 7.5 MB. So no module of the fourth result is
heavy by itself.

Measured, Mac: the 24 modules, one per `lake build`: 24 of 24 with exit code 0, 183 s in all, largest resident
size of one Lean process 6.06 GB; the 15 modules of `Work/FourierE8`, one per `lake build`: 15 of 15 with exit code
0, 168 s in all, largest resident size 6.15 GB; each of the three axiom lines ends with
`depends on axioms: [propext, Classical.choice, Quot.sound]`; `lake build --no-build` of the 16 roots (the 11
published roots and the 5 new ones) ended with `All targets up-to-date (9304 jobs).` and exit code 0.

Measured, Mac, a second run: the same 39 modules were built from the files of this revision alone, in a copy
of the complete build of the published tree, with the commands above as printed. The loop: 24 of 24 with
exit code 0, 146 s in all, largest resident size 6.05 GB. The Fourier command: `Build completed successfully
(9227 jobs).`, exit code 0, 47 s, 6.16 GB. The same three axiom lines, and `lake build --no-build` of the 16
roots again ended with `All targets up-to-date (9304 jobs).`.

Measured, Linux: 39 of 39 `lake build` calls ended with `Build completed successfully` and exit code 0;
1185 s in all (476 s for the 24 modules, 709 s for the 15); largest process 6.69 GB; the three axiom lines as
on the Mac. The machine was shared with other jobs (load average 16 to 20 on 8 threads), so its times are
upper bounds.

### 10.2 `tools/Compare.lean`

This step compares each challenge module with its solution module, using this repository's own tool
(section 4). It shows three things: the solution proves the same statement as a kernel term, every definition
the statement depends on is identical, and only the permitted axioms are used.

    lake env lean --run tools/Compare.lean Work.GCert.Data.ChallengeB2Ge8x \
        Work.GCert.Data.SolutionB2Ge8x OAI.PowerSaving.WHT.wht_main_block_B2Ge8x
    lake env lean --run tools/Compare.lean Work.FourierE8.UniformFourierChallenge \
        Work.FourierE8.Main OAI.PowerSaving.transform_mainE8 OAI.PowerSaving.convolution_mainE8

Expected, exit code 0, for each theorem:

    OK: statement is the same kernel term in both
    OK: every statement dependency is identical in the solution
    OK: only permitted axioms

    RESULT: PASS

Measured, Mac: the first command printed
`solution proof reaches 36821 constants; axioms: #[Classical.choice, Quot.sound, propext]` and
`RESULT: PASS`, exit code 0, 107 s, largest resident size 6.61 GB. The second printed 41331 constants and
41354 constants with the same axioms and `RESULT: PASS`, exit code 0, 103 s, 7.32 GB. The footprint, sampled
every four seconds, peaked at about 9.8 GB in both runs. Both went to the end on the 16 GB Mac; the two runs
of section 2 that were stopped at 9.3 and 9.7 GB were stopped by the auditing agents' own memory rule, not
by the machine. Measured, Linux: the same lines and counts, exit code 0; 99 s and 15.7 GB for the first command,
150 s and 16.0 GB for the second.

### 10.3 The official comparator

This step makes the same comparison with the official comparator (github.com/leanprover/comparator;
section 5). It shows that each solution proves the same statement as its challenge, uses only the permitted
axioms, and is accepted by the Lean kernel.

    lake env comparator comparator/B2Ge8x.json
    lake env comparator comparator/UniformFourierE8.json

Expected last lines of each run, exit code 0, as in section 5:

    Lean default kernel accepts the solution
    Your solution is okay!

Measured, Mac (the stand-in of section 5, `enable_nanoda` false): for both configurations
`Lean default kernel accepts the solution` and `Your solution is okay!`, exit code 0; 85 s and 6.05 GB for
`comparator/B2Ge8x.json`, 89 s and 6.09 GB for `comparator/UniformFourierE8.json`; 171 modules replayed, none
built. Linux: the same two lines for both configurations, exit code 0; 406 s and 6.80 GB, 432 s and 6.68 GB;
171 `Replayed` lines and no `Built` line in each log. A memory figure in GB for the Linux machine, here and in
the rest of this file, is the kbytes figure that GNU `time` printed, divided by a million; a kbyte there is
1024 bytes, so in units of 10^9 bytes each figure is 2.4 % higher. The sandbox there was a stand-in built on
Linux namespaces, not `landrun`, because that kernel has no Landlock; the stand-in was not tested on its own
in these runs. This is not a run with the real sandbox. The times are from a machine under load.

### 10.4 Regenerate the generated files, and check the certificate in Python (no Lean needed)

This step has two parts. The first shows that the generated files in the repository are what the generators
write. The second runs two Python checkers on the certificate.

**Part 1: regenerate.**

    python3 tools/gx/regen_check_e8.py
    python3 tools/fourier/regen_check_e8.py

Expected last lines, exit code 0. Both were obtained on the files of this revision:

    RESULT: REPRODUCED (25 repository files reproduced byte for byte, 0 different, 3 generated files not in the repository)
    RESULT: REPRODUCED (14 repository files reproduced byte for byte, 0 different, 0 generated files not in the repository)

What the two lines count:

* The 25 files are the 24 Lean modules and `comparator/B2Ge8x.json`.
* The three generated files that are not in the repository are the eight-decimal variant: `ChallengeB2Ge8`,
  `SolutionB2Ge8` and its comparator configuration. It is the same for the second result.
* The 14 files are the 11 copies of OpenAI's files, `Goal.lean`, the challenge and
  `comparator/UniformFourierE8.json`.

These scripts do not generate `Work/FourierE8/Seam.lean` and `Axioms.lean`. Those two are the files of the
same name in `Work/Fourier/`, with names and numbers replaced (ORIGIN.md, last section).

**Part 2: check the certificate.**

    python3 -B tools/gx/refcheck.py tools/certificate/gcert1-e8-r783.json.gz
    python3 -B tools/gx/gxdry.py tools/certificate/gcert1-e8-r783.json.gz

Expected, among the lines printed:

    ACCEPTED by gx.check1: h=9 v=120 R=783 N=9039 gates A=773 B=1220; ...
    ... MIRROR ACCEPTS (label side E0 E1 E2 E6, scalar side E3 E4 E5) ...
    unit of B_2: m=45 W=1263 D=120 tally=56715 ...
    10 decimals, fill 1 - 2^-40: whole-block 8762479 (evaluated 8762482 holds / 8762483 fails); ...

On the Mac both checkers answer in under a second: about 0.2 s for `refcheck.py` and about 0.25 s for
`gxdry.py`, with Python 3.9.6.

The words "evaluated 8762482 holds / 8762483 fails" report the
rate inequality as Python evaluates it: it holds at 8762482 and fails at 8762483. The theorem that Lean checks
states 8762479, which is three units lower. The three larger figures, 8762480, 8762481 and 8762482, are
neither proved in Lean nor claimed.

### 10.5 Negative controls: what was run, and on which design

A negative control gives the checks an input that is wrong on purpose; the checks must refuse it. Words used
below:

- a unit is the finite network that the saving is computed from;
- a helper is an extra array that holds partial sums;
- a frame is the subspace of labels where an array stands; it only grows;
- the scatter is the point of the circuit where a few helpers, the totals (8 here), are added into all
  targets;
- an in-place pair forms a + b and a - b on the two arrays that held a and b.

The controls, in four groups:

- **In Python, on the certificate of this revision.** A checking agent made seven damaged copies of it:

  1. a wrong coefficient in an in-place pair;
  2. the compensating read of a re-used helper removed;
  3. two helpers folded into one;
  4. a total declared one dimension lower;
  5. two labels exchanged;
  6. a source added with 1 in place of 1/2;
  7. three blocks booked as one cheaper block.

  `refcheck.py` and `gxdry.py` refuse all seven. So do two programs that read only the certificate: a
  stand-alone replay, and the checking agent's own program. All four accept an unchanged copy.

  One more copy carries the circuit to another family of 120 labels. `refcheck.py` and `gxdry.py` accept it,
  with the same figure. The two other programs test the family, and they refuse it. So the two published
  checkers, like the theorem, do not say which family the labels are.
- **In Lean, on two earlier units of the same day, not on this certificate.** The units at 7558959 and at
  8320202 were built in Lean during the work. Their Lean files are not included.

  - The first unit, on the Mac. An auditing agent ran damaged inputs through Lean. Five damaged inputs and two
    raised figures were refused, and the intact control was accepted. The official comparator accepted the
    right pair and refused both challenges with a larger saving. `tools/Compare.lean` could not be run to its
    end in that audit, because of memory (section 2).
  - The second unit, on Linux. Another agent ran 21 builds. It had written down the outcome of each
    beforehand: 13 had to fail and failed, 8 had to build and built. `tools/Compare.lean` passed the true pair
    and failed both challenges with a larger saving.
- **One damaged record goes through Lean, in both of those audits.** The certificate has a list of totals.
  Each entry has a helper and a frame number. Change the frame number of one entry and nothing else. Then:

  - the Python checkers refuse the file, and the generator writes nothing;
  - the Lean text with that one number changed builds to the end, to the unchanged theorem.

  No Lean check was found that reads that number: the frame of a total is taken from the replayed state.
  Lean refuses a version that tries to profit from the wrong frame. For nine more fields of the certificate, the auditing
  agent's text search of the Lean sources found no place that reads them. That is a search by pattern, not a
  test: a use under another name would escape it, and none of the nine was changed and built. The theorem is
  not affected. The point is that "Lean accepted
  the instance" says nothing about those fields of the JSON file.
- **In Lean, on the certificate of this revision.** The agents ran these controls on the Mac on 2026-10-10,
  in a separate copy of the build tree. The damaged files are not in this repository. The file names below
  say where each change sits; the agents made the changes in copies of those files under other module names.

  - **Twelve wrong Lean modules were refused.** None of them left a compiled file.

    - Six had one wrong coefficient each in the data module, `Work/GCert/Data/Gen/E8/S0.lean`. Three were
      in in-place pairs on a helper that holds halves. One was in the addition of a source into a helper
      with the coefficient 1/2. One was in an addition between two helpers at a frame of dimension 8, after
      the scatter. One was in the scatter itself. Each text differs from the accepted text in one token.
      Lean's scalar check refused all six: the theorem `scal0_ok` in `Work/GCert/Data/Gen/E8Scal0.lean`.
    - Two had a wrong step on the label side. In one, a helper is sent to a line that is not inside its
      next frame. In the other, a helper that climbs from a line to the full frame in one block has a wrong
      rank (the rank of a block is the number of dimensions it climbs; 7 was written for 8). Lean's label
      check refused both: the theorem `seg0_ok` in `Work/GCert/Data/Gen/E8Seg0.lean`.
      It built a legal twin of the first, in which the same helper goes to a line inside that frame.
    - Three claimed a larger saving. `Work/GCert/Data/RateB2Ge8x.lean` did not build with every `8762479`
      replaced by `8762483`, and did not build with `8800000`. A copy of the solution module that states
      8762483 and cites the real theorem failed with a type mismatch.
    - The twelfth tried to profit from the gap described below.
  - **Accepted, as it must be:** the unchanged text under other names. It built to the theorem at
    `1 - 8762479/10^10`, with the axioms `propext`, `Classical.choice` and `Quot.sound`.
  - **`tools/Compare.lean` printed `RESULT: FAIL`** for a challenge at 8762483 against the real solution
    module, with the line `FAIL: statement (type) differs between challenge and solution`.
  - **The Python side.** `refcheck.py` and `gxdry.py` refused every damaged certificate, and the published
    generator `gxgen.py` wrote no Lean text for any of them. So the damaged Lean text came from an
    unpublished driver, or from an edit of generated text. The agent that ran the driver describes it as the
    generator's own functions with the Python test switched off.
  - **The known gap is present here too.** Lean does not use the frame number stored in an entry of the
    list of totals. That list is `ret` in the certificate and `ret_0` on line 30 of `S0.lean`. With the
    first entry changed from `(504, 1)` to `(504, 122)` and nothing else, all 24 generated modules built to
    the unchanged theorem, with the three axioms. `refcheck.py` and `gxdry.py` refuse the same change in the
    certificate, with the line `REJECTED: E2: retained slot 504 is not at its frame`. A run showed this for
    the first entry. For the other seven it rests on a reading of the Lean sources, which did not cover
    every function.

    What the gap means: that field is checked by the Python checkers only. The theorem that Lean proves is
    unchanged, and the entries of the published certificate are correct. The twelfth refused module was the
    attempt to profit from it: all 8 totals declared one dimension lower, and priced so. Lean built the 13
    modules before it and refused the check `cen_ok` (line 28 of `Work/GCert/Data/E8Cert.lean`). That check
    compares the declared dimension with the frame where the totals stand in the replayed circuit.

  Limits of these controls:

  - **Coverage.** They touch 3 of the 21 in-place pairs on a helper that holds halves, 1 of the 4 additions
    with 1/2, 1 of the 5 additions at dimension 8, and 1 of the 2 climbs from a line to the full frame.
    None is aimed at the table depth 10, the digit width 7 or the largest digit 36 of this certificate.
  - **What a refusal says.** A refusal by Lean says only that a check is false, not which test inside it
    failed. That the intended test failed rests on two things: the one-token difference from text that
    builds, and a Python replay of the same check.
  - **The raised figures.** That the rate lemma does not build at 8762483 shows that this proof does not go
    through there. That the inequality itself fails at 8762483 is a separate theorem of this repository:
    `fold_B2Ge8x_sharp`, line 124 of `Work/GCert/Data/RateB2Ge8x.lean`.
  - **Not covered.** There was no control on the modules of `Work/FourierE8`, none on Linux, and the
    official comparator was not run on a wrong pair. One kernel.
  - **Who checked it.** One agent ran the controls. A second agent audited that record. It compared the
    damaged certificates and the Lean text that was built with the originals, and it read the logs. It ran
    no Lean and did not repeat the Python runs. It found no wrong statement in the record. By the time
    stamps in the first agent's file, the expected outcome of each control was written down before its run;
    those stamps are that agent's own. No person has read either record.

### 10.6 The published third revision on Linux

On 2026-10-10 the published third revision was cloned on the Linux machine and built from nothing with the
three build commands of section 3 (commit `f010392`). Mathlib was not compiled from source: its compiled files came from its
cache server. All 342 Lean modules of the repository were compiled from source: 252 under `Work/`, the 89
files of OpenAI under `OAI/`, and one under `WHTCheck/`. Then the four commands of section 4 were run on that
build: the three listed at its start and the one for the Fourier pair.

The build:

* 44 `lake build` calls. Each ended with `Build completed successfully`, exit code 0, and no line starting
  with `error`.
* Total time 5644 s (94 minutes): 758 s for the first command, 4781 s for the loop, 105 s for the Fourier
  command.
* Largest single process: 10.5 GB (`Gen.P193Scal0`, `Gen.P193Scal1`).
* "declaration uses `sorry`" appeared 5 times, in the 4 challenge modules and nowhere else.
* The axioms of all five published theorems were `[propext, Classical.choice, Quot.sound]`.
* The slowest modules were `Gen.P193Seg3` and `Gen.P193Seg4`: 757 s and 1032 s, against 156 s and 186 s on
  the Mac.

`tools/Compare.lean` on that build:

* It printed `RESULT: PASS`, exit code 0, for all four commands of section 4.
* Section 4 gives an expected constant count for three of the four commands, and those three counts were
  printed: 38678 for the second result, 36438 and 33452 for the two theorems of the first. For the fourth
  command, the Fourier pair of the third result, section 4 gives no expected count; the run printed 43188
  and 43211.
* Time: 128 s and 155 s for the second and third results, 47 s and 45 s for the two of the first.
* Largest process: 15.7 to 16.0 GB (section 2).

The comparator on that build: point 14 of section 8.

### 10.7 Limits particular to the fourth result

1. **No run inside the comparator's real sandbox.** On macOS the sandbox is a stand-in. On the Linux machine
   the kernel was built without Landlock. Section 8, point 1, has both.
2. **One kernel.** `enable_nanoda` is false in both new configurations, so no second, independent kernel
   checked the proofs (section 8, point 2). `leanchecker` was not run.
3. **The modules were compiled before the comparator ran**, as in section 8, point 3. No comparator
   run on the fourth result was made in a tree where nothing had been compiled.
4. **No build of the whole fourth revision from an empty directory.** Such a build starts with nothing
   compiled and compiles every module of the revision in one run. None has been made, on either machine.
   What was made is this. On the Mac the 39 modules were compiled from source in a tree that
   held the compiled modules of the earlier results. That tree also held the source and compiled files of the
   two earlier units of the same day, which are not in this repository; no new module imports them. The 39
   modules were therefore built a second time on the Mac from the files of this revision alone, in a copy of
   the complete build of the published tree, with the same outcome (10.1). For that second build the files were exported from git, and
   the tree held none of the working files of the day. On Linux the 39 modules were compiled from source on
   top of that machine's own build from nothing of the published tree (10.6). So on Linux every module of the
   revision was compiled from source on that machine, in two steps and not in one run. The revision has 391
   files that Lean reads: the Lean sources, the comparator configurations, `lakefile.lean`,
   `lake-manifest.json` and `lean-toolchain`. By sha256, all 391 are identical to the files in the trees
   where the builds and checks were run, on the Mac and on Linux.
5. **Negative controls in Lean, and fields that Lean does not read.** Negative controls in Lean were run
   on two earlier units and on this certificate (10.5). Those on this certificate were run by one agent, on
   the Mac. A second agent audited the record and ran no Lean. They cover one to three instances of each new
   kind of step, not all of them. None is on the modules of `Work/FourierE8`, none was run on Linux, and
   none uses the official comparator. One field of the certificate is
   shown by a test not to be read by Lean: the frame number of a total. For it, the Python checkers refuse a
   changed value and Lean builds all the same. On this certificate the test changed the first of the eight
   totals; for the other seven there is a reading of the Lean sources only. For nine more fields the evidence is weaker: a text search of
   the Lean sources found no place that reads them, and a search by pattern can miss a use under another
   name (10.5).
6. **The search programs are not checked by anything.** `tools/e8/` shows how the certificate was found. The
   theorem rests on the certificate as the kernel checks it, not on those programs. That the certificate is
   what those programs produce was reproduced: `python3 -B tools/e8/run.py`, run on the files of this
   revision, starts from the recorded hierarchy (the output of the solver stage), writes the JSON of the
   certificate byte for byte and ends with `RESULT: REPRODUCED` (76 to 98 s on the Mac, five runs). With
   `--pool --solve` it also rebuilds the pool and re-runs the solver, which needs NumPy and SciPy; two runs
   gave the recorded hierarchy and the same bytes (95 s and 107 s).
7. **One day.** The certificate was found, built and checked on 2026-10-10. The measurements are single runs.

## 11. The results built on the community's Lean tree (2026-10-13)

This section is how to check the results of 13 October 2026 yourself. They are not built in this repository's
own Lean project. They are laid over Ryan Shea's Lean tree, the folder `lean/complex588` of his repository
`fourier-transform-below-nlogn`, which runs on the Lean engine of this repository. Each result is a *package*
under `community/`: a folder `overlay/` that mirrors `lean/complex588/`, the list of its new modules in build
order (`MODULES.txt`), a file `result.json` with its figures, theorem names and the status of every check,
`SHA256SUMS`, and the logs of the runs.

| package | what is proved | saving | status |
|---|---|---|---|
| `community/u8xmhh` (the headline) | the Fourier transform of every length, and convolution: OpenAI's statement with `1 - 16197705/10^10`; the community's U8 word with exit frames, mutable sources and same-run hosting | 1.6198e-3 | checked in Lean on the second Colab machine (11.2), comparator accepted (11.3); the three phases of the leaderboard's verifier re-created offline with the agents' scripts, in the stand-in sandbox; statement re-checked and modules rebuilt by two independent agents |
| `community/c588mhh` | the same statement with `1 - 15699305/10^10`: Kothe's circuit with mutable sources and same-run hosting | 1.5699e-3 | checked in Lean (11.2d), comparator accepted (11.3), rebuilt by an independent agent; the leaderboard's verifier not run |
| `community/c588m` | the same statement with `1 - 13552606/10^10`: the same circuit with mutable sources, no hosting | 1.3553e-3 | checked in Lean (11.2a); the three phases of the leaderboard's verifier re-created offline with the agents' scripts, in the stand-in sandbox (11.6, point 4); rebuilt by an independent agent |
| `community/u8xm` | the same statement with `1 - 13868856/10^10`: the community's U8 word with exit frames and mutable sources | 1.3869e-3 | checked in Lean, comparator accepted (11.2c); the agents' second build and two independent agents' checks not shown |
| `community/u8x` | the same statement with `1 - 12704023/10^10`: the U8 word with exit frames | 1.2704e-3 | checked in Lean, comparator accepted (11.2c); the agents' second build and two independent agents' checks not shown |
| `community/c588-repack` | the same statement with `1 - 12383291/10^10`, on Shea's certificate with a new bank table | 1.2383e-3 | checked in Lean, comparator accepted; an independent agent's re-check of the built tree (its comparison one unit stronger not shown) |
| `community/wht-c588mhh` | the Walsh-Hadamard transform, this repository's statement, with `1 - 15699306/10^10`, on the kernel of `community/c588mhh` | 1.5699e-3 | checked in Lean, comparator accepted (11.2b); no independent rebuild shown |
| `community/wht-c588` | the Walsh-Hadamard transform, this repository's statement, with `1 - 12331686/10^10`, `1 - 12383292/10^10` and, on the kernel of `community/c588m`, `1 - 13552607/10^10` | 1.2332e-3, 1.2383e-3 and 1.3553e-3 | checked in Lean, comparator accepted; an independent agent's rebuild of the first (not shown for the second, not run for the third) |

`community/current.json` names the headline package. To replace it by another one, see 11.5. `u8xmhh`,
`c588mhh`, `u8xm` and `u8x` have a `README.txt`; every package's `result.json` lists every check with its status. A check marked
"pass" gives the path of its log inside the package; every other check is marked "not shown" (it ran, but its log
is not in the package or has no verdict) or "not run". The static checks of 11.5 (check 9) refuse a package that
breaks this rule.

### 11.1 What you need

* A Linux machine the size of a Colab VM: the agents used 8 logical processors and 51 GiB of memory. One run
  of `tools/Compare.lean` took 14.9 to 16.1 GB (resident size by `/usr/bin/time -v`); the largest module took
  9.6 GB.
* Shea's tree, built. His repository pins Lean v4.34.1 and the same Mathlib as this repository. Building his
  tree from source took the agents 25 minutes; Mathlib's compiled files came from its cache server:

      git clone https://github.com/shea256/fourier-transform-below-nlogn shea
      cd shea && git checkout --detach 52394849c7ecf3412d39b925d75fc7fdb98994d8 && cd lean/complex588
      lake exe cache get
      lake build Work.FourierC588.Main Work.FourierC588.UniformFourierChallenge \
        Work.FourierC588Neg.UniformFourierChallenge Work.GCert.Data.C588Controls Work.FourierC588.Axioms

  `community/base-logs/` holds the logs of that build and the scripts the agents used (`scripts/cb-setup.sh`,
  `scripts/cb-build.sh`).
* Work on a copy of the built tree for each package (`cp -a shea shea-c588m`), so that one package never
  meets another.

### 11.2 Build and check the main result

The main result, `community/u8xmhh`, puts the hosted network of `community/c588mhh` (11.2d) on the certificate
of `community/u8xm` (11.2c): the community's U8 word with exit frames and mutable sources, with most leftover
helper pieces stored in data registers that are idle during the same run, instead of in banks (310 banks
instead of 3,771). Its twelve network modules of `Work/XHost` are the same files, byte for byte, as in
`community/c588mhh`. New are the host table and bank table of the U8 word (2,016 hosts and 4,092 hosted pieces
in each of the 8 runs), which the kernel checker `hostSchedK` reads (`Work/XHost/Data`), the glue
`Work/XHost/U8XMHH.lean`, and the Fourier chain `Work/FourierU8XMHH`, with its seam pointed at
`XHost.engine_of_ginv_host`. From the root of this repository, with `T` the folder `lean/complex588` of your
copy, after adding the `XExit` line of 11.2c to `$T/lakefile.lean`:

    bash community/scripts/apply-overlay.sh community/u8xmhh $T
    bash community/scripts/build-and-check.sh community/u8xmhh $T

The first command lays `community/u8xm` first (the package's `result.json` names it under "requires"), then the
34 files of this package (32 modules and 2 comparator configurations); no file of the tree changes except
`lakefile.lean`. The second builds the 32 modules in the order of `MODULES.txt` and runs the checks listed in
11.2a, against the statement `1 - 16197705/10^10` and the one unit stronger, `16197706`.

The agents' run is in `community/u8xmhh/logs/colab/`, on the second Colab machine, in a copy of a tree in which
they had built `community/u8xm` a second time and checked it (`setup.out`: the 96 files of that overlay equal
to its checksums, every target up to date, and that build's comparator verdicts; its other logs are not in the
package). `build.out` has the six Lake steps (the rate certificate, the checker, the host table, the network,
the glue, the Fourier chain), each ending without error, 15 minutes in all and 8.6 GB at most; then the axioms
(`axioms.log`, `axioms-glue.log`: `propext`, `Classical.choice`, `Quot.sound`, and `propext` alone for the
host-table check), the source scan (`sorry-grep.txt`), the check that both statement files equal the agents' own
derivation from OpenAI's file, a negative control (a host table with one remainder off by one, for which Lean's
kernel proves that `hostSchedK` returns false: `hsneg.time`, `ax/HsU8XMHHNeg.lean`), and the two comparisons
(`compare-pos.log`: `RESULT: PASS`, 86 s, 15.8 GB; `compare-neg.log`: `RESULT: FAIL` on
`decimalExponentU8XMHH`). The comparator's runs are in `cmp/` (11.3), and the offline re-creation of the
leaderboard's verifier, done as in 11.6, point 4, on this package, is in `bos-H.out` and `bos-H/`: the Lean
check accepted it (68 minutes, 8.4 GB, every module of the packaged submission compiled from source), the axiom
audit printed the three standard axioms, and nanoda and Lean's kernel accepted it (14 minutes).

Two independent agents then checked it. One, on the same machine in its own copy of the tree, compared both
statement files with OpenAI's file, checked the theorems with plain `lean` as an unprivileged user, ran
`tools/Compare.lean` both ways and replayed the comparator, which accepted it and rejected one unit stronger
(`logs/skeptic-statement/`). The other compiled the 32 modules again on the first Colab machine from the two
archives the agents had uploaded (their hashes checked), one module per Lake call, over its own copy of the
agents' first build of `community/u8xm` there, with no compiled hosting or Fourier file before it: all compiled,
`RESULT: PASS`, one unit stronger `RESULT: FAIL`, standard axioms; `uses.log` shows that both theorems reach
`XHost.engine_of_ginv_host` and `hostSchedK` (`logs/skeptic-rebuild/`). It did not run the comparator. The
`report.txt` in each folder is the agent's written verdict.

### 11.2a The package with mutable sources only

From the root of this repository:

    bash community/scripts/apply-overlay.sh community/c588m $T
    bash community/scripts/build-and-check.sh community/c588m $T

The first command checks the package against its `SHA256SUMS` and copies its 72 files into the tree. Every
file is new: it stops before copying anything if a file of the tree would change. The second builds the 70
new modules one at a time, in the order of `MODULES.txt`, then runs the checks and prints one line per check:

* `compare: RESULT: PASS`: `tools/Compare.lean` finds the theorem proved identical to the statement. The
  statement is OpenAI's `UniformFourier.lean` with only the exponent and seven names changed;
* `compare-negative`: the same comparison against the statement one unit stronger (`13552607`) fails on the
  exponent, as it must;
* `axioms`: the two theorems and the program behind them use only `propext`, `Classical.choice` and
  `Quot.sound`;
* `source scan`: no `sorry`, `admit`, `axiom` or `native_decide` in the new files, except the two `sorry` of
  each statement file; no file of Shea's tree changed.

It prints `ALL CHECKS PASSED` at the end. The agents' run of these steps is in `community/c588m/logs/colab/`
(`mb-build.out`, `axioms.log`, `compare-fourier.log`, `compare-negative.log`, `source-scan.txt`). The 70
modules took 29 minutes, 36 with the axioms and the two comparisons.

The module `Work.GCert.Data.C588MControls` is a set of negative controls inside Lean: a certificate with one
sign changed, the classic check (which forbids writing the sources) and a bank table with an owner used twice are each
rejected by the kernel.

### 11.2b The Walsh-Hadamard corollaries (`community/wht-c588` and `community/wht-c588mhh`)

The theorem `OAI.PowerSaving.WHT.wht_main_block_C588Mx` (saving 1.3552607e-3) is one line on the kernel program
of `community/c588m`, so build that package first, at least up to `Work.FourierC588M.Seam`. Its files are in
`community/wht-c588/overlay/` (folders `Work/WHTC588M`, `Work/WHTC588MNeg`, and the three `WHTC588M*.json`
comparator configurations). The statement file `Work/WHTC588M/Challenge.lean` is this repository's
`Work/GCert/Data/ChallengeB2Ge8x.lean` with the theorem name and the exponent changed. Then, in `$T`:

    cp -R <this repository>/community/wht-c588/overlay/Work/WHTC588M <this repository>/community/wht-c588/overlay/Work/WHTC588MNeg Work/
    cp <this repository>/community/wht-c588/overlay/comparator/WHTC588M*.json comparator/
    lake build Work.WHTC588M.Challenge Work.WHTC588MNeg.Challenge Work.WHTC588M.Main Work.WHTC588M.Solution
    lake env lean Work/WHTC588M/Axioms.lean
    lake env lean --run <this repository>/tools/Compare.lean Work.WHTC588M.Challenge Work.WHTC588M.Solution OAI.PowerSaving.WHT.wht_main_block_C588Mx
    lake env lean --run <this repository>/tools/Compare.lean Work.WHTC588MNeg.Challenge Work.WHTC588M.Solution OAI.PowerSaving.WHT.wht_main_block_C588Mx

The first comparison must print `RESULT: PASS`, the second (one unit stronger, `13552608`) `RESULT: FAIL`. The
agents' run, on the second Colab machine, used `community/wht-c588/scripts-as-run/c588m/rwht-run.sh`; its logs
are in `community/wht-c588/logs/colab/c588m/` (`run.out` prints `ALL CHECKS PASSED`; `cmp.out` has the
comparator runs: both kernels accept, the configuration one unit stronger is rejected). The two other
Walsh-Hadamard corollaries of `community/wht-c588` are built the same way on their own kernels
(`WHTC588` on Shea's unchanged tree, `WHTC588Q` after `community/c588-repack`).

The corollary on the kernel of `community/c588mhh`, `OAI.PowerSaving.WHT.wht_main_block_C588MHHx` (saving
1.5699306e-3), is in `community/wht-c588mhh`, built the same way after `community/c588m` up to
`Work.FourierC588M.Seam` and the 17 files of `community/c588mhh` that it needs (16 of `Work/XHost` and
`Work/FourierC588MHH/Seam.lean`, listed with their checksums in
`community/wht-c588mhh/scripts-as-run/hosting-sources.MANIFEST.sha256`). The agents' run, on the second Colab
machine, used `community/wht-c588mhh/scripts-as-run/rfwht-run.sh`: 21 modules compiled with plain `lean`, then
the axioms, both comparisons and the comparator (`logs/colab/run.out` prints `ALL CHECKS PASSED`; `cmp.out`
has the comparator runs, accepted with both configurations and rejected one unit stronger, `15699307`). No
Walsh-Hadamard corollary has been built on the kernel of the main result, `community/u8xmhh` (16197706).

### 11.2c The packages on the U8 word

`community/u8x` and `community/u8xm` are laid over Shea's tree directly, in a copy of their own. Both contain the
five modules of `XExit` (exit frames), a Lean library of its own: first add the line

    lean_lib XExit where globs := #[`XExit.+]

to `$T/lakefile.lean` (the agents' line: `logs/colab-xexit/lakefile-tail.txt` of either package); this is the
only change to a file of Shea's tree. Then run `apply-overlay.sh` and `build-and-check.sh` as in 11.2a with the
package name. The agents' runs, in each package's `logs/colab/`:

* `u8xm` (saving 1.3868856e-3): one Lake build of 9,305 jobs (`lake.out`), then `checks2.out`: standard axioms,
  `RESULT: PASS`, the statement one unit stronger fails, the controls of `Work.GCert.Data.U8XMControls` proved;
  comparator: `cmp/kcomb-U8XM-pos.log` accepted, `cmp/kcomb-U8XM-neg.log` rejected.
* `u8x` (saving 1.2704023e-3): the 73 modules one at a time (`direct.out`), then `checks.out`: standard axioms,
  `RESULT: PASS`, the stronger statement fails; comparator: `cmp/klower-U8X-pos.log` accepted (its own Lake build
  of 9,289 jobs included), `cmp/klower-U8X-neg.log` rejected.

In both runs the five `XExit` modules had been compiled earlier in the same tree, from the same sources
(`logs/colab-xexit/`), and Lake replayed them. The agents built both packages a second time on the second
Colab machine, and two independent agents checked them; those logs are not in the packages ("not shown"),
except the lines of `community/u8xmhh/logs/colab/setup.out` that show the comparator's verdicts on the second
build of `u8xm`. The main result (11.2) is laid over `u8xm` and needs the same line in the lakefile.

### 11.2d Kothe's circuit with hosting (`community/c588mhh`)

This package, `community/c588mhh`, extends `community/c588m` (11.2a): the same certificate of Kothe's
circuit with mutable sources, with most leftover helper pieces stored in data registers that are idle during the
same run, instead of in banks (111 banks instead of 4,032). This is DaysSky's idle-axis hosting. Its Lean proof
is `Work/XHost`: Mamano's nine-stage network from Shea's tree with hosting added (`HDef`, `HStage`, `HGeom`,
`HCert` copy his `GNetDef`, `GNetStage`, `GGeom`, `GNetCert`; `Frames` adapts his bank lemmas), the kernel
checker `hostSchedK` of the host and bank tables with its soundness lemma `HostSched.ofTab_fitsR`, and the end
theorem `XHost.engine_of_ginv_host`, adapted from his `GX.engine_of_ginv_401`, so the Fourier chain runs on it
unchanged except for its seam (the first file of the chain, the only one that names the end theorem). A control instance, `C588MH`, runs the hosted network with an empty host table
and must give exactly the figure of `community/c588m`. From the root of this repository, with `T` the folder
`lean/complex588` of your copy:

    bash community/scripts/apply-overlay.sh community/c588mhh $T
    bash community/scripts/build-and-check.sh community/c588mhh $T

The first command lays `community/c588m` first (the package's `result.json` names it under "requires"), then the
53 files of this package (49 modules and 4 comparator configurations); every file is new. The second builds the
49 modules in the order of `MODULES.txt` and runs the checks listed in 11.2a, against the statement
`1 - 15699305/10^10` and the one unit stronger, `15699306`.

The agents' run is in `community/c588mhh/logs/colab/`. The modules were compiled with plain `lean`, one at a
time, in rounds while the agents fixed errors (`r2` to `r14`; the earlier rounds hold the failures); the last
compile of every module ended without error (`r14.summary` and the earlier summaries). Then:
`final-C588MHH-compare-pos.log` (`RESULT: PASS`), `final-C588MHH-compare-neg.log` (`RESULT: FAIL` on
`decimalExponentC588MHH`), `r14-Work.FourierC588MHH.Axioms.log` (`propext`, `Classical.choice`, `Quot.sound`),
and `final-C588MH.out` for the control (`RESULT: PASS` at `1 - 13552606/10^10`). An independent agent then
copied the sources afresh, compiled the 32 modules of the package (all but the control) in import order into
its own build folder, and ran `tools/Compare.lean` and the axioms on its build: all modules exit 0,
`RESULT: PASS`, the stronger statement fails, standard axioms (`community/c588mhh/logs/rebuild-skeptic/`). It did
not run the comparator.

Two lemmas proved before the network was built are in `community/c588mhh/supporting/`: L1, which restates what
Shea's tree already proves (`B401.gauged_step`, `B401.gauged_Wn`, `GX.entW_spec`, `GX.entW_form`,
`GX.gauge_step`), and L2, which moves a helper's leftover piece into its host by a relabelling of addresses.
They were compiled on the Mac, and no module of the overlay imports them: the result does the same step with
`host_adapt` and `host_pad` (`Work/XHost/Frames.lean`), adapted from `bank_adapt` and `bank_pad` of Shea's
`Work/Bank401/Adapt.lean`. Two independent agents rebuilt L1 and L2 on the Mac from their own copies, with
the files of Shea's tree they need taken from his tree (`supporting/L1/skeptic/`, `supporting/L2/skeptic/`;
`report.txt` in each is the agent's written verdict).

### 11.3 The official comparator

    chown -R 65534 $T/.lake
    CC=<folder with the comparator's programs> bash community/scripts/run-comparator.sh pos $T comparator/fourier-c588m.json
    CC=... bash community/scripts/run-comparator.sh neg $T comparator/fourier-c588m-neg.json

The first must end with `Your solution is okay!` after both kernels (nanoda and Lean's) accept; the second
must be rejected. `community/base-logs/scripts/cb-tools.sh` builds the comparator, `lean4export`, nanoda and
landrun at the pins of the leaderboard's verifier. The comparator's own sandbox, landrun, needs Landlock, which
the Colab VM's kernel lacks. So the agents ran it inside a stand-in, `community/base-logs/scripts/jail2.sh`: new
mount, network and IPC namespaces, every mount read-only except the tree's `.lake`, the user nobody, no
capabilities, and Unix sockets refused (`nounix.c`). It is not landrun. For the main result the
configurations are `comparator/fourier-u8xmhh.json` and `comparator/fourier-u8xmhh-neg.json`; the agents' runs
are `community/u8xmhh/logs/colab/cmp/fsb-U8XMHH-pos.log` (accepted, 15 minutes, 7.0 GB) and
`fsb-U8XMHH-neg.log` (rejected). Its Lake build found the compiled files up to date (480 compiled files before
and after, `cmp.out`), so it checked them with both kernels without compiling them again; the re-creation of the
leaderboard's verifier and the independent rebuild of 11.2 compiled them from source. For `community/c588mhh`
the configurations are `comparator/fourier-c588mhh.json` and `comparator/fourier-c588mhh-neg.json`; the agents'
runs are `community/c588mhh/logs/colab/final-C588MHH-cmp-pos.out` (accepted, 14 minutes, 8.6 GB; its Lake build
compiled the hosted modules again from source) and `final-C588MHH-cmp-neg.out` (rejected). A first attempt of both
stopped after one second with "Child exited with 1" (end of `final-C588MHH.out`); the cause is not in the log, and
the runs above started a minute later. The control `C588MH` was accepted and its negative rejected
(`cmp-C588MH.out`). For `community/c588m`: `community/c588m/logs/colab/cmp-A1.out` (accepted, 9 minutes) and
`cmp-A2neg.out` (rejected).

After the `chown`, do not run `lake` as root in that tree: git refuses the package clones owned by nobody,
Lake takes that for a changed address, deletes `.lake/packages` and fetches it again without its compiled
files. Copy `.lake/packages` back from the built tree if that happens.

### 11.4 Regenerate the package (no Lean needed)

    bash tools/community/c588m-gen/REPRODUCE.sh

It regenerates the certificate and all 147 Lean and comparator files of `community/c588m` (overlay and
variants) from the inputs in this repository and compares them byte for byte with the published ones. It ends
with `REGENERATION CHECK PASSED` (a few seconds on the Mac). Its steps also re-run the exact Python mirrors of
every kernel check; their outputs are in `community/c588m/preflight/`.

    SHEA=<Shea's clone> bash tools/community/mutable/REPRODUCE.sh

re-derives, from Shea's generated Lean data, his certificate with the source copies, the same circuit with
the copies folded back (mutable sources), the exact checks and controls of that circuit, the annealed bank
table and the two engine files made by text substitution, and compares them with the published ones.
`SHEA=<Shea's clone> VERIFIER=<a clone of bigospeedrun/verifier> bash tools/community/repack/REPRODUCE.sh` does the
same for `community/c588-repack`. Its integer programme needs HiGHS in a `venv/` inside `tools/community/repack`
(`python3 -m venv venv; venv/bin/pip install highspy`, as the script's header says).

### 11.5 Replacing the certificate

    tools/swap-certificate.sh --dry-run <package folder or .tar.gz>
    tools/swap-certificate.sh <package folder or .tar.gz>
    tools/swap-certificate.sh --tree $T --comparator <package>        # on the Linux machine

The script first runs the static checks of the new package as it arrived (`tools/community/package_check.py
check <package> .`), numbered as in its output: 1, the checksums; 2, the statement equal to OpenAI's up to the
exponent and the names, once the modification notice of lines 3 to 12 is removed; 2b, that notice present in
both statement files; 3, the negative statement one unit stronger; 4, the exponent of `Goal.lean`; 5, the
comparator configurations; 6, the module list; 7, the source scan; 8, the figures of `result.json` against its
numerators; 9, every check of `result.json` either "pass" with a log file inside the package, or "not shown" or
"not run". The packages of the Walsh-Hadamard corollaries get checks 1, 3 and 5 to 9 for each of their
theorems. `package_check.py same <file> <file>` compares two statement files with their notices removed.
Then it installs the package under `community/` and keeps the old one (`--remove-old` removes it). With `--tree`
it builds the package in Shea's tree and runs 11.2; with `--comparator` also 11.3; the statuses and log paths go
into the package's `result.json`. It rewrites the package's `SHA256SUMS` and the repository's `MANIFEST.sha256`
and prints the list of checks with their status and log.

It changes the text files (README.md, VERIFY.md, NOTICE, ORIGIN.md, RELATED-WORK.md and `notes/`) and
`community/current.json` only if every check of the new package is "pass" and names a log file inside the
package. Otherwise it changes neither and lists the checks that block it. When it does change them, it replaces
the figures and names of the old package by those of the new one with whole-word matches (`community/c588m` does
not match `community/c588m-gen`), never touches a line about the Walsh-Hadamard corollary (its figures come from
a separate build), and leaves alone the notes that start with a marker comment because they record the figures
of their day (`notes/fifth-result-story.md`, `notes/open-directions-3.md`). It changes no prose: the title, the
credits, the check tables and the sentences that describe how the certificate was made are mine to read and
edit. It never commits. I moved the headline from `community/c588mhh` to `community/u8xmhh` by hand, running the
script only with `--dry-run`: it would also have renamed the many places that name `community/c588mhh` itself,
which stays.

### 11.6 Limits particular to these results

1. **Not this repository's Lean project.** `lake build` at the root of this repository does not build
   `community/`. The packages need Shea's tree, which holds Nil Mamano's gauged extension of the engine and
   the Fourier chain that the packages extend.
2. **No real sandbox.** The comparator ran with the stand-in of 11.3, never with landrun, for every package.
3. **Colab machines for the results.** Every Lean run of the results of this section was on Colab VMs of one
   kind: Linux x86-64, 8 logical processors, 51 GiB, kernel without Landlock, no systemd. The main result,
   `community/u8xmhh`, was built and checked on the second machine, and its independent rebuild ran on the
   first (11.2). `community/c588mhh` and its independent rebuild ran on the first machine, as did most of the
   other packages; the independent rebuild of `community/c588m` and the Walsh-Hadamard corollaries on the
   kernels of `community/c588m` and `community/c588mhh` (11.2b) ran on the second. On the Mac ran only the two
   supporting lemmas of `community/c588mhh` and their two independent rebuilds (11.2d), which no hosting
   result imports; a negative control of `hostSchedK` for `community/c588mhh` (a host table with one remainder
   off by one, rejected) on an earlier version of the data files, whose log is not in the package ("not
   shown"); and the first run of `hostSchedK` on the tables of `community/u8xmhh` and on a copy with one
   remainder off by one (`community/u8xmhh/logs/generation/`). The packages were built in copies of a tree that
   already held Shea's compiled modules (and, for `community/c588mhh`, those of `community/c588m`; for
   `community/u8xmhh`, those of a second build of `community/u8xm`). The exceptions are the re-creations of the
   leaderboard's verifier (point 4 and 11.2), which compiled every module of the packaged submission from
   source (Mathlib came already compiled), and the independent rebuild of `community/u8xmhh`, which compiled
   its 32 modules from source over the agents' first build of `community/u8xm`.
4. **Two later checks of `community/c588m`, both passed, with their limits.** The agents re-created the
   leaderboard's verifier (bigospeedrun.com) offline, on a submission packaged like Nil Mamano's
   (`community/scripts/bos/`). Its own controller needs bwrap and systemd-run, which the Colab VMs lack, so it
   never ran. Instead `community/scripts/bos/mb-bos-run.sh` ran the verifier's three phases with the agents'
   scripts, each inside the stand-in sandbox of 11.3 with `LEAN_NUM_THREADS=1`, and all three passed: the Lean
   check ("Your solution is okay!"), the axiom audit (`propext`, `Classical.choice`, `Quot.sound`) and nanoda
   (both kernels accept). Logs: `community/c588m/logs/colab/bos/`. An independent agent also rebuilt the 70
   modules on a second Colab machine against that machine's own build of Shea's tree, and `tools/Compare.lean`
   passed on its build and failed one unit stronger (`community/c588m/logs/skeptic-rebuild/`); it did not run
   the comparator.
5. **The leaderboard's verifier was not run on the newer packages.** It was never run online. Its three
   phases were re-created offline, as in point 4, for `community/u8xmhh` (11.2) as well; not for
   `community/c588mhh`, `community/u8xm` and `community/u8x`. Independent agents' rebuilds with their logs cover
   `community/u8xmhh` (11.2) and `community/c588mhh` (11.2d); none is shown for the two U8 packages under the
   main result.
6. **Variants not built.** `community/c588m/variants/` holds two other bank tables for the same certificate:
   `cyclic` (1.3513e-3) and `packed` (1.3566e-3). Their exact Python checks pass; Lean has not built them.

## 12. The lever and the sixth result (`package/`, 2026-10-14 to 2026-10-20)

This section is how to check the packages under `package/` yourself: the sixth result, `package/u8x0hq`, and the
recursion lever. Like the results of section 11 they are laid over Ryan Shea's tree; unlike them, they come with
their own scripts in `scripts/` (adapted from `community/scripts/`, with fixes found when they were tested on a
fresh tree on 19 October).

| package | what is proved | saving | status |
|---|---|---|---|
| `package/leverF` (the headline) | the Fourier transform of every length, and convolution: OpenAI's statement with `1 - 17806935/10^10`; the community's #3945 with every head kept, and the recursion lever with the full cut | 1.7807e-3 | checked in Lean, comparator accepted, the leaderboard's verifier re-created offline on the merged form (450 files) and on the unmerged one; rebuilt by an independent agent; this section's scripts run on it from a fresh clone of this repository on 20 October ([scripts/logs-2026-10-19-fresh-clone/](scripts/logs-2026-10-19-fresh-clone/README.md)); the package's `result.json` lists each check with its log |
| `package/leverB` (the half cut on the same table) | the Fourier transform of every length, and convolution: OpenAI's statement with `1 - 17581520/10^10`; the community's #3945 with every head kept, and the recursion lever | 1.7582e-3 | checked in Lean, comparator accepted, the leaderboard's verifier re-created offline on the merged form (449 files) and on the unmerged one; rebuilt by an independent agent; this section's scripts run on it from a fresh clone of this repository on 19 and 20 October ([scripts/logs-2026-10-19-fresh-clone/](scripts/logs-2026-10-19-fresh-clone/README.md)); the package's `result.json` lists each check with its log |
| `package/leverA` | the Fourier transform of every length, and convolution: OpenAI's statement with `1 - 17564707/10^10`; the community's #2979 with every head kept, and the recursion lever | 1.7565e-3 | checked in Lean, comparator accepted, the leaderboard's verifier re-created offline on the merged form (449 files); rebuilt by an independent agent; this section's scripts run on it from a fresh clone of this repository on 19 October ([scripts/logs-2026-10-19-fresh-clone/](scripts/logs-2026-10-19-fresh-clone/README.md)); the package's `result.json` lists each check with its log |
| `package/lever2979` | the Fourier transform of every length, and convolution: OpenAI's statement with `1 - 17556473/10^10`; the community's #2979 as published, with the recursion lever | 1.7556e-3 | checked in Lean, comparator accepted, the leaderboard's verifier re-created offline on the merged form (449 files) and on the unmerged one; rebuilt by an independent agent (the package's `result.json` lists each check with its log) |
| `package/lever` | the Fourier transform of every length, and convolution: OpenAI's statement with `1 - 17501109/10^10`; the community's #3945 as published, with the recursion lever | 1.7501e-3 | checked in Lean, comparator accepted, the leaderboard's verifier re-created offline on the unmerged form (702 files); rebuilt by an independent agent; this section's scripts run on it from a fresh tree on 19 October (the package's `result.json` lists each check with its log) |
| `package/u8x0hq` | the Fourier transform of every length, and convolution: OpenAI's statement with `1 - 16256629/10^10`; the community's U8 word with exit frames, mutable sources without centre copies and same-run hosting | 1.6257e-3 | checked in Lean, comparator accepted, the leaderboard's verifier re-created offline; rebuilt by an independent agent; this section's scripts run on it from a fresh tree on 19 October |

Every lever package requires `package/u8x0hq` (its `result.json`, "requires"), and `scripts/apply-overlay.sh`
lays that first. Each package has `overlay/` (files under `Work/` and `comparator/` only, so the tree's
`lakefile.lean` gains only the two library lines of `package/u8x0hq`), `MODULES.txt` (the new modules in build
order), `result.json` (figures, theorem names, Lean controls and the status of every check with its log),
`README.txt`, `SHA256SUMS`, `logs/`, and, for the lever packages, `controls/` and `generators/`.

### 12.1 What you need

- A Linux machine the size of a Colab VM: the agents used 8 logical processors and about 50 GB of memory. One run
  of `tools/Compare.lean` took about 16 GB and one run of the official comparator on a lever package about 15 GB.
  Allow 20 to 40 GB of disk for one tree with its compiled files (Shea's built tree alone is about 10 GB).
- Lean v4.34.1 through elan, and Mathlib at `d13f23b723b8a846827a245b89c10fc7d3f11612` (Shea's tree pins both).
- Python 3.9 or newer for the scripts (standard library only).

### 12.2 Build Shea's tree

As in 11.1:

    git clone https://github.com/shea256/fourier-transform-below-nlogn shea
    cd shea && git checkout --detach 52394849c7ecf3412d39b925d75fc7fdb98994d8 && cd lean/complex588
    lake exe cache get
    lake build Work.FourierC588.Main Work.FourierC588.UniformFourierChallenge \
      Work.FourierC588Neg.UniformFourierChallenge Work.GCert.Data.C588Controls Work.FourierC588.Axioms

`scripts/base/cb-setup.sh` and `cb-build.sh` are the scripts the agents ran (paths under `/content/work`);
`cb-build.sh` also runs Shea's own Compare control, which must pass before anything is laid over the tree. Below,
`T` is the folder `lean/complex588` of a copy of this built tree. Use one copy per package (`cp -a shea
shea-lever`), so that one package never meets another.

### 12.3 Lay a package over the tree, build it, and check it

From the root of this repository, with `P` one of the folders of the table above:

    bash scripts/apply-overlay.sh package/P $T
    bash scripts/build-and-check.sh package/u8x0hq $T      # optional: the sixth result, 1.6256629e-3
    bash scripts/build-and-check.sh package/P $T

The first command lays `package/u8x0hq` first, then `package/P`. Both are checked against their `SHA256SUMS`; no
file of the tree is overwritten (the script refuses if one would be), and `lakefile.lean` gains the library lines
that each `result.json` lists under `lakefile_add` (`lean_lib XExit` and `lean_lib CSX`, from `package/u8x0hq`).

`build-and-check.sh` builds the modules in the order of `MODULES.txt`, one `lake build` each, waiting while
`MAX_LEAN` Lean processes are running (default 2) and, if `MIN_FREE_MB` is set, while less memory than that is
available. One module at a time, the several hundred modules of a lever package take hours. It is faster to build
them all at once first, which lets Lake compile independent modules side by side, and then run the script, whose
build steps then only confirm that each module is up to date:

    M=$(grep -v '^#' package/P/MODULES.txt); (cd $T && lake build $M)

Then it:

1. prints the axioms of the theorems (the file named `axioms_file` in `result.json`, and a file it writes with
   the same imports for every name under `axiom_targets`);
2. derives both statement files from OpenAI's file (`third-party/openai-math/`) with `scripts/derive_challenges.sh`
   and compares them byte for byte with the package's;
3. runs the Lean controls that `result.json` lists (each must compile, except those marked `"expect": "fails"`,
   which must fail with a Lean error: one head too many, a table chunk one too high);
4. runs `tools/Compare.lean` against the statement (must print `RESULT: PASS`) and against the statement one unit
   stronger (must fail, on the package's `decimalExponent`);
5. scans the package's files for `sorry`, `admit`, `axiom`, `native_decide`, `implemented_by`, `extern`, `unsafe`
   and `opaque` outside comments (the only allowed hits are the two `sorry` lines of each statement file, which
   OpenAI's file has too), and checks that no tracked file of the tree changed except `lakefile.lean`.

It ends with `ALL CHECKS PASSED`. The agents' runs, with machine, times and peak memory, are in each package's
`logs/` and `result.json`. `controls/` of a lever package holds further controls compiled against the built tree:
files that must fail to compile (one head too many, a table chunk one too high), and side checks.

### 12.4 The official comparator

    export CC=<folder with bin/comparator, bin/lean4export, bin/nanoda_bin> JAIL=<folder with jail2.sh, nounix, fresh-bin/lake>
    bash scripts/run-comparator.sh pos $T <the positive configuration of result.json, under comparator/>
    bash scripts/run-comparator.sh neg $T <the negative configuration>

The first must end with `Your solution is okay!` after both kernels accept; the second must be rejected
(`Const does not match ... decimalExponent<suffix>`). `scripts/base/cb-tools.sh` builds the comparator,
`lean4export`, nanoda and landrun at the pins of the leaderboard's verifier. On a lever package the agents'
runs took between about 45 minutes and an hour and a half, and about 15 GB: the kernel phases check every kernel decision of the table of hosts and
of both words again.

The comparator's own sandbox, landrun, needs Landlock, which the Colab VM's kernel lacks. So the agents ran it
inside a stand-in, `scripts/cc-jail/jail2.sh`: new mount, network and IPC namespaces, every mount read-only except
the tree's `.lake` and a private temporary folder, the user nobody, no capabilities, and Unix sockets refused by
the seccomp helper `nounix` (`nounix.c`). It is not landrun. After the `chown`, do not run `lake` as root in that
tree: git refuses the package clones owned by nobody, Lake takes that for a changed address and fetches
`.lake/packages` again without its compiled files. (The script itself runs `chown -R 65534` on the tree's `.lake`.)

### 12.5 The leaderboard's verifier, re-created offline

`scripts/bos/` and `scripts/lever-run/` hold the scripts the agents used (their paths are those of the Colab
machine). They package the proof as a submission folder in the layout of the entries already merged on
[bigospeedrun/verifier](https://github.com/bigospeedrun/verifier) (the import closure of the solution module
copied under `Submission/`, only local import paths changed, and an adapter `Race.bound`), and run the verifier's
three phases at its commit `c92b9b6`: the Lean check, which compiles the packaged proof from source, the audit of
the axioms, and the nanoda kernel. The verifier admits at most 500 files, and a lever package's closure has about 580 to
700 modules, so the packaged proof merges small generated modules into larger ones; the Lean check of the
packaged proof took the agents between about one and two and a half hours (2 h 07 for the headline), and about 10 GB. They run in the same stand-in sandbox as 12.4; the
verifier's own controller is not run.

### 12.6 Re-run the lever's searches and prices (no Lean)

    bash tools/reproduce-search.sh <a new folder> <Shea's tree>/lean/complex588/Work/Bridge401/Gates.lean \
        package/u8x0hq/overlay/Work/XHost

This re-runs every program of `tools/lever/` on the community's files in `third-party/` and on `data/` (the head
test, its control and recount, the prices, the conversion of #3945's words with seven mutants each, the head test in
this project's network, and the T74 regression), and compares each output with `results/`. It ends with
`ALL EQUAL`. It took 14 minutes and 0.97 GB on a Mac (`results/reproduced/`). What each program does, and its
command on its own: [tools/lever/README.md](tools/lever/README.md). The head tests and prices of #2979 and of the
tables with every head kept are in `research/record-track/programs/` and in each package's `fill/checks/`.

To regenerate the Lean data of a lever package from the community's files and compare it with the package:

    bash package/P/generators/REGENERATE.sh <a new folder>

### 12.7 The community's own checks of the circuits (optional)

The circuits' files are in `third-party/` (#3945) and in the generators of `package/lever2979` and
`package/leverA` (#2979; `NOTICE`, section 10.c); the community's check programs are not included as a whole. To re-run them, check out CrocSwap/integer-mult-bounds on a Linux machine with
Python 3.12:

- #3945: pull request #3945 at `91f81d131252db704f945891da1ea54e29d46b8b`, and in
  `research/staged544-terminal-gap-fill/`:

      python3 native/run_native.py inputs/NEW-base.json.gz <out-base>     # and NEW-extra: both end in PASS_COMPLETE_T74_NATIVE_REPLAY
      python3 admission/full_replay.py --work <empty folder> --jobs 3       # PASS_TERMINAL_WIDTH5_COMPLETE_CONDITIONAL_RELEASE

  The first writes the canonical native word that `data/native/` holds ([data/README.md](data/README.md) gives
  the file and content hashes). The second is the community's complete stated verification of #3945; it took the
  agents 1 hour 8 minutes and 6.2 GB, and rejected all 119 of its controls.
- #2979: the main branch at `20ecf289`, `python3 scripts/research_package.py --package research/rank7-joint-endframes`
  (`PASS_REGISTERED_RESEARCH_PACKAGE`) and, in that folder, its complete replay `reproduce.py`
  (`PASS_RANK7_JOINT_ENDFRAMES_COMPLETE_REPLAY`), with the author's recorded versions of numpy, scipy and ortools.
- The tables with every head kept: each package's `fill/` README gives the commands that made the table and the
  community checks that were run on it.

### 12.8 Limits particular to these packages

- Every Lean run of these packages was on Colab machines (Linux). None was on a Mac.
- The comparator and the re-created verifier ran in the stand-in sandbox of 12.4, not in landrun or the
  verifier's own controller, and nothing was submitted to the leaderboard.
- The Lean data of the circuits is generated by programs from the community's files; Lean checks that data, so an
  error in a generator could make the build fail but could not make a false statement pass.
- The times and memory figures are those of the agents' runs on machines that were also running other jobs.

## 13. The programs of `research/` (no Lean)

Each section of [research/](research/README.md) has a `programs/` folder, grouped by finding; each program folder
has a README with what it shows, the command, its time and memory, and what it should print. Python 3.9 or later;
numpy for most programs and scipy for some; sympy, python-flint, Singular, msolve or a C compiler where a README
says so. The heavy ones ran on Colab machines and the light ones on an Apple M4 Mac mini, as each README says.
Two warnings apply everywhere: msolve 0.6.5 reads a product written `x*x` as `x` (write powers as `x^k` and
cross-check with Singular), and `ledger.py` certifies a word only when the number of digit columns is 1 more than
a multiple of 4 ([research/README.md](research/README.md), "Running the programs").

    python3 -I -B research/make_index.py --check      # index.json and open-problems.md agree with the sections

## Sources

For section 11: the U8 word is pull request #869 (Jeremy Kothe) of CrocSwap/integer-mult-bounds, on
his #799, which carries Chafik Boukhalfa's frame changes of #741, and on icekylinx's #645, rebuilt by the agents
at its commit `407d52ff15b2217b6e22bce44728b4c3907b6515` (`community/base-logs/u8-word/`); the community's
owner maps for it and its bit supplier are not used (README.md, section 9, names them); Kothe's circuit is its
pull request #588 and DaysSky's idle-axis hosting its pull request #442 (README.md, section 9, "For the fifth
result", has the other credits); the offline
re-creations of the leaderboard's verifier used bigospeedrun/verifier at `c92b9b6`;
Shea's repository `github.com/shea256/fourier-transform-below-nlogn` at
`52394849c7ecf3412d39b925d75fc7fdb98994d8` (its `lean/` folder is the same at `0e61af6`); Lean v4.34.1;
Mathlib `d13f23b723b8a846827a245b89c10fc7d3f11612`; comparator `d03acab`, lean4export `076e8e5`, nanoda
`3a24072`, landrun `811cfff` (the pins of `bigospeedrun/verifier`, `profiles/tools.json`).

For sections 12 and 13: the community's pull requests #3945 (head `91f81d13`) and #2979 (main
`20ecf289`); Shea's repository at `52394849c7ecf3412d39b925d75fc7fdb98994d8`; bigospeedrun/verifier at `c92b9b6`
(re-created) and `1e009d7` (its file limit); README.md, section 9, "For the lever, the sixth result and the
research", has the rest.
