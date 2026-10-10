# How to rebuild this repository and re-run the checks

This file says how to rebuild the proofs and re-run the checks yourself. To check the newest result, the
fourth, start at section 10. It sends you to sections 2 to 5 for what you need and for the first commands.

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

`MANIFEST.sha256` lists every file of the repository except itself.

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
   certificate at the same time and may need more than 16 GB.
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
    of the third). One run was still going when this text was written: the comparator with
    `comparator/UniformFourier.json` in a tree where nothing had been compiled before. It is the third attempt
    in that tree, and each attempt continued the build of the one before. It has no verdict here.
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
171 modules replayed, none built. The sandbox there was a stand-in built on Linux namespaces, not `landrun`,
because that kernel has no Landlock; the stand-in was not tested on its own in these runs. This is not a run
with the real sandbox. The times are from a machine under load.

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
