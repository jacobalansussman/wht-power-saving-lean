# Where every file comes from

This file was written for the first publication (commit `024f763`, 194 Lean modules). Its tables describe the
files of that publication, which are unchanged apart from the documentation files. The section "Files added in the second revision" lists
what the second revision (2026-10-09) added, and the last section, "Files added in the third revision", what
the third (2026-10-09) added.

"openai/math" below means the repository github.com/openai/math at commit
`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb` (Apache-2.0), directory `lean/`.
"This project" means: not from openai/math (see the README for who wrote it and how).

| path | files | origin |
|---|---|---|
| `OAI/Computability/` | 89 Lean files | **openai/math, unmodified.** Byte-identical to `lean/OAI/Computability/...` of the commit above (the "family 130" development: `FourierTransform/` 88 files, `FourierCircuit/Core.lean`). openai/math has 139 files in these two directories; only the 89 that the theorems of this repository import, directly or indirectly, are here. |
| `WHTCheck/Solution.lean` | 1 Lean file | This project. The Walsh-Hadamard statement (`wht`, `WHTProgram`, `WHTGoal`) and a proof of it at the exponent of openai/math (`wht_main`), on top of the OpenAI development. |
| `Work/` | 104 Lean files | This project. Written directly as Lean text (not produced by a script), except the 14 generated files listed below. Some files copy or generalise proofs of openai/math and say so in their header comments (for example `Work/Scratch/Engine.lean`, `Work/Block/*.lean`). |
| `Work/CarrierCheck/ChallengeB2Ke16x.lean`, `ChallengeB2Ke16xR.lean` | 2 of the 104 | This project, with a block from openai/math: lines 22-270 of each file are lines 9-257 of `lean/ComparatorChallenges/UniformFourier.lean` (the RAM cost model), verbatim. They import only Mathlib. |
| `comparator/` | 2 JSON files | This project. Configurations for the `comparator` tool (generated, see below). |
| `tools/Compare.lean` | 1 Lean script | This project. A second, independent check of "same statement, permitted axioms only". Not imported by any module. |
| `tools/certificate/c2-combine-best-h16.json` | 1 JSON file | This project. The carrier certificate (output of a computer search) from which the Lean data module is generated. sha256 `1325bd8a6ccf1bfa2b773d9846c2e7dc8c120c4860c2b4317e109de3e3a20d51`. |
| `tools/gen/`, `tools/regen_check.py`, `tools/ccheck.py` | Python scripts, 2 text templates | This project. Generators of the 14 generated Lean files and of `comparator/*.json`; a reference checker of the certificate that does not use Lean. The template `tools/gen/templates/ChallengeB2L16x.lean.txt` contains the same block from openai/math as the two challenge files (its lines 22-270). |
| `tools/landrun-macos-shim.sh` | 1 shell script | This project. The stand-in for the Linux sandbox `landrun` with which the comparator was run on macOS (see VERIFY.md, Limits). |
| `lakefile.lean`, `lean-toolchain`, `lake-manifest.json` | 3 files | `lean-toolchain` is byte-identical to openai/math `lean/lean-toolchain`. `lakefile.lean` keeps the package name, the package option and the Mathlib revision of openai/math `lean/lakefile.lean` and drops every other dependency. `lake-manifest.json` lists Mathlib and its 8 dependencies at the revisions of openai/math `lean/lake-manifest.json`. |
| `README.md`, `RELATED-WORK.md`, `NOTICE`, `ORIGIN.md`, `VERIFY.md`, `MANIFEST.sha256`, `.gitignore` | 7 files | This project. `RELATED-WORK.md` quotes other people's pull requests and reproduces two BibTeX entries of openai/math. |
| `LICENSE` | 1 file | The Apache License 2.0, byte-identical to `LICENSE` of openai/math (sha256 `c71d239df91726fc519c6eb72d318ec65820627232b2f796219e87dcf35d0ab4`). |

No file under `OAI/` is modified. Since the third revision, modified copies of 11 of them are under
`Work/Fourier/` (last section). Two files of this project contain copies of upstream proofs that were then
changed (`Work/Scratch/Engine.lean`, `Work/Block/Engine.lean`; their comments mark the places, and `NOTICE`
lists them). The namespace `OAI.PowerSaving` and the package name `OAI` are those of
openai/math; they are kept because the new modules extend that development. Their use here does not mean that
OpenAI wrote, checked or endorsed the files of this project.

## How to check the first row yourself

    git clone https://github.com/openai/math openai-math
    git -C openai-math checkout fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb
    for f in $(find OAI -name '*.lean'); do cmp "$f" "openai-math/lean/$f" || echo "DIFFERENT $f"; done
    cmp lean-toolchain openai-math/lean/lean-toolchain
    # the cost model inside the two challenge files (all three lines print the same hash, 39866ff8...87411):
    sed -n 22,270p Work/CarrierCheck/ChallengeB2Ke16x.lean  | shasum -a 256
    sed -n 22,270p Work/CarrierCheck/ChallengeB2Ke16xR.lean | shasum -a 256
    sed -n 9,257p  openai-math/lean/ComparatorChallenges/UniformFourier.lean | shasum -a 256

(`sha256sum` in place of `shasum -a 256` on Linux.) The loop prints nothing when all 89 files are identical.
This was run for this tree against a checkout of that commit: 89 of 89 identical, the three hashes equal.

## Generated files

Made by `tools/gen/ccgen.py` from the certificate (6 files):
`Work/CarrierCheck/Gen/Cr2h16.lean` (the data, 3.9 MB), `Gen/Cr2h16Lab.lean`, `Gen/Cr2h16Shape.lean`,
`Gen/Cr2h16Scal.lean`, `Gen/Cr2h16Hist.lean` (the kernel checks), `Work/CarrierCheck/Cr2h16Cert.lean`.

Made by `tools/gen/crate.py` from `Cr2h16Cert.lean` and the two templates (8 Lean files, 2 JSON files):
`Work/CarrierCheck/RateB2Ke16.lean`, `RateB2Ke16x.lean`, `InstB2Ke16.lean`, `B2Ke16.lean`,
`ChallengeB2Ke16x.lean`, `SolutionB2Ke16x.lean`, `ChallengeB2Ke16xR.lean`, `SolutionB2Ke16xR.lean`,
`comparator/B2Ke16x.json`, `comparator/B2Ke16xR.json`.

`python3 tools/regen_check.py` regenerates all 16 into a temporary directory and compares them byte for byte
with the files of this repository (see VERIFY.md). The generators are not trusted: Lean checks every file.

## Path names inside comments

The Lean files were written in a working tree with another layout, and their comments still use its path names.
Lean files are not edited for the release (their bytes are the ones that were checked), so:

| name in a comment | in this repository |
|---|---|
| `checks/wht12/carrier-check/py/ccgen.py` | `tools/gen/ccgen.py` |
| `checks/wht12/carrier-check/py/crate.py` | `tools/gen/crate.py` |
| `checks/wht9/bridge-rate/py` (`foldlib`, `foldgen`, `b2inst`) | `tools/gen/foldlib.py`, `foldgen.py`, `b2inst.py` (only the definitions that `crate.py` uses) |
| `checks/wht9/bridge-rate/pack.py` (header of the challenge files) | not here; the challenge files are made by `pack` in `tools/gen/crate.py` from `tools/gen/templates/ChallengeB2L16x.lean.txt` |
| `checks/Compare.lean` | `tools/Compare.lean` |
| `Work.BridgeRate.ChallengeB2L16x`, `Work.BridgeRate.SolutionB2L16x` | `tools/gen/templates/*.lean.txt` (text templates of an earlier package, exponent `1 - 4059928/10^10`; not built, not a result of this repository) |
| `Corollaries/...`, `WHTCheck/Challenge.lean`, `Work/Linked/...`, `Work/Fold/...`, other `checks/...` paths, `PLAN.md`, `PLAN2.md` | not in this repository (earlier packages, notes and experiments of the working tree) |
| `(key: ...)`, `(agent key: ...)` | names of the work packages in which a file was written |
| "Scratch file", "feasibility study" (headers of `WHTCheck/Solution.lean`, `Work/Scratch/Engine.lean`) | how these two files started; both are in the import closure of the theorems and are checked with them |
| "upstream" | openai/math |

**The header comment of the two challenge files** (the files to read for the statement) is inherited from an
earlier package and is out of date in three ways. It names `checks/wht9/bridge-rate/pack.py` as the generator;
the files are made by `tools/gen/crate.py`. It describes the file by its differences from
`Corollaries/WHTTwoStage/Challenge.lean`, an earlier statement file of this project with the same text and the
exponent `1 - 66451/10^12`, which is not here. And it does not say that lines 22-270 are the cost model of
openai/math; they are (first table above, and `NOTICE`). The doc comment before the theorem refers to
`WHTCheck/Challenge.lean`, also not here; the statement `WHTGoal` that it means is defined in
`WHTCheck/Solution.lean`. None of these comments affects what Lean checks.

**The usage line of `tools/Compare.lean`** says `checks/Compare.lean`; the file is `tools/Compare.lean` here
(VERIFY.md, section 4, has the commands).

**The field `info` of the certificate** `tools/certificate/c2-combine-best-h16.json` holds notes of the search
that produced it, with names of work packages and files of the working tree (`copt-mech, wht13`,
`copt-design`, `/shared/c2-search-state-5.json`). They point to nothing in this repository, and no script here
reads the field. It was left in place because the sha256 of the file is the one quoted above and in the audits.

## Files added in the second revision

Nothing of the first publication was changed except the documentation files, `MANIFEST.sha256` and `.gitignore` (one more entry).

| path | files | origin |
|---|---|---|
| `Work/GFrame/` | 54 Lean files | This project. Written directly as Lean text. The general frame lemma and the generalised engine. |
| `Work/GCert/Chain/`, `Labels/`, `Scalar/` | 11 + 12 + 15 Lean files | This project. Written directly as Lean text, except five files of `Chain/` (`NetDef`, `NetStage`, `NetCert`, `Stage`, `Geom`), which a script of the working tree made from the files of the first result that their headers name (`Work/Bridge/...`, `Work/BridgeGeom/Net.lean`). That script is not in this repository. The five files are ordinary Lean proofs and are checked like all others. |
| `Work/GCert/Data/` | 39 Lean files | This project. `Raw.lean`, `Count.lean`, `Price.lean`, `Glue.lean` were written directly. The other 35 are generated from the certificate by `tools/gx/gxgen.py` and `tools/gx/gxrate.py`; `python3 tools/gx/regen_check.py` lists them and compares them byte for byte. The 25 files under `Gen/` hold the certificate data and its kernel checks; that data is derived from data files of CrocSwap/integer-mult-bounds (`NOTICE`, section 6). |
| `Work/GCert/Data/ChallengeB2Gp193x.lean` | 1 of the 39 | This project, with the block from openai/math: its lines 22-270 are lines 9-257 of `lean/ComparatorChallenges/UniformFourier.lean`, as in the two challenge files of the first publication (`sed -n 22,270p` of it prints the same sha256, 39866ff8...87411). Its header comment is inherited in the same way as theirs (below) and still names `checks/wht9/bridge-rate/pack.py`; the file is made by `tools/gx/gxrate.py`. |
| `Work/Combine/Split.lean`, `SplitDef.lean` | 2 Lean files | This project. Written directly; imported by the new checker. |
| `comparator/B2Gp193x.json` | 1 JSON file | This project. Generated by `tools/gx/gxrate.py`. |
| `tools/gx/` (15 files), `tools/certificate/gcert1-p11-pr193.json.gz` | Python scripts, 1 certificate | See `tools/gx/ORIGIN.md`. The circuit in the certificate is the design of the authors of CrocSwap/integer-mult-bounds, rebuilt from their published data (`NOTICE`, section 6; README, "Whose circuit this is"). |
| `tools/whatif_copies.py`, `notes/scratch-copies.md` | 1 Python script, 1 note | This project. The what-if of the README, section 7. No part of any proof. |
| `third-party/CrocSwap-integer-mult-bounds/NOTICE` | 1 file | A byte copy of the file `NOTICE` of CrocSwap/integer-mult-bounds at commit `187e1010ac8b259af8e9b5166f68b64bc27b4b47` (the head of its pull request #193); sha256 `c9850418622aea8d1382a230306c940947b5da818f4854b4b272ed60cf239f6c`. Not written by this project. |

Path names in the comments of the new Lean files:

| name in a comment | in this repository |
|---|---|
| `checks/wht26/gx-data/py/gxgen.py`, `gxrate.py` and the other scripts of that directory | `tools/gx/` |
| `checks/wht26/gx-chain/py/netgen.py` | not in this repository (the script behind the five files of `Work/GCert/Chain/` named above) |
| `checks/wht26/shared/gcert-interface.md` | not in this repository (a working specification of the certificate format) |
| `(key: gx-data)`, `(key: gx-chain)`, `(key: gx-labels)`, `(key: gx-scalar)`, `(key: eng-labels)`, `(key: eng-integrate)`, `(key: combine)` | names of the work packages in which a file was written |

## Files added in the third revision

Nothing of the earlier publications was changed except the documentation files and `MANIFEST.sha256`.

| path | files | origin |
|---|---|---|
| `Work/Fourier/Seam.lean` | 1 Lean file | This project. Written directly as Lean text. The kernel program of the generalised engine in the form of OpenAI's theorem `hills_program`, and the facts about the work bound that OpenAI's reduction uses. It repeats the statement of `hills_program` (`TensorProgram.lean`) and the definition of `hills` (`TensorSaving.lean`) with the two names replaced; its header says so. |
| the 11 copied proof files under `Work/Fourier/` (table below) | 11 Lean files | **openai/math, modified.** Made by `tools/fourier/mkchain.py` from the unmodified files under `OAI/Computability/FourierTransform/`: every declaration has the letter `Z` in its name (as a suffix of its own name or of its parent's), `alpha` and `hills` are replaced by `alphaZ` and `hillsZ`, the definitions of the statement by those of `Goal.lean`, and the imports point to the copies. Two further edits were needed, and the first comment of each file names its own: in the copy of `Main` the same constant stands once as a literal number, 2/10^11, and is replaced there too (the two final theorems of that copy also received a comment of their own); and in the copy of `SynchronizedAlgorithm` the name `canopy` is written out as `Grove.canopy` in 5 places, because in the larger environment of the copy the short name would mean another declaration. The script also put the letter on three names inside OpenAI's comments; the comments are otherwise OpenAI's and describe OpenAI's network. Each file says this in its first comment (`NOTICE`, section 7). |
| `Work/Fourier/Goal.lean` | 1 Lean file | **openai/math, modified.** Five definitions of `OAI/Computability/FourierTransform/Goal.lean` with the exponent replaced and the suffix `Z`; made by the same script. |
| `Work/Fourier/UniformFourierChallenge.lean` | 1 Lean file | **openai/math, modified.** `lean/ComparatorChallenges/UniformFourier.lean` with a comment added, the exponent replaced and the suffix `Z` on seven names; made by the same script. Imports only Mathlib. |
| `comparator/UniformFourier.json` | 1 JSON file | This project, on the pattern of `lean/ComparatorChallenges/UniformFourier.json` of openai/math: the same fields, with the module names and theorem names of this repository. |
| `third-party/openai-math/lean/ComparatorChallenges/UniformFourier.lean` | 1 Lean file | **openai/math, unmodified.** Byte-identical to `lean/ComparatorChallenges/UniformFourier.lean` of the commit above (336 lines, sha256 `4c571f275bd7506c3ccd1e6d78f2fca77e07112c969c1c62ad8359d61544ec74`). Not part of the Lean build: `tools/fourier/mkchain.py` reads it to make the challenge, and `diff` shows the changes (VERIFY.md, section 1). With the checkout of "How to check the first row yourself": `cmp third-party/openai-math/lean/ComparatorChallenges/UniformFourier.lean openai-math/lean/ComparatorChallenges/UniformFourier.lean`. |
| `Work/Fourier/Axioms.lean` | 1 Lean file | This project. It proves nothing: when it is built it prints the axioms of the final theorems and of `hillsZ_program` (VERIFY.md, section 3). |
| `tools/fourier/` | 2 Python scripts | This project. `mkchain.py` makes the copies; `regen_check.py` makes the eleven copies, `Goal.lean`, the challenge and the comparator configuration again from the unmodified files under `OAI/` and `third-party/openai-math/` and compares bytes (VERIFY.md, section 6). |

| File of this repository | Derived from (openai/math at fd4aeeb) | Lines: original, copy | Identical lines | Kind of change | Change notice in the file |
|---|---|---|---|---|---|
| `Work/Fourier/ArbitraryLength.lean` | `lean/OAI/Computability/FourierTransform/ArbitraryLength.lean` | 224, 239 | 159 | identifiers renamed on 64 lines (18 names: 15 with the suffix Z, 3 with Z inserted); 1 import line replaced by 1; 15 comment lines added; no other change to the code | yes |
| `Work/Fourier/Asymptotics.lean` | `lean/OAI/Computability/FourierTransform/Asymptotics.lean` | 187, 202 | 112 | identifiers renamed on 74 lines (32 names: 26 with the suffix Z, 6 with Z inserted); 1 import line replaced by 1; 15 comment lines added; no other change to the code | yes |
| `Work/Fourier/Axioms.lean` | none (new file; nearest upstream file `lean/OAI/Computability/FourierCircuit/Core.lean`, similarity 0.06) | -, 13 | - | not a copy | - |
| `Work/Fourier/ConvolutionProgram.lean` | `lean/OAI/Computability/FourierTransform/ConvolutionProgram.lean` | 100, 115 | 77 | identifiers renamed on 22 lines (10 names: 9 with the suffix Z, 1 with Z inserted); 1 import line replaced by 1; 15 comment lines added; no other change to the code | yes |
| `Work/Fourier/Goal.lean` | parts of `lean/OAI/Computability/FourierTransform/Goal.lean` | -, 35 | 6 | written for this repository, with upstream lines: of its 7 code lines, 6 are lines of `lean/OAI/Computability/FourierTransform/Goal.lean` (original lines 35, 75-77, 79, 81): unchanged: its lines 25; with identifiers renamed: its lines 21, 24, 26, 28, 30 | yes |
| `Work/Fourier/Main.lean` | `lean/OAI/Computability/FourierTransform/Main.lean` | 151, 176 | 90 | identifiers renamed on 59 lines (40 names: 35 with the suffix Z, 5 with Z inserted); 1 import line replaced by 2; 24 comment lines added; other code lines changed: 1 of the original, 1 of the copy (copy line 59) | yes |
| `Work/Fourier/Seam.lean` | parts of `lean/OAI/Computability/FourierTransform/Main.lean`, `lean/OAI/Computability/FourierTransform/TensorProgram.lean`, `lean/OAI/Computability/FourierTransform/TensorSaving.lean` | -, 245 | 19 | written for this repository, with upstream lines: of its 122 code lines, 15 are lines of `lean/OAI/Computability/FourierTransform/TensorProgram.lean` (original lines 36-37, 95-96): unchanged: its lines 93-94, 113-114, 141-142, 163-164, 185-187, 194-196; with identifiers renamed: its lines 184; 3 are lines of `lean/OAI/Computability/FourierTransform/TensorSaving.lean` (original lines 83-84, 114): with identifiers renamed: its lines 45-46, 59; 1 is a line of `lean/OAI/Computability/FourierTransform/Main.lean` (original line 60): with identifiers renamed: its lines 213 | yes |
| `Work/Fourier/SectorAlgorithm.lean` | `lean/OAI/Computability/FourierTransform/SectorAlgorithm.lean` | 224, 239 | 162 | identifiers renamed on 61 lines (17 names: 13 with the suffix Z, 4 with Z inserted); 1 import line replaced by 1; 15 comment lines added; no other change to the code | yes |
| `Work/Fourier/SynchronizedAlgorithm.lean` | `lean/OAI/Computability/FourierTransform/SynchronizedAlgorithm.lean` | 126, 142 | 94 | identifiers renamed on 31 lines (12 names: 10 with the suffix Z, 1 with Z inserted, 1 other), the other being `canopy` to `Grove.canopy` (5 lines); 1 import line replaced by 1; 16 comment lines added; no other change to the code | yes |
| `Work/Fourier/TransformProgram.lean` | `lean/OAI/Computability/FourierTransform/TransformProgram.lean` | 146, 161 | 94 | identifiers renamed on 51 lines (21 names: 19 with the suffix Z, 2 with Z inserted); 1 import line replaced by 1; 15 comment lines added; no other change to the code | yes |
| `Work/Fourier/UniformBounds.lean` | `lean/OAI/Computability/FourierTransform/UniformBounds.lean` | 183, 198 | 115 | identifiers renamed on 67 lines (29 names: 25 with the suffix Z, 4 with Z inserted); 1 import line replaced by 1; 15 comment lines added; no other change to the code | yes |
| `Work/Fourier/UniformFourierChallenge.lean` | `lean/ComparatorChallenges/UniformFourier.lean` | 336, 348 | 328 | identifiers renamed on 7 lines (7 names: 7 with the suffix Z); 12 comment lines added; other code lines changed: 1 of the original, 1 of the copy (copy line 281) | yes |
| `Work/Fourier/WorkingCompiler.lean` | `lean/OAI/Computability/FourierTransform/WorkingCompiler.lean` | 153, 168 | 104 | identifiers renamed on 48 lines (14 names: 14 with the suffix Z); 1 import line replaced by 1; 15 comment lines added; no other change to the code | yes |
| `Work/Fourier/WorkingPreparation.lean` | `lean/OAI/Computability/FourierTransform/WorkingPreparation.lean` | 212, 227 | 156 | identifiers renamed on 55 lines (30 names: 28 with the suffix Z, 2 with Z inserted); 1 import line replaced by 1; 15 comment lines added; no other change to the code | yes |
| `Work/Fourier/WorkingTransform.lean` | `lean/OAI/Computability/FourierTransform/WorkingTransform.lean` | 119, 134 | 91 | identifiers renamed on 27 lines (9 names: 8 with the suffix Z, 1 with Z inserted); 1 import line replaced by 1; 15 comment lines added; no other change to the code | yes |

How the table was made: each file was compared line by line with the upstream file (`diff`-style alignment). "Identical lines" are the same bytes. "Identifiers renamed" counts the lines that are equal to the upstream line except for identifiers; "suffix Z" is a name with the letter `Z` appended (`gleamT` to `gleamTZ`), "Z inserted" a name with `Z` inside it or inside one component of a dotted name (`hills_pos` to `hillsZ_pos`, `Breezy.of_moves` to `BreezyZ.of_moves`), and any other substitution is written out. Every line that differs for another reason is counted as an import line, a comment line or an "other code line".

## Files added in the fourth revision (pull request, 2026-10-10)

Contributed by Chafik Boukhalfa (Anthropic Claude assistance). Nothing of the earlier publications was changed
except the documentation files (README.md, VERIFY.md, this file) and `MANIFEST.sha256`.

| path | files | origin |
|---|---|---|
| `tools/certificate/gcert1-p11-pr233-flow.json.gz` | 1 JSON file | The circuit of pull request #233 of CrocSwap/integer-mult-bounds (its author's reuse pairing on the frames of #168 and the word of #193/#194) as one explicit program in the format gcert/1: output of `tools/emit/gcert_emit.py` on the data that the outside repository's scripts write (`tools/emit/ORIGIN.md`). The outside repository's pull request #256 regenerates it and compares bytes. |
| `tools/emit/gcert_emit.py`, `tools/emit/ORIGIN.md` | 1 Python script, 1 text | The contributor's. Reads the outside data as JSON; imports nothing of the outside repository. |
| `Work/GCert/Data/Gen/P233/S0.lean` to `S8.lean`, `Gen/P233.lean`, `Gen/P233Par.lean`, `Gen/P233Tab0-3.lean`, `Gen/P233Seg0-5.lean`, `Gen/P233SegF.lean`, `Gen/P233End.lean`, `Gen/P233Hist.lean`, `Gen/P233YChk.lean`, `Gen/P233Scal0.lean`, `Gen/P233Scal1.lean`, `InstB2Gp233.lean`, `RateB2Gp233.lean`, `RateB2Gp233x.lean`, `ChallengeB2Gp233.lean`, `ChallengeB2Gp233x.lean`, `P233Price.lean`, `P233Cert.lean`, `P233ScalOK.lean`, `B2Gp233.lean`, `B2Gp233Main.lean`, `SolutionB2Gp233.lean`, `SolutionB2Gp233x.lean`; `comparator/B2Gp233.json`, `comparator/B2Gp233x.json` | 34 Lean files, 2 JSON files | **Generated by this project's unchanged generators** `tools/gx/gxgen.py` and `tools/gx/gxrate.py` from the certificate above (`tools/gx/regen_check.py tools/certificate/gcert1-p11-pr233-flow.json.gz P233 B2Gp233` reproduces all 41 generated files). The two challenge files contain the same block of openai/math as `ChallengeB2Gp193x.lean`. |
| `Work/Fourier233/Seam.lean` | 1 Lean file | `Work/Fourier/Seam.lean` with `B2Gp193`/`P193` replaced by `B2Gp233`/`P233`, `7474547`/`7474546` by `7547361`/`7547360`, `alphaZ`/`hillsZ` by `alphaY`/`hillsY` and the paths by those of this directory; nothing else. |
| the 11 copied proof files, `Goal.lean`, `UniformFourierChallenge.lean` under `Work/Fourier233/` | 13 Lean files | **openai/math, modified**, exactly as the files of the third revision under `Work/Fourier/` (table above) with the letter `Y` in place of `Z`, the seam `Work.Fourier233.Seam` and the exponent `1 - 7547360/(10^(10:ℕ))`; made by `tools/fourier/mkchain233.py` (`tools/fourier/regen_check233.py` reproduces them). |
| `Work/Fourier233/Axioms.lean` | 1 Lean file | `Work/Fourier/Axioms.lean` with the names of this chain. Proves nothing. |
| `comparator/UniformFourier233.json` | 1 JSON file | Made by `tools/fourier/mkchain233.py`, as `comparator/UniformFourier.json`. |
| `tools/fourier/mkchain233.py`, `tools/fourier/regen_check233.py` | 2 Python scripts | `mkchain.py` and `regen_check.py` with their constants changed (`diff` them). |
