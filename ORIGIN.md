# Where every file comes from

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

No file of openai/math is modified. Two files of this project contain copies of upstream proofs that were then
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
