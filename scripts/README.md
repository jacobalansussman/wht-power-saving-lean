# Scripts

| script | what | used in |
|---|---|---|
| `apply-overlay.sh PACKAGE TREE` | lays a package (and the packages it requires) over Shea's `lean/complex588`, after checking its SHA256SUMS; refuses to overwrite a file of the tree; appends the `lakefile_add` lines of its `result.json` | VERIFY.md 12.3 |
| `build-and-check.sh PACKAGE TREE` | builds the package's modules and runs its checks: axioms, the statement files re-derived from OpenAI's, the Lean controls (each compiles, or fails where `result.json` says it must), `tools/Compare.lean` both ways, the source scan; ends in `ALL CHECKS PASSED` | VERIFY.md 12.3 |
| `run-comparator.sh LABEL TREE CONFIG` | one run of the official comparator, in the stand-in sandbox | VERIFY.md 12.4 |
| `derive_challenges.sh REF NUM SUFFIX` | OpenAI's statement with the exponent `1 - NUM/10^10` and the suffix on its seven names | build-and-check.sh |
| `base/` | the scripts the agents ran to set up a Colab machine: elan, Shea's tree and Mathlib's cache (`cb-setup.sh`), its build with his Compare control (`cb-build.sh`), the comparator tools at the verifier's pins (`cb-tools.sh`, `cb-cmp.sh`); paths under `/content/work` | VERIFY.md 12.2, 12.4 |
| `cc-jail/` | the stand-in sandbox: `jail2.sh` and `nounix.c` (a seccomp helper that refuses Unix sockets), as used | VERIFY.md 12.4 |
| `logs-2026-10-19-fresh-clone/` | the run of `apply-overlay.sh`, `build-and-check.sh` and `run-comparator.sh` from a fresh clone of this repository on a fresh copy of Shea's tree (19 October), with its logs | VERIFY.md 12 |
| `bos/` | the offline re-creation of the leaderboard's verifier (`fspb-bos-prep.sh`, `fspb-bos-run.sh`, `fspb_submission.py`) and the comparator recipe (`fspb-cmp.sh`), as run for the earlier result; their paths are those of the Colab machine | VERIFY.md 12.4, 12.5 |

`apply-overlay.sh`, `build-and-check.sh` and `run-comparator.sh` are adapted from the scripts of the same names
that I wrote for the packages of 13 October (`community/scripts/`); the first lines of each say what changed. The
scripts in `base/`, `cc-jail/` and `bos/` are kept exactly as they were run.
