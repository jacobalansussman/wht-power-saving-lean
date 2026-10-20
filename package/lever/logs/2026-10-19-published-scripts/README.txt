The repository's published scripts on package/lever, 19 October 2026, Colab C3 (8 processors, 50 GB), on a fresh
copy of Ryan Shea's built tree at 52394849 (scripts/lever-run/if-pubtest.sh).

  pubtest.out            the whole run: the copy of the tree, apply-overlay.sh (it laid package/u8x0hq first, as
                         result.json "requires", then package/lever), the laid tree against the hash list of the tree
                         checked on C4 (equal except lakefile.lean's comment lines: laid-vs-checked.diff), one parallel
                         `lake build` of every module (prebuild.log, prebuild.time), then build-and-check.sh
  pubtest/apply/         apply-overlay.sh's lists of the files it laid
  pubtest/check/         build-and-check.sh's step logs: one per module (each confirms the module is up to date),
                         axioms, the statement files derived from OpenAI's, the two Lean controls, Compare positive
                         and negative, the source scan
  pubtest/eval-fixed.out the script's final evaluation re-run on the same logs (bac_eval.py is that block of
                         build-and-check.sh, unchanged), with package/lever/result.json's corrected tree_change

On the first run every Lean check passed, and one rule failed: lakefile.lean carried the two library lines that
package/u8x0hq adds, which package/lever's result.json did not yet quote under tree_change. The result.json in
this repository quotes them; with it, the evaluation ends with ALL CHECKS PASSED.
