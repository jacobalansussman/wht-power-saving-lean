README-r-wht (r-wht, 13 Oct 2026). Walsh-Hadamard corollary at the C588M kernel figure, built on Colab #2.

Result: OAI.PowerSaving.WHT.wht_main_block_C588Mx
    ∃ solve W, WHTProgram solve W ∧ WHTTimeBoundsAt (1 - 13552607/(10:ℝ)^10) W
  = saving 1.3552607e-3 for the Walsh-Hadamard transform (one unit above the DFT figure 1.3552606e-3, as with
  C588: kernel 12331686, DFT 12331685).  Proof: one term, BlockWHT.wht_main_of_engine (our Work/Block/WHT.lean,
  sha256 e00c1415..., byte-identical to wht-power-saving-lean 9c94857) applied to SeamC588M.engine_C588Mx
  (m-gen package c588m-colab-pkg.tar.gz sha256 6554fd38..., variant annealed).

Tree: /content/work/r-wht/shea/lean/complex588 on Colab #2 (cp -a of /content/work/shea at 5239484; tracked
  changes 0 before and after; overlay = 72 m-gen files + 8 WHT files, all new).  Only the C588M modules up to
  Work.FourierC588M.Seam (55 of 70) were built: the WHT corollary needs nothing after the Seam.
  NOTE: .lake is now owned by nobody (for the comparator). Do NOT run `lake` as root in this tree (m-build's
  gotcha: Lake deletes and re-clones .lake/packages).

Checks (run.out, summary at its end: ALL CHECKS PASSED):
  * build: 55 C588M modules + Challenge, Neg Challenge, Main, Solution, one module per `lake build`, all rc=0;
    `lake build --no-build` rc=0 (NN-*.log/.time, wht-build-*.log/.time, wht-nobuild.log)
  * axioms (axioms.log): both theorems on [propext, Classical.choice, Quot.sound], no sorryAx
  * sorry/escape grep (sorry-grep.txt): Main/Solution only a comment mention; the only code `sorry` is the
    statement proof in the two Challenge files
  * Compare.lean (= published tools/Compare.lean, sha256 e48d6d3b...) vs Work.WHTC588M.Challenge: RESULT: PASS
    (compare-pos.log); negative control one unit stronger (13552608): rc=1, "FAIL: statement (type) differs",
    RESULT: FAIL (compare-neg.log)
  * official comparator (cmp.out, cmp/), c-base's cb-cmp.sh procedure with jail2.sh (stand-in for landrun, no
    Landlock on Colab), private copies under /content/work/r-wht/cc (only the paths changed; jail probed first:
    runs as nobody, read-only outside .lake and /tmp, AF_UNIX refused, no network):
      WHTC588M.json (Lean kernel): rc=0 "Lean default kernel accepts the solution / Your solution is okay!" 288 s, 7.0 GB
      WHTC588M-nanoda.json: rc=0 "nanoda kernel accepts ... Lean default kernel accepts ... okay!" 305 s, 7.0 GB
      WHTC588MNeg.json (13552608): rc=1 "Challenge and solution theorem statement do not match" 61 s
    413 oleans before and after each run, no module rebuilt; config sha256 = MANIFEST.sha256.

Challenge = OUR published Work/GCert/Data/ChallengeB2Ge8x.lean (9c94857, sha256 a59db9d5...) with exactly:
  header comment (module names + a derivation note), theorem name B2Ge8x -> C588Mx, exponent 8762479 -> 13552607.
  gen_wht.py asserts that undoing the name and exponent makes every line from `namespace OAI` on byte-identical
  to the published file (3 changed body lines: doc comment exponent, theorem name, statement exponent).

Files here: src/ (the 8 generated files + MANIFEST.sha256 + PARAMS.json), scripts/ (gen_wht.py, kernel-note-C588M.txt,
  rwht-run.sh, rwht-swap.sh, cc/), applied.txt (m-gen overlay record), cp.log.

SWAP (for the ~1.39e-3 or ~1.5785e-3 certificate): scripts/rwht-swap.sh <TAG> <Seam module> <engine theorem>
  <NUM> <kernel-note.txt> <remote package dir> <variant>  -- regenerates the files (gen_wht.py), copies them,
  makes a fresh tree copy and runs rwht-run.sh (build to the Seam, WHT build, all checks) in tmux on Colab #2.
  Sanity: gen_wht.py with C588 / SeamC588.engine_C588x / 12331686 reproduces c-wht's Challenge, Neg, Solution,
  Axioms and comparator configs byte for byte (apart from c-wht's agent label in the derivation note).
