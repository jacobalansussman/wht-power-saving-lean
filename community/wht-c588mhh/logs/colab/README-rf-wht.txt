README-rf-wht (the build of this corollary, 13 Oct 2026). Walsh-Hadamard corollary at the C588MHH kernel figure, built on Colab #2.

Result: OAI.PowerSaving.WHT.wht_main_block_C588MHHx
    ∃ solve W, WHTProgram solve W ∧ WHTTimeBoundsAt (1 - 15699306/(10:ℝ)^10) W
  = saving 1.5699306e-3 for the Walsh-Hadamard transform: one unit above the DFT figure of the same kernel,
  1.5699305e-3, as with C588 (kernel 12331686, DFT 12331685) and C588M (kernel 13552607, DFT 13552606).
  Proof: one term, BlockWHT.wht_main_of_engine (Work/Block/WHT.lean, sha256 e00c1415..., byte-identical to
  wht-power-saving-lean 9c94857) applied to SeamC588MHH.engine_C588MHHx (the hosted build's Work/FourierC588MHH/Seam.lean, community/c588mhh:
  XHost.engine_of_ginv_host with the same-run host table, 111 banks, W 7279, fill parameter s = 40).
  Same method as the corollaries at C588M and C588 (community/wht-c588): gen_wht.py, unchanged (sha256 a9596571...), with
  --tag C588MHH --seam Work.FourierC588MHH.Seam --engine SeamC588MHH.engine_C588MHHx --num 15699306.

Tree: /content/work/rf-wht/shea/lean/complex588 on Colab #2 = a copy (tar, 15 s) of the C588M corollary's tree
  /content/work/r-wht/shea (Shea's tree at 5239484; the c588m overlay, annealed; the C588M modules built up to
  Work.FourierC588M.Seam, 55 of 70; 413 oleans) WITHOUT .lake/packages, which is a symlink to that tree's (read only).
  Tracked changes: 0 before the overlay, after it, and after every run.
  Overlay, new files only (a file that exists with other bytes stops the run): the 17 files of the hosted build the corollary
  needs (16 of Work/XHost + Work/FourierC588MHH/Seam.lean, byte-identical to the hosted build's sources and to the hashes of
  its colab-sources.sha256 (community/c588mhh/data/); list: scripts/hosting-sources.MANIFEST.sha256) and the 8 generated WHT files
  (scripts/wht-files.MANIFEST.sha256).  The Test/ modules and the control instance C588MH were not needed or built.

Build (run.out, build.summary, b-<module>.log/.time): plain `lean` v4.34.1, -DautoImplicit=false (the lakefile's only
  option), LEAN_PATH = the tree's build dir + the Mathlib package build dirs; one process per module, two lanes, each
  start waiting for >= 14 GB MemAvailable.  21 modules, all EXIT 0, no error:
  17 hosting modules (HostK .. C588MHH, FourierC588MHH.Seam; HRoles 1:50, H588MSched 2:03, max RSS 8.6 GB),
  Challenge, Neg Challenge, Main, Solution (about 10 s each, 6.9 GB).
  Attempts: run.out has three.  Attempt 1 built every module of both lanes (all EXIT 0) but stopped because the lane
  function returned the status of its last `[ $ec -ne 0 ]` test (fixed: `return 0`).  Attempt 2 resumed, built the
  remaining 8 modules (all EXIT 0), passed axioms and greps, but Compare.lean could not start its `lean` subprocess
  (toolchain bin not on PATH; fixed).  Attempt 3 resumed with nothing left to build and ran every check.
  The compare-*.log, axioms*.log and sorry-grep*.txt files here are those of attempt 3.  The rfwht-run.sh kept with
  the package is the attempt-3 file (sha256 dc9919fa...); attempts 1 and 2 lacked the two lines that carry a comment
  naming the attempt.

Checks (attempt 3; summary at the end of run.out: ALL CHECKS PASSED):
  * axioms (axioms.log): wht_main_C588MHHx and wht_main_block_C588MHHx on [propext, Classical.choice, Quot.sound];
    (axioms-seam.log): SeamC588MHH.engine_C588MHHx, XHost.engine_of_ginv_host, XHost.HostSched.ofTab_fitsR on the
    same three; the engine's statement is printed there (envelope (1 - 15699306 / 10 ^ 10)).
  * sorry/escape grep (sorry-grep.txt): Main/Solution only a comment mention; the only code `sorry` is the statement
    proof of the two Challenge files.  (sorry-grep-overlay.txt): Work/XHost (without Test/) and the Seam: comment
    mentions ("No `sorry`") only.
  * Compare.lean (/content/work/ours/Compare.lean = published tools/Compare.lean = the tree's tools/Compare.lean,
    sha256 e48d6d3b...) vs Work.WHTC588MHH.Challenge: rc 0, "OK: statement is the same kernel term in both",
    "OK: only permitted axioms", RESULT: PASS (compare-pos.log, 80 s, 15.9 GB); negative control one unit stronger
    (15699307) vs the same solution: rc 1, "FAIL: statement (type) differs", RESULT: FAIL (compare-neg.log, 46 s).

Challenge = OUR published Work/GCert/Data/ChallengeB2Ge8x.lean (9c94857, sha256 a59db9d5...) with exactly:
  header comment (module names + a derivation note), theorem name B2Ge8x -> C588MHHx, exponent 8762479 -> 15699306.
  gen_wht.py asserts that undoing the name and exponent makes every line from `namespace OAI` on byte-identical to the
  published file (3 changed body lines: doc comment exponent, theorem name, statement exponent); a diff on the Mac
  shows the same 3 lines.

Official comparator (cmp.out, cmp/): the C588M corollary's procedure (the base build's cb-cmp.sh with jail2.sh, a stand-in for landrun:
  no Landlock on Colab, no systemd), copies under /content/work/rf-wht/cc with only the paths changed (nounix built
  from the same nounix.c: sha256 1e6528c3..., identical to the C588M corollary's binary).  The tree's .lake was chowned to 65534 first
  (the packages symlink not followed).  One run at a time, each start waiting for >= 16 GB MemAvailable:
    WHTC588MHH-nanoda.json: rc=0 "nanoda kernel accepts the solution", "Lean default kernel accepts the solution",
      "Your solution is okay!"  721 s, 8.8 GB.  Its own `lake build` (as nobody, inside the jail) found no Lake traces for
      the plain-lean oleans and rebuilt the 17 hosting modules, Main and Solution (and the Challenge), all built, nothing
      of the C588M base rebuilt; Lake ran up to three of these compiles at once during its first few minutes (seen in
      samples taken by hand; MemAvailable 33 GB or more in every sample), the only time more than two heavy Lean
      processes of this run were alive together.
    WHTC588MHHNeg.json (15699307): rc=1 "Challenge and solution theorem statement do not match:
      'OAI.PowerSaving.WHT.wht_main_block_C588MHHx'"  283 s.
    WHTC588MHH.json (Lean kernel only): rc=0 "Lean default kernel accepts the solution", "Your solution is okay!"
      327 s, 6.8 GB; nothing rebuilt.
  434 oleans before and after the second and third runs; git: 65 status lines (untracked overlay only) before and after
  every run, 0 tracked changes.  This is the sandbox stand-in, not the official bigospeedrun verifier, which has not
  been run on this result.
  NOTE: the tree's .lake is now owned by nobody.  Do NOT run `lake` as root in /content/work/rf-wht/shea.

Files here: run.out (all three attempts), build.summary, b-<module>.log/.time, laneA.out/laneB.out (those of the last
  attempt; attempt 1's lane output is in run.out), axioms.log/.time, axioms-seam.log/.time, sorry-grep.txt,
  sorry-grep-overlay.txt, compare-pos/neg.log/.time, cmp.out, cmp/.  Scripts and generated files:
  the agents' working folder (scripts/, ovl/, src-C588MHH/) and this package (scripts-as-run/).
