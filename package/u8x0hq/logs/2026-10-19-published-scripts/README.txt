The published scripts of this repository, run on this package on 19 October 2026 (Colab, 8 processors, about
50 GB), on a fresh copy of Ryan Shea's built tree at 5239484:

  test.out, apply/   scripts/apply-overlay.sh package/u8x0hq: the 257 files laid, the two library lines appended
                     to lakefile.lean; then scripts/build-and-check.sh, which built the U8X0HQ targets from source
                     (Work.XHost.U8X0HQ and its dependencies in 734 s; build logs in check/). Every Lean check of that
                     run passed (axioms, statement files, both Lean controls, Compare PASS and one unit stronger
                     FAIL); two of its source-scan rules were too narrow for this package and failed.
  check3/, test3.out the same build-and-check.sh after the fixes (the scan rules, and resolving its own folder
                     before changing directory), started with relative paths in the built tree: ALL CHECKS PASSED.
  cmp/, cmp-*.out    scripts/run-comparator.sh in the stand-in sandbox: positive accepted by nanoda and by Lean's
                     kernel ("Your solution is okay!", 895 s, 7.0 GB at most); one unit stronger rejected ("Const
                     does not match ... decimalExponentU8X0HQ").
