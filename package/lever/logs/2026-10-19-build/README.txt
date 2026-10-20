Logs of the L3945 build and checks, 19 October 2026, Colab C4 (8 processors, 50 GB), by the integration agent i-final.
  build.out, build-a/, build-b/   lake build of the whole closure in a copy of the U8X0HQ stack: phase A (the light
                                  kernel-data modules, 8 Lean threads, 165 jobs), phase B (everything else, 3 threads,
                                  1085 s, at most 9.6 GB)
  tree/                           the copy of the stack is up to date before the build
  final.out                       the fixed final run (i-skel's isk-final.sh): install, build, axioms, greps,
                                  challenges re-derived from OpenAI's statement, Compare.lean positive / negative
  l3945/                          its step logs: axioms.log, axioms-extra.log, derived-*.lean, compare-pos.log,
                                  compare-neg.log, cmp/ (official comparator positive / negative)
  L3945-cmp.out                   the official comparator, summary
  L3945-bos.out, bos-L3945*/      the bigospeedrun verifier re-created offline (lean, audit, nanoda; negative)
  neg/                            the must-fail negatives (57,313 heads; residual bound 6) and the side check
                                  headBins3945_card = 57,312
  controls/                       the controls of result.json, compiled against the final tree
  controls/Rate3945Sharp.*        a first attempt at the rate control that failed on a wrong object path (the file alone
                                  could not find the module Work.Mixed.Rate3945); Rate3945SharpCtl.* replaced it
  tree-equality.txt               Shea's base + package/u8x0hq + package/lever, laid in that order, equal the files of the
                                  tree that was built and checked, byte for byte, except lakefile.lean's comment lines
