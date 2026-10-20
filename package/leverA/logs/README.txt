Logs of the L2979H build and checks, 19 October 2026 (Colab machines: 8 processors, 50 GB each).

2026-10-19-build/   Colab C1
  tree/, cp.out                   the copy of the stack (aA-data's lake-built tree) is up to date before the build
  reuse.out                       Lake's outputs of 107 word-stage modules copied from the #2979 lever build's tree,
                                  only for modules whose source is byte-identical; Lake re-checked each trace
  build.out, build-a/, build-b/   lake build of the whole closure: phase A (the host-table data, 41 jobs), phase B
                                  (the engine, the chain and the side checks, 9,511 jobs); the first one-thread attempt,
                                  stopped after 2 minutes, is kept as build-1-onethread.out, build-b-1-onethread/
  final.out                       the final run (isk-final.sh): install, build, axioms, greps, challenges re-derived from
                                  OpenAI's statement, Compare.lean positive / negative
  l2979h/                         its step logs: axioms*.log, sorry/forbidden greps, derived-*.lean, compare-pos.log,
                                  compare-neg.log, cmp/ (official comparator positive / negative)
  L2979H-cmp.out                  the official comparator, summary
  ctl.out, controls/              the controls of result.json, compiled against the final tree
  L2979H-bos.out, bos-L2979H*/    a second copy of the bigospeedrun re-creation (2 Lean threads)
  ident.out, fast-L2979H/         the identity check: every theorem of the built tree against the merged stage
2026-10-19-bos-C4/  Colab C4
  bos4.out, L2979H-bos.out        the bigospeedrun verifier re-created offline on the merged stage (lean, audit, nanoda)
  bos-L2979H/, bos-L2979Hneg/     and its negative at one grid step higher
regenerate/         the Mac: generators/REGENERATE.sh, its comparison table and log
