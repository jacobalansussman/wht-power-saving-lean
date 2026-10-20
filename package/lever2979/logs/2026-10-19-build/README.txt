Logs of the L2979 build and checks, 19 October 2026, Colab C1 (8 processors, 50 GB), by the integration agent x-final.
  build.out, build-a/, build-b/   lake build of the whole closure in a copy of the U8X0HQ stack: phase A (the host
                                  table's kernel data, 8 Lean threads, 41 jobs), phase B (everything else, 3 threads,
                                  9,511 jobs, 1,204 s, at most 9.5 GB)
  tree/                           the copy of the stack is up to date before the build
  final.out                       the fixed final run (isk-final.sh): install, build, axioms, greps, challenges
                                  re-derived from OpenAI's statement, Compare.lean positive / negative
  l2979/                          its step logs: axioms.log, axioms-extra.log, derived-*.lean, compare-pos.log,
                                  compare-neg.log, the grep results, cmp/ (official comparator positive / negative)
  L2979-cmp.out                   the official comparator, summary
  L2979-bos.out, bos-L2979/,      the bigospeedrun verifier re-created offline on the MERGED stage (449 files, under
  bos-L2979neg/                   the 500-file limit): lean, audit, nanoda; negative one grid step higher
  ident.out, fast-L2979/          the identity check of the merge: the merged stage and the built tree hold the same
                                  constants and the same proof terms (the two constant dumps, 8 MB, are not included)
The unmerged re-creation (Colab C5) is in ../2026-10-19-bos-unmerged/, the controls in ../controls/, the replay of the
generators (on a laptop) in ../regenerate/.
