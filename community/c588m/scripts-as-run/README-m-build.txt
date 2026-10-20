README-m-build (m-build, 13 Oct 2026). Copy-free #588 (C588M, annealed bank table) on top of c-base's built tree.

Tree: /content/work/mb/shea/lean/complex588  (cp -a of /content/work/shea at 5239484; m-gen overlay "annealed"
      applied = 72 new untracked files; no tracked file changed). ALL 70 new modules are built.
Result there (logs /content/work/logs/m-build/):
  * transform_mainC588M, convolution_mainC588M: axioms [propext, Classical.choice, Quot.sound] (axioms.log)
  * tools/Compare.lean (= our Compare.lean, same sha256) vs OpenAI's statement at 1 - 13552606/10^10: RESULT: PASS
    negative one unit stronger (13552607): RESULT: FAIL, 'decimalExponentC588M: definition body differs'
  * official comparator (c-base's cb-cmp.sh procedure, jail2.sh stand-in): nanoda + Lean kernel accept,
    'Your solution is okay!' (cmp-A1.out); negative config rejected (cmp-A2neg.out)

GOTCHA (cost me the tree's packages once, repaired): after `chown -R 65534 .lake` for the comparator, do NOT run
`lake ...` as root in that tree: git refuses the nobody-owned package repos ("dubious ownership"), Lake reads that as
"URL has changed", DELETES .lake/packages/* and re-clones them from GitHub without build outputs. Repair used:
  rm -rf .lake/packages && cp -a /content/work/shea/lean/complex588/.lake/packages .lake/packages
Current ownership: .lake/build nobody, .lake/packages root (root lake works; for another comparator run as
nobody, chown -R 65534 .lake again first, and after that no root lake).

Offline bigospeedrun/verifier re-creation for C588M: stage /content/work/mb/bos (packaged like nmamano's
submissions/dft/nmamano-shea-c588; bound 6776303/5000000000), logs /content/work/logs/m-build/bos/.
