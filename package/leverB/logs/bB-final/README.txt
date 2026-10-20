logs/bB-final/: the integration run of 19 Oct 2026 on two cloud machines (C3 and C2, 8 cores and 50 GB each)

C3/final.out          the final run: the tree (a copy of the tree in which logs/bB-data/ built every file), the
                      overlay installed identical, lake build up to date, then the checks of C3/l3945h/ (axioms,
                      source greps, the challenge files derived from OpenAI's statement, Compare.lean positive and
                      negative)
C3/L3945H-cmp.out     the official comparator in a stand-in sandbox (uid 65534, no network), positive and negative;
                      the full positive log is C3/l3945h/cmp/
C3/L3945H-bos.out     the leaderboard verifier re-created offline on the merged stage (the import closure has 703
                      local modules, over the pull-request limit of 500; merged to 447): Lean check from source,
                      axiom audit, nanoda, then the bound one grid step higher; C3/bos-L3945H/ has each phase's log
C3/ident.out          the merged stage against the built tree, constant by constant (C3/fast-L3945H/; the two
                      constant dumps are gzipped)
C3/ctl.out            the release controls of ../../controls/ (C3/controls/)
C2/bosU.out           the same verifier on the unmerged 703-module stage, compiled from sources only (the form for
                      the verifier's issue route, which has no file limit)
*-1-*.out             first attempts kept as they ran (prep-1-dupcheckbug.out: a bug in my own duplicate check)
