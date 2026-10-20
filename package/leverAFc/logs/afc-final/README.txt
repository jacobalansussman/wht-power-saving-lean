logs/afc-final/: the integration run of 20 Oct 2026 on three cloud machines (C4, C5 and C2, 8 cores and 50 GB each)

C4/final.out          the final run: the tree (a copy of the tree in which logs/afc-data/ and logs/afc-cut/ built every
                      file), the overlay and the stage installed identical, lake build up to date, then the checks of
                      C4/l2979fc/ (axioms, source greps, the challenge files derived from OpenAI's statement,
                      Compare.lean positive and negative), the side modules (C4/side/) and the official comparator
C4/L2979Fc-cmp.out    the official comparator in a stand-in sandbox (uid 65534, no network), positive and negative;
                      the full positive log is C4/l2979fc/cmp/
C4/ctl.out            the release controls of ../../controls/ (C4/controls/)
C5/bosM.out           the leaderboard verifier re-created offline on the merged stage (the import closure has 577
                      local modules, over the pull-request limit of 500; merged to fewer): Lean check from source,
                      axiom audit, nanoda, then the bound one grid step higher; C5/bos-L2979Fc*/ has each phase's log
C5/identC5.out        a second build from Ryan Shea's tree with package/u8x0hq and this package's overlay/ alone
                      (C5/l2979fci/), its axioms, then the merged stage against that tree, constant by constant
                      (C5/fast-L2979Fc/; the two constant dumps are gzipped)
C2/bosU.out           the same verifier on the unmerged 577-module stage, compiled from sources only (the form for
                      the verifier's issue route, which has no file limit), and a grep of every packaged file
                      (C2/closure-grep/)
*/markers.txt         the marker files of each machine and the hashes of every file installed there
