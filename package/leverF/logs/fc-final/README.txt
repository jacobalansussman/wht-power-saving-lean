logs/fc-final/: the integration run of 19-20 Oct 2026 on three cloud machines (C3, C5 and C2, 8 cores and 50 GB each)

C3/final.out          the final run: the tree (a copy of the tree in which logs/fc-engine/ built every file), the
                      overlay installed identical, lake build up to date, then the checks of C3/l3945f/ (axioms,
                      source greps, the challenge files derived from OpenAI's statement, Compare.lean positive and
                      negative), the side modules (C3/side/) and the official comparator
C3/L3945F-cmp.out     the official comparator in a stand-in sandbox (uid 65534, no network), positive and negative;
                      the full positive log is C3/l3945f/cmp/
C3/ctl.out            the release controls of ../../controls/ (C3/controls/)
C5/bosM.out           the leaderboard verifier re-created offline on the merged stage (the import closure has 707
                      local modules, over the pull-request limit of 500; merged to fewer): Lean check from source,
                      axiom audit, nanoda, then the bound one grid step higher; C5/bos-L3945F*/ has each phase's log
C5/identC5.out        a second build from Ryan Shea's tree with package/u8x0hq and this package's overlay/ alone
                      (C5/l3945fi/), its axioms, then the merged stage against that tree, constant by constant
                      (C5/fast-L3945F/; the two constant dumps are gzipped)
C2/bosU.out           the same verifier on the unmerged 707-module stage, compiled from sources only (the form for
                      the verifier's issue route, which has no file limit)
*/markers.txt         the marker files of each machine and the hashes of every file installed there
*-1-*.out             first attempts kept as they ran (bosM-1-metafiles.out, bosU-1-metafiles.out: the overlay
                      tarball's three top-level list files stopped the install; only Work/ and comparator/ are
                      installed now)
