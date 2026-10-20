The scripts the agents ran for package/lever on 19 October 2026, as run (paths are those of the Colab machine).

  isk-lib.sh          shared helpers: paths, the memory-aware wait (at most 3 large Lean processes), timed steps
  isk-tree.sh         a private copy of the built stack of package/u8x0hq, checked up to date
  isk-install.sh      lays new files into that copy; refuses any file that would change an existing one
  if-build.sh         the build: the light kernel-data modules first (8 Lean threads), then everything else
  isk-final.sh        the final run: install, build, then isk-checks.sh, then the comparator and the verifier
                      re-creation side by side
  isk-checks.sh       axioms of the final theorems and of the engine, the seam and the rate; the source scan; both
                      statement files re-derived from OpenAI's and compared byte for byte; tools/Compare.lean
                      positive (must pass) and one unit stronger (must fail)
  isk-cmp.sh          the official comparator, positive and negative, in the stand-in sandbox of scripts/cc-jail/
  isk-bos.sh          the leaderboard's verifier re-created offline: packages the proof as a submission folder
  isk_submission.py   (after checking that the same packaging of Shea's tree reproduces Nil Mamano's merged entry
                      byte for byte), then the verifier's three phases, and the bound one grid step higher
  bos-comment-L3945.txt  the comment block of the packaged Submission.lean
  isk-build.sh, isk-control.sh, isk-hyp.sh, isk-skel.sh   helpers of the same set used in the trial runs
  if-srcpack.sh       the tarball of every source file of the final tree, with its hash list
  if-pubtest.sh       the run of scripts/apply-overlay.sh and scripts/build-and-check.sh on package/lever, on a
                      fresh copy of Shea's tree

They do the same steps as scripts/build-and-check.sh and scripts/run-comparator.sh, with the machine's paths.
