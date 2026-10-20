package/u8x0hq: the stack underneath, and its own result 1.6256629e-3

What it is
  The Lean files I built between 13 and 14 October 2026 on Ryan Shea's tree, before the lever. They prove
  OpenAI's statement (the Fourier transform of every length, and convolution) with the exponent
  1 - 16256629/10^10, a saving of 1.6256629e-3, on the community's U8 circuit (Jeremy Kothe's #869 on his #799,
  with Chafik Boukhalfa's frame changes of #741 and icekylinx's gauge44 construction of #645). That result was
  never published, because it did not beat the community's own figure. package/lever is laid over this package
  and reuses its engine extensions.

What the 257 files are (paths relative to lean/complex588)
  XExit/ (5), Work/KLow/ (5)      exit frames: a helper may stop at a smaller frame, and the undo step charges only
                                  for the climb it needs
  Work/XMut/ (7), Work/KMX/ (5)   mutable sources: a run may overwrite its inputs if the undo step restores them;
                                  KMX joins this with exit frames
  CSX/ (9)                        the central totals read from the sources themselves, so no source copies remain
  Work/XHost/ (24)                same-run hosting: the nine-stage network of Nil Mamano's files with helper storage
                                  kept in arrays that are idle in the same run, the kernel checker hostSchedK and its
                                  soundness, and the host and bank tables of the instances
  Work/GCert/ (111)               the U8 circuit's certificate data, generated, and its checks
  Work/BlockAccounting/Moment4    the quartic bound term_le4 used by the rate facts
  Work/Fourier*/ (80)             the Fourier chains of five instances on the way: U8XM (1.3869e-3), U8XMHH
                                  (1.6198e-3), U8X0, U8X0HH and U8X0HQ (1.6256629e-3, this package's result), each
                                  with the statement one unit stronger in the folder ending in Neg
  comparator/ (10)                the configurations of the official comparator for those instances
  meta/SOURCES.txt says from which of the agents' packages each file comes; meta/SHA256SUMS.overlay is the hash
  list of the overlay as built.

The tree change
  lakefile.lean gains two library lines (meta/lakefile.lean.ours):
      lean_lib XExit where globs := #[`XExit.+]
      lean_lib CSX where globs := #[`CSX.+]
  scripts/apply-overlay.sh appends them.

Checks (all passed; result.json lists each one with its log)
  14 October 2026, on a Colab machine: lake build (299 s, 8.9 GB at most); the three standard axioms on both
  theorems and on the kernel program; the host-table checker proved, by Lean's kernel, to reject a table with one
  remainder changed; both statement files OpenAI's statement with only the exponent changed
  (16256629 and 16256630) and a suffix added to seven names;
  tools/Compare.lean RESULT: PASS (112 s, 16.0 GB) and, one unit stronger, RESULT: FAIL on decimalExponentU8X0HQ;
  the official comparator in the stand-in sandbox: nanoda and Lean's kernel accept, the stronger statement is
  rejected; the leaderboard's verifier (bigospeedrun/verifier at c92b9b6) re-created offline on the proof packaged
  as a submission: the Lean check compiled it from source (78 minutes, 9.1 GB), the axiom audit and nanoda passed.
  An agent whose only job was to refute it rebuilt every module from a fresh copy of Shea's tree on another
  machine, re-ran the generators (byte-identical; they refuse one unit higher) and did not refute it.
  19 October 2026: rebuilt on all five machines of this round before the lever build, two of them with every
  check above except the comparator (logs/2026-10-19-rebuilds/); and the repository's own scripts were run on
  this package from a fresh copy of Shea's tree: apply-overlay.sh, build-and-check.sh (ALL CHECKS PASSED) and
  run-comparator.sh (accepted; the stronger statement rejected) (logs/2026-10-19-published-scripts/).

Credits
  The circuit is the community's (above). The engine is mine; the network files are adapted from Nil Mamano's,
  as the first lines of each file of Work/XHost say; the entrance gauges and banks are Dugongue's, the layout
  Zhihao Chen's, the bank patterns lydakis's and the starting layout DaysSky's. DaysSky proposed hosting in idle
  arrays and Zhihao Chen wrote the argument for mutable sources, in the community's terms; Jeremy Kothe and
  DreamingOfClouds wrote the community's lemmas on where a helper may end. The Lean proofs of exit frames,
  mutable sources, centre reading and hosting are new; the agents (Claude) wrote them under my direction.

Note
  Comments in seven files (Work/XHost/U8X0HH.lean, U8X0HQ.lean, U8XMHH.lean, Work/XHost/Data/RateU8X0HHx.lean,
  RateU8X0HQx.lean, RateU8XMHHx.lean and XExit/Pack.lean) describe the U8 word as "Kothe #869/#870, Boukhalfa
  #890, Capalbo #942", an earlier and inexact list that names the community's owner maps (#870, #890, #942) among
  its makers. Those owner maps are not used: these files take from the U8 circuit only its additions and frames.
  The credits above are the right ones; the files are left exactly as they were compiled and checked.

Known cosmetic issues (left as compiled, because the files are checked byte for byte): the comment above the
Fourier theorem in each Main.lean calls it a sibling of a Walsh-Hadamard theorem wht_main_block_... that does
not exist (the correct text: it follows from the package's kernel program in Seam.lean); the header of each
Goal.lean cites OpenAI's lines 75-78, 80 and 82 where the definitions are at 75-77, 79 and 81.
