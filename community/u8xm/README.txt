community/u8xm: the community's U8 word with exit frames and mutable sources
============================================================================

What it proves. The Fourier transform of every length, and convolution, with z = 1 - 13868856/10^10 (saving
1.3868856e-3), on OpenAI's statement with only the exponent and seven names changed (suffix U8XM).

What it is. The U8 word of CrocSwap/integer-mult-bounds (Jeremy Kothe's #869, on his #799, which carries Chafik
Boukhalfa's frame changes of #741, and on icekylinx's gauge44 construction #645; README.md, section 9, gives the
rest of the line and what the community's figures for this word use that this package does not) in Ryan Shea's
Lean tree, with two changes: 126 helpers end at exit frames of dimension 4 or 8 instead of the full frame (as in
community/u8x), and Shea's source copies are folded back into their sources (mutable sources, as in
community/c588m), except the 5 copies of the sources that host a split centre. Work/KMX is the new Lean that
combines the two. Its bank table has 3,771 banks.

How to lay it. On a built copy T of Shea's lean/complex588 at 52394849c7ecf3412d39b925d75fc7fdb98994d8:
add the line   lean_lib XExit where globs := #[`XExit.+]   to T/lakefile.lean (the five XExit modules are a
library of their own; logs/colab-xexit/lakefile-tail.txt), then

    bash community/scripts/apply-overlay.sh community/u8xm $T
    bash community/scripts/build-and-check.sh community/u8xm $T

Folders. overlay/ and MODULES.txt (every module of the overlay, in import order); logs/colab/ the agents' run on
the first Colab machine (lake.out: the build; checks2.out: axioms, Compare.lean positive and negative, controls;
cmp/: the official comparator, positive accepted and negative rejected; run.sh, run2.sh: the scripts of the
run); logs/colab-xexit/ the earlier build of the five XExit modules in the same tree; logs/agent-notes/ the
agents' build notes and the module list as run; scripts-as-run/ the package scripts as they were on the Colab;
logs/static-check.txt the static checks on the Mac. result.json lists every check with its log; the checks
without a log are marked "not shown" or "not run".

Credit comments, left as they were compiled or written. Line 9 of overlay/XExit/Pack.lean names the community's
line as "#799/#869/#890"; #890 (Chafik Boukhalfa's owner packing) leaves the word unchanged, and this package does
not use it. The agents' notes in logs/agent-notes/BUILD.md credit the word to "#869/#870, Chafik Boukhalfa #890,
Capalbo #942" and the owners to Hao Yan; the word is Kothe's #869, and the community's owner maps (#870 with
Hao Yan's #858, #890, #942) are not used. The right credits are in result.json and in README.md.

Limits. The comparator ran inside the stand-in sandbox jail2.sh, never with landrun. The leaderboard's own
verifier has not been run on this package. No independent rebuild of this package is shown.
