community/u8xmhh: the community's U8 word with exit frames, mutable sources and same-run hosting
==============================================================================================

What it proves. The Fourier transform of every length, and convolution, in time O(n (log n)^z) with
z = 1 - 16197705/10^10, a saving of 1.6197705e-3. The statement is OpenAI's UniformFourier.lean with only the
exponent and seven names changed (suffix U8XMHH). The theorems are OAI.PowerSaving.transform_mainU8XMHH and
OAI.PowerSaving.convolution_mainU8XMHH.

What it is. Two packages of this folder put together, each unchanged:
  * the circuit and its certificate are those of community/u8xm: the U8 word of CrocSwap/integer-mult-bounds
    (Jeremy Kothe's #869, on his #799, which carries Chafik Boukhalfa's frame changes of #741, and on
    icekylinx's gauge44 construction #645; README.md, section 9, gives the rest of the line and what the
    community's figures for this word use that this package does not) in Ryan Shea's Lean tree, with 126 helpers
    ending at exit frames and with Shea's source copies folded back into their sources (mutable sources), except
    the 5 copies of the sources that host a split centre. In community/u8xm every leftover helper piece is
    stored in a bank: 3,771 banks;
  * the hosted network of community/c588mhh (Work/XHost): most leftover pieces are stored instead in data
    registers that sit idle during the same run of the network, an idea DaysSky proposed for the community's
    circuits (pull request #442 of CrocSwap/integer-mult-bounds). community/c588mhh/README.txt says what each
    module of the network proves and which files of Shea's tree (from Nil Mamano's submission) it copies or
    adapts.
What is new here is the instance: a host table for the U8 circuit (2,016 hosts, 4,092 hosted pieces in each of the
8 runs, the same table in every run), the 310 banks that remain (data/), the Lean data modules that state them
and are checked by the kernel checker hostSchedK (Work/XHost/Data), the glue Work/XHost/U8XMHH.lean, and the
Fourier chain Work/FourierU8XMHH with its seam pointed at XHost.engine_of_ginv_host. The twelve network modules
in overlay/Work/XHost are byte-identical to those of community/c588mhh; they are repeated here so that this
package needs only community/u8xm under it.

Wrong comments, left as they were compiled.
  * Comments in Work/FourierU8XMHH/Seam.lean (its header and the comment on its main theorem), carried over by
    the generator from the seam of community/u8xm, name GX.engine_of_ginv_401, the unit list blocksU8XMx, the
    rate fact fold_U8XMx, the scalar half GXD.U8XM.scalOK and a file Work/GCert/Data/U8XMHH.lean. The code
    calls XHost.engine_of_ginv_host and fold_U8XMHHx, and the instance is in Work/XHost/U8XMHH.lean; the
    independent rebuild's Uses.lean (logs/skeptic-rebuild/logs/uses.log) shows that the proofs reach
    fold_U8XMHHx and not fold_U8XMx.
  * The header of Work/XHost/HEnd.lean is the one community/c588mhh/README.txt describes as wrong in part.
  * The credit for the U8 word in two comments is not exact: Work/XHost/U8XMHH.lean, line 10 ("the U8 word of
    Kothe #869/#870, Boukhalfa #890, Capalbo #942") and Work/XHost/Data/RateU8XMHHx.lean, line 5 ("72 copies
    of the U8 word (#869/#870/#890/#942)"), both written by generators/lean-gen/gen_hu8xm.py, which says the
    same in its lines 2, 124 and 161. The word is Kothe's #869; #870, #890 and #942 add owner maps for it,
    which this package does not use, and leave the word unchanged. The same earlier list is in the scripts as
    run (scripts-as-run/REPRODUCE-build.sh; scripts-as-run/build/fsb_submission.py and the file it wrote for the
    verifier re-creation, logs/colab/bos-H-stage/Submission.lean, whose header also credits Hao Yan for owners)
    and in a key of data/words.json. The right credits are in result.json and in README.md.

How to lay it. It needs community/u8xm first (result.json: "requires"), and community/u8xm needs one line added to
the tree's lakefile.lean (community/u8xm/README.txt). From the root of this repository, with T the folder
lean/complex588 of a built copy of Shea's tree at 52394849c7ecf3412d39b925d75fc7fdb98994d8:

    echo 'lean_lib XExit where globs := #[`XExit.+]' >> $T/lakefile.lean
    bash community/scripts/apply-overlay.sh community/u8xmhh $T     # applies community/u8xm, then this package
    bash community/scripts/build-and-check.sh community/u8xmhh $T

The package adds 34 files (32 modules in MODULES.txt and 2 comparator configurations); no file of the tree
changes except lakefile.lean.

Folders.
  overlay/           the files laid over lean/complex588
  data/              the host table of one run (HU8XM-hosttab.json), the banks (hu8xm_sr.padsa.json, with the slot
                     list .slots.gz that audit_slots.py reads), the price summary of the generator (U8XMHH-price.json:
                     W 7478, 310 banks, kernel 16197706, DFT 16197705), the prices of every variant tried that day
                     with the controls (PRICES.json, controls-price.json), the inputs of the generators (words.json,
                     U8XM-price.json, u8xmhh-profile.json, seam_U8XMHH.json), the provenance of the Fourier files,
                     and the checksums of the two archives uploaded to the Colab machine and of their files as built
                     (colab-sources.sha256)
  generators/        the programs that chose the host table and wrote the Lean files (generators/README.txt)
  logs/colab/        the build on the second Colab machine, unedited: setup.out (the copy of the u8xm tree and the
                     hash checks), build.out (the six staged builds b1 to b6, axioms, source scan, challenge
                     derivation, the host-table negative control, Compare.lean), the logs of each step, cmp/ and
                     cmp.out (the official comparator), bos-H.out, bos-H/ and bos-H-stage/ (the offline re-creation
                     of the leaderboard's verifier), ax/ (the axiom check and the negative-control table). The two
                     statement files the build derived again from OpenAI's file (derived-*.lean) are left out: they
                     equal the as-built statement files.
  logs/mac-preflight/  the Python checks on the Mac before the build
  logs/generation/   the logs of the pricing and of the generators on the Mac, and a first run of hostSchedK on the
                     Mac on the generated tables (true) and on a copy with one remainder off by one (false)
  logs/skeptic-statement/  the independent check of the statement and its replay of the comparator (report.txt,
                     logs/, scripts/)
  logs/skeptic-rebuild/  the independent rebuild on the first Colab machine (report.txt, RESULT.json, logs/: one
                     log per module, build.out, compare-pos.log, compare-neg.log, uses.log; scripts/)
  logs/static-check.txt  tools/community/package_check.py check, run on the Mac
  scripts-as-run/    build/: the Colab scripts as run (fsb-setup.sh, fsb-build.sh, fsb-cmp.sh, fsb-chain.sh, the
                     verifier re-creation fsb-bos-*.sh and fsb_submission.py, the stand-in sandbox jail2.sh);
                     REPRODUCE-price.sh and REPRODUCE-build.sh: the whole sequence of each step as it ran;
                     fsb_rate_check.py: the rate check of the preflight
Home-folder paths in the logs and scripts are shortened to ~ and <verify>.

Limits.
  * The comparator ran inside the stand-in sandbox jail2.sh (community/base-logs/scripts/), never with landrun:
    the Colab kernel lacks Landlock. Its Lake build found the compiled files up to date, so it checked them with
    both kernels without compiling them again; the verifier re-creation and the independent rebuild compiled the
    modules from source.
  * The leaderboard's own verifier (bigospeedrun.com) has not been run on this package online. The offline
    re-creation used the agents' scripts, not the verifier's own controller, which needs bwrap and systemd-run.
  * Every Lean run was on Colab machines of one kind; none on the Mac except the first run of hostSchedK.
  * The build compiled over a tree that already held the compiled modules of Shea's tree and of community/u8xm.
  * The logs of the agents' second build of community/u8xm, over which this one was compiled, are not in this
    package; community/u8xm holds the logs of its first build.
