community/c588mhh: Kothe's #588 circuit with mutable sources and same-run hosting
====================================================================================

What it proves. The Fourier transform of every length, and convolution, in time O(n (log n)^z) with
z = 1 - 15699305/10^10, a saving of 1.5699305e-3. The statement is OpenAI's UniformFourier.lean with only the
exponent and seven names changed (suffix C588MHH). The theorems are OAI.PowerSaving.transform_mainC588MHH and
OAI.PowerSaving.convolution_mainC588MHH.

What is new. The circuit and its certificate are those of community/c588m: Jeremy Kothe's #588 circuit in Ryan
Shea's Lean tree, with Shea's 112 source copies removed (mutable sources). In community/c588m every leftover
helper piece is stored in a bank: 4,032 banks. Here most of them are stored instead in data registers that sit
idle during the same run of the network, an idea DaysSky proposed for the community's circuits (pull request
#442 of CrocSwap/integer-mult-bounds). Only 111 banks remain. Work/XHost is the Lean proof that this is sound:
  * the role maps of a hosted stage and their injectivity (HRoles), the hosted invocation (HDef, HStage, HCert:
    copies of Nil Mamano's GNetDef, GNetStage and GNetCert in Shea's tree, with hosting added), the frames of a
    host (Frames: host_adapt, host_pad, slot_frames and host_frames, adapted from Mamano's bank_adapt, bank_pad
    and bank_frames in Work/Bank401);
  * host schedules (HSched), and the kernel checker hostSchedK with its soundness lemma
    HostSched.ofTab_fitsR (HSchedK), which turns the published host and bank tables into a proved schedule;
    the checker is modelled on Mamano's schedK and copies its list lemmas;
  * the geometry of the idle climbs of the nine-stage network (HGeom, a copy of the instance part of Mamano's
    GGeom with hosting added) and the end theorem XHost.engine_of_ginv_host (HEnd, adapted from Mamano's
    Work/Bridge401/End.lean), which has the same form as his GX.engine_of_ginv_401, so the Fourier chain of
    community/c588m runs on it unchanged except for its seam (Work/FourierC588MHH);
  * the instance C588MHH (the c588m certificate with the hosting table) and a control, C588MH, the same network
    with an empty host table, which must give exactly the c588m figure 1 - 13552606/10^10 (it does).
One header comment is wrong in part. That of HEnd.lean (written by generators/lean-gen/gen_hgeom.py) says that
Shea's Work/Bridge401/End.lean was written "after Jacob Sussman's Work/Bridge377/End.lean". No such file of mine
exists, in this repository or in Shea's tree: his End.lean says only that it follows a file of that name, and
that the engine, the invocation package and the stage geometry it uses are mine. The comment is left as it was
compiled.
Two lemmas proved before the network was built are in supporting/: L1 restates in the form the plan needed what
was already in Shea's tree (B401.gauged_step and gauged_Wn, GX.entW_spec and entW_form, GX.gauge_step, built on
the entrance gauges of Dugongue's #401, Mamano's files and core_nested of this repository's engine); L2 moves a
helper's leftover piece into a host by a free relabelling of addresses. They were compiled on the Mac and are
not imported by the modules of overlay/: the result does that step with host_adapt (Frames.lean), adapted from
Mamano's bank_adapt. Two independent agents rebuilt them on the Mac (supporting/L1/skeptic, supporting/L2/skeptic).

How to lay it. It needs community/c588m first (result.json: "requires"). From the root of this repository, with
T the folder lean/complex588 of a built copy of Shea's tree at 52394849c7ecf3412d39b925d75fc7fdb98994d8:

    bash community/scripts/apply-overlay.sh community/c588mhh $T     # applies community/c588m, then this package
    bash community/scripts/build-and-check.sh community/c588mhh $T

The package adds 53 files (49 modules in MODULES.txt and 4 comparator configurations); no file of the tree
changes.

Folders.
  overlay/           the files laid over lean/complex588
  data/              the host table (h588_c0_samerun.repacked2.json, from H588-hosttab.json), the price summary of the
                     generator (C588MHH-price.json: W 7279, kernel 15699306, DFT 15699305), the provenance of the
                     Fourier files, and the checksums of the sources uploaded to the Colab (colab-sources.sha256;
                     colab-sources-NOTE.txt lists the ten uploaded files that are not in this package)
  generators/        lean-gen/: the programs that wrote the generated Lean files (gen_hnet.py, gen_hgeom.py,
                     gen_h588m.py, mk_fourier_h.py). Of the twelve network modules of Work/XHost, gen_hnet.py and
                     gen_hgeom.py write six (HRoles, HDef, HStage, HCert, HGeom, HEnd); the other six (HostK,
                     Frames, HSched, HBank, HSchedK, HPriceK) have no generator here, and lean-gen/hsched_text.lean
                     is a copy of HSched.lean. design-tools/: the pricing and the integer programme that chose the
                     host table (host_price.py, host_repack2.py, route_price.py, ...). Two lines of code differ from the copies
                     that ran: the agents' working-folder path in gen_h588m.py is read from the variable BC (the
                     comment on that line points to generators/README.txt, which this package does not have: the
                     working folder is the agents' own), and a local sys.path line of host_price.py is a comment.
                     Besides these, the agents' labels were removed from the docstrings and comments of nine files
                     (gen_hgeom.py, gen_hnet.py, mk_fourier_h.py, gen_h588m.py, emit_hosttab.py, host_price.py,
                     host_repack.py, host_repack2.py, price_hist.py); no other line changed. Regenerating the files
                     from them was not re-run for this package.
  logs/colab/        the agents' build on the first Colab machine, unedited: r<round>-<module>.log and
                     r<round>.summary for rounds r2 to r14 (modules were compiled again after each fix; earlier
                     rounds hold failures), final-*.out and final-*-compare-*.log (Compare.lean), final-*-cmp-*.out
                     and cmp-*.out (the official comparator), cmp/ (its timings)
  logs/rebuild-skeptic/  the independent rebuild: a.summary (32 modules, all exit 0), cmp-pos.log, cmp-neg.log,
                     ax.log, compare.sha
  logs/static-check.txt  tools/community/package_check.py check, run on the Mac
  scripts-as-run/    build/: the Colab scripts as run (hl3-mk.sh compiles with plain lean, one module at a time);
                     rebuild-skeptic/: the rebuild's scripts
  supporting/        L1 and L2 with their Mac logs (home-folder paths shortened to ~ and <verify>)

Limits.
  * The comparator ran inside the stand-in sandbox jail2.sh (community/base-logs/scripts/), never with landrun:
    the Colab kernel lacks Landlock.
  * The leaderboard's own verifier (bigospeedrun.com) has not been run on this package, offline or online.
  * Every Lean run was on Colab machines of one kind; none on the Mac except the supporting lemmas and the
    negative control below.
  * The build compiled over a tree that already held the compiled modules of Shea's tree and of community/c588m.
  * A Lean negative control of hostSchedK (one remainder off by one, rejected) ran on the Mac on an earlier split
    of the data modules; it is not in this package.
