# This repository's scripts, run from a fresh clone (19 October 2026)

To check that `scripts/` works for someone who starts from nothing but this repository, the agents ran it the
way [VERIFY.md](../../VERIFY.md), sections 12.3 and 12.4, says, on a Colab machine (Linux, 8 processors, about
50 GB of memory): a fresh clone of a commit of this repository (sent to the machine as a Git bundle), a fresh copy of
Ryan Shea's built tree at commit `5239484`, then

    bash scripts/apply-overlay.sh package/leverA $T          # lays package/u8x0hq first, then package/leverA
    (cd $T && lake build $(grep -v '^#' package/leverA/MODULES.txt))
    bash scripts/build-and-check.sh package/leverA $T
    bash scripts/build-and-check.sh package/u8x0hq $T
    bash scripts/run-comparator.sh neg $T comparator/fourier-l2979h-neg.json
    bash scripts/run-comparator.sh pos $T comparator/fourier-l2979h.json

`rf-fresh.sh` is the script that ran these steps; `rf-fresh-leverA.out` is its output, and `leverA/` holds the logs
of each step (the per-module build logs and the build output inside the comparator runs are left out). The cloned commit was a test commit made for this run; its
`package/`, `third-party/` and every script of `scripts/` are byte for byte those of this revision (only
`scripts/README.md` and this folder were added to afterwards).

Results for `package/leverA` (the lever on #2979 with every head kept, 1.7564707e-3): the overlay laid with no file
of the tree overwritten; one parallel build of every module (9,507 jobs, 25 minutes, 9.6 GB at most);
`build-and-check.sh` ended `ALL CHECKS PASSED` for `package/leverA` and for `package/u8x0hq`; the official comparator
accepted the proof with both kernels ("Your solution is okay!", 26 minutes, 15.1 GB) and rejected the statement one
unit stronger ("Const does not match ... decimalExponentL2979H").

Preparing this run, the agents found three faults in `build-and-check.sh` with the lever packages, which would have ended in
`SOME CHECK FAILED` although every Lean check passed: it required every Lean control to compile, while the lever
packages also list controls that must fail; it did not accept the two lakefile lines that `package/u8x0hq` adds
when a package requires it; and it looked for the axioms of names that the package's axioms file does not print.
They were fixed before the run above (the script's first lines say what changed), and the fixed script is the one
in this revision.

The same run on `package/leverB` (the lever on #3945 with every head kept, 1.7581520e-3), from the same test commit,
on another Colab machine: the overlay laid the same way; one parallel build of every module (9,629 jobs, 58 minutes,
9.6 GB at most); `build-and-check.sh` ended `ALL CHECKS PASSED` for `package/leverB` and for `package/u8x0hq` (its
steps waited about an hour for the machine, which another job was using); the official comparator accepted the
proof with both kernels ("Your solution is okay!", 88 minutes, 15.0 GB) and rejected the statement one unit
stronger ("Const does not match ... decimalExponentL3945H"). `rf-fresh-leverB.out` is its output and `leverB/` its
logs, kept the same way (the comparator summaries keep the first lines and the verdict lines of each run).

The same run on `package/leverF` (the full cut on #3945 with every head kept, 1.7806935e-3, the headline), on 20
October (UTC), from a second test commit with `package/leverF` added. That commit differs from this revision only
in documentation and in `package/leverF`'s records (its logs, `README.txt`, the checks listed in `result.json`,
`SHA256SUMS`, and the import path of one pricing program); its scripts and the overlay, controls and module list
of `package/leverF`, which are what the scripts read, are byte for byte those of this revision. The script was
`rf-fresh2.sh`, which is `rf-fresh.sh` with the bundle as a setting, a wait for free memory before the comparator,
and the optional `package/u8x0hq` step left out (it had passed twice above). The overlay laid with no file of the
tree overwritten; one parallel build of every module (9,633 jobs, 48 minutes, 9.6 GB at most); `build-and-check.sh`
ended `ALL CHECKS PASSED` for `package/leverF`; the official comparator accepted the proof with both kernels ("Your
solution is okay!", 76 minutes, 16.1 GB) and rejected the statement one unit stronger ("Const does not match ...
decimalExponentL3945F"). `rf-fresh-leverF.out` is its output and `leverF/` its logs.
