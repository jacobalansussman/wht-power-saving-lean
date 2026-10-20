BASE-README.txt  (the build of Shea's base tree, 13 Oct 2026; a Google Colab machine, Ubuntu 24.04, 8 logical CPUs on 4 cores, 50 GB, no swap)
=====================================================================================================================

WHAT IS HERE
  Ryan Shea's lean/complex588 project (github.com/shea256/fourier-transform-below-nlogn, proof commit
  52394849c7ecf3412d39b925d75fc7fdb98994d8; lean/ is byte-identical at repo head 0e61af6), BUILT from source on this
  machine: all 354 Lean modules compiled, Mathlib d13f23b from its cache server, Lean 4.34.1.
  Theorems: OAI.PowerSaving.transform_mainC588 / convolution_mainC588 (DFT of every length, saving 12331685/10^10),
  RAM.hillsC588_program (kernel interface, 12331686/10^10). Credit: Shea (instantiation), Kothe (#588 circuit and
  its community predecessors), Mamano (gauged #401 Lean extension), Dugongue + jacklightChen (#401/#377 geometry),
  Sussman (engine/checkers), OpenAI (model + Fourier reduction).

ENVIRONMENT (every shell and every tmux command)
  . /content/work/env.sh          # ELAN_HOME=/content/work/elan, PATH, MATHLIB_CACHE_DIR=/content/work/mathlib-cache
  elan 4.2.4; lean 4.34.1 (commit 5045d005); /usr/bin/time installed (apt package "time").

PATHS
  /content/work/shea                          the git clone (detached at 5239484). Treat as READ-ONLY.
  /content/work/shea/lean/complex588          the Lake project (lakefile.lean, .lake/build 2.0 GB, .lake/packages 7.7 GB)
  /content/work/ours/Compare.lean             our tools/Compare.lean (sha256 e48d6d3b...dae9; byte-identical to Shea's tools/Compare.lean)
  /content/work/ref/                          OpenAI's original ComparatorChallenges/UniformFourier.lean at fd4aeeb, 3014888, 1581d4c
                                              (all three sha256 4c571f27...ec74) and my own derivations of the C588 challenges
  /content/work/cc/bin/                       comparator d03acab, lean4export 076e8e5, nanoda_bin 3a24072, landrun 811cfff (0.1.18), nounix
  /content/work/cc/jail2.sh                   the sandbox STAND-IN (see COMPARATOR below)
  /content/work/cc/t-c588                     a copy of /content/work/shea used by the comparator (.lake owned by uid 65534)
  /content/work/bos/                          offline re-creation of bigospeedrun/verifier's check of its C588 submission
  /content/work/bos-verifier                  bigospeedrun/verifier clone at c92b9b6 (read only; none of its scripts were run)
  /content/work/logs/                         all logs (cb-*.log, cb-build/, cb-cmp/, cb-bos/); copies on the Mac under
                                              community/base-logs/ (this repository)
  /content/work/rust, /content/work/cc/goroot toolchains used only to build nanoda and landrun

HOW TO OVERLAY (do not edit the base tree)
  cp -a /content/work/shea /content/work/<label>        # about 10 GB; Lake accepts the copy and rebuilds nothing
  cd /content/work/<label>/lean/complex588
  # put new or changed files under Work/... (any module under Work/ is in the lean_lib "Work", glob Work.+;
  # no lakefile change needed). Then rebuild only what changed and its dependents:
  . /content/work/env.sh
  LEAN_NUM_THREADS=2 lake build Work.<Your>.Main <other roots>
  lake build --no-build <roots>        # must print "All targets up-to-date"
  Lake decides by source hashes (traces), so only modified modules and modules importing them recompile.
  The full base build took 1515 s with LEAN_NUM_THREADS=2 (2 Lean processes, each 6-7 GB resident, mostly mapped
  Mathlib files); the largest single process was 9.7 GB. Typical module 5-15 s; slowest C588Sched 42 s, Bridge401.Gates 40 s,
  C588EChk 39 s. RULE: at most 2 heavy Lean processes at once on this machine; run long jobs inside a fresh tmux server:
    tmux -L <name> new -d -s <name> '<cmd> > /content/work/logs/<name>.log 2>&1; echo EXIT $? >> /content/work/logs/<name>.log'

CHECKS (all ran on the base tree; logs in /content/work/logs/cb-build/)
  1. axioms:     lake env lean Work/FourierC588/Axioms.lean
                 -> the three theorems depend on [propext, Classical.choice, Quot.sound]   (8 s, 7.0 GB)
  2. sorry grep: grep -rnw sorry --include=*.lean OAI Work WHTCheck | grep -v '`sorry`'
                 -> only Work/FourierC588/UniformFourierChallenge.lean:327,331 and the same lines of
                    Work/FourierC588Neg/UniformFourierChallenge.lean (the challenge statements). Build log: 4 warnings
                    "declaration uses `sorry`", all in those two challenge files. No axiom declarations, no native_decide,
                    implemented_by, extern, unsafe, ofReduceBool, skipKernelTC. (log: /content/work/logs/cb-source-grep.log)
  3. statement check (our Compare.lean; Shea's is byte-identical):
     lake env lean --run /content/work/ours/Compare.lean Work.FourierC588.UniformFourierChallenge Work.FourierC588.Main \
        OAI.PowerSaving.transform_mainC588 OAI.PowerSaving.convolution_mainC588          -> RESULT: PASS (67-119 s, 16.1 GB)
     negative control (exponent one unit stronger, 12331686):
     ... Work.FourierC588Neg.UniformFourierChallenge ...   -> RESULT: FAIL, exit 1 ("decimalExponentC588: definition body differs")
     The challenge file Work/FourierC588/UniformFourierChallenge.lean is byte-identical to MY OWN derivation from OpenAI's
     original (exponent 1 - 1/10^13 -> 1 - 12331685/10^10, seven names + "C588"): see /content/work/ref/mine-12331685.lean.
     For a new exponent N: derive the challenge the same way (N = numerator over 10^10, SUF = your name suffix):
       perl -Mutf8 -CSD -pe "s{1 - 1/\(10\^\(13:ℕ\)\)}{1 - N/(10^(10:ℕ))}; s/\b(decimalExponent|decimalTime|TimeBounds|DFTGoal|ConvGoal|transform_main|convolution_main)\b/\$1SUF/g" \
         /content/work/ref/UniformFourier-fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb.lean > Work/<Dir>/UniformFourierChallenge.lean )
     and always also run a one-unit-stronger negative control, which must FAIL.
  4. Shea's Python checks: python3 scripts/check_provenance.py -> "PASS: 425 hashes ..."; python3 scripts/regenerate.py ->
     "PASS: 63 new Lean modules regenerated byte for byte" (1.5 s). Python here is 3.13 (his README says 3.11+).

COMPARATOR (official tool, NOT the official sandbox; say so in any report)
  This kernel has no Landlock (landrun refuses in strict mode; with --best-effort, which the comparator uses, it enforces
  nothing) and there is no systemd. Every comparator run therefore goes through /content/work/cc/jail2.sh: new mount/net/IPC
  namespaces, all mounts read-only except <project>/.lake, /dev and a private /tmp (bind of a scratch dir, as the bigospeedrun
  verifier does), uid 65534, no capabilities, no_new_privs, socket(AF_UNIX) blocked by seccomp. Home-made, reviewed by
  nobody. jail.sh of 10 Oct lacked the writable /tmp; with nanoda enabled the comparator writes temp files and failed
  ("read-only file system") under it.
  Run:   /content/work/cb-cmp.sh <label> <project dir whose .lake is owned by 65534> <config.json>
         e.g. /content/work/cb-cmp.sh A1b /content/work/cc/t-c588/lean/complex588 comparator/fourier.json
  (chown -R 65534:65534 <project>/.lake first; log /content/work/logs/cb-cmp/<label>.log, time in <label>.time)
  Shea's config comparator/fourier.json has enable_nanoda: true (nanoda + Lean kernel).
  Results on the base (13 Oct, all through jail2.sh):
    A1b  comparator/fourier.json unchanged, pre-built tree t-c588: "nanoda kernel accepts the solution",
         "Lean default kernel accepts the solution", "Your solution is okay!", exit 0, 389 s, largest process 7.0 GB.
    A2   same with challenge Work.FourierC588Neg.UniformFourierChallenge (one unit stronger): exit 1,
         "Const does not match between challenge and target 'OAI.PowerSaving.decimalExponentC588'", 84 s.
    A3   comparator/fourier.json in a copy with NO module of the repository compiled (t-fresh): the comparator compiled
         the 351 modules of the closure itself as uid 65534 inside the jail, one Lean process at a time (lake adapter
         /content/work/cc/fresh-bin/lake, put first on PATH inside `lake env`, because landrun passes on only
         PATH/HOME/LEAN_ABORT_ON_PANIC and Lake would otherwise start 8 jobs); both kernels accept, exit 0, 3655 s,
         largest process 9.7 GB. Script: /content/work/cc/start-A3.sh.
  For an overlay tree, a from-nothing comparator run needs the same adapter; for a pre-built tree use cb-cmp.sh.

BIGOSPEEDRUN VERIFIER
  Its own controller cannot run here (needs sudo systemd-run with cgroup limits, apparmor, bwrap, a GitHub Actions event).
  /content/work/bos is a re-creation of its prepare()+verify() core: stage = the 351 .lean files of
  submissions/dft/nmamano-shea-c588 (verifier c92b9b6), VerifierChallenge.lean = challenges/dft/UniformFourier233.lean
  (sha256 3a646acd..., = OpenAI fd4aeeb + Y names + header comment) + the Race.bound target at 2466337/2000000000
  (= 12331685/10^10), its lakefile.toml/manifest; phases lean, audit, nanoda through jail2.sh:
     /content/work/cb-bos-run.sh [lean] [audit] [nanoda]     (log /content/work/logs/cb-bos/)
  The comparator compiles all submission modules itself there (thread count from /content/work/bos/control/threads).
  Result 13 Oct: lean phase "Lean default kernel accepts the solution" / "Your solution is okay!", exit 0, 3610 s
  (351 modules compiled inside the jail, 1 thread); audit "'Race.bound' depends on axioms: [propext, Classical.choice,
  Quot.sound]"; nanoda phase "nanoda kernel accepts the solution" + Lean kernel accepts, exit 0, 502 s.
  To check a NEW submission folder this way: rerun cb-bos-prep.sh with SUB and the bound changed (it deletes /content/work/bos).

DISK  /content/work/cc/t-c588, t-fresh and /content/work/bos are ~10 GB each and can be deleted if space is needed;
      never delete /content/work/shea, /content/work/elan, /content/work/mathlib-cache.
