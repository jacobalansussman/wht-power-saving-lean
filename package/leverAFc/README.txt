package/leverAFc: the recursion lever (full cut) on the community's #2979 words with a re-allocated fill, saving 1.7887045e-3

What it is
  The new Lean files of this result: OpenAI's statement (the Fourier transform of every length, and convolution)
  with the exponent 1 - 17887045/10^10, that is a power saving of 1.7887045e-3. Laid over Ryan Shea's tree after
  package/u8x0hq (result.json, "requires"); every file is under Work/ or comparator/, so the lakefile does not
  change. Final theorems: OAI.PowerSaving.transform_mainL2979Fc and OAI.PowerSaving.convolution_mainL2979Fc, in
  Work.FourierL2979Fc.Main. Their axioms are propext, Classical.choice and Quot.sound.

The circuit, the fill and the lever
  #2979 is Marius-Paul Dumitrean's circuit (rank-7 joint end frames, on Jeremy Kothe's #2957 and Mason Hale's #2917
  words), merged on the community repository as research; its own figure is 1.7290340e-3. It runs 74 copies of a
  9-call pattern: 666 calls of two provider words, 66,304 arrays per gauge copy, 271,728 hosted roles. The fill says
  which call hosts each role.
  This package keeps #2979's words, clock, call order, dispatch and providers, and uses a new fill (fill/). I made it
  for the full cut rather than for #2979's own figure: with #3815's allocation model (Erwan Secousses), I solved for
  the fill whose profile is cheapest AFTER the cut, and placed the roles with #3815's own placement program, made
  tolerant. The new fill leaves every bin of children 8 and 10 empty, and no other block idle at rank 48.
  - In my lowering, its 16,576 rank-48 blocks (8,288 in child 8 and 8,288 in child 10) are all literal heads: each
    is the first operation on its own array and acts on that array alone.
  - The other idle blocks have ranks 1, 6, 7, 17, 19 and 23. The paid profile of one gauge copy is
    1: 1,630,639; 2: 470,196; 3: 154,512; 4: 71,928; 5: 11,988; 6: 5; 7: 3; 8: 5,328; 17: 8; 19: 1; 23: 25,905;
    48: 16,576.
  - Without the lever this fill is worse than #2979's own: its uniform figure is 1.6897030e-3. It pays off only with
    the full cut.
  - The placement needed 9 spare roles: 8 roles of width 6 at call 553 and 1 of width 4 at call 0 had no planned slot
    left, and went to free places that obey every rule of the community's checks. The community's own owner engine
    for #2979 (from #2839) and its other fill checks accept the fill as a fill (fill/checks/hub/).
  The full cut removes each head entirely, so its array leaves the recursion with target 73 - 48 = 25 instead of 73.
  No other block has rank 48, so nothing has to be split. Every paid rank is then at most 23, below the smaller
  target 25, so the mixed-target recursion applies with the two targets 73 (49,728 arrays per gauge copy) and 25
  (16,576). The rate equation of the cut profile has its root at 1.78870456430e-3; the Lean figure 1.7887045e-3 is
  the grid point just below it, and the rate fact is false one grid step higher (controls/Rate2979FcSharpCtl.lean).
  The figure is 0.45 % above the full cut on #3945's words with every head kept (package/leverF, 1.7806935e-3) and
  3.45 % above #2979's own figure.

What the files are (paths relative to lean/complex588)
  Work/Mixed/ (Flow, Phi, LemmaA, Bound, TwoTarget)   the mixed-target recursion: weights of the class tree, exact
                                                     balance, sources and exits, and the cost bound for every exponent
                                                     above the rate point
  Work/Mixed/Prog/, Work/MixedCost/                  the lockstep tree of calls as a program in OpenAI's machine
                                                     model, its correctness and its cost
  Work/Mixed/ (SeamConst, Rate2979Fc, Hist2979Fc)    the seam with a constant; the rate fact of the full cut at the
                                                     figure; the histogram from the unit profile to the cut profile
  Work/Mixed/HeadM1                                  the head lemma, the full cut of a head (cut_heads_full) and the
                                                     split of a block (splitWord, here splitting nothing)
  Work/Mixed/Bridge/ (Core, CoreF, Inst2979Fc,       the generic bridge and the full-cut bridge: from a heads-first
  HFTest2979FcW)                                     network to the recursion's hypotheses with the heads removed;
                                                     its instance on this fill; the same on the warp network of the
                                                     #2979 provider words, with the fact that a head's 48 directions
                                                     are the 48 lines that the cut removes
  Work/KLowW/, Work/Warp/, Work/PWGen/               the warp checker, what a warp means, the word helpers, and the
                                                     pair of #2979 provider words (S2979Pair)
  Work/GCert/Data/ (S2979B*, S2979E*, Gen/S2979*)    both provider words as checked Lean data (generated)
  Work/XHostD/, Work/XHostR/,                        the hosted network with 74 copies and hosting across copies,
  Work/XHostT/ (T*, K*, HF/, HFW2979/)               its warp variant, the host-table checker, the heads-first
                                                     certificate on the warp network (Work/XHost/ itself comes with
                                                     package/u8x0hq)
  Work/XHostT/Data/H2979Fc/, Data/D2979/,            the host table of this fill (all 271,728 roles), the dispatch,
  Inst/S2979Fc, S2979FcFacts                         the schedule and its facts: 16,576 heads, no other rank-48
                                                     block, the idle price (generated data and their glue)
  Work/L2979Fc/ (Stmt, EngineOfGeomW, Prov, Engine)  the engine statement, the network's data, and the one engine
                                                     theorem the chain imports
  Work/FourierL2979Fc/, Work/FourierL2979FcNeg/      the Fourier chain at this exponent; the two challenge files are
                                                     OpenAI's statement with only the exponent and names changed
  Work/BlockAccounting/Sharp.lean                    this project's rate lemmas, needed by the rate files
  comparator/fourier-l2979fc.json, -neg.json         the official comparator's configurations
  MODULES.txt lists all 280 modules in build order; ORIGINS.tsv says which part of the work wrote each file.

Checks
  [pass]    build (20 Oct 2026, Colab C4; first built by afc-data and afc-cut): a copy of the tree in which afc-data
            and afc-cut built every file (logs/afc-data/, logs/afc-cut/); the 282 overlay files and the 20 stage files
            installed identical (0 new); lake build 9503 jobs, all up to date (logs/afc-final/C4/final.out)
  [pass]    build from the release overlay (Colab C5): Ryan Shea's tree with package/u8x0hq (prebuilt) plus this
            package's overlay/ (282 files, hash-checked), on a second machine: lake build 9503 jobs completed, 38 min;
            the final theorems' axioms are propext, Classical.choice, Quot.sound
            (logs/afc-final/C5/identC5/axioms-fresh.log) (logs/afc-final/C5/identC5.out)
  [pass]    axioms: 9 theorems (the two final theorems and the chain's seven named steps; logs also
            l2979fc/axioms-extra.log): propext, Classical.choice, Quot.sound (logs/afc-final/C4/l2979fc/axioms.log)
  [pass]    source greps: the 278 Lean files of this package's overlay other than the two challenge files (which hold
            OpenAI's own two sorry lines each): 0 lines with sorry or admit, 0 with axiom, native_decide, unsafe,
            implemented_by or extern (sorry-grep.txt) (logs/afc-final/C4/l2979fc/forbidden-grep.txt)
  [pass]    source greps over the whole closure: every one of the 577 Lean files of the import closure as packaged for
            the verifier (Ryan Shea's, package/u8x0hq's and this package's; the challenge files are not part of it): 0
            lines with sorry or admit, 0 with axiom, native_decide, unsafe, implemented_by or extern
            (C2/closure-grep/) (logs/afc-final/C2/closure-grep/sorry-grep.txt)
  [pass]    statement files derived from OpenAI's: both challenge files byte-identical to the derivation (17887045 and
            17887046) (logs/afc-final/C4/l2979fc/derived-FourierL2979Fc.lean)
  [pass]    Compare.lean positive: RESULT: PASS, 77 s, 15.9 GB (logs/afc-final/C4/l2979fc/compare-pos.log)
  [pass]    Compare.lean negative: RESULT: FAIL on decimalExponentL2979Fc (one grid step higher, 1 - 17887046/10^10)
            (logs/afc-final/C4/l2979fc/compare-neg.log)
  [pass]    official comparator positive (stand-in sandbox): nanoda and Lean's kernel accept; 28 min, 14.6 GB
            (l2979fc/cmp/l2979fc-pos.log) (logs/afc-final/C4/L2979Fc-cmp.out)
  [pass]    official comparator negative: rejected: Const does not match (decimalExponentL2979Fc); 2 min
            (logs/afc-final/C4/L2979Fc-cmp.out)
  [pass]    leaderboard verifier re-created offline (c92b9b6), merged stage (C5): the import closure has 577 local
            modules, over the pull-request limit of 500; merged to 447 modules (449 files with Submission.lean and
            manifest.json), static check passed; Lean check from source 1 h 03 min (10.0 GB), axiom audit standard (2
            min), nanoda 37 min; the bound one grid step higher is rejected (logs/afc-final/C5/bosM.out)
  [pass]    merged stage = built tree, constant by constant: h-bos's fast identity check on C5: the 21,665
            non-auxiliary constants of the closure are identical in the merged stage and in the tree built from the
            release overlay, and so are the statements and proof terms of all 9,601 theorems (the adapter module
            aside) (logs/afc-final/C5/identC5.out)
  [pass]    leaderboard verifier re-created offline, unmerged stage (C2, sources only): the 577-module form for the
            verifier's issue route (no file limit), compiled from sources only: Lean check 1 h 12 min, audit standard,
            nanoda 37 min; the bound one grid step higher is rejected (logs/afc-final/C2/bosU.out)
  [pass]    release controls: Heads2979FcNeg, Rate2979FcSharpCtl and Ctl2979Fc compile with their axioms; NegHeads,
            NegGlob, NegF1 and the perturbed pieces NegSome119, NegSome46, NegSome334 fail (Lean's decide refutes
            each); the unperturbed NegNone* compile (logs/afc-final/C4/ctl.out)
  [pass]    side modules (controls/side/): R5Joint2979Fc (the rate fact against the recursion's rule, separation 23 <
            25), FullCut2979Fc (the full-cut word from the head lemma; the split keeps the histogram),
            Rate2979FcSharp, and the schedule's side checks S2979FcHB, S2979FcHFW, S2979FcUnit, built in the same
            tree; standard axioms or fewer (the build log prints each theorem's axioms)
            (logs/afc-final/C4/side/build.log)
  [pass]    generators re-run (generators/REGENERATE.sh): 179 generated overlay files made again byte-identical; each
            generator refuses one step too far (logs/regenerate/regenerate.out)
  [pass]    price of the full cut, independent routes: the full cut prices at 1.78870456430e-3 by three independent
            routes (exact rational enclosures on the 1e-24 grid, a 60-digit root, 70-digit interval arithmetic) and by
            the Lean rate generator's scan, whose largest grid numerator is 17887045 (17887046 refused)
            (fill/checks/price/hx-price-NOTE.md)
  [pass]    the community's checks of #2979's line on this fill: the community's own checks of #2979's line,
            unmodified, on this fill: #2979's PR2839 owner engine with its controls, verify.py, the complete replay's
            fresh natives and the admission's owner, lift, rebind, discovery and outer-ledger steps, and the release's
            verify.py pass; it is rejected only by the steps that ask for a higher uniform figure than the one it
            improves on (its uniform figure is 1.6897030e-3) (fill/checks/hub/hx-admit-NOTE.md)
  [pass]    gate: the cut fact on these words: the 48-line cut fact on #2979's words (a head's 48 directions are the
            48 lines the cut removes) compiled on its own first, on the design-A tree, with the heads-first
            certificate and the full-cut bundle at this fill's numbers (logs/afc-gate/GATE.md)
  result.json lists each check with its log under logs/.

Controls (controls/, compiled against the built tree; logs/afc-final/C4/controls/)
  scripts/build-and-check.sh of the repository runs the five that need only this package's overlay
  (result.json, lean_controls: Heads2979FcNeg, Rate2979FcSharpCtl, NegHeads, NegGlob, NegF1). The other seven need
  more in the tree (result.json, lean_controls_extra, says what): Ctl2979Fc needs package/leverA's files, and the
  perturbed/ controls need the two perturbed tables of side/ built first. All twelve were compiled on the build
  machine, whose tree had both.
  Heads2979FcNeg.lean    Lean's kernel decides that the host table does NOT have 16,577 heads, and that NOT every
                         residual is at least 5 (the claims one step too far); compiles, axioms propext.
  Rate2979FcSharpCtl.lean the full cut's rate fact is false at 1 - 17887046/10^10 and at every larger saving;
                         compiles, standard axioms.
  Ctl2979Fc.lean         the full-cut profile identity is false with one head fewer, and false for the unit profile
                         of the every-head-kept fill on #2979 (package/leverA): that profile with its heads removed
                         still has 11,614 rank-48 blocks, above the target 25, while this one has none; compiles,
                         axioms propext. It imports package/leverA's Work/Mixed/Bridge/Inst2979H.lean for that
                         comparison, so it compiles in a tree that also holds package/leverA's files.
  NegHeads.lean, NegGlob.lean, NegF1.lean   must FAIL to compile (16,577 heads; one table chunk count one too high;
                         one idle rank-48 block outside children 8 and 10): they do.
  perturbed/             the host table with a load moved (NegSome119) and with a window moved (NegSome46,
                         NegSome334) must FAIL to compile; the same pieces unmoved (NegNone*) compile.
  side/                  side modules built in the same tree, never imported by the final theorems (paths as in the
                         tree): Rate2979FcSharp (the rate control as a module) and RateMono; R5Joint2979Fc (the rate
                         fact against the recursion's rule, the separation 23 < 25); FullCut2979Fc (the full-cut word
                         composed from the head lemma, as a cross-check of the bridge; the split keeps the
                         histogram); S2979FcHB (the head bins of the schedule are 16,576), S2979FcHFW (the heads-first
                         certificate at this schedule), S2979FcUnit (the unit price on the plain network) and the
                         generic modules they import; the two perturbed tables as Lean modules.

The fill (fill/)
  FILL.json.gz (sha256 7fc766c6...) and PROFILE.json are the fill in #2979's format; ORDER.json and DISPATCH.json
  are the same order and dispatch as #2979's (equal as JSON; the Lean data were made from #2979's own files).
  made-by/ holds the programs, the class solution and the record of the run that made the fill; checks/ holds its
  head counts and price by independent evaluators (eval/, price/) and the run of the community's own checks on it
  (hub/). fill/README.txt says more.

Generators
  generators/REGENERATE.sh makes every generated file of this package again from the community's #2979 files and
  fill/, and compares with SHA256SUMS; generators/README.txt says what makes what. It needs Python 3.11 or earlier
  (KNOWN-ISSUES.txt, item 11).

Known cosmetic issues
  KNOWN-ISSUES.txt lists them (comments that name earlier builds, one profile in several files, a comment above the
  Fourier theorem that names a Walsh-Hadamard theorem that does not exist, line numbers of OpenAI's Goal.lean cited one
  off, and the like). None affects a theorem, a check or a figure. One affects a command: generators/REGENERATE.sh
  needs Python 3.11 or earlier (item 11 there, with the workaround).

Credits
  The circuit: #2979 by Marius-Paul Dumitrean (rank-7 joint end frames, merged as research), on Jeremy Kothe's
  #2957/#2948/#2966 rewrites and Mason Hale's #2917 words, T74 chronology and allocation checker; with Jeremy
  Diviney (#2908), Dugongue (#2969, and #2863/#2839: the owner engine that checks this fill), RossRichardson (#2918),
  Yiyang Lu (#2857/#2848), Tekalign, icekylinx, themilkman216 and Gareth Davidson (#1995/#2013); the warps and dirty
  borrowing of #2226 (Cinna Davis) and the helper width exchange of #2312 (Jeremy Kothe). The fill was made with the
  allocation model and placement program of #3815 (Erwan Secousses / ChaloupeMusic). Under it: Ryan Shea's Lean tree
  (generators, bank and seam pipeline), Nil Mamano's gauged Lean extension and record files, Dugongue's #401 and
  jacklightChen's #377 banks, and OpenAI's model, statement and Fourier reduction. Mine: the mixed-target recursion,
  the re-allocated fill, the lever with its full cut, and the Lean that connects them to the circuit.
