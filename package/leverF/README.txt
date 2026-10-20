package/leverF: the recursion lever (full cut) on the community's #3945 circuit with every head kept, saving 1.7806935e-3

What it is
  The new Lean files of this result: OpenAI's statement (the Fourier transform of every length, and convolution)
  with the exponent 1 - 17806935/10^10, that is a power saving of 1.7806935e-3. Laid over Ryan Shea's tree after
  package/u8x0hq (result.json, "requires"); every file is under Work/ or comparator/, so the lakefile does not
  change. Final theorems: OAI.PowerSaving.transform_mainL3945F and OAI.PowerSaving.convolution_mainL3945F, in
  Work.FourierL3945F.Main. Their axioms are propext, Classical.choice and Quot.sound.

The circuit, the fill and the lever
  #3945 is Dugongue's circuit (staged544 T360 terminal-gap fill) on the community's staged544 line; its own figure
  is 1.7306632e-3. It runs 360 copies of a 9-call pattern: 3,240 calls of two provider words, 322,560 arrays per
  gauge copy, 1,321,920 hosted roles. The paid profile of one gauge copy has 137,952 blocks of rank 48.
  This package uses the same fill as package/leverB: I re-allocated the hosted roles so that every block in
  children 8 and 10 is empty (fill/), which leaves the words, the call order, the dispatch, the template, the plan
  and the paid profile of #3945 unchanged. In my lowering, 80,640 of the rank-48 blocks (40,320 + 40,320, all of
  children 8 and 10) are then literal heads: the first operation on their own array, acting on that array alone.
  The other 57,312 (children 12 and 13) are whole idle blocks of rank 48 that host no role.
  The full cut does two things:
  - it removes each head entirely, so its array leaves the recursion with target 73 - 48 = 25 instead of 73;
  - it splits each of the other 57,312 rank-48 blocks into two blocks of rank 24 on the same role. A block of rank
    24 + 24 is the product of its first 24 and its last 24 directions, so the split changes no output.
  Every paid rank is then at most 24, below the smaller target 25, so the mixed-target recursion applies with the
  two targets 73 (241,920 arrays per gauge copy) and 25 (80,640). The rate equation of the full-cut profile has its
  root at 1.78069354058e-3; the Lean figure 1.7806935e-3 is the grid point just below it, and the rate fact is
  false one grid step higher (controls/Rate3945FSharpCtl.lean). The half cut on the same fill (package/leverB) gives
  1.7581520e-3; the full cut is 1.28 % above it and 2.89 % above #3945's own figure. The price was found by four
  independent routes (fill/checks/full-cut/).

What the files are (paths relative to lean/complex588)
  Work/Mixed/ (Flow, Phi, LemmaA, Bound, TwoTarget)   the mixed-target recursion: weights of the class tree, exact
                                                     balance, sources and exits, and the cost bound for every exponent
                                                     above the rate point
  Work/Mixed/Prog/, Work/MixedCost/                  the lockstep tree of calls as a program in OpenAI's machine
                                                     model, its correctness and its cost
  Work/Mixed/ (SeamConst, Rate3945F)                 the seam with a constant; the rate fact of the full cut at the
                                                     figure
  Work/Mixed/HeadM1                                  the head lemma, the full cut of a head (cut_heads_full) and the
                                                     split of a block (splitWord, gwalk_splitMove; the file's
                                                     split_block is not used by the proof)
  Work/Mixed/Bridge/CoreF, Inst3945F, HFTest3945FW   the full-cut bridge: from a heads-first network to the
                                                     recursion's hypotheses with the heads removed and the other
                                                     rank-48 blocks split 24 + 24; its instance on this fill; the
                                                     same on the warp network of the #3945 provider words
  Work/Mixed/Bridge/ (Core, HeadsFirst, Inst,        the generic bridge and the heads-first network of this fill,
  HFTest, Inst3945H, HeadsFirst3945H, HFTest3945HW)  shared with package/leverB
  Work/KLowW/, Work/Warp/, Work/PWGen/               the warp checker, what a warp means, the word helpers, and the
                                                     pair of #3945 provider words (S5Pair)
  Work/GCert/Data/ (S5B*, S5E*, Gen/S5*)             both provider words as checked Lean data (generated)
  Work/XHostD/, Work/XHostR/,                        the hosted network with 360 copies and hosting across copies,
  Work/XHostT/ (T*, K*, HF/, HFW2979/)               its warp variant, the host-table checker, the heads-first
                                                     certificate on the warp network (Work/XHost/ itself comes with
                                                     package/u8x0hq)
  Work/XHostT/Data/S5H/, Data/D3945/, Inst/S5H*      the host table of this fill (all 1,321,920 roles), the dispatch,
                                                     the schedule and its facts: 80,640 heads, the unit price
                                                     (generated data and their glue, as in package/leverB)
  Work/L3945H/Prov                                   the network's data (providers, dispatch, windows, residuals,
                                                     unit price), shared with package/leverB
  Work/L3945F/ (Stmt, EngineOfGeomW, Engine)         the engine statement and the one engine theorem the chain
                                                     imports
  Work/FourierL3945F/, Work/FourierL3945FNeg/        the Fourier chain at this exponent; the two challenge files are
                                                     OpenAI's statement with only the exponent and names changed
  Work/BlockAccounting/Sharp.lean                    this project's rate lemmas, needed by the rate files
  Work/Mixed/Rate3945.lean, Rate3815.lean,           rate modules of earlier builds that the generic bridge files
  Rate3945H.lean                                     import; not used for this figure
  comparator/fourier-l3945f.json, -neg.json          the official comparator's configurations
  MODULES.txt lists all 410 modules in build order; ORIGINS.tsv says which part of the work wrote each file.

Checks
  [pass]    build (19 Oct 2026, Colab C3; first built by fc-engine): a copy of the tree in which fc-engine built every
            file (logs/fc-engine/); the 567 + 20 overlay files installed identical (0 new); lake build 9633 jobs, all
            up to date (logs/fc-final/C3/final.out)
  [pass]    build from the release overlay (Colab C5): Ryan Shea's tree with package/u8x0hq (prebuilt) plus this
            package's overlay/ (412 files, hash-checked), on a second machine: lake build 9633 jobs completed, 1 h 08
            min; the final theorems' axioms are propext, Classical.choice, Quot.sound
            (logs/fc-final/C5/identC5/axioms-fresh.log) (logs/fc-final/C5/identC5.out)
  [pass]    axioms: 9 theorems (the two final theorems and the chain's seven named steps; logs also
            l3945f/axioms-extra.log): propext, Classical.choice, Quot.sound (logs/fc-final/C3/l3945f/axioms.log)
  [pass]    source greps: 583 files of the closure that are not Ryan Shea's: 0 lines with sorry or admit (OpenAI's two
            challenge statements aside), 0 with axiom, native_decide, unsafe, implemented_by or extern
            (sorry-grep.txt) (logs/fc-final/C3/l3945f/forbidden-grep.txt)
  [pass]    statement files derived from OpenAI's: both challenge files byte-identical to the derivation (17806935 and
            17806936) (logs/fc-final/C3/l3945f/derived-FourierL3945F.lean)
  [pass]    Compare.lean positive: RESULT: PASS, 56 s, 16.1 GB (logs/fc-final/C3/l3945f/compare-pos.log)
  [pass]    Compare.lean negative: RESULT: FAIL on decimalExponentL3945F (one grid step higher, 1 - 17806936/10^10)
            (logs/fc-final/C3/l3945f/compare-neg.log)
  [pass]    official comparator positive (stand-in sandbox): nanoda and Lean's kernel accept; 1 h 06 min, 14.9 GB
            (l3945f/cmp/l3945f-pos.log) (logs/fc-final/C3/L3945F-cmp.out)
  [pass]    official comparator negative: rejected: Const does not match (decimalExponentL3945F); 4 min
            (logs/fc-final/C3/L3945F-cmp.out)
  [pass]    leaderboard verifier re-created offline (c92b9b6), merged stage (C5): the import closure has 707 local
            modules, over the pull-request limit of 500; merged to 448 modules (450 files with Submission.lean and
            manifest.json), static check passed; Lean check from source 2 h 07 min (10.0 GB), axiom audit standard (3
            min), nanoda 1 h 28 min; the bound one grid step higher is rejected (logs/fc-final/C5/bosM.out)
  [pass]    merged stage = built tree, constant by constant: h-bos's fast identity check on C5: the 29,540
            non-auxiliary constants of the closure are identical in the merged stage and in the tree built from the
            release overlay, and so are the statements and proof terms of all 12,330 theorems (the adapter module
            aside) (logs/fc-final/C5/identC5.out)
  [pass]    leaderboard verifier re-created offline, unmerged stage (C2, sources only): the 707-module form for the
            verifier's issue route (no file limit), compiled from sources only: Lean check 2 h 03 min, audit standard,
            nanoda 1 h 25 min; the bound one grid step higher is rejected (logs/fc-final/C2/bosU.out)
  [pass]    release controls: Heads3945HNeg, Rate3945FSharpCtl and Ctl3945F compile with their axioms; NegHeads,
            NegHeadsOld, NegRes, NegGlob fail (Lean's decide refutes each) (logs/fc-final/C3/ctl.out)
  [pass]    side modules (controls/side/): fc-rate's R5Joint3945F (the rate fact against the recursion's rule,
            separation 24 < 25), Hist3945F (unit to cut profile) and FullCut3945F (the full-cut word from the head
            lemma) built in the same tree; standard axioms (the build log prints each theorem's axioms)
            (logs/fc-final/C3/side/build.log)
  [pass]    generators re-run (generators/REGENERATE.sh): 192 generated overlay files made again byte-identical; each
            generator refuses one step too far (logs/regenerate/regenerate.out)
  [pass]    price of the full cut, independent routes: the full cut prices at 1.7806935e-3 by four routes (the lever
            pricer; exact rational interval arithmetic; the Lean data generator's scan; the community moment engine's
            enclosures); one grid step higher is refused by three of them; the profile was recounted from the fill
            (fill/checks/full-cut/, logs/fc-gate/GATE.md) (logs/fc-gate/GATE.md)
  [pass]    the community's complete admission replay of #3945 on this fill:
            PASS_TERMINAL_WIDTH5_COMPLETE_CONDITIONAL_RELEASE with fresh native runs; 119 of 119 strict controls
            rejected (logs/ah-gate/GATE.md; the same fill as package/leverB) (logs/ah-gate/frB/fullB.out)
  result.json lists each check with its log under logs/.

Controls (controls/, compiled against the built tree; logs/fc-final/C3/controls/)
  Heads3945HNeg.lean     Lean's kernel decides that the host table does NOT have 80,641 heads, and that NOT every
                         residual is at least 6 (the claims one step too far); compiles, axioms propext.
  Rate3945FSharpCtl.lean the full cut's rate fact is false at 1 - 17806936/10^10 and at every larger saving;
                         compiles, standard axioms.
  Ctl3945F.lean          the full-cut profile identity is false with #3945's published head count (57,312), false
                         with one head more (80,641), and false without the split (a rank-48 block would remain,
                         above the target 25); compiles, axioms propext.
  NegHeads.lean, NegHeadsOld.lean, NegRes.lean, NegGlob.lean   must FAIL to compile (80,641 heads; #3945's own
                         count of 57,312 on this table; a width bound of 6; one table chunk count one too high):
                         they do.
  side/                  side modules built in the same tree, never imported by the final theorems:
                         Rate3945FSharp (the rate control as a module) and RateMono; R5Joint3945F (the rate fact
                         against the recursion's rule, the separation 24 < 25); Hist3945F (the histogram from the
                         unit profile to the cut profile, its rank mass and block count); FullCut3945F (the full-cut
                         word composed from the head lemma, as a cross-check of the bridge); S5HHB (the head bins of
                         the schedule are 80,640).
  side/perturbed/        the host table with a load moved and with a window moved as Lean parts: the kernel refutes
                         each moved piece (package/leverB's, the same table).

The fill (fill/)
  The same files as in package/leverB: FILL.json.gz (sha256 9c06677d...) and PROFILE.json are the fill in #3945's
  format; old-FILL.json.gz and old-PROFILE.json are the same allocation at the #3815 stage, which #3945's offset
  step reads; made-by/ holds the programs and the record of the run that made them; checks/ holds the head counts
  and prices of the fill measured by two independent evaluators. checks/full-cut/ adds the prices of the full cut
  by four routes (the lever pricer, my own exact rational interval code, the Lean data generator's scan, and the
  community moment engine's own enclosures), with the programs (tools/), and a recount of the profile from the fill.

Generators
  generators/REGENERATE.sh makes every generated file of this package again from the community's #3945 owner folder
  and fill/, and compares with SHA256SUMS; generators/README.txt says what makes what.

Known cosmetic issues
  KNOWN-ISSUES.txt lists them (names kept from package/leverB, unused rate modules in the closure, and the like).
  None affects a theorem, a check or a figure.

Credits
  The circuit: #3945 by Dugongue (staged544 T360 terminal-gap fill), on the lineage T74 words, #2979, #3133/#3148,
  #3605 (Ross Richardson), #3642, #3682 (zxv), #3734 (dataisfire: the staged544 words and owner interface),
  #3802/#3875 (hexapode) and #3815 (Erwan Secousses / ChaloupeMusic: the clock and the fill construction whose
  programs made this fill); the physical owner predicates of #3638, the warps and dirty borrowing of #2226 (Cinna
  Davis), Sagar Thakrar (#3687), roger-sato (#3693), Jeremy Kothe, Tekalign, icekylinx, DreamingOfClouds and
  chafreaky. Under it: Ryan Shea's Lean tree (generators, bank and seam pipeline), Nil Mamano's gauged Lean extension
  and record files, Dugongue's #401 and jacklightChen's #377 banks, and OpenAI's model, statement and Fourier
  reduction. Mine: the mixed-target recursion, the every-head-kept fill, the lever with its full cut, and the Lean
  that connects them to the circuit.
