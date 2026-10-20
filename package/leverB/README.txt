package/leverB: the recursion lever on the community's #3945 circuit with every head kept, saving 1.7581520e-3

What it is
  The new Lean files of this result: OpenAI's statement (the Fourier transform of every length, and convolution)
  with the exponent 1 - 17581520/10^10, that is a power saving of 1.7581520e-3. Laid over Ryan Shea's tree after
  package/u8x0hq (result.json, "requires"); every file is under Work/ or comparator/, so the lakefile does not
  change. Final theorems: OAI.PowerSaving.transform_mainL3945H and OAI.PowerSaving.convolution_mainL3945H, in
  Work.FourierL3945H.Main. Their axioms are propext, Classical.choice and Quot.sound.

The circuit, the fill and the lever
  #3945 is Dugongue's circuit (staged544 T360 terminal-gap fill) on the community's staged544 line; its own figure
  is 1.7306632e-3. It runs 360 copies of a 9-call pattern: 3,240 calls of two provider words, 322,560 arrays per
  gauge copy, 1,321,920 hosted roles. The community's construction places each copy's empty rank-48 blocks either in
  children 8 and 10 or in children 12 and 13, at the same price; as published, 57,312 of them sit in children 8 and 10.
  I re-allocated the hosted roles so that every child-8/10 block is empty (fill/: the community's own LP and placement
  programs of #3815 with a head bonus in the LP objective, an exact MIP on the generated columns, then #3945's offset
  step). The words, the call order, the dispatch, the template, the plan and the paid profile are #3945's, unchanged;
  #3945's own owner check passes on the new fill, and so does the community's complete admission replay. In my lowering
  all 80,640 child-8/10 rank-48 blocks (40,320 + 40,320, the most there can be: 224 per copy) are now literal heads:
  the first operation on their own array, acting on that array alone. I cut each head by 24 of its 48 axes, so it
  becomes a block of rank 24 that leaves its array with target 49; every paid rank stays at most 48, below 49, so the
  mixed-target recursion applies with the two targets 73 and 49. The rate equation of the cut profile has its root at
  1.75815201676e-3; the Lean figure 1.7581520e-3 is the grid point just below it, and the rate fact is false one grid
  step higher (controls/Rate3945HSharpCtl.lean).

What the files are (paths relative to lean/complex588)
  Work/Mixed/ (Flow, Phi, LemmaA, Bound, TwoTarget)   the mixed-target recursion: weights of the class tree, exact
                                                     balance, sources and exits, and the cost bound for every exponent
                                                     above the rate point
  Work/Mixed/Prog/, Work/MixedCost/                  the lockstep tree of calls as a program in OpenAI's machine
                                                     model, its correctness and its cost
  Work/Mixed/ (SeamConst, Rate3945H)                 the seam with a constant; the rate fact at the figure
  Work/Mixed/HeadM1                                  the head lemma and the half cut (24 + 24)
  Work/Mixed/Bridge/                                 from a heads-first network to the recursion's hypotheses
                                                     (Core, HeadsFirst, Inst, HFTest; the instance of this fill
                                                     Inst3945H, HeadsFirst3945H, HFTest3945HW)
  Work/KLowW/, Work/Warp/, Work/PWGen/               the warp checker, what a warp means, the word helpers, and the
                                                     pair of #3945 provider words (S5Pair)
  Work/GCert/Data/ (S5B*, S5E*, Gen/S5*)             both provider words as checked Lean data (generated)
  Work/XHostD/, Work/XHostR/,                        the hosted network with 360 copies and hosting across copies,
  Work/XHostT/ (T*, K*, HF/, HFW2979/)               its warp variant, the host-table checker, the heads-first
                                                     certificate on the warp network (Work/XHost/ itself comes with
                                                     package/u8x0hq)
  Work/XHostT/Data/S5H/, Data/D3945/, Inst/S5H*      the host table of this fill (all 1,321,920 roles), the dispatch,
                                                     the schedule and its facts: 80,640 heads, the unit price
                                                     (generated data and their glue)
  Work/L3945H/ (Stmt, EngineOfGeomW, Prov, Engine)   the engine statement, the network's data, and the one engine
                                                     theorem the chain imports
  Work/FourierL3945H/, Work/FourierL3945HNeg/        the Fourier chain at this exponent; the two challenge files are
                                                     OpenAI's statement with only the exponent and names changed
  Work/BlockAccounting/Sharp.lean                    this project's rate lemmas, needed by the rate files
  Work/Mixed/Rate3945.lean, Rate3815.lean            rate modules of the earlier #3945 build that the generic bridge
                                                     files import; not used for this figure
  comparator/fourier-l3945h.json, -neg.json          the official comparator's configurations
  MODULES.txt lists all 406 modules in build order; ORIGINS.tsv says which part of the work wrote each file.

Checks
  [pass]    build (19 Oct 2026, Colab C3; first built by bB-data, logs/bB-data/): a copy of bB-data's built tree; the
            561 + 20 overlay files installed identical (0 new); lake build 9629 jobs, all up to date
            (logs/bB-final/C3/final.out)
  [pass]    axioms: 9 theorems (the two final theorems and the chain's seven named steps; logs also
            l3945h/axioms-extra.log): propext, Classical.choice, Quot.sound (logs/bB-final/C3/l3945h/axioms.log)
  [pass]    source greps: 577 files of the closure that are not Ryan Shea's: 0 lines with sorry or admit (OpenAI's two
            challenge statements aside), 0 with axiom, native_decide, unsafe, implemented_by or extern
            (sorry-grep.txt) (logs/bB-final/C3/l3945h/forbidden-grep.txt)
  [pass]    statement files derived from OpenAI's: both challenge files byte-identical to the derivation (17581520 and
            17581521) (logs/bB-final/C3/l3945h/derived-FourierL3945H.lean)
  [pass]    Compare.lean positive: RESULT: PASS, 51 s, 16.0 GB (logs/bB-final/C3/l3945h/compare-pos.log)
  [pass]    Compare.lean negative: RESULT: FAIL on decimalExponentL3945H (one grid step higher, 1 - 17581521/10^10)
            (logs/bB-final/C3/l3945h/compare-neg.log)
  [pass]    official comparator positive (stand-in sandbox): nanoda and Lean's kernel accept; 1 h 05 min, 15.1 GB
            (l3945h/cmp/l3945h-pos.log) (logs/bB-final/C3/L3945H-cmp.out)
  [pass]    official comparator negative: rejected: Const does not match (decimalExponentL3945H); 3 min
            (logs/bB-final/C3/L3945H-cmp.out)
  [pass]    leaderboard verifier re-created offline (c92b9b6), merged stage: the import closure has 703 local modules,
            over the pull-request limit of 500; merged to 447 modules (449 files with Submission.lean and
            manifest.json), static check passed; Lean check from source 1 h 27 min (10.1 GB), axiom audit standard (1
            min), nanoda 1 h 01 min; the bound one grid step higher is rejected (logs/bB-final/C3/L3945H-bos.out)
  [pass]    merged stage = built tree, constant by constant: h-bos's fast identity check: the 29,512 non-auxiliary
            constants of the closure are identical in the merged stage and in the built tree, and so are the
            statements and proof terms of all 12,307 theorems (the adapter module aside) (logs/bB-final/C3/ident.out)
  [pass]    leaderboard verifier re-created offline, unmerged stage (C2, sources only): the 703-module form for the
            verifier's issue route (no file limit), compiled from sources only: Lean check 2 h 29 min, audit standard,
            nanoda 1 h 30 min; the bound one grid step higher is rejected (logs/bB-final/C2/bosU.out)
  [pass]    release controls: Heads3945HNeg and Rate3945HSharpCtl compile with their axioms; NegHeads, NegHeadsOld,
            NegRes, NegGlob fail (Lean's decide refutes each) (logs/bB-final/C3/ctl.out)
  [pass]    generators re-run (generators/REGENERATE.sh): 193 generated overlay files made again byte-identical; each
            generator refuses one step too far (logs/regenerate/regenerate.out)
  [pass]    the community's complete admission replay of #3945 on this fill:
            PASS_TERMINAL_WIDTH5_COMPLETE_CONDITIONAL_RELEASE with fresh native runs; 119 of 119 strict controls
            rejected (logs/ah-gate/GATE.md) (logs/ah-gate/frB/fullB.out)
  result.json lists each check with its log under logs/.

Controls (controls/, compiled against the built tree; logs/bB-final/C3/controls/)
  Heads3945HNeg.lean     Lean's kernel decides that the host table does NOT have 80,641 heads, and that NOT every
                         residual is at least 6 (the claims one step too far); compiles, axioms propext.
  Rate3945HSharpCtl.lean the rate fact is false at 1 - 17581521/10^10 and at every larger saving; compiles,
                         standard axioms.
  NegHeads.lean, NegHeadsOld.lean, NegRes.lean, NegGlob.lean   must FAIL to compile (80,641 heads; #3945's own
                         count of 57,312 on this table; a width bound of 6; one table chunk count one too high):
                         they do.
  side/                  four side checks built with the tree: Rate3945HSharp and RateMono (the rate control as two
                         modules), R5Joint3945H (the rate fact against the recursion's rule), S5HHB (the head bins
                         of the schedule are 80,640).
  side/perturbed/        the host table with a load moved (S5Hn) and with a window moved (S5Hw) as Lean parts: the
                         kernel refutes each moved piece, while its "= none" twin compiles (logs/bB-data/neg/).

The fill (fill/)
  FILL.json.gz (sha256 9c06677d...) and PROFILE.json are the final fill in #3945's format; old-FILL.json.gz
  (0e61fce2...) and old-PROFILE.json are the same allocation at the #3815 stage, which #3945's offset step reads.
  made-by/ holds the programs and the record of the run that made them (made-by/README.txt); checks/ holds the
  head counts and prices of the fill measured by two independent evaluators.

Generators
  generators/REGENERATE.sh makes every fill-dependent generated file of this package again from the community's #3945
  owner folder and fill/, and compares with SHA256SUMS; generators/README.txt says what makes what and what is the
  same as in the #3945 build that this package extends.

Credits
  The circuit: #3945 by Dugongue (staged544 T360 terminal-gap fill), on the lineage T74 words, #2979, #3133/#3148,
  #3605 (Ross Richardson), #3642, #3682 (zxv), #3734 (dataisfire: the staged544 words and owner interface),
  #3802/#3875 (hexapode) and #3815 (Erwan Secousses / ChaloupeMusic: the clock and the fill construction whose
  programs made this fill); the physical owner predicates of #3638, the warps and dirty borrowing of #2226 (Cinna
  Davis), Sagar Thakrar (#3687), roger-sato (#3693), Jeremy Kothe, Tekalign, icekylinx, DreamingOfClouds and
  chafreaky. Under it: Ryan Shea's Lean tree (generators, bank and seam pipeline), Nil Mamano's gauged Lean extension
  and record files, Dugongue's #401 and jacklightChen's #377 banks, and OpenAI's model, statement and Fourier
  reduction. Mine: the mixed-target recursion, the every-head-kept fill, the lever, and the Lean that connects them
  to the circuit.

Known cosmetic issues (left as compiled, because the files are checked byte for byte): the comment above the
Fourier theorem in each Main.lean calls it a sibling of a Walsh-Hadamard theorem wht_main_block_... that does
not exist (the correct text: it follows from the package's kernel program in Seam.lean); the header of each
Goal.lean cites OpenAI's lines 75-78, 80 and 82 where the definitions are at 75-77, 79 and 81. Line 15 of Work/FourierL3945H/Seam.lean says 57,312 heads; it should say 80,640.
