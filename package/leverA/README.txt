package/leverA: the recursion lever on the community's #2979 circuit with every head kept, saving 1.7564707e-3

What it is
  The new Lean files of this result: OpenAI's statement (the Fourier transform of every length, and convolution)
  with the exponent 1 - 17564707/10^10, that is a power saving of 1.7564707e-3. Laid over Ryan Shea's tree after
  package/u8x0hq (result.json, "requires"); every file is under Work/ or comparator/, so the lakefile does not
  change. Final theorems: OAI.PowerSaving.transform_mainL2979H and OAI.PowerSaving.convolution_mainL2979H, in
  Work.FourierL2979H.Main. Their axioms are propext, Classical.choice and Quot.sound.

The circuit, the fill and the lever
  #2979 is Marius-Paul Dumitrean's circuit (rank-7 joint end frames, on Jeremy Kothe's #2957 and Mason Hale's #2917
  words), merged on the community repository as research; its own figure is 1.7290340e-3. It runs 74 copies of a
  9-call pattern: 666 calls of two provider words, 66,304 arrays per gauge copy, 271,728 hosted roles.
  The fill says which call hosts each role. In #2979's allocation model a rank-48 bin that is left empty can sit in
  children 8/10, where it is a literal head (the first operation on its own array, acting on that array alone), or in
  children 12/13, where it never is; both cost the same, and #2979's solver split them arbitrarily (16,086 heads). I
  re-ran #2979's own fill construction with one extra constraint: every child-8/10 bin is empty (fill/README.txt).
  The new fill keeps the circuit's paid profile, idle mix and price exactly, and has 16,576 literal heads
  (8,288 + 8,288), the most possible. Only the fill changes: the words, the clock, the dispatch and the providers
  are #2979's, and the community's own checker and complete admission replay pass on it (fill/checks/).
  I cut each head by 24 of its 48 axes, so it becomes a block of rank 24 that leaves its array with target 49;
  every paid rank stays at most 48, below 49, so the mixed-target recursion applies with the two targets 73 and 49.
  The rate equation of the cut profile has its root at 1.75647080068e-3. The Lean figure 1.7564707e-3 is a grid point
  below it: the chain puts its rate point 0.02 grid steps above the figure, and 1.7564708e-3 lies only 0.0068 grid
  steps below the root, so its rate point would be above the root. For the same reason the rate fact is still true
  one grid step higher; it is false from two grid steps higher on (1 - 17564709/10^10, controls/Rate2979HSharpCtl.lean).
  The statement one grid step higher (1 - 17564708/10^10) is the negative challenge, which the checks reject.

What the files are (paths relative to lean/complex588)
  Work/Mixed/ (Flow, Phi, LemmaA, Bound, TwoTarget)   the mixed-target recursion: weights of the class tree, exact
                                                     balance, sources and exits, and the cost bound for every exponent
                                                     above the rate point
  Work/Mixed/Prog/, Work/MixedCost/                  the lockstep tree of calls as a program in OpenAI's machine
                                                     model, its correctness and its cost
  Work/Mixed/ (SeamConst, Rate2979H)                 the seam with a constant; the rate fact at the figure
  Work/Mixed/HeadM1                                  the head lemma and the half cut (24 + 24)
  Work/Mixed/Bridge/                                 from a heads-first network to the recursion's hypotheses
                                                     (Core, HeadsFirst, Inst; this instance's Inst2979H,
                                                     HeadsFirst2979H, HFTest2979H, HFTest2979HW)
  Work/KLowW/, Work/Warp/, Work/PWGen/               the warp checker, what a warp means, the word helpers, and the
                                                     pair of #2979 provider words (S2979Pair)
  Work/GCert/Data/ (S2979B*, S2979E*, Gen/S2979*)    both provider words as checked Lean data (generated)
  Work/XHostD/, Work/XHostR/,                        the hosted network with 74 copies and hosting across copies,
  Work/XHostT/ (T*, K*, HF/, HFW2979/)               its warp variant, the host-table checker, the heads-first
                                                     certificate on the warp network (Work/XHost/ itself comes with
                                                     package/u8x0hq)
  Work/XHostT/Data/H2979H/, Data/D2979/,             the host table of the new fill (all 271,728 roles), the dispatch,
  Inst/S2979H*                                       the schedule and its facts: 16,576 heads (generated data and
                                                     their glue)
  Work/L2979H/ (Stmt, EngineOfGeomW, Prov, Engine)   the engine statement, the network's data, and the one engine
                                                     theorem the chain imports
  Work/FourierL2979H/, Work/FourierL2979HNeg/        the Fourier chain at this exponent; the two challenge files are
                                                     OpenAI's statement with only the exponent and names changed
  Work/BlockAccounting/Sharp.lean                    this project's rate lemmas, needed by the rate files
  Work/Mixed/Rate3945.lean, Rate3815.lean            rate modules of the #3945 build that the generic bridge files
                                                     import; not used for this figure
  comparator/fourier-l2979h.json, -neg.json          the official comparator's configurations
  MODULES.txt lists all 284 modules in build order; ORIGINS.tsv says which part of the work wrote each file.

Checks
  All on Colab machines (8 processors, 50 GB), 19 October 2026; nothing was compiled elsewhere.
  - The whole closure builds with lake in a copy of the U8X0HQ stack (9,511 jobs, at most 7.7 GB), with no
    "declaration uses sorry" warning. #print axioms of the two final theorems and of seven key lemmas of the chain
    (the engine theorem, the seam, the program envelope, the rate point): propext, Classical.choice, Quot.sound.
  - Source greps over the 282 Lean files of the closure (the two challenge files hold OpenAI's own two sorry lines
    each, by design) and the four side checks: no sorry or admit, no axiom declaration, no native_decide,
    implemented_by, extern, unsafe or other kernel-bypassing construct.
  - The two challenge files are OpenAI's UniformFourier statement with only the exponent and the names changed: each
    is byte-identical to a fresh derivation from OpenAI's file, and they differ from each other only in the numerator.
  - tools/Compare.lean: the statement at 1 - 17564707/10^10 PASSES; one grid step higher, at 1 - 17564708/10^10, it
    FAILS ("decimalExponentL2979H: definition body differs").
  - The official comparator (nanoda and Lean's own kernel, in a jail as an unprivileged user) accepts the solution at
    17564707 ("Your solution is okay!", 47 minutes) and rejects it at 17564708 ("Const does not match between
    challenge and target 'OAI.PowerSaving.decimalExponentL2979H'").
  - The bigospeedrun leaderboard verifier, re-created offline (verifier clone at c92b9b6): the import closure is 581
    modules, so a pull request would hold 583 files, over the verifier's 500-file limit. The merge tool of the
    release (merging modules of one level of the import chain into one file, each original module verbatim exactly
    once) brings it to 447 modules, 449 files. The verifier compiles that merged stage from source and accepts it at
    the bound 17564707/10^10 (lean phase, 64 minutes; audit: Race.bound depends on propext, Classical.choice and
    Quot.sound only; nanoda accepts it too), and rejects the bound 17564708/10^10 (the adapter does not compile). A
    second machine built the same merged stage again: the same verdicts. An identity check compares every constant of
    the built tree with the merged stage: the 9,584 theorems have the same statements and identical proof terms, and
    both final theorems have the same axioms.
  - Controls (below) behave as stated; generators/REGENERATE.sh makes every generated file again byte for byte.
  - The fill: #2979's own allocation checker and its complete admission replay (fresh native runs) pass on it
    (fill/README.txt).
  result.json lists each check with its log under logs/.

Controls (controls/, compiled against the built tree; logs/2026-10-19-build/controls/)
  Heads2979HNeg.lean      Lean's kernel decides that the host table does NOT have 16,577 heads, and that NOT every
                          residual is at least 5 (the claims one step too far); compiles, axioms propext.
  Rate2979HSharpCtl.lean  the rate fact is false at 1 - 17564709/10^10 and at every larger saving; compiles,
                          standard axioms.
  NegHeads.lean, NegGlob.lean   must FAIL to compile (16,577 heads; one table chunk count one too high): they do.
  side/                   four side checks built with the tree: Rate2979HSharp, RateMono, R5Joint2979H (the rate
                          fact against the recursion's rule), S2979HHB (the head bins of the schedule are 16,576).

Generators
  generators/REGENERATE.sh makes every generated file of this package again from the community's #2979 files and
  the new fill, and compares with SHA256SUMS and controls/; generators/README.txt says what makes what and what is
  hand-written. fill/ holds the fill and the commands that made it.

Credits
  The circuit: #2979 by Marius-Paul Dumitrean, on Jeremy Kothe's #2957/#2948/#2966 rewrites and Mason Hale's #2917
  words, T74 chronology and allocation checker; with Jeremy Diviney (#2908), Dugongue (#2969, #2863/#2839),
  RossRichardson (#2918), Yiyang Lu (#2857/#2848), Tekalign, icekylinx, themilkman216 and Gareth Davidson
  (#1995/#2013); the warps and dirty borrowing of #2226 (Cinna Davis) and the helper width exchange of #2312
  (Jeremy Kothe). The fill is made with #2979's own fill17.py and double.py. Under it: Ryan Shea's Lean tree
  (generators, bank and seam pipeline), Nil Mamano's gauged Lean extension and record files, Dugongue's #401 and
  jacklightChen's #377 banks, and OpenAI's model, statement and Fourier reduction. Mine: the mixed-target recursion,
  the every-head-kept fill, the lever, and the Lean that connects them to the circuit.

Known cosmetic issues (left as compiled, because the files are checked byte for byte): the comment above the
Fourier theorem in each Main.lean calls it a sibling of a Walsh-Hadamard theorem wht_main_block_... that does
not exist (the correct text: it follows from the package's kernel program in Seam.lean); the header of each
Goal.lean cites OpenAI's lines 75-78, 80 and 82 where the definitions are at 75-77, 79 and 81.
