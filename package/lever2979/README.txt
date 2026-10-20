package/lever2979: the recursion lever on the community's #2979 circuit, saving 1.7556473e-3

What it is
  The new Lean files of this result: OpenAI's statement (the Fourier transform of every length, and convolution)
  with the exponent 1 - 17556473/10^10, that is a power saving of 1.7556473e-3. Laid over Ryan Shea's tree after
  package/u8x0hq (result.json, "requires"); every file is under Work/ or comparator/, so the lakefile does not
  change. Final theorems: OAI.PowerSaving.transform_mainL2979 and OAI.PowerSaving.convolution_mainL2979, in
  Work.FourierL2979.Main. Their axioms are propext, Classical.choice and Quot.sound.

The circuit and the lever
  #2979 is Marius-Paul Dumitrean's circuit (rank-7 joint end frames, on Jeremy Kothe's #2957 and Mason Hale's #2917
  words), merged on the community repository as research; its own figure is 1.7290340e-3. It runs 74 copies of a
  9-call pattern: 666 calls of two provider words, 66,304 arrays per gauge copy, 271,728 hosted roles. In my
  lowering 16,086 of its rank-48 blocks are literal heads: the first operation on their own array, acting on that
  array alone. I cut each head by 24 of its 48 axes, so it becomes a block of rank 24 that leaves its array with
  target 49; every paid rank stays at most 48, below 49, so the mixed-target recursion applies with the two targets
  73 and 49. The rate equation of the cut profile has its root at 1.75564731529e-3; the Lean figure 1.7556473e-3 is
  the grid point just below it, and the rate fact is false one grid step higher (controls/Rate2979SharpCtl.lean).

What the files are (paths relative to lean/complex588)
  Work/Mixed/ (Flow, Phi, LemmaA, Bound, TwoTarget)   the mixed-target recursion: weights of the class tree, exact
                                                     balance, sources and exits, and the cost bound for every exponent
                                                     above the rate point
  Work/Mixed/Prog/, Work/MixedCost/                  the lockstep tree of calls as a program in OpenAI's machine
                                                     model, its correctness and its cost
  Work/Mixed/ (SeamConst, Rate2979)                  the seam with a constant; the rate fact at the figure
  Work/Mixed/HeadM1                                  the head lemma and the half cut (24 + 24)
  Work/Mixed/Bridge/                                 from a heads-first network to the recursion's hypotheses
                                                     (Core, HeadsFirst, Inst; the #2979 instance Inst2979,
                                                     HeadsFirst2979, HFTest2979, HFTest2979W)
  Work/KLowW/, Work/Warp/, Work/PWGen/               the warp checker, what a warp means, the word helpers, and the
                                                     pair of #2979 provider words (S2979Pair)
  Work/GCert/Data/ (S2979B*, S2979E*, Gen/S2979*)    both provider words as checked Lean data (generated)
  Work/XHostD/, Work/XHostR/,                        the hosted network with 74 copies and hosting across copies,
  Work/XHostT/ (T*, K*, HF/, HFW2979/)               its warp variant, the host-table checker, the heads-first
                                                     certificate on the warp network (Work/XHost/ itself comes with
                                                     package/u8x0hq)
  Work/XHostT/Data/H2979/, Data/D2979/, Inst/S2979*  the host table (all 271,728 roles), the dispatch, the schedule
                                                     and its facts: 16,086 heads, the unit price (generated data and
                                                     their glue)
  Work/L2979/ (Stmt, EngineOfGeomW, Prov, Engine)    the engine statement, the network's data, and the one engine
                                                     theorem the chain imports
  Work/FourierL2979/, Work/FourierL2979Neg/          the Fourier chain at this exponent; the two challenge files are
                                                     OpenAI's statement with only the exponent and names changed
  Work/BlockAccounting/Sharp.lean                    this project's rate lemmas, needed by the rate files
  Work/Mixed/Rate3945.lean, Rate3815.lean            rate modules of the #3945 build that the generic bridge files
                                                     import; not used for this figure
  comparator/fourier-l2979.json, -neg.json           the official comparator's configurations
  MODULES.txt lists all 284 modules in build order; ORIGINS.tsv says which part of the work wrote each file.

Checks
  All of these ran on 19 October 2026 on Google Colab machines (8 processors, 50 GB); the generators' replay ran on a
  laptop.
  Build              the whole closure builds in a copy of the U8X0HQ stack (9,511 jobs, 20 minutes, at most 9.5 GB).
  Axioms             transform_mainL2979, convolution_mainL2979, the engine theorem and the seam depend on propext,
                     Classical.choice and Quot.sound only.
  Source greps       no sorry or admit, and no axiom, native_decide, unsafe, implemented_by or extern, in the 282 new
                     Lean files (OpenAI's two challenge statements, which contain sorry by design, are not counted).
  Statement          both challenge files are OpenAI's UniformFourier statement with only the exponent and the names
                     changed, byte for byte.
  Compare.lean       PASS at 1 - 17556473/10^10; FAIL one grid step higher, at 1 - 17556474/10^10.
  Official comparator  accepted (nanoda and Lean kernels); one grid step higher rejected.
  bigospeedrun       the verifier re-created offline (commit c92b9b6, profile direct). The proof needs 581 modules,
                     583 files in a pull request, over the verifier's 500-file limit; merged into 447 modules
                     (449 files), it is accepted by all three phases (Lean kernel, axiom audit, nanoda) and rejected
                     one grid step higher. A check of the merge finds the same 9,584 theorems, with identical proof
                     terms, in the merged stage and in the tree checked above. The unmerged stage also passes all
                     three phases.
  result.json lists each check with its log under logs/.

Controls (controls/, compiled against the built tree; logs/controls/)
  Heads2979Neg.lean      Lean's kernel decides that the host table does NOT have 16,087 heads, and that NOT every
                         residual is at least 5 (the claims one step too far); compiles, axioms propext.
  Rate2979SharpCtl.lean  the rate fact is false at 1 - 17556474/10^10 and at every larger saving; compiles,
                         standard axioms.
  NegHeads.lean, NegGlob.lean   must FAIL to compile (16,087 heads; one table chunk count one too high): they do.
  side/                  four side checks built with the tree: Rate2979Sharp, RateMono, R5Joint2979 (the rate fact
                         against the recursion's rule), S2979HB (the head bins of the schedule are 16,086).

Generators
  generators/REGENERATE.sh makes every generated file of this package again from the community's #2979 files and
  compares with SHA256SUMS; generators/README.txt says what makes what and what is hand-written.

Credits
  The circuit: #2979 by Marius-Paul Dumitrean, on Jeremy Kothe's #2957/#2948/#2966 rewrites and Mason Hale's #2917
  words, T74 chronology and allocation checker; with Jeremy Diviney (#2908), Dugongue (#2969, #2863/#2839),
  RossRichardson (#2918), Yiyang Lu (#2857/#2848), Tekalign, icekylinx, themilkman216 and Gareth Davidson
  (#1995/#2013); the warps and dirty borrowing of #2226 (Cinna Davis) and the helper width exchange of #2312
  (Jeremy Kothe). Under it: Ryan Shea's Lean tree (generators, bank and seam pipeline), Nil Mamano's gauged Lean
  extension and record files, Dugongue's #401 and jacklightChen's #377 banks, and OpenAI's model, statement and
  Fourier reduction. Mine: the mixed-target recursion, the lever, and the Lean that connects them to the circuit.

Known cosmetic issues (left as compiled, because the files are checked byte for byte): the comment above the
Fourier theorem in each Main.lean calls it a sibling of a Walsh-Hadamard theorem wht_main_block_... that does
not exist (the correct text: it follows from the package's kernel program in Seam.lean); the header of each
Goal.lean cites OpenAI's lines 75-78, 80 and 82 where the definitions are at 75-77, 79 and 81.
