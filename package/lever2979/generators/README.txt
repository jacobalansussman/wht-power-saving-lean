generators/: how the generated Lean files of this package are made again

Run (Python 3, standard library only; about 3 minutes and under 0.3 GB on a laptop):

    bash generators/REGENERATE.sh /some/new/folder

It runs every program below from tree/ on the community's #2979 files (tree/wht45/lever/p-alt/dl/main2979/
research/rank7-joint-endframes/, read from the community repository's main branch at bb4bccbf) and the two vendor
pins TEMPLATE.json.gz and PLAN.json.gz, writes everything into the new folder, and compares the result byte for byte
with ../SHA256SUMS. Every generator that writes a figure or a count is also run one unit too high and must refuse.
On 19 Oct 2026 it ended in "REGENERATE: ALL AS EXPECTED": 174 overlay files came out byte-identical (116 written by
the generators, 58 hand-written checker and helper files that the word-stage tool copies in), none differed.
The record of that run is ../logs/regenerate/.

tree/ holds the programs at the same relative paths as in my workspace, so each runs unchanged. The only change is
one path-only patch (absolute paths become $GEN_ROOT, which REGENERATE.sh sets to tree/); PATCHES.txt shows it line
by line, with the hashes before and after. tree.sha256 lists every file of tree/.

What makes what (paths of the outputs relative to lean/complex588):

  the two provider words (gauge44 = WORD-base, extra = WORD-extra), as checked Lean data
    p-lower/tools/transport.py, lower_w.py      the community's words in our lowering (lowered words byte-identical
                                                to the ones the Lean build used: 8f2db138..., 163466d5...)
    x-words/tools/mk_info.py, pack_ceil2.py,     Work/GCert/Data/Gen/S2979B*, S2979E* (labels, scalars, histogram,
    wgen2.py (with wgen.py), gen_global_w.py,    residuals, windows, fingerprints), Work/GCert/Data/S2979*Glue401,
    glue2.lean.tmpl, k-combine/tools/            S2979*Price, InstS2979*, RateS2979*x; refuses the standalone figure
    make_sched_x.py, generate_checks_x.py        one unit higher (--dft 14611616)
    x-words/tools/mk_stage.py, mk_stage_pair.py  the word stage: the generated files plus the hand-written checker
                                                (Work/KLowW), warp meaning (Work/Warp), helpers (Work/PWGen) and the
                                                pair module Work/PWGen/S2979Pair.lean from tree/.../x-words/stage/
  the host table (who hosts every one of the 271,728 roles), its Lean parts and the dispatch
    x-tables/tools/mk_resid2979.py,              the table from the community's allocation (FILL, ORDER, DISPATCH,
    own_to_host2979.py                          PROFILE, PROVIDERS) and our lowering's residuals
    x-tables/tools/heads2979.py                  the 16,086 literal heads (8,064 + 8,022); one head more must fail
    x-tables/tools/mk_parts.py                   Work/XHostT/Data/H2979/{Calls, Order, P0..P26, Whole, GlobC0..5,
                                                Glob}; refuses a perturbed table
    x-tables/tools/mk_disp2979.py                Work/XHostT/Data/D2979/Disp.lean (416 extra, 250 gauge44 calls);
                                                refuses 417 / 249
  the figure, the statement chain, the warp heads-first files
    e-r45m1/tools/mk_rate_mixed.py               Work/Mixed/Rate2979.lean (and Rate2979Sharp.lean, shipped in
                                                ../controls/): the rate fact at 1 - 17556473/10^10; refuses 17556474
    x-inst/tools/mk_l2979.py (with               Work/FourierL2979/*, Work/FourierL2979Neg/*, Work/L2979/Stmt.lean,
    mk_fourier_mixed.py)                         Work/Mixed/SeamConst.lean, comparator/*.json; the two challenge files
                                                are OpenAI's statement with only the exponent and names changed;
                                                refuses --num 17556474
    x-inst/tools/mk_hfw2979.py                   Work/XHostT/HFW2979/{Def, Stage, Cert, ToHF, Heads}.lean: the
                                                heads-first certificate moved to the network with warp providers
                                                (inputs pinned in x-inst/tools/hfw2979-inputs.sha256)

Not made again here (../logs/regenerate/compare.tsv lists them, 112 files):
  - hand-written Lean: the recursion engine (Work/Mixed/{Flow, Phi, LemmaA, Bound, TwoTarget}, Work/Mixed/Prog/*,
    Work/MixedCost/*), the bridge (Work/Mixed/Bridge/*), the head lemma (Work/Mixed/HeadM1), the table glue and
    counters (Work/XHostT/K*, Data/H2979/{Heads, DispLink}, Inst/S2979*), the #2979 instance layer
    (Work/XHostT/HF/{Count2979C, Inst2979}, HFW2979/{GeomA4W, GeomInstW, Disp2979W}, Work/L2979/*);
  - generic generated files of the network that do not depend on the word (Work/XHostD/*, Work/XHostT/T*.lean and
    T*W.lean, Work/XHostT/HF/{Def, Stage, Cert, ToHF}): made by the hosted-network generators of the earlier
    #3945 build (gen_xnet.py, gen_tgeom.py, mk_provw.py, mk_warp.py, mk_xneth.py, mk_framesk.py, gen_hnet.py) and
    used here byte for byte;
  - Work/Mixed/Rate3945.lean and Rate3815.lean: rate modules of the #3945 build that the generic bridge files
    import; not used by this result's figure.
