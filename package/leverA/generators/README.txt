generators/: how the generated Lean files of this package are made again

Run (Python 3, standard library only; about 3.5 minutes and under 0.3 GB on a laptop):

    bash generators/REGENERATE.sh /some/new/folder

It runs every program below from tree/ on the community's #2979 files (tree/wht45/lever/p-alt/dl/main2979/
research/rank7-joint-endframes/, read from the community repository's main branch at bb4bccbf), the two vendor pins
TEMPLATE.json.gz and PLAN.json.gz, and design A's fill (tree/wht45/allheads/aA-data/in/alloc2979H/FILL.json.gz, the
same file as ../fill/FILL.json.gz). It writes everything into the new folder and compares the result byte for byte
with ../SHA256SUMS and with the controls in ../controls/. Every generator that writes a figure or a count is also run
one unit too high and must refuse. On 19 Oct 2026 it ended in "REGENERATE: ALL AS EXPECTED": 185 overlay files came
out byte-identical (written by the generators, or copied in by the word-stage tool), none differed, and the four
controls and three side checks came out byte-identical too. The record of that run is ../logs/regenerate/.

tree/ holds the programs at the same relative paths as in my workspace, so each runs unchanged. The changes are one
path-only patch (absolute paths become $GEN_ROOT, which REGENERATE.sh sets to tree/) and one sha256 pin that follows
it (aA-data's mk_rate2979H.py pins e-r45m1's mk_rate_mixed.py, whose path line the patch changes); PATCHES.txt shows
both line by line, with the hashes before and after. tree.sha256 lists every file of tree/.

What makes what (paths of the outputs relative to lean/complex588)

  the two provider words (gauge44 = WORD-base, extra = WORD-extra), as checked Lean data: as for the #2979 lever
  package (the words do not depend on the fill)
    p-lower/tools/transport.py, lower_w.py      the community's words in our lowering (lowered words byte-identical
                                                to the ones the Lean build used: 8f2db138..., 163466d5...)
    x-words/tools/mk_info.py, pack_ceil2.py,     Work/GCert/Data/Gen/S2979B*, S2979E*, Work/GCert/Data/S2979*Glue401,
    wgen2.py, gen_global_w.py, glue2.lean.tmpl,  S2979*Price, InstS2979*, RateS2979*x; refuses the standalone figure
    k-combine/tools/make_sched_x.py,             one unit higher (--dft 14611616)
    generate_checks_x.py
    x-words/tools/mk_stage.py, mk_stage_pair.py  the word stage (generated files plus the hand-written checker
                                                Work/KLowW, warp meaning Work/Warp, helpers Work/PWGen and the pair
                                                module Work/PWGen/S2979Pair.lean from tree/.../x-words/stage/)
  the host table of design A (who hosts every one of the 271,728 roles), its Lean parts and the dispatch
    x-tables/tools/mk_resid2979.py               our lowering's residuals (= x-tables', checked)
    x-tables/tools/own_to_host2979.py            the table of #2979's PUBLISHED fill (only for the comparison below)
    aA-data/tools/own_to_host2979H.py            the table of design A's fill (x-tables' converter with only the
                                                FILL pin changed); full profile = #2979's PROFILE.json
    aA-final/tools/fill_indep2979H.py            the comparison of the two tables: every call field (provider, time,
                                                order, windows, residuals) is the same; only role rows, offsets and
                                                remainders change
    aA-data/tools/heads2979.py                   the 16,576 literal heads (8,288 + 8,288); one head more must fail
    aA-data/tools/mk_parts.py                    Work/XHostT/Data/H2979H/{Calls, Order, P0..P26, Whole, GlobC0..5,
                                                Glob} and controls/NegGlob.lean; refuses a perturbed table
    x-tables/tools/mk_disp2979.py                Work/XHostT/Data/D2979/Disp.lean (416 extra, 250 gauge44 calls; the
                                                same file as for #2979's own fill); refuses 417 / 249
  the figure, the instance modules, the statement chain, the warp heads-first files, the engine theorem
    aA-data/tools/mk_profile2979H.py             the cut profile (16,576 heads cut to rank 24; 11,614 other rank-48)
    aA-data/tools/mk_rate2979H.py (e-r45m1's     Work/Mixed/Rate2979H.lean: the rate fact at 1 - 17564707/10^10;
    mk_rate_mixed.py, imported unchanged)        refuses 17564708; and controls/side/.../Rate2979HSharp.lean (the
                                                rate fact false at 1 - 17564709/10^10, see ../README.txt)
    aA-data/tools/mk_inst2979H.py                a CHECKED renaming of the hand-written #2979 instance modules
                                                (Data/H2979H/{Heads, DispLink}, Inst/S2979H*, Mixed/Bridge/*2979H,
                                                Mixed/R5Joint2979H, L2979H/EngineOfGeomW, controls/NegHeads.lean):
                                                it takes the new numbers from the data above and refuses unless the
                                                remainder histogram, residual sums, dispatch and call fields are
                                                #2979's
    aA-data/tools/mk_l2979H.py (with             Work/FourierL2979H/*, Work/FourierL2979HNeg/*, Work/L2979H/Stmt.lean,
    mk_fourier_mixed.py)                         Work/Mixed/SeamConst.lean, comparator/*.json; the two challenge files
                                                are OpenAI's statement with only the exponent and names changed;
                                                refuses --num 17564708
    x-inst/tools/mk_hfw2979.py                   Work/XHostT/HFW2979/{Def, Stage, Cert, ToHF, Heads}.lean (the
                                                heads-first certificate on the warp network; fill-independent)
    aA-final/tools/mk_prov2979H.py               Work/L2979H/{Prov, Engine}.lean: the #2979 network's data and the one
                                                engine theorem, from the #2979 lever package's two files with the
                                                design-A names (refuses any fill-dependent name left unrenamed)
    aA-final/tools/mk_controls2979H.py           controls/{Heads2979HNeg, Rate2979HSharpCtl}.lean

Not made again here (../logs/regenerate/compare.tsv lists them, 101 files):
  - hand-written Lean: the recursion engine (Work/Mixed/{Flow, Phi, LemmaA, Bound, TwoTarget}, Work/Mixed/Prog/*,
    Work/MixedCost/*), the generic bridge (Work/Mixed/Bridge/{Core, HeadsFirst, Inst}), the head lemma
    (Work/Mixed/HeadM1), the table glue and counters (Work/XHostT/K*, TEnd), the #2979 instance layer that does not
    depend on the fill (Work/XHostT/HF/{Count2979C, Inst2979}, HFW2979/{GeomA4W, GeomInstW, Disp2979W});
  - generic generated files of the network that do not depend on the word or the fill (Work/XHostD/*,
    Work/XHostT/T*.lean and T*W.lean, Work/XHostT/HF/{Def, Stage, Cert, ToHF, Heads, GeomA4, GeomInst}): made by
    the hosted-network generators of the earlier #3945 build and used here byte for byte;
  - Work/Mixed/Rate3945.lean and Rate3815.lean: rate modules of the #3945 build that the generic bridge files
    import; not used by this result's figure.

The fill itself is not made again here: it comes from a CP-SAT search (OR-Tools, about 5 minutes); ../fill/README.txt
has the commands and the record of that run.
