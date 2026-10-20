generators/: how the generated Lean files of this package are made again

Run (Python 3.11 or earlier, standard library only; a few minutes and under 0.5 GB on a laptop):

    bash generators/REGENERATE.sh /some/new/folder

It runs every program below from tree/ on the community's #2979 files (tree/wht45/lever/p-alt/dl/main2979/
research/rank7-joint-endframes/, read from the community repository's main branch at bb4bccbf), the two vendor pins
TEMPLATE.json.gz and PLAN.json.gz, and this package's fill (tree/wht45/afc/afc-data/in/alloc2979Fc/FILL.json.gz and
PROFILE.json, the same files as ../fill/FILL.json.gz and ../fill/PROFILE.json; ORDER, DISPATCH and PROVIDERS are
#2979's own bytes). It writes everything into the new folder and compares the result byte for byte with ../SHA256SUMS
and with the controls in ../controls/. Every generator that writes a figure, a count or a table is also run one step
too far and must refuse. The record of the run on 20 Oct 2026 is ../logs/regenerate/.

Python 3.11 or earlier. x-words/tools/pack_ceil2.py chooses a bank pattern by the smallest sum of squares of floats,
and Python 3.12 changed how sum() adds floats (compensated summation), so near-ties break differently: on Python 3.12
or later, 36 row tables of the two provider words (Work/GCert/Data/Gen/S2979BRows0..17.lean and S2979ERows0..17.lean)
come out different and the script ends "REGENERATE: SOMETHING FAILED". Those files are the ones Lean compiled and
checked; only their regeneration depends on the Python version. With Python 3.12 or later installed, make a virtual
environment with an older Python and put it first on PATH, for example:

    python3.11 -m venv /tmp/py311 && PATH=/tmp/py311/bin:$PATH bash generators/REGENERATE.sh /some/new/folder

The two lowering reports in tree/wht45/lever2979/x-gate/out/ (lower-base.json, lower-extra.json) name the machine
paths of the run that made them; those paths were shortened for publication. tree.sha256 gives their published bytes,
and REGENERATE.sh puts their original sha256 (f219f479..., 5e5d7945...) back into the residual file it makes, since
that file records them; PATCHES.txt says the same.

tree/ holds the programs at the same relative paths as in my workspace, so each runs unchanged. It is package/leverA's
generators tree (the same files; leverA's PATCHES.txt is the first part of PATCHES.txt here) with this package's
programs added. The changes are the path-only patch (absolute paths become $GEN_ROOT, which REGENERATE.sh sets to
tree/) and the sha256 pins that follow it (here: mk_l2979Fc.py pins its copy of mk_rate_mixed.py, whose path line the
patch changes); PATCHES.txt shows each change line by line, with the hashes before and after. tree.sha256 lists every
file of tree/.

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
  the host table of this fill (who hosts every one of the 271,728 roles), its Lean parts and the dispatch
    x-tables/tools/mk_resid2979.py               our lowering's residuals (= x-tables', checked)
    x-tables/tools/own_to_host2979.py            the table of #2979's PUBLISHED fill (only for the comparison below)
    afc-data/tools/own_to_host2979Fc.py          the table of this fill (x-tables' converter with only the FILL and
                                                PROFILE pins changed); its sha256 is checked (a4baa4e6...)
    afc-data/tools/fill_indep2979Fc.py           the comparison of the two tables: every call field (provider, time,
                                                order, windows, residuals) is the same; only role rows, offsets and
                                                remainders change
    afc-data/tools/heads2979Fc.py                the 16,576 literal heads (8,288 + 8,288) and no other rank-48 block
                                                (x-tables' heads2979.py; one comparison fixed so that an expected
                                                count of 0 is accepted); one head more, or one other rank-48 block,
                                                must fail
    afc-data/tools/mk_parts.py                   (p-hostck's program, unchanged) Work/XHostT/Data/H2979Fc/{Calls,
                                                Order, P0..P26, Whole, GlobC0..5, Glob} and controls/NegGlob.lean;
                                                refuses a table with a load moved or a window moved; with --emit-neg
                                                the perturbed tables (controls/side/.../H2979Fcn*/Calls.lean) and the
                                                controls/perturbed/ files
    afc-data/tools/mk_disp2979.py                Work/XHostT/Data/D2979/Disp.lean (416 extra, 250 gauge44 calls; the
                                                same file as for #2979's own fill); refuses 417 / 249
  the instance modules, the figure, the statement chain, the warp heads-first files, the engine's data, the controls
    afc-data/tools/mk_inst2979Fc.py              a CHECKED renaming of the hand-written #2979 instance modules
                                                (Data/H2979Fc/{Heads, DispLink}, Inst/S2979Fc, S2979FcFacts, and the
                                                side checks S2979FcHB, S2979FcHFW, S2979FcUnit; controls/NegHeads.lean,
                                                NegF1.lean): it takes the new numbers (heads, idle remainders, unit
                                                profile) from the data above and refuses unless 16,576 heads, no other
                                                rank-48 block, 666 calls x the local word + the remainders = PROFILE,
                                                and residual sums, dispatch and call fields are #2979's
    afc-cut/tools/mk_profile2979Fc.py            the full-cut profile (the 16,576 heads removed, target 25; nothing
                                                split; refuses any other profile, and any rank 24..47 or above 48)
    afc-cut/tools/mk_rate_mixed.py               (e-r45m1's generator; only the header text option differs)
                                                Work/Mixed/Rate2979Fc.lean: the rate fact at 1 - 17887045/10^10;
                                                refuses 17887046; and controls/side/.../Rate2979FcSharp.lean (the rate
                                                fact false one grid step higher)
    afc-cut/tools/mk_l2979Fc.py (with            Work/FourierL2979Fc/*, Work/FourierL2979FcNeg/*, Work/L2979Fc/Stmt.lean,
    mk_fourier_mixed.py)                         Work/Mixed/SeamConst.lean, comparator/*.json; the two challenge files
                                                are OpenAI's statement with only the exponent and names changed;
                                                refuses --num 17887046 (final and hypothesis modes)
    x-inst/tools/mk_hfw2979.py                   Work/XHostT/HFW2979/{Def, Stage, Cert, ToHF, Heads}.lean (the
                                                heads-first certificate on the warp network; fill-independent)
    afc-cut/tools/mk_prov2979Fc.py               Work/L2979Fc/Prov.lean: the #2979 network's data, from package/
                                                leverA's Prov.lean with this fill's names and idle price (read from the
                                                regenerated S2979FcFacts.lean and checked against PROFILE.json)
    afc-final/tools/mk_controls2979Fc.py         controls/{Heads2979FcNeg, Rate2979FcSharpCtl}.lean and copies of the
                                                other controls (Ctl2979Fc, NegHeads, NegGlob, NegF1, perturbed/)

Not made again here (../logs/regenerate/compare.tsv lists them, 103 files):
  - hand-written Lean: this project's rate lemmas (Work/BlockAccounting/Sharp), the recursion engine (Work/Mixed/{Flow, Phi, LemmaA, Bound, TwoTarget}, Work/Mixed/Prog/*,
    Work/MixedCost/*), the generic bridges (Work/Mixed/Bridge/{Core, CoreF}), the head lemma (Work/Mixed/HeadM1), this
    instance's bridge and bookkeeping (Work/Mixed/Bridge/{Inst2979Fc, HFTest2979FcW}, Work/Mixed/Hist2979Fc,
    Work/L2979Fc/{EngineOfGeomW, Engine}), the table glue and counters (Work/XHostT/K*, TEnd), the #2979 instance
    layer that does not depend on the fill (Work/XHostT/HF/{Count2979C, Inst2979}, HFW2979/{GeomA4W, GeomInstW,
    Disp2979W}); the control Ctl2979Fc;
  - generic generated files of the network that do not depend on the word or the fill (Work/XHostD/*, Work/XHostR/*,
    Work/XHostT/T*.lean and T*W.lean, Work/XHostT/HF/{Def, Stage, Cert, ToHF, Heads, GeomA4, GeomInst}): made by
    the hosted-network generators of the earlier #3945 build and used here byte for byte.

The fill itself is not made again here: it comes from an LP and MIP solve and #3815's placement program on a cloud
machine (about 1 minute); ../fill/made-by/ has the programs, the class solution and the record of that run.
