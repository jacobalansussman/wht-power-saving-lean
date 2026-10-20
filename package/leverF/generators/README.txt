generators/: how the generated Lean files of this package are made again

Run (Python 3, standard library only; about 4 minutes and under 1 GB on a laptop):

    bash generators/REGENERATE.sh /some/new/folder

It lays the design's fill (../fill/FILL.json.gz, old-FILL.json.gz, PROFILE.json) over the community's #3945 owner
folder (tree/wht45/lever/g-heads/dl/pr3945/owner, read from the community repository at PR head 91f81d13) and checks
the result against owner-design.sha256, runs every program below from tree/, writes everything into the new folder,
and compares the result byte for byte with ../SHA256SUMS and with ../controls/. The fill and every fill-dependent
file are the ones of the half-cut package on the same fill (leverB); only the rate fact and the statement chain
belong to the full cut. Every generator that writes a figure,
a count or a table is also run one step too far and must refuse. The record of the run of 19 Oct 2026 is
../logs/regenerate/.

tree/ holds the programs at the same relative paths as in my workspace, so each runs unchanged. The only change is a
path-only patch (absolute paths become $GEN_ROOT, which REGENERATE.sh sets to tree/); PATCHES.txt shows it line by
line, with the hashes before and after. tree.sha256 lists every file of tree/.

What makes what (paths of the outputs relative to lean/complex588)

  the head gate in my lowering
    bB-data/tools/g0_staged.py               (g-heads' program, unchanged) runs the fill through my hosted network
                                             of the two #3945 words: 80,640 literal heads (40,320 + 40,320), the
                                             other 57,312 rank-48 blocks not heads, no window violation, paid
                                             profile = #3945's PROFILE.json
  the host table of the fill (who hosts each of the 1,321,920 roles) and its Lean parts
    bB-data/tools/own_to_host3945H.py        p-own's converter (the fill's pin and its expected count of roles
                                             hosted across copies, 1,161,423, are the only changes); the table's
                                             sha256 is checked
    bB-data/tools/heads3945H.py              p-own's head list with the design's expected counts (80,640 heads,
                                             57,312 others), cross-checked three ways
    bB-data/tools/mk_parts.py                (p-hostck's program, unchanged) Work/XHostT/Data/S5H/{Calls, Order,
                                             P0..P129, Whole, GlobC0..26, Glob}; refuses a table with a load moved
                                             or a window moved; its NegGlob.lean is ../controls/NegGlob.lean
  the dispatch, which does not depend on the fill
    lever/p-disp/tools/mk_disp3945.py        Work/XHostT/Data/D3945/Disp.lean, byte-identical to the earlier
                                             #3945 build's (1,672 extra and 1,568 gauge44 calls); refuses 1,673 /
                                             1,567
  the figure
    fullcut/fc-rate/tools/mk_rate_mixed.py   (e-r45m1's program; one change: a profile may carry a description,
                                             so the Lean documentation says "full cut") Work/Mixed/Rate3945F.lean:
                                             the rate fact at 1 - 17806935/10^10 for the full-cut profile
                                             (profile3945F.json); refuses 17806936; its Rate3945FSharp.lean is the
                                             second half of ../controls/Rate3945FSharpCtl.lean
    bB-data/tools/mk_rate_mixed.py           (e-r45m1's program, unchanged) Work/Mixed/Rate3945H.lean, the half
                                             cut's rate fact: the generic bridge files of this closure import it;
                                             it is not used for this figure; refuses 17581521
  the hand-written files that name the fill's numbers
    bB-data/tools/mk_3945H.py                Work/XHostT/Data/S5H/{Heads, DispLink}, Work/XHostT/Inst/{S5H,
                                             S5HFacts}, Work/Mixed/Bridge/{Inst3945H, HeadsFirst3945H, HFTest3945HW},
                                             Work/L3945H/{EngineOfGeomW, Prov, Engine} (only Prov is imported here;
                                             the other two are compared with the half-cut package's
                                             designB-only.sha256): the files of the earlier
                                             #3945 build (and two of the #2979 build) that compiled there, with
                                             exact text replacements, each of which must occur the stated number of
                                             times; any leftover old count or name makes it refuse. Its negatives
                                             are ../controls/NegHeads.lean, NegHeadsOld.lean, NegRes.lean
  the statement chain
    fullcut/fc-engine/tools/mk_l3945.py      (i-skel's program; one change: --word replaces the description of the
    (with mk_fourier_mixed.py,               cut word in one comment line of the generated Seam.lean) Work/FourierL3945F/*,
    profile3945F.json, word3945F.txt)        Work/FourierL3945FNeg/*, Work/L3945F/Stmt.lean, Work/Mixed/SeamConst.lean,
                                             comparator/*.json; the two challenge files are OpenAI's statement with
                                             only the exponent and names changed; refuses --num 17806936

Not made again here (the comparison lists them, about 220 files): everything that does not depend on the fill, used
byte for byte as the earlier #3945 build compiled it (ORIGINS.tsv names the part of the work that wrote each file):
  - hand-written Lean of the full cut: Work/Mixed/Bridge/{CoreF, Inst3945F, HFTest3945FW}, Work/L3945F/{EngineOfGeomW,
    Engine} (fc-engine);
  - hand-written Lean: the recursion engine (Work/Mixed/{Flow, Phi, LemmaA, Bound, TwoTarget}, Work/Mixed/Prog/*,
    Work/MixedCost/*), the generic bridge (Work/Mixed/Bridge/{Core, HeadsFirst, Inst, HFTest}), the head lemma
    (Work/Mixed/HeadM1), the table glue and counters (Work/XHostT/K*), the warp heads-first certificate
    (Work/XHostT/HFW2979/*, generic in the instance despite its name);
  - the two #3945 provider words as Lean data (Work/GCert/Data/Gen/S5*, the word stage of p-wgen) and the generic
    generated files of the hosted network (Work/XHostD/*, Work/XHostT/T*.lean, Work/XHostT/HF/*): made by the word
    and network generators of the earlier #3945 build from #3945's words, which this fill does not change.
