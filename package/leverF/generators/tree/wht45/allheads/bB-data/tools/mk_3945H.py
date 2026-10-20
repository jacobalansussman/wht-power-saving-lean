"""bB-data (wht45/allheads, 19 Oct 2026): the hand-written fill-dependent Lean files of the #3945 lever instance,
re-instantiated for design B (#3945 with every child-8/10 bin a literal head; suffix 3945H / S5H / L3945H).

Each output file is a source file of the #3945 round (p-disp, p-a3, i-final) or of the #2979 round (x-inst), compiled
there, with a list of EXACT text replacements.  Every replacement must occur exactly the stated number of times, and
after all replacements no forbidden leftover (old head count, old namespace, old instance name) may remain; otherwise
the tool REFUSES (exit 3) and writes nothing.  The numbers of design B:
  heads (A4 bins with remainder 48, children 8/10) 80,640 = 40,320 + 40,320; F1 bins with remainder 48 (children
  12/13, not heads) 57,312; remainders {3: 23,328, 5: 4,320, 48: 137,952} (unchanged); unit profile = #3945's;
  cut profile cut3945H = {1: 8041680, 2: 2306880, 3: 729648, 4: 353160, 5: 69120, 8: 25920, 24: 80640, 48: 57312};
  arrays per gauge copy W0 = 322,560, cut Hd = 80,640; DFT figure N = 17581520 (e-r45m1's R5 scan, N + 1 refused):
  seam exponent 1 - 1758152001/10^12.
usage: python3 -I -B mk_3945H.py OUTDIR [--neg-heads]
  --neg-heads: also write the NEGATIVE control files (must FAIL to compile) under OUTDIR/Work/XHostT/Data/S5Hneg/.
"""
import os, sys, hashlib

OUT = sys.argv[1]
L = (os.environ['GEN_ROOT'] + '/wht45/lever')
L2 = (os.environ['GEN_ROOT'] + '/wht45/lever2979')
ENG = os.path.join(L, 'i-final/out/a1/ENG')
IF = os.path.join(L, 'i-final/lean')
XI = os.path.join(L2, 'x-inst/lean')

BY = 'bB-data (wht45/allheads, 19 Oct 2026, design B: #3945 with every child-8/10 bin a literal head)'

jobs = []   # (src, dst, [(old, new, count)], [forbidden leftovers])

# ---------------------------------------------------------------- 1. DispLink (p-disp) on the S5H call tree
jobs.append((os.path.join(L, 'p-disp/lean/Work/XHostT/Data/S5/DispLink.lean'), 'Work/XHostT/Data/S5H/DispLink.lean', [
    ('import Work.XHostT.Data.S5.Calls', 'import Work.XHostT.Data.S5H.Calls', 1),
    ('p-hostck\'s call tree `DS5.ct` (from p-own\'s #3945 host table, `Work/XHostT/Data/S5/Calls.lean`)',
     'p-hostck\'s call tree `DS5H.ct` (design B: p-own\'s converter on the all-heads fill, `Work/XHostT/Data/S5H/Calls.lean`)', 1),
    ('Written by p-disp (wht45/lever, 19 Oct 2026);', 'Written by p-disp (wht45/lever, 19 Oct 2026); instance S5H by ' + BY + ';', 1),
    ('K.DS5\n', 'K.DS5H\n', 2),
    ('K.DS5.', 'K.DS5H.', 2),
], ['DS5.', 'DS5\n', 'Data.S5.']))

# ---------------------------------------------------------------- 2. the schedule (p-a3's Inst/S5.lean)
jobs.append((os.path.join(ENG, 'Work/XHostT/Inst/S5.lean'), 'Work/XHostT/Inst/S5H.lean', [
    ('import Work.XHostT.Data.S5.Whole', 'import Work.XHostT.Data.S5H.Whole', 1),
    ('# (key: bridge-host-cross) The host table of hub PR #3945 as a Lean cross-replica host schedule',
     '# (key: bridge-host-cross) The host table of design B (hub PR #3945, all-heads fill) as a Lean cross-replica host schedule', 'H'),
    ('p-hostck\'s whole-table theorem `DS5.whole`', 'p-hostck\'s whole-table theorem `DS5H.whole` (design B: the all-heads fill of round wht45/reopt)', 'H'),
    ('Written by p-a3 (wht45/lever, 19 Oct 2026);', 'Written by p-a3 (wht45/lever, 19 Oct 2026); instance S5H by ' + BY + ';', 'H'),
    ('the cross-replica host schedule of #3945.', 'the cross-replica host schedule of design B (#3945, all-heads fill).', 'H'),
    ('K.DS5\n', 'K.DS5H\n', 1),
    ('sched3945', 'sched3945H', 6),
], ['DS5.', 'DS5\n', 'Data.S5.', 'sched3945 ', 'sched3945_']))

# ---------------------------------------------------------------- 3. head data (i-final's Data/S5/Heads.lean)
HEADS_SRC = os.path.join(IF, 'Work/XHostT/Data/S5/Heads.lean')
jobs.append((HEADS_SRC, 'Work/XHostT/Data/S5H/Heads.lean', [
    ('import Work.XHostT.Data.S5.Calls', 'import Work.XHostT.Data.S5H.Calls', 1),
    ('/-! i-final (wht45/lever): per-kind data facts of #3945\'s host table S5 (p-hostck\'s data)',
     '/-! ' + BY + ': per-kind data facts of the host table S5H (p-hostck\'s generator on the all-heads fill)', 1),
    ('the number of A4 bins (kind 4) with remainder 48 (the 57,312 literal\nheads: 28,672 of child 8 + 28,640 of child 10), the F1 bins with remainder 48 (80,640, not heads)',
     'the number of A4 bins (kind 4) with remainder 48 (the 80,640 literal\nheads: 40,320 of child 8 + 40,320 of child 10), the F1 bins with remainder 48 (57,312, not heads)', 1),
    ('Same text as x-tables\'\nData/H2979/Heads.lean with the instance changed.', 'i-final\'s Data/S5/Heads.lean\nwith the instance and the counts changed.', 1),
    ('K.DS5\n', 'K.DS5H\n', 2),
    ('K.DS5.', 'K.DS5H.', 3),
    ('/-- the A4 bins with remainder 48 (all of them unhosted: the 57,312 literal heads). -/\ntheorem heads_cnt : cntK rt 48 4 360 = 57312 := by decide +kernel',
     '/-- the A4 bins with remainder 48 (all of them unhosted: the 80,640 literal heads). -/\ntheorem heads_cnt : cntK rt 48 4 360 = 80640 := by decide +kernel', 1),
    ('theorem f1_cnt : cntK rt 48 6 360 = 80640 := by decide +kernel', 'theorem f1_cnt : cntK rt 48 6 360 = 57312 := by decide +kernel', 1),
], ['DS5.', 'DS5\n', 'Data.S5.', '28,672', '= 57312 := by decide +kernel\ntheorem f1', '6 360 = 80640']))

# ---------------------------------------------------------------- 4. schedule facts (i-final's Inst/S5Facts.lean)
jobs.append((os.path.join(IF, 'Work/XHostT/Inst/S5Facts.lean'), 'Work/XHostT/Inst/S5HFacts.lean', [
    ('import Work.XHostT.Inst.S5\n', 'import Work.XHostT.Inst.S5H\n', 1),
    ('import Work.XHostT.Data.S5.Glob', 'import Work.XHostT.Data.S5H.Glob', 1),
    ('import Work.XHostT.Data.S5.Heads', 'import Work.XHostT.Data.S5H.Heads', 1),
    ('import Work.XHostT.Data.S5.DispLink', 'import Work.XHostT.Data.S5H.DispLink', 1),
    ('# (key: bridge-host-cross) The data facts of #3945\'s schedule that the bridge needs',
     '# (key: bridge-host-cross) The data facts of design B\'s schedule (#3945, all-heads fill) that the bridge needs', 'H'),
    ('For p-a3\'s `sched3945` (Inst/S5.lean, the cross-replica host schedule of hub PR #3945 from p-hostck\'s checked table):',
     'For `sched3945H` (Inst/S5H.lean, the cross-replica host schedule of design B from p-hostck\'s checked table S5H):', 'H'),
    ('57,312 hosts of kind A4 have remainder 48', '80,640 hosts of kind A4 have remainder 48', 1),
    ('The proofs are x-tables\' (wht45/lever2979, Inst/S2979Facts.lean) with the instance changed.\ni-final (wht45/lever, 19 Oct 2026);',
     'The proofs are i-final\'s Inst/S5Facts.lean (wht45/lever) with the instance and the head counts changed.\n' + BY + ';', 'H'),
    ('K.DS5\n', 'K.DS5H\n', 1),
    ('/-- **the head count of #3945**: 57,312 hosts of kind A4 (both sides) have remainder 48. -/',
     '/-- **the head count of design B**: 80,640 hosts of kind A4 (both sides) have remainder 48. -/', 1),
    ('sched3945.hrem p.1 4 p.2 = 48} = 57312 := by', 'sched3945H.hrem p.1 4 p.2 = 48} = 80640 := by', 1),
    ('remAt rt (2016 * p.1.val + 224 * (4 : Fin 9).val + (if p.2.1 then 112 else 0) + p.2.2.val) = 48} = 57312',
     'remAt rt (2016 * p.1.val + 224 * (4 : Fin 9).val + (if p.2.1 then 112 else 0) + p.2.2.val) = 48} = 80640', 1),
    ('sched3945.hrem p.1 6 p.2 = 48} = 80640 := by', 'sched3945H.hrem p.1 6 p.2 = 48} = 57312 := by', 1),
    ('remAt rt (2016 * p.1.val + 224 * (6 : Fin 9).val + (if p.2.1 then 112 else 0) + p.2.2.val) = 48} = 80640',
     'remAt rt (2016 * p.1.val + 224 * (6 : Fin 9).val + (if p.2.1 then 112 else 0) + p.2.2.val) = 48} = 57312', 1),
    ('/-- **the price of the remaining climbs of #3945\'s host chains.** -/', '/-- **the price of the remaining climbs of design B\'s host chains** (the same idle profile as #3945). -/', 1),
    ('sched3945_hcost', 'sched3945H_hcost', 3),
    ('sched3945_heads', 'sched3945H_heads', 3),
    ('sched3945_f1', 'sched3945H_f1', 2),
    ('sched3945_hrem', 'sched3945H_hrem', 2),
    ('hcostT sched3945 ', 'hcostT sched3945H ', 2),
    ('resAt3945_ge', 'resAt3945H_ge', 2),
    ('resK3945_disp', 'resK3945H_disp', 4),
    ('resK3945_ge5', 'resK3945H_ge5', 3),
    ('resK3945_pos', 'resK3945H_pos', 3),
    ('resK3945_eq0', 'resK3945H_eq0', 3),
    ('card_live3945', 'card_live3945H', 3),
    ('/-- **the live roles of #3945\'s network**', '/-- **the live roles of design B\'s network** (the same as #3945\'s)', 1),
], ['DS5.', 'DS5\n', 'Data.S5.', 'Inst.S5\n', 'sched3945 ', 'sched3945.', 'sched3945_', '57,312 hosts', '3945_']))

# ---------------------------------------------------------------- 5. head bins (i-final's Inst/S5HB.lean)
jobs.append((os.path.join(IF, 'Work/XHostT/Inst/S5HB.lean'), 'Work/XHostT/Inst/S5HHB.lean', [
    ('import Work.XHostT.Inst.S5Facts', 'import Work.XHostT.Inst.S5HFacts', 1),
    ('# (key: heads-first) #3945\'s head bins: `headBins sched3945` has 57,312 elements',
     '# (key: heads-first) design B\'s head bins: `headBins sched3945H` has 80,640 elements', 1),
    ('p-hfirst\'s `headBins` (A4 bins with remaining climb 48 in which no slot is hosted) on `sched3945`: with positive\nresiduals (`resK3945_pos`) p-hfirst\'s `headBins_eq` reduces it to the A4 bins with remainder 48, counted by\n`sched3945_heads` (x-tables\' S2979HB.lean with the instance changed).\ni-final (wht45/lever, 19 Oct 2026);',
     'p-hfirst\'s `headBins` (A4 bins with remaining climb 48 in which no slot is hosted) on `sched3945H`: with positive\nresiduals (`resK3945H_pos`) p-hfirst\'s `headBins_eq` reduces it to the A4 bins with remainder 48, counted by\n`sched3945H_heads` (i-final\'s S5HB.lean with the instance and the count changed).\n' + BY + ';', 1),
    ('K.DS5\n', 'K.DS5H\n', 1),
    ('/-- **the head bins of #3945**: 57,312. -/\ntheorem headBins3945_card : (headBins sched3945).card = 57312 := by',
     '/-- **the head bins of design B**: 80,640. -/\ntheorem headBins3945H_card : (headBins sched3945H).card = 80640 := by', 1),
    ('rw [headBins_eq sched3945 (resK ct resT 360 408) sched3945_fits resK3945_pos]',
     'rw [headBins_eq sched3945H (resK ct resT 360 408) sched3945H_fits resK3945H_pos]', 1),
    ('exact sched3945_heads', 'exact sched3945H_heads', 1),
    ('Inst.headBins3945_card', 'Inst.headBins3945H_card', 1),
], ['DS5.', 'DS5\n', 'sched3945 ', 'sched3945_', '57312', '3945_card']))

# ---------------------------------------------------------------- 6. the bridge instance (x-inst's Inst2979.lean, #3945H numbers)
def bridge_num(inst):
    return [
        ('2979', '3945H', None),
        ('16086', '80640', None),
        ('66304', '322560', None),
        ('1755647301', '1758152001', None),
        ('1.7556473e-3', '1.7581520e-3', None),
    ] + ([('16087', '80641', None), ('N = 17556473', 'N = 17581520', None)] if inst else [])

INST_HDR_OLD = """# x-inst: the #2979 instance of p-bridge's bridge (round wht45/lever2979)

p-bridge's `Inst.lean` (#3945) with the numbers of hub PR #2979 (mdumitrean, on Jeremy Kothe's #2957 and Mason Hale's
#2917 words; T74, 74 replicas, 666 calls).  Target rule: `card α = 73`, heads of rank `48 = 24 + 24` cut by 24 lines
(target 49), rank-24 remnants, every paid rank at most 48 (`48 < 49`), `W0 = 66304` arrays per gauge copy
(74 replicas x 896) of which `Hd = 16086` are cut, the cut profile `cut2979` of R5 (`Work.Mixed.Rate2979`).

* `unit2979`: the expected UNIT profile of the #2979 network (p-alt's measured full profile = the hub's
  `allocation/PROFILE.json`; the 16,086 literal heads still rank 48), and `hcut2979`: unit profile -> cut profile
  (checked by `decide`);"""
INST_HDR_NEW = """# bB-data: the design-B instance (#3945H) of p-bridge's bridge (round wht45/allheads)

x-inst's `Inst2979.lean` (itself p-bridge's `Inst.lean` with other numbers) with the numbers of design B: hub PR #3945
(Dugongue's staged544 T360 terminal-gap fill, 360 replicas, 3,240 calls) with the all-heads fill of round wht45/reopt.
Target rule: `card α = 73`, heads of rank `48 = 24 + 24` cut by 24 lines (target 49), rank-24 remnants, every paid
rank at most 48 (`48 < 49`), `W0 = 322560` arrays per gauge copy (360 replicas x 896) of which `Hd = 80640` are cut,
the cut profile `cut3945H` of R5 (`Work.Mixed.Rate3945H`).

* `unit3945H`: the expected UNIT profile of the design-B network (= #3945's PROFILE.json, unchanged by the new fill;
  the 80,640 literal heads still rank 48; `unit3945H_eq : unit3945H = unit3945`), and `hcut3945H`: unit profile ->
  cut profile (checked by `decide`);"""
jobs.append((os.path.join(XI, 'Work/Mixed/Bridge/Inst2979.lean'), 'Work/Mixed/Bridge/Inst3945H.lean', [
    (INST_HDR_OLD, INST_HDR_NEW, 'H'),
    ('Written by x-inst (wht45/lever2979, 19 Oct 2026) on p-bridge\'s Inst.lean;', 'Written by x-inst (wht45/lever2979) on p-bridge\'s Inst.lean; instance 3945H by ' + BY + ';', 'H'),
    ('import Work.Mixed.Bridge.Core\nimport Work.Mixed.Rate2979', 'import Work.Mixed.Bridge.Core\nimport Work.Mixed.Bridge.Inst\nimport Work.Mixed.Rate3945H', 1),
    ('/-- Expected unit profile of the #2979 network (rank, count): the 16,086 literal heads are rank 48. -/\ndef unit2979 : List (ℕ × ℕ) :=\n  [(1, 1645680), (2, 470196), (3, 159468), (4, 71928), (5, 13716), (8, 5328), (48, 28190)]',
     '/-- Expected unit profile of the design-B network (rank, count): #3945\'s profile; the 80,640 literal heads are rank 48. -/\ndef unit3945H : List (ℕ × ℕ) :=\n  [(1, 8041680), (2, 2306880), (3, 729648), (4, 353160), (5, 69120), (8, 25920), (48, 137952)]\n\n/-- the unit profile of design B is #3945\'s (only the fill changed, not the profile). -/\ntheorem unit3945H_eq : unit3945H = unit3945 := rfl', 1),
    ('/-- Unit profile -> R5\'s cut profile for #2979: 16,086 heads of rank 48 become rank-24 remnants. -/',
     '/-- Unit profile -> R5\'s cut profile for design B: 80,640 heads of rank 48 become rank-24 remnants. -/', 1),
    ('/-- **(f) for #2979**:', '/-- **(f) for design B (#3945H)**:', 1),
    ('/-- **THE ENGINE THEOREM FOR #2979 (half cut)**', '/-- **THE ENGINE THEOREM FOR DESIGN B, #3945H (half cut)**', 1),
] + [(a, b, c) for a, b, c in bridge_num(inst=True)], ['2979', '16086', '66304', '1755647', '57312', '1750110']))

HF_HDR_OLD = """# x-inst: the #2979 bridge in the exact shape of p-hfirst's `XNet.heads_first` conclusion (round wht45/lever2979)

p-bridge's `HeadsFirst.lean` with the #2979 instance"""
HF_HDR_NEW = """# bB-data: the design-B (#3945H) bridge in the exact shape of p-hfirst's `XNet.heads_first` conclusion (round wht45/allheads)

x-inst's `HeadsFirst2979.lean` (p-bridge's `HeadsFirst.lean` with other numbers) with the design-B instance"""
jobs.append((os.path.join(XI, 'Work/Mixed/Bridge/HeadsFirst2979.lean'), 'Work/Mixed/Bridge/HeadsFirst3945H.lean', [
    (HF_HDR_OLD, HF_HDR_NEW, 'H'),
    ('Written by x-inst (wht45/lever2979, 19 Oct 2026) on p-bridge\'s file;', 'Written by x-inst (wht45/lever2979) on p-bridge\'s file; instance 3945H by ' + BY + ';', 'H'),
    ('/-- **THE ENGINE THEOREM FOR #2979 FROM p-hfirst\'s CONCLUSION.**', '/-- **THE ENGINE THEOREM FOR DESIGN B (#3945H) FROM p-hfirst\'s CONCLUSION.**', 1),
    ('/-- **The engine theorem for #2979 with the cut data given directly**', '/-- **The engine theorem for design B (#3945H) with the cut data given directly**', 1),
    ('/-- (f) for #2979 from', '/-- (f) for design B (#3945H) from', 1),
] + [(a, b, c) for a, b, c in bridge_num(inst=False)], ['2979', '16086', '66304', '1755647', '57312', '1750110']))

# ---------------------------------------------------------------- 7. the warp-network engine theorem (i-final's HFTest3945W.lean)
HFT_HDR_OLD = """# i-final: the #3945 engine theorem on the WARP network (the network of the real #3945 provider words)

The #3945 provider words are windowed packages (`KLW.GInvEW .. okW`, p-wgen), so the real network is p-a3's
`xnetGeomW`.  Its heads-first certificate (`XHostTW.xnetGeomW_heads_first`: x-inst's generic warp port of p-hfirst's
heads-first stage, Work/XHostT/HFW2979/, written for wht45/lever2979 and generic in the instance) is fed into
p-bridge's `engine3945_of_cut`, exactly as p-bridge's `HFTest.engine3945_geom` does for the plain network
(x-inst's HFTest2979W.lean with the #3945 numbers):
* `engine3945_geomW`: any head set `hb` (never hosted, remaining climb 48), data facts `hU` (unit price
  `unitPriceTW`), `hcut`, `hLlen`, `hσ`, `hα`;
* `engine3945_geomW_auto`: 360 replicas, classes `Fin 112`, the canonical head set `headBins S`, the head count in
  Finset form, positive residuals, any unit profile `U` with `hcut`;
* `engine3945_geomW_unit`: the same with `U = unit3945` (`hcut3945`).
i-final (wht45/lever, 19 Oct 2026); prepared with assistance from Claude.  No `sorry`."""
HFT_HDR_NEW = """# bB-data: the design-B (#3945H) engine theorem on the WARP network (the network of the real #3945 provider words)

i-final's `HFTest3945W.lean` with the numbers of design B (#3945's words, order and dispatch with the all-heads fill:
80,640 literal heads, figure 1.7581520e-3).  The provider words are windowed packages (`KLW.GInvEW .. okW`, p-wgen),
so the real network is p-a3's `xnetGeomW`; x-inst's generic warp heads-first certificate
(`XHostTW.xnetGeomW_heads_first`) is fed into `engine3945H_of_cut` (Work/Mixed/Bridge/HeadsFirst3945H.lean):
* `engine3945H_geomW`: any head set `hb` (never hosted, remaining climb 48), data facts `hU` (unit price
  `unitPriceTW`), `hcut`, `hLlen`, `hσ`, `hα`;
* `engine3945H_geomW_auto`: 360 replicas, classes `Fin 112`, the canonical head set `headBins S`, the head count in
  Finset form (80,640), positive residuals, any unit profile `U` with `hcut`;
* `engine3945H_geomW_unit`: the same with `U = unit3945H` (`hcut3945H`).
""" + BY + """; prepared with assistance from Claude.  No `sorry`."""
jobs.append((os.path.join(IF, 'Work/Mixed/Bridge/HFTest3945W.lean'), 'Work/Mixed/Bridge/HFTest3945HW.lean', [
    (HFT_HDR_OLD, HFT_HDR_NEW, 'H'),
    ('import Work.Mixed.Bridge.HFTest\n', 'import Work.Mixed.Bridge.HFTest\nimport Work.Mixed.Bridge.HeadsFirst3945H\n', 1),
    ('namespace HFTest3945W', 'namespace HFTest3945HW', 1),
    ('end HFTest3945W', 'end HFTest3945HW', 1),
    ('/-- **The #3945 engine theorem from the geometric WARP network.** -/', '/-- **The design-B (#3945H) engine theorem from the geometric WARP network.** -/', 1),
    ('/-- **The #3945 engine theorem from the warp network and the table data**', '/-- **The design-B (#3945H) engine theorem from the warp network and the table data**', 1),
    ('/-- The same with the expected unit profile `unit3945` (`hcut3945`). -/', '/-- The same with the expected unit profile `unit3945H` (`hcut3945H`). -/', 1),
    ('engine3945_geomW', 'engine3945H_geomW', 8),
    ('engine3945_of_cut', 'engine3945H_of_cut', 1),
    ('HFTest3945W.', 'HFTest3945HW.', 3),
    ('57312', '80640', 7),
    ('MX.cut3945 ', 'MX.cut3945H ', 2),
    ('1750110901', '1758152001', 3),
    ('unit3945 hU hcut3945 m', 'unit3945H hU hcut3945H m', 1),
    ('listCost φ unit3945)', 'listCost φ unit3945H)', 1),
], ['57312', '1750110', 'cut3945 ', 'engine3945_', 'HFTest3945W', 'unit3945)', 'unit3945 ']))

# ---------------------------------------------------------------- 8. integration drafts: L3945H adapter, providers, engine (i-final's)
jobs.append((os.path.join(IF, 'Work/L3945/EngineOfGeomW.lean'), 'Work/L3945H/EngineOfGeomW.lean', [
    ('import Work.L3945.Stmt', 'import Work.L3945H.Stmt', 1),
    ('import Work.Mixed.Bridge.HFTest3945W', 'import Work.Mixed.Bridge.HFTest3945HW', 1),
    ('# i-final: the adapter from the #3945 warp network to the seam\'s engine statement `MX.L3945.EngineStmt`',
     '# bB-data: the adapter from the design-B (#3945H) warp network to the seam\'s engine statement `MX.L3945H.EngineStmt`', 'H'),
    ('`MX.L3945.engine_of_geomW`', '`MX.L3945H.engine_of_geomW`', 'H'),
    ('`1 - 1750110901/10^12` (DFT figure 1.7501109e-3), i.e. `MX.L3945.EngineStmt.{0}`',
     '`1 - 1758152001/10^12` (DFT figure 1.7581520e-3), i.e. `MX.L3945H.EngineStmt.{0}`', 'H'),
    ('`Work/L3945/Stmt.lean`', '`Work/L3945H/Stmt.lean` (suffix L3945H)', 'H'),
    ('(x-inst\'s L2979/EngineOfGeomW.lean with the #3945\nnumbers.)  i-final (wht45/lever, 19 Oct 2026);',
     '(i-final\'s L3945/EngineOfGeomW.lean with the design-B\nnumbers.)  ' + BY + ';', 'H'),
    ('namespace L3945\n', 'namespace L3945H\n', 1),
    ('end L3945\n', 'end L3945H\n', 1),
    ('/-- **The engine statement of L3945 from the #3945 warp network and its data.** -/',
     '/-- **The engine statement of L3945H from the design-B warp network and its data.** -/', 1),
    ('S.hrem p.1 4 p.2 = 48)).card = 57312)', 'S.hrem p.1 4 p.2 = 48)).card = 80640)', 1),
    ('listCost φ PB.unit3945)', 'listCost φ PB.unit3945H)', 1),
    ('PB.HFTest3945W.engine3945_geomW_unit', 'PB.HFTest3945HW.engine3945H_geomW_unit', 1),
    ('MX.L3945.engine_of_geomW\n', 'MX.L3945H.engine_of_geomW\n', 1),
], ['57312', '1750110', 'L3945.', 'L3945\n', 'unit3945)', 'HFTest3945W']))

PROV_HDR_OLD = """# i-final: the data of the #3945 warp network (providers, dispatch, window facts, residuals, unit price)"""
PROV_HDR_NEW = """# bB-data: the data of the design-B (#3945H) warp network (providers, dispatch, window facts, residuals, unit price)

i-final's `Work/L3945/Prov.lean` with the schedule and host table of design B (`Inst.sched3945H`, `K.DS5H`): the
providers, dispatch and window facts do not depend on the fill (same words, order and dispatch); the residual check,
`hfit`, `hpos` and the unit price are restated on the S5H table.  Original description:"""
jobs.append((os.path.join(IF, 'Work/L3945/Prov.lean'), 'Work/L3945H/Prov.lean', [
    (PROV_HDR_OLD, PROV_HDR_NEW, 'H'),
    ('p-hostck\'s residual table `DS5.resT`', 'p-hostck\'s residual table `DS5H.resT`', 'H'),
    ('`listCost φ unit3945`', '`listCost φ unit3945H`', 'H'),
    ('import Work.XHostT.Inst.S5Facts', 'import Work.XHostT.Inst.S5HFacts', 1),
    ('import Work.Mixed.Bridge.Inst\n', 'import Work.Mixed.Bridge.Inst3945H\n', 1),
    ('i-final (wht45/lever, 19 Oct 2026); prepared', 'i-final (wht45/lever, 19 Oct 2026); instance 3945H by ' + BY + '; prepared', 'H'),
    ('namespace L3945\n', 'namespace L3945H\n', 1),
    ('end L3945\n', 'end L3945H\n', 1),
    ('K.DS5.', 'K.DS5H.', 8),
    ('Inst.sched3945', 'Inst.sched3945H', 5),
    ('`sched3945_fits`', '`sched3945H_fits`', 2),
    ('Inst.resK3945_disp', 'Inst.resK3945H_disp', 1),
    ('Inst.resK3945_pos', 'Inst.resK3945H_pos', 1),
    ('PB.unit3945', 'PB.unit3945H', 3),
    ('MX.L3945.', 'MX.L3945H.', 6),
    ('/-- **the unit price of the #3945 warp network is its profile** `unit3945`. -/',
     '/-- **the unit price of the design-B warp network is its profile** `unit3945H` (= #3945\'s). -/', 'H'),
], ['DS5.', 'sched3945 ', 'sched3945.', 'sched3945_', 'unit3945 ', 'unit3945\n', 'L3945.', 'L3945\n', 'S5Facts', 'Bridge.Inst\n']))

jobs.append((os.path.join(IF, 'Work/L3945/Engine.lean'), 'Work/L3945H/Engine.lean', [
    ('import Work.L3945.EngineOfGeomW\nimport Work.L3945.Prov', 'import Work.L3945H.EngineOfGeomW\nimport Work.L3945H.Prov', 1),
    ('# i-final: THE engine theorem of the #3945 lever instance L3945', '# bB-data: THE engine theorem of the design-B lever instance L3945H (#3945 with the all-heads fill)', 1),
    ('`MX.L3945.engine : MX.L3945.EngineStmt.{0}`: OpenAI\'s `hills_program` statement at the seam exponent\n`1 - 1750110901/10^12` (DFT figure 1.7501109e-3), from the #3945 warp network',
     '`MX.L3945H.engine : MX.L3945H.EngineStmt.{0}`: OpenAI\'s `hills_program` statement at the seam exponent\n`1 - 1758152001/10^12` (DFT figure 1.7581520e-3), from the #3945 warp network with the all-heads fill of round\nwht45/reopt (design B)', 1),
    ('with the lever (half cut of the 57,312 literal rank-48 heads', 'with the lever (half cut of the 80,640 literal rank-48 heads', 1),
    ('This is the one theorem i-skel\'s generated chain (`Work.FourierL3945.Seam`) imports.\ni-final (wht45/lever, 19 Oct 2026);',
     'This is the one theorem the generated chain (`Work.FourierL3945H.Seam`) imports.\ni-final\'s Engine.lean (wht45/lever) with the design-B instance; ' + BY + ';', 1),
    ('namespace L3945\n', 'namespace L3945H\n', 1),
    ('end L3945\n', 'end L3945H\n', 1),
    ('/-- **THE ENGINE THEOREM OF L3945.** -/', '/-- **THE ENGINE THEOREM OF L3945H.** -/', 1),
    ('XHostT.Inst.sched3945 hfit3945W hpos3945W\n    XHostT.Inst.sched3945_heads hU3945W',
     'XHostT.Inst.sched3945H hfit3945W hpos3945W\n    XHostT.Inst.sched3945H_heads hU3945W', 1),
    ('MX.L3945.engine\n', 'MX.L3945H.engine\n', 1),
], ['57,312', '1750110', 'L3945.', 'L3945\n', 'sched3945 ', 'sched3945_']))

# ---------------------------------------------------------------- negatives (must FAIL to compile)
NEG = []
if '--neg-heads' in sys.argv:
    NEG.append((os.path.join(IF, 'Work/XHostT/Data/S5neg/NegHeads.lean'), 'Work/XHostT/Data/S5Hneg/NegHeads.lean', [
        ('import Work.XHostT.Data.S5.Calls', 'import Work.XHostT.Data.S5H.Calls', 1),
        ('/-! i-final NEGATIVE control (must FAIL to compile): one head more than #3945\'s table has. -/',
         '/-! ' + BY + ': NEGATIVE control (must FAIL to compile): one head more than the S5H table has. -/', 1),
        ('K.DS5\n', 'K.DS5H\n', 2),
        ('theorem heads_cnt_plus1 : cntK rt 48 4 360 = 57313 := by decide +kernel', 'theorem heads_cnt_plus1 : cntK rt 48 4 360 = 80641 := by decide +kernel', 1),
    ], ['DS5\n', 'Data.S5.', '57313']))
    NEG.append((os.path.join(IF, 'Work/XHostT/Data/S5neg/NegRes.lean'), 'Work/XHostT/Data/S5Hneg/NegRes.lean', [
        ('import Work.XHostT.Data.S5.Calls', 'import Work.XHostT.Data.S5H.Calls', 1),
        ('/-! i-final NEGATIVE control (must FAIL to compile): a residual bound one too high (#3945\'s smallest width is 5). -/',
         '/-! ' + BY + ': NEGATIVE control (must FAIL to compile): a residual bound one too high (the smallest width is 5). -/', 1),
        ('K.DS5\n', 'K.DS5H\n', 2),
    ], ['DS5\n', 'Data.S5.']))
    # the head count of the OLD fill must fail on the new table (the data really changed)
    NEG.append((os.path.join(IF, 'Work/XHostT/Data/S5neg/NegHeads.lean'), 'Work/XHostT/Data/S5Hneg/NegHeadsOld.lean', [
        ('import Work.XHostT.Data.S5.Calls', 'import Work.XHostT.Data.S5H.Calls', 1),
        ('/-! i-final NEGATIVE control (must FAIL to compile): one head more than #3945\'s table has. -/',
         '/-! ' + BY + ': NEGATIVE control (must FAIL to compile): the head count of #3945\'s OWN fill (57,312) on the S5H table. -/', 1),
        ('K.DS5\n', 'K.DS5H\n', 2),
        ('theorem heads_cnt_plus1 : cntK rt 48 4 360 = 57313 := by decide +kernel', 'theorem heads_cnt_old : cntK rt 48 4 360 = 57312 := by decide +kernel', 1),
    ], ['DS5\n', 'Data.S5.', '57313']))

outs = []
fail = []
for src, dst, reps, forb in jobs + NEG:
    t = open(src, encoding='utf-8').read()
    held = []
    # header / prose replacements (count 'H'): taken out first (placeholders), so that the instance replacements and
    # the leftover check below see only the Lean text and the other prose; put back at the end.
    for old, new, cnt in reps:
        if cnt != 'H':
            continue
        n = t.count(old)
        if n != 1:
            fail.append('%s: header %r occurs %d times, expected 1' % (dst, old[:70], n))
            continue
        ph = '\x00HDR%d\x00' % len(held)
        t = t.replace(old, ph); held.append((ph, new))
    for old, new, cnt in reps:
        if cnt == 'H':
            continue
        n = t.count(old)
        if cnt is not None and n != cnt:
            fail.append('%s: %r occurs %d times, expected %d' % (dst, old[:70], n, cnt))
            continue
        if cnt is None and n == 0:
            fail.append('%s: %r does not occur' % (dst, old[:70]))
        t = t.replace(old, new)
    for f in forb:
        if f in t:
            i = t.index(f)
            fail.append('%s: forbidden leftover %r near %r' % (dst, f, t[max(0, i - 60):i + 40]))
    for ph, new in held:
        t = t.replace(ph, new)
    if 'sorry' in t.replace('No `sorry`', '').replace('no `sorry`', ''):
        fail.append('%s: contains sorry' % dst)
    outs.append((src, dst, t))
if fail:
    sys.stderr.write('REFUSED (nothing written):\n  ' + '\n  '.join(fail) + '\n')
    sys.exit(3)
for src, dst, t in outs:
    p = os.path.join(OUT, dst)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, 'w', encoding='utf-8').write(t)
    print('%s  %s  <- %s (%s)' % (hashlib.sha256(t.encode()).hexdigest()[:16], dst, os.path.relpath(src, (os.environ['GEN_ROOT'] + '/wht45')),
                                  hashlib.sha256(open(src, 'rb').read()).hexdigest()[:16]))
