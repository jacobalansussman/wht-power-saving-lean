"""fs-price (from h-L3-build's gen_h588m.py): the instance U8XMHH = k-combine's U8XM certificate (GXD.U8XM: the U8 word of
Kothe #869/#870, Boukhalfa #890, Capalbo #942; hub exit frames; x-mutable's mutable sources; 5 centre-source copies kept)
through h-L3-build's HOSTED network with a same-run host table (host_price.py --samerun 1 on U8XM's 509 residual widths,
unhosted residuals re-packed by s-stack's host_repack2.py).
Writes the data modules (host rows, bank rows, pads, remainders, checked by XHost.hostSchedK), the rate certificate
(Shea's/Sussman's foldgen at s = 40), the unit price, and a JSON summary.  Every Lean check is mirrored here first.
usage: python3 -I gen_hu8xm.py OUTDIR HOSTTAB.json REPACKED2.json"""
import sys, os, re, json, gzip
from collections import Counter, defaultdict
BC = os.environ.get('BC', '.')  # as run: the agents' working folder; see generators/README.txt
sys.path.insert(0, BC + '/m-gen/toolkit/tools/gen')
import foldlib as fl
import foldgen as fg
OUT = sys.argv[1]
NAME = 'U8XMHH'
CERT = 'U8XM'
s = open(BC + '/k-combine/out/global-U8XM/Work/GCert/Data/Gen/U8XMSched.lean').read()
res = [int(x) for x in re.search(r'def schedRes : List Nat := \[([^\]]*)\]', s).group(1).split(',')]
R = len(res); assert R == 509
HT, RP = sys.argv[2], sys.argv[3]
ht = json.load(open(HT))
rp = json.load(open(RP))
price = json.load(open(BC + '/k-combine/out/price/U8XM.json'))
nT = 112; NH = 18 * nT
assert ht['nT'] == nT and len(ht['rem']) == NH and Counter(ht['res']) == Counter(res)
# emitter slot e (width order) -> Lean slot q (same width, same rank among that width)
qs = defaultdict(list)
for q, w in enumerate(res):
    qs[w].append(q)
seen = Counter(); perm = []
for e, w in enumerate(ht['res']):
    perm.append(qs[w][seen[w]]); seen[w] += 1
assert sorted(perm) == list(range(R))
# host rows
hr = [[2016 + q for q in range(R)] for _ in range(9)]
for l in range(9):
    for e, v in enumerate(ht['rows'][l]):
        if v:
            hr[l][perm[e]] = v - 1
rem = list(ht['rem'])
# bank rows from the re-packed bank items (width, call k)
banks = rp['bank_items']; nb = len(banks); assert nb == rp['banks']
free = {}
for k in range(72):
    l = k % 9
    d = defaultdict(list)
    for q in range(R):
        if hr[l][q] >= 2016:
            d[res[q]].append(q)
    free[k] = d
br = [[nb + q for q in range(R)] for _ in range(72)]
used = defaultdict(Counter)
for b, items in enumerate(banks):
    for w, k in items:
        q = free[k][w][used[k][w]]; used[k][w] += 1
        assert br[k][q] == nb + q
        br[k][q] = b
for k in range(72):
    for w, lst in free[k].items():
        assert used[k][w] == len(lst), ('unplaced unhosted slots', k, w)
pad = [73 - sum(w for w, _ in items) for items in banks]
assert Counter(pad) == Counter({int(a): b for a, b in rp['pads'].items()})
# ---- mirror of XHost.hostSchedK
def inWin(c, l):
    return [l < 2, l == 2, l < 4, l == 4, l < 6, l == 6, 3 <= l < 9, 5 <= l < 9, 7 <= l < 9][c] if c < 9 else False
rank = [16, 8, 32, 8, 48, 8, 48, 32, 16]
assert 0 < nb and len(hr) == 9 and all(len(r) == R and len(set(r)) == R and max(r) < 2016 + R for r in hr)
assert len(br) == 72 and all(len(r) == R and len(set(r)) == R and max(r) < nb + R for r in br)
for k in range(72):
    for q in range(R):
        assert (hr[k % 9][q] < 2016) == (not br[k][q] < nb)
for l in range(9):
    for e in hr[l]:
        assert e >= 2016 or inWin(e // 224, l)
assert all(p < 73 for p in pad) and len(rem) == NH and all(rem[i] <= rank[i // 224] for i in range(NH))
assert 72 * sum(res) < 2 ** 20
bl = Counter(); hl = Counter()
for k in range(72):
    for q in range(R):
        if br[k][q] < nb: bl[br[k][q]] += res[q]
for l in range(9):
    for q in range(R):
        if hr[l][q] < 2016: hl[hr[l][q]] += res[q]
assert all(bl[b] == 73 - pad[b] for b in range(nb))
assert all(hl[i] == rank[i // 224] - rem[i] for i in range(NH))
# ---- the unit: 72 invocations, 8 runs of remaining climbs, one pad per bank
lp = Counter({int(a): b for a, b in price['local_profile'].items()})
hist = Counter({r: 72 * n for r, n in lp.items()})
remc = Counter(x for x in rem if x)
for r, n in remc.items(): hist[r] += 8 * n
for p in pad:
    if p: hist[p] += 1
W = 64 * nT + nb
assert dict(hist) == {int(a): b for a, b in rp['hist'].items()}, (dict(hist), rp['hist'])
assert W == rp['W'] and 73 * W - sum(r * n for r, n in hist.items()) == 2560
u = dict(m=73, W=W, D=2560, H=dict(sorted(hist.items())))
A = rp['kernel_A']
assert fl.plan_block(u, A, 10, 40) is not None and fl.plan_block(u, A + 1, 10, 40) is None
print(NAME + ': kernel numerator', A, 'DFT', A - 1, 'W', W, 'banks', nb, 'hosted entries/run', sum(1 for l in range(9) for e in hr[l] if e < 2016))
# ---- Lean texts
D = os.path.join(OUT, 'Work/XHost/Data'); os.makedirs(D, exist_ok=True)
hdr = '/-! fs-price gen_hu8xm.py (after h-L3-build gen_h588m.py): data of `GXD.%s` (U8 + exit frames + mutable sources, 5 centre copies, + same-run hosting). -/\n' % NAME
parts = []
rtop = max(rem); ptop = max(pad)
rc = Counter(rem); pc = Counter(pad)
t = 'import Work.XHost.HSchedK\nimport Work.GCert.Data.Gen.U8XMSched\nimport Work.GCert.Data.Count\n' + hdr
t += 'set_option maxRecDepth 100000\nset_option maxHeartbeats 4000000\nnamespace GXD.%s\nopen GXD OAI.PowerSaving OAI.PowerSaving.XHost\n' % NAME
t += 'noncomputable def hrowsL : List (List Nat) := %s\n' % hr
for i in range(9):
    t += 'noncomputable def brows%d : List (List Nat) := %s\n' % (i, br[8 * i:8 * i + 8])
t += 'noncomputable def browsL : List (List Nat) := %s\n' % ' ++ '.join('brows%d' % i for i in range(9))
t += 'noncomputable def hpadL : List Nat := %s\nnoncomputable def hremL : List Nat := %s\n' % (pad, rem)
t += 'theorem hs_ok : hostSchedK %d %d hrowsL browsL hpadL hremL GXD.U8XM.schedRes = true := by decide +kernel\n' % (R, nb)
t += 'theorem pad_len : hpadL.length = %d := by decide +kernel\n' % nb
t += 'theorem rem_len : hremL.length = 2016 := by decide +kernel\n'
t += 'theorem remBound : wsumK (SSC.wGt %d) hremL 0 = 0 := by decide +kernel\n' % rtop
t += ''.join('theorem rankR%d_ok : wsumK (SSC.wEq %d) hremL 0 = %d := by decide +kernel\n' % (r, r, rc[r]) for r in range(rtop + 1))
t += 'theorem padBound : wsumK (SSC.wGt %d) hpadL 0 = 0 := by decide +kernel\n' % ptop
t += ''.join('theorem rankP%d_ok : wsumK (SSC.wEq %d) hpadL 0 = %d := by decide +kernel\n' % (r, r, pc[r]) for r in range(ptop + 1))
t += 'end GXD.%s\n#print axioms GXD.%s.hs_ok\n' % (NAME, NAME)
open(os.path.join(D, 'HU8XMSched.lean'), 'w').write(t)
body, plan = fg.block_text(NAME + 'x', u, A, 10, 40)
rate = 'import Work.BlockAccounting.Moment\nset_option maxRecDepth 4000\nnamespace OAI.PowerSaving.BlockAccounting\nopen Finset\n'
rate += fg.list_text(NAME + 'x', u, '72 copies of the U8 word (#869/#870/#890/#942) with hub exit frames, mutable sources and 5 '
                     'centre-source copies, #401 banks for the helpers that are not hosted, same-run idle-axis hosting (DaysSky #442) '
                     'of the others; fs-price.')
rate += '\n' + body + '\nend OAI.PowerSaving.BlockAccounting\n#print axioms OAI.PowerSaving.BlockAccounting.fold_%sx\n' % NAME
open(os.path.join(D, 'Rate%sx.lean' % NAME), 'w').write(rate)
remp = ' + '.join('%d * bcost φ %d' % (n, r) for r, n in sorted(rc.items()) if r)
padp = ' + '.join('%d * bcost φ %d' % (n, r) for r, n in sorted(pc.items()) if r)
inst = f'''import Work.GCert.Data.InstU8XM
import Work.XHost.Data.Rate{NAME}x
namespace OAI.PowerSaving.BridgeRate.{NAME}
open Finset Binary Matrix RAM CB FoldRate
noncomputable section
/-- the remaining climbs of ONE run (2,016 hosts, ranks `rem`). -/
def remPrice (φ : ℕ → ℝ) : ℝ := {remp}
def padPrice (φ : ℕ → ℝ) : ℝ := {padp}
def unitPrice (φ : ℕ → ℝ) : ℝ := 72 * BridgeRate.U8XM.invPrice φ + 8 * remPrice φ + padPrice φ
theorem unit_cost_x (φ : ℕ → ℝ) : unitPrice φ = listCost φ BlockAccounting.blocks{NAME}x := by
  unfold unitPrice BridgeRate.U8XM.invPrice remPrice padPrice
  simp only [listCost, BlockAccounting.blocks{NAME}x, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  norm_num [bcost]
  try ring
end
end OAI.PowerSaving.BridgeRate.{NAME}
#print axioms OAI.PowerSaving.BridgeRate.{NAME}.unit_cost_x
'''
open(os.path.join(D, 'Inst%s.lean' % NAME), 'w').write(inst)
rsum = ', '.join('rankR%d_ok' % r for r in range(rtop + 1))
psum = ', '.join('rankP%d_ok' % r for r in range(ptop + 1))
glue = f'''import Work.XHost.HBank
import Work.XHost.HPriceK
import Work.XHost.Data.HU8XMSched
import Work.XHost.Data.Inst{NAME}
import Work.GCert.Data.U8XMGlue401

/-!
# (key: bridge-host) the instance `{NAME}`: k-combine's U8XM certificate `GXD.U8XM` through the HOSTED network

Same gauged invocation package (`GXD.U8XM.ginv = ginvX.toE'`: the U8 word of Kothe #869/#870, Boukhalfa #890, Capalbo #942,
hub exit frames, x-mutable's mutable sources with 5 centre-source copies, Mamano's gauges);
host schedule: h-L3-design's same-run host table (`host_price.py --samerun 1`, after x-hosting / DaysSky's hosting_core
layout) with the {nb} banks of s-stack's re-pack, checked by `XHost.hostSchedK`.  Unit: W = {W}, kernel numerator {A}.
Written by fs-price/tools/gen_hu8xm.py (after h-L3-build's gen_h588m.py; wht43, 13 Oct 2026), prepared with assistance from Claude.  No `sorry`.
-/

set_option maxRecDepth 100000
namespace GXD.{NAME}
open GXD OAI.PowerSaving OAI.PowerSaving.CB OAI.PowerSaving.XHost
noncomputable section

/-- the host schedule (checked table) -/
def hsched : HostSched (Fin GXD.U8XM.raw.R) (Fin GXD.U8XM.raw.v) :=
  HostSched.ofTab {R} {nb} hrowsL browsL hpadL hremL GXD.U8XM.schedRes hs_ok

theorem hfits : hsched.fits GXD.U8XM.ginv.ent := by
  unfold HostSched.fits
  have e : (fun q : Fin GXD.U8XM.raw.R => Fintype.card (Fin GXD.U8XM.aux.p.h) - (GXD.U8XM.ginv.ent q).dim)
      = (fun q => GXD.U8XM.schedRes.getD q.val 0) := funext GXD.U8XM.hent
  rw [e]
  exact HostSched.ofTab_fitsR {R} {nb} hrowsL browsL hpadL hremL GXD.U8XM.schedRes hs_ok

theorem hwidth : 64 * Fintype.card (Fin GXD.U8XM.raw.v) + hsched.nb = {W} := by
  rw [Fintype.card_fin]
  rfl

theorem hrem_price (φ : ℕ → ℝ) : hcostS hsched 0 φ + hcostS hsched 1 φ + hcostS hsched 2 φ + hcostS hsched 3 φ
    + hcostS hsched 4 φ + hcostS hsched 5 φ + hcostS hsched 6 φ + hcostS hsched 7 φ + hcostS hsched 8 φ
      = rcost φ hremL := by
  rw [← sum_fin9 (fun c => hcostS hsched c φ)]
  unfold hcostS
  refine (remTot (fun i => bcost φ (hremL.getD i 0))).trans ?_
  rw [← rem_len]
  exact rcost_range φ hremL

theorem hpad_price (φ : ℕ → ℝ) : ∑ b, bcost φ (hsched.pad b) = rcost φ hpadL :=
  rcost_fin φ hpadL {nb} pad_len

/-- the hosted unit price of `{NAME}` -/
theorem hprice (φ : ℕ → ℝ) : unitPriceHost GXD.U8XM.ginv hsched φ = BridgeRate.{NAME}.unitPrice φ := by
  unfold unitPriceHost
  rw [GXD.U8XM.ginv_cost, hrem_price, hpad_price, GXD.rcost_wsumK φ hremL {rtop} remBound,
    GXD.rcost_wsumK φ hpadL {ptop} padBound]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, {rsum}, {psum}]
  unfold BridgeRate.{NAME}.unitPrice BridgeRate.{NAME}.remPrice BridgeRate.{NAME}.padPrice BridgeRate.U8XM.invPrice
    GXD.U8XM.invPrice
  norm_num [bcost]
  try ring

end
end GXD.{NAME}

#print axioms GXD.{NAME}.hfits
#print axioms GXD.{NAME}.hprice
'''
open(os.path.join(OUT, 'Work/XHost/%s.lean' % NAME), 'w').write(glue)
json.dump(dict(name=NAME, kernel_numerator=A, dft_numerator=A - 1, W=W, banks=nb, hist={str(k): v for k, v in sorted(hist.items())},
               pads=dict(Counter(pad)), rem_hist={str(k): v for k, v in sorted(rc.items())}, slack=plan['slack_rel'],
               source=dict(hosttab=HT, repack=RP, price='k-combine/out/price/U8XM.json', sched='k-combine/out/global-U8XM/Work/GCert/Data/Gen/U8XMSched.lean')), open(os.path.join(OUT, '..', 'out', NAME + '.json'), 'w'), indent=1)
print('written; rate slack', plan['slack_rel'])
