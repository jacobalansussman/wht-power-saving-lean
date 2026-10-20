"""h-L3-design: emit the per-run host table of a SAME-RUN hosting result (host_price.py --samerun 1) as a Lean data
module checked by XHost.hostK (`decide +kernel`), plus a JSON copy.

Slot order q: widths descending, then helper index h (the residual of a width-w helper is w; the checker only needs
res q).  Host index i = (c * 2 + side) * nT + t, c in CH order (A2 B2 A3 B3 A4 B4 F1 F2 F3), side = inst // nT,
t = inst % nT.  Entry 0 = not hosted (bank).  Asserts: same table in every run, one host per (stage, slot).
usage: python3 -I emit_hosttab.py WORDS.json RESULT.json NAME OUT.lean OUT.json [--runs-check 1]"""
import sys, json, gzip
from collections import Counter, defaultdict
CH = ['A2', 'B2', 'A3', 'B3', 'A4', 'B4', 'F1', 'F2', 'F3']
RANK = dict(A2=16, B2=8, A3=32, B3=8, A4=48, B4=8, F1=48, F2=32, F3=16)
words, resf, name, outl, outj = sys.argv[1:6]
w = json.load(open(words))[name]; res = json.load(open(resf))
assert res.get('samerun') == 1 and res['copies'] == 0
wid = Counter({int(a): b for a, b in w['widths'].items()})
slots_order = [(ww, h) for ww in sorted(wid, reverse=True) for h in range(wid[ww])]
qof = {s: q for q, s in enumerate(slots_order)}
R = len(slots_order); nT = 112
slots = json.load(gzip.open(resf + '.slots.gz', 'rt'))
per = defaultdict(set)
for k, ww, h, (rp, ch), inst, off in slots:
    assert k // 9 == rp, 'cross-run slot'
    per[rp].add((k % 9, ww, h, ch, inst))
assert all(per[r] == per[0] for r in range(8)), 'table differs between runs'
rows = [[0] * R for _ in range(9)]
load = Counter()
for l, ww, h, ch, inst in sorted(per[0]):
    c = CH.index(ch); side, t = divmod(inst, nT)
    i = (c * 2 + side) * nT + t
    q = qof[(ww, h)]
    assert rows[l][q] == 0
    rows[l][q] = i + 1
    load[i] += ww
rem = []
for i in range(18 * nT):
    c = i // (2 * nT)
    rem.append(RANK[CH[c]] - load[i])
    assert rem[-1] >= 0
resl = [ww for ww, h in slots_order]
# cross-check against the result's child remainders (all 8 replicas)
chk = {ch: Counter() for ch in CH}
for i, m in enumerate(rem):
    chk[CH[i // (2 * nT)]][m] += 8
for ch in CH:
    want = Counter({int(a): b for a, b in res['child_remainders'][ch].items()})
    assert chk[ch] == want, (ch, chk[ch], want)
hosted = sum(1 for row in rows for e in row if e)
json.dump(dict(name=name, nT=nT, R=R, rows=rows, rem=rem, res=resl, hosted_per_run=hosted), open(outj, 'w'))
mod = outl.split('/')[-1].replace('.lean', '')
L = []
L.append('import XHost.HostK\n')
L.append('/-! Host table of one run (same in all 8 runs) for `%s`, from %s (h-L3-design, emit_hosttab.py).\n'
         '%d slots, %d hosted entries per run, %d hosts. Data only; checked by `hostK` with `decide +kernel`. -/\n' % (name, resf.split('/')[-1], R, hosted, 18 * nT))
L.append('namespace OAI.PowerSaving.XHost.%s\n' % mod)
for l in range(9):
    L.append('def row%d : List Nat := [%s]\n' % (l, ', '.join(map(str, rows[l]))))
L.append('def rows : List (List Nat) := [%s]\n' % ', '.join('row%d' % l for l in range(9)))
L.append('def rem : List Nat := [%s]\n' % ', '.join(map(str, rem)))
L.append('def res : List Nat := [%s]\n' % ', '.join(map(str, resl)))
L.append('theorem host_ok : hostK %d %d rows rem res = true := by decide +kernel\n' % (nT, R))
hist = Counter(rem)
for m in sorted(hist):
    L.append('theorem hist_%d : hkHist %d rem = %d := by decide +kernel\n' % (m, m, hist[m]))
L.append('end OAI.PowerSaving.XHost.%s\n' % mod)
open(outl, 'w').write('\n'.join(L))
print(f'{name}: R={R} hosted/run={hosted} hosts={18*nT} rem-hist={dict(sorted(hist.items()))} -> {outl}')
