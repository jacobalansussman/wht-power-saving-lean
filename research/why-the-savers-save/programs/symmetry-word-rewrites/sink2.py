"""Maximal sink set: scanning main+tail backwards, a gate joins S if both registers end full and it commutes with
every later gate not in S (later S gates keep their relative order in the cleanup prefix)."""
import sys, copy
sys.path.insert(0, '.')
import fw, search, reorder2
from sink import reframe
w = fw.load(sys.argv[1] if len(sys.argv) > 1 else 'K8-FULL-WORD.json.gz')
W = search.Word(w); O = reorder2.OW2(W)
S = []; keep = []
for g in reversed(O.order):
    if g == O.CUT: keep.append(g); continue
    if O.sec0[g] != 'T' and all(W.end[r] == W.full for r in O.regs[g]) and all(O.commute(g, h) for h in keep):
        S.append(g)
    else: keep.append(g)
S.reverse()
print('sink set size', len(S), 'by section', {s: sum(1 for g in S if O.sec[g] == s) for s in 'AB'})
full = list(fw.basis([1 << i for i in range(9)]))
w2 = copy.deepcopy(w)
Aset = {W.seq[g][1] for g in S if O.sec[g] == 'A'}; Bset = {W.seq[g][1] for g in S if O.sec[g] == 'B'}
moved = [ (w['A'] if O.sec[g] == 'A' else w['B'])[W.seq[g][1]][:-1] + [full] for g in S]
w2['A'] = [e for i, e in enumerate(w['A']) if i not in Aset]
w2['B'] = [e for i, e in enumerate(w['B']) if i not in Bset]
m = len(moved); w2['cleanup'] = moved + w['cleanup']
w2['early_cleanup'] = [[i + m, U] for i, U in w['early_cleanup']]
w3, info = reframe(w2)
print(info)
if w3:
    print('profile', w3['report']['local_profile']); fw.save(w3, 'W-sinkmax.json.gz')
