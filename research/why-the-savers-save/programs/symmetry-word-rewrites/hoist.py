"""Maximal hoist set: scanning A forwards, a gate joins H if every register it touches starts at the empty frame and
it commutes with every earlier gate not in H; H moves to the front of A (all at frame ()), which is free."""
import sys, copy
sys.path.insert(0, '.')
import fw, search, reorder2
from sink import reframe
w = fw.load(sys.argv[1] if len(sys.argv) > 1 else 'K8-FULL-WORD.json.gz')
W = search.Word(w); O = reorder2.OW2(W)
H = []; keep = []
for g in O.order:
    if g == O.CUT: break
    if all(W.start[r] == () for r in O.regs[g]) and all(O.commute(g, h) for h in keep):
        H.append(g)
    else: keep.append(g)
already = sum(1 for g in H if W.f[g] == ())
print('hoist set', len(H), 'already at ()', already)
w2 = copy.deepcopy(w)
Hs = set(W.seq[g][1] for g in H)
w2['A'] = [w['A'][i][:-1] + [[]] for i in sorted(Hs)] + [e for i, e in enumerate(w['A']) if i not in Hs]
w3, info = reframe(w2)
print(info)
if w3: print('profile', w3['report']['local_profile']); fw.save(w3, 'W-hoist.json.gz')
