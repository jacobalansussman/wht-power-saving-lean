"""mincert.py -- o3-finish2: smallest BFS-prefix certificate of each dead off-line class (exact), and the depth /
distance-to-K profile of the junk rows it uses.  Diagnostic only."""
import sys, os, json, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import kfD
from kfD import *

def min_cert(S, start, withE, maxrows=100):
    G = S.G
    old = kfD.CHECKS
    for n in range(0, maxrows + 1):
        kfD.CHECKS = (n,) if n > 0 else (1,)
        r = dead_test(S, start, withE, max(n, 1))
        if r[0] == "dead":
            kfD.CHECKS = old
            rows = r[3]
            prof = collections.Counter((l, S.w.depth(g), S.w.dist(g)) for (l, g) in rows)
            return n, r[2], dict(prof)
    kfD.CHECKS = old
    return None

recs = {}
for l in open(sys.argv[1]):
    r = json.loads(l); recs[(r['run'].replace('.jsonl', ''), r['it'])] = r
for l in open(sys.argv[2]):
    o = json.loads(l)
    if o.get('how') != 'dead':
        continue
    rec = recs[(o['run'], o['it'])]
    S = Sys(rec); G = S.G
    # rebuild the direct classes and test the dead ones
    seen = set(); k = 0
    cl = []
    for t in range(2, len(S.gens)):
        if any(G.frm(gid) == S.Kfrm for (ll, gid) in S.gens[t]):
            cl.append(('span', t, None))
    for s in range(2, len(S.gens)):
        for (ll, gid) in S.gens[s]:
            d = decomp1(G.mul(G.K, np.conj(G.MATS[gid]).T))
            if d is None: continue
            key = (s, G.frame_key(G.mul(G.D1, d[1])))
            if key in seen: continue
            seen.add(key); cl.append(('dir', s, d[1]))
    for kk, sts in o['class_status'].items():
        c = cl[int(kk)]
        if c[0] != 'dir' or sts[0][0] != 'dead':
            continue
        res = min_cert(S, ('E', c[1], c[2]), True)
        print(o['run'], o['it'], 'class', kk, 'D1 on u%d' % (c[1] - 1), 'min certificate rows', res[0] if res else None,
              'cols', res[1] if res else None, 'profile (input, depth, dist)', res[2] if res else None, flush=True)
