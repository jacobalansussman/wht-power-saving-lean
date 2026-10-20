import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import poolsuff as PS
from poolsuff import *
H = FT.hdata()
want = {(x.split(':')[0], int(x.split(':')[1])) for x in sys.argv[2].split(',')}
for l in open(sys.argv[1]):
    rec = json.loads(l); key = (rec['run'].replace('.jsonl', ''), rec['it'])
    if key not in want: continue
    t0 = time.time()
    Hp, Hq, HisR, Hm, Hkeyidx = H
    S = Sys(rec); G = S.G
    kgen = [t for t in range(2, len(S.gens)) if any(G.frm(gid) == S.Kfrm for (l2, gid) in S.gens[t])]
    R0 = {}
    for t in kgen:
        r = dead_test(S, ('S', t, np.eye(16, dtype=complex)), False, 40000)
        R0[t] = (r[0], len(r[3] or []))
    w = FT.W4(1); G2 = w.G; gens = FT.build(w, rec['prog']); ctx = FT.Ctx(w, gens)
    pool = FT.e_classes(w, gens, Hp, Hq, HisR)
    print(key, 'kgens', kgen, 'R0', R0, 'pool classes', len({k for (_, _, k) in pool}), round(time.time() - t0, 1), 's', flush=True)
    poolkeys = {k for (_, _, k) in pool}
    threats = {}
    for t in kgen:
        rr = dead_test(S, ('S', t, np.eye(16, dtype=complex)), False, 40000)
        R0m = [(l2, S.G.MATS[gid].copy()) for (l2, gid) in rr[3]]
        cterms = [G2.MATS[gid] for (l2, gid) in gens[t] if G2.frm(gid) == ctx.Kfrm]
        for cK in cterms:
            ci = np.conj(G2.mul(cK, np.conj(G2.K).T)).T
            for (l2, Mr) in R0m:
                for s2 in range(len(gens)):
                    for (l3, g2) in gens[s2]:
                        if l3 != l2: continue
                        h = G2.mul(Mr, np.conj(G2.MATS[g2]).T)
                        if decomp1(h) is None: continue
                        E = G2.mul(ci, h)
                        k = PS.pclass_key(E, Hp, Hq, HisR, Hm)
                        threats.setdefault(k, 0); threats[k] += 1
    print('   threat classes', len(threats), 'of which pool', sum(1 for k in threats if k in poolkeys), 'K-frame terms', len(cterms), flush=True)
