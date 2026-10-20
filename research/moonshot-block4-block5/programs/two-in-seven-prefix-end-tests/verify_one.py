"""verify_one.py -- o3-finish: take a prefix with M_D1 != 0, extract one exact solution, and check the identity
    sum_j c_j D1 b_j gen_ij + sum_k d_k b_k gen_ik + (free terms) = sum_s t_s (b_s K) x_l(s)
formally (term by term) and on random complex arrays at f = 1..5 (g -> g^(x)f).  usage: verify_one.py RUN IT"""
import sys, json
args = sys.argv[1:]; sys.argv = ['x']
import fulltest as F, fulltestG as FG, epscheck as EC
from fulltest import *
from fractions import Fraction as Fr
recs = {}
for l in open('all_prefixes.jsonl'):
    r = json.loads(l); recs[(r['run'].replace('.jsonl', ''), r['it'])] = r
rec = recs[(args[0], int(args[1]))]
w = W4(1); G = w.G
gens = build(w, rec['prog'])
ctx = FG.CtxG(w, gens)
d1 = {0: [], 1: []}
for g in gens:
    for (l, gid) in g:
        if w.dist(gid) == 1 and gid not in d1[l]: d1[l].append(gid)
Ki = np.conj(G.K).T
done = False
for l in range(2):
    for gid2 in d1[l]:
        h = G.mul(G.K, np.conj(G.MATS[gid2]).T)
        d = FG.decomp(G, h)
        if d is None: continue
        sid = G.intern(G.mul(np.conj(d[0]).T, G.K))
        rows, cols, A = EC.component_system(ctx, G.D1, (l, sid))
        ns = EC.nullspace(A, len(rows), len(cols))
        if not ns: continue
        v = ns[0]; ic = {c: k for k, c in cols.items()}
        # formal identity: sum over non-T columns of coef * column terms == sum over T columns of t * seed row, up to free rows
        lhs = {}
        for c, x in enumerate(v):
            if x == 0: continue
            k = ic[c]
            if k[0] == 'T':
                continue
            kind, i, b = k
            P2 = G.MATS[b] if kind == 'S' else G.mul(G.D1, G.MATS[b])
            for (l2, g2), co in gens[i].items():
                t = (l2, G.intern(G.mul(P2, G.MATS[g2])))
                lhs[t] = lhs.get(t, 0) + x * co
        rhs = {}
        for c, x in enumerate(v):
            k = ic[c]
            if k[0] == 'T' and x != 0:
                rhs[k[1]] = rhs.get(k[1], 0) + x
        diff = dict(lhs)
        for t, x in rhs.items():
            diff[t] = diff.get(t, 0) - x
        bad = [(t, x) for t, x in diff.items() if x != 0 and not (G.is_free(t[1]) or G.is_monomial(G.mul(np.conj(G.D1).T, G.MATS[t[1]])))]
        eps = [sum(x for t, x in rhs.items() if t[0] == ll) for ll in range(2)]
        print('array', l, 'solution: E-columns', sum(1 for c, x in enumerate(v) if x and ic[c][0] == 'E'),
              'S-columns', sum(1 for c, x in enumerate(v) if x and ic[c][0] == 'S'), 'seeds', len(rhs), 'eps', [str(e) for e in eps],
              'formal identity (non-free leftover terms):', len(bad))
        # raw arrays: the leftover free / D1-free terms are free images of inputs (or D1 of them): include them on both sides
        res = {}
        rng = np.random.default_rng(3)
        for f in range(1, 6):
            Nf = 16 ** f
            xs = [rng.normal(size=Nf) + 1j * rng.normal(size=Nf) for _ in range(2)]
            vecs = list(xs)
            for mv in rec['prog']:
                cc = np.zeros(Nf, dtype=complex)
                for (nu, de, a, j) in mv['cs']:
                    vv = vecs[j] if a is None else apply_tensor(G.mono(*a), vecs[j], f)
                    cc = cc + (nu / de) * vv
                vecs.append(apply_tensor(G.mul(w.D1, G.mono(*mv['b'])), cc, f))
            L = np.zeros(Nf, dtype=complex)
            for c, x in enumerate(v):
                if x == 0 or ic[c][0] == 'T': continue
                kind, i, b = ic[c]
                P2 = G.MATS[b] if kind == 'S' else G.mul(G.D1, G.MATS[b])
                L += float(x) * apply_tensor(P2, vecs[i], f)
            R = np.zeros(Nf, dtype=complex)
            for t, x in diff.items():
                if x != 0:
                    R += float(x) * apply_tensor(G.MATS[t[1]], xs[t[0]], f)
            for t, x in rhs.items():
                R += float(x) * apply_tensor(G.MATS[t[1]], xs[t[0]], f)
            res[f] = float(np.abs(L - R).max() / max(1e-300, np.abs(R).max()))
        print('   raw arrays f = 1..5, max relative error:', {k: float('%.1e' % e) for k, e in res.items()})
        sd = [(t[0], G.mono_of(G.mul(G.MATS[t[1]], Ki)), str(x)) for t, x in rhs.items()]
        print('   seed letters b (row = b K x_l) and coefficients:', [(a, list(map(int, s)), list(map(int, q)), x) for (a, (s, q), x) in sd])
        done = True; break
    if done: break
