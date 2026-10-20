"""epscheck.py -- o3-finish: the pair test for prefixes whose M_D1(B') is NONZERO (fulltestG.py verdict 'nonzero').
A 2-in-7 needs rows v6 in M_E6, v7 in M_E7 forming an invertible 2 x 2 matrix A over Q[B'].  The augmentation
eps : Q[B'] -> Q (sum of coefficients) is a ring map, so eps(A) must be invertible.  Every M_E (E one move) is
b M_D1, and M_D1 is the Q[Lambda]-span of the solutions of the representatives' components; eps(b v) = eps(v).  So
if  dim eps(M_D1) <= 1  (eps applied to the two coordinates), eps(A) is singular for every choice: no 2-in-7.
This script recomputes each representative's component EXACTLY (Fractions; these components are small), takes the
T-projection of its nullspace, and reports the span of (eps_0, eps_1) over all representatives.
usage: python3 epscheck.py PREFIXES.jsonl RUN IT [RUN IT ...]"""
import sys, json, collections
from fractions import Fraction as Fr
args = sys.argv[1:]; sys.argv = ['x']
import fulltest as F, fulltestG as FG
from fulltest import *

def nullspace(A, nr, nc):
    M = [[Fr(0)] * nc for _ in range(nr)]
    for (r, c), v in A.items():
        M[r][c] += v
    piv = []; pr = 0
    for c in range(nc):
        p = next((r for r in range(pr, nr) if M[r][c] != 0), None)
        if p is None: continue
        M[pr], M[p] = M[p], M[pr]
        inv = 1 / M[pr][c]; M[pr] = [x * inv for x in M[pr]]
        for r in range(nr):
            if r != pr and M[r][c] != 0:
                f = M[r][c]; M[r] = [x - f * y for x, y in zip(M[r], M[pr])]
        piv.append(c); pr += 1
    free = [c for c in range(nc) if c not in piv]
    out = []
    for fc in free:
        v = [Fr(0)] * nc; v[fc] = Fr(1)
        for i, pc in enumerate(piv):
            v[pc] = -M[i][fc]
        out.append(v)
    return out

def component_system(ctx, E, start, maxrows=20000):
    G = ctx.G; gens = ctx.gens
    Einv = np.conj(E).T
    rows = {start: 0}; cols = {('T', start): 0}; A = {(0, 0): Fr(-1)}
    q = collections.deque([start]); seen = set()
    while q:
        t = q.popleft()
        if t in seen: continue
        seen.add(t)
        if len(rows) > maxrows: return None
        lt, gt = t; Mg = G.MATS[gt]
        cands = []
        for (i, gterm) in ctx.index.get((G.frm(gt), lt), ()):
            cands.append(('S', i, G.mul(Mg, np.conj(G.MATS[gterm]).T), None))
        hm = G.mul(Einv, Mg)
        for (i, gterm) in ctx.index.get((G.frame_key(hm), lt), ()):
            cands.append(('E', i, G.mul(hm, np.conj(G.MATS[gterm]).T), E))
        for (kind, i, bmat, pre) in cands:
            if not G.is_monomial(bmat): continue
            ck = (kind, i, G.intern(bmat))
            if ck in cols: continue
            c = len(cols); cols[ck] = c
            P2 = bmat if pre is None else G.mul(pre, bmat)
            for (l2, g2), co in gens[i].items():
                g3 = G.intern(G.mul(P2, G.MATS[g2]))
                if G.is_free(g3) or G.is_monomial(G.mul(Einv, G.MATS[g3])): continue
                t3 = (l2, g3)
                if t3 not in rows:
                    rows[t3] = len(rows)
                    if ctx.seed_of(g3):
                        cols[('T', t3)] = len(cols); A[(rows[t3], cols[('T', t3)])] = Fr(-1)
                    q.append(t3)
                A[(rows[t3], c)] = A.get((rows[t3], c), Fr(0)) + Fr(co)
    return rows, cols, A

def eps_span(rec):
    w = W4(1); G = w.G
    gens = build(w, rec["prog"])
    ctx = FG.CtxG(w, gens)
    d1 = {0: [], 1: []}
    for g in gens:
        for (l, gid) in g:
            if w.dist(gid) == 1 and gid not in d1[l]: d1[l].append(gid)
    vecs = []; info = []
    for l in range(2):
        for gid2 in d1[l]:
            h = G.mul(G.K, np.conj(G.MATS[gid2]).T)
            if not (G.rowsupp(h) == 2).all(): continue
            d = FG.decomp(G, h)
            if d is None: continue
            sid = G.intern(G.mul(np.conj(d[0]).T, G.K))
            cs = component_system(ctx, G.D1, (l, sid))
            if cs is None:
                info.append((l, 'too big')); continue
            rows, cols, A = cs
            ns = nullspace(A, len(rows), len(cols))
            tcols = [(k[1], c) for k, c in cols.items() if k[0] == 'T']
            for v in ns:
                e = [sum(v[c] for (row, c) in tcols if row[0] == ll) for ll in range(2)]
                vecs.append(e)
            info.append((l, len(rows), len(cols), len(tcols), len(ns)))
    # rank of the eps vectors
    r = 0; basis = []
    for e in vecs:
        e = list(e)
        for b in basis:
            k = next(i for i in range(2) if b[i] != 0)
            if e[k] != 0:
                f = e[k] / b[k]; e = [x - f * y for x, y in zip(e, b)]
        if any(x != 0 for x in e):
            basis.append(e)
    return len(basis), [[str(x) for x in e] for e in vecs[:8]], info

if __name__ == "__main__":
    recs = {}
    for l in open(args[0]):
        r = json.loads(l); recs[(r['run'].replace('.jsonl', ''), r['it'])] = r
    for i in range(1, len(args), 2):
        k = (args[i], int(args[i + 1]))
        d, ex, info = eps_span(recs[k])
        print(k, "dim eps(M_D1) =", d, "-> no invertible pair (no 2-in-7)" if d <= 1 else "-> eps does not decide", "| eps vectors", ex, "| components", info, flush=True)
