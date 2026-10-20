#!/usr/bin/env python3
"""cx-mirror / chaincost.py (own code, stdlib).  First-order price of cleaning (and, by the same count, gathering) one
total's content against its 56 x's when no drop-free chain covers the support.

Model (optimistic on purpose): one register holds the total's content sigma_k at frame 0 after the U-drop and removes
it chain by chain; a chain is a set of x's it can meet without dropping (for every twisted mirror, N is symmetric
and a chain is an independent set of {N = 1}); the x's frames are taken as small as possible (only their own lines).
Before chain C_{i+1} the register must drop from its frame F to F & (meet of v_b^perp over b in C_{i+1}); the drop
costs the codimension.  Greedy: repeatedly pick, among maximum chains of the remaining support (exact maximum
independent set, ties broken at random over restarts), the one with the cheapest drop.  Cases: the natural mirror
and the best twisted mirror of twist3.json.  Reports drops per total and per call (x2 for gather and cleanup)."""
import os, json, random, importlib.util
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('tw', os.path.join(here, 'twist_lib.py')); tw = importlib.util.module_from_spec(spec); spec.loader.exec_module(tw)
U, Sk, v, dot, apply, I9 = tw.U, tw.Sk, tw.v, tw.dot, tw.apply, tw.I9
S10 = 1 << 9; J = (1 << 9) - 1
def red(vs):
    b = []
    for x in vs:
        for y in b: x = min(x, x ^ y)
        if x: b.append(x)
    return b
def dim(vs): return len(red(vs))
def meet_perp(F, ws):
    """dim of {f in span F : f . w = 0 for all w in ws}: = dim F - rank of the matrix [f_i . w_j]"""
    Fb = red(F)
    rows = [sum(dot(f, w) << j for j, w in enumerate(ws)) for f in Fb]
    return len(Fb) - dim(rows), Fb
def sub_perp(F, ws):
    """a basis of the subspace of span F orthogonal to every w in ws"""
    Fb = red(F); out = []
    # Gaussian elimination on the functionals
    vecs = [(f, sum(dot(f, w) << j for j, w in enumerate(ws))) for f in Fb]
    piv = []
    for i in range(len(vecs)):
        f, s = vecs[i]
        for (pf, ps) in piv:
            if s & (ps & -ps): f ^= pf; s ^= ps
        if s: piv.append((f, s))
        else: out.append(f)
        # keep pivots reduced
    return red(out)
def mis_all(nodes, adj, rnd, k=40):
    """several maximum independent sets (randomised branch and bound order)"""
    res = []
    for _ in range(k):
        order = sorted(nodes, key=lambda x: (-len(adj[x] & set(nodes)), rnd.random()))
        best = []
        def bb(Rr, Pp):
            nonlocal best
            if len(Rr) + len(Pp) <= len(best): return
            if not Pp: best = list(Rr); return
            x = Pp[0]
            bb(Rr + [x], [y for y in Pp[1:] if y not in adj[x]])
            bb(Rr, Pp[1:])
        bb([], order); res.append(best)
    return res
def cost_case(A, label, restarts=6, seed=0):
    rnd = random.Random(seed)
    Au = [apply(A, u) for u in U]
    V = [S10 | (J ^ x) for x in Au]
    adj = {a: {b for b in range(v) if a != b and dot(U[b], Au[a]) == 0} for a in range(v)}   # N = 1 + u_b . A u_a
    per_total = []
    for S in Sk:
        bestc = None
        for _ in range(restarts):
            rest = list(S); F = []; cost = 0; nch = 0
            while rest:
                cands = mis_all(rest, adj, rnd, k=12)
                def drop_of(C):
                    return len(red(F)) - len(sub_perp(F, [V[b] for b in C])) if F else 0
                C = min(cands, key=lambda c: (drop_of(c), rnd.random()))
                d = drop_of(C); cost += d; nch += 1
                F = sub_perp(F, [V[b] for b in C]) + [U[b] for b in C]
                rest = [x for x in rest if x not in C]
            if bestc is None or cost < bestc[0]: bestc = (cost, nch)
        per_total.append(bestc)
    tot = sum(c for c, _ in per_total)
    print('%-28s chains per total %s | cleanup drops per total %s | per call: cleanup %d, gather + cleanup %d (room 34)'
          % (label, [n for _, n in per_total], [c for c, _ in per_total], tot, 2 * tot))
    return dict(case=label, chains=[n for _, n in per_total], drops=[c for c, _ in per_total], per_call_cleanup=tot, per_call_both=2 * tot)
out = [cost_case(I9, 'natural mirror')]
t3 = json.load(open(os.path.join(here, 'twist3.json')))
out.append(cost_case(t3['best']['A'], 'best sampled twisted mirror'))
json.dump(out, open(os.path.join(here, 'chaincost.json'), 'w'), indent=1)
