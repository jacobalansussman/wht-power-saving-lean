#!/usr/bin/env python3
"""heightfilter.py (p6-prove, 15 Oct 2026): exact search filter from Theorem C (PROOF.md section 3).

A word CANNOT SAVE if there is a height system: integers h_l(x) (input l, frame x of Gamma_m) with
  h_l(L0) = 0, h_l(L_T) = m, |h_l(x) - h_l(y)| <= d(x, y)   (then h_l extends to a 1-Lipschitz function on Gamma_m),
such that after every paid move all terms (l, x) of the moved register have one common height.
This is a system of difference constraints; it is feasible iff the constraint graph has no negative cycle
(Bellman-Ford), and then an integer solution exists.  Exact and polynomial.

Interface:
  feasible(F, outputs, W) -> (bool, heights or None)   outputs = list of sets {(l, x)}, x = frame index of Frames(m)
  word_outputs(E, F, W, R, word) -> outputs             for words in ledger.py format (true-model terms)
"""
import sys
import os as _os; sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "ledger"))
import numpy as np

def feasible(F, outputs, W):
    nodes = {("Z",)}
    for l in range(W): nodes |= {(l, F.L0), (l, F.LT)}
    for out in outputs: nodes |= set(out)
    nodes = sorted(nodes, key=repr); idx = {v: i for i, v in enumerate(nodes)}
    E = []   # (u, v, w): h(v) - h(u) <= w
    z = idx[("Z",)]
    def eq(u, v, c):   # h(v) - h(u) = c
        E.append((u, v, c)); E.append((v, u, -c))
    for l in range(W):
        eq(z, idx[(l, F.L0)], 0); eq(z, idx[(l, F.LT)], F.m)
    byl = {}
    for v in nodes:
        if v[0] != "Z": byl.setdefault(v[0], []).append(v)
    for l, vs in byl.items():
        for a in vs:
            for b in vs:
                if a != b: E.append((idx[a], idx[b], F.d(a[1], b[1])))
    for out in outputs:
        out = sorted(out, key=repr)
        for a, b in zip(out, out[1:]): eq(idx[a], idx[b], 0)
    n = len(nodes); dist = [0] * n
    for it in range(n + 1):
        ch = False
        for (u, v, w) in E:
            if dist[u] + w < dist[v]: dist[v] = dist[u] + w; ch = True
        if not ch: break
    else:
        return False, None
    if ch: return False, None
    h = {v: dist[idx[v]] - dist[z] for v in nodes if v[0] != "Z"}
    return True, h

def best_bound(F, outputs, W):
    """quantitative Theorem C: max over height systems (h_l(L0) = 0, 1-Lipschitz, every move output homogeneous)
    of sum_l h_l(L_T); M >= this for every unit with these move outputs.  The pointwise-maximal solution of a
    system of difference constraints is the shortest-path distance from the zero node, so the maximum of the sum
    is the sum of the distances to the nodes (l, L_T).  Returns None if even h = 0 ... never infeasible (h = 0 works)."""
    nodes = {("Z",)}
    for l in range(W): nodes |= {(l, F.L0), (l, F.LT)}
    for out in outputs: nodes |= set(out)
    nodes = sorted(nodes, key=repr); idx = {v: i for i, v in enumerate(nodes)}
    E = []
    z = idx[("Z",)]
    def eq(u, v, c): E.append((u, v, c)); E.append((v, u, -c))
    for l in range(W): eq(z, idx[(l, F.L0)], 0)
    byl = {}
    for v in nodes:
        if v[0] != "Z": byl.setdefault(v[0], []).append(v)
    for l, vs in byl.items():
        for a in vs:
            for b in vs:
                if a != b: E.append((idx[a], idx[b], F.d(a[1], b[1])))
    for out in outputs:
        out = sorted(out, key=repr)
        for a, b in zip(out, out[1:]): eq(idx[a], idx[b], 0)
    INF = 10 ** 9
    n = len(nodes); dist = [INF] * n; dist[z] = 0
    for it in range(n + 1):
        ch = False
        for (u, v, w) in E:
            if dist[u] < INF and dist[u] + w < dist[v]: dist[v] = dist[u] + w; ch = True
        if not ch: break
    assert not ch, "negative cycle impossible: h = 0 is feasible"
    return sum(min(dist[idx[(l, F.LT)]], F.m) for l in range(W))

def position(E, F, gid):
    """frame lambda(g) = gbar^-1 L0 as a Frames index; syms columns = images of basis labels, a-bits first"""
    from ledger import inv2
    S = inv2(np.array(E.syms[gid]) % 2)
    m = E.m
    vecs = []
    for i in range(m):                      # L0 basis (0 | e_i) = label index m + i
        col = S[:, m + i] % 2
        a = sum(int(col[k]) << k for k in range(m)); b = sum(int(col[m + k]) << k for k in range(m))
        vecs.append(a | (b << m))
    return F.idx[F.span(vecs)]

def word_outputs(E, F, W, R, word):
    from t_chain import tstep
    from ledger import ONE
    C = [dict() for _ in range(R)]
    for r in range(W): C[r] = {(r, E.ID, 0): ONE}
    outs = []
    for st in word:
        C = tstep(E, C, st)
        if st[0] == "move":
            outs.append({(l, position(E, F, g)) for (l, g, s) in C[st[1]]})
    return outs

if __name__ == "__main__":
    from ledger import Engine
    from frames import Frames
    import checks
    for m in (2, 3):
        E = Engine(m); F = Frames(m)
        assert position(E, F, E.ID) == F.L0 and position(E, F, E.Kid) == F.LT, "frame convention"
        cases = [("plain", 2, 2, checks.plain(m, 2))]
        w = [("move", 0, 1), ("mix", 0, [(0, 1, None), (1, 1, None)])]
        for j in range(1, m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
        w += [("mix", 0, [(0, 1, None), (1, -1, None)]), ("move", 1, 1)]
        cases.append(("gaugeX", 2, 2, w))
        s = E.addr([2, 1] + [1 << j for j in range(2, m)])
        w = [("mix", 0, [(0, 1, None), (1, 1, None)]), ("free", 1, s), ("mix", 1, [(1, 1, None), (2, 1, None)]),
             ("mix", 2, [(2, 1, None), (0, 1, None)])] + checks.plain(m, 3)
        w += [("mix", 2, [(2, 1, None), (0, -1, None)]), ("mix", 1, [(1, 1, None), (2, -1, None)]), ("free", 1, s),
              ("mix", 0, [(0, 1, None), (1, -1, None)])]
        cases.append(("HOLO3", 3, 3, w))
        if m >= 3:
            from fractions import Fraction as Fr
            s3 = E.addr([2, 4, 1] + [1 << j for j in range(3, m)]); s32 = E.addr([4, 1, 2] + [1 << j for j in range(3, m)])
            w = [("mix", 0, [(0, 1, None), (1, 1, None)]), ("mix", 1, [(0, 1, None), (1, -1, None), (1, -1, s3)])]
            for j in range(m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
            h = Fr(1, 2)
            w += [("mix", 1, [(0, h, None), (1, -h, None), (0, -h, s3), (1, h, s3), (0, h, s32), (1, -h, s32)]),
                  ("mix", 0, [(0, 1, None), (1, -1, None)])]
            cases.append(("C%dO3" % (2 * m), 2, 2, w))
        # a crossing-type word: input 1 first moved alone along a line, then joined (inhomogeneous in raw levels)
        for name, W, R, w in cases:
            outs = word_outputs(E, F, W, R, w)
            ok, h = feasible(F, outs, W)
            bb = best_bound(F, outs, W)
            print("m=%d %-7s moves=%d  height system exists: %s  -> %s;  best Theorem-C bound M >= %d (W*m = %d)"
                  % (m, name, len(outs), ok, "cannot save (Theorem C)" if ok else "not excluded by Theorem C",
                     bb, W * m))
