#!/usr/bin/env python3
"""vmod.py (v-p6-refute): test Lemma G (p6-prove) and the permissive position-level closure against EXPLICIT modules
V in C[G]^2, with my own code.  Random single-term plays; at every step:
  S_true = every single-term vector of V (singles e_(l,g); two-input pairs e_(1,g1) + c e_(2,g2), forced or free c),
           computed from an annihilator basis of V (columns zero / parallel);
  S_pred = what the position game offers: singles at available positions; pairs with both ends available (free c);
           pairs (g1, g2) in one LOCKED tree component with g2 = g1 * (product of twists along the tree path).
Reports: missed = S_true minus S_pred (must be empty), over = S_pred minus S_true (allowed only in cycle components).
Usage: vmod.py GAME PLAYS STEPS SEED"""
import sys, random
import numpy as np
import warnings; warnings.filterwarnings('ignore')
MODE = 'full'
import vgroups as VG

def make(name):
    return {"D7": lambda: VG.dihedral(7), "D9": lambda: VG.dihedral(9), "B3": lambda: VG.hyperoct(3),
            "B4": lambda: VG.hyperoct(4), "S4Z3": VG.s4z3, "J42": VG.johnson42, "H23": lambda: VG.hamming(2),
            "H33": lambda: VG.hamming(3), "D5": lambda: VG.dihedral(5)}[name]()

class PosState:
    def __init__(self, g):
        self.g = g; self.A = [{g.raw}, {g.raw}]; self.edges = []      # edges: (t, x)
    def comps(self):
        g = self.g; N = g.N
        par = {}
        def f(a):
            while par.get(a, a) != a: a = par[a]
            return a
        for (t, x, _l) in self.edges:
            u, w = x, N + g.ract[t][x]
            ru, rw = f(u), f(w)
            if ru != rw: par[ru] = rw
        V = {}; E = {}
        for (t, x, _l) in self.edges:
            r = f(x); V.setdefault(r, set()).update({x, N + g.ract[t][x]}); E[r] = E.get(r, 0) + 1
        return f, V, E
    def close(self):
        g = self.g; N = g.N
        while True:
            f, V, E = self.comps()
            ch = False; cyc = set()
            for r, vs in V.items():
                av = any((v < N and v in self.A[0]) or (v >= N and v - N in self.A[1]) for v in vs)
                if E[r] >= len(vs): cyc.add(r)
                if MODE == 'nocycle' and not av: continue
                if not av and E[r] >= len(vs):
                    self.cycun = getattr(self, 'cycun', set()) | {(0, v) if v < N else (1, v - N) for v in vs}
                if av or E[r] >= len(vs):
                    for v in vs:
                        if v < N:
                            if v not in self.A[0]: self.A[0].add(v); ch = True
                        else:
                            if v - N not in self.A[1]: self.A[1].add(v - N); ch = True
            ne = [(t, x, l) for (t, x, l) in self.edges if x not in self.A[0]]
            if len(ne) != len(self.edges): ch = True
            self.edges = ne
            if not ch: return
    def predicted(self):
        """singles set {(l, elem)}, free pairs (both available), forced pairs {(g1, g2)}"""
        g = self.g; N = g.N
        singles = {(l, e) for l in (0, 1) for e in range(g.n) if g.pos[e] in self.A[l]}
        adj = {}
        for (t, x, _l) in self.edges:
            u, w = x, N + g.ract[t][x]
            adj.setdefault(u, []).append((w, t, +1)); adj.setdefault(w, []).append((u, t, -1))
        forced = set()
        for e1 in range(g.n):
            x = g.pos[e1]
            if x not in adj: continue
            # DFS on the tree, carrying the cover element
            seen = {x: e1}; st = [x]
            while st:
                v = st.pop(); ev = seen[v]
                for (nv, t, sgn) in adj[v]:
                    if nv in seen: continue
                    ne = g.M[ev][t] if sgn > 0 else g.M[ev][g.I[t]]
                    seen[nv] = ne; st.append(nv)
                    if nv >= N and (MODE != 'nocomp' or v == x): forced.add((e1, ne))
        return singles, forced

def annihilator(Q, tol=1e-9):
    # Q: columns orthonormal basis of V in C^D; return rows spanning V-perp
    D = Q.shape[0]
    if Q.shape[1] == 0: return np.eye(D, dtype=complex)
    U, s, Vh = np.linalg.svd(Q, full_matrices=True)
    r = int((s > tol).sum())
    return U[:, r:].conj().T

def true_sets(g, Q, tol=1e-7):
    n = g.n
    Nm = annihilator(Q)
    if Nm.shape[0] == 0:
        return {(l, e) for l in (0, 1) for e in range(n)}, set(), set()
    nrm = np.linalg.norm(Nm, axis=0)
    singles = {(l, e) for l in (0, 1) for e in range(n) if nrm[l * n + e] < tol}
    nz1 = [e for e in range(n) if nrm[e] >= tol]; nz2 = [e for e in range(n) if nrm[n + e] >= tol]
    forced = set()
    if nz1 and nz2:
        A = Nm[:, nz1] / nrm[nz1]; B = Nm[:, [n + e for e in nz2]] / nrm[[n + e for e in nz2]]
        Gm = np.abs(A.conj().T @ B)
        for i, j in zip(*np.where(Gm > 1 - 1e-7)): forced.add((nz1[i], nz2[j]))
    return singles, forced, Nm

def leftmul(g, a, v):
    n = g.n; w = np.zeros_like(v)
    for l in (0, 1):
        for e in range(n):
            if v[l * n + e] != 0: w[l * n + g.M[a][e]] = v[l * n + e]
    return w

def add_gen(g, Q, v):
    vs = [leftmul(g, b, v) for b in g.P]
    Mx = np.column_stack([Q] + [x[:, None] for x in vs]) if Q.shape[1] else np.column_stack(vs)
    U, s, _ = np.linalg.svd(Mx, full_matrices=False)
    return U[:, : int((s > 1e-9).sum())]

def run(name, plays, steps, seed):
    rng = random.Random(seed); nrng = np.random.default_rng(seed)
    g = make(name); n = g.n; D = 2 * n
    tot = dict(checks=0, missed=0, over=0, over_noncycle=0, pairs_true=0, forced_true=0)
    for pl in range(plays):
        Q = np.zeros((D, 0), dtype=complex)
        for l in (0, 1):
            v = np.zeros(D, dtype=complex); v[l * n + 0] = 1
            Q = add_gen(g, Q, v)
        ps = PosState(g); ps.close()
        for st in range(steps):
            sT, fT, _ = true_sets(g, Q)
            sP, fP = ps.predicted()
            # compare singles
            tot["checks"] += 1
            miss_s = sT - sP; over_s = sP - sT
            # forced pairs: true pairs between non-available ends must be predicted
            miss_f = fT - fP; over_f = fP - fT
            miss_f = {(a, b) for (a, b) in miss_f if not ((0, a) in sP and (1, b) in sP)}
            if miss_s or miss_f:
                tot["missed"] += 1
                print("MISSED", name, pl, st, sorted(miss_s)[:5], sorted(miss_f)[:5], flush=True)
            if over_s or over_f:
                tot["over"] += 1
                # is every over-predicted single in a component with a cycle? (checked at position level)
                cu = getattr(ps, 'cycun', set())
                if any((l, g.pos[e]) not in cu for (l, e) in over_s) or any((0, g.pos[a]) not in cu for (a, b) in over_f):
                    tot["over_noncycle"] += 1
            tot["forced_true"] += len(fT)
            # choose a random single-term vector of V (true), biased to pairs
            choices = []
            av1 = [e for (l, e) in sT if l == 0]; av2 = [e for (l, e) in sT if l == 1]
            if fT and rng.random() < 0.6:
                e1, e2 = rng.choice(sorted(fT))
                # the forced ratio: solve from annihilator
                _, _, Nm = true_sets(g, Q)
                a_, b_ = Nm[:, e1], Nm[:, n + e2]
                c = -(np.vdot(b_, a_) / np.vdot(b_, b_))
                v = np.zeros(D, dtype=complex); v[e1] = 1; v[n + e2] = c
                kind = ("pair", e1, e2)
            elif av1 and av2 and rng.random() < 0.5:
                e1, e2 = rng.choice(av1), rng.choice(av2)
                c = complex(nrng.normal(), nrng.normal())
                v = np.zeros(D, dtype=complex); v[e1] = 1; v[n + e2] = c
                kind = ("created", e1, e2)
            else:
                l = rng.randrange(2); e = rng.choice(av1 if l == 0 else av2)
                v = np.zeros(D, dtype=complex); v[l * n + e] = 1
                kind = ("single", l, e)
            assert np.linalg.norm(annihilator(Q) @ v) < 1e-6, "chosen vector not in V"
            a = rng.choice(g.H)
            # bias upward: retry a few times for a move that raises some level
            w = leftmul(g, a, v)
            Q = add_gen(g, Q, w)
            if kind[0] == "single":
                ps.A[kind[1]].add(g.pos[g.M[a][kind[2]]])
            else:
                e1, e2 = kind[1], kind[2]
                t = g.M[g.I[e1]][e2]
                # lineage: a direct copy of an existing locked edge keeps its lineage, else a new lineage
                lin = None
                if kind[0] == "pair":
                    for (tt, xx, ll) in ps.edges:
                        if tt == t and xx == g.pos[e1]: lin = ll; break
                if lin is None: ps.nlin = getattr(ps, 'nlin', 0) + 1; lin = ps.nlin
                key = (t, g.pos[g.M[a][e1]], lin)
                if key not in ps.edges: ps.edges.append(key)
            ps.close()
    print(name, "plays", plays, "steps", steps, tot, flush=True)
    return tot

if __name__ == "__main__":
    if len(sys.argv) > 5: MODE = sys.argv[5]
    run(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
