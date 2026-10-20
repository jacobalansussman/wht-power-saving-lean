# w3cut.py -- session 2. Undirected cut LP for a FIXED layout of opposite pairs on Gamma_m (sl-combinatorics
# Lemma 4.2 form: moves crossing Q in either direction >= pairs split by Q), by constraint generation with an exact
# rational certificate extracted from the dual (a cut packing), plus an analysis of which cuts the packing uses, and a
# RESTRICTED-family packing LP (balls, residue cuts F_p, pulled-back balls) whose value is a valid lower bound too.
#
#   cut value(D) = min sum_e c_e  s.t.  c(delta(Q)) >= sep_D(Q) for all Q
#                = max sum_Q lam_Q sep_D(Q)  s.t.  lam >= 0, sum_{Q : e in delta(Q)} lam_Q <= 1 for every edge e.
# Any lam feasible for the right-hand LP (over any family of Q) is a lower bound on the moves of every
# coherent unit with layout D. The certificate is checked in exact rationals (fractions.Fraction).
import sys, random, time, json, itertools
from fractions import Fraction
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
from scipy.sparse.csgraph import maximum_flow, breadth_first_order
import dpg

SCALE = 10 ** 6

class G:
    def __init__(self, m):
        self.m = m
        verts, idx, adj, dist = dpg.build(m)
        self.verts, self.idx, self.adj = verts, idx, adj
        self.n = len(verts)
        self.D = np.array(dist, dtype=np.int16)
        self.edges = [(u, v) for u in range(self.n) for v in adj[u] if u < v]
        self.eu = np.array([u for u, v in self.edges]); self.ev = np.array([v for u, v in self.edges])
        self.opp = [np.nonzero(self.D[s] == m)[0].tolist() for s in range(self.n)]

def sep_of(inQ, pairs):
    return sum(1 for (s, t) in pairs if inQ[s] != inQ[t])

def full_cut_lp(g, pairs, max_rounds=400, verbose=False):
    m = g.m; n = g.n
    X = sorted({v for p in pairs for v in p})
    pats = []
    for mask in range(0, 1 << (len(X) - 1)):
        A = {X[0]} | {X[i + 1] for i in range(len(X) - 1) if mask >> i & 1}
        if len(A) == len(X): continue
        s_ = sum(1 for (s, t) in pairs if (s in A) != (t in A))
        if s_ > 0: pats.append((sorted(A), [x for x in X if x not in A], s_))
    rows = []; rhs = []; Qs = []
    def addQ(inQ):
        rows.append((inQ[g.eu] != inQ[g.ev]).astype(float)); rhs.append(sep_of(inQ, pairs)); Qs.append(inQ.copy())
    for x in X:
        for k in range(m):
            addQ(g.D[x] <= k)
    E = len(g.edges); S_, T_ = n, n + 1
    for rnd in range(max_rounds):
        A = sp.csr_matrix(np.array(rows))
        res = linprog(np.ones(E), A_ub=-A, b_ub=-np.array(rhs, dtype=float), bounds=(0, None), method='highs')
        c = res.x; val = res.fun
        cap = np.floor(c * SCALE).astype(np.int64)
        added = 0
        for (Aside, Bside, s_) in pats:
            big = SCALE * 1000
            r_ = np.concatenate([g.eu, g.ev, np.full(len(Aside), S_), np.array(Bside)])
            c_ = np.concatenate([g.ev, g.eu, np.array(Aside), np.full(len(Bside), T_)])
            w_ = np.concatenate([cap, cap, np.full(len(Aside), big), np.full(len(Bside), big)])
            Cap = sp.csr_matrix((w_.astype(np.int32), (r_, c_)), shape=(n + 2, n + 2)); Cap.sum_duplicates()
            fl = maximum_flow(Cap, S_, T_)
            if fl.flow_value < s_ * SCALE - 1:
                F = fl.flow.tocsr()
                R = (Cap - F).tocsr(); R.data = (R.data > 0).astype(np.int8); R.eliminate_zeros()
                reach = breadth_first_order(R, S_, directed=True, return_predecessors=False)
                inQ = np.zeros(n, dtype=bool); inQ[[v for v in reach.tolist() if v < n]] = True
                row = (inQ[g.eu] != inQ[g.ev]).astype(float)
                if row @ c < sep_of(inQ, pairs) - 1e-7:
                    addQ(inQ); added += 1
        if verbose: print(f"   round {rnd}: LP {val:.6f}, cuts {len(rows)}, added {added}", flush=True)
        if added == 0:
            lam = -np.array(res.ineqlin.marginals)
            return val, lam, Qs, rhs, rnd
    return None, None, Qs, rhs, max_rounds

def exact_packing(g, lam, Qs, pairs, den=10 ** 6):
    """round lam to rationals, scale so every edge load <= 1, return exact value (Fraction) and the support."""
    sup = [i for i in range(len(lam)) if lam[i] > 1e-9]
    fr = {i: Fraction(int(round(lam[i] * den)), den) for i in sup}
    fr = {i: f for i, f in fr.items() if f > 0}
    load = [Fraction(0)] * len(g.edges)
    for i, f in fr.items():
        cut = np.nonzero(Qs[i][g.eu] != Qs[i][g.ev])[0]
        for e in cut.tolist(): load[e] += f
    mx = max(load) if load else Fraction(1)
    scale = Fraction(1) if mx <= 1 else 1 / mx
    val = sum(f * sep_of(Qs[i], pairs) for i, f in fr.items()) * scale
    return val, {i: f * scale for i, f in fr.items()}

def classify(g, inQ):
    """name the cut if it is a ball B(A,k) (or its complement), a residue F_p, or a pulled-back ball."""
    m = g.m
    for side, Q in (('', inQ), ('~', ~inQ)):
        sz = int(Q.sum())
        # ball: Q = {d(., A) <= k}
        for k in range(m):
            cand = np.nonzero(Q)[0]
            for A in cand.tolist():
                if np.array_equal(g.D[A] <= k, Q): return f"{side}ball(k={k})"
        # residue F_p
    return 'other'

def families(g, which=('ball', 'res', 'pull')):
    """candidate cuts as boolean vectors."""
    m, n = g.m, g.n
    out = []; names = []
    if 'ball' in which:
        for A in range(n):
            for k in range(m):
                out.append(g.D[A] <= k); names.append(('ball', A, k))
    vecs = range(1, 1 << (2 * m))
    mem = None
    if 'res' in which or 'pull' in which:
        mem = np.zeros((1 << (2 * m), n), dtype=bool)
        for j, L in enumerate(g.verts):
            for x in L: mem[x, j] = True
    if 'res' in which:
        for p in vecs:
            out.append(mem[p].copy()); names.append(('res', p))
    if 'pull' in which:
        # f(L) = d(L, M) - [p not in L] for Lagrangians M containing p: the pull-back of d(., M/p) along the
        # projection from p (1-Lipschitz on Gamma_m); level sets {f <= k}, k = 0..m-2
        for p in vecs:
            Ms = np.nonzero(mem[p])[0].tolist()
            notp = ~mem[p]
            for M in Ms:
                f = g.D[M].astype(int) - notp.astype(int)
                for k in range(0, m - 1):
                    out.append(f <= k); names.append(('pull', p, M, k))
    return out, names

def restricted_lp(g, pairs, fam):
    Q = np.array(fam)                     # F x n
    cutmat = (Q[:, g.eu] != Q[:, g.ev])   # F x E
    sepv = np.array([sep_of(q, pairs) for q in fam], dtype=float)
    keep = np.nonzero(sepv > 0)[0]
    A = sp.csr_matrix(cutmat[keep].T.astype(float))   # E x F'
    res = linprog(-sepv[keep], A_ub=A, b_ub=np.ones(len(g.edges)), bounds=(0, None), method='highs')
    return -res.fun, keep, res.x

def rand_3k2(g, rng, W=3):
    while True:
        P = []; used = set()
        for _ in range(W):
            for _t in range(100):
                s = rng.randrange(g.n); t = rng.choice(g.opp[s])
                if s not in used and t not in used: break
            P.append((s, t)); used |= {s, t}
        if len(used) == 2 * W: return P

if __name__ == '__main__':
    m = int(sys.argv[1]); trials = int(sys.argv[2]); seed = int(sys.argv[3]); W = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    mode = sys.argv[5] if len(sys.argv) > 5 else 'both'
    t0 = time.time(); g = G(m); rng = random.Random(seed)
    print(f"Gamma_{m}: {g.n} frames, {len(g.edges)} edges ({time.time()-t0:.0f}s)", flush=True)
    fam = names = None
    if mode in ('both', 'restricted'):
        fam, names = families(g); print(f"family: {len(fam)} cuts", flush=True)
    out = []
    for tr in range(trials):
        pairs = rand_3k2(g, rng, W)
        rec = dict(pairs=pairs)
        if mode in ('both', 'restricted'):
            v, keep, x = restricted_lp(g, pairs, fam)
            used = [names[keep[i]][0] for i in range(len(x)) if x[i] > 1e-9]
            rec['restricted'] = v; rec['restricted_use'] = {k: used.count(k) for k in set(used)}
        if mode in ('both', 'full'):
            val, lam, Qs, rhs, rnd = full_cut_lp(g, pairs)
            ev, fr = exact_packing(g, lam, Qs, pairs)
            kinds = [classify(g, Qs[i]) for i in fr]
            rec.update(full=val, exact=str(ev), exact_float=float(ev), rounds=rnd, support=len(fr),
                       kinds={k: kinds.count(k) for k in set(kinds)})
        out.append(rec)
        print(f"trial {tr}: {pairs} " + " ".join(f"{k}={rec[k]}" for k in rec if k != 'pairs') + f" ({time.time()-t0:.0f}s)", flush=True)
    json.dump(out, open(f"w3cut_m{m}_W{W}_s{seed}_{mode}.json", "w"), default=str)
