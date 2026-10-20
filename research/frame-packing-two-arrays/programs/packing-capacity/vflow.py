#!/usr/bin/env python3
"""vflow.py (fp-capacity): exact check of the virtual live-class flow theorem (PROOF.md Theorem A) on known units.

For a word, a gauge h (one Sp(2m,2) element per input) and exact contents (sl-combinatorics' ledger engine, imported
read-only), every register's terms are grouped into FRAME CLASSES (label g h_l mod Paulis and scalars).  A class is
LIVE if, in the class-level data-flow graph, it reaches a delivering register's final class.  Each paid move on a
register with k live classes is k VIRTUAL MOVES, one per live class, each along an arc of the frame graph Gamma_m
(position lambda^-1 L0 before -> after).  Checks, all exact (integers / exact LP over rationals not needed):
  raw gauge (h = 1): directed max-flow from L0 to L1 = K L0 with arc capacity = #virtual live moves is >= W, and
      sum over virtual moves of w = max(0, Delta f)/2  >=  W m,  f = d(L0, .) - d(L1, .);
      hence  Delta = Wm - M  <=  P_w := sum_moves (sum_live w - 1)  <=  P_live := sum_moves (k_live - 1)^+.
  W = 2, any gauge: undirected two-commodity fractional flow S_a -> T_a (value 1 each) inside the virtual
      capacities exists (scipy LP; Hu's theorem says the cut condition suffices), so M + P_live(h) >= 2m.
Negative control: the same flows with ONE unit of capacity per move (packing ignored).
usage: python3 vflow.py [seed] [n_random]"""
import sys, os, random, json, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
SLC = os.path.join(HERE, "..", "ledger")
sys.path.insert(0, SLC)
from ledger import *                     # Engine, run, gauges, GQ, Fr, mul2, inv2, key2, plain_word, random_sp ...
import ledger
from cutconst import setup
from collections import deque
import numpy as np

_step = ledger._step

def lag_of(E, lam):
    """bitset over vector indices of the Lagrangian lam^-1 L0 (as in t_cutcheck.py)"""
    m = E.m; Li = inv2(lam)
    cols = [Li[:, m + j] for j in range(m)]
    out = 0
    for mask in range(1 << m):
        v = np.zeros(2 * m, dtype=np.int64)
        for j in range(m):
            if (mask >> j) & 1: v ^= cols[j]
        idx = 0
        for k in range(2 * m):
            if v[k] & 1: idx |= 1 << k
        out |= 1 << idx
    return out

def keymat(k, m): return np.frombuffer(k, dtype=np.uint8).reshape(2 * m, 2 * m).astype(np.int64)

def class_of(E, gid, l, gauge): return key2(mul2(E.syms[gid], gauge[l]))

def classes(E, c, gauge):
    out = {}
    for (l, gid) in c: out.setdefault(class_of(E, gid, l, gauge), []).append((l, gid))
    return out

def free_sym(E, name): return E.syms[E.letters[name][0]]

def virtual_moves(E, W, R, word, gauge, start=None, deliver=None):
    """replay; return list of moves, each (step index, register, [ (class key, live?) ... ]) and the unit verdict"""
    m = E.m
    if start is None: start = {r: r for r in range(W)}
    if deliver is None: deliver = dict(start)
    C = [dict() for _ in range(R)]
    for r, l in start.items(): C[r] = {(l, E.ID): ONE}
    snaps = [C]
    for st in word:
        C = _step(E, C, st); snaps.append(C)
    # unit check
    res = run(E, W, R, word, start=start, deliver=deliver)
    T = len(word)
    # backward liveness on (t, r, class)
    live = [None] * (T + 1)
    live[T] = set()
    for r, l in deliver.items():
        cl = classes(E, snaps[T][r], gauge)
        for k in cl: live[T].add((r, k))
    for t in range(T - 1, -1, -1):
        st = word[t]; cur = [classes(E, c, gauge) for c in snaps[t]]
        L = set()
        kind = st[0]
        touched = set()
        if kind in ("move", "free"):
            r = st[1]; touched.add(r)
            name = ("D%d" % st[2]) if kind == "move" else st[2]
            S = free_sym(E, name)
            for k in cur[r]:
                k2 = key2(mul2(S, keymat(k, m)))
                if (r, k2) in live[t + 1]: L.add((r, k))
        elif kind == "mix":
            r, terms = st[1], st[2]; touched.add(r)
            for (sreg, mu, name) in terms:
                S = None if name is None else free_sym(E, name)
                for k in cur[sreg]:
                    k2 = k if S is None else key2(mul2(S, keymat(k, m)))
                    if (r, k2) in live[t + 1]: L.add((sreg, k))
        elif kind == "gate":
            regs, A = st[1], st[2]
            for i, r in enumerate(regs):
                touched.add(r)
                for j, s in enumerate(regs):
                    if gq(A[i][j]).iszero(): continue
                    for k in cur[s]:
                        if (r, k) in live[t + 1]: L.add((s, k))
        for r in range(R):
            if r in touched: continue
            for k in cur[r]:
                if (r, k) in live[t + 1]: L.add((r, k))
        live[t] = L
    moves = []
    for t, st in enumerate(word):
        if st[0] != "move": continue
        r = st[1]; cur = classes(E, snaps[t][r], gauge)
        S = free_sym(E, "D%d" % st[2])
        lst = []
        for k in cur:
            k2 = key2(mul2(S, keymat(k, m)))
            lst.append((k, k2, (r, k) in live[t]))
        moves.append((t, r, lst))
    return moves, res

def maxflow(nodes, cap, s, t):
    """Edmonds-Karp on a directed multigraph given as dict (u, v) -> integer capacity"""
    adj = {}
    resid = {}
    for (u, v), c in cap.items():
        resid[(u, v)] = resid.get((u, v), 0) + c
        resid.setdefault((v, u), 0)
        adj.setdefault(u, set()).add(v); adj.setdefault(v, set()).add(u)
    flow = 0
    while True:
        par = {s: None}; dq = deque([s])
        while dq and t not in par:
            u = dq.popleft()
            for v in adj.get(u, ()):
                if v not in par and resid[(u, v)] > 0: par[v] = u; dq.append(v)
        if t not in par: return flow
        b = None; v = t
        while par[v] is not None:
            u = par[v]; b = resid[(u, v)] if b is None else min(b, resid[(u, v)]); v = u
        v = t
        while par[v] is not None:
            u = par[v]; resid[(u, v)] -= b; resid[(v, u)] += b; v = u
        flow += b

def two_commodity_feasible(cap_und, S, T):
    """LP: undirected capacities cap_und[(u,v)] (u<v); commodities a: S[a] -> T[a], value 1 each. Feasible?"""
    from scipy.optimize import linprog
    edges = sorted(cap_und)
    nodes = sorted({u for e in edges for u in e} | set(S) | set(T))
    idx = {v: i for i, v in enumerate(nodes)}
    K = len(S); E2 = len(edges)
    # variables: for each commodity a, each edge e, two directions: x[a, e, d]
    nv = K * E2 * 2
    def var(a, e, d): return (a * E2 + e) * 2 + d
    Aeq = []; beq = []
    for a in range(K):
        for v in nodes:
            row = np.zeros(nv)
            for e, (u, w) in enumerate(edges):
                # direction 0: u -> w, direction 1: w -> u
                if u == v: row[var(a, e, 0)] -= 1; row[var(a, e, 1)] += 1
                if w == v: row[var(a, e, 0)] += 1; row[var(a, e, 1)] -= 1
            rhs = 0
            if v == S[a]: rhs -= 1
            if v == T[a]: rhs += 1
            Aeq.append(row); beq.append(rhs)
    Aub = []; bub = []
    for e in range(E2):
        row = np.zeros(nv)
        for a in range(K):
            row[var(a, e, 0)] = 1; row[var(a, e, 1)] = 1
        Aub.append(row); bub.append(cap_und[edges[e]])
    if not edges:
        return all(S[a] == T[a] for a in range(K))
    r = linprog(np.zeros(nv), A_ub=np.array(Aub), b_ub=np.array(bub), A_eq=np.array(Aeq), b_eq=np.array(beq),
                bounds=[(0, None)] * nv, method="highs")
    return r.status == 0

def analyse(E, W, R, word, gauge, labs_index, Dm, start=None, deliver=None, raw=True):
    m = E.m
    moves, res = virtual_moves(E, W, R, word, gauge, start, deliver)
    if not res["unit"]: return dict(unit=False)
    pos = lambda k: labs_index[lag_of(E, keymat(k, m))]
    S = [labs_index[lag_of(E, gauge[a])] for a in range(W)]
    Tt = [labs_index[lag_of(E, mul2(E.Ksym, gauge[a]))] for a in range(W)]
    cap = {}; cap1 = {}
    P_all = 0; P_live = 0; P_w = 0.0; Mv = 0; idle_dead = 0; wsum = 0.0; asc_sum = 0; P_asc = 0
    f = None
    if raw:
        L0, L1 = S[0], Tt[0]
        f = Dm[L0] - Dm[L1]
    for (t, r, lst) in moves:
        k_all = len(lst); kl = [x for x in lst if x[2]]
        P_all += max(0, k_all - 1); P_live += max(0, len(kl) - 1); Mv += len(kl)
        if not kl: idle_dead += 1
        ww = 0.0; aa = 0
        for (k, k2, _) in kl:
            u, v = pos(k), pos(k2)
            assert Dm[u, v] == 1, "a virtual move is not an edge"
            cap[(u, v)] = cap.get((u, v), 0) + 1
            if raw:
                ww += max(0, int(f[v] - f[u])) / 2.0
                if Dm[L0, v] > Dm[L0, u]: aa += 1
        if raw: asc_sum += aa; P_asc += aa - 1
        if kl:   # negative control: one unit of capacity per move, on the first live class's arc
            k, k2, _ = sorted(kl)[0]; u, v = pos(k), pos(k2)
            cap1[(u, v)] = cap1.get((u, v), 0) + 1
        wsum += ww; P_w += ww - 1.0
    M = res["M"]; Delta = W * m - M
    out = dict(unit=True, M=M, Wm=W * m, Delta=Delta, P_all=P_all, P_live=P_live, P_w=(P_w if raw else None), Mv=Mv,
               dead_moves=idle_dead, wsum=(wsum if raw else None))
    if raw:
        fl = maxflow(None, cap, S[0], Tt[0]); fl1 = maxflow(None, cap1, S[0], Tt[0])
        out.update(flow=fl, flow_one_per_move=fl1, ok_flow=(fl >= W), ok_w=(wsum >= W * m - 1e-9),
                   ok_bound=(Delta <= P_w + 1e-9 and P_w <= P_live + 1e-9), live_ascents=asc_sum, P_asc=P_asc,
                   ok_asc=(asc_sum >= W * m and Delta <= P_asc))
    if W == 2:
        und = {}; und1 = {}
        for (u, v), c in cap.items():
            e = (min(u, v), max(u, v)); und[e] = und.get(e, 0) + c
        for (u, v), c in cap1.items():
            e = (min(u, v), max(u, v)); und1[e] = und1.get(e, 0) + c
        out.update(hu_feasible=two_commodity_feasible(und, S, Tt), hu_one_per_move=two_commodity_feasible(und1, S, Tt),
                   ok_W2=(M + P_live >= 2 * m))
    return out

def pstar_W2(E, R, word, labs_index, Dm):
    """exact min over Sp-gauges of P_live for W = 2: P_live depends only on the relative gauge
    delta = h_2 h_1^-1 (classes of input 1 and input 2 merge iff g_1 h_1 = g_2 h_2), so it suffices to try h_1 = 1
    and h_2 in {raw, one generic element} u {sym(g_2)^-1 sym(g_1) : terms (g_1, 1), (g_2, 2) of one register at any
    time}: every relative gauge outside that set merges no classes at all, exactly like the generic one."""
    m = E.m
    C = [dict() for _ in range(R)]
    for r in range(2): C[r] = {(r, E.ID): ONE}
    cands = {}
    def collect(C):
        for c in C:
            t1 = [gid for (l, gid) in c if l == 0]; t2 = [gid for (l, gid) in c if l == 1]
            for a in t1:
                for b in t2:
                    d = mul2(inv2(E.syms[b]), E.syms[a]); cands[key2(d)] = d
    collect(C)
    for st in word:
        C = _step(E, C, st); collect(C)
    I = np.eye(2 * m, dtype=np.int64)
    rng = random.Random(7)
    while True:                    # a generic relative gauge: merges no classes at any time
        gen = random_sp(E, rng, steps=20)
        if key2(gen) not in cands: break
    best = None; P_raw = None; Delta = None
    for kk, d in [("raw", I), ("generic", gen)] + list(cands.items()):
        o = analyse(E, 2, R, word, [I, d], labs_index, Dm, raw=False)
        if kk == "raw": P_raw = o["P_live"]; Delta = o["Delta"]
        if best is None or o["P_live"] < best: best = o["P_live"]
    return dict(Pstar=best, P_raw=P_raw, Delta=Delta, ncand=len(cands) + 2, ok=(Delta <= best))

def cases_for(E, rng):
    m = E.m; out = []
    for W in (1, 2, 3): out.append(("plain W=%d" % W, W, W, plain_word(m, W, rng)))
    A = [[1, 1], [1, -1]]; Ai = [[Fr(1, 2), Fr(1, 2)], [Fr(1, 2), Fr(-1, 2)]]
    out.append(("pair trick", 2, 2, [("gate", [0, 1], A)] + plain_word(m, 2, rng) + [("gate", [0, 1], Ai)]))
    w = [("mix", 2, [(0, 1, None)]), ("mix", 3, [(1, 1, None)])]
    w += [("move", 2, 1 << j) for j in range(m)] + [("move", 3, 1 << j) for j in range(m)]
    w += [("mix", 0, [(2, 1, None)]), ("mix", 1, [(3, 1, None)])]
    out.append(("copy-move-copy", 2, 4, w))
    if m >= 2:
        w = [("move", 0, 1), ("mix", 0, [(0, 1, None), (1, 1, None)])]
        for j in range(1, m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
        w += [("mix", 0, [(0, 1, None), (1, -1, None)]), ("move", 1, 1)]
        out.append(("gauge-exchange", 2, 2, w))
        w = [("mix", 2, [(1, 1, None)]), ("mix", 0, [(0, 1, None), (1, 1, None)]), ("move", 0, 1),
             ("mix", 0, [(0, 1, None), (2, -1, None)]), ("move", 0, 2), ("move", 2, 2),
             ("mix", 0, [(0, 1, None), (2, 1, None)]), ("move", 1, 1), ("move", 1, 2),
             ("mix", 0, [(0, 1, None), (1, -1, None)])]
        for j in range(2, m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
        out.append(("packing gadget", 2, 3, w))
    if m >= 3:
        s = E.addr([2, 4, 1] + [1 << j for j in range(3, m)]); s2 = E.addr([4, 1, 2] + [1 << j for j in range(3, m)])
        h = Fr(1, 2)
        w = [("mix", 0, [(0, 1, None), (1, 1, None)]), ("mix", 1, [(0, 1, None), (1, -1, None), (1, -1, s)])]
        for j in range(m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
        w += [("mix", 1, [(0, h, None), (1, -h, None), (0, -h, s), (1, h, s), (0, h, s2), (1, -h, s2)]),
              ("mix", 0, [(0, 1, None), (1, -1, None)])]
        out.append(("order-3 control C%dO3" % (2 * m), 2, 2, w))
        w = [("mix", 0, [(0, 1, None), (1, 1, None)]), ("mix", 1, [(0, 1, None), (1, -2, None)])]
        for j in range(m): w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
        w += [("mix", 1, [(0, h, None), (1, -h, None)]), ("mix", 0, [(0, 1, None), (1, -1, None)])]
        out.append(("+-1 halving control C%dH" % (2 * m), 2, 2, w))
    return out

def random_order3(E, rng):
    m = E.m
    a, b, c = rng.sample(range(m), 3)
    cols = [1 << j for j in range(m)]; cols[a], cols[b], cols[c] = 1 << b, 1 << c, 1 << a
    s = E.addr(cols)
    cols2 = [1 << j for j in range(m)]; cols2[a], cols2[b], cols2[c] = 1 << c, 1 << a, 1 << b
    s2 = E.addr(cols2)
    W = rng.choice([2, 3]); w = []
    if W == 2:
        w += [("mix", 0, [(0, 1, None), (1, 1, None)]), ("mix", 1, [(0, 1, None), (1, -1, None), (1, -1, s)])]
    else:
        w += [("mix", 0, [(0, 1, None), (1, 1, None), (2, 1, None)]),
              ("mix", 1, [(0, 1, None), (1, -1, None), (1, -1, s)]),
              ("mix", 2, [(0, 1, None), (2, -1, None), (2, -1, s)])]
    moves = []
    for r in range(W):
        order = list(range(m)); rng.shuffle(order); moves += [("move", r, 1 << j) for j in order]
    rng.shuffle(moves); w += moves
    h = Fr(1, 2)
    inv1 = lambda r: [(0, h, None), (r, -h, None), (0, -h, s), (r, h, s), (0, h, s2), (r, -h, s2)]
    if W == 2: w += [("mix", 1, inv1(1)), ("mix", 0, [(0, 1, None), (1, -1, None)])]
    else: w += [("mix", 1, inv1(1)), ("mix", 2, inv1(2)), ("mix", 0, [(0, 1, None), (1, -1, None), (2, -1, None)])]
    return ("random order-3 packer W=%d" % W, W, W, w)

def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    nrand = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    rng = random.Random(seed)
    allrows = []; bad = []
    for m in (2, 3, 4):
        E = Engine(m)
        labs, L0bits, Dm = setup(m)
        labs_index = {x: i for i, x in enumerate(labs)}
        cs = cases_for(E, rng)
        if m >= 3: cs += [random_order3(E, rng) for _ in range(nrand)]
        for (name, W, R, w) in cs:
            gauges = [("raw", gauge_trivial(E, W), True)]
            if W == 2:
                gc, why = find_coherent_gauge(E, W, R, w)
                if gc is not None: gauges.append(("coherent", gc, False))
                gauges += [("random%d" % i, [random_sp(E, rng) for _ in range(W)], False) for i in range(2)]
            for (gname, g, raw) in gauges:
                o = analyse(E, W, R, w, g, labs_index, Dm, raw=raw)
                o.update(case=name, m=m, W=W, gauge=gname)
                allrows.append(o)
                oks = [o.get(k, True) for k in ("ok_flow", "ok_w", "ok_bound", "ok_asc", "ok_W2", "hu_feasible")]
                if not o["unit"] or not all(oks): bad.append(o)
                print("m=%d %-28s %-9s M=%s Wm=%s Delta=%s P_all=%s P_live=%s P_w=%s P_asc=%s Mv=%s flow=%s flow1=%s hu=%s hu1=%s"
                      % (m, name, gname, o.get("M"), o.get("Wm"), o.get("Delta"), o.get("P_all"), o.get("P_live"),
                         o.get("P_w"), o.get("P_asc"), o.get("Mv"), o.get("flow"), o.get("flow_one_per_move"),
                         o.get("hu_feasible"), o.get("hu_one_per_move")), flush=True)
            if W == 2:
                ps = pstar_W2(E, R, w, labs_index, Dm)
                o2 = dict(case=name, m=m, W=W, gauge="P*-search", **ps); allrows.append(o2)
                if not ps["ok"]: bad.append(o2)
                print("m=%d %-28s P* (min over all relative gauges) = %s over %d candidate gauges (raw %s); Delta=%s"
                      % (m, name, ps["Pstar"], ps["ncand"], ps["P_raw"], ps["Delta"]), flush=True)
    neg = [(o["case"], o["m"], o["gauge"]) for o in allrows if o.get("flow_one_per_move") is not None and o["flow_one_per_move"] < o["W"]]
    neg2 = [(o["case"], o["m"], o["gauge"]) for o in allrows if o.get("hu_one_per_move") is False]
    print("\nSUMMARY: %d rows; violations of Theorem A: %d; negative control (one capacity per move) fails max-flow in %d raw rows: %s; fails two-commodity LP in %d rows: %s"
          % (len(allrows), len(bad), len(neg), sorted(set(neg)), len(neg2), sorted(set(neg2))))
    with open(os.path.join(HERE, "vflow_s%d.json" % seed), "w") as fh:
        json.dump(dict(rows=allrows, bad=bad), fh, indent=1, default=str)

if __name__ == "__main__":
    main()
