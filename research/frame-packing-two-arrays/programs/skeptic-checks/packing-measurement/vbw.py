#!/usr/bin/env python3
"""v-fp-measure: independent replay of the BACKWARD call (x -= y) of a gcert/1 word, read as the transposed inverse
of the forward call F = g8 . W' . g1 (W' = A, scatter, B, final climbs; g1 = "y -= K s" at the start, g8 = "s :=
Fss^-1 (s - J x)" at the end), gate by gate in the same time order and at the same frames.  Written from the
statement in fp-measure NOTE.md / bw_ledger.py's docstring only (bw_ledger.py and fm_ledger.py not imported).
Frame model as in vfm.py.  Clean helpers (slots start empty).
usage: python3 -I vbw.py WORD.gcert.json.gz [B1|B2]"""
import sys, json, gzip
from fractions import Fraction as Q
from collections import Counter


def rank(vecs):
    rows = []
    for x in vecs:
        for b in rows:
            x = min(x, x ^ b)
        if x:
            rows.append(x)
            rows.sort(reverse=True)
    return len(rows)


def add(d, src, q):
    for key, u in src.items():
        w = d.get(key, Q(0)) + q * u
        if w:
            d[key] = w
        else:
            d.pop(key, None)


def scat_rows(c):
    sc, h, ports = c["scat"], c["h"], c["ports"]
    if "table" in sc:
        return sc["table"]
    return [[(k,) + tuple(sc["inside"] if ports[t] >> k & 1 else sc["outside"]) for k in range(h)]
            for t in range(c["v"])]


def forward_symbolic(c):
    v, R = c["v"], c["R"]
    n = 2 * v + R
    reg = [dict() for _ in range(n)]
    for t in range(v):
        reg[t] = {("x", t): Q(1)}
        reg[v + t] = {("y", t): Q(1)}
    for r in range(2 * v, n):
        reg[r] = {("d", r): Q(1)}
    def run(gs):
        for g in gs:
            if g[0] == "out":
                src = dict(reg[g[2]])
                for t, a, b in g[3]:
                    add(reg[t], src, Q(a, b))
            else:
                for s, a, b in g[3]:
                    add(reg[g[2]], dict(reg[s]), Q(a, b))
    run(c["A"])
    ret = sorted(c["ret"])
    tot = [dict(reg[s]) for _, s, _ in ret]
    for t, row in enumerate(scat_rows(c)):
        for k, a, b in row:
            add(reg[v + t], tot[k], Q(a, b))
    run(c["B"])
    for t in range(v):
        assert reg[t] == {("x", t): Q(1)}
    K = [{} for _ in range(v)]
    for t in range(v):
        rest = dict(reg[v + t])
        assert rest.pop(("y", t)) == 1 and rest.pop(("x", t)) == 1
        for key, u in rest.items():
            assert key[0] == "d", key
            K[t][key[1]] = u
    J, F = {}, {}
    for r in range(2 * v, n):
        for key, u in reg[r].items():
            if key[0] == "x":
                J.setdefault(r, {})[key[1]] = u
            elif key[0] == "d":
                F.setdefault(r, {})[key[1]] = u
            else:
                raise AssertionError(("y in slot", r, key))
    return K, J, F


def main():
    path = sys.argv[1]
    init = sys.argv[2] if len(sys.argv) > 2 else "B1"
    c = json.load(gzip.open(path, "rt"))
    h, v, R = c["h"], c["v"], c["R"]
    frames = c["frames"]
    dim = [rank(f) for f in frames]
    ports = c["ports"]
    n = 2 * v + R
    ret = sorted(c["ret"])
    nc = len(ret)
    K, J, F = forward_symbolic(c)
    where = list(c["start"]) + [None] * nc
    cont = [dict() for _ in range(n + nc)]
    for t in range(v):
        if init == "B1":
            cont[t] = {("x", t): Q(1)}; cont[v + t] = {("y", t): Q(1)}
        elif init == "B2":
            cont[t] = {("x", t): Q(1)}; cont[v + t] = {("y", t): Q(1), ("x", t): Q(1)}
    x0 = [dict(cont[t]) for t in range(v)]
    y0 = [dict(cont[v + t]) for t in range(v)]
    def twist(lab):
        return ports[lab[1]] if lab[0] == "x" else 0
    capc = {}
    def cap(f, s):
        if s == 0:
            return 0
        if (f, s) not in capc:
            capc[f, s] = 1 if rank(frames[f] + [s]) == dim[f] else 0
        return capc[f, s]
    def lev(lab, f):
        s = twist(lab)
        return dim[f] + (1 if s else 0) - 2 * cap(f, s)
    def L(r):
        return {("Y",) if k[0] == "y" else ("x", k[1]) for k in cont[r]}
    def psi(r):
        return sum(lev(l, where[r]) for l in L(r)) if where[r] is not None else 0
    S = Counter(); kP = Counter(); kM = Counter(); kD = Counter(); kI = Counter()
    mx = [0]
    def kind(r):
        return "x" if r < v else "y" if r < 2 * v else "s" if r < n else "c"
    def climb(r, f):
        o = where[r]
        if o == f:
            return
        assert dim[f] > dim[o] and rank(frames[f] + frames[o]) == dim[f], "not nested"
        k = dim[f] - dim[o]
        labs = L(r)
        S["M"] += k; kM[kind(r)] += k
        if not labs:
            S["I0"] += k; kI[kind(r)] += k
        else:
            S["P"] += k * (len(labs) - 1); kP[kind(r)] += k * (len(labs) - 1)
            mx[0] = max(mx[0], len(labs))
            for l in labs:
                d = cap(f, twist(l)) - cap(o, twist(l))
                S["D"] += d; kD[kind(r)] += d
        where[r] = f
    def mix(targets, op):
        before = sum(psi(r) for r in targets)
        op()
        S["Bmix"] -= sum(psi(r) for r in targets) - before
    psi0 = sum(psi(r) for r in range(n))
    # g1^T: s += K^T y   (at the start frames)
    tg = sorted({r for t in range(v) for r in K[t]})
    for r in tg:
        assert where[r] == 0, ("g1 slot not at frame 0", r, where[r])
    for t in range(v):
        if K[t]:
            assert where[v + t] == 0
    def g1T():
        snap = [dict(cont[v + t]) for t in range(v)]
        for t in range(v):
            for r, u in K[t].items():
                add(cont[r], snap[t], u)
    mix(tg, g1T)
    def gate(g):
        f = g[1]
        rs = [g[2]] + [x[0] for x in g[3]] + (list(g[4]) if g[0] == "out" else [])
        for r in rs:
            climb(r, f)
        if g[0] == "out":          # forward: t += k src for each t  ->  src -= sum k t   (old t values)
            s0 = g[2]
            def op():
                snap = {t: dict(cont[t]) for t, _, _ in g[3]}
                for t, a, b in g[3]:
                    add(cont[s0], snap[t], -Q(a, b))
            mix([s0], op)
        else:                      # forward: t += sum k_s s  ->  each s -= k_s t   (old t value)
            t0 = g[2]
            def op():
                snap = dict(cont[t0])
                for s, a, b in g[3]:
                    add(cont[s], snap, -Q(a, b))
            mix([s for s, _, _ in g[3]], op)
    for g in c["A"]:
        gate(g)
    rows = scat_rows(c)
    for k, (_, s, f) in enumerate(ret):
        assert where[s] == f
        cr = n + k
        where[cr] = 0
        def make(cr=cr, k=k):
            for t in range(v):
                assert where[v + t] == 0
                for kk, a, b in rows[t]:
                    if kk == k:
                        add(cont[cr], dict(cont[v + t]), -Q(a, b))
        mix([cr], make)
        climb(cr, f)
        def into(cr=cr, s=s):
            add(cont[s], dict(cont[cr]), Q(1))
            cont[cr] = {}
        mix([s, cr], into)
    for g in c["B"]:
        gate(g)
    for r in range(n):
        climb(r, c["final"][r])
    # g8^T: x += J^T s, s := Fss^T s  (at the end frames)
    def g8T():
        snap = {r: dict(cont[r]) for r in range(2 * v, n)}
        for r, row in J.items():
            for t, u in row.items():
                add(cont[t], snap[r], u)
        for r in range(2 * v, n):
            cont[r] = {}
        for r, row in F.items():          # forward s'_r = sum_q F[r][q] s_q ; transpose: s'_q += F[r][q] s_r
            for q, u in row.items():
                add(cont[q], snap[r], u)
    mix(list(range(n)), g8T)
    psiE = sum(psi(r) for r in range(n + nc))
    def minus(a, b):
        o = dict(a)
        for kk, u in b.items():
            o[kk] = o.get(kk, 0) - u
        return {kk: u for kk, u in o.items() if u}
    okx = all(cont[t] == minus(x0[t], y0[t]) for t in range(v))
    oky = all(cont[v + t] == y0[t] for t in range(v))
    oks = all(not cont[r] for r in range(2 * v, n + nc))
    deliv = sum(psi(r) for r in range(2 * v))
    ident = psiE - psi0 + S["Bmix"] + 2 * S["D"] - S["P"] + S["I0"] == S["M"]
    B = S["Bmix"] + psiE - deliv
    print(json.dumps(dict(word=path.split("/")[-1], init=init, M=S["M"], P=S["P"], B=str(B), D=S["D"], I0=S["I0"],
                          maxlab=mx[0], psi0=psi0, deliv=deliv, x_ok=okx, y_ok=oky, s_restored=oks, identity=ident,
                          balance=str(S["P"] - B - 2 * S["D"] - S["I0"]), P_per_move=S["P"] / S["M"],
                          P_by_kind=dict(kP), M_by_kind=dict(kM), D_by_kind=dict(kD), I0_by_kind=dict(kI),
                          nK=sum(len(x) for x in K), nJ=sum(len(x) for x in J.values()))))


if __name__ == "__main__":
    main()
