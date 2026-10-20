#!/usr/bin/env python3
"""v-fp-measure: independent replay of one forward call of a gcert/1 word in the frame model stated by fp-measure
(NOTE.md section 3), written from gxcore.py's semantics only (fm_ledger.py not imported, not copied).

Model (as claimed):
- register at subspace frame U (a subspace of F_2^h); input with twist S sits at raw level dim U + dim S - 2 dim(U^S);
- x_t has twist <p_t>; y_t has twist 0; all y share ONE label (same twist); distinct x_t distinct labels;
- a climb U -> U' (nested) is dim U' - dim U paid moves; label S descends dim(U'^S) - dim(U^S) times;
- centre: each retained total is copied at its frame f, the copy descends f -> 0 in dim f moves, then the scatter
  rows add the copies into the y at frame 0 (copies left as scratch).
Ledger: P = sum over moves of (#labels - 1) on non-empty registers; I0 = moves on empty registers; D descents;
Bmix = -(sum of Psi changes at gates); the segment identity M = Psi_end - Psi0 + Bmix + 2D + Z - P + I0 is checked.
usage: python3 -I vpos.py WORD.gcert.json.gz [F1|F2|F3]"""
import sys, json, gzip
from fractions import Fraction as Q
from collections import Counter


def rank(vecs):
    rows, r = [], 0
    for x in vecs:
        for b in rows:
            x = min(x, x ^ b)
        if x:
            rows.append(x)
            rows.sort(reverse=True)
    return len(rows)


def main():
    path = sys.argv[1]
    init = sys.argv[2] if len(sys.argv) > 2 else "F1"
    c = json.load(gzip.open(path, "rt"))
    h, v, R = c["h"], c["v"], c["R"]
    frames = c["frames"]
    dim = [rank(f) for f in frames]
    assert all(dim[i] == len(f) for i, f in enumerate(frames))
    ports = c["ports"]
    n = 2 * v + R
    ret = sorted(c["ret"])
    nc = len(ret)
    where = list(c["start"]) + [None] * nc
    # twist of a label: 0 for y, port for x_t; label id: ('Y',) or ('x', t)
    def twist(lab):
        return ports[lab[1]] if lab[0] == "x" else 0
    capc = {}
    def cap(f, s):           # dim(U ^ <s>) for a line s (or 0)
        if s == 0:
            return 0
        k = (f, s)
        if k not in capc:
            capc[k] = 1 if rank(frames[f] + [s]) == dim[f] else 0
        return capc[k]
    posc = {}
    def rref(vs):
        rows = []
        for x in vs:
            for b in rows:
                x = min(x, x ^ b)
            if x:
                rows = [min(r, r ^ x) if (r >> (x.bit_length() - 1)) & 1 else r for r in rows]
                rows.append(x)
                rows.sort(reverse=True)
        return tuple(rows)
    def position(lab, f):
        key = (lab if lab[0] == "x" else ("Y",), f)
        if key not in posc:
            U = frames[f]
            perp = [a for a in range(1, 1 << h) if all(bin(a & u).count("1") % 2 == 0 for u in U)]
            base = [(u << h) | u for u in U] + [a for a in rref(perp)]
            if lab[0] == "x":
                p = ports[lab[1]]
                base = [w ^ ((p << h) if bin(p & (w & ((1 << h) - 1))).count("1") & 1 else 0) for w in base]
            r = rref(base)
            assert len(r) == h, "not Lagrangian"
            posc[key] = r
        return posc[key]
    def lev(lab, f):
        s = twist(lab)
        return dim[f] + (1 if s else 0) - 2 * cap(f, s)
    cont = [dict() for _ in range(n + nc)]
    def L(r):
        return {("Y",) if k[0] == "y" else ("x", k[1]) for k in cont[r]}
    for t in range(v):
        if init == "F1":
            cont[t] = {("x", t): Q(1)}; cont[v + t] = {("y", t): Q(1)}
        elif init == "F2":
            cont[t] = {("x", t): Q(1), ("y", t): Q(-1)}; cont[v + t] = {("y", t): Q(1)}
        elif init == "F3":
            cont[t] = {("y", t): Q(-1)}; cont[v + t] = {("y", t): Q(1), ("x", t): Q(1)}
    x0 = [dict(cont[t]) for t in range(v)]
    y0 = [dict(cont[v + t]) for t in range(v)]
    def psi(r):
        return sum(lev(l, where[r]) for l in L(r)) if where[r] is not None else 0
    S = Counter()
    kindP = Counter(); kindM = Counter()
    maxl = [0]
    def kind(r):
        return "x" if r < v else "y" if r < 2 * v else "s" if r < n else "c"
    def climb(r, f, down=False):
        o = where[r]
        if o == f:
            return
        if not down:
            assert dim[f] > dim[o] and rank(frames[f] + frames[o]) == dim[f], "not nested"
            k = dim[f] - dim[o]
        else:
            assert f == 0
            k = dim[o]
        labs = L(r)
        S["M"] += k; kindM[kind(r)] += k
        if not labs:
            S["I0"] += k
        else:
            S["P"] += k * (len(labs) - 1); kindP[kind(r)] += k * (len(labs) - 1)
            S["Ppos"] += k * (len({position(l, o) for l in labs}) - 1)
            maxl[0] = max(maxl[0], len(labs))
            for l in labs:
                s = twist(l)
                if not down:
                    d = cap(f, s) - cap(o, s)
                else:
                    # level goes lev(o) -> dim S over k moves; ascents - descents = dim S - lev(o)
                    ds = 1 if s else 0
                    a_minus_d = ds - lev(l, o)
                    d = (k - a_minus_d) // 2
                    assert (k - a_minus_d) % 2 == 0
                assert 0 <= d <= k
                S["D"] += d
        where[r] = f
    def add(t, src, q):
        d = cont[t]
        for key, u in src.items():
            w = d.get(key, Q(0)) + q * u
            if w:
                d[key] = w
            else:
                d.pop(key, None)
    def gate_mix(targets, ops):
        before = sum(psi(r) for r in targets)
        ops()
        S["Bmix"] -= sum(psi(r) for r in targets) - before
    psi0 = sum(psi(r) for r in range(n))
    def do(g):
        f = g[1]
        rs = [g[2]] + [x[0] for x in g[3]] + (list(g[4]) if g[0] == "out" else [])
        for r in rs:
            climb(r, f)
        if g[0] == "out":
            src = dict(cont[g[2]])
            gate_mix([x[0] for x in g[3]], lambda: [add(t, src, Q(a, b)) for t, a, b in g[3]])
        else:
            def op():
                for s, a, b in g[3]:
                    add(g[2], dict(cont[s]), Q(a, b))
            gate_mix([g[2]], op)
    for g in c["A"]:
        do(g)
    cr = []
    for k, (_, s, f) in enumerate(ret):
        assert where[s] == f
        r = n + k
        where[r] = f
        gate_mix([r], lambda: add(r, dict(cont[s]), Q(1)))
        climb(r, 0, down=True)
        cr.append(r)
    sc = c["scat"]
    if "table" in sc:
        rows = sc["table"]
    else:
        rows = [[(k,) + tuple(sc["inside"] if ports[t] >> k & 1 else sc["outside"]) for k in range(h)] for t in range(v)]
    for t in range(v):
        assert where[v + t] == 0
        def op(t=t):
            for k, a, b in rows[t]:
                add(v + t, dict(cont[cr[k]]), Q(a, b))
        gate_mix([v + t], op)
    for g in c["B"]:
        do(g)
    for r in range(n):
        climb(r, c["final"][r])
    psiE = sum(psi(r) for r in range(n + nc))
    ok_x = all(cont[t] == x0[t] for t in range(v))
    def plus(a, b):
        o = dict(a)
        for k2, u in b.items():
            o[k2] = o.get(k2, 0) + u
        return {k2: u for k2, u in o.items() if u}
    ok_y = all(cont[v + t] == plus(y0[t], x0[t]) for t in range(v))
    deliv = sum(psi(r) for r in range(2 * v))
    ident = (psiE - psi0 + S["Bmix"] + 2 * S["D"] - S["P"] + S["I0"]) == S["M"]
    B = S["Bmix"] + psiE - deliv
    out = dict(word=path.split("/")[-1], init=init, M=S["M"], N=c["N"], P=S["P"], D=S["D"], I0=S["I0"],
               Bmix=str(S["Bmix"]), B=str(B), psi0=psi0, psiE=psiE, deliv=deliv, maxlab=maxl[0],
               balance=str(S["P"] - B - 2 * S["D"] - S["I0"]), deliv_minus_psi0_minus_M=deliv - psi0 - S["M"],
               x_ok=ok_x, y_ok=ok_y, identity=ident, P_per_move=S["P"] / S["M"], BoverP=float(B) / S["P"] if S["P"] else None,
               P_by_kind=dict(kindP), Ppos=S["Ppos"], Ppos_per_move=S["Ppos"] / S["M"], M_by_kind=dict(kindM))
    print(json.dumps(out))


if __name__ == "__main__":
    main()
