#!/usr/bin/env python3
"""fp-measure: the raw-gauge ledger of the BACKWARD call (x -= y) of a gcert/1 local word, in the frame model.

The backward reading (as an-build read it in our Lean: GCert/Chain/Phased.lean, Carrier/Live.lean) is the transposed
inverse of every gate of the forward call, applied in the same time order, with the same frames, including the two
dirt gates of the forward call:
  g1 = "y -= K s" at the start frames, g8 = "slots := rows of the inverse" (s := Fss^-1 (s - J x)) at the end frames.
Here the forward word W' (A, scatter, B, final climbs) maps (x, y, s) to (x, y + x + K s, J x + Fss s); K, J, Fss are
read off an exact forward replay with symbolic dirt in every slot.  Then F = g8 W' g1 = (x, y + x, s) and the
backward word is (F^-1)^T = (x - y, y, s), which this program checks exactly.
Transposed inverses: "t += k src" -> "src -= k t";  g1 -> "s += K^T y";  g8 -> "x += J^T s, s := Fss^T s";
the scatter "y_t += c_tk s_k" (a total read across frames, charged dim f_k) -> a centre register created at frame 0
holding -sum_t c_tk y_t, climbing 0 -> f_k (the same dim f_k moves), added into s_k at f_k and consumed.
(In the forward ledger, fm_ledger.py, the centre is a copy of s_k that descends f_k -> 0.)
usage: python3 -I bw_ledger.py WORD.gcert.json.gz NAME [init ...]   init: B1 (x, y fresh: sources y carry one twist),
       B2 (x, y + x: the second shear of a swap)."""
import sys, json, gzip
from fractions import Fraction as Fr
from collections import Counter, defaultdict

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from fm_ledger import load, span_has


_FM = {}


def forward_matrices(D):
    if id(D) in _FM:
        return _FM[id(D)]
    _FM[id(D)] = _forward_matrices(D)
    return _FM[id(D)]


def _forward_matrices(D):
    """exact forward replay of W' with symbolic x, y and dirt in every slot: returns K[t] (dict slot -> coef),
    J[r] (dict t -> coef), Fss[r] (dict slot -> coef)."""
    v, R = D["v"], D["R"]
    n = 2 * v + R
    reg = [dict() for _ in range(n)]
    for t in range(v):
        reg[t] = {("x", t): Fr(1)}
        reg[v + t] = {("y", t): Fr(1)}
    for r in range(2 * v, n):
        reg[r] = {("d", r): Fr(1)}
    def axpy(t, src, k):
        d = reg[t]
        for kk, u in src.items():
            w = d.get(kk, Fr(0)) + k * u
            if w:
                d[kk] = w
            else:
                d.pop(kk, None)
    def run(gs):
        for g in gs:
            if g[0] == "out":
                src = dict(reg[g[2]])
                for t, a, b in g[3]:
                    axpy(t, src, Fr(a, b))
            else:
                for s, a, b in g[3]:
                    axpy(g[2], dict(reg[s]), Fr(a, b))
    run(D["A"])
    tot = {k: dict(reg[s]) for k, s, f in D["ret"]}
    for t, row in enumerate(scat(D)):
        for k, a, b in row:
            axpy(v + t, tot[k], Fr(a, b))
    run(D["B"])
    K = [{kk[1]: c for kk, c in reg[v + t].items() if kk[0] == "d"} for t in range(v)]
    for t in range(v):
        rest = {kk: c for kk, c in reg[v + t].items() if kk[0] != "d"}
        assert rest == {("x", t): 1, ("y", t): 1}, "forward y_t is not y_t + x_t + K s"
        assert reg[t] == {("x", t): 1}
    J = {r: {kk[1]: c for kk, c in reg[r].items() if kk[0] == "x"} for r in range(2 * v, n)}
    Fss = {r: {kk[1]: c for kk, c in reg[r].items() if kk[0] == "d"} for r in range(2 * v, n)}
    assert all(kk[0] in "xd" for r in range(2 * v, n) for kk in reg[r]), "y reached a slot"
    return K, J, Fss


def scat(D):
    v, h, ports = D["v"], D["h"], D["ports"]
    if "table" in D["scat"]:
        return D["scat"]["table"]
    ins, outs = D["scat"]["inside"], D["scat"]["outside"]
    return [[(k,) + tuple(ins if ports[t] >> k & 1 else outs) for k in range(h)] for t in range(v)]


def run_bw(D, init="B1", dirt=False, dirt_one=False):
    h, v, R = D["h"], D["v"], D["R"]
    fr, ports = D["frames"], D["ports"]
    dim = [len(f) for f in fr]
    n = 2 * v + R
    rets = sorted(D["ret"])
    nreg = n + len(rets)
    K, J, Fss = forward_matrices(D)
    cur = list(D["start"]) + [None] * len(rets)
    final = list(D["final"])
    reg = [dict() for _ in range(nreg)]
    for t in range(v):
        if init == "B1":
            reg[t] = {("x", t): Fr(1)}
            reg[v + t] = {("y", t): Fr(1)}
        elif init == "B2":
            reg[t] = {("x", t): Fr(1)}
            reg[v + t] = {("y", t): Fr(1), ("x", t): Fr(1)}
        else:
            raise ValueError(init)
    if dirt:
        for r in range(2 * v, n):
            reg[r] = {("d", r): Fr(1)}
    x0 = [dict(reg[t]) for t in range(v)]
    y0 = [dict(reg[v + t]) for t in range(v)]
    s0 = [dict(reg[r]) for r in range(2 * v, n)]
    def cls(key):
        k, i = key
        return ("x", i) if k == "x" else ("Y",) if k == "y" else (("D",) if dirt_one else ("d", i))
    memb = {}
    def inU(p, f):
        if (p, f) not in memb:
            memb[(p, f)] = span_has(fr[f], p)
        return memb[(p, f)]
    def level(c, f):
        if c[0] == "x":
            return dim[f] + 1 - 2 * (1 if inU(ports[c[1]], f) else 0)
        return dim[f]
    def classes(r):
        return set(cls(k) for k in reg[r])
    def psi(r):
        return sum(level(c, cur[r]) for c in classes(r)) if reg[r] else 0
    L = dict(M=0, A=0, D=0, Z=0, P=0, I0=0, Bmix=Fr(0))
    by = defaultdict(Counter)
    hist = Counter()
    maxlab = [0]
    def kind(r):
        return "x" if r < v else "y" if r < 2 * v else "s" if r < n else "c"
    def move(r, f_new):
        f_old = cur[r]
        if f_old == f_new:
            return
        k = dim[f_new] - dim[f_old]
        assert k > 0 and all(span_has(fr[f_new], u) for u in fr[f_old]), "not a nested climb"
        cl = classes(r)
        L["M"] += k
        by[kind(r)]["M"] += k
        if not cl:
            L["I0"] += k; by[kind(r)]["I0"] += k; hist[0] += k
        else:
            L["P"] += k * (len(cl) - 1); by[kind(r)]["P"] += k * (len(cl) - 1)
            hist[len(cl)] += k; maxlab[0] = max(maxlab[0], len(cl))
            for c in cl:
                dd = ((1 if inU(ports[c[1]], f_new) else 0) - (1 if inU(ports[c[1]], f_old) else 0)) if c[0] == "x" else 0
                L["D"] += dd; L["A"] += k - dd; by[kind(r)]["D"] += dd
        cur[r] = f_new
    def mix(targets, fn):
        before = sum(psi(r) for r in targets)
        fn()
        for r in targets:
            reg[r] = {kk: c for kk, c in reg[r].items() if c != 0}
        L["Bmix"] -= sum(psi(r) for r in targets) - before
    def axpy(t, src, k):
        d = reg[t]
        for kk, u in src.items():
            d[kk] = d.get(kk, Fr(0)) + k * u
    psi0 = sum(psi(r) for r in range(n))
    # g1^T at the start frames: s += K^T y
    KT = defaultdict(list)
    for t in range(v):
        for s, c in K[t].items():
            KT[s].append((t, c))
    for s, lst in KT.items():
        assert cur[s] == cur[v + lst[0][0]] == 0, "g1 joins a slot that is not at frame 0"
    mix(list(KT), lambda: [axpy(s, dict(reg[v + t]), c) for s, lst in KT.items() for t, c in lst])
    def regs(g):
        return [g[2]] + [t for t, _, _ in g[3]] + (g[4] if g[0] == "out" else [])
    def gate_T(g):
        f = g[1]
        for r in regs(g):
            move(r, f)
        if g[0] == "out":        # forward: t += k src for each target -> src -= sum k t
            snap = {t: dict(reg[t]) for t, _, _ in g[3]}
            mix([g[2]], lambda: [axpy(g[2], snap[t], -Fr(a, b)) for t, a, b in g[3]])
        else:                    # forward: tgt += sum k s -> s -= k tgt for each s
            tg = dict(reg[g[2]])
            mix([s for s, _, _ in g[3]], lambda: [axpy(s, tg, -Fr(a, b)) for s, a, b in g[3]])
    for g in D["A"]:
        gate_T(g)
    sc = scat(D)
    for k, (kk, s, f) in enumerate(rets):
        assert cur[s] == f
        c = n + k
        cur[c] = 0
        assert all(cur[v + t] == 0 for t in range(v))
        mix([c], lambda: [axpy(c, dict(reg[v + t]), -Fr(a, b)) for t in range(v) for k2, a, b in sc[t] if k2 == kk])
        move(c, f)
        mix([s, c], lambda: (axpy(s, dict(reg[c]), Fr(1)), reg[c].clear()))
    for g in D["B"]:
        gate_T(g)
    for r in range(n):
        move(r, final[r])
    # g8^T at the end frames: x += J^T s (sum over slots r of J[r][t] s_r), then s := Fss^T s
    xs = [t for t in range(v)]
    snap = {r: dict(reg[r]) for r in range(2 * v, n)}
    def g8():
        for r, row in J.items():
            for t, c in row.items():
                axpy(t, snap[r], c)
        new = {r: dict() for r in range(2 * v, n)}
        for r, row in Fss.items():
            for s, c in row.items():                 # (Fss^T s)_s = sum_r Fss[r][s] s_r
                for kk, u in snap[r].items():
                    new[s][kk] = new[s].get(kk, Fr(0)) + c * u
        for r in range(2 * v, n):
            reg[r] = new[r]
    mix(list(range(v)) + list(range(2 * v, n)), g8)
    psi_end = sum(psi(r) for r in range(nreg))
    def minus(a, b):
        o = dict(a)
        for kk, u in b.items():
            o[kk] = o.get(kk, 0) - u
        return {kk: u for kk, u in o.items() if u != 0}
    okx = all(reg[t] == minus(x0[t], y0[t]) for t in range(v))
    oky = all(reg[v + t] == y0[t] for t in range(v))
    oks = all(reg[2 * v + i] == s0[i] for i in range(R))
    rhs = psi_end - psi0 + L["Bmix"] + 2 * L["D"] + L["Z"] - L["P"] + L["I0"]
    assert rhs == L["M"]
    deliv = sum(level(c, final[t]) for t in range(v) for c in classes(t)) + \
        sum(level(c, final[v + t]) for t in range(v) for c in classes(v + t))
    out = dict(M=L["M"], P=L["P"], D=L["D"], Z=L["Z"], I0=L["I0"], A=L["A"], Bmix=str(L["Bmix"]), psi0=psi0,
               psi_end=psi_end, deliv=deliv, B=str(L["Bmix"] + psi_end - deliv), maxlab=maxlab[0],
               hist=dict(sorted(hist.items())), by={k: dict(c) for k, c in by.items()},
               x_is_x_minus_y=okx, y_kept=oky, slots_restored=oks)
    return out


def main():
    p, name = sys.argv[1], sys.argv[2]
    inits = [a for a in sys.argv[3:] if not a.startswith("models=")] or ["B1", "B2"]
    mods = [a[7:] for a in sys.argv[3:] if a.startswith("models=")]
    models = {"clean": (False, False), "dirty": (True, False), "dirty-one": (True, True)}
    combos = [models[x] for x in (mods[0].split(",") if mods else ["clean", "dirty", "dirty-one"])]
    D = load(p)
    res = {}
    for ini in inits:
        for dirt, one in combos:
            r = run_bw(D, ini, dirt, one)
            tag = ini + ("/dirty-one" if one else "/dirty" if dirt else "/clean")
            res[tag] = r
            M = r["M"]
            print("%s backward %s: M=%d P=%d (%.2f per move) B=%s D=%d I0=%d maxlab=%d | x = x - y %s, y kept %s, slots"
                  " restored %s" % (name, tag, M, r["P"], r["P"] / M, r["B"], r["D"], r["I0"], r["maxlab"],
                                    r["x_is_x_minus_y"], r["y_kept"], r["slots_restored"]))
            for k in "xysc":
                b = r["by"].get(k, {})
                if b:
                    print("    %s: M=%d P=%d D=%d I0=%d P/M=%.2f" % (k, b.get("M", 0), b.get("P", 0), b.get("D", 0),
                                                                 b.get("I0", 0), b.get("P", 0) / max(1, b.get("M", 0))))
    json.dump(res, open("out/bw_%s.json" % name, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
