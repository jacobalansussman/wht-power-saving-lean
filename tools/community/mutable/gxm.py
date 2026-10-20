#!/usr/bin/env python3
"""x-mutable: exact checks of a gcert/1 word under the MUTABLE-SOURCE interface (own code; reuses gx.py helpers
imported by path).  Everything is exact rational arithmetic.

  hist(c)            the block histogram (climbs), exactly as gx.check1 counts it (E1 nesting asserted)
  kinds(c)           gate kinds (reader class, written class) per phase
  structure(c, mut)  E0/E1/E2 + kind rules: classic (gx E3) or mutable (no gate reads a target into x or s;
                     phase A names no target)
  hid_hx(c)          x-column replay: hid (y_t receives exactly x_t) and whether x is restored (classic hx)
  hy(c)              y-column replay: the y columns of the product are unit columns
  word tests         forward / backward words with the DERIVED gates (dirt gate y -= K s at frame 0, undo
                     z := (M^-1 z) at the full frame) on random rational states: (x,y,s) -> (x, y+x, s) and
                     (x,y,s) -> (x-y, y, s)
"""
import random, importlib.util, sys
from fractions import Fraction as Fr
from collections import Counter

def load_gx(path):
    spec = importlib.util.spec_from_file_location("gx", path)
    gx = importlib.util.module_from_spec(spec); spec.loader.exec_module(gx); return gx

def cls_of(c):
    v = c["v"]
    return lambda r: "x" if r < v else "y" if r < 2 * v else "s"

def adds(g):
    """elementary adds of a gate: list of (written, read, coefficient)"""
    if g[0] == "out":
        return [(t, g[2], Fr(a, b)) for t, a, b in g[3]]
    if g[0] == "in":
        return [(g[2], s, Fr(a, b)) for s, a, b in g[3]]
    return []

def gate_regs(g):
    if g[0] == "neg": return [g[2]]
    if g[0] == "out": return [g[2]] + [t for t, _, _ in g[3]] + list(g[4])
    return [g[2]] + [s for s, _, _ in g[3]]

def hist(c, gx):
    h, v, R = c["h"], c["v"], c["R"]
    fr = [tuple(f) for f in c["frames"]]; dim = [len(f) for f in fr]
    cur = list(c["start"]); H = Counter(); per = {k: Counter() for k in "xysc"}; cl = cls_of(c)
    steps = set()
    def climb(r, f):
        o = cur[r]
        if o != f:
            if (o, f) not in steps:
                assert dim[o] < dim[f] and gx.inside(fr[o], fr[f]), ("E1 not nested", r, o, f)
                steps.add((o, f))
            H[dim[f] - dim[o]] += 1; per[cl(r)][dim[f] - dim[o]] += 1; cur[r] = f
    for ph in ("A", "B"):
        if ph == "B":
            for k, s, f in c["ret"]:
                assert cur[s] == f, "E2: retained slot not at its frame"
                H[dim[f]] += 1; per["c"][dim[f]] += 1
        for g in c[ph]:
            for r in gate_regs(g): climb(r, g[1])
    for r in range(2 * v + R): climb(r, c["final"][r])
    return dict(sorted(H.items())), {k: dict(sorted(w.items())) for k, w in per.items()}

def kinds(c):
    cl = cls_of(c); out = {}
    for ph in ("A", "B"):
        out[ph] = Counter("%s->%s" % (cl(s), cl(t)) for g in c[ph] for t, s, _ in adds(g))
    return out

def structure(c, mutable):
    """returns list of violations (empty = pass)"""
    v, h = c["v"], c["h"]; cl = cls_of(c); bad = []
    for ph in ("A", "B"):
        for i, g in enumerate(c[ph]):
            regs = gate_regs(g)
            if len(set(regs)) != len(regs): bad.append("E3 repeated register %s%d" % (ph, i))
            for t, s, _ in adds(g):
                kd = (cl(s), cl(t))
                if mutable:
                    ok = (kd[0] != "y" or kd[1] == "y") and (ph == "B" or "y" not in kd)
                else:
                    ok = kd in {("x", "s"), ("s", "s"), ("s", "y"), ("x", "y"), ("x", "x"), ("y", "y")} and \
                         (ph == "B" or kd in {("x", "s"), ("s", "s")})
                if not ok: bad.append("kind %s->%s at %s%d" % (kd + (ph, i)))
    for k, s, f in c["ret"]:
        if cl(s) != "s": bad.append("E2 retained total %d is not a slot (register %d)" % (k, s))
    for t in range(v):
        if c["final"][t] != 1: bad.append("E1 x final")
    return bad

# ---------------------------------------------------------------- linear simulation on full states
def n_regs(c): return 2 * c["v"] + c["R"]

def scat_rows(c):
    return [[(k, Fr(a, b)) for k, a, b in row] for row in c["scat"]["table"]]

def run_gates(gates, u, inverse=False):
    seq = reversed(gates) if inverse else gates
    for g in seq:
        if g[0] == "neg":
            u[g[2]] = -u[g[2]]; continue
        ad = adds(g)
        if g[0] == "out":     # all reads are of the source register, which no add of this gate writes
            for t, s, co in ad: u[t] += (-co if inverse else co) * u[s]
        else:                 # 'in': t += sum co*s; sources distinct from t
            acc = sum(co * u[s] for t, s, co in ad)
            u[g[2]] += -acc if inverse else acc
    return u

def run_gates_T(gates, u, inverse=False):
    """apply (product)^T (inverse=False) or (product^-1)^T (inverse=True) to u"""
    seq = gates if inverse else list(reversed(gates))
    for g in seq:
        if g[0] == "neg":
            u[g[2]] = -u[g[2]]; continue
        ad = adds(g)
        if g[0] == "out":     # transposed: source s += co * t (for each t); all t distinct from s
            u[g[2]] += sum((-co if inverse else co) * u[t] for t, s, co in ad)
        else:
            for t, s, co in ad: u[s] += (-co if inverse else co) * u[t]
    return u

def scatter(c, u, inverse=False, transpose=False):
    v = c["v"]; rows = scat_rows(c); ret = {k: s for k, s, _ in c["ret"]}
    if not transpose:
        tot = {k: u[s] for k, s in ret.items()}
        for t, row in enumerate(rows):
            add = sum(co * tot[k] for k, co in row)
            u[v + t] += -add if inverse else add
    else:
        acc = Counter()
        for t, row in enumerate(rows):
            for k, co in row: acc[k] += co * u[v + t]
        for k, s in ret.items(): u[s] += -acc[k] if inverse else acc[k]
    return u

def M(c, u):      return run_gates(c["B"], scatter(c, run_gates(c["A"], list(u))))
def Minv(c, u):   return run_gates(c["A"], scatter(c, run_gates(c["B"], list(u), True), True), True)
def MT(c, u):     return run_gates_T(c["A"], scatter(c, run_gates_T(c["B"], list(u)), transpose=True))

def forward_word(c, x, y, s, undo="z"):
    v = c["v"]; n = n_regs(c); Y = range(v, 2 * v)
    zs = [Fr(0)] * (2 * v) + list(s)
    Ks = M(c, zs)                                   # K s (y rows, slot columns): x = y = 0
    u = list(x) + [y[t] - Ks[v + t] for t in range(v)] + list(s)      # dirt gate g1 at frame 0
    u = M(c, u)                                     # A, scatter, B
    w = Minv(c, u)
    if undo == "z":   # derived undo on x AND slots (mutable sources)
        u = [w[i] if not (v <= i < 2 * v) else u[i] for i in range(n)]
    else:             # the classic engine's g8: slots only
        u = [w[i] if i >= 2 * v else u[i] for i in range(n)]
    return u

def backward_word(c, x, y, s, undo="z"):
    v = c["v"]; n = n_regs(c)
    u = list(x) + list(y) + list(s)
    # (g1i M)^T: slots += (M^T (0, y, 0)) restricted to slots
    m = MT(c, [Fr(0)] * v + list(y) + [Fr(0)] * c["R"])
    u = [u[i] + (m[i] if i >= 2 * v else 0) for i in range(n)]
    u = run_gates_T(c["A"], u, inverse=True)        # MAi^T
    u = scatter(c, u, inverse=True, transpose=True) # scatMi^T
    u = run_gates_T(c["B"], u, inverse=True)        # MBi^T
    # (h8 M)^T: z := (M^T (z, 0)) on z (both variants read only the undone rows)
    # h8^T = P_keep + M^T P_undo  (undo rows: x and slots for the mutable engine, slots only for the classic one)
    keep = (lambda i: v <= i < 2 * v) if undo == "z" else (lambda i: i < 2 * v)
    zz = [Fr(0) if keep(i) else u[i] for i in range(n)]
    m = MT(c, zz)
    return [(u[i] if keep(i) else 0) + m[i] for i in range(n)]

def word_tests(c, trials=2, seed=1, undo="z"):
    rnd = random.Random(seed); v, R = c["v"], c["R"]; res = []
    for _ in range(trials):
        x = [Fr(rnd.randint(-10**6, 10**6)) for _ in range(v)]
        y = [Fr(rnd.randint(-10**6, 10**6)) for _ in range(v)]
        s = [Fr(rnd.randint(-10**6, 10**6)) for _ in range(R)]
        f = forward_word(c, x, y, s, undo)
        ok_f = f == x + [y[t] + x[t] for t in range(v)] + s
        b = backward_word(c, x, y, s, undo)
        ok_b = b == [x[t] - y[t] for t in range(v)] + y + s
        res.append((ok_f, ok_b))
    return res

def hid_hx(c):
    """x-column replay (columns of M on x inputs): hid and classic hx(x-columns)"""
    v = c["v"]; n = n_regs(c); hid = True; hxcols = True
    for t in range(v):
        u = [Fr(0)] * n; u[t] = Fr(1); u = M(c, u)
        for S in range(v):
            if u[v + S] != (1 if S == t else 0): hid = False
        for r in range(v):
            if u[r] != (1 if r == t else 0): hxcols = False
    return hid, hxcols

def hy(c):
    v = c["v"]; n = n_regs(c)
    for t in range(v):
        u = [Fr(0)] * n; u[v + t] = Fr(1); u = M(c, u)
        if any(u[i] != (1 if i == v + t else 0) for i in range(n)): return False
    return True

def hx_full(c):
    """classic hx on ALL columns: the x rows of M are unit rows (row replay through M^T)"""
    v = c["v"]; n = n_regs(c)
    for t in range(v):
        u = [Fr(0)] * n; u[t] = Fr(1); u = MT(c, u)
        if any(u[i] != (1 if i == t else 0) for i in range(n)): return False
    return True
