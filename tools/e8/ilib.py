"""ilib.py -- load a saved hierarchy (JSON, see hio.py), compile it, prune target stops, evaluate."""
import sys, math
sys.dont_write_bytecode = True
from lib import *
import comp, hio
from collections import Counter, defaultdict
class H: pass
def load(path="hierarchy.json"):
    gd = hio.load_h(path)
    h = H(); h.content, h.frame, h.kids, h.need, h.twin = gd["content"], gd["frame"], gd["kids"], gd["need"], dict(gd.get("twin", {}))
    h.need = [dict(n) for n in h.need]
    return h
def save(h, path):
    hio.dump_h(dict(content=h.content, frame=h.frame, kids=h.kids, need=h.need, twin=h.twin), path)
NEED = [need(S) for S in range(V)]
def match(h, u, S):
    """coef c such that c * content[u] is part of need(S), else None"""
    nd = NEED[S]; c = None
    for T, x in h.content[u].items():
        y = nd.get(T)
        if y is None: return None
        r = y / x
        if c is None: c = r
        elif c != r: return None
    return c
def weight(c):
    inv = e8.inv_hist(c)
    return 5 * sum(n * r * math.log(45.0 / r) for r, n in inv.items())
def build(h, ylow="all", drop=frozenset(), name="e8 search"):
    st = {}
    d = comp.compile_design(h, name=name, ylow=ylow, stats=st, drop=drop)
    c = e8.build_cert(d)
    return c, st
def prune_low(h, name="e8 search", verbose=False):
    """the pruning pass: drop target stops that do not pay"""
    def trial(drop):
        c, st = build(h, "all", drop, name)
        return weight(c), c, st
    drop = set()
    w, c, st = trial(drop)
    low0 = list(st["lowset"]); byS = {}
    for S, u in low0: byS.setdefault(S, []).append(u)
    for S in sorted(byS):
        d2 = drop | set((S, u) for u in byS[S])
        w2, c2, st2 = trial(d2)
        if w2 < w - 1e-9: drop, w, c, st = d2, w2, c2, st2
    for S, u in low0:
        if (S, u) in drop: continue
        d2 = drop | {(S, u)}
        w2, c2, st2 = trial(d2)
        if w2 < w - 1e-9: drop, w, c, st = d2, w2, c2, st2
    for S in sorted(byS):
        d2 = drop - set((S, u) for u in byS[S])
        if d2 != drop:
            w2, c2, st2 = trial(d2)
            if w2 < w - 1e-9: drop, w, c, st = d2, w2, c2, st2
    return c, st, drop, w
def figure(c):
    q = e8.quick(c)
    ok = q["accepted_ref"] and q["accepted_mirror"]
    return (q["whole_block"] if ok else 0), q
