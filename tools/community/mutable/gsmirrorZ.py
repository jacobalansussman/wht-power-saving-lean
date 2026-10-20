#!/usr/bin/env python3
"""x-mutable: Python mirror of GS.scalCheckZ / GS.yCheckZ (XMut/ScalarZ.lean): Sussman's gsmirror.replay
(r5 tools/gx/gsmirror.py, own code) with the mutable-source kinds okAZ / okBZ and WITHOUT the x-final test (X3).
Searches register units (ux, us, uy) for which both packed replays pass; prints the digit width sw they need.
usage: gsmirrorZ.py <gcert1.json.gz>"""
import sys, json, gzip
from math import gcd

class Fail(Exception): pass

def adds_of(g):
    if g[0] == "out": return [(t, g[2], a, b) for t, a, b in g[3]]
    return [(g[2], s, a, b) for s, a, b in g[3]]

def replay(c, units, mode, lo=0, n=None):
    v, R = c["v"], c["R"]; v2, nr = 2 * v, 2 * v + R
    n = v if n is None else n
    ux, us, uy = units
    unit = lambda r: ux if r < v else uy if r < v2 else us
    P = [dict() for _ in range(nr)]; N = [dict() for _ in range(nr)]; mx = [max(ux, uy)]
    if mode == "x":
        for t in range(lo, lo + n): P[t][t] = ux
    elif mode == "e":                       # E7 (eCheckZ): sources = the gauged slot registers gl = lo (a list)
        for T, r in enumerate(lo): P[r][T] = us
    else:
        for t in range(lo, lo + n): P[v + t][t] = uy
    def add(t, s, a, b, ok):
        if not ok or t == s or not (t < nr and s < nr) or b <= 0:
            raise Fail("kind / bounds: %d += (%d/%d) %d" % (t, a, b, s))
        A, B = abs(a) * unit(t), b * unit(s); g = gcd(A, B); wa, wb = A // g, B // g
        sp, sn = (P[s], N[s]) if a >= 0 else (N[s], P[s])
        if wb > 1 and (any(x % wb for x in sp.values()) or any(x % wb for x in sn.values())):
            raise Fail("digit not divisible by %d (%d += %d/%d %d)" % (wb, t, a, b, s))
        if wa == 0: return
        for dst, src in ((P[t], sp), (N[t], sn)):
            for k, x in src.items():
                if x == 0: continue
                y = dst.get(k, 0) + wa * (x // wb); dst[k] = y
                if y > mx[0]: mx[0] = y
    okAZ = lambda t, s: (t < v or v2 <= t) and (s < v or v2 <= s)
    okS = lambda t, s: v <= t < v2 and v2 <= s < nr
    okBZ = lambda t, s: s < v or v2 <= s or (v <= t < v2)
    if mode in ("x", "e"):
        for g in c["A"]:
            for t, s, a, b in adds_of(g): add(t, s, a, b, okAZ(t, s))
        rows = c["scat"]["table"]
        for t, row in enumerate(rows):
            for k, a, b in row:
                s = c["ret"][k][1]; add(v + t, s, a, b, okS(v + t, s))
    for g in c["B"]:
        for t, s, a, b in adds_of(g):
            if mode in ("x", "e") or v <= s < v2: add(t, s, a, b, okBZ(t, s))
    net = lambda r: {k: x - N[r].get(k, 0) for k, x in P[r].items() if x != N[r].get(k, 0)}
    for t in range(v):
        if mode == "e":
            if net(v + t): raise Fail("E7: a gauged slot reaches y %d" % t)
            continue
        want = {t: uy} if lo <= t < lo + n else {}
        if net(v + t) != want: raise Fail("identity at y %d" % t)
    return mx[0]

def main(path):
    c = json.loads(gzip.open(path, "rt").read()); v = c["v"]
    res = []
    for ux in (1, 2, 4, 8, 16):
        for us in (1, 2, 4, 8, 16):
            for uy in (1, 2, 4, 8, 16):
                try:
                    m1 = max(replay(c, (ux, us, uy), "x", lo, 28) for lo in range(0, v, 28))
                    m2 = replay(c, (ux, us, uy), "y")
                except Fail as e:
                    continue
                res.append((max(m1, m2).bit_length() + 1, ux, us, uy, m1, m2))
    res.sort()
    gl = [r for r, f in c["ext"]]
    eres = []
    for ux in (1, 4, 8, 16, 32):
        for us in (1, 2, 4, 8):
            try:
                m = max(replay(c, (ux, us, 8), "e", gl[i:i + 80]) for i in range(0, len(gl), 80))
                eres.append((m.bit_length() + 1, ux, us, m))
            except Fail as e:
                pass
    eres.sort()
    print(json.dumps({"cert": path, "passing_units": len(res), "best": res[:5], "E7_batches_of_80_pass_(sw,ux,us,maxdigit)": eres[:4]}))

if __name__ == "__main__":
    main(sys.argv[1])
