#!/usr/bin/env python3
"""m-gen (13 Oct 2026): EXACT Python mirror of the kernel terms of the mutable-source scalar checks

    GS.scalCheckZ c p lo n,  GS.yCheckZ c p lo n   (x-mutable, XMut/ScalarZ.lean)
    GS.eCheckZ c p gl                              (x-mutable, XMut/ECheckZ.lean)

and, for controls, of the classic GS.scalCheck / GS.yCheck / GS.eCheck (Jacob Sussman's Work/GCert/Scalar/Def.lean
and Seg.lean; Nil Mamano's Work/Gauge/ECheck.lean).  Function by function the same Boolean computation on Python
integers: gPre (parameter guards, mask tm, K = 64 + sw n, m2 = 2^K), gInit / gInitL (tagged first contents), gAdd
(kind test, t /= s, register range of the trie, gW reduced digit weight a/b with its three guards, gDiv exact
division tested by b masked additions, addNC masked additions, gEnc bound on the positive part), gScat (row
count = v, total index < number of retained totals), gFin (packed digit equality, tags dropped), incK.
The Lean kinds are copied literally:
    okAZK  t < nr, s < nr, t not a target, s not a target           (phase A, sources mutable)
    okBZK  t < nr, s < nr, s not a target or t a target             (phase B, sources mutable)
    okSK   target <- slot                                           (scatter)
    okAK / okBK the classic kinds (controls only).
Input: a gcert/1 JSON (A, B as ["out", frame, s, [[t, a, b], ..], ex] with signed numerators a), exactly the data the
generated Lean `raw` holds (the generator's Emit writes these lists; see roundtrip.py for the parse-back check).
Nothing here is trusted by Lean; it predicts the value of `decide +kernel`.
"""
import sys, json, gzip, time, argparse
from math import gcd

M32 = 4294967296


class No(Exception):
    """the kernel term evaluates to false"""


def guard(b, why):
    if not b:
        raise No(why)


def unitK(ux, us, uy, v, v2, r):
    return ux if r + 1 <= v else (uy if r + 1 <= v2 else us)


# ---- kinds (Lean text: Def.lean okAK/okSK/okBK, ScalarZ.lean okAZK/okBZK) -------------------------------
def okAK(v, v2, nr, t, s):
    return t + 1 <= nr and s + 1 <= nr and v2 <= t and (s + 1 <= v or v2 <= s)


def okSK(v, v2, nr, t, s):
    return v <= t and t + 1 <= v2 and v2 <= s and s + 1 <= nr


def okBK(v, v2, nr, t, s):
    return t + 1 <= nr and s + 1 <= nr and (v <= t or s + 1 <= v) and \
        (s + 1 <= v or v2 <= s or (v <= t and t + 1 <= v2))


def okAZK(v, v2, nr, t, s):
    return t + 1 <= nr and s + 1 <= nr and (t + 1 <= v or v2 <= t) and (s + 1 <= v or v2 <= s)


def okBZK(v, v2, nr, t, s):
    return t + 1 <= nr and s + 1 <= nr and (s + 1 <= v or v2 <= s or (v <= t and t + 1 <= v2))


def selY(v, v2, s):
    return v <= s and s + 1 <= v2


def rep(sw, n):
    x = 0
    for _ in range(n):
        x = (x << sw) + 1
    return x


class Run:
    """one evaluation: parameters from gPre and the trie as a dict register -> leaf (unset = 0)"""

    def __init__(self, c, p, lo, n):
        d, sw, ux, us, uy = p
        self.c, self.p = c, p
        v, R = c["v"], c["R"]
        v2, nr, hi = 2 * v, 2 * v + R, lo + n
        half = 2 ** (sw - 1) if sw >= 1 else 1        # Nat.sub sw 1 truncates at 0
        K = 64 + sw * n
        guard(1 <= ux and 1 <= us and 1 <= uy and 1 <= sw, "gPre: units / sw")
        guard(nr <= 2 ** d, "gPre: registers exceed 2^d")
        guard(hi <= v, "gPre: lo + n > v")
        guard(ux + 1 <= half and uy + 1 <= half, "gPre: ux / uy >= 2^(sw-1)")
        guard(ux + 1 <= 2147483648 and uy + 1 <= 2147483648, "gPre: ux / uy >= 2^31")
        guard(all(v2 <= s and s + 1 <= nr for _, s, _ in c["ret"]), "gPre: retOKK")
        self.tm = ((rep(sw, n) << (sw - 1 if sw >= 1 else 0)) << 64) | (2 ** 63)
        self.v, self.v2, self.nr, self.hi, self.K, self.m2, self.d = v, v2, nr, hi, K, 2 ** K, d
        self.sw, self.ux, self.us, self.uy, self.lo, self.n = sw, ux, us, uy, lo, n
        self.tr = {}
        self.adds = 0

    def key_ok(self, r):
        return 0 <= r < 2 ** self.d

    def addNC(self, a, n, acc):
        for _ in range(n):
            acc = acc + a
            guard(acc & self.tm == 0, "addNC: digit reaches 2^(sw-1) (or tag bit 63)")
        return acc

    def gDiv(self, b, X):
        if b == 1:
            return X
        q = X // b
        guard(q & self.tm == 0, "gDiv: quotient mask")
        r = self.addNC(q, b, 0)
        guard(r == X, "gDiv: not divisible by %d" % b)
        return q

    def gAdd(self, t, s, co, ok):
        guard(ok, "kind test fails: %d += %d" % (t, s))
        guard(t != s, "t = s")
        guard(self.key_ok(s) and self.key_ok(t), "register outside the trie")
        ng, num, den = co
        A = num * unitK(self.ux, self.us, self.uy, self.v, self.v2, t)
        B = den * unitK(self.ux, self.us, self.uy, self.v, self.v2, s)
        g = gcd(A, B)
        a = A // g if g else 0
        b = B // g if g else 0
        guard(1 <= B and 1 <= b and a * B == A * b, "gW guards")
        vs, vt = self.tr.get(s, 0), self.tr.get(t, 0)
        ps, ns, pt, nt = vs % self.m2, vs >> self.K, vt % self.m2, vt >> self.K
        qp = self.gDiv(b, ps)
        qn = self.gDiv(b, ns)
        if not ng:
            pt2 = self.addNC(qp, a, pt)
            nt2 = self.addNC(qn, a, nt)
        else:
            pt2 = self.addNC(qn, a, pt)
            nt2 = self.addNC(qp, a, nt)
        guard(pt2 + 1 <= self.m2, "gEnc: positive part too long")
        self.tr[t] = pt2 + (nt2 << self.K)
        self.adds += 1

    def gGates(self, gates, okf, sel):
        for g in gates:
            if g[0] == "out":
                s = g[2]
                if sel(s):
                    for t, k in g[3]:
                        self.gAdd(t, s, k, okf(t, s))
            else:
                t = g[2]
                for s, k in g[3]:
                    if sel(s):
                        self.gAdd(t, s, k, okf(t, s))

    def gScat(self, rows, ret, okf):
        rl = len(ret)
        t = self.v
        for row in rows:
            for k, co in row:
                guard(k + 1 <= rl, "gScat: total index")
                s = ret[k][0]
                self.gAdd(t, s, co, okf(t, s))
            t += 1
        guard(t == self.v2, "gScat: number of rows is not v")

    def gInit(self, u, off, n):
        for i in range(n):
            r = off + i
            guard(self.key_ok(r), "gInit: register outside the trie")    # Trie.set on a deeper key: never here
            self.tr[r] = (u << (64 + self.sw * i)) + u * ((2654435761 * (i + 1)) % M32)

    def gInitL(self, u, gl):
        for i, r in enumerate(gl):
            self.tr[r] = (u << (64 + self.sw * i)) + u * ((2654435761 * (i + 1)) % M32)

    def gFin(self, u, base, lo, hi):
        for S in range(self.v):
            r = base + S
            guard(self.key_ok(r), "gFin: register outside the trie")
            x = self.tr.get(r, 0)
            e = (u << (self.sw * (S - lo))) if (lo <= S and S + 1 <= hi) else 0
            if ((x % self.m2) >> 64) != ((x >> self.K) >> 64) + e:
                raise No("gFin: register %d (class base %d) does not hold exactly its source" % (r, base))

    def items(self, okA, okB):
        c, v, v2, nr = self.c, self.v, self.v2, self.nr
        self.gGates(c["A"], lambda t, s: okA(v, v2, nr, t, s), lambda s: True)
        self.gScat(c["rows"], c["retL"], lambda t, s: okSK(v, v2, nr, t, s))
        self.gGates(c["B"], lambda t, s: okB(v, v2, nr, t, s), lambda s: True)


def lean_data(c):
    """the Lean-facing fields: coefficients (ng, num, den), ret as a list (register, frame) in total order"""
    co = lambda a, b: (a < 0, abs(a), b)
    gate = lambda g: ("out", g[1], g[2], [(t, co(a, b)) for t, a, b in g[3]]) if g[0] == "out" else \
        ("in", g[1], g[2], [(s, co(a, b)) for s, a, b in g[3]])
    assert all(g[0] in ("out", "in") for g in c["A"] + c["B"])
    return dict(v=c["v"], R=c["R"], A=[gate(g) for g in c["A"]], B=[gate(g) for g in c["B"]],
                rows=[[(k, co(a, b)) for k, a, b in row] for row in c["scat"]["table"]],
                retL=[(s, f) for _, s, f in sorted(c["ret"])], ret=c["ret"])


def scalCheck(c, p, lo, n, mutable=True):
    R = Run(c, p, lo, n)
    R.gInit(R.ux, lo, n)
    R.items(okAZK if mutable else okAK, okBZK if mutable else okBK)
    if not mutable:
        R.gFin(R.ux, 0, lo, R.hi)
    R.gFin(R.uy, R.v, lo, R.hi)
    return R.adds


def yCheck(c, p, lo, n, mutable=True):
    R = Run(c, p, lo, n)
    R.gInit(R.uy, R.v + lo, n)
    okB = okBZK if mutable else okBK
    v, v2, nr = R.v, R.v2, R.nr
    R.gGates(c["B"], lambda t, s: okB(v, v2, nr, t, s), lambda s: selY(v, v2, s))
    R.gFin(R.uy, R.v, lo, R.hi)
    return R.adds


def eCheck(c, p, gl, mutable=True):
    R = Run(c, p, 0, len(gl))
    lo = R.v2
    for r in gl:                                       # incK nr gl v2
        guard(lo <= r and r + 1 <= R.nr, "incK: gauge list not increasing / not slots")
        lo = r + 1
    guard(R.us + 1 <= 2 ** (R.sw - 1 if R.sw >= 1 else 0), "eCheck: us >= 2^(sw-1)")
    guard(R.us + 1 <= 2147483648, "eCheck: us >= 2^31")
    R.gInitL(R.us, gl)
    R.items(okAZK if mutable else okAK, okBZK if mutable else okBK)
    R.gFin(R.uy, R.v, 0, 0)
    return R.adds


def evaluate(f, *a, **k):
    t0 = time.time()
    try:
        n = f(*a, **k)
        return dict(value=True, adds=n, seconds=round(time.time() - t0, 2))
    except No as e:
        return dict(value=False, why=str(e), seconds=round(time.time() - t0, 2))


def suite(c, spar, nb, espar, eb, mutable=True):
    """the evaluations of the build: scalar batches, y check, E7 batches"""
    v = c["v"]
    step = -(-v // nb)
    out = dict(spar=spar, espar=espar, scal=[], ychk=None, echk=[])
    for lo in range(0, v, step):
        n = min(step, v - lo)
        out["scal"].append(dict(lo=lo, n=n, **evaluate(scalCheck, c, spar, lo, n, mutable)))
    out["ychk"] = dict(lo=0, n=v, **evaluate(yCheck, c, spar, 0, v, mutable))
    gl = [r for r, _ in c["ext"]]
    for i in range(0, len(gl), eb):
        out["echk"].append(dict(batch=i // eb, size=len(gl[i:i + eb]), **evaluate(eCheck, c, espar, gl[i:i + eb], mutable)))
    out["all_true"] = all(x["value"] for x in out["scal"] + out["echk"] + [out["ychk"]])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cert")
    ap.add_argument("--spar", default="10,21,4,8,8")
    ap.add_argument("--batches", type=int, default=4)
    ap.add_argument("--espar", default="10,8,4,8,8")
    ap.add_argument("--ebatch", type=int, default=80)
    ap.add_argument("--classic", action="store_true", help="the classic checks (controls)")
    ap.add_argument("--margin", action="store_true", help="also evaluate at sw - 1 (expected false)")
    a = ap.parse_args()
    c0 = json.loads(gzip.open(a.cert, "rt").read())
    c = lean_data(c0)
    c["ext"] = c0["ext"]
    sp = tuple(map(int, a.spar.split(",")))
    ep = tuple(map(int, a.espar.split(",")))
    res = dict(cert=a.cert, mutable=not a.classic, main=suite(c, sp, a.batches, ep, a.ebatch, not a.classic))
    if a.margin:
        res["sw_minus_1"] = suite(c, (sp[0], sp[1] - 1) + sp[2:], a.batches, (ep[0], ep[1] - 1) + ep[2:], a.ebatch,
                                  not a.classic)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
