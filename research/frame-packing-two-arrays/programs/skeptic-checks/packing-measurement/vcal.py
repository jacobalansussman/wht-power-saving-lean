#!/usr/bin/env python3
"""v-fp-measure: own exact engine (Gaussian-rational 2^m x 2^m matrices, no code of sl-combinatorics or fp-measure)
for the small calibration units: raw-gauge ledger P, B, D, Z, I0 with frame labels (Sp(2m,2) images), exactness
of the unit (every delivering register = unit scalar * K y_pi).  Words transcribed from the definitions in
sl-combinatorics t_known.py (main()), the o3 control of PROOF.md Example 3.4.
usage: python3 -I vcal.py"""
from fractions import Fraction as Fr
import json, itertools

ZERO = (Fr(0), Fr(0))


def cm(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def ca_(a, b):
    return (a[0] + b[0], a[1] + b[1])


def cinv(a):
    d = a[0] * a[0] + a[1] * a[1]
    return (a[0] / d, -a[1] / d)


def cconj(a):
    return (a[0], -a[1])


class Eng:
    def __init__(self, m):
        self.m, self.n = m, 1 << m
        n = self.n
        self.I = [[(Fr(1), Fr(0)) if i == j else ZERO for j in range(n)] for i in range(n)]
        ca, cb = (Fr(1, 2), Fr(1, 2)), (Fr(1, 2), Fr(-1, 2))
        self.dirc = {}
        for z in range(1, n):
            M = [[ZERO] * n for _ in range(n)]
            for x in range(n):
                M[x][x] = ca_(M[x][x], ca)
                M[x ^ z][x] = ca_(M[x ^ z][x], cb)
            self.dirc[z] = M
        K = self.I
        for j in range(m):
            K = self.mul(self.dirc[1 << j], K)
        self.K = K

    def mul(self, A, B):
        n = self.n
        out = [[ZERO] * n for _ in range(n)]
        for i in range(n):
            Ai = A[i]
            for k in range(n):
                a = Ai[k]
                if a == ZERO:
                    continue
                Bk = B[k]
                row = out[i]
                for j in range(n):
                    if Bk[j] != ZERO:
                        row[j] = ca_(row[j], cm(a, Bk[j]))
        return out

    def dag(self, A):
        n = self.n
        return [[cconj(A[j][i]) for j in range(n)] for i in range(n)]

    def addr(self, images):          # address map x -> A x, A given by images of e_j
        n = self.n
        M = [[ZERO] * n for _ in range(n)]
        for x in range(n):
            y = 0
            for j in range(self.m):
                if x >> j & 1:
                    y ^= images[j]
            M[y][x] = (Fr(1), Fr(0))
        return M

    def norm(self, g):
        """(key, scalar) with g = scalar * normalised(g)"""
        for i in range(self.n):
            for j in range(self.n):
                if g[i][j] != ZERO:
                    s = g[i][j]
                    si = cinv(s)
                    key = tuple(tuple(cm(e, si) for e in row) for row in g)
                    return key, s
        raise ValueError

    def pauli(self, a, b):
        n = self.n
        M = [[ZERO] * n for _ in range(n)]
        for x in range(n):
            sgn = -1 if bin(b & x).count("1") & 1 else 1
            M[x ^ a][x] = (Fr(sgn), Fr(0))
        return M

    def symp(self, g):
        """Sp(2m,2) image: images (a|b) of X^{e_i} and Z^{e_i}"""
        gd = self.dag(g)
        out = []
        for kind in "XZ":
            for i in range(self.m):
                P = self.pauli(1 << i, 0) if kind == "X" else self.pauli(0, 1 << i)
                C = self.mul(self.mul(g, P), gd)
                col = [C[r][0] for r in range(self.n)]
                a = [r for r in range(self.n) if col[r] != ZERO]
                assert len(a) == 1
                a = a[0]
                b = 0
                for j in range(self.m):
                    x = 1 << j
                    ratio = cm(C[x ^ a][x], cinv(C[a][0]))
                    assert ratio in ((Fr(1), Fr(0)), (Fr(-1), Fr(0)))
                    if ratio[0] < 0:
                        b |= x
                out.append((a, b))
        return tuple(out)

    def level(self, sp):
        zs = sp[self.m:]
        rows = []
        for a, _ in zs:
            x = a
            for r in rows:
                x = min(x, x ^ r)
            if x:
                rows.append(x)
                rows.sort(reverse=True)
        return len(rows)


def run(E, W, R, word):
    """contents: register -> {(key, input): [g, coef]}"""
    reg = [dict() for _ in range(R)]
    for a in range(W):
        k, s = E.norm(E.I)
        reg[a][(k, a)] = [E.I, (Fr(1), Fr(0))]
    spc = {}
    def labels(r):
        out = {}
        for (k, l), (g, c) in reg[r].items():
            if k not in spc:
                sp = E.symp(g)
                spc[k] = (sp, E.level(sp))
            out[spc[k][0]] = spc[k][1]
        return out
    def psi(r):
        return sum(labels(r).values())
    L = dict(M=0, P=0, D=0, Z=0, A=0, I0=0, Bmix=0)
    def put(d, g, c, l):
        k, s = E.norm(g)
        key = (k, l)
        c = cm(c, s)
        g0 = [list(row) for row in k]
        if key in d:
            nc = ca_(d[key][1], c)
            if nc == ZERO:
                del d[key]
            else:
                d[key][1] = nc
        else:
            d[key] = [g0, c]
    for st in word:
        if st[0] == "move":
            _, r, z = st
            L["M"] += 1
            before = labels(r)
            if not before:
                L["I0"] += 1
            new = {}
            for (k, l), (g, c) in reg[r].items():
                put(new, E.mul(E.dirc[z], g), c, l)
            reg[r] = new
            after = labels(r)
            assert len(after) == len(before)
            if before:
                L["P"] += len(before) - 1
                # pair labels by left multiplication: compute per term
                for (k, l), (g, c) in list(reg[r].items()):
                    pass
                # level changes per label: map each old label to new via the letter's symplectic image
                gz = E.dirc[z]
                seen = set()
                for (k, l), (g, c) in new.items():
                    sp_new = spc[k][0]
                    if sp_new in seen:
                        continue
                    seen.add(sp_new)
                    # old label = dir z^-1 * g ; level(old)
                    gold = E.mul(E.dag(gz), g)
                    ko, _ = E.norm(gold)
                    if ko not in spc:
                        sp = E.symp(gold); spc[ko] = (sp, E.level(sp))
                    d = spc[k][1] - spc[ko][1]
                    if d > 0: L["A"] += 1
                    elif d < 0: L["D"] += 1
                    else: L["Z"] += 1
        elif st[0] == "mix":
            _, r, terms = st
            before = psi(r)
            new = {}
            for s, mu, letter in terms:
                mu = (Fr(mu), Fr(0))
                for (k, l), (g, c) in reg[s].items():
                    gg = g if letter is None else E.mul(letter, g)
                    put(new, gg, cm(mu, c), l)
            reg[r] = new
            L["Bmix"] -= psi(r) - before
    # exactness: delivering registers 0..W-1 each hold one term proportional to K y_pi
    kK, sK = E.norm(E.K)
    perm, ok = [], True
    for r in range(W):
        items = list(reg[r].items())
        if len(items) != 1 or items[0][0][0] != kK:
            ok = False
            break
        (k, l), (g, c) = items[0]
        u = cm(c, cinv(sK))     # content = c * normalised(K) = u * K
        perm.append(l)
        if u[0] * u[0] + u[1] * u[1] != 1:
            ok = False
    ok = ok and sorted(perm) == list(range(W))
    psiE = sum(psi(r) for r in range(R))
    B = L["Bmix"] + psiE - W * E.m
    Delta = W * E.m - L["M"]
    return dict(L, B=B, Delta=Delta, unit=ok, identity=(Delta == L["P"] - B - 2 * L["D"] - L["Z"] - L["I0"]))


def words(E):
    m = E.m
    out = []
    def plain(W):
        w = []
        for r in range(W):
            w += [("move", r, 1 << j) for j in range(m)]
        return w
    for W in (1, 2, 3):
        out.append(("plain W=%d" % W, W, W, plain(W)))
    if m >= 2:
        w = [("move", 0, 1), ("mix", 0, [(0, 1, None), (1, 1, None)])]
        for j in range(1, m):
            w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
        w += [("mix", 0, [(0, 1, None), (1, -1, None)]), ("move", 1, 1)]
        out.append(("gauge exchange", 2, 2, w))
        w = [("mix", 2, [(1, 1, None)]), ("mix", 0, [(0, 1, None), (1, 1, None)]), ("move", 0, 1),
             ("mix", 0, [(0, 1, None), (2, -1, None)]), ("move", 0, 2), ("move", 2, 2),
             ("mix", 0, [(0, 1, None), (2, 1, None)]), ("move", 1, 1), ("move", 1, 2),
             ("mix", 0, [(0, 1, None), (1, -1, None)])]
        for j in range(2, m):
            w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
        out.append(("packing gadget", 2, 3, w))
    if m >= 3:
        s = E.addr([2, 4, 1] + [1 << j for j in range(3, m)])
        s2 = E.mul(s, s)
        h = Fr(1, 2)
        w = [("mix", 0, [(0, 1, None), (1, 1, None)]), ("mix", 1, [(0, 1, None), (1, -1, None), (1, -1, s)])]   # sequential: r1 = (x1 + x2) - x2 - s x2
        for j in range(m):
            w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
        # K x2 = (1 - s + s^2)/2 (r0 - r1) = (1 - s + s^2)(1 + s)/2 K x2;  then K x1 = r0 - K x2
        w += [("mix", 1, [(0, h, None), (1, -h, None), (0, -h, s), (1, h, s), (0, h, s2), (1, -h, s2)]),
              ("mix", 0, [(0, 1, None), (1, -1, None)])]
        out.append(("order-3 control C%dO3" % (2 * m), 2, 2, w))
        w = [("mix", 0, [(0, 1, None), (1, 1, None)]), ("mix", 1, [(0, 1, None), (1, -2, None)])]   # r1 = x1 - x2
        for j in range(m):
            w += [("move", 0, 1 << j), ("move", 1, 1 << j)]
        w += [("mix", 1, [(0, h, None), (1, -h, None)]), ("mix", 0, [(0, 1, None), (1, -1, None)])]
        out.append(("+-1 halving C%dH" % (2 * m), 2, 2, w))
    return out


def main():
    res = []
    for m in (2, 3, 4):
        E = Eng(m)
        for name, W, R, w in words(E):
            r = run(E, W, R, w)
            r.update(case=name, m=m, W=W)
            res.append(r)
            print("m=%d %-24s W=%d M=%2d Delta=%2d P=%d B=%d D=%d Z=%d I0=%d unit=%s identity=%s"
                  % (m, name, W, r["M"], r["Delta"], r["P"], r["B"], r["D"], r["Z"], r["I0"], r["unit"], r["identity"]))
    json.dump(res, open("out/vcal.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
