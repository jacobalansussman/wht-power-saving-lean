"""v-xc3-s1s2-compute: residue facts, own code.  Symbols at lambda of order 5 (F16) and order 7 (F8) read by the same
definition-based regularity routine as vcheck.py (deep rows of T_f(q) M banded Toeplitz), with q now allowed to contain
1+S+S^2 factors.  Checks: (F16) all symbols preserve H(u,v) = sum u_x v_z^4 + u_z v_x^4; orbit of L0 = 325 planes;
(F8) symbols of E, Gm, G3b do not scale the symplectic form; orbit of L0 = 4745 planes."""
import sys, itertools
src = open("vcheck.py").read().split('print("== (A) regularity')[0]
exec(src)
class GF:
    def __init__(self, k, mod):
        self.k, self.mod, self.size = k, mod, 1 << k
        self.mul_t = [[self._mul(a, b) for b in range(self.size)] for a in range(self.size)]
        self.inv_t = [0] + [next(b for b in range(1, self.size) if self.mul_t[a][b] == 1) for a in range(1, self.size)]
    def _mul(self, a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1; a <<= 1
            if a >> self.k: a ^= self.mod
        return r
    def pw(self, a, e):
        e %= (self.size - 1); r = 1
        for _ in range(e): r = self.mul_t[r][a]
        return r
def order(F, a):
    e, x = 1, a
    while x != 1: x = F.mul_t[x][a]; e += 1
    return e
def run(F, lam, tag):
    lp = lambda d: F.pw(lam, d)
    def qe(q):
        s = 0
        for d, c in q.items():
            if c: s ^= lp(d)
        return s
    base = [{0: 1, 1: 1}, {0: 1, -1: 1}, {0: 1, 1: 1, 2: 1}, {0: 1, -1: 1, -2: 1}]
    cands = []
    for ex in itertools.product(range(3), repeat=4):
        q = {0: 1}
        for b, e in zip(base, ex):
            for _ in range(e): q = qmul(q, b)
        cands.append((sum(ex), ex, q))
    cands.sort(key=lambda t: t[0])
    letters = {
        "chi": lambda f: label_from_addr(addr_poly({(0, 0): [0, 1], (1, 1): [0, 1]}, f), f),
        "E": lambda f: label_from_addr(addr_poly({(0, 0): [0], (0, 1): [1], (1, 1): [0]}, f), f),
        "Gm": lambda f: label_from_addr(addr_poly({(0, 0): [0], (0, 1): [1], (1, 0): [1], (1, 1): [0, 1, 2]}, f), f),
        "G3": lambda f: label_from_addr(addr_poly({(0, 0): [0, 1, 2], (1, 1): [0, 1, 2]}, f), f),
        "G3b": lambda f: label_from_addr(addr_poly({(0, 0): [0, 1], (0, 1): [1], (1, 0): [1], (1, 1): [0]}, f), f),
        "R3": lambda f: label_from_addr(addr_poly({(0, 0): [0, 3], (1, 1): [0, 3]}, f), f),
        "D1": LETTERS["D1"][1], "Dcz": LETTERS["Dcz"][1], "CZch": LETTERS["CZch"][1], "CNOT": LETTERS["CNOT"][1],
        "SWAPQ": LETTERS["SWAPQ"][1], "H0": LETTERS["H0"][1], "H1": LETTERS["H1"][1], "HH": LETTERS["HH"][1],
        "P": LETTERS["P"][1]}
    SY = {}
    for nm, fn in letters.items():
        got = None
        for _, ex, q in cands:
            if qe(q) == 0: continue
            ns = []
            for f in (23, 31):
                n = deep_toeplitz(fn(f), q, f)
                if n is None: break
                ns.append(n)
            if len(ns) == 2 and set(ns[0]) == set(ns[1]) and all(np.array_equal(ns[0][d], ns[1][d]) for d in ns[0]):
                got = (ex, q, ns[0]); break
        if got is None: print(f"  [{tag}] {nm}: NOT regular"); continue
        ex, q, n = got
        S = [[0] * 4 for _ in range(4)]
        for d, blk in n.items():
            for i in range(4):
                for k in range(4):
                    if blk[i, k]: S[i][k] ^= lp(d)
        iq = F.inv_t[qe(q)]
        SY[nm] = [[F.mul_t[iq][x] for x in row] for row in S]
        print(f"  [{tag}] {nm}: regular, q exponents (1+S, 1+S^-1, 1+S+S^2, rev) = {ex}")
    return SY
def apply(F, s, v):
    return tuple(_xor([F.mul_t[s[i][k]][v[k]] for k in range(4)]) for i in range(4))
def _xor(xs):
    r = 0
    for x in xs: r ^= x
    return r
def rref(F, rows):
    rows = [list(r) for r in rows]; out = []; c = 0
    for c in range(4):
        p = next((i for i, r in enumerate(rows) if r[c]), None)
        if p is None: continue
        r = rows.pop(p); iv = F.inv_t[r[c]]; r = [F.mul_t[iv][x] for x in r]
        rows = [[x ^ F.mul_t[rr[c]][y] for x, y in zip(rr, r)] for rr in rows]
        out = [[x ^ F.mul_t[o[c]][y] for x, y in zip(o, r)] for o in out]
        out.append(r)
    return tuple(tuple(r) for r in out if any(r))
def orbit(F, SY, L0):
    gens = list(SY.values())
    seen = {L0}; frontier = [L0]
    while frontier:
        nxt = []
        for P in frontier:
            for s in gens:
                Q = rref(F, [apply(F, s, v) for v in P])
                if Q not in seen: seen.add(Q); nxt.append(Q)
        frontier = nxt
    return seen
ok = True
F16 = GF(4, 0b10011)
lam5 = next(a for a in range(2, 16) if order(F16, a) == 5)
print(f"F16: lambda = {lam5} of order 5; lambda^-1 = lambda^4 = Frob^2(lambda): {F16.pw(lam5, -1) == F16.pw(lam5, 4)}")
SY16 = run(F16, lam5, "F16")
fr4 = lambda a: F16.pw(a, 4) if a else 0
def Hh(u, v): return F16.mul_t[u[0]][fr4(v[2])] ^ F16.mul_t[u[1]][fr4(v[3])] ^ F16.mul_t[u[2]][fr4(v[0])] ^ F16.mul_t[u[3]][fr4(v[1])]
BAS = [tuple(int(i == t) for i in range(4)) for t in range(4)]
herm = {nm: all(Hh(apply(F16, s, u), apply(F16, s, v)) == Hh(u, v) for u in BAS for v in BAS) for nm, s in SY16.items()}
print("F16: every symbol preserves the Hermitian form (x -> x^4):", all(herm.values()), herm)
ok &= all(herm.values()) and {"chi", "G3", "G3b", "R3", "E", "Gm"} <= set(SY16)
L0 = rref(F16, [(0, 0, 1, 0), (0, 0, 0, 1)])
orb = orbit(F16, SY16, L0)
print(f"F16: orbit of L0 under the letter symbols: {len(orb)} planes (claim 325)")
ok &= len(orb) == 325
iso = all(Hh(u, v) == 0 for P in orb for u in P for v in P)
print("F16: every orbit plane is totally H-isotropic (a generator of H(3,16)):", iso); ok &= iso
F8 = GF(3, 0b1011)
lam7 = 2
print(f"F8: lambda = {lam7} of order {order(F8, lam7)}")
SY8 = run(F8, lam7, "F8")
def om(u, v): return F8.mul_t[u[0]][v[2]] ^ F8.mul_t[u[1]][v[3]] ^ F8.mul_t[u[2]][v[0]] ^ F8.mul_t[u[3]][v[1]]
for nm in ("chi", "G3", "E", "Gm", "G3b"):
    s = SY8[nm]
    scal = set()
    for u in BAS:
        for v in BAS:
            a, b = om(apply(F8, s, u), apply(F8, s, v)), om(u, v)
            if b: scal.add(F8.mul_t[a][F8.inv_t[b]])
            elif a: scal.add("nonzero-on-zero")
    print(f"F8: symbol of {nm} scales omega by {scal} ({'a scalar' if len(scal) == 1 and 'nonzero-on-zero' not in scal else 'NOT a scalar multiple'})")
L08 = rref(F8, [(0, 0, 1, 0), (0, 0, 0, 1)])
orb8 = orbit(F8, SY8, L08)
print(f"F8: orbit of L0: {len(orb8)} planes (claim 4745 = all planes of F8^4)")
ok &= len(orb8) == 4745
print("ALL PASS" if ok else "SOME FAIL")
