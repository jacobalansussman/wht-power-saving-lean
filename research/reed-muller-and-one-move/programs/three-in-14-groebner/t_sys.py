"""t_sys.py: checks of gb2.system (the equation generator), exact arithmetic.
(1) For random rational A, B, H and a random order, every polynomial equals the matching entry of
    A_D (I + H_D)^-1 B_D^T - diag(P) computed by Gaussian elimination (the closed form of rm2-solve's inst.py).
(2) For plain controls and random orders, the plain program (XOR chains: a = b = e_owner, H = 2 between moves of one
    array in time order) makes every polynomial vanish."""
import json, random, re, itertools
from fractions import Fraction as Fr
from f2q import Field, Config, popc, all_frames
from gb2 import supports_in, orders, Ring, system

F5 = Field(5); certs = json.load(open("certs_min.json")); rng = random.Random(5)


def inv(Mx):
    n = len(Mx); A = [row[:] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(Mx)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[p] = A[p], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


def evalp(p, val):
    t = Fr(0)
    for mono, c in p.items():
        x = Fr(c)
        for v in mono: x *= val[v]
        t += x
    return t


bad1 = 0; n1 = 0
for c in [x for x in certs if x["tag"] == "mine-b5"][:6]:
    cfg = Config(F5, c["Ts"], c["moves"]); U = [i for i in range(cfg.M) if c["cert"]["U"] >> i & 1]
    sup = supports_in(cfg, U); R = Ring(cfg.W, U); o = rng.choice(list(orders(cfg, U)))
    val = [Fr(rng.randint(-5, 5), rng.randint(1, 4)) for _ in R.names]
    polys = system(cfg, R, o, sup); k = 0
    pos = {i: t for t, i in enumerate(o)}
    for d, p in sorted(sup.items()):
        D = sorted((i for i in o if d >> i & 1), key=lambda i: pos[i])
        Hm = [[Fr(int(a == b)) + (val[R.idx[f"h{a}_{b}"]] if pos[a] < pos[b] else 0) for b in D] for a in D]
        N = inv(Hm)
        for l in range(cfg.W):
            for l2 in range(cfg.W):
                X = sum(val[R.idx[f"a{l}_{D[s]}"]] * N[s][t] * val[R.idx[f"b{l2}_{D[t]}"]]
                        for s in range(len(D)) for t in range(len(D)))
                X -= int(l == l2 and (p >> l & 1))
                n1 += 1
                if evalp(polys[k], val) != X: bad1 += 1
                k += 1
print(f"(1) closed form vs generated polynomials: {n1} entries compared, {bad1} mismatches")

bad2 = 0; n2 = 0
for c in [x for x in certs if x["tag"] == "mine-b5-plain"][:8]:
    cfg = Config(F5, c["Ts"], c["moves"])
    frs = [all_frames(F5, k) for k in c["ks"]]
    for ch in itertools.product(*frs):
        if sorted(z for f in ch for z in f) == sorted(cfg.moves): break
    owner = [None] * cfg.M; used = set()
    for l, f in enumerate(ch):
        for z in f:
            i = next(i for i in range(cfg.M) if cfg.moves[i] == z and i not in used); used.add(i); owner[i] = l
    U = rng.sample(range(cfg.M), 6); sup = supports_in(cfg, U); R = Ring(cfg.W, U)
    for o in rng.sample(list(orders(cfg, U)), 3):
        pos = {i: t for t, i in enumerate(o)}; val = [None] * len(R.names)
        for nm, k in R.idx.items():
            if nm[0] in "ab":
                l, i = map(int, nm[1:].split("_")); val[k] = Fr(int(owner[i] == l))
            else:
                i, j = map(int, nm[1:].split("_")); val[k] = Fr(2 if owner[i] == owner[j] else 0)
        for p in system(cfg, R, o, sup):
            n2 += 1; bad2 += evalp(p, val) != 0
print(f"(2) plain programs substituted: {n2} polynomials, {bad2} nonzero")
