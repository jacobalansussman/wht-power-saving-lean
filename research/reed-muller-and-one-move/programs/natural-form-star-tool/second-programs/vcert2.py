"""(construction part copied from vstar2.py) (rm2-exact, independent of xbase / gen_tables / scanrc / starrec*): an independent re-derivation of the
exclusions made by starrec2 at block 5 (one move per nonzero direction; array l owed T_l = mult(alpha^(2 k_l)),
alpha a root of x^5 + x^2 + 1).  It does not read the C program's certificates; it re-proves each exclusion itself.

Own constructions: GF(32) by shift-and-reduce; the coordinate basis (argument) checked self-dual; directions as field
elements x = sum_j z_j f_j; forms decoded from their 15 coefficients by definition; parity of Q at T_k computed as
sum_j Q(alpha^k f_j); affine flats built as cosets of F_2-spans of field elements.
Facts (proved; see NOTE.md): T10 rank E[D_Q,D_Q] = |D_Q| - r_Q;  STAR (g-2) rank E[Psi] <= sum_B rank E[B] - 2|Psi|
for the g flats B of one dimension more through a flat A (Psi = D_F n A, E[D_F,D_F] nonsingular);  CLQ E[Y,Y] = 0
=> |Y| <= W;  MAT rank E[D,D] <= largest matching avoiding known zero cells.
usage: vstar2.py <basis f1,..,f5> <maxd> <rounds> < lines "k1,k2,..."     prints one verdict per line.
"""
import sys
import numpy as np

POLY = 0b100101
def gmul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a & 32: a ^= POLY
    return r
def tr(a):
    s, y = 0, a
    for _ in range(5): s ^= y; y = gmul(y, y)
    return s
basis = [int(x) for x in sys.argv[1].split(",")]

for i in range(5):
    for j in range(5): assert tr(gmul(basis[i], basis[j])) == (i == j)
pw = [1]
for _ in range(30): pw.append(gmul(pw[-1], 2))
assert gmul(pw[-1], 2) == 1 and len(set(pw)) == 31
z_of = {}                                 # field element -> direction index z (1..31)
for z in range(1, 32):
    x = 0
    for j in range(5):
        if z >> j & 1: x ^= basis[j]
    z_of[x] = z
def emask(elems):                          # set of nonzero field elements -> direction bit mask (bit z-1)
    m = 0
    for x in elems:
        if x: m |= 1 << (z_of[x] - 1)
    return m
PAIRS = [(j, k) for j in range(5) for k in range(j, 5)]
NF = 1 << 15
supp = np.zeros(NF, dtype=np.int64); size = np.zeros(NF, dtype=np.int64); par = np.zeros((NF, 31), dtype=np.int64)
inv = {z: x for x, z in z_of.items()}
for q in range(NF):
    S = set()
    for x, z in z_of.items():
        v = 0
        for p, (j, k) in enumerate(PAIRS):
            if q >> p & 1: v ^= (z >> j & 1) & (z >> k & 1)
        if v: S.add(x)
    supp[q] = emask(S); size[q] = len(S)
    for k in range(31): par[q, k] = sum(1 for f in basis if gmul(pw[k], f) in S) % 2
# affine flats of F_32 by dimension
def span(gs):
    V = {0}
    for g in gs: V |= {v ^ g for v in V}
    return frozenset(V)
subsp = {0: {frozenset([0])}}
for d in range(1, 5):
    subsp[d] = set()
    for V in subsp[d - 1]:
        for g in range(1, 32):
            if g not in V: subsp[d].add(span(list(V) + [g]))
flats = {}
for d in range(5):
    L = set()
    for V in subsp[d]:
        for a in range(32): L.add(frozenset(a ^ v for v in V))
    flats[d] = sorted(L, key=lambda f: sorted(f))
assert [len(flats[d]) for d in range(5)] == [32, 496, 1240, 620, 62]
fmask = {d: [emask(A) for A in flats[d]] for d in range(5)}
upl = {}
for d in range(4):
    idx = {B: i for i, B in enumerate(flats[d + 1])}
    upl[d] = [[idx[B] for B in flats[d + 1] if A <= B] for A in flats[d]]
    assert all(len(u) == (1 << (5 - d)) - 1 for u in upl[d])
popc = lambda x: bin(x).count("1")


def maxmatch(D, Z):
    rows = [x for x in range(1, 32) if D >> (x - 1) & 1]
    match = {}
    def aug(x, seen):
        for y in rows:
            if (x, y) in Z or y in seen: continue
            seen.add(y)
            if y not in match or aug(match[y], seen): match[y] = x; return True
        return False
    return sum(1 for x in rows if aug(x, set()))

_cache = {}
def round1(F, s, key):
    """round-1 star recursion for the even form F: returns (contradiction?, set of zero-block masks)"""
    if (key, F) in _cache: return _cache[(key, F)]
    DF = int(supp[F]); dirs = [z for z in range(1, 32) if DF >> (z - 1) & 1]; n = len(dirs)
    bits = ((supp[:, None] >> (np.array(dirs) - 1)[None, :]) & 1)
    cidx = bits @ (1 << np.arange(n))
    best = np.full(1 << n, 10 ** 6, dtype=np.int64)
    sel = np.arange(NF) != F; sel[0] = False
    np.minimum.at(best, cidx[sel], s[sel])
    for b in range(n):
        v = best.reshape(-1, 2, 1 << b)
        v[:, 0, :] = np.minimum(v[:, 0, :], v[:, 1, :])
    def comp(msk): return sum(1 << i for i, z in enumerate(dirs) if msk >> (z - 1) & 1)
    u = {}; zb = set(); contra = False
    for d in range(5):
        u[d] = [min(popc(A & DF), int(best[comp(A & DF)])) for A in fmask[d]]
        for i, A in enumerate(fmask[d]):
            if A & DF and u[d][i] == 0: zb.add(A & DF)
    for d in range(3, -1, -1):
        g = (1 << (5 - d)) - 1
        for i, A in enumerate(fmask[d]):
            S = A & DF; nS = popc(S)
            if nS == 0: continue
            tot = sum(u[d + 1][j] for j in upl[d][i]); val = tot - 2 * nS
            if val < 0: contra = True; break
            u[d][i] = min(u[d][i], val // (g - 2), min(u[d + 1][j] for j in upl[d][i]))
            if u[d][i] == 0: zb.add(S)
        if contra: break
    _cache[(key, F)] = (contra, zb)
    return contra, zb

nok = nbad = nskip = 0
for line in sys.stdin:
    p = line.split()
    if not p or p[0] != "S" or len(p) < 3 or p[2] not in ("C", "M", "Z"): continue
    ks = [int(x) for x in p[1].split(",")]; key = p[1]
    cnt = np.zeros(31, dtype=np.int64)
    for k in ks: cnt[k % 31] += 1
    W = int(cnt.sum()); r = par @ cnt; s = size - r
    assert (s[1:] >= 0).all()
    typ, rnd = p[2], int(p[3])
    if rnd != 1: nskip += 1; continue
    ok = True
    if typ == "Z":
        F = int(p[4]); assert r[F] == 0
        contra, _ = round1(F, s, key); ok = contra
    else:
        if typ == "C":
            q = int(p[4]); Y = int(p[6].split("=")[1]); cells = p[7:]
        else:
            Q = int(p[4]); mmc = int(p[5]); sQ = int(p[6]); cells = p[7:]
        Z = set()
        for c in cells:
            xy, src = c.split(":"); x, y = map(int, xy.split(",")); f = list(map(int, src.split("/")))
            F, S, rr = f[0], f[1], f[2]; dd = f[3] if len(f) > 3 else 0; ii = f[4] if len(f) > 4 else 0
            if rr != 1 or r[F] != 0: ok = False; break
            if not (S >> (x - 1) & 1 and S >> (y - 1) & 1): ok = False; break
            if dd == -1:
                # quadratic star block: Psi = D_F n D_G1 n D_G2, sets D_F n D_G1, D_F n D_G2, D_F n D_(F+G1+G2)
                G1, G2 = ii // 65536, ii % 65536; G3 = F ^ G1 ^ G2; DF = int(supp[F])
                P = DF & int(supp[G1]) & int(supp[G2])
                if P != S or G1 in (0, F) or G2 in (0, F): ok = False; break
                sig = 0
                for G in (G1, G2, G3):
                    Sg = DF & int(supp[G]); b = popc(Sg)
                    if G not in (0, F): b = min(b, int(s[G]))
                    sig += b
                if sig - 2 * popc(P) > 0: ok = False; break
            else:
                contra, zb = round1(F, s, key)
                if S not in zb: ok = False; break
            Z.add((x, y))
        if ok and typ == "C":
            ys = [x for x in range(1, 32) if Y >> (x - 1) & 1]
            ok = len(ys) > W and all((x, y) in Z for x in ys for y in ys)
        elif ok:
            assert s[Q] == sQ
            ok = maxmatch(int(supp[Q]), Z) < s[Q]
    if ok: nok += 1
    else: nbad += 1; print("NOT VERIFIED", line[:200].strip())
    sys.stdout.flush()
print("round-1 certificates verified %d, not verified %d, skipped (round >= 2) %d" % (nok, nbad, nskip))
