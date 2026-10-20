"""vf.py (v-rm2-exact): skeptic's own GF(32) / forms / canonical-target parities, written from scratch.
Field: GF(2)[x]/(x^5+x^2+1) by carry-less multiplication (no log tables). Directions z in 1..31, bit j = coordinate
on basis element f_j. Forms Q(z) = sum_{j<=k} q_jk z_j z_k, index bits in lexicographic (j,k) j<=k order."""
import numpy as np, struct, sys
POLY = 0b100101; m = 5; n = 31
def gmul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a & 32: a ^= POLY
    return r
def gpow(a, e):
    r = 1
    for _ in range(e % n if a else e): r = gmul(r, a)
    return r
def tr(x):
    s, y = 0, x
    for _ in range(m): s ^= y; y = gmul(y, y)
    assert s in (0, 1); return s
ALPHA = 2
assert len({gpow(ALPHA, e) for e in range(31)}) == 31  # primitive
BASIS = [3, 5, 12, 17, 26]
for i in range(m):
    for j in range(m):
        assert tr(gmul(BASIS[i], BASIS[j])) == (i == j), "basis not self-dual"
def elem(z): # direction -> field element
    x = 0
    for j in range(m):
        if z >> j & 1: x ^= BASIS[j]
    return x
def coord(x): return sum(tr(gmul(x, BASIS[j])) << j for j in range(m))
for z in range(32): assert coord(elem(z)) == z
PAIRS = [(j, k) for j in range(m) for k in range(j, m)]
NF = 1 << len(PAIRS)
def Qval(q, z):
    s = 0
    for p, (j, k) in enumerate(PAIRS):
        if q >> p & 1: s ^= (z >> j) & (z >> k) & 1
    return s
# support masks, bit z-1 <-> direction z
VALS = np.zeros((NF, 31), dtype=np.uint8)
Z = np.arange(1, 32)
for p, (j, k) in enumerate(PAIRS):
    mon = ((Z >> j) & (Z >> k) & 1).astype(np.uint8)
    qs = (np.arange(NF) >> p) & 1
    VALS ^= (qs[:, None].astype(np.uint8) * mon[None, :])
MASK = (VALS.astype(np.int64) << np.arange(31)[None, :]).sum(axis=1)
WT = VALS.sum(axis=1).astype(np.int64)
def Tmat(k):
    beta = gpow(ALPHA, 2 * k)
    return [[tr(gmul(beta, gmul(BASIS[i], BASIS[j]))) for j in range(m)] for i in range(m)]
# parity <Q, T_k> by pairing with the multiplication matrix
QB = ((np.arange(NF)[:, None] >> np.arange(len(PAIRS))[None, :]) & 1).astype(np.int64)
PAR = np.zeros((NF, 31), dtype=np.uint8)
for k in range(31):
    T = Tmat(k)
    for i in range(m):
        for j in range(m): assert T[i][j] == T[j][i]
    tv = np.array([T[j][kk] for (j, kk) in PAIRS], dtype=np.int64)
    PAR[:, k] = (QB @ tv) & 1
def check_frames():
    # second route: frame g_j = alpha^k f_j, parity = sum_j Q(coord(g_j))
    rng = np.random.default_rng(1)
    for q in list(rng.integers(0, NF, 400)) + [1, 2, 17170]:
        for k in range(31):
            g = [coord(gmul(gpow(ALPHA, k), f)) for f in BASIS]
            s = 0
            for z in g: s ^= Qval(int(q), z)
            assert s == PAR[q, k], (q, k)
    # and Tr(b alpha^k) with b = linear part: Q(x) restricted... check F3 via linear-part classes
    return True
if __name__ == "__main__":
    check_frames(); print("frames OK")
    # compare with rm2-exact tables5.bin
    fn = sys.argv[1]
    with open(fn, 'rb') as f:
        nf = struct.unpack('<i', f.read(4))[0]
        mask = np.frombuffer(f.read(4 * nf), dtype='<u4').astype(np.int64)
        wt = np.frombuffer(f.read(4 * nf), dtype='<i4').astype(np.int64)
        vid = np.frombuffer(f.read(4 * nf), dtype='<i4').astype(np.int64)
        nv = struct.unpack('<i', f.read(4))[0]
        pk = np.frombuffer(f.read(4 * nv), dtype='<u4').astype(np.int64)
    assert nf == NF
    print("mask equal:", np.array_equal(mask, MASK), "wt equal:", np.array_equal(wt, WT))
    mypat = (PAR.astype(np.int64) << np.arange(31)[None, :]).sum(axis=1)
    print("parity patterns equal:", np.array_equal(pk[vid], mypat), "classes", len(set(mypat.tolist())))
    # weights of the even class
    ev = np.nonzero(mypat == 0)[0]
    print("even forms", len(ev), "weights", sorted(set(WT[ev].tolist())), [int((WT[ev] == w).sum()) for w in sorted(set(WT[ev].tolist()))])
