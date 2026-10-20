# meta-spectrum liftsample.py: sampled distance profile (abar_1..abar_m) of the b56-methods "start lift" at any m,
# without listing the frames (needed at m = 6: 4,922,775 frames).
# Model (b56-methods gen.py): n = 2m+2 points, INF = n-1; frames = self-dual codes C (dim m+1, ONE in C);
# frame distance d(C, C') = m + 1 - dim(C cap C').  Lift rule copied from gen.py: weight-2 words {a,b} of C avoiding INF
# give rows e_a - e_b; the other rows are the 0/1 RREF basis (gen.py's rref) of the words avoiding INF and the pair
# points, signed alternately.  delta(C, C') = rank_Q(stacked rows) - m (exact, fraction-free integer elimination).
# Uniform pairs at distance d: U uniform (distance d' from a fixed frame L0 with probability k_d'/N, then uniform at that
# distance), V uniform at distance d from U: W = uniform (m+1-d)-dim subspace of U containing ONE, then
# V = W + span(f_i + sum_j S_ij e_j), S a uniform symmetric d x d matrix over F_2 (all Lagrangians V with V cap U = W).
# usage: python3 -I liftsample.py m pairs_per_distance seed
import sys, random
m = int(sys.argv[1]); NP = int(sys.argv[2]); rng = random.Random(int(sys.argv[3]))
n = 2 * m + 2; INF = n - 1; ONE = (1 << n) - 1; dcol = n - 1

def rref(vecs):  # gen.py
    basis = []
    for v in vecs:
        for p, b in basis:
            if v >> p & 1: v ^= b
        if v:
            p = v.bit_length() - 1
            basis = [(q, (b ^ v) if (b >> p & 1) else b) for q, b in basis]
            basis.append((p, v))
    basis.sort(reverse=True)
    return tuple(b for p, b in basis)
def span(basis):
    out = [0]
    for b in basis: out += [x ^ b for x in out]
    return out
def dot(a, b): return bin(a & b).count('1') & 1

def solve(rows, rhs):  # find x (n bits) with dot(rows[k], x) = rhs[k]; rows independent
    piv = []
    for e, r in zip(rows, rhs):
        for (pe, pr, pb) in piv:
            if e >> pb & 1: e ^= pe; r ^= pr
        assert e, 'dependent'
        pb = e.bit_length() - 1
        piv = [((qe ^ e, qr ^ r, qb) if (qe >> pb & 1) else (qe, qr, qb)) for (qe, qr, qb) in piv]
        piv.append((e, r, pb))
    x = 0
    for (pe, pr, pb) in piv:
        if pr: x |= 1 << pb
    return x

def rand_at(U, d):
    """U: code basis (tuple, rref, contains ONE in span). Return a uniform frame at distance d from U."""
    if d == 0: return U
    Uw = [w for w in span(U)]
    # W: uniform subspace of U containing ONE of dimension m+1-d
    W = [ONE]
    while len(W) < m + 1 - d:
        v = rng.choice(Uw)
        if len(rref(W + [v])) > len(W): W.append(v)
    # complete W to U: e_1..e_d
    E = []
    while len(W) + len(E) < m + 1:
        v = rng.choice(Uw)
        if len(rref(W + E + [v])) > len(W) + len(E): E.append(v)
    # g_j in W^perp with <e_i, g_j> = delta_ij (W contains ONE, so g_j has even weight)
    F = []
    for j in range(d):
        g = solve(W + E, [0] * len(W) + [1 if i == j else 0 for i in range(d)])
        for k in range(j):
            if dot(g, F[k]): g ^= E[k]
        F.append(g)
    for i in range(d):
        for j in range(d): assert dot(E[i], F[j]) == (i == j) and dot(F[i], F[j]) == 0
    S = [[0] * d for _ in range(d)]
    for i in range(d):
        for j in range(i, d):
            S[i][j] = S[j][i] = rng.randrange(2)
    V = list(W)
    for i in range(d):
        v = F[i]
        for j in range(d):
            if S[i][j]: v ^= E[j]
        V.append(v)
    V = rref(V); assert len(V) == m + 1
    return V

def lift_rows(C):  # gen.py rule
    sp = span(C)
    w2 = [w for w in sp if bin(w).count('1') == 2]
    rows = []; cov = 0
    for w in w2: cov |= w
    for w in w2:
        if w >> INF & 1: continue
        a, b = [k for k in range(n) if w >> k & 1]
        r = [0] * dcol; r[a] = 1; r[b] = -1; rows.append(r)
    if cov != ONE:
        restL = rref([w for w in sp if not (w >> INF & 1) and not (w & cov)])
        for w in restL:
            pts = [k for k in range(dcol) if w >> k & 1]
            r = [0] * dcol
            for t, p in enumerate(pts): r[p] = 1 if t % 2 == 0 else -1
            rows.append(r)
    assert len(rows) == m, (len(rows), m)
    return rows

def rank_q(A):  # fraction-free Bareiss, exact
    A = [row[:] for row in A]; R = len(A); Cc = len(A[0]); rk = 0; prev = 1
    for c in range(Cc):
        if rk == R: break
        p = next((r for r in range(rk, R) if A[r][c]), None)
        if p is None: continue
        A[rk], A[p] = A[p], A[rk]
        for r in range(rk + 1, R):
            for k in range(c + 1, Cc): A[r][k] = (A[r][k] * A[rk][c] - A[r][c] * A[rk][k]) // prev
            A[r][c] = 0
        prev = A[rk][c]; rk += 1
    return rk
def rank2(A):
    rows = [sum((x & 1) << k for k, x in enumerate(r)) for r in A]
    return len(rref(rows))

def qbinom(a, b):
    num = den = 1
    for i in range(b): num *= (2 ** (a - i) - 1); den *= (2 ** (i + 1) - 1)
    return num // den
ks = [2 ** (i * (i + 1) // 2) * qbinom(m, i) for i in range(m + 1)]; N = sum(ks)
L0 = rref([ONE] + [(1 << (2 * i)) | (1 << (2 * i + 1)) for i in range(m + 1)])
cache = {}
def lr(C):
    if C not in cache:
        if len(cache) > 200000: cache.clear()
        cache[C] = lift_rows(C)
    return cache[C]
print(f'm={m}: frames {N}, sphere sizes {ks}; {NP} uniform pairs per distance, seed {sys.argv[3]}', flush=True)
for d in range(1, m + 1):
    s = s2 = 0; hist = {}
    for t in range(NP):
        dd = rng.choices(range(m + 1), weights=ks)[0]
        U = rand_at(L0, dd); V = rand_at(U, d)
        A = lr(U) + lr(V)
        assert rank2(A) - m == d
        q = rank_q(A) - m; s += q; s2 += q * q; hist[q] = hist.get(q, 0) + 1
    mu = s / NP; se = ((s2 / NP - mu * mu) / (NP - 1)) ** 0.5
    print(f'd={d} pairs {NP} sumdelta {s} mean {mu:.6f} stderr {se:.6f} hist {sorted(hist.items())}', flush=True)
