"""v-rm-theory: own check of rm-theory 5b 'block 5' (mt_scan5z / t_tightcode claims) on MY OWN random draws.
Block 5, one move per direction (31 moves), W distinct canonical Singer targets T_k = sigma^k sigma^kT (own GF(32),
own self-dual basis). For each draw: tight forms? (8.8 criterion: some nonzero b with #{l: Tr(b alpha^k_l)=1} = 8),
deletion, cliques (C), zero block (Z). Claim to test: every set that has a tight form is excluded (by C or Z);
the criterion agrees with the direct tight-form count.
usage: vk_b5.py W ndraws seed
"""
import itertools, random, sys
import numpy as np

W, ND, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
m = 5; q = 32; PRIM = 0b100101
def mul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a & q: a ^= PRIM
    return r
def power(a, e):
    r = 1
    for _ in range(e): r = mul(r, a)
    return r
def tr(a):
    s = 0; x = a
    for _ in range(m): s ^= x; x = mul(x, x)
    return s
TR = [tr(a) for a in range(q)]
f = next(c for c in itertools.combinations(range(1, q), m)
         if all(TR[mul(a, b)] == (1 if a == b else 0) for a in c for b in c))
def coord(a): return [TR[mul(a, x)] for x in f]
def mulmat(beta):
    M = np.zeros((m, m), dtype=np.int64)
    for j in range(m): M[:, j] = coord(mul(beta, f[j]))
    return M
sig = mulmat(2)
Tk = []; S = np.eye(m, dtype=np.int64)
for k in range(q - 1):
    Tk.append((S @ S.T) % 2); S = (sig @ S) % 2
pairs = [(i, j) for i in range(m) for j in range(i, m)]
nf = 1 << len(pairs)
Qbits = np.array([[(Q >> t) & 1 for t in range(len(pairs))] for Q in range(nf)], dtype=np.int64)
pts = list(range(1, q))
X = np.array([[(z >> j) & 1 for j in range(m)] for z in pts], dtype=np.int64)
mon = np.array([X[:, i] * X[:, j] for (i, j) in pairs]).T
Qtab = ((Qbits @ mon.T) % 2).astype(np.int8)         # nf x 31 ; one move per direction
D = Qtab.sum(1).astype(np.int16)
Pk = np.array([((Qbits @ np.array([T[i, j] for (i, j) in pairs])) % 2) for T in Tk]).astype(np.int16)  # 31 x nf

def max_clique(cands, adj):
    best = 0
    def rec(s, Pset):
        nonlocal best
        if s + len(Pset) <= best: return
        if not Pset: best = max(best, s); return
        while Pset:
            if s + len(Pset) <= best: return
            v = Pset[0]
            rec(s + 1, [u for u in Pset[1:] if adj[v][u]])
            Pset = Pset[1:]
    rec(0, list(cands))
    return best

rng = random.Random(seed)
stats = {}
agree = 0
for d in range(ND):
    ks = sorted(rng.sample(range(q - 1), W))
    r = Pk[ks].sum(0)
    # 8.8 criterion: some nonzero b with #{l : Tr(b * alpha^(k_l)) = 1} == 8   (copy-k parity = Tr(b alpha^k), checked in vk_field)
    crit = any(sum(TR[mul(b, power(2, k))] for k in ks) == 8 for b in range(1, q))
    if (r > D).any(): v = 'deletion'; tight = None
    else:
        tight = np.where((D == r) & (r >= 1))[0]
        if len(tight) == 0: v = 'no tight form'
        else:
            live = [z for z in range(31) if Qtab[tight, z].any()]
            adj = {a: {b: bool((Qtab[tight, a] & Qtab[tight, b]).any()) for b in live} for a in live}
            c = max_clique(live, adj)
            if c > W: v = 'C'
            else:
                v = 'passes C,Z'
                for F in np.where(D + r < 2 * c)[0]:
                    sub = [z for z in live if Qtab[F, z]]
                    if sub and D[F] + r[F] < 2 * max_clique(sub, adj): v = 'Z'; break
    has_tight = None if tight is None else (len(tight) > 0)
    if has_tight is not None: agree += (has_tight == crit)
    stats[v] = stats.get(v, 0) + 1
print(f"block 5, W={W}, {ND} draws (seed {seed}): {stats}; criterion 8.8 agrees with direct tight count in {agree} of "
      f"{ND - stats.get('deletion', 0)} deletion-passing draws")
