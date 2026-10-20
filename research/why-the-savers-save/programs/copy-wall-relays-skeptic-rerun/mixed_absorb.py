"""v-R1-recheck: does the '>= 17 descents' step hold?  Mixed relay decoders: relays to columns in D (descended
registers) are exempt from closure; all other relays respect the transitive closure.  Absorb ONLY the
closure-respecting relays (c^S = B[S] - sum_{S relays T, T not in D} B[S,T] c^T, layer order) and check that the
result is a kfree-type decoder: same row space, diagonal 1, unknown nonzeros only in D columns (earlier layers).
If so, #descended columns >= the kfree relayed-column bound.  python3 -I mixed_absorb.py [trials]"""
import random, sys, itertools
sys.dont_write_bytecode = True
def rank_p(rows, p):
    rows = [list(r) for r in rows]; rk = 0; n = len(rows[0]) if rows else 0
    for c in range(n):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] % p), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], p - 2, p); rows[rk] = [x * inv % p for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % p:
                f = rows[i][c]; rows[i] = [(a - f * b) % p for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk
rng = random.Random(7)
trials = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
stats = [0, 0, 0, 0, 0]
for t in range(trials):
    v = rng.randrange(5, 13); p = rng.choice([2, 3, 5]); L = rng.choice([2, 3, 4, v])
    dens = rng.choice([0.3, 0.5, 0.7])
    K = [[False] * v for _ in range(v)]
    for a, b in itertools.combinations(range(v), 2):
        if rng.random() < dens: K[a][b] = K[b][a] = True
    lay = [rng.randrange(L) for _ in range(v)]
    D = {x for x in range(v) if rng.random() < 0.3}
    order = sorted(range(v), key=lambda x: lay[x])
    B = [[0] * v for _ in range(v)]; F = {}
    for S in order:
        B[S][S] = 1; F[S] = set()
        for T in range(v):
            if T == S: continue
            if K[S][T]:
                if rng.random() < 0.5: B[S][T] = rng.randrange(1, p); F[S].add(T)
            elif lay[T] < lay[S] and rng.random() < 0.6:
                if T in D:   # descended register: footprint dropped; G_U contributes labels known to S and to T
                    B[S][T] = rng.randrange(1, p)
                    F[S] |= {x for x in range(v) if K[S][x] and K[T][x] and rng.random() < 0.3}
                elif F[T] <= {x for x in range(v) if K[S][x]}:
                    B[S][T] = rng.randrange(1, p); F[S] |= F[T]
    C = {}
    for S in order:
        c = list(B[S])
        for T in range(v):
            if T != S and not K[S][T] and B[S][T] and T not in D:
                f = B[S][T]; c = [(a - f * b) % p for a, b in zip(c, C[T])]
        C[S] = c
    diag = all(C[S][S] % p == 1 for S in range(v))
    unk_ok = all(C[S][T] % p == 0 for S in range(v) for T in range(v) if T != S and not K[S][T] and T not in D)
    lay_ok = all(lay[T] < lay[S] for S in range(v) for T in D if T != S and not K[S][T] and C[S][T] % p)
    span = rank_p([B[i] for i in range(v)] + [C[i] for i in range(v)], p) == rank_p(B, p)
    stats[0] += 1; stats[1] += diag; stats[2] += unk_ok; stats[3] += lay_ok; stats[4] += span
print('instances %d  diag=1 %d  unknown nonzeros only in D columns %d  D-relays respect layers %d  rows in rowspace(B) %d' % tuple(stats))
