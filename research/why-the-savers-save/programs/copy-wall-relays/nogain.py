"""nogain.py (d-relay) -- the relay no-gain lemma, checked mechanically on random instances.

Physical statement (frame game: equal-frame additions H10, frames climb H9, sources at their labels / sandwich H7,
content leaves a frame only through a paid copy child).  Invariant: every register at frame F holds
    (combination of copy-dropped contents) + (sources x_T with a_T in F) + (y-content).
Base: a source at <a_T>; a target at 0 holding y.  Step: an addition at frame F joins two registers at F, both
satisfying the invariant for F; climbing F -> F' (F inside F') keeps it; a drop is a paid copy.  A relay (S reads
T's register at G) is such an addition, so it changes nothing.  At the end target S sits inside a_S^perp holding
y_S + x_S, hence  x_S in span(dropped contents) + span{x_T : a_T in a_S^perp} = the STATIC fit of the dropped
totals.  So relays (any order, any number of rounds, second-round reads included) never lower the static minrank
or the star floor: the dropped totals already fit without them.  Any field, any unit with non-isotropic labels.

Algebraic statement checked here (the SAT model of relaysat.py, transitive closure):
B[S,S] = 1; B[S,T] != 0 with S not knowing T only if layer(T) < layer(S) ('S relays T');
F(T) = {T' known to T : B[T,T'] != 0} u U_{T relays T''} F(T'');  S relays T  =>  F(T) inside known(S).
Then c^T = B[T] - sum_{T relays T''} B[T,T''] c^{T''} (layer order) lies in rowspace(B), c^T[T] = 1, and
supp(c^T) - {T} lies in F(T), inside known(T): a STATIC fit inside rowspace(B).  Proof: induction on the layer;
c^{T''}[T] = 0 since T is not in known(T) (non-isotropic), and c^{T'''}[T''] = 0 for T'' not known to T because
supp(c^{T'''}) - {T'''} lies in F(T''') inside known(T).
With only the one-level clause (S relays T => B[T,T'] = 0 for T' known to T, not to S) the lemma can fail when
L >= 3; this script shows instances.  usage: python3 -I nogain.py [trials]"""
import random, sys, itertools
sys.dont_write_bytecode = True

def rank_p(rows, p):
    rows = [list(r) for r in rows]; rk = 0; n = len(rows[0]) if rows else 0
    for c in range(n):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] % p), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], p - 2, p)
        rows[rk] = [x * inv % p for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % p:
                f = rows[i][c]; rows[i] = [(a - f * b) % p for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk

def random_knows(v, dens, rng):
    K = [[False] * v for _ in range(v)]
    for a, b in itertools.combinations(range(v), 2):
        if rng.random() < dens: K[a][b] = K[b][a] = True
    return K

def random_B(K, lay, p, rng, mode, prelay=0.6, pknown=0.5):
    """rows in layer order; random entries on allowed cells respecting the closure mode"""
    v = len(K); order = sorted(range(v), key=lambda t: lay[t]); B = [[0] * v for _ in range(v)]; F = {}
    for S in order:
        B[S][S] = 1
        rel = []
        for T in range(v):
            if T == S: continue
            if K[S][T]:
                if rng.random() < pknown: B[S][T] = rng.randrange(1, p)
            elif lay[T] < lay[S] and rng.random() < prelay:
                if mode == 'trans': ok = F[T] <= {x for x in range(v) if K[S][x]}
                else: ok = all(not B[T][x] for x in range(v) if x != T and K[T][x] and not K[S][x])
                if ok: B[S][T] = rng.randrange(1, p); rel.append(T)
        F[S] = {x for x in range(v) if x != S and K[S][x] and B[S][x]}
        for T in rel: F[S] |= F[T]
    return B

def c_vectors(B, K, lay, p):
    v = len(B); C = {}
    for T in sorted(range(v), key=lambda t: lay[t]):
        c = list(B[T])
        for T2 in range(v):
            if T2 != T and not K[T][T2] and B[T][T2]:
                f = B[T][T2]; c = [(a - f * b) % p for a, b in zip(c, C[T2])]
        C[T] = c
    return C

def static_ok(C, K, p):
    v = len(C)
    return all(C[S][S] % p == 1 and all(C[S][T] % p == 0 for T in range(v) if T != S and not K[S][T]) for S in range(v))

def main(trials):
    rng = random.Random(20261013)
    stats = {}
    for t in range(trials):
        v = rng.randrange(5, 13); p = rng.choice([2, 3, 5]); L = rng.choice([2, 3, 4, v])
        K = random_knows(v, rng.choice([0.3, 0.5, 0.7]), rng); lay = [rng.randrange(L) for _ in range(v)]
        for mode in ('trans', 'one'):
            B = random_B(K, lay, p, rng, mode)
            used = sum(1 for S in range(v) for T in range(v) if S != T and not K[S][T] and B[S][T])
            C = c_vectors(B, K, lay, p)
            ok = static_ok(C, K, p)
            inspan = rank_p([B[i] for i in range(v)] + [C[i] for i in range(v)], p) == rank_p(B, p)
            key = (mode, 'relays used' if used else 'no relay')
            s = stats.setdefault(key, [0, 0, 0]); s[0] += 1; s[1] += ok; s[2] += inspan
    for key in sorted(stats):
        n, ok, sp = stats[key]
        print('closure=%-5s %-12s instances %5d   c^T static fit %5d   c^T in rowspace(B) %5d' % (key + (n, ok, sp)))

if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 3000)
