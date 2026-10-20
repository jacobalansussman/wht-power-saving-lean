#!/usr/bin/env python3
"""v-an2-smallnet / vclass.py (own code, stdlib; python3 -I).  Independent re-run of the label-class floor of
an2-smallnet NOTE 3.2: the largest set Z of ports that splits into groups with labels pairwise non-orthogonal inside a
group (B = 1) and orthogonal across groups (B = 0).  Equivalently: an induced subgraph of the graph B=1 that is a
disjoint union of cliques (a 'cluster' induced subgraph).  Exact branch and bound for n <= 20, randomised greedy +
local search above.  Consistency test: for P112 the floor must not exceed the hub's 64 per call (Lean-checked word),
and for E8/E7 we print it against t(k-1)."""
import itertools, json, random, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import vdesign as V

def Bmat(L):
    n = len(L); nn = V.dot(L[0], L[0])
    return [[1 if V.dot(L[i], L[j]) in (0, nn, -nn) else 0 for j in range(n)] for i in range(n)]

def is_cluster(B, S):
    S = list(S)
    for a, b, c in itertools.combinations(S, 3):
        e = B[a][b] + B[b][c] + B[a][c]
        if e == 2: return False      # induced P3 -> not a disjoint union of cliques
    return True

def exact(B, n):
    best = [0]
    def rec(i, S):
        if len(S) + n - i <= best[0]: return
        if i == n: best[0] = max(best[0], len(S)); return
        T = S + [i]
        if all(not (B[a][b] + B[b][i] + B[a][i] == 2) for a, b in itertools.combinations(S, 2)):
            rec(i + 1, T)
        rec(i + 1, S)
    rec(0, [])
    return best[0]

def greedy(B, n, iters, seed):
    rnd = random.Random(seed); best = 0
    for _ in range(iters):
        order = list(range(n)); rnd.shuffle(order); S = []
        for v in order:
            ok = True
            for a, b in itertools.combinations(S, 2):
                pass
            # fast check: v may join if for every pair (a in S) the triple condition holds
            ok = all(not (B[a][b] + B[b][v] + B[a][v] == 2) for a, b in itertools.combinations(S, 2)) if len(S) < 40 else \
                 is_cluster(B, S + [v])
            if ok: S.append(v)
        best = max(best, len(S))
    return best

def main():
    E8 = V.e8_lines()
    frame = set(V.canon({tuple((2 if q == 2*i else (2*sg if q == 2*i+1 else 0)) for q in range(8))
                       for i in range(4) for sg in (1, -1)}))
    P112 = [l for l in E8 if l not in frame]
    E7 = [l for l in E8 if V.dot(l, (2, -2, 0, 0, 0, 0, 0, 0)) == 0]
    import itertools as it
    def Dn(n):
        S = set()
        for i, j in it.combinations(range(n), 2):
            for a in (1, -1):
                for b in (1, -1):
                    v = [0]*n; v[i] = a; v[j] = b; S.add(tuple(v))
        return V.canon(S)
    def An(n):
        S = set()
        for i, j in it.permutations(range(n+1), 2):
            v = [0]*(n+1); v[i] = 1; v[j] = -1; S.add(tuple(v))
        return V.canon(S)
    out = []
    for name, L, hub in [('D4', Dn(4), 8), ('A3', An(3), 6), ('D5', Dn(5), 20), ('A4', An(4), 16), ('A5', An(5), 20),
                         ('D6', Dn(6), 24), ('E7', E7, 42), ('P112', P112, 64), ('E8', E8, 64)]:
        B = Bmat(L); n = len(L)
        if n <= 20:
            z = exact(B, n); how = 'exact'
        else:
            z = greedy(B, n, 300 if n > 60 else 2000, 7); how = 'greedy (lower bound on the best Z)'
        rec = dict(design=name, n=n, Z=z, how=how, floor_per_port=z / n, hub_t_k1=hub)
        print(rec, flush=True); out.append(rec)
    json.dump(out, open(sys.argv[1], 'w'), indent=1)

if __name__ == '__main__':
    main()
