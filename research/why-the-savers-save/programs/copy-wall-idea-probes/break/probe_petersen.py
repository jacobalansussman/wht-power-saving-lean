"""probe_petersen.py (g-break, H2) -- exact minrank over F_q of the D5 bit pattern reduced to its 10 supports
(the Petersen pattern: B[A,B] = 0 when supports share exactly one index, B[A,A] != 0, free when disjoint).
The D5 bit unit (20 lines) has the same minrank (both signs of a support share a row; they know each other).
Natural F_2 decoder (1+G) mod 2 has rank 5; clique bound alpha = 4.  Question: rank 4 over some field?"""
import itertools, sys
pairs = list(itertools.combinations(range(5), 2))
n = 10
known = [[b for b in range(n) if not set(pairs[a]) & set(pairs[b])] for a in range(n)]
def reduce_q(x, basis, q):
    # basis: list of (pivot, vec) in echelon form, vec[pivot] = 1
    x = list(x)
    for p, b in basis:
        if x[p]:
            f = x[p]
            x = [(u - f * w) % q for u, w in zip(x, b)]
    return x
def add_q(x, basis, q):
    p = next(i for i in range(n) if x[i])
    inv = pow(x[p], q - 2, q)
    x = [u * inv % q for u in x]
    nb = []
    for pp, b in basis:
        if b[p]:
            f = b[p]; b = [(u - f * w) % q for u, w in zip(b, x)]
        nb.append((pp, b))
    return nb + [(p, x)]
def minrank_le(r, q):
    sols = [0]
    def dfs(k, basis):
        if k == n: return True
        kn = known[k]
        for vals in itertools.product(range(q), repeat=len(kn)):
            x = [0] * n; x[k] = 1
            for T, c in zip(kn, vals): x[T] = c
            y = reduce_q(x, basis, q)
            if not any(y):
                if dfs(k + 1, basis): return True
            elif len(basis) < r:
                if dfs(k + 1, add_q(y, basis, q)): return True
        return False
    return dfs(0, [])
for q in (2, 3, 5, 7, 11, 13):
    print('Petersen/D5 pattern over F_%d: minrank <= 4: %s ; <= 5: %s' % (q, minrank_le(4, q), minrank_le(5, q)))
    sys.stdout.flush()
