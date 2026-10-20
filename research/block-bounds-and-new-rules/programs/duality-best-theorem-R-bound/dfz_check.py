# dfz_check.py -- sanity test of the DFZ Theorem 8.6 encoding used in lp3.py / lpm.py (pf-duality, 14 Oct 2026).
# Pure Python. Ranks of sums of subspaces over GF(p), p odd (the theorem holds for every field with 2 invertible,
# arXiv:1311.4601 q.txt:1775-2039), so a violation here would mean the encoding is wrong.  Also checks that the
# GF(2) Fano arrangement of the paper's Theorem 8.7 violates it (the encoding is not vacuous).
import random, sys

def rank_mod(rows, p):
    rows = [r[:] for r in rows]; rk = 0; ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] % p), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], p - 2, p) if p > 2 else 1
        rows[rk] = [(v * inv) % p for v in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % p:
                f = rows[i][c]; rows[i] = [(a - f * b) % p for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk

def dfz_value(sp, p):
    """sp: dict letter -> list of row vectors (spanning set). returns RHS - LHS of (65)."""
    def H(*ls):
        rows = [v for l in ls for v in sp[l]]
        return rank_mod(rows, p) if rows else 0
    def Hc(x, ys): return H(x, *ys) - H(*ys)
    lhs = 2 * H('A') + H('B') + 2 * H('C')
    rhs = (H('W') + H('X') + H('Y') + H('Z') + 2 * Hc('A', ['Z', 'Y']) + Hc('B', ['X', 'Z']) + 2 * Hc('C', ['A', 'X'])
           + 3 * Hc('X', ['W', 'Y']) + 3 * Hc('Z', ['W', 'C']) + 5 * Hc('W', ['A', 'B']) + 5 * Hc('Y', ['B', 'C'])
           + 5 * (H('A') + H('B') + H('C') - H('A', 'B', 'C')))
    return rhs - lhs

L = ['A', 'B', 'C', 'W', 'X', 'Y', 'Z']
fano = dict(A=[[1, 0, 0]], B=[[0, 1, 0]], C=[[0, 0, 1]], W=[[1, 1, 0]], X=[[1, 0, 1]], Y=[[0, 1, 1]], Z=[[1, 1, 1]])
print("GF(2) Fano arrangement (DFZ Thm 8.7): RHS-LHS =", dfz_value(fano, 2), "(must be < 0)")
print("same vectors over GF(3): RHS-LHS =", dfz_value(fano, 3), "(must be >= 0)")
rng = random.Random(1)
worst = None; n = 0
for p in (3, 5, 7, 1000003):
    for trial in range(4000 if p < 100 else 1500):
        dim = rng.randint(2, 6)
        sp = {}
        mode = trial % 3
        base = [[rng.randrange(p) for _ in range(dim)] for _ in range(3)]
        for l in L:
            k = rng.randint(0, 2)
            vecs = [[rng.randrange(p) for _ in range(dim)] for _ in range(k)]
            if mode >= 1:   # near-Fano: combinations of three base vectors with small coefficients
                vecs = [[sum(rng.choice([0, 1, 1, p - 1]) * b[i] for b in base) % p for i in range(dim)] for _ in range(max(k, 1))]
            sp[l] = vecs if vecs else [[0] * dim]
        v = dfz_value(sp, p); n += 1
        if worst is None or v < worst: worst = v
        assert v >= 0, ("VIOLATION", p, sp, v)
print("random arrangements over GF(3), GF(5), GF(7), GF(1000003):", n, "trials, min RHS-LHS =", worst, "(no violation)")
