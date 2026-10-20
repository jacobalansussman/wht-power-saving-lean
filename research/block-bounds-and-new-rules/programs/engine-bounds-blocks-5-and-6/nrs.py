# pf-induction, 13 Oct 2026.  CLOSED FORM of the projection counts N_r(s) (PROOF.md, Lemma 2), checked
# (a) against the exact table of b56-methods / v-b56-methods (m = 2..6, read from their output files), and
# (b) against a fresh brute-force count written here (every point p, a RANDOM pair of Lagrangians at
#     distance r built from a random symplectic basis, m = 2..5; all r).
# Usage: python3 -I nrs.py
import re, ast, random, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
W16 = os.path.join(HERE, '..', '..', 'wave16')

def N_closed(m, r):
    """number of nonzero p in F_2^{2m} with d(pi_p L, pi_p M) = s, for d(L, M) = r (1 <= r <= m)."""
    out = {}
    def add(s, c):
        if c: out[s] = out.get(s, 0) + c
    add(r, (2 ** (m - r) - 1) + (4 ** m - 2 ** (m + r)))            # p in L cap M, or p outside L + M
    add(r - 1, 2 * (2 ** m - 2 ** (m - r)) + (2 ** r - 1) * 2 ** (r - 1) * 2 ** (m - r))   # p in one of L, M; or p = l + m', <l,m'> = 1
    add(r - 2, (2 ** r - 1) * (2 ** (r - 1) - 1) * 2 ** (m - r))      # p = l + m', <l,m'> = 0
    assert sum(out.values()) == 4 ** m - 1
    return dict(sorted(out.items()))

# (a) the exact table
tab = {}
for fn in [os.path.join(W16, 'b56-methods', 'out_project.txt'), os.path.join(W16, 'v-b56-methods', 'out_projcheck.txt')]:
    for line in open(fn):
        mm = re.match(r"m=(\d+) r=(\d+): (?:N_r\(s\) = )?(\{[^}]*\})", line)
        if mm:
            key = (int(mm.group(1)), int(mm.group(2)))
            val = {int(k): v for k, v in ast.literal_eval(mm.group(3)).items()}
            if key in tab: assert tab[key] == val, (fn, key)
            tab[key] = val
bad = 0
for (m, r), val in sorted(tab.items()):
    if N_closed(m, r) != val:
        bad += 1; print("MISMATCH table", m, r, val, N_closed(m, r))
print(f"(a) closed form vs exact table: {len(tab)} (m, r) cases (m = 2..6), mismatches {bad}")

# (b) brute force with random symplectic bases
def form(m, a, b):
    M = (1 << m) - 1
    return (bin((a & M) & (b >> m)).count('1') + bin((b & M) & (a >> m)).count('1')) & 1
def span(B):
    out = {0}
    for b in B:
        if b not in out: out |= {x ^ b for x in out}
    return out
def rand_symp_basis(m, rng):
    # random symplectic basis e_1..e_m, f_1..f_m by symplectic Gram-Schmidt on random vectors
    while True:
        es, fs = [], []
        ok = True
        rest = list(range(1, 1 << (2 * m))); rng.shuffle(rest)
        pool = rest
        def proj(v):
            for e, f in zip(es, fs):
                if form(m, v, f): v ^= e
                if form(m, v, e): v ^= f
            return v
        for i in range(m):
            e = next((proj(v) for v in pool if proj(v) != 0), None)
            if e is None: ok = False; break
            f = next((proj(v) for v in pool if form(m, proj(v), e) == 1), None)
            if f is None: ok = False; break
            es.append(e); fs.append(f)
        if ok: return es, fs
def dist(m, A, B):
    return m - (len(A & B).bit_length() - 1)
def proj_lag(m, A, p):
    S = [v for v in A if form(m, v, p) == 0]
    return span(S + [p])
rng = random.Random(20261013)
bad = 0; cases = 0
for m in range(2, 6):
    for r in range(1, m + 1):
        for trial in range(2 if m < 5 else 1):
            es, fs = rand_symp_basis(m, rng)
            L = span(es)
            M = span(es[:m - r] + fs[m - r:])     # shares e_1..e_{m-r}
            assert dist(m, L, M) == r
            cnt = {}
            for p in range(1, 1 << (2 * m)):
                A, B = proj_lag(m, L, p), proj_lag(m, M, p)
                s = (m - 1) - (len(A & B).bit_length() - 1 - 1)
                cnt[s] = cnt.get(s, 0) + 1
            cases += 1
            if dict(sorted(cnt.items())) != N_closed(m, r):
                bad += 1; print("MISMATCH brute", m, r, cnt, N_closed(m, r))
print(f"(b) closed form vs fresh brute force (random symplectic frames): {cases} cases, m = 2..5, mismatches {bad}")
for m in range(2, 9):
    print(f"m={m}: " + "; ".join(f"r={r}: {N_closed(m, r)}" for r in (1, m - 1, m) if r >= 1))
