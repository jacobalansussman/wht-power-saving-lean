"""patterns.py (d-relay) -- side-information patterns for the relay study.  knows[S][T] (S != T) = S can cancel x_T
for free (bit side: lines orthogonal).  Standard library only."""
import itertools, sys
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-floor-checker'))

def from_gram(G):
    n = len(G)
    return [[S != T and G[S][T] == 0 for T in range(n)] for S in range(n)]

def c5():
    v = 5
    return [[(t - s) % v in (1, v - 1) for t in range(v)] for s in range(v)]

def petersen():
    pairs = list(itertools.combinations(range(5), 2))
    # D5 bit pattern reduced to supports: S knows T iff supports disjoint (g-break/probe_petersen.py)
    return [[a != b and not set(pairs[a]) & set(pairs[b]) for b in range(10)] for a in range(10)]

def lines_of(roots):
    out, seen = [], set()
    for r in roots:
        k = tuple(r); m = tuple(-x for x in r)
        if k in seen or m in seen: continue
        seen.add(k); out.append(list(r))
    return out

def dn(n):
    rs = []
    for i, j in itertools.combinations(range(n), 2):
        for s in (1, -1):
            v = [0] * n; v[i] = 1; v[j] = s; rs.append(v)
    L = lines_of(rs)
    return [[sum(a * b for a, b in zip(x, y)) for y in L] for x in L]

def e8_ours():
    from copyfloor import e8_lines
    return e8_lines()

def e6_gram():
    # same construction as g-break/probe_h2.py: lines of E8 orthogonal to line 0 and to the first line at 60 deg
    L, G = e8_ours()
    j60 = next(j for j in range(120) if G[0][j] != 0 and j != 0)
    e6 = [i for i in range(120) if i not in (0, j60) and G[i][0] == 0 and G[i][j60] == 0]
    return [[G[i][j] for j in e6] for i in e6]

def hamming_e8_hub():
    """Pure-Python re-implementation of the hub's s2/e8units.hamming_e8 (same generator, same order) so that the
    hub's P111 numbering 39,52,57,71,76,78,103,109,114 applies.  G = L L^T / 2."""
    gen = [[1, 1, 1, 1, 0, 0, 0, 0], [0, 0, 1, 1, 1, 1, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1], [0, 1, 0, 1, 0, 1, 0, 1]]
    code = set()
    for c in itertools.product([0, 1], repeat=4):
        code.add(tuple(sum(c[k] * gen[k][i] for k in range(4)) % 2 for i in range(8)))
    w4 = sorted(c for c in code if sum(c) == 4)
    assert len(w4) == 14
    L = [tuple(2 if i == a else 0 for i in range(8)) for a in range(8)]
    for c in w4:
        idx = [i for i in range(8) if c[i]]
        for sg in itertools.product([1, -1], repeat=3):
            x = [0] * 8; x[idx[0]] = 1
            for k in range(3): x[idx[k + 1]] = sg[k]
            L.append(tuple(x))
    G = [[sum(a * b for a, b in zip(x, y)) // 2 for y in L] for x in L]
    assert len(L) == 120 and all(G[i][i] == 2 for i in range(120))
    return L, G

HUB_P111 = [39, 52, 57, 71, 76, 78, 103, 109, 114]

def instance(name):
    """returns (knows, gram or None, label)"""
    if name == 'C5': return c5(), None
    if name == 'petersen': return petersen(), None
    if name == 'D5':
        G = dn(5); return from_gram(G), G
    if name == 'E6':
        G = e6_gram(); return from_gram(G), G
    if name == 'P120':
        L, G = e8_ours(); return from_gram(G), G
    if name == 'P120hub':
        L, G = hamming_e8_hub(); return from_gram(G), G
    if name == 'P111':
        L, G0 = hamming_e8_hub(); live = [i for i in range(120) if i not in HUB_P111]
        G = [[G0[a][b] for b in live] for a in live]; return from_gram(G), G
    if name == 'P111ours':
        L, G0 = e8_ours(); P = [8, 9, 10, 12, 16, 17, 24, 32, 64]; live = [i for i in range(120) if i not in P]
        G = [[G0[a][b] for b in live] for a in live]; return from_gram(G), G
    if name.startswith('sub:'):          # sub:out/shrink_L2.json = sub-pattern of P111 (indices into its live list)
        import json
        d = json.load(open(name[4:])); k, G = instance('P111'); ports = d['ports']
        return [[k[a][b] for b in ports] for a in ports], [[G[a][b] for b in ports] for a in ports]
    raise SystemExit('unknown instance ' + name)

def max_clique_nonknow(knows, size=None):
    """a largest set of pairwise non-knowing ports (clique bound for the bit pattern)"""
    n = len(knows)
    adj = [sum(1 << T for T in range(n) if T != S and not knows[S][T] and not knows[T][S]) for S in range(n)]
    best = []
    def rec(cur, cand):
        nonlocal best
        if len(cur) + bin(cand).count('1') <= len(best): return
        if size and len(best) >= size: return
        if not cand:
            if len(cur) > len(best): best = list(cur)
            return
        i = cand.bit_length() - 1
        rec(cur + [i], cand & adj[i])
        rec(cur, cand & ~(1 << i))
    rec([], (1 << n) - 1)
    return best
