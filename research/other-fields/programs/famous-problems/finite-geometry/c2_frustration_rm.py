# c2: "frustrated sign patterns = words outside RM(1,k)" (ms-fano section 3), re-derived as an exact
# linear-algebra statement and checked, plus the facts that make RM(1,k) a famous object at larger k.
# Pieces p_v = x1 + s(v) Z^v x2 (v in S, a subset of F_2^k; Z^v free Pauli relabellings).
# Separation: g_v in Q[F_2^k] with sum_v g_v = 1 and sum_v s(v) g_v Z^v = 0  (then sum g_v p_v = x1).
# Exact Gaussian elimination over Q. Claim: solvable iff s restricted to S is not +-(an affine character).
import itertools, random
from fractions import Fraction
def solvable(k, S, s):
    N = 1 << k
    idx = {(v, w): i for i, (v, w) in enumerate(itertools.product(S, range(N)))}   # a_{v,w}: coeff of Z^w in g_v
    rows = []
    for w in range(N):   # sum_v a_{v,w} = [w == 0]
        r = [Fraction(0)] * (len(idx) + 1)
        for v in S: r[idx[(v, w)]] = Fraction(1)
        r[-1] = Fraction(1 if w == 0 else 0); rows.append(r)
    for w in range(N):   # coefficient of Z^w in sum_v s(v) g_v Z^v: sum_v s(v) a_{v, w^v} = 0
        r = [Fraction(0)] * (len(idx) + 1)
        for v in S: r[idx[(v, w ^ v)]] += s[v]
        rows.append(r)
    # eliminate; inconsistent iff a row 0 = nonzero
    m, ncol, piv = rows, len(idx), 0
    for c in range(ncol):
        p = next((i for i in range(piv, len(m)) if m[i][c] != 0), None)
        if p is None: continue
        m[piv], m[p] = m[p], m[piv]
        inv = 1 / m[piv][c]; m[piv] = [x * inv for x in m[piv]]
        for i in range(len(m)):
            if i != piv and m[i][c] != 0:
                f = m[i][c]; m[i] = [a - f * b for a, b in zip(m[i], m[piv])]
        piv += 1
    return all(any(x != 0 for x in r[:-1]) or r[-1] == 0 for r in m)
def affine_on(k, S, s):   # is s|S = c * (-1)^{a.v + b} for some a, b?
    return any(all(s[v] == (1 - 2 * ((bin(a & v).count("1") + b) & 1)) * s[S[0]] * (1 - 2 * ((bin(a & S[0]).count("1") + b) & 1)) for v in S)
               for a in range(1 << k) for b in (0,))
def dist_RM1(k, f):   # f: list of 0/1 of length 2^k; distance to RM(1,k)
    N = 1 << k; best = N
    for a in range(N):
        for b in (0, 1):
            d = sum(1 for v in range(N) if f[v] != ((bin(a & v).count("1") + b) & 1)); best = min(best, d)
    return best
for k in (2, 3):
    N = 1 << k; S = list(range(N)); agree = 0; nfr = 0
    for bits in range(1 << N):
        s = {v: (-1 if bits >> v & 1 else 1) for v in S}
        sol = solvable(k, S, s); aff = affine_on(k, S, s); agree += (sol == (not aff)); nfr += sol
    print(f"k={k}: all {1<<N} sign patterns on F_2^{k}: 'separable' == 'not +-affine' for {agree}; separable {nfr}; RM(1,{k}) cosets {(1<<N)//(1<<(k+1))} (frustrated classes {(1<<N)//(1<<(k+1)) - 1})")
# minimal frustrated supports, k = 3 and 4
for k in (3, 4):
    N = 1 << k; found = None
    for size in range(1, N + 1):
        for S in itertools.combinations(range(N), size):
            for bits in range(1 << size):
                s = {v: (-1 if bits >> i & 1 else 1) for i, v in enumerate(S)}
                if not affine_on(k, list(S), s):
                    found = (size, S, bits); break
            if found: break
        if found: break
    size, S, bits = found
    x = 0
    for v in S: x ^= v
    print(f"k={k}: smallest frustrated support has {size} points, e.g. {S} (XOR of the points = {x}: an affine plane), pattern bits {bits}; exact solve agrees: {solvable(k, list(S), {v: (-1 if bits >> i & 1 else 1) for i, v in enumerate(S)})}")
random.seed(1)
k = 4; N = 16; ok = 0; T = 60
for t in range(T):
    size = random.randint(4, 7); S = sorted(random.sample(range(N), size)); s = {v: random.choice((1, -1)) for v in S}
    ok += (solvable(k, S, s) == (not affine_on(k, S, s)))
print(f"k=4 random supports: criterion agrees with exact solve in {ok}/{T}")
# covering radius of RM(1,k) = max distance to RM(1,k) (k <= 4 by brute force; famous open for odd k >= 9)
for k in (2, 3, 4):
    N = 1 << k
    print(f"covering radius of RM(1,{k}) = max nonlinearity on F_2^{k}: {max(dist_RM1(k, [(b >> v) & 1 for v in range(N)]) for b in range(1 << N))}")
