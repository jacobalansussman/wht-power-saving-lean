# v-iso-combi check v2: does the separation criterion ("separable iff s is not +-chi_a on S") survive when the
# separating coefficients g_v may be ANY combination of Pauli letters (X and Z), i.e. any N x N matrix?
# Pieces p_v = x1 + s(v) Z^v x2. Need g_v (N x N) with sum_v g_v = I and sum_v s(v) g_v Z^v = 0.
# Brute-force linear solve (no use of the column decoupling), over GF(p), p = 2^31 - 1. Small: seconds, < 30 MB.
import itertools, random
P = 2147483647
wt = lambda x: bin(x).count("1")
def solvable(k, S, s):
    N = 1 << k
    var = {(v, r, c): i for i, (v, r, c) in enumerate(itertools.product(S, range(N), range(N)))}
    rows = []
    for r in range(N):
        for c in range(N):
            row = [0] * (len(var) + 1)
            for v in S: row[var[(v, r, c)]] = 1
            row[-1] = 1 if r == c else 0; rows.append(row)          # (sum g_v)_{rc} = I_{rc}
            row = [0] * (len(var) + 1)
            for v in S: row[var[(v, r, c)]] = (s[v] * (-1) ** wt(c & v)) % P   # (g_v Z^v)_{rc} = (g_v)_{rc} (-1)^{c.v}
            rows.append(row)
    piv = 0; ncol = len(var)
    for col in range(ncol):
        pr = next((i for i in range(piv, len(rows)) if rows[i][col]), None)
        if pr is None: continue
        rows[piv], rows[pr] = rows[pr], rows[piv]
        inv = pow(rows[piv][col], P - 2, P); rows[piv] = [x * inv % P for x in rows[piv]]
        for i in range(len(rows)):
            if i != piv and rows[i][col]:
                f = rows[i][col]; rows[i] = [(a - f * b) % P for a, b in zip(rows[i], rows[piv])]
        piv += 1
    return all(any(r[:-1]) or r[-1] == 0 for r in rows)
def pm_affine(k, S, s):
    return any(len({s[v] * (-1) ** wt(a & v) for v in S}) == 1 for a in range(1 << k))
for k, trials in ((2, None), (3, 24)):
    N = 1 << k; agree = tot = 0
    cases = []
    if trials is None:
        for size in range(1, N + 1):
            for S in itertools.combinations(range(N), size):
                for bits in range(1 << size):
                    cases.append((list(S), {v: (-1 if bits >> i & 1 else 1) for i, v in enumerate(S)}))
    else:
        random.seed(7)
        for _ in range(trials):
            size = random.choice((4, 5, 8)); S = sorted(random.sample(range(N), size))
            cases.append((S, {v: random.choice((1, -1)) for v in S}))
    for S, s in cases:
        tot += 1; agree += (solvable(k, S, s) == (not pm_affine(k, S, s)))
    print(f"k={k}: full matrix algebra (X and Z letters): criterion agrees with brute-force solve in {agree}/{tot} cases")
