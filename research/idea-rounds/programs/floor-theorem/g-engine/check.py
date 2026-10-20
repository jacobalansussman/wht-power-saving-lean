#!/usr/bin/python3
"""check.py (agent g-engine, skeptic of f-engine, 10 Oct 2026).  Own code, nothing imported from floor.py.
Exact arithmetic in the prime field F_p, p = 998244353 (p = 1 mod 4), i := 3^((p-1)/4), so i*i = -1.
Every statement tested is an identity or a dimension count that the written proof claims over any field
with 1/2 and i, so a failure in F_p would be a failure of the proof.  Pure Python, no numpy.
NOT Lean.  Nobody has read this.
"""
import sys, random, itertools
P = 998244353
I = pow(3, (P - 1) // 4, P)
assert I * I % P == P - 1
INV2 = pow(2, P - 2, P)
CA = (1 + I) * INV2 % P          # ca = (1+i)/2   BinaryFrames.lean:66
CB = (1 - I) * INV2 % P          # cb = (1-i)/2   BinaryFrames.lean:67

def zeros(r, c): return [[0] * c for _ in range(r)]
def eye(n): return [[1 if i == j else 0 for j in range(n)] for i in range(n)]
def mul(A, B):
    n, k, m = len(A), len(B), len(B[0]); Bt = list(zip(*B))
    return [[sum(a * b for a, b in zip(row, col)) % P for col in Bt] for row in A]
def add(A, B): return [[(x + y) % P for x, y in zip(r, s)] for r, s in zip(A, B)]
def scal(c, A): return [[c * x % P for x in r] for r in A]
def kron(A, B):
    return [[a * b % P for a in ra for b in rb] for ra in A for rb in B]
def rank(rows):
    M = [r[:] for r in rows if any(r)]; rk = 0
    if not M: return 0
    nc = len(M[0])
    for c in range(nc):
        p = next((r for r in range(rk, len(M)) if M[r][c]), None)
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]; inv = pow(M[rk][c], P - 2, P)
        M[rk] = [x * inv % P for x in M[rk]]
        for r in range(len(M)):
            if r != rk and M[r][c]:
                t = M[r][c]; M[r] = [(x - t * y) % P for x, y in zip(M[r], M[rk])]
        rk += 1
        if rk == len(M): break
    return rk
def cols(A): return [list(c) for c in zip(*A)]

# ---- the Lean matrices, typed again from the files ----
def shift(m, z):                       # shift z i j = 1 iff i+z = j     BinaryFrames.lean:137-138
    n = 1 << m; A = zeros(n, n)
    for i in range(n): A[i][i ^ z] = 1
    return A
def dirm(m, z): return add(scal(CA, eye(1 << m)), scal(CB, shift(m, z)))     # BinaryFrames.lean:152
def littleC(): return [[CA, CB], [CB, CA]]                                    # BinaryFrames.lean:68
def kernel(m):
    K = [[1]]
    for _ in range(m): K = kron(K, littleC())
    return K
def indep(zs):
    span = {0}
    for z in zs:
        if z in span: return False
        span |= {s ^ z for s in span}
    return True
def block(m, zs):
    M = eye(1 << m)
    for z in zs: M = mul(M, dirm(m, z))
    return M
def delta(n):                          # basis vectors e_x (x) e_x of U (x) U, as rows
    return [[1 if j == x * n + x else 0 for j in range(n * n)] for x in range(n)]
def Fpair(g):                          # F(g(x)g) = dim(Delta + (g(x)g)Delta) - dim Delta
    n = len(g); D = delta(n); G = kron(g, g)
    img = cols(mul(G, cols(D)))
    return rank(D + img) - n

def part1():
    print("PART 1  Lemma W, own transcription, exact in F_p")
    for m in (1, 2, 3):
        n = 1 << m
        K = kernel(m)
        assert K == block(m, [1 << j for j in range(m)]), "kernel is not the product of the m directions"
        for k in range(1, m + 1):
            vals = set()
            for zs in itertools.product(range(1, n), repeat=k):
                if indep(zs): vals.add(Fpair(block(m, zs)))
            print(f"  m={m} rank {k}{' (kernel rank)' if k == m else ''}: F over all independent tuples = {sorted(vals)}; predicted {n - (n >> k)}")
    # m = 4: a few tuples only (64 s budget)
    m = 4; n = 16; random.seed(7)
    for k in (1, 3, 4):
        vals = set()
        for _ in range(3):
            while True:
                zs = [random.randrange(1, n) for _ in range(k)]
                if indep(zs): break
            vals.add(Fpair(block(m, zs)))
        print(f"  m=4 rank {k}: 3 random tuples: F = {sorted(vals)}; predicted {n - (n >> k)}")

# ---- the potential, class P(m,2) widened: any mixing (also singular), any b with at most one non-zero entry per
#      column (also singular), any paid g (also singular) ----
def rnd(): return random.randrange(P)
def col_monomial(n, singular):
    A = zeros(n, n)
    if singular: tgt = [random.randrange(n) for _ in range(n)]
    else: tgt = list(range(n)); random.shuffle(tgt)
    for x in range(n):
        A[tgt[x]][x] = random.randrange(1, P) if not (singular and random.random() < 0.2) else 0
    return A
def phi(S, R, W, n):
    """S[r] = matrix (n*n rows) x (W*n*n cols) of array r.  A = image of Delta^W.  Potential = dim((A+Z)/Z)."""
    nn = n * n; rows = []
    ycols = [l * nn + x * n + x for l in range(W) for x in range(n)]
    off = [x * n + y for x in range(n) for y in range(n) if x != y]
    # generators of A as vectors over all roles, kept only on the off-diagonal coordinates (quotient by Z)
    for c in ycols:
        v = []
        for r in range(R): v += [S[r][o][c] for o in off]
        rows.append(v)
    return rank(rows)
def run_random(m, R, W, nwords, nlet, planted=None, seed=1):
    random.seed(seed); n = 1 << m; nn = n * n
    viol = 0; rises = 0; total = 0; maxrise = 0
    for _ in range(nwords):
        S = [zeros(nn, W * nn) for _ in range(R)]
        for l in range(W):
            for j in range(nn): S[l][j][l * nn + j] = 1
        cur = phi(S, R, W, n)
        for _ in range(nlet):
            t = random.random(); total += 1
            if planted is not None and t < 0.15:
                r = random.randrange(R); G = kron(planted, planted); S[r] = mul(G, S[r]); cap = 0
            elif t < 0.35:      # mixing, sometimes singular
                g = [[rnd() if random.random() < 0.6 else 0 for _ in range(R)] for _ in range(R)]
                S = [[[sum(g[r][j] * S[j][a][b] for j in range(R)) % P for b in range(W * nn)] for a in range(nn)] for r in range(R)]
                cap = 0
            elif t < 0.6:       # free one-array step
                r = random.randrange(R); b = col_monomial(n, random.random() < 0.3)
                S[r] = mul(kron(b, b), S[r]); cap = 0
            else:               # paid step: arbitrary matrix, sometimes singular, sometimes a Lean block
                r = random.randrange(R); u = random.random()
                if u < 0.4:
                    while True:
                        k = random.randrange(1, m + 1); zs = [random.randrange(1, n) for _ in range(k)]
                        if indep(zs): break
                    g = block(m, zs)
                else:
                    g = [[rnd() if random.random() < 0.7 else 0 for _ in range(n)] for _ in range(n)]
                    if u > 0.85: g[0] = [0] * n
                cap = Fpair(g); S[r] = mul(kron(g, g), S[r])
            new = phi(S, R, W, n)
            if new > cur + cap: viol += 1
            if new > cur: rises += 1; maxrise = max(maxrise, new - cur)
            cur = new
    return total, rises, viol, maxrise

def part2():
    print("PART 2  the potential lemma on random words of the widened class (steps 2-4 of the proof)")
    for (m, R, W, nw, nl) in ((1, 3, 2, 300, 20), (2, 3, 2, 60, 16)):
        tot, rises, viol, mx = run_random(m, R, W, nw, nl, seed=11 + m)
        print(f"  m={m}, {R} arrays ({W} live): {tot} letters, potential rose at {rises} (largest rise {mx}), letters breaking the lemma: {viol}")
    # CONTROL: plant a non-monomial matrix as a 'free' letter; the test must fire
    for (m, R, W, nw, nl) in ((1, 3, 2, 100, 20), (2, 3, 2, 30, 16)):
        tot, rises, viol, mx = run_random(m, R, W, nw, nl, planted=dirm(m, 1), seed=5)
        print(f"  CONTROL m={m}: dir(e1) planted as a free letter: letters breaking the lemma: {viol} (must be > 0)")

# ---- T2 widened: arrays X_r (x) U (x) U, any linear map on the X parts free, paid step  Pi (x) (g(x)g) + (1-Pi) (x) 1 ----
def part3():
    print("PART 3  T2 (coverage form): weight of a paid step applied to a part of an array")
    random.seed(3)
    for m in (1, 2):
        n = 1 << m; nn = n * n
        for dx in (2, 3):
            bad = 0; trials = 0
            for _ in range(12 if m == 1 else 4):
                # one array X (x) U (x) U, Y = X (x) Delta.  Pi = projection of rank rho on X (coordinate or oblique)
                rho = random.randrange(1, dx)
                Bm = [[rnd() for _ in range(dx)] for _ in range(dx)]
                while rank(Bm) < dx: Bm = [[rnd() for _ in range(dx)] for _ in range(dx)]
                # inverse by Gauss-Jordan
                aug = [Bm[i] + eye(dx)[i] for i in range(dx)]
                for c in range(dx):
                    p = next(r for r in range(c, dx) if aug[r][c]); aug[c], aug[p] = aug[p], aug[c]
                    inv = pow(aug[c][c], P - 2, P); aug[c] = [x * inv % P for x in aug[c]]
                    for r in range(dx):
                        if r != c and aug[r][c]:
                            t = aug[r][c]; aug[r] = [(x - t * y) % P for x, y in zip(aug[r], aug[c])]
                Bi = [row[dx:] for row in aug]
                Dg = [[1 if (i == j and i < rho) else 0 for j in range(dx)] for i in range(dx)]
                Pi = mul(mul(Bm, Dg), Bi) if random.random() < 0.5 else Dg
                Q = add(eye(dx), scal(P - 1, Pi))
                g = [[rnd() for _ in range(n)] for _ in range(n)]
                if random.random() < 0.5:
                    while True:
                        k = random.randrange(1, m + 1); zs = [random.randrange(1, n) for _ in range(k)]
                        if indep(zs): break
                    g = block(m, zs)
                Ax = [[rnd() for _ in range(dx)] for _ in range(dx)]      # anything on the other columns
                L = add(kron(mul(mul(Pi, Ax), Pi), kron(g, g)), kron(Q, eye(nn)))
                Y = [[1 if j == a * nn + x * n + x else 0 for j in range(dx * nn)] for a in range(dx) for x in range(n)]
                LY = cols(mul(L, cols(Y)))
                Fgot = rank(Y + LY) - dx * n
                trials += 1
                if Fgot > rho * Fpair(g): bad += 1
            print(f"  m={m}, other columns of dimension {dx}: {trials} random controlled steps, rise above (rank of selector) x F(g): {bad}")
    # CONTROL / the hatch: a non-monomial matrix on ONE of the two columns only
    m = 1; n = 2
    D = delta(n); L = kron(littleC(), eye(n)); LD = cols(mul(L, cols(D)))
    print(f"  CONTROL: C on one column only (cheap, O(n)): Delta rises by {rank(D + LD) - n} (not kept, as the hypothesis says)")
    # the missed hatch: an address permutation ACROSS two columns (x,y) -> (x, x+y): monomial, one pass, cost O(n)
    A = zeros(n * n, n * n)
    for x in range(n):
        for y in range(n): A[x * n + (x ^ y)][x * n + y] = 1
    AD = cols(mul(A, cols(D)))
    print(f"  cross-column address permutation (x,y)->(x,x+y): Delta rises by {rank(D + AD) - n} (a monomial step that does NOT keep Delta)")
    sw = zeros(n * n, n * n)
    for x in range(n):
        for y in range(n): sw[y * n + x][x * n + y] = 1
    SD = cols(mul(sw, cols(D)))
    print(f"  column swap (x,y)->(y,x): Delta rises by {rank(D + SD) - n}")

# ---- T1 end to end on right words: equality cases and the smallest sizes ----
def part4():
    print("PART 4  the two sides of T1 at the smallest sizes")
    for m in (1, 2, 3, 4):
        n = 1 << m
        need = n - 1
        print(f"  m={m}: kernel needs {need} of {n}; one block of rank m-1 gives {n - 2 if m > 1 else 0}; so per array at least "
              f"{'no block exists (rank<1): the kernel cannot be made at all' if m == 1 else f'{need}/{n-2} = {need/(n-2):.4f} blocks'}")
    # f = 4 certificate from the even-weight code, one bit: ratio dim(Y ^ gY)/dim Y for g = C on the bit
    C4 = kron(kron(littleC(), littleC()), kron(littleC(), littleC()))
    Y = [[1 if j == s else 0 for j in range(16)] for s in range(16) if bin(s).count("1") % 2 == 0]
    CY = cols(mul(C4, cols(Y)))
    inter = 2 * len(Y) - rank(Y + CY)
    print(f"  even-weight code on 4 columns, one bit: dim Y = {len(Y)}, dim(Y meet C^4 Y) = {inter}: ratio {inter}/{len(Y)} (Delta has 1/2)")
    C3 = kron(kron(littleC(), littleC()), littleC())
    Y = [[1 if j == s else 0 for j in range(8)] for s in (0, 7)]
    CY = cols(mul(C3, cols(Y)))
    print(f"  repetition code on 3 columns, one bit: dim Y = 2, dim(Y meet C^3 Y) = {4 - rank(Y + CY)}: ratio {4 - rank(Y + CY)}/2")

if __name__ == "__main__":
    for a in sys.argv[1:]:
        {"1": part1, "2": part2, "3": part3, "4": part4}.get(a, lambda: None)()

def part5():
    print("PART 5  T1 against the published block lists (arithmetic only; W and m from the Lean data files)")
    from fractions import Fraction as Fr
    lists = {"B2Gp193 (m=110, W=14692)": (110, 14692, [(1, 145890), (2, 61990), (3, 49615), (4, 23685), (5, 9750), (6, 8025), (7, 3530), (8, 7655), (9, 2145), (10, 4650), (11, 1770), (12, 5985), (13, 635), (14, 3710), (15, 9395), (16, 440), (17, 405), (18, 13265), (19, 835), (20, 110), (21, 2640), (42, 2640), (46, 2640)])}
    for name, (m, W, L) in lists.items():
        lhs = sum(c * (1 - Fr(1, 2 ** r)) for r, c in L); rhs = W * (1 - Fr(1, 2 ** m)); B = sum(c for r, c in L)
        print(f"  {name}: blocks {B} = {B / W:.3f} per live array; sum of weights {float(lhs):.1f} against needed {float(rhs):.1f}: "
              f"{'holds' if lhs >= rhs else 'VIOLATED'}, ratio {float(lhs / rhs):.2f}; total rank {sum(r * c for r, c in L)} against 2W(1-2^-m) = {float(2 * rhs):.0f}")

def part2b():
    print("PART 2b  the same potential test at m = 3 (few words)")
    tot, rises, viol, mx = run_random(3, 3, 2, 6, 10, seed=21)
    print(f"  m=3, 3 arrays (2 live): {tot} letters, potential rose at {rises} (largest rise {mx}), letters breaking the lemma: {viol}")

if __name__ == "__main__" and "5" in sys.argv[1:]: part5()
if __name__ == "__main__" and "2b" in sys.argv[1:]: part2b()
