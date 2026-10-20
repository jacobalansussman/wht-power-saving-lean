"""padic.py -- mw-xw2tight (14 Oct 2026).  Exact 2-adic lattice potentials for the XW game at block 2.

O = Z_(2)[i], uniformiser pi = 1 - i (pi^2 = -2i, 2 = i pi^2).  Gaussian integers are stored as pairs of Python ints
(exact; no modular truncation).  A lattice is given by an integral basis matrix B (columns) of full rank.

  smith_vals(M)            : pi-valuations of the elementary divisors of an integral matrix M (rows x cols, full row rank)
  lam(g, B)                : length((g Lam + Lam)/Lam) for Lam = column span of B, g = pi^-s * G with G integral

Addresses at block 2 with f columns: x = sum_j col_j 4^j, col_j in 0..3, bit r of col_j = row r.  Bit index b = 2j + r.
"""
import itertools
from fractions import Fraction

# ---------------------------------------------------------------- Gaussian integers as (a, b) = a + b i
def gmul(x, y): return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def gadd(x, y): return (x[0]+y[0], x[1]+y[1])
def gsub(x, y): return (x[0]-y[0], x[1]-y[1])
def gneg(x): return (-x[0], -x[1])
ZERO, ONE, I = (0, 0), (1, 0), (0, 1)
PI = (1, -1)
IPOW = [(1, 0), (0, 1), (-1, 0), (0, -1)]

def vpi(x):
    """pi-adic valuation of a Gaussian integer (None for 0)"""
    a, b = x
    if a == 0 and b == 0: return None
    v = 0
    while True:
        if (a - b) % 2:          # a + b i with a, b of different parity: a unit
            return v
        # divisible by pi = 1 - i:  (a + b i)/(1 - i) = (a + b i)(1 + i)/2 = ((a - b) + (a + b) i)/2
        a, b = (a - b) // 2, (a + b) // 2
        v += 1

def divpi(x, k=1):
    a, b = x
    for _ in range(k):
        assert (a - b) % 2 == 0
        a, b = (a - b) // 2, (a + b) // 2
    return (a, b)

def unit_inv_times(u, x, prec):
    """x * u^-1 modulo pi^prec-ish: u a unit; we return exact Gaussian rational? -- avoided; see smith_vals"""
    raise NotImplementedError

# ---------------------------------------------------------------- Smith valuations by elimination over Z[i] localised at pi
def smith_vals(M):
    """M: list of rows, each a list of Gaussian ints. Returns sorted list of pi-valuations of elementary divisors
    (length = rank).  Exact: works over Z[i] with the trick that dividing by a unit u = a+bi is replaced by
    multiplying the other row by u and the pivot row by the entry (both rows stay integral; the operations
    are invertible over O because u is a unit and we only ever scale by units)."""
    M = [list(r) for r in M]
    nr = len(M); nc = len(M[0]) if nr else 0
    vals = []
    rows = list(range(nr)); cols = list(range(nc))
    # strip common pi-powers is not needed; we do valuation pivoting
    while rows and cols:
        best = None
        for r in rows:
            Mr = M[r]
            for c in cols:
                x = Mr[c]
                if x[0] == 0 and x[1] == 0: continue
                v = vpi(x)
                if best is None or v < best[0]:
                    best = (v, r, c)
                    if v == 0: break
            if best is not None and best[0] == 0: break
        if best is None: break
        v, pr, pc = best
        vals.append(v)
        p = M[pr][pc]
        # p = pi^v * u ; for each other row r with entry e at pc (vpi(e) >= v): row_r <- u*row_r - (e/pi^v) * row_p_unitpart
        pu = divpi(p, v)                          # unit part
        for r in rows:
            if r == pr: continue
            e = M[r][pc]
            if e == ZERO: continue
            q = divpi(e, v)                       # integral since vpi(e) >= v
            Mr = M[r]; Mp = M[pr]
            # new row = pu * row_r - q * row_p  (multiplying a row by the unit pu is invertible over O)
            M[r] = [gsub(gmul(pu, Mr[c]), gmul(q, Mp[c])) for c in range(nc)]
        # column elimination is not needed for the valuations: after clearing the pivot column, remove pivot row/col
        rows.remove(pr); cols.remove(pc)
        # but entries in the pivot row in other columns could be needed? No: removing row pr after clearing column pc
        # leaves the module of the remaining rows (projection); valuations of the remaining block are correct because
        # the pivot has minimal valuation (standard DVR Smith argument).
    return sorted(vals)


# ---------------------------------------------------------------- block-2 letters at f columns
def nbits(f): return 2 * f
def bit(x, b): return (x >> b) & 1

class Mono:
    """(g v)[sig[x]] = i^ph[x] v[x]"""
    def __init__(self, sig, ph): self.sig = list(sig); self.ph = [p % 4 for p in ph]
    def then(self, other):  # other after self
        return Mono([other.sig[s] for s in self.sig], [(self.ph[x] + other.ph[self.sig[x]]) % 4 for x in range(len(self.sig))])
    def inv(self):
        n = len(self.sig); sig = [0]*n; ph = [0]*n
        for x in range(n): sig[self.sig[x]] = x; ph[self.sig[x]] = (-self.ph[x]) % 4
        return Mono(sig, ph)
    def dense(self):
        n = len(self.sig); G = [[ZERO]*n for _ in range(n)]
        for x in range(n): G[self.sig[x]][x] = IPOW[self.ph[x]]
        return G

def linmap(f, A):
    """address permutation x -> A x over F_2 (A: list of rows, n x n bits)"""
    n = 2 * f; N = 1 << n
    sig = []
    for x in range(N):
        y = 0
        for i in range(n):
            s = 0
            for j in range(n):
                if A[i][j]: s ^= bit(x, j)
            y |= s << i
        sig.append(y)
    return Mono(sig, [0]*N)

def phase(f, q):
    N = 1 << (2 * f)
    return Mono(list(range(N)), [q(x) % 4 for x in range(N)])

def shift(f, mask):
    N = 1 << (2 * f)
    return Mono([x ^ mask for x in range(N)], [0]*N)

def ident(n): return [[1 if i == j else 0 for j in range(n)] for i in range(n)]

def xw_generators(f, cs=True, chain=True, boundary=True):
    """monomial generators of the XW free group at block 2, f columns (boundary = first/last column letters)"""
    n = 2 * f; gens = {}
    b = lambda r, j: 2 * j + r
    # column-uniform GL(2,2): row swap, row0 ^= row1, row1 ^= row0
    for name, (a00, a01, a10, a11) in {"swap": (0, 1, 1, 0), "C01": (1, 1, 0, 1), "C10": (1, 0, 1, 1)}.items():
        A = [[0]*n for _ in range(n)]
        for j in range(f):
            A[b(0, j)][b(0, j)] = a00; A[b(0, j)][b(1, j)] = a01
            A[b(1, j)][b(0, j)] = a10; A[b(1, j)][b(1, j)] = a11
        gens[name] = linmap(f, A)
    # column-uniform Paulis and phases
    for r in range(2):
        gens["X%d" % r] = shift(f, sum(1 << b(r, j) for j in range(f)))
        gens["S%d" % r] = phase(f, lambda x, r=r: sum(bit(x, b(r, j)) for j in range(f)))
    gens["CZu"] = phase(f, lambda x: 2 * sum(bit(x, b(0, j)) * bit(x, b(1, j)) for j in range(f)))
    if cs:
        gens["CSu"] = phase(f, lambda x: sum(bit(x, b(0, j)) * bit(x, b(1, j)) for j in range(f)))
    if chain:
        for name, (rt, rs) in {"chi0": (0, 0), "chi1": (1, 1), "E": (0, 1), "Em": (1, 0)}.items():
            A = ident(n)
            for j in range(1, f): A[b(rt, j)][b(rs, j - 1)] ^= 1
            gens[name] = linmap(f, A)
        if f > 3:            # R3 = row ^= row << 3 on both rows
            A = ident(n)
            for j in range(3, f):
                for r in range(2): A[b(r, j)][b(r, j - 3)] ^= 1
            gens["R3"] = linmap(f, A)
        for r in range(2):   # parity phases P_r = i^(parity of row r)
            gens["P%d" % r] = phase(f, lambda x, r=r: sum(bit(x, b(r, j)) for j in range(f)) % 2)
    if boundary:
        for jj, tag in ((0, "f"), (f - 1, "l")):
            for r in range(2):
                gens["X%s%d" % (tag, r)] = shift(f, 1 << b(r, jj))
                gens["S%s%d" % (tag, r)] = phase(f, lambda x, r=r, jj=jj: bit(x, b(r, jj)))
            gens["CZ%s" % tag] = phase(f, lambda x, jj=jj: 2 * bit(x, b(0, jj)) * bit(x, b(1, jj)))
            for name, (a00, a01, a10, a11) in {"swap": (0, 1, 1, 0), "C01": (1, 1, 0, 1), "C10": (1, 0, 1, 1)}.items():
                A = ident(n)
                A[b(0, jj)][b(0, jj)] = a00; A[b(0, jj)][b(1, jj)] = a01
                A[b(1, jj)][b(0, jj)] = a10; A[b(1, jj)][b(1, jj)] = a11
                gens["%s%s" % (name, tag)] = linmap(f, A)
            if cs:
                gens["CS%s" % tag] = phase(f, lambda x, jj=jj: bit(x, b(0, jj)) * bit(x, b(1, jj)))
    return gens

# ---------------------------------------------------------------- paid moves and kernel as pi^-s * integral matrix
def matmul(A, B):
    n = len(A); m = len(B[0]); k = len(B)
    C = [[ZERO]*m for _ in range(n)]
    for i in range(n):
        Ai = A[i]; Ci = C[i]
        for t in range(k):
            a = Ai[t]
            if a == ZERO: continue
            Bt = B[t]
            for j in range(m):
                bt = Bt[j]
                if bt != ZERO:
                    Ci[j] = gadd(Ci[j], gmul(a, bt))
    return C

def one_minus_iX(f, mask):
    """integral matrix 1 - i X^mask"""
    N = 1 << (2 * f)
    M = [[ZERO]*N for _ in range(N)]
    for x in range(N):
        M[x][x] = gadd(M[x][x], ONE)
        M[x ^ mask][x] = gadd(M[x ^ mask][x], (0, -1))
    return M

def paid(f, z):
    """D_z = pi^-f prod_j (1 - i X^z_j); returns (s, G) with D_z = pi^-s G"""
    N = 1 << (2 * f)
    G = ident_g(N)
    for j in range(f):
        G = matmul(one_minus_iX(f, z << (2 * j)), G)
    return f, G

def kernel(f):
    N = 1 << (2 * f)
    G = ident_g(N)
    for bb in range(2 * f):
        G = matmul(one_minus_iX(f, 1 << bb), G)
    return 2 * f, G

def ident_g(N): return [[ONE if i == j else ZERO for j in range(N)] for i in range(N)]

def lam(g, B):
    """g = (s, G) meaning pi^-s G; B basis (N x N) of Lam. length((g Lam + Lam)/Lam)."""
    s, G = g
    GB = matmul(G, B)
    pisB = [[gmul(x, pipow(s)) for x in row] for row in B]
    top = smith_vals([pisB[i] + GB[i] for i in range(len(B))])
    bot = smith_vals(pisB)
    return sum(bot) - sum(top)

_pp = {}
def pipow(k):
    if k not in _pp:
        x = ONE
        for _ in range(k): x = gmul(x, PI)
        _pp[k] = x
    return _pp[k]

def mono_as_g(m): return (0, m.dense())

# ---------------------------------------------------------------- Reed-Muller graded lattices
def monomials(n, r):
    return [S for d in range(0, r + 1) for S in itertools.combinations(range(n), d)]

def rm_lattice(f, rs):
    """Lam = sum_k pi^k * Zspan{0/1 monomials of degree <= r_k} + pi^len(rs) L ; returns a square basis (via Smith-free
    generating set reduced to N columns by Hermite-like elimination)."""
    n = 2 * f; N = 1 << n
    gensv = []
    for k, r in enumerate(rs):
        for S in monomials(n, min(r, n)) if r >= 0 else []:
            v = [pipow(k) if all(bit(x, b) for b in S) else ZERO for x in range(N)]
            gensv.append(v)
    K = len(rs)
    for x in range(N):
        v = [ZERO]*N; v[x] = pipow(K); gensv.append(v)
    return reduce_basis(gensv, N)

def reduce_basis(vecs, N):
    """column-style Hermite reduction over O: return N vectors spanning the same O-module (full rank assumed)."""
    # put vectors as rows; do row reduction with valuation pivots on columns, producing a triangular basis
    rows = [list(v) for v in vecs]
    basis = []
    for c in range(N):
        best = None
        for idx, r in enumerate(rows):
            x = r[c]
            if x == ZERO: continue
            v = vpi(x)
            if best is None or v < best[0]: best = (v, idx)
        if best is None: raise ValueError("not full rank")
        v, idx = best
        p = rows.pop(idx); pu = divpi(p[c], v)
        new = []
        for r in rows:
            e = r[c]
            if e == ZERO: new.append(r); continue
            q = divpi(e, v)
            new.append([gsub(gmul(pu, r[t]), gmul(q, p[t])) for t in range(N)])
        rows = new
        basis.append(p)
    # basis as columns
    return [[basis[j][i] for j in range(N)] for i in range(N)]
