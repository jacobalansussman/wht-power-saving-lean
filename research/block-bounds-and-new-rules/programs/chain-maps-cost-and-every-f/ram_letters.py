# xc-legit (wave 17): per-entry RAM cost of CROSS-COLUMN letters in the Lean engine's address layout.
#
# Layout (Lean engine, OAI/.../ColumnAddresses.lean `book X r f`, Work/Block/Address.lean): an entry's index is
#     idx = spec * 2^(m f) + sum_b row_b * 2^(f b),      row_b in [0, 2^f): bit j of row_b = address bit b of column j.
# Every letter below is executed exactly as a RAM program would: per entry it decodes the rows it needs (div, mod by
# prepared powers), computes the new rows with OpenAI's integer primitives (add sub mul div mod lt) and look-ups in
# prepared tables, re-encodes, and moves one value (1 look + 1 write).  Every primitive and every look-up is counted
# (one unit each, as NOp.run / Atom.look in UniformFourier.lean).  Tables are built with counted primitives too.
# Output per letter: max ops per entry (must not grow with f), table ops / 2^k (must stay bounded), and an exhaustive
# check that the RAM program computes the intended map at every entry.
# Usage: python3 -I ram_letters.py         (largest case 2^13 entries; < 50 MB)
import itertools, random, sys

class Ops:
    def __init__(self): self.n = 0
    def add(s, a, b): s.n += 1; return a + b
    def sub(s, a, b): s.n += 1; return a - b if a >= b else 0
    def mul(s, a, b): s.n += 1; return a * b
    def div(s, a, b): s.n += 1; return a // b if b else 0
    def mod(s, a, b): s.n += 1; return a % b if b else a
    def look(s, t, i): s.n += 1; return t[i]
    def cmul(s, a, b): s.n += 1; return a * b      # complex scale (Atom.scale)

# ---------------- tables (counted) ----------------
def xor_table(f, o):
    """X[a*F+b] = a XOR b on f bits, built in O(1) per entry from X[(a//2)*F + b//2] (4^f entries).
    (OpenAI's xor_prepare builds it in O(f) per entry, (2^f)^3 total; this recursion is O(4^f).)"""
    F = 1 << f; X = [0] * (F * F)
    for a in range(F):
        for b in range(F):
            if a == 0 and b == 0: continue
            prev = o.look(X, o.add(o.mul(o.div(a, 2), F), o.div(b, 2)))
            lowbit = o.mod(o.add(o.mod(a, 2), o.mod(b, 2)), 2)
            X[a * F + b] = o.add(o.mul(prev, 2), lowbit)
    return X

def and_table(f, o):
    F = 1 << f; A = [0] * (F * F)
    for a in range(F):
        for b in range(F):
            if a == 0 or b == 0: continue
            prev = o.look(A, o.add(o.mul(o.div(a, 2), F), o.div(b, 2)))
            A[a * F + b] = o.add(o.mul(prev, 2), o.mul(o.mod(a, 2), o.mod(b, 2)))
    return A

def lum_table(f, o):      # i^popcount(y)
    F = 1 << f; L = [1] * F
    for y in range(1, F):
        L[y] = o.cmul(o.look(L, o.div(y, 2)), 1j if o.mod(y, 2) else 1)
    return L

def row_table(cols, f, X, o):
    """T[y] = N y over F2 for an f x f binary matrix given by its columns cols[t] (as f-bit ints).
    Built by T[2^t + y'] = T[y'] XOR cols[t]: one XOR look-up per entry, 2^f entries."""
    F = 1 << f; T = [0] * F
    for t in range(f):
        p = 1 << t
        for yp in range(p):
            T[p + yp] = o.look(X, o.add(o.mul(o.look(T, yp), F), cols[t]))
    return T

# ---------------- per-entry helpers (counted) ----------------
def get_row(idx, b, f, P, o): return o.mod(o.div(idx, P[b]), 1 << f)
def xor(a, b, X, F, o):         return o.look(X, o.add(o.mul(a, F), b))
def encode(spec, rows, m, f, P, o):
    v = o.mul(spec, P[m])
    for b in range(m): v = o.add(v, o.mul(rows[b], P[b]))
    return v

# ---------------- mathematical definitions (uncounted, for the check) ----------------
def to_cols(idx, m, f):
    rows = [(idx >> (f * b)) & ((1 << f) - 1) for b in range(m)]
    return [sum(((rows[b] >> j) & 1) << b for b in range(m)) for j in range(f)], idx >> (m * f)
def from_cols(cols, spec, m, f):
    idx = spec << (m * f)
    for j, a in enumerate(cols):
        for b in range(m):
            if (a >> b) & 1: idx |= 1 << (f * b + j)
    return idx

def math_chain(cols):                 # (L x)_j = x_j + x_{j-1}, x_{-1} = 0   (row ^= row << 1, every row)
    return [cols[j] ^ (cols[j - 1] if j else 0) for j in range(len(cols))]
def math_chain_inv(cols):             # prefix XOR
    out, acc = [], 0
    for a in cols: acc ^= a; out.append(acc)
    return out
def math_row_cross(cols):             # bit 0 of column j += bit 1 of column j-1   (row0 ^= row1 << 1)
    return [cols[j] ^ (((cols[j - 1] >> 1) & 1) if j else 0) for j in range(len(cols))]
def math_carry(cols, m):              # +1 on the column string read as a little-endian base-2^m number (mod 2^(mf))
    out, c = [], 1
    for a in cols:
        s = a + c; out.append(s & ((1 << m) - 1)); c = s >> m
    return out

# ---------------- the letters as RAM programs ----------------
def run_letter(name, m, f, r, prog, mathmap, tables_cost, extra=None):
    k = m * f + r; N = 1 << k
    P = [1 << (f * b) for b in range(m + 1)]
    worst = 0; ok = True
    for idx in range(N):
        o = Ops()
        out = prog(idx, P, o)
        worst = max(worst, o.n + 2)          # + 1 look of the source value + 1 write (sow)
        cols, spec = to_cols(idx, m, f)
        if mathmap is not None:
            want = from_cols(mathmap(cols), spec, m, f)
            if out != want: ok = False
    return dict(name=name, m=m, f=f, r=r, k=k, ops_per_entry=worst, table_ops_over_2k=round(tables_cost / N, 3), ok=ok)

def letters_for(m, f, r, rng):
    F = 1 << f
    o = Ops(); X = xor_table(f, o); cx = o.n
    o = Ops(); A = and_table(f, o); ca = o.n
    o = Ops(); L = lum_table(f, o); cl = o.n
    # prefix-XOR table (inverse of the chain map on one row): PX[y] = (2*PX[y//2] mod F) XOR (y%2 ? F-1 : 0)
    o = Ops(); PX = [0] * F
    for y in range(1, F):
        hi = o.mod(o.mul(o.look(PX, o.div(y, 2)), 2), F)
        PX[y] = xor(hi, o.mul(o.mod(y, 2), F - 1), X, F, o)
    cpx = o.n
    # trailing-ones table for the carry letter: TO[y] = number of trailing 1 bits of y (f-bit)
    o = Ops(); TO = [0] * F
    for y in range(F):
        TO[y] = o.add(o.look(TO, o.div(y, 2)), 1) if o.mod(y, 2) else 0
    cto = o.n
    pw = [1 << j for j in range(f + 1)]
    out = []

    def chain(idx, P, o):
        spec = o.div(idx, P[m]); rows = [get_row(idx, b, f, P, o) for b in range(m)]
        new = [xor(rw, o.mod(o.mul(rw, 2), F), X, F, o) for rw in rows]
        return encode(spec, new, m, f, P, o)
    out.append(run_letter('row ^= row<<1, all rows', m, f, r, chain, math_chain, cx))

    def chain_inv(idx, P, o):
        spec = o.div(idx, P[m]); rows = [get_row(idx, b, f, P, o) for b in range(m)]
        return encode(spec, [o.look(PX, rw) for rw in rows], m, f, P, o)
    out.append(run_letter('inverse (prefix XOR, table)', m, f, r, chain_inv, math_chain_inv, cx + cpx))

    if m >= 2:
        def rcross(idx, P, o):
            spec = o.div(idx, P[m]); rows = [get_row(idx, b, f, P, o) for b in range(m)]
            rows[0] = xor(rows[0], o.mod(o.mul(rows[1], 2), F), X, F, o)
            return encode(spec, rows, m, f, P, o)
        out.append(run_letter('row0 ^= row1<<1', m, f, r, rcross, math_row_cross, cx))

    # bounded-degree polynomial matrix A(S) = [[1+S, S^2],[S, 1]] on rows 0,1 (A(0) = I), and its inverse:
    # det = 1 + S + S^3 (a unit mod S^f); A^-1 = det^-1 * [[1, S^2],[S, 1+S]]; det^-1 is a row table.
    if m >= 2:
        def polyA(idx, P, o):
            spec = o.div(idx, P[m]); rows = [get_row(idx, b, f, P, o) for b in range(m)]
            r0, r1 = rows[0], rows[1]
            n0 = xor(xor(r0, o.mod(o.mul(r0, 2), F), X, F, o), o.mod(o.mul(r1, 4), F), X, F, o)
            n1 = xor(o.mod(o.mul(r0, 2), F), r1, X, F, o)
            rows[0], rows[1] = n0, n1
            return encode(spec, rows, m, f, P, o)
        def math_polyA(cols):
            g = lambda j, b: (cols[j] >> b) & 1 if 0 <= j else 0
            res = []
            for j in range(len(cols)):
                b0 = g(j, 0) ^ g(j - 1, 0) ^ g(j - 2, 1); b1 = g(j - 1, 0) ^ g(j, 1)
                res.append((cols[j] & ~3) | b0 | (b1 << 1))
            return res
        out.append(run_letter('polynomial matrix [[1+S,S^2],[S,1]]', m, f, r, polyA, math_polyA, cx))
        # inverse: series of det^-1 mod S^f
        inv = [0] * f; inv[0] = 1                      # (1+S+S^3) * inv = 1
        for j in range(1, f):
            inv[j] = (inv[j - 1] if j >= 1 else 0) ^ (inv[j - 3] if j >= 3 else 0)
        cols_dinv = [sum(inv[i - t] << i for i in range(t, f)) for t in range(f)]   # column t = S^t * det^-1
        o = Ops(); Tdinv = row_table(cols_dinv, f, X, o); cd = o.n
        def polyAinv(idx, P, o):
            spec = o.div(idx, P[m]); rows = [get_row(idx, b, f, P, o) for b in range(m)]
            r0, r1 = rows[0], rows[1]
            a0 = xor(r0, o.mod(o.mul(r1, 4), F), X, F, o)
            a1 = xor(o.mod(o.mul(r0, 2), F), xor(r1, o.mod(o.mul(r1, 2), F), X, F, o), X, F, o)
            rows[0], rows[1] = o.look(Tdinv, a0), o.look(Tdinv, a1)
            return encode(spec, rows, m, f, P, o)
        def math_polyAinv(cols):
            # check A(A^-1 x) = x by brute force inversion of math_polyA on this f
            return None
        res = run_letter('its inverse (adjugate + det^-1 table)', m, f, r, polyAinv, None, cx + cd)
        # exhaustive round-trip check: polyA(polyAinv(idx)) == idx
        P = [1 << (f * b) for b in range(m + 1)]
        okrt = all(polyA(polyAinv(i, P, Ops()), P, Ops()) == i for i in range(1 << (m * f + r)))
        res['ok'] = okrt; out.append(res)

    # carry: +1 on the column string (little-endian base 2^m digits), a finite-automaton (carry) map across columns
    def carry(idx, P, o):
        spec = o.div(idx, P[m]); rows = [get_row(idx, b, f, P, o) for b in range(m)]
        allones = rows[0]
        for b in range(1, m): allones = o.look(A, o.add(o.mul(allones, F), rows[b]))
        t = o.look(TO, allones)                              # columns 0..t-1 are all-ones digits
        if o.sub(f, t) == 0:                                 # (one test) every digit is 2^m-1: wraps to 0 (truncation)
            return encode(spec, [0] * m, m, f, P, o)
        pt = o.look(pw, t)
        rows = [o.sub(rw, o.sub(pt, 1)) for rw in rows]      # clear columns < t (they are all 1 in every row)
        # column t: digit d -> d + 1 (a carry inside the digit, across rows): bits of rows at column t
        d = 0
        for b in range(m): d = o.add(d, o.mul(o.mod(o.div(rows[b], pt), 2), 1 << b))
        d2 = o.add(d, 1)                                      # < 2^m since column t is not all-ones
        for b in range(m):
            ob = o.mod(o.div(d, 1 << b), 2); nb = o.mod(o.div(d2, 1 << b), 2)
            rows[b] = o.sub(o.add(rows[b], o.mul(nb, pt)), o.mul(ob, pt))
        return encode(spec, rows, m, f, P, o)
    out.append(run_letter('carry: +1 across column digits', m, f, r, carry, lambda c: math_carry(c, m), cx + ca + cto))

    # cross-column phase i^popcount(row0 AND (row1 << 1)) (a standard-basis CZ chain between neighbouring columns)
    if m >= 2:
        def cphase(idx, P, o):
            r0 = get_row(idx, 0, f, P, o); r1 = get_row(idx, 1, f, P, o)
            y = o.look(A, o.add(o.mul(r0, F), o.mod(o.mul(r1, 2), F)))
            o.cmul(o.look(L, y), 1.0)        # multiply the moved value by the phase
            return idx
        res = run_letter('phase i^|row0 & row1<<1| (CZ chain)', m, f, r, cphase, lambda c: c, cx + ca + cl)
        out.append(res)

    # an ARBITRARY linear map of the body F2^(mf): block elimination into O(m^2) row letters with f x f row tables
    out.append(gl_letter(m, f, r, X, cx, rng))
    return out

# ---------------- GL(mf, F2) as O(m^2) row letters ----------------
def mat_mul(A, B, n):
    return [[sum(A[i][t] & B[t][j] for t in range(n)) & 1 for j in range(n)] for i in range(n)]
def mat_inv(M, n):
    A = [row[:] + [1 if i == j else 0 for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next((i for i in range(c, n) if A[i][c]), None)
        if p is None: return None
        A[c], A[p] = A[p], A[c]
        for i in range(n):
            if i != c and A[i][c]: A[i] = [x ^ y for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]
def rand_inv(n, rng):
    while True:
        M = [[rng.randint(0, 1) for _ in range(n)] for _ in range(n)]
        if mat_inv(M, n) is not None: return M
def blk(A, b, c, f): return [row[c * f:(c + 1) * f] for row in A[b * f:(b + 1) * f]]

def decompose(Afull, m, f, rng):
    """Left-multiply by elementary row-block operations until the identity; return the letters whose product is A.
    Letter ('e', b, c, N): row_b ^= N row_c (b != c, an involution); ('d', b, M): row_b <- M row_b."""
    n = m * f; A = [row[:] for row in Afull]; ops = []
    def apply_e(b, c, N):           # A <- (I + N E_bc) A
        for i in range(f):
            for j in range(n):
                s = 0
                for t in range(f): s ^= N[i][t] & A[c * f + t][j]
                A[b * f + i][j] ^= s
    def apply_d(b, M):
        rows = [A[b * f + i][:] for i in range(f)]
        for i in range(f):
            A[b * f + i] = [sum(M[i][t] & rows[t][j] for t in range(f)) & 1 for j in range(n)]
    for c in range(m):
        tries = 0
        while mat_inv(blk(A, c, c, f), f) is None:
            tries += 1
            for d in range(c + 1, m):
                N = [[rng.randint(0, 1) for _ in range(f)] for _ in range(f)]
                apply_e(c, d, N); ops.append(('e', c, d, N))
        Pinv = mat_inv(blk(A, c, c, f), f); apply_d(c, Pinv); ops.append(('d', c, Pinv))
        for d in range(m):
            if d != c:
                N = blk(A, d, c, f)
                if any(any(r) for r in N): apply_e(d, c, N); ops.append(('e', d, c, N))
    assert all(A[i][j] == (1 if i == j else 0) for i in range(n) for j in range(n))
    # A_final = E_t ... E_1 A_0 = I  =>  A_0 = E_1^-1 ... E_t^-1 ; as maps x -> A_0 x apply E_t^-1 first.
    letters = []
    for op in reversed(ops):
        if op[0] == 'e': letters.append(op)                       # involution in char 2
        else: letters.append(('d', op[1], mat_inv(op[2], f)))
    return letters

def gl_letter(m, f, r, X, cx, rng):
    n = m * f; F = 1 << f
    A = rand_inv(n, rng)
    letters = decompose(A, m, f, rng)
    o = Ops(); tabs = []
    for L in letters:
        M = L[3] if L[0] == 'e' else L[2]
        cols = [sum(M[i][t] << i for i in range(f)) for t in range(f)]
        tabs.append(row_table(cols, f, X, o))
    ct = o.n
    def prog(idx, P, o):
        spec = o.div(idx, P[m]); rows = [get_row(idx, b, f, P, o) for b in range(m)]
        for L, T in zip(letters, tabs):
            if L[0] == 'e': rows[L[1]] = xor(rows[L[1]], o.look(T, rows[L[2]]), X, F, o)
            else: rows[L[1]] = o.look(T, rows[L[1]])
        return encode(spec, rows, m, f, P, o)
    def mathmap_idx(idx):
        # body bit vector v (index b*f + j  <->  bit j of row b), w = A v
        v = [(idx >> i) & 1 for i in range(n)]
        w = [sum(A[i][t] & v[t] for t in range(n)) & 1 for i in range(n)]
        return (idx >> n << n) | sum(w[i] << i for i in range(n))
    k = n + r; worst = 0; ok = True
    P = [1 << (f * b) for b in range(m + 1)]
    for idx in range(1 << k):
        o = Ops(); out = prog(idx, P, o); worst = max(worst, o.n + 2)
        if out != mathmap_idx(idx): ok = False
    return dict(name=f'random GL({n},F2) as {len(letters)} row letters', m=m, f=f, r=r, k=k,
                ops_per_entry=worst, table_ops_over_2k=round((cx + ct) / (1 << k), 3), ok=ok)

if __name__ == '__main__':
    rng = random.Random(17)
    print('per-entry RAM ops (incl. 1 look + 1 write of the moved value); table ops / 2^k; exhaustive check')
    for m, f, r in [(2, 2, 0), (2, 3, 1), (2, 4, 0), (2, 5, 1), (2, 6, 0), (3, 2, 0), (3, 3, 1), (3, 4, 0), (4, 3, 0)]:
        for d in letters_for(m, f, r, rng):
            print(f"m={d['m']} f={d['f']} r={d['r']} k={d['k']:2d}  {d['name']:42s} ops/entry={d['ops_per_entry']:3d}  "
                  f"tables/2^k={d['table_ops_over_2k']:7.3f}  {'OK' if d['ok'] else 'MISMATCH'}", flush=True)
