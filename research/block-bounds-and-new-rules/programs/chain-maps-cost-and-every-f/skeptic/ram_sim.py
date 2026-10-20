# xc-v-legit (wave 17 skeptic): independent re-check of the RAM cost of cross-column letters.
# Written from scratch (not a copy of xc-legit/ram_letters.py).  Semantics follow UniformFourier.lean namespace RAM:
#   * integer primitives add/sub(truncated)/mul/div/mod/lt and `look` cost 1 each;
#   * tapes are IMMUTABLE: a tape is only built by `tab len body` or `sow len count body`, whose body sees the
#     environment and the index but NOT the tape under construction (Code.tab / Code.sow in UniformFourier.lean).
#     Bill.tab/sow charge 1 + len plus, per write, the body's work and O(1) step overhead (counted here as 2);
#   * a table that needs its own earlier entries (dynamic programming) is therefore built in STAGES: stage j is a
#     fresh tab that reads stage j-1 (a `loop` whose accumulator is a tape).  This is stricter than xc-legit, whose
#     xor_table reads the array it is filling.
# For each letter: max body work per entry (must not grow with f), table work / 2^k (must stay bounded), and an
# exhaustive check against a bit-level reference on the column representation.
# Usage: python3 -I ram_sim.py      (largest: m=2 f=7 r=1, k=15, 32768 entries; < 60 MB)
import random, sys

class RAM:
    def __init__(s): s.n = 0
    def add(s, a, b): s.n += 1; return a + b
    def sub(s, a, b): s.n += 1; return a - b if a >= b else 0
    def mul(s, a, b): s.n += 1; return a * b
    def div(s, a, b): s.n += 1; return a // b if b else 0
    def mod(s, a, b): s.n += 1; return a % b if b else a
    def lt(s, a, b):  s.n += 1; return 1 if a < b else 0
    def look(s, t, i):
        s.n += 1
        assert 0 <= i < len(t), (i, len(t))
        return t[i]
    def neg(s, v):    s.n += 1; return -v          # Atom.sub from Atom.cz: 2 ops really; constant anyway
    def tab(s, length, body):
        s.n += 1 + length; out = []
        worst = 0
        for i in range(length):
            n0 = s.n; out.append(body(i)); s.n += 2; worst = max(worst, s.n - n0)
        s.last_worst = worst
        return tuple(out)
    def sow(s, length, count, body, blank=0):
        s.n += 1 + length; arr = [blank] * length; worst = 0
        for i in range(count):
            n0 = s.n; pos, val = body(i); s.n += 2; worst = max(worst, s.n - n0)
            assert 0 <= pos < length; arr[pos] = val
        s.last_worst = worst
        return tuple(arr)

# ---------------- tables, built in stages ----------------
def stage_powers(R, f):
    """2^0..2^f as loop accumulators (O(f) once); returned as a plain list for the driver (each value is a register)."""
    p, out = 1, [1]
    for _ in range(f): p = R.mul(p, 2); out.append(p)
    return out

def xor_table(R, f, pw):
    X = R.tab(1, lambda i: 0)                       # X_0: 0 xor 0
    for j in range(1, f + 1):
        Pj, Pj1 = pw[j], pw[j - 1]
        prev = X
        def body(i, prev=prev, Pj=Pj, Pj1=Pj1):
            a = R.div(i, Pj); b = R.mod(i, Pj)
            hi = R.look(prev, R.add(R.mul(R.div(a, 2), Pj1), R.div(b, 2)))
            lo = R.mod(R.add(R.mod(a, 2), R.mod(b, 2)), 2)
            return R.add(R.mul(hi, 2), lo)
        X = R.tab(Pj * Pj, body)
    return X

def prefix_xor_table(R, f, pw):
    """PX[y]: bit j = y_0 ^ ... ^ y_j  (the inverse of row ^= row<<1 on one row)."""
    T = R.tab(1, lambda i: 0)
    for j in range(1, f + 1):
        Pj1 = pw[j - 1]; Pj2 = pw[j - 2] if j >= 2 else 0; prev = T
        def body(y, prev=prev, Pj1=Pj1, Pj2=Pj2, j=j):
            lo = R.mod(y, Pj1); p = R.look(prev, lo)
            top = R.mod(R.div(p, Pj2), 2) if j >= 2 else 0
            bit = R.mod(R.add(top, R.div(y, Pj1)), 2)
            return R.add(p, R.mul(bit, Pj1))
        T = R.tab(pw[j], body)
    return T

def chain_parity_table(R, f, pw):
    """Q[y] = sum_j y_j y_{j+1} mod 2  (standard-basis CZ chain along one row)."""
    T = R.tab(1, lambda i: 0)
    for j in range(1, f + 1):
        Pj1 = pw[j - 1]; Pj2 = pw[j - 2] if j >= 2 else 0; prev = T
        def body(y, prev=prev, Pj1=Pj1, Pj2=Pj2, j=j):
            lo = R.mod(y, Pj1); q = R.look(prev, lo)
            if j >= 2: q = R.add(q, R.mul(R.div(y, Pj1), R.mod(R.div(y, Pj2), 2)))
            return R.mod(q, 2)
        T = R.tab(pw[j], body)
    return T

# ---------------- reference (bit level, column picture) ----------------
def rows_of(idx, m, f): return [(idx >> (f * b)) & ((1 << f) - 1) for b in range(m)], idx >> (m * f)
def idx_of(rows, spec, m, f):
    v = spec << (m * f)
    for b, r in enumerate(rows): v |= r << (f * b)
    return v
def cols_of(rows, m, f): return [sum(((rows[b] >> j) & 1) << b for b in range(m)) for j in range(f)]
def rows_from_cols(cols, m, f): return [sum(((cols[j] >> b) & 1) << j for j in range(f)) for b in range(m)]

def ref_chain(idx, m, f):        # column j += column j-1 (j >= 1)
    rows, spec = rows_of(idx, m, f); c = cols_of(rows, m, f)
    c2 = [c[j] ^ (c[j - 1] if j else 0) for j in range(f)]
    return idx_of(rows_from_cols(c2, m, f), spec, m, f)
def ref_cross(idx, m, f):        # bit 0 of column j += bit 1 of column j-1
    rows, spec = rows_of(idx, m, f); c = cols_of(rows, m, f)
    c2 = [c[j] ^ (((c[j - 1] >> 1) & 1) if j else 0) for j in range(f)]
    return idx_of(rows_from_cols(c2, m, f), spec, m, f)

def run(m, f, r, rng):
    k = m * f + r; N = 1 << k; F = 1 << f
    R = RAM()
    pw = stage_powers(R, max(f, 1) * (m + 1) + r)
    PB = [pw[f * b] for b in range(m + 1)]
    n0 = R.n
    X = xor_table(R, f, pw); nX = R.n - n0
    n0 = R.n; PX = prefix_xor_table(R, f, pw); nPX = R.n - n0
    n0 = R.n; Q = chain_parity_table(R, f, pw); nQ = R.n - n0
    data = tuple(rng.randrange(-99, 100) for _ in range(N))
    res = []

    def rows_get(i):
        return R.div(i, PB[m]), [R.mod(R.div(i, PB[b]), F) for b in range(m)]
    def enc(spec, rows):
        v = R.mul(spec, PB[m])
        for b in range(m): v = R.add(v, R.mul(rows[b], PB[b]))
        return v

    # L1 forward: row ^= row<<1 on every row, as a sow (each source entry writes its value at its new address)
    def b_chain(i):
        spec, rows = rows_get(i)
        new = [R.look(X, R.add(R.mul(rw, F), R.mod(R.mul(rw, 2), F))) for rw in rows]
        return enc(spec, new), R.look(data, i)
    n0 = R.n; out1 = R.sow(N, N, b_chain); w1 = R.last_worst; c1 = R.n - n0
    ok1 = all(out1[ref_chain(i, m, f)] == data[i] for i in range(N))
    res.append(('row ^= row<<1 (all rows), sow', w1, c1 / N, nX / N, ok1))

    # the SAME letter as a tab: each output entry i fetches from its source A^-1(i) = prefix XOR of its rows
    # (this is where the inverse map's table is used); must equal the sow result.  Then P^-1 as a tab
    # (fetch from A(i) with the XOR table) must give the input back.
    def b_tabP(i):
        spec, rows = rows_get(i)
        return R.look(data, enc(spec, [R.look(PX, rw) for rw in rows]))
    n0 = R.n; viaTab = R.tab(N, b_tabP); w2 = R.last_worst; c2 = R.n - n0
    def b_back(i):
        spec, rows = rows_get(i)
        return R.look(viaTab, enc(spec, [R.look(X, R.add(R.mul(rw, F), R.mod(R.mul(rw, 2), F))) for rw in rows]))
    back = R.tab(N, b_back)
    ok2 = (viaTab == out1) and (back == data)
    res.append(('same letter as tab via prefix-XOR table', w2, c2 / N, (nX + nPX) / N, ok2))

    # L2: row0 ^= row1<<1 (an involution)
    if m >= 2:
        def b_cross(i):
            spec, rows = rows_get(i)
            rows[0] = R.look(X, R.add(R.mul(rows[0], F), R.mod(R.mul(rows[1], 2), F)))
            return enc(spec, rows), R.look(data, i)
        n0 = R.n; out3 = R.sow(N, N, b_cross); w3 = R.last_worst; c3 = R.n - n0
        ok3 = all(out3[ref_cross(i, m, f)] == data[i] for i in range(N))
        res.append(('row0 ^= row1<<1, sow', w3, c3 / N, nX / N, ok3))

    # L3: phase (-1)^{sum_j row0_j row0_{j+1}} (standard-basis CZ chain on row 0), a tab
    def b_phase(i):
        q = R.look(Q, R.mod(i, F)); v = R.look(data, i)
        R.n += 1                                     # the ifz test
        return v if q == 0 else R.neg(v)
    n0 = R.n; out4 = R.tab(N, b_phase); w4 = R.last_worst; c4 = R.n - n0
    def refq(y): return bin(y & (y >> 1)).count('1') % 2
    ok4 = all(out4[i] == (data[i] if refq(i % F) == 0 else -data[i]) for i in range(N))
    res.append(('phase (-1)^CZchain(row0), tab', w4, c4 / N, nQ / N, ok4))
    return k, res

rng = random.Random(20261013)
cases = [(2, f, r) for f in range(1, 8) for r in (0, 1)] + [(3, f, r) for f in range(1, 5) for r in (0, 2)] + \
        [(4, f, 0) for f in range(1, 4)]
allok = True
print("m f r  k  letter                                  max work/entry  pass work/2^k  table work/2^k  exhaustive")
for (m, f, r) in cases:
    k, res = run(m, f, r, rng)
    for (name, w, c, t, ok) in res:
        allok &= ok
        print(f"{m} {f} {r} {k:2d}  {name:40s} {w:5d}           {c:7.2f}        {t:7.3f}        {'OK' if ok else 'FAIL'}")
print("ALL OK" if allok else "SOME FAILED")
