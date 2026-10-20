# mytab.py (v-r3-classes): my own GF(2^m) field, quadratic forms and parities, written from the definitions only.
# Moves = nonzero field elements z (1..2^m-1, polynomial basis bits); D_Q = {z : Q(z) = 1}.
# Target T_k = mult(alpha^(2k)); frame = alpha^k f_j for a self-dual basis f_j (Tr(f_i f_j) = delta_ij);
# parity p_k(Q) = sum_j Q(alpha^k f_j) mod 2.
# Output (binary, little-endian uint32): m, nforms, then per form: mask, parity bits over k = 0..2^m-2.
import sys, itertools, struct
POLY = {3: 0b1011, 4: 0b10011, 5: 0b100101}

def run(m, out, alt=False):
    P = POLY[m]; n = 1 << m
    def mul(a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1; a <<= 1
            if a & n: a ^= P
        return r
    def tr(x):
        t, y = 0, x
        for _ in range(m): t ^= y; y = mul(y, y)
        assert t in (0, 1); return t
    # self-dual bases
    bases = []
    for B in itertools.combinations(range(1, n), m):
        if all(tr(mul(B[i], B[j])) == (1 if i == j else 0) for i in range(m) for j in range(m)):
            bases.append(B)
    assert bases
    f = bases[-1] if alt else bases[0]
    al = [1]
    for k in range(1, n - 1): al.append(mul(al[-1], 2))
    assert len(set(al)) == n - 1
    pairs = [(i, j) for i in range(m) for j in range(i + 1, m)]
    nq = 1 << (m + len(pairs))
    def Qv(c, x):
        v = 0
        for i in range(m):
            if (c >> i) & 1 and (x >> i) & 1: v ^= 1
        for t, (i, j) in enumerate(pairs):
            if (c >> (m + t)) & 1 and (x >> i) & 1 and (x >> j) & 1: v ^= 1
        return v
    rows = []
    for c in range(nq):
        val = [Qv(c, x) for x in range(n)]
        mask = sum(1 << x for x in range(1, n) if val[x])
        par = 0
        for k in range(n - 1):
            p = 0
            for j in range(m): p ^= val[mul(al[k], f[j])]
            par |= p << k
        rows.append((mask, par))
    with open(out, 'wb') as fo:
        fo.write(struct.pack('<II', m, nq))
        for mask, par in rows: fo.write(struct.pack('<II', mask, par))
    return rows, bases, f

if __name__ == '__main__':
    m = int(sys.argv[1]); out = sys.argv[2]
    rows, bases, f = run(m, out)
    print('m', m, 'self-dual bases', len(bases), 'used', f)
    if len(bases) > 1:
        rows2, _, f2 = run(m, out + '.alt', alt=True)
        print('alt basis', f2, 'parities identical:', all(a[1] == b[1] for a, b in zip(rows, rows2)))
