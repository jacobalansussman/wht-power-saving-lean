# ut-normal relcheck.py (wave 17): exact checks behind PROOF.md section 5 (m = 2, f = 1..6).
# Operators are stored as monomials: (perm[x], phase-exponent mod 4 at x)  meaning  e_x -> i^ph[x] e_perm[x].
# Everything is exact (integers mod 4); no floats, memory < 10 MB.
import itertools, sys
m = 2
def run(f):
    k = m * f; N = 1 << k
    mask = (1 << f) - 1
    def rows(x): return [sum(((x >> (m*j+i)) & 1) << j for j in range(f)) for i in range(m)]
    def frm(r): return sum(((r[i] >> j) & 1) << (m*j+i) for i in range(m) for j in range(f))
    def pmul(p, q):
        out = 0
        for j in range(f):
            if (p >> j) & 1: out ^= q << j
        return out & mask
    pc = lambda v: bin(v).count('1')
    I = (list(range(N)), [0]*N)
    def mono_perm(fun): return ([fun(x) for x in range(N)], [0]*N)
    def mono_diag(fun): return (list(range(N)), [fun(x) % 4 for x in range(N)])
    def mul(A, B):          # (A B) e_x = A (i^b e_{pb}) = i^(b + a[pb]) e_{pa[pb]}
        pa, aa = A; pb, bb = B
        return ([pa[pb[x]] for x in range(N)], [(bb[x] + aa[pb[x]]) % 4 for x in range(N)])
    def inv(A):
        pa, aa = A; p = [0]*N; a = [0]*N
        for x in range(N): p[pa[x]] = x; a[pa[x]] = (-aa[x]) % 4
        return (p, a)
    def G(Gm):
        return mono_perm(lambda x: frm([ (lambda r: __import__('functools').reduce(lambda s, l: s ^ pmul(Gm[i][l], r[l]), range(m), 0))(rows(x)) for i in range(m)]))
    chi = G([[0b11, 0], [0, 0b11]])            # row ^= row << 1 on both rows
    E   = G([[1, 0b10], [0, 1]])               # row0 ^= row1 << 1
    Dcz = mono_diag(lambda x: 2 * pc(rows(x)[0] & rows(x)[1]))     # uniform CZ between the two rows
    D1  = mono_diag(lambda x: pc(rows(x)[1]))                      # uniform S on row 1
    A = mul(mul(mul(inv(E), Dcz), E), inv(Dcz))
    B = mul(mul(mul(inv(chi), D1), chi), inv(mul(D1, D1)))
    wS = mul(A, inv(B))
    Slast = mono_diag(lambda x: (rows(x)[1] >> (f-1)) & 1)
    P = mul(mul(chi, Slast), inv(chi))
    Ppar = mono_diag(lambda x: pc(rows(x)[1]) & 1)                 # i^(parity of row 1)
    Zstr = mono_diag(lambda x: 2 * pc(rows(x)[1]))                 # Z on row 1 in every column (uniform)
    ok1 = wS == Slast
    ok2 = P == Ppar
    ok3 = mul(Ppar, Ppar) == Zstr
    # relation P = a I + b Z_str with a = (1+i)/2, b = (1-i)/2: check pointwise on the diagonal (all diagonal ops)
    # i^p = a + b (-1)^p for p in {0,1}: p=0 -> a+b = 1 ; p=1 -> a-b = i.  Holds for every x since parity in {0,1}.
    ok4 = all(Ppar[1][x] in (0, 1) and Zstr[1][x] == 2 * Ppar[1][x] % 4 for x in range(N))
    # normalized trace of P: (1/N) sum i^parity = (1+i)/2
    cnt = [Ppar[1].count(e) for e in range(4)]
    tr = (cnt[0] - cnt[2], cnt[1] - cnt[3])   # real, imag parts times N
    ok5 = (2 * tr[0] == N and 2 * tr[1] == N)
    # f = 1 collapse: chi and E are the identity
    ok6 = (f != 1) or (chi == I and E == I)
    return ok1, ok2, ok3, ok4, ok5, ok6
for f in range(1, 7):
    r = run(f)
    print(f"f={f}: wS = S_last {r[0]};  chi S_last chi^-1 = i^parity(row1) {r[1]};  (i^parity)^2 = uniform Z string {r[2]};"
          f"  i^parity = a + b Zstr pointwise {r[3]};  tau(i^parity) = (1+i)/2 {r[4]};  f=1: chi = E = I {r[5]}")
# No-go for the f = 1 collapse composed with a tensor power: theta_t(S) versus a + b theta_t(Z), on t copies of one bit.
# theta_t(S) = diag i^{wt(y)}, a + b theta_t(Z) = diag (a + b (-1)^{wt(y)}) = diag i^{wt(y) mod 2}.
for t in range(1, 6):
    diffs = sum(1 for y in range(1 << t) if (bin(y).count('1') % 4) != (bin(y).count('1') % 2))
    print(f"t={t}: entries where theta_t(S) != a + b theta_t(Z): {diffs} of {1 << t}  -> collapse {'IS' if diffs == 0 else 'is NOT'} a representation of the relation")
