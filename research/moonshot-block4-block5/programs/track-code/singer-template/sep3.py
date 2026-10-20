"""sep3.py (mw2-singer4): the order-3 analogue of step 1 (v-ms-fano's sigma = cyclic shift of three address bits).
Entries of Z[i][Z3] with support <= 3 and coefficients in {+-1, +-i}: which are invertible over C (all 3 characters
nonzero) and singular mod 2 (F_2[Z3] = F_2 x F_4: kill F2 = trivial character mod 2, kill F4 = the primitive cube roots
of unity mod 2)?  Exact inverse denominators via the 3x3 circulant (Fractions).  usage: python3 sep3.py"""
import itertools, cmath
from fractions import Fraction as Fr
W3 = cmath.exp(2j * cmath.pi / 3)
UN = [1, -1, 1j, -1j]

def chars(a):
    return [sum(a[k] * W3 ** (s * k) for k in range(3)) for s in range(3)]

def mod2(a):
    return [int(round(z.real) + round(z.imag)) % 2 for z in a]

def kills(a):
    f = mod2(a); out = []
    if sum(f) % 2 == 0:
        out.append('F2')
    # F4 = F2[w]/(w^2+w+1): f(w) = f0 + f1 w + f2 w^2 = (f0 + f2) + (f1 + f2) w
    if (f[0] + f[2]) % 2 == 0 and (f[1] + f[2]) % 2 == 0:
        out.append('F4')
    return tuple(out)

def inv_den(a):
    """exact inverse in Q(i)[Z3] by solving the 6x6 real circulant system; returns lcm of denominators"""
    # represent as real 6x6: unknown b = (re0..2, im0..2), a*b = 1
    n = 3
    M = [[Fr(0)] * 6 for _ in range(6)]
    for i in range(n):
        for j in range(n):
            z = a[(i - j) % n]
            re, im = Fr(int(round(z.real))), Fr(int(round(z.imag)))
            M[i][j] += re; M[i][j + 3] += -im
            M[i + 3][j] += im; M[i + 3][j + 3] += re
    rhs = [Fr(1), Fr(0), Fr(0), Fr(0), Fr(0), Fr(0)]
    A = [row[:] + [rhs[r]] for r, row in enumerate(M)]
    for c in range(6):
        p = next(r for r in range(c, 6) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for r in range(6):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    from math import lcm
    d = 1
    for r in range(6):
        d = lcm(d, A[r][6].denominator)
    return d, [A[r][6] for r in range(3)]

seen = {}
for supp in (1, 2, 3):
    for S in itertools.combinations(range(3), supp):
        for cs in itertools.product(UN, repeat=supp):
            if cs[0] != 1:
                continue
            a = [0, 0, 0]
            for k, c in zip(S, cs):
                a[k] = c
            ch = chars(a)
            inv = all(abs(z) > 1e-9 for z in ch)
            kl = kills(a)
            if inv and kl:
                d, b = inv_den(a)
                key = (supp, kl)
                seen.setdefault(key, []).append((a, d, b))
for key in sorted(seen):
    print('support %d, kills %s: %d entries; e.g.' % (key[0], key[1], len(seen[key])))
    for a, d, b in seen[key][:4]:
        print('   a =', a, ' inverse denominator', d, ' real part of inverse coefficients', [str(x) for x in b])
