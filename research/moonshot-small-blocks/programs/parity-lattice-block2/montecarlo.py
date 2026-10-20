"""montecarlo.py -- a fourth, completely different check of lambda(K) = 10 and lambda(D) = 6 on the parity lattice at f = 2:
lambda(g) = log2 [Lam : Lam cap g^-1 Lam], and for v uniform in Lam (mod a deep sublattice) P(g v in Lam) = 2^-lambda(g).
Membership in Lam_par is tested from its DEFINITION (entries congruent mod pi, sum = 0 mod 2), not from any basis."""
import random, math
rng = random.Random(2026)
N = 16
def gmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def vpi_ge(x, k):
    """is the Gaussian integer x divisible by pi^k ?"""
    a, b = x
    for _ in range(k):
        if (a - b) % 2: return False
        a, b = (a - b) // 2, (a + b) // 2
    return True
def divpi(x, k):
    a, b = x
    for _ in range(k): a, b = (a - b) // 2, (a + b) // 2
    return (a, b)
def in_par(v):
    """v: list of Gaussian integers. Lam_par membership by definition."""
    v0 = v[0]
    if not all(vpi_ge((x[0] - v0[0], x[1] - v0[1]), 1) for x in v): return False
    s = (sum(x[0] for x in v), sum(x[1] for x in v))
    return s[0] % 2 == 0 and s[1] % 2 == 0
def one_minus_iX(v, mask):
    return [(v[x][0] + v[x ^ mask][1], v[x][1] - v[x ^ mask][0]) for x in range(N)]   # v - i v[x^mask]
def apply_scaled(v, masks):
    for m in masks: v = one_minus_iX(v, m)
    return v
def in_par_after(v, masks, s):
    """is pi^-s prod(1 - i X_m) v in Lam_par ?"""
    w = apply_scaled(v, masks)
    if not all(vpi_ge(x, s) for x in w): return False
    return in_par([divpi(x, s) for x in w])
def rand_par():
    """uniform-ish random element of Lam_par modulo pi^8 L: a*1 + pi*w + pi^2*z with w of even sum (mod pi^?)"""
    a = (rng.randrange(256), rng.randrange(256))
    w = [(rng.randrange(256), rng.randrange(256)) for _ in range(N)]
    # force sum w = 0 mod pi: adjust w[0]
    s = (sum(x[0] for x in w), sum(x[1] for x in w))
    if (s[0] - s[1]) % 2: w[0] = (w[0][0] + 1, w[0][1])
    z = [(rng.randrange(256), rng.randrange(256)) for _ in range(N)]
    pi = (1, -1); pi2 = (0, -2)
    v = [(a[0] + gmul(pi, w[x])[0] + gmul(pi2, z[x])[0], a[1] + gmul(pi, w[x])[1] + gmul(pi2, z[x])[1]) for x in range(N)]
    assert in_par(v)
    return v
T = 400000
K_masks = [1 << b for b in range(4)]
D_masks = {z: [z, z << 2] for z in (1, 2, 3)}
hitK = 0; hitD = {1: 0, 2: 0, 3: 0}
for t in range(T):
    v = rand_par()
    if in_par_after(v, K_masks, 4): hitK += 1
    for z in (1, 2, 3):
        if in_par_after(v, D_masks[z], 2): hitD[z] += 1
print("samples", T)
print("K : fraction %.6f -> lambda ~ %.3f (exact 10 => %.6f)" % (hitK / T, -math.log2(hitK / T), 2 ** -10))
for z in (1, 2, 3):
    print("D_%d: fraction %.6f -> lambda ~ %.3f (exact 6 => %.6f)" % (z, hitD[z] / T, -math.log2(hitD[z] / T), 2 ** -6))
