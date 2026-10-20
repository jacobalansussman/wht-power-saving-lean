"""r-swaps/hubcheck.py -- the two failing conditions of the permutation door, on the hub's real bit geometry
(#589 nine-stage GEOMETRY.json, 111 ports, h = 8, a = 7, rho15 rational; copy in ../s-fix/, sha256 checked).
Exact arithmetic modulo the prime 2^61-1 (no rho15 denominator is divisible by it: the hub's denominators are tiny).
  (1) frame part: for orthogonal pairs (S, T), the S-local end frame K_T (+) E vs K_S (+) E, mapped to global frames by
      rho_S: rank distance (the return leg a swapped 0-bank content must pay).
  (2) route part: boundary-transfer condition rho(1,S) rho(1,T)^-1 rho(2,T) == rho(2,S) on the window U (+) E1 (+) E2
      (a content of S crossing the 1 -> 2 boundary inside T's owner lands on its own route).
usage: python3 -I hubcheck.py [pairs=60]"""
import json, sys, hashlib, random
from fractions import Fraction as F
P = (1 << 61) - 1
path = '../../data/hub-inputs/GEOMETRY-nine-bit589.json'
sha = hashlib.sha256(open(path, 'rb').read()).hexdigest()
assert sha == '0058c491d839142aa99ef86a5ec697ce99857632228531512ea3fa74fb24d734', sha
g = json.load(open(path)); ports = g['ports']; V = len(ports); h = g['h']; a = h - 1; n = h + 2 * a
def modp(x):
    x = F(x); return x.numerator % P * pow(x.denominator % P, P - 2, P) % P
def mm(A, B): return [[sum(A[i][k] * B[k][j] for k in range(len(B))) % P for j in range(len(B[0]))] for i in range(len(A))]
def eye(k): return [[int(i == j) for j in range(k)] for i in range(k)]
def rank(rows):
    M = [r[:] for r in rows]; rk = 0
    for c in range(len(M[0]) if M else 0):
        p = next((r for r in range(rk, len(M)) if M[r][c] % P), None)
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]; iv = pow(M[rk][c], P - 2, P); M[rk] = [x * iv % P for x in M[rk]]
        for r in range(len(M)):
            if r != rk and M[r][c] % P:
                f = M[r][c]; M[r] = [(x - f * y) % P for x, y in zip(M[r], M[rk])]
        rk += 1
    return rk
def route(t, s):                       # rho15 placed on U (+) E_s inside the window U (+) E1 (+) E2
    R15 = [[modp(x) for x in row] for row in ports[t]['rho15']]; ax = list(range(h)) + list(range(h + a * (s - 1), h + a * s))
    R = eye(n)
    for i, r in enumerate(ax):
        for j, c in enumerate(ax): R[r][c] = R15[i][j]
    return R
q = [[F(x) for x in p['q']] for p in ports]
dot = lambda u, v: sum(x * y for x, y in zip(u, v))
pairs = [(s, t) for s in range(V) for t in range(V) if s != t and dot(q[s], q[t]) == 0]
random.seed(1); N = int(sys.argv[1]) if len(sys.argv) > 1 else 60; sample = random.sample(pairs, N)
def perp_cols(t):                      # basis of q_t^perp in U (standard dot, as the hub's projector identities use)
    qq = [modp(x) for x in ports[t]['q']]; j = next(k for k in range(h) if qq[k]); out = []
    for k in range(h):
        if k == j: continue
        v = [0] * h; v[k] = 1; v[j] = (-qq[k] * pow(qq[j], P - 2, P)) % P; out.append(v)
    return out
def local_K_plus_E(t, s):              # columns of K_t (+) E_s in window coordinates
    cols = [v + [0] * (2 * a) for v in perp_cols(t)]
    for j in range(h + a * (s - 1), h + a * s): cols.append([int(k == j) for k in range(n)])
    return cols
def apply(R, cols): return [[sum(R[i][k] * c[k] for k in range(n)) % P for i in range(n)] for c in cols]
def dist(A, B): return rank(A) + rank(B) - 2 * (rank(A) + rank(B) - rank(A + B))
d_hist = {}; cond = 0; same_pair_ctrl = 0
for S, T in sample:
    R1S = route(S, 1)
    own = apply(R1S, local_K_plus_E(S, 1)); foreign = apply(R1S, local_K_plus_E(T, 1))
    d = dist(own, foreign); d_hist[d] = d_hist.get(d, 0) + 1
    lhs = mm(mm(R1S, route(T, 1)), route(T, 2))        # rho^-1 = rho (hub: rho^2 = I)
    cond += lhs == route(S, 2)
S0 = sample[0][0]
same_pair_ctrl = mm(mm(route(S0, 1), route(S0, 1)), route(S0, 2)) == route(S0, 2)
print('geometry sha256 ok; ports %d, orthogonal ordered pairs %d, sampled %d' % (V, len(pairs), N))
print('(1) global distance between a swapped 0-bank content\'s end frame (S-local K_T) and its own (S-local K_S): %s'
      % d_hist)
print('(2) boundary-transfer condition rho(1,S) rho(1,T) rho(2,T) == rho(2,S): %d of %d sampled pairs '
      '(control S = T: %s)' % (cond, N, same_pair_ctrl))
