# Cut (Grassmann-distance) lower bound for in-place words computing T^(x)3, T=[[1,1],[0,1]].
# Registers 0..7 start with x_c (c = 3-bit subset mask).  Register c ends with a scaled target row pi(c),
# row S of Z = T^(x)3 being sum_{U superset of S} x_U.
# LEMMA [proved in NOTE]: for a register set X, only calls (i <- a x_i + b x_j) with i in X, j notin X change
# V_X = span(rows of X), and each moves it to an adjacent subspace (meet of dim |X|-1).  Hence
#     #{calls i<-j : i in X, j notin X}  >=  rank Z[pi(X), complement(X)].
# Here: for every pi (40320), the LP/IP  min sum n_ij  s.t.  in-cut_n(X) >= r_pi(X) for all X.
import itertools, sys, time
import numpy as np
from fractions import Fraction
from scipy.optimize import linprog, milp, LinearConstraint, Bounds

N = 8
Z = [[1 if (S & ~U) == 0 else 0 for U in range(N)] for S in range(N)]

def rank_exact(rows):
    M = [[Fraction(v) for v in r] for r in rows]
    r = 0
    if not M: return 0
    cols = len(M[0])
    for c in range(cols):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r

# rank table: rowmask R (set of subsets), colmask C
RK = [[0] * 256 for _ in range(256)]
for R in range(256):
    rs = [S for S in range(N) if R >> S & 1]
    for C in range(256):
        cs = [U for U in range(N) if C >> U & 1]
        RK[R][C] = rank_exact([[Z[S][U] for U in cs] for S in rs]) if rs and cs else 0

arcs = [(i, j) for i in range(N) for j in range(N) if i != j]   # call: register i <- a x_i + b x_j
Xs = list(range(1, 255))
A = np.zeros((len(Xs), len(arcs)))
for a, X in enumerate(Xs):
    for k, (i, j) in enumerate(arcs):
        if (X >> i & 1) and not (X >> j & 1):
            A[a, k] = 1

def rvec(pi):
    out = []
    for X in Xs:
        R = 0
        for c in range(N):
            if X >> c & 1: R |= 1 << pi[c]
        out.append(RK[R][255 & ~X])
    return np.array(out, dtype=float)

def lp(pi):
    r = rvec(pi)
    res = linprog(np.ones(len(arcs)), A_ub=-A, b_ub=-r, bounds=(0, None), method="highs")
    return res.fun, r

def ip(pi, r=None):
    if r is None: r = rvec(pi)
    res = milp(np.ones(len(arcs)), constraints=LinearConstraint(A, lb=r, ub=np.inf),
               integrality=np.ones(len(arcs)), bounds=Bounds(0, np.inf))
    return res.fun, res.x

if __name__ == "__main__":
    t0 = time.time()
    ident = tuple(range(N))
    v, r = lp(ident)
    print("identity pi: LP =", v, " IP =", ip(ident)[0])
    # all pi, modulo the S3 bit symmetry (conjugation)
    bitperms = list(itertools.permutations(range(3)))
    def act(g, S):
        return sum(((S >> k) & 1) << g[k] for k in range(3))
    seen = set(); reps = []
    for pi in itertools.permutations(range(N)):
        if pi in seen: continue
        orb = set()
        for g in bitperms:
            gi = [act(g, S) for S in range(N)]
            ginv = [0] * N
            for S in range(N): ginv[gi[S]] = S
            q = tuple(gi[pi[ginv[c]]] for c in range(N))
            orb.add(q)
        seen |= orb; reps.append(pi)
    print("pi classes:", len(reps), flush=True)
    hist = {}; low = []
    for pi in reps:
        v, r = lp(pi)
        key = round(v, 6)
        hist[key] = hist.get(key, 0) + 1
        if v < 11 + 1e-6: low.append((v, pi))
    print("LP value histogram:", sorted(hist.items()))
    print("classes with LP <= 11:", len(low), " time", round(time.time() - t0, 1), flush=True)
    low.sort()
    ipl = {}
    worst = []
    for v, pi in low:
        w, x = ip(pi)
        ipl[round(w)] = ipl.get(round(w), 0) + 1
        if w < 11.5: worst.append((w, pi, x))
    print("IP histogram over LP<=11 classes:", sorted(ipl.items()))
    print("classes with IP <= 11:", len(worst))
    for w, pi, x in worst[:20]:
        used = [(arcs[k], int(round(x[k]))) for k in range(len(arcs)) if x[k] > 0.5]
        print("  IP", w, "pi", pi, "arcs", used)
    print("time", round(time.time() - t0, 1))
