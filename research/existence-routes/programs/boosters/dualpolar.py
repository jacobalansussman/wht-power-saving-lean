# e3-boost dualpolar.py: exact intersection numbers p^d_ij of the frame graph at block M (the dual polar graph of
# Sp(2M, 2): frames = Lagrangian subspaces, distance = M - dim of the intersection), and the DISTANCE-SET certificates
#   V = Q^(frames), Y_U = span{ e_W : dist(U, W) in S }       (S any subset of {0..M})
#   f(rank-d block) = |S(U)| - |S(U) cap S(V)| for dist(U, V) = d  =  sum_{i in S} k_i - sum_{i, j in S} p^d_ij.
# Valid for Theorem R [derived here]: every engine letter maps frames to frames preserving distance, the free letters
# fix L_0, and in a distance-regular graph |S(U) cap S(V)| depends only on dist(U, V).  S = {0} is the frame count;
# S = {0, 1} at M = 3 is certificate E of w-bounds 5.4 (12, 12, 15).
# Intersection array of the dual polar graph C_M(2) (Brouwer-Cohen-Neumaier 9.4): b_i = 2^(i+1) (2^(M-i) - 1),
# c_i = 2^i - 1.  Checked below against the sphere sizes 15, 135, 2295, 75735, 4922775 (M = 2..6).
from fractions import Fraction as F

def array(M):
    b = [2 ** (i + 1) * (2 ** (M - i) - 1) for i in range(M)] + [0]
    c = [0] + [2 ** i - 1 for i in range(1, M + 1)]
    k = b[0]
    a = [k - b[i] - c[i] for i in range(M + 1)]
    return a, b, c

def spheres(M):
    a, b, c = array(M)
    ks = [1]
    for i in range(1, M + 1): ks.append(ks[-1] * b[i - 1] // c[i])
    return ks

def polymul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q): r[i + j] += x * y
    return r

def intersection_numbers(M):
    """p[d][i][j] = #{W : dist(U,W) = i, dist(V,W) = j} for dist(U,V) = d."""
    a, b, c = array(M)
    # v_i(x): A_i = v_i(A_1);  x v_i = b_{i-1} v_{i-1} + a_i v_i + c_{i+1} v_{i+1}
    v = [[F(1)], [F(0), F(1)]]
    for i in range(1, M):
        nxt = polymul([F(0), F(1)], v[i])
        nxt = [x - a[i] * (v[i][t] if t < len(v[i]) else 0) - b[i - 1] * (v[i - 1][t] if t < len(v[i - 1]) else 0)
               for t, x in enumerate(nxt)]
        v.append([x / c[i + 1] for x in nxt])
    # minimal polynomial: x v_M - a_M v_M - b_{M-1} v_{M-1}
    mp = polymul([F(0), F(1)], v[M])
    mp = [x - a[M] * (v[M][t] if t < len(v[M]) else 0) - b[M - 1] * (v[M - 1][t] if t < len(v[M - 1]) else 0)
          for t, x in enumerate(mp)]
    def reduce(p):
        p = list(p)
        while len(p) > M + 1:
            lead = p[-1] / mp[-1]; s = len(p) - len(mp)
            for t in range(len(mp)): p[s + t] -= lead * mp[t]
            p.pop()
        return p + [F(0)] * (M + 1 - len(p))
    def coords(p):  # express p (deg <= M) in the basis v_0..v_M (v_i has degree i, leading coeff 1/(c_2...c_i))
        p = list(p); out = [F(0)] * (M + 1)
        for i in range(M, -1, -1):
            lead = p[i] / v[i][i]; out[i] = lead
            for t in range(i + 1): p[t] -= lead * v[i][t]
        return out
    ks = spheres(M)
    P = [[[0] * (M + 1) for _ in range(M + 1)] for _ in range(M + 1)]
    for i in range(M + 1):
        for j in range(M + 1):
            co = coords(reduce(polymul(v[i], v[j])))   # A_i A_j = sum_d p^d_ij A_d
            for d in range(M + 1):
                assert co[d].denominator == 1
                P[d][i][j] = int(co[d])
    return P, ks

def dist_cert(M, S, P=None, ks=None):
    if P is None: P, ks = intersection_numbers(M)
    tot = sum(ks[i] for i in S)
    return [tot - sum(P[d][i][j] for i in S for j in S) for d in range(1, M + 1)]

if __name__ == "__main__":
    for M in range(2, 7):
        print(M, spheres(M), sum(spheres(M)))
    P, ks = intersection_numbers(3)
    print("M=3 S={0}:", dist_cert(3, [0], P, ks), " S={0,1} (E: 12,12,15):", dist_cert(3, [0, 1], P, ks))
    # sanity: sum_j p^d_ij = k_i
    for M in range(2, 9):
        P, ks = intersection_numbers(M)
        ok = all(sum(P[d][i][j] for j in range(M + 1)) == ks[i] for d in range(M + 1) for i in range(M + 1))
        print(f"M={M}: row sums = sphere sizes: {ok}; p^1_11 = a_1 = {P[1][1][1]}")
