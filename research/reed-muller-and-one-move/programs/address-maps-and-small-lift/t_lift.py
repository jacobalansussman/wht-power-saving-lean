"""t_lift.py: exact check of the "boundary lift" (rm-theory Theorem B): a multi-target commutative-slice program whose
array l has target T_l = h_l h_l^T (h_l in H) becomes a unit of C[Sym_m x| H] with the same moves, by relabelling
input l with [h_l] before the program and multiplying output l by [h_l^-1] after it.

Formal group algebra C[Sym_m x| H] with H = <sigma> the Singer cycle (self-dual coordinates, so the trace placement):
elements are dicts {(S, k): coeff}, S a symmetric matrix as a bitmask of its m(m+1)/2 upper entries, k in Z_n;
(S, k)(S', k') = (S + sigma^k S' sigma^kT, k + k').
Programs checked (block 3, n = 7), each as a register machine run on formal inputs:
 P1  2 arrays, targets I and T = sigma sigma^T, each by its own frame (plain multi-target 2 in 6);
 P2  the same as P1 but array 2 is served by the 4-move non-frame decomposition of T (U_3-type 1 + 4 shape);
 P3  a deliberately WRONG lift (array 2 relabelled by [sigma^2] instead of [sigma]) must fail.
"""
import itertools
from gf import GF
from t_field import mons


def run(m=3):
    F = GF(m); n = F.n; q = F.q
    fb = F.self_dual_basis(1)
    co = {x: F.coords(x, fb) for x in range(q)}
    ML = mons(m)
    def sym_bits(S):
        v = 0
        for idx, (i, j) in enumerate(ML):
            if S[i][j]: v |= 1 << idx
        return v
    def bits_sym(b):
        S = [[0] * m for _ in range(m)]
        for idx, (i, j) in enumerate(ML):
            if (b >> idx) & 1: S[i][j] = S[j][i] = 1
        return S
    def mat(a):   # matrix of x -> a x in coordinates (columns = images of f_j)
        cols = [co[F.mul(a, f)] for f in fb]
        return [[cols[j][i] for j in range(m)] for i in range(m)]
    Sig = [mat(F.alpha(k)) for k in range(n)]
    def conj(k, b):   # sigma^k S sigma^kT
        A = Sig[k % n]; S = bits_sym(b)
        R = [[sum(A[i][s] & S[s][t] & A[j][t] for s in range(m) for t in range(m)) & 1 for j in range(m)] for i in range(m)]
        return sym_bits(R)
    def gmul(g, h):
        (S, k), (S2, k2) = g, h
        return (S ^ conj(k, S2), (k + k2) % n)
    def amul(x, y):
        out = {}
        for g, a in x.items():
            for h, b in y.items():
                p = gmul(g, h); out[p] = out.get(p, 0) + a * b
        return {g: c for g, c in out.items() if c != 0}
    def elt(S=0, k=0, c=1): return {(S, k % n): c}
    def zzT(z):
        v = co[z]; return sym_bits([[v[i] & v[j] for j in range(m)] for i in range(m)])
    I_bits = sym_bits([[int(i == j) for j in range(m)] for i in range(m)])
    # a slice program: product of moves (one array at a time; sharing does not matter for the lift check)
    def slice_prog(dirs, x):
        for z in dirs: x = amul(elt(zzT(z)), x)
        return x
    one = elt()
    results = []
    # P1: array 1 target I via frame f_j; array 2 target T = sigma sigma^T via frame alpha f_j
    T2 = conj(1, I_bits)
    dirs1 = list(fb)
    dirs2 = [F.mul(F.alpha(1), f) for f in fb]
    for label, dirs, hk, T in [("P1 array 1 (target I, frame f_j)", dirs1, 0, I_bits),
                               ("P1 array 2 (target sigma sigma^T, frame alpha f_j)", dirs2, 1, T2)]:
        x = elt(0, hk)                       # input relabelled by [h]
        y = slice_prog(dirs, x)              # slice program
        y = amul(elt(0, -hk), y)             # output times [h^-1]
        ok = (y == elt(I_bits, 0))
        # the slice program alone delivers [T]
        alone = slice_prog(dirs, one) == elt(T, 0)
        results.append((label, alone, ok))
    # P2: array 2 via a 4-direction decomposition of T2 that is not a frame: find 4 directions with sum zz^T = T2
    found = None
    for quad in itertools.combinations(range(1, q), 4):
        b = 0
        for z in quad: b ^= zzT(z)
        if b == T2: found = quad; break
    x = elt(0, 1); y = amul(elt(0, -1), slice_prog(list(found), x))
    results.append((f"P2 array 2 via 4 directions {found}", slice_prog(list(found), one) == elt(T2, 0), y == elt(I_bits, 0)))
    # P3: wrong relabelling
    x = elt(0, 2); y = amul(elt(0, -2), slice_prog(dirs2, x))
    results.append(("P3 wrong lift ([sigma^2] for target sigma sigma^T): must FAIL", True, y == elt(I_bits, 0)))
    for lab, alone, ok in results:
        print(f"{lab}: slice program delivers its target: {alone}; lifted word delivers [I]: {ok}")
    return results


if __name__ == "__main__":
    run(3); run(4)
