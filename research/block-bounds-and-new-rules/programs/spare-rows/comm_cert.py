# comm_cert.py (d-G3-SPECTATOR/more-spectators, wave 17, 13 Oct 2026).
# The Walsh-diagonal ("commutative", OpenAI-frame) game and a PERFECT char-0 certificate for rank-2 and rank-3 targets
# B_k (x) I_(m-k), any number of spectator rows.
#
# Facts checked here (exact integer / Gaussian-integer arithmetic, no floats in the decisive checks):
#  (1) dir(z) = J diag(i^[z.x]) J  (J = normalised Walsh), so every word in dir(z), shifts R_w and the scalar i is
#      J diag(i^q) J with q : F_2^m -> Z/4; q -> B_q (B_q(x,y) = (q(x+y)-q(x)-q(y)+q(0))/2 mod 2, a symmetric F_2
#      form) is a homomorphism onto Sym_m(F_2) with dir(z) -> z z^T, R_w -> 0, i -> 0.  B_k (x) I -> M_K = sum_{j<k} e_j e_j^T.
#  (2) phi = (Q_1..Q_k): Sym_m -> F_2^k linear, phi(z z^T) = (Q_j(z)), with Q_j quadratic, pairwise disjoint supports,
#      sum_{j<k} Q_j(e_i) ... = 1 for every j.  Then phi(zz^T) is 0 or a unit vector and phi(M_K) = (1,..,1).
#  (3) Certificate: V = (+)_j Q[F_2^k] (G acts through phi by translation), Y_0 = (+)_j span{delta_x : x_j = 0}.
#      f(g) = dim(Y_0 + gY_0) - dim Y_0 = 2^(k-1) * (Hamming weight of phi(g)).  Rank-1 move: <= 2^(k-1); target: k 2^(k-1).
#      Theorem R (wht35/wider/w-bounds.md section 2) then gives  moves >= k W.
#  (4) The same phi is kept by address maps A whose active rows are a permutation of the active coordinates
#      ((Az)_j = z_pi(j) for j < k; spectator rows arbitrary): checked on all such A for small m.
# Usage: python3 -I comm_cert.py
import itertools, random, sys
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import numpy as np

def popc(x): return bin(x).count('1')

# ---------- (1) Walsh-diagonal form of dir(z), exact -------------------------------------------------------------
def check_dir_walsh(m):
    n = 2 ** m
    H = np.array([[(-1) ** popc(s & x) for x in range(n)] for s in range(n)], dtype=np.int64)
    ok = True
    for z in range(1, n):
        # 2 dir(z) = (1+i) I + (1-i) R_z ; represent Gaussian integers as pairs (re, im)
        Re = np.eye(n, dtype=np.int64); Im = np.eye(n, dtype=np.int64)
        R = np.zeros((n, n), dtype=np.int64)
        for x in range(n): R[x ^ z, x] = 1
        Re = Re + R; Im = Im - R
        # H (2 dir) H = 2 * n * diag(i^[z.s])
        A = H @ Re @ H; B = H @ Im @ H
        for s in range(n):
            e = popc(z & s) & 1
            want = (2 * n, 0) if e == 0 else (0, 2 * n)
            if (A[s, s], B[s, s]) != want: ok = False
        offd = A - np.diag(np.diag(A)), B - np.diag(np.diag(B))
        if np.any(offd[0]) or np.any(offd[1]): ok = False
    return ok

# ---------- Z/4 exponents and their bilinear forms -----------------------------------------------------------------
def q_dir(z, m):  return [popc(z & x) & 1 for x in range(2 ** m)]          # [z.x] as integer 0/1
def q_shift(w, m): return [2 * (popc(w & x) & 1) for x in range(2 ** m)]  # (-1)^{w.x} = i^{2 w.x}
def q_add(q1, q2): return [(a + b) % 4 for a, b in zip(q1, q2)]
def bil(q, m):
    n = 2 ** m; M = np.zeros((m, m), dtype=np.int64)
    for a in range(m):
        for b in range(m):
            x, y = 1 << a, 1 << b
            v = (q[x ^ y] - q[x] - q[y] + q[0]) % 4
            assert v % 2 == 0
            M[a, b] = (v // 2) % 2
    return M
def bil_full_ok(q, m):
    """B_q bilinear on ALL pairs (sanity: q is a Z/4 quadratic form)"""
    n = 2 ** m; M = bil(q, m)
    for x in range(n):
        for y in range(n):
            v = (q[x ^ y] - q[x] - q[y] + q[0]) % 4
            bx = [(x >> a) & 1 for a in range(m)]; by = [(y >> a) & 1 for a in range(m)]
            if v != 2 * (int(np.array(bx) @ M @ np.array(by)) % 2): return False
    return True

# ---------- (2) the quadratic forms ---------------------------------------------------------------------------------
def Q_list(k):
    """k = 2: Q1 = z1(1+z2), Q2 = z2(1+z1)  (supports {10}, {01} on the active bits).
    k = 3: Q_j = z_j + z_j z_{j+1} + z_j z_{j+2} + z_{j+1} z_{j+2}  (support {e_j, e_{j+1}+e_{j+2}} on the active bits;
    the three supports are disjoint and permuted by EVERY permutation of the active bits)."""
    if k == 2: return [lambda z: (z & 1) * (1 - ((z >> 1) & 1)), lambda z: ((z >> 1) & 1) * (1 - (z & 1))]
    if k == 3:
        def mk(j):
            a, b, c = j, (j + 1) % 3, (j + 2) % 3
            return lambda z: (((z >> a) & 1) + ((z >> a) & 1) * ((z >> b) & 1) + ((z >> a) & 1) * ((z >> c) & 1)
                              + ((z >> b) & 1) * ((z >> c) & 1)) % 2
        return [mk(j) for j in range(3)]
    raise ValueError
def phi_of_sym(M, k):
    """linear functional version: lambda_j(M) = sum_i a_ji M_ii + sum_{i<l} b_jil M_il, read off Q_j"""
    Qs = Q_list(k); out = []
    m = M.shape[0]
    for Q in Qs:
        v = 0
        for i in range(m):
            a = Q(1 << i)
            v ^= a & int(M[i, i])
            for l in range(i + 1, m):
                b = (Q((1 << i) | (1 << l)) - Q(1 << i) - Q(1 << l)) % 2
                v ^= b & int(M[i, l])
        out.append(v)
    return tuple(out)

def main():
    print("(1) dir(z) = J diag(i^[z.x]) J exactly:", {m: check_dir_walsh(m) for m in (2, 3, 4, 5)})
    for k in (2, 3):
        for m in range(k, 11):
            n = 2 ** m
            vals = set()
            for z in range(1, n):
                zz = np.array([[((z >> a) & 1) * ((z >> b) & 1) for b in range(m)] for a in range(m)])
                vals.add(phi_of_sym(zz, k))
            MK = np.diag([1] * k + [0] * (m - k))
            tgt = phi_of_sym(MK, k)
            unit = all(sum(v) <= 1 for v in vals)
            if m <= 6 or m == 10:
                print(f"(2) k={k} m={m}: phi(zz^T) over all {n-1} z: {sorted(vals)}  all weight<=1: {unit};  phi(M_K) = {tgt}")
            assert unit and tgt == (1,) * k
    # (1') homomorphism on random words, exact: q of a word = sum of q of letters; B_q = sum of z z^T; phi additive
    rng = random.Random(1)
    for m in (3, 4, 5):
        bad = 0
        for trial in range(200):
            q = [0] * (2 ** m); Msum = np.zeros((m, m), dtype=np.int64); kind = []
            for _ in range(rng.randrange(1, 12)):
                t = rng.random()
                if t < 0.6:
                    z = rng.randrange(1, 2 ** m); q = q_add(q, q_dir(z, m))
                    Msum = (Msum + np.array([[((z >> a) & 1) * ((z >> b) & 1) for b in range(m)] for a in range(m)])) % 2
                elif t < 0.9:
                    q = q_add(q, q_shift(rng.randrange(2 ** m), m))
                else:
                    q = [(v + 1) % 4 for v in q]           # scalar i
            if not bil_full_ok(q, m) or not np.array_equal(bil(q, m), Msum): bad += 1
        print(f"(1') m={m}: 200 random words in dir/shift/i: exponent is a Z/4 quadratic form with B_q = sum z z^T: {bad == 0}")
    # targets
    for m in (3, 4, 5):
        for k in (2, 3):
            q = [0] * (2 ** m)
            for j in range(k): q = q_add(q, q_dir(1 << j, m))
            print(f"    m={m} k={k}: B_k (x) I has B_q = M_K: {np.array_equal(bil(q, m), np.diag([1]*k + [0]*(m-k)))}", end='')
        print()
    # (3) certificate numbers, computed in the representation (dimensions of subspaces of Q^(k 2^k))
    for k in (2, 3):
        N = 2 ** k
        def rep(v):          # v in F_2^k (tuple) -> permutation matrix on (+)_j Q[F_2^k]
            P = np.zeros((k * N, k * N), dtype=np.int64)
            s = sum(b << j for j, b in enumerate(v))
            for j in range(k):
                for x in range(N): P[j * N + (x ^ s), j * N + x] = 1
            return P
        Y = np.zeros((k * N, k * N // 2), dtype=np.int64); c = 0
        for j in range(k):
            for x in range(N):
                if not (x >> j) & 1: Y[j * N + x, c] = 1; c += 1
        def f(v):
            g = rep(v)
            return int(np.linalg.matrix_rank(np.hstack([Y, g @ Y]))) - int(np.linalg.matrix_rank(Y))
        unitf = {tuple(int(i == j) for i in range(k)): f(tuple(int(i == j) for i in range(k))) for j in range(k)}
        print(f"(3) k={k}: dim V = {k*N}, dim Y_0 = {Y.shape[1]}; f(0) = {f((0,)*k)}, f(unit vectors) = {unitf}, "
              f"f(target (1..1)) = {f((1,)*k)}  ->  moves >= {f((1,)*k)}/{max(unitf.values())} W = {f((1,)*k)/max(unitf.values())} W")
    # (4) address maps that keep phi: active rows a permutation of active coordinates, spectator rows arbitrary
    for k in (2, 3):
        for m in (k + 1, k + 2):
            n = 2 ** m; okA = 0; tot = 0; cnt = {}
            perms = list(itertools.permutations(range(k)))
            for _ in range(300):
                # random A: rows j<k = e_pi(j); rows >= k random, A invertible
                while True:
                    pi = rng.choice(perms)
                    rows = [1 << pi[j] for j in range(k)] + [rng.randrange(n) for _ in range(m - k)]
                    # invertible over F_2?
                    r = rows[:]; rank = 0
                    for bit in range(m):
                        piv = next((i for i in range(rank, m) if (r[i] >> bit) & 1), None)
                        if piv is None: continue
                        r[rank], r[piv] = r[piv], r[rank]
                        for i in range(m):
                            if i != rank and (r[i] >> bit) & 1: r[i] ^= r[rank]
                        rank += 1
                    if rank == m: break
                def Az(z): return sum((popc(rows[j] & z) & 1) << j for j in range(m))
                Qs = Q_list(k)
                # phi(A zz^T A^T) = (Q_j(Az)); compare as multisets up to the permutation pi of the k copies
                # need ONE permutation sigma of the k copies with Q_j(Az) = Q_sigma(j)(z) for all z (phi o A = sigma o phi)
                good = any(all(Qs[j](Az(z)) == Qs[sg[j]](z) for z in range(1, n) for j in range(k))
                           for sg in itertools.permutations(range(k)))
                tot += 1; okA += good
                cnt[good] = cnt.get(good, 0) + 1
            print(f"(4) k={k} m={m}: random address maps with active rows = a permutation pi of the active bits: "
                  f"phi o A = sigma o phi in {okA}/{tot}")

if __name__ == '__main__':
    main()

# ---------- (5) m-independent floor in the FULL game: permutation module on Pauli labels ------------------------------
# V = Q[F_2^(2m)] (labels (x|w) of X^x Z^w; every Clifford letter permutes labels through Sp(2m,2)); Y_0 = Q[L_0],
# L_0 = {(0|w)}.  Shifts, phases i^{z.x+c} and address maps fix L_0; dir(z) acts as the transvection along (z|0).
# f(g) = |g L_0 minus L_0|.  Bound = f(B_k (x) I) / f(dir e1) = (2^m - 2^(m-k)) / 2^(m-1) = 2 - 2^(1-k), every m.
def label_floor(m):
    n = 2 ** m
    def T(z, lab):          # transvection along (z|0): (x|w) -> (x + (z.w) z | w)
        x, w = lab
        return (x ^ (z if popc(z & w) & 1 else 0), w)
    def phase(z, lab):      # i^{z.x}: (x|w) -> (x | w + (z.x) z)   (symplectic image of the diagonal phase)
        x, w = lab
        return (x, w ^ (z if popc(z & x) & 1 else 0))
    def addr(A, lab):       # x -> A x on addresses: (x|w) -> (Ax | A^-T w); here A = elementary x_i ^= x_j
        (i, j), (x, w) = A, lab
        x2 = x ^ (((x >> j) & 1) << i)
        w2 = w ^ (((w >> i) & 1) << j)       # A^-T for the elementary transvection
        return (x2, w2)
    L0 = {(0, w) for w in range(n)}
    keep = all({phase(z, l) for l in L0} == L0 for z in range(n)) and \
           all({addr((i, j), l) for l in L0} == L0 for i in range(m) for j in range(m) if i != j)
    caps = {len({T(z, l) for l in L0} - L0) for z in range(1, n)}
    out = {}
    for k in (2, 3):
        img = L0
        for i in range(k): img = {T(1 << i, l) for l in img}
        out[k] = len(img - L0)
    return keep, caps, out
for m in (2, 3, 4, 5, 6):
    keep, caps, out = label_floor(m)
    print(f"(5) m={m}: free letters keep L_0: {keep}; f(dir z) over all z: {caps}; f(B_2 (x) I) = {out[2]} -> {out[2]/max(caps):.4f}; "
          f"f(B_3 (x) I) = {out[3]} -> {out[3]/max(caps):.4f}")
