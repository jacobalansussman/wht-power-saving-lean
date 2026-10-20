# an-openai: exact check of the frame identities the OpenAI unit rests on, at small m (pure Python; all entries
# are Gaussian integers over 2^m, so Python complex floats are exact here and we compare with tolerance 1e-12):
#   C_z^2 = R_z ; C_z^{-1} = R_z C_z ; Phi_F = C^{⊗m} with Phi_U = J diag(i^{q_U}) J, q_U(x) = wt(P_U x) mod 4 ;
#   Phi_{<u,w>} = Phi_<u> Phi_<w> for orthogonal unit vectors u, w ; Phi_{<u,w>} Phi_<u>^{-1} in {C_w, C_w^{-1}} ;
#   exceptional route ratio Phi_{U^perp} Phi_U^{-1} = R_u C^{⊗m} for a unit-norm u (U = <u>).
import random
def popc(x): return bin(x).count("1")
def dot(x, y): return popc(x & y) & 1
def mm(A, B):
    n = len(A); Bt = list(zip(*B))
    return [[sum(A[i][k] * Bt[j][k] for k in range(n)) for j in range(n)] for i in range(n)]
def close(A, B): return all(abs(A[i][j] - B[i][j]) < 1e-12 for i in range(len(A)) for j in range(len(A)))
def ct(A): return [[A[j][i].conjugate() for j in range(len(A))] for i in range(len(A))]
def run(m, trials=25, seed=3):
    rnd = random.Random(seed); n = 1 << m
    J = [[(-1) ** dot(x, y) / 2 ** (m / 2) for y in range(n)] for x in range(n)]
    I = [[1.0 + 0j if i == j else 0j for j in range(n)] for i in range(n)]
    def R(z): return [[1.0 + 0j if y == x ^ z else 0j for y in range(n)] for x in range(n)]
    a, b = (1 + 1j) / 2, (1 - 1j) / 2
    def Cz(z):
        Rz = R(z); return [[a * I[i][j] + b * Rz[i][j] for j in range(n)] for i in range(n)]
    Cm = [[1.0 + 0j]]
    for _ in range(m):  # kron with C
        C1 = [[a, b], [b, a]]
        Cm = [[Cm[i // 2][j // 2] * C1[i % 2][j % 2] for j in range(2 * len(Cm))] for i in range(2 * len(Cm))]
    # NB: kron ordering puts the first factor on the high bit; C^{⊗m} is symmetric under bit order anyway.
    def proj(Ub):
        k = len(Ub)
        # P_U x = sum_i lam_i b_i where G lam = (x.b_j)_j ; solve over F_2 by brute force on k <= m
        G = [[dot(p, q) for q in Ub] for p in Ub]
        table = {}
        for mask in range(1 << k):
            vec = tuple(sum(G[j][i] for i in range(k) if (mask >> i) & 1) % 2 for j in range(k))
            table[vec] = mask
        assert len(table) == 1 << k, "degenerate"
        def P(x):
            mask = table[tuple(dot(x, q) for q in Ub)]
            out = 0
            for i in range(k):
                if (mask >> i) & 1: out ^= Ub[i]
            return out
        return P
    def Phi(Ub, sign=1):
        if not Ub: return I
        P = proj(Ub)
        D = [[(1j ** ((sign * popc(P(x))) % 4)) if x == y else 0j for y in range(n)] for x in range(n)]
        return mm(mm(J, D), J)
    ok = True
    for z in range(1, n):
        C = Cz(z); ok &= close(mm(C, C), R(z)); ok &= close(ct(C), mm(R(z), C))
    full = [1 << i for i in range(m)]
    ok &= close(Phi(full), Cm)
    cnt = 0
    for _ in range(trials):
        u = rnd.choice([x for x in range(1, n) if popc(x) % 2 == 1])
        piv = {}; basis = []
        for y in range(n):
            if dot(y, u): continue
            yy = y
            while yy and (yy.bit_length() - 1) in piv: yy ^= piv[yy.bit_length() - 1]
            if yy: piv[yy.bit_length() - 1] = yy; basis.append(y)
        ok &= close(mm(Phi(basis), Phi([u], -1)), mm(R(u), Cm)); cnt += 1
        ws = [w for w in range(1, n) if popc(w) % 2 == 1 and dot(w, u) == 0]
        if ws:
            w = rnd.choice(ws)
            ok &= close(Phi([u, w]), mm(Phi([u]), Phi([w])))
            E = mm(Phi([u, w]), Phi([u], -1))
            ok &= (close(E, Cz(w)) or close(E, ct(Cz(w))))
    print("m=%d: all frame identities hold: %s (%d random unit vectors u)" % (m, bool(ok), cnt))
for m in (2, 3, 4, 5):
    run(m)
