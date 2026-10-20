# t10check.py (v-w8-rmexact): end-to-end check of Theorem MT + Theorem 10' on a genuine unit that USES address maps.
# Physical unit (1 in 4, H = Singer Z_15): r = [1]; for j = 1..4: r <- P_{g_j} [e1 e1^T] P_{g_j}^-1 r with g_j e1 = e_j.
# Push it through Phi (rows k of the n x n matrix; [S] acts on row k by [k^-1 S k^-T]; P_g sends row k to row g k),
# record the Q-independent transfer data (B, Gamma, Atil) symbolically, then for EVERY form Q check
#   (i)  Atil_D (I + 2 Gamma_D)^-1 B_D = O_Q   (the program serves virtual array k with target k^-1 k^-T), and
#   (ii) rank E[D, D] = |D| - n_odd(Q),  E = I + 2 Gamma - B Atil   (Theorem 10'), ranks mod two primes.
import sys
sys.argv = ['x']
exec(open('vchain.py').read().split("if __name__ == '__main__':")[0])
def gmul(a, b):
    r = 0
    for i in range(4):
        if (b >> i) & 1: r ^= a << i
    for d in (6, 5, 4):
        if (r >> d) & 1: r ^= 0b10011 << (d - 4)
    return r
cols = [gmul(2, 1 << j) for j in range(4)]
A = tuple(tuple((cols[j] >> i) & 1 for j in range(4)) for i in range(4))
H = cyclic(A); n = len(H); idx = {h: i for i, h in enumerate(H)}
inv = {h: next(g for g in H if mul(h, g) == ID) for h in H}
gj = [next(g for g in H if act(g, 1) == 1 << j) for j in range(4)]
# registers: one physical register; virtual rows k (indexed by H element); each row: dict symbol -> coeff
rows = {h: {('in', idx[h]): 1} for h in H}
moves = []   # (direction u, row-vector snapshot)
def perm(g):  # left multiplication by P_g: new row (g k) = old row k
    global rows
    rows = {mul(g, k): rows[k] for k in H}
for j in range(4):
    perm(inv[gj[j]])
    z = 1   # direction e1
    for k in H:   # [z z^T] on row k acts as move in direction k^-1 z
        u = act(inv[k], z)
        t = len(moves); moves.append((u, dict(rows[k])))
        rows[k] = dict(rows[k]); rows[k][('w', t)] = rows[k].get(('w', t), 0) + 1
    perm(gj[j])
P = len(moves); N = n
B = [[mv[1].get(('in', i), 0) for i in range(N)] for mv in moves]
Gam = [[moves[t][1].get(('w', s), 0) for s in range(P)] for t in range(P)]
Ain = [[rows[k].get(('in', i), 0) for i in range(N)] for k in H]
At = [[rows[k].get(('w', t), 0) for t in range(P)] for k in H]
assert Ain == [[int(i == idx[k]) for i in range(N)] for k in H], 'input part not identity'
assert all(Gam[t][s] == 0 for t in range(P) for s in range(t, P)), 'Gamma not strictly lower'
targets = [mul(k, tp(k)) for k in H]   # virtual array k (row k) target k^-1 k^-T
def rank_mod(M, p):
    M = [[x % p for x in r] for r in M]; r = 0; nr = len(M); nc = len(M[0]) if M else 0
    for c in range(nc):
        piv = next((i for i in range(r, nr) if M[i][c]), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]; iv = pow(M[r][c], p - 2, p)
        M[r] = [x * iv % p for x in M[r]]
        for i in range(nr):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        r += 1
    return r
def solve_mod(L, Bm, p):   # L lower unitriangular (mod p): solve L X = Bm
    X = []
    for i in range(len(L)):
        row = [Bm[i][c] % p for c in range(len(Bm[0]))]
        for s in range(i):
            if L[i][s]:
                row = [(a - L[i][s] * b) % p for a, b in zip(row, X[s])]
        X.append(row)
    return X
bad1 = bad2 = 0; p1, p2 = 2147483647, 1000000007
for Q in range(1, NF):
    D = [t for t in range(P) if qval(Q, moves[t][0])]
    odd = [pairing(Q, targets[k]) for k in range(N)]
    nodd = sum(odd)
    if not D:
        if nodd: bad1 += 1
        continue
    L = [[int(a == b) + 2 * Gam[a][b] for b in D] for a in D]
    X = solve_mod(L, [B[t] for t in D], p1)     # L^-1 B_D
    AtD = [[At[k][t] for t in D] for k in range(N)]
    prod = [[sum(AtD[k][i] * X[i][c] for i in range(len(D))) % p1 for c in range(N)] for k in range(N)]
    if prod != [[odd[k] if k == c else 0 for c in range(N)] for k in range(N)]: bad1 += 1
    E = [[int(a == b) + 2 * Gam[a][b] - sum(B[a][k] * At[k][b] for k in range(N)) for b in D] for a in D]
    r1, r2 = rank_mod(E, p1), rank_mod(E, p2)
    if not (r1 == r2 == len(D) - nodd): bad2 += 1
print('virtual moves', P, 'virtual arrays', N, 'distinct targets', len(set(targets)))
print('forms where the transfer identity fails:', bad1, ' forms where rank E[D,D] != |D| - n_odd:', bad2, ' (of', NF - 1, ')')
