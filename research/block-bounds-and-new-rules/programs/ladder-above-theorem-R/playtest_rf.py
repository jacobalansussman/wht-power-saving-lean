# pf-ladder: bookkeeping play test of Theorem RF / Theorem D in the coset module Q[G/F], engine game, m = 2.
# Builds G (46080 elements) and its 15 cosets gF by 0-1 BFS (F = elements with 0 directions), then plays random words
# (random mixing of all registers, random free letters, random directions on random registers) with W = 2 live
# registers and 2 scratch registers, and checks after every paid move:
#   Phi_0 rises by <= 1 (f = 1),  Phi_1(t+1) <= Phi_0(t) (shift step),  Phi_j non-decreasing, Phi_1 <= Phi_0;
# and for the plain schedule (2 directions per array) that Phi_1 at the end equals W (dem_(m-1) = 1 per array).
# Contents are tracked by their action on the coset [F]: for each register, a W x 15 rational matrix.
import itertools, random
from fractions import Fraction as Fr
from collections import deque
import sys
sys.path.insert(0, '.')
from ladder_check import gmul, dir2, perm_mat, diag_mat, I4

m, n, W = 2, 4, 2
def red(M, k):
    while k > 0 and all((e.real % 2 == 0 and e.imag % 2 == 0) for r in M for e in r):
        M = tuple(tuple(e / 2 for e in r) for r in M); k -= 1
    return (M, k)
free = [perm_mat(n, [x ^ z for x in range(n)]) for z in (1, 2)]
free.append(perm_mat(n, [((x & 1) << 1) | (x >> 1) for x in range(n)]))
free.append(perm_mat(n, [x ^ ((x & 1) << 1) for x in range(n)]))
for z in range(1, n): free.append(diag_mat([I4[bin(z & x).count('1') % 4] for x in range(n)]))
free.append(diag_mat([1j] * n))
dirs = [dir2(m, z) for z in range(1, n)]
start = red(perm_mat(n, list(range(n))), 0)
dist = {start: 0}; dq = deque([start])
while dq:
    key = dq.popleft(); M, k = key; dcur = dist[key]
    for f in free:
        nk = red(gmul(f, M), k)
        if dist.get(nk, 10**9) > dcur: dist[nk] = dcur; dq.appendleft(nk)
    for D in dirs:
        nk = red(gmul(D, M), k + 1)
        if dist.get(nk, 10**9) > dcur + 1: dist[nk] = dcur + 1; dq.append(nk)
Fset = [g for g, d in dist.items() if d == 0]
# cosets gF: canonical id = frozenset of the coset's keys (computed via right multiplication by F)
def rmul(key, fkey):
    (A, a), (B, b) = key, fkey
    return red(gmul(A, B), a + b)
coset_of = {}; reps = []
for g in dist:
    if g in coset_of: continue
    cid = len(reps); reps.append(g)
    for f in Fset: coset_of[rmul(g, f)] = cid
NC = len(reps)
level = [dist[reps[c]] for c in range(NC)]          # least number of directions in the coset = frame distance
def act(letter_key, cid):                           # left multiplication of a coset by a letter
    return coset_of[rmul(letter_key, reps[cid])]
fkeys = [red(f, 0) for f in free]; dkeys = [red(D, 1) for D in dirs]
Kkey = red(gmul(dir2(m, 1), dir2(m, 2)), 2)
print(f"|G| = {len(dist)}, |F| = {len(Fset)}, cosets = {NC}, levels {sorted(set(level))}, K's coset level = {level[coset_of[Kkey]]}")

def rank(rows):                                     # exact rank over Q of a list of vectors (lists of Fractions)
    rows = [list(r) for r in rows if any(r)]
    rk = 0; ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] != 0), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][c] != 0:
                fct = rows[i][c] / rows[rk][c]
                rows[i] = [x - fct * y for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk

def Phi(history, j):
    # S_j = { v in Q^W : sum_l v_l * a_l (a vector over cosets) has no weight on cosets of level > j, for all contents a }
    # Phi_j = rank of the linear conditions: for each content a and each coset c with level > j, the row (a[l][c])_l.
    rows = []
    for a in history:
        for c in range(NC):
            if level[c] > j: rows.append([a[l][c] for l in range(W)])
    return rank(rows) if rows else 0

def play(rng, plain=False):
    R = 4
    regs = [[[Fr(0)] * NC for _ in range(W)] for _ in range(R)]
    Fid = coset_of[start]
    for l in range(W): regs[l][l][Fid] = Fr(1)
    history = [r for r in regs[:W]]
    phis = [(Phi(history, 0), Phi(history, 1))]
    ok = True
    script = [(l, z) for l in range(W) for z in (0, 1)] if plain else None
    steps = len(script) if plain else rng.randrange(1, 6)
    for s in range(steps):
        if not plain:
            for _ in range(rng.randrange(0, 3)):           # random mixing (rational matrix)
                g = [[Fr(rng.randrange(-2, 3), rng.randrange(1, 3)) for _ in range(R)] for _ in range(R)]
                regs = [[[sum(g[i][k] * regs[k][l][c] for k in range(R)) for c in range(NC)] for l in range(W)] for i in range(R)]
            for _ in range(rng.randrange(0, 3)):           # random free letters on random registers
                i = rng.randrange(R); fk = rng.choice(fkeys)
                new = [[Fr(0)] * NC for _ in range(W)]
                for l in range(W):
                    for c in range(NC):
                        if regs[i][l][c]: new[l][act(fk, c)] += regs[i][l][c]
                regs[i] = new
            i = rng.randrange(R); dk = rng.choice(dkeys)
        else:
            i, z = script[s]; dk = dkeys[[0, 1][z]]
        new = [[Fr(0)] * NC for _ in range(W)]
        for l in range(W):
            for c in range(NC):
                if regs[i][l][c]: new[l][act(dk, c)] += regs[i][l][c]
        regs[i] = new
        history = history + [r for r in regs]              # every content so far generates M_t (with free images)
        # free images of contents keep level sets, so generators suffice; include all current registers
        p0, p1 = Phi(history, 0), Phi(history, 1)
        q0, q1 = phis[-1]
        if not (p0 <= q0 + 1 and p1 <= q0 and p0 >= q0 and p1 >= q1 and p1 <= p0): ok = False
        phis.append((p0, p1))
    return ok, phis, regs

rng = random.Random(7)
bad = 0
for trial in range(60):
    ok, phis, _ = play(rng)
    bad += (not ok)
print(f"60 random words: violations of (Phi_0 step <= 1, Phi_1(t+1) <= Phi_0(t), monotone, Phi_1 <= Phi_0): {bad}")
ok, phis, regs = play(rng, plain=True)
print(f"plain schedule (2 directions per array): ok = {ok}; (Phi_0, Phi_1) after each move: {phis}")
print(f"  end: Phi_1 = {phis[-1][1]} (Theorem D needs W = {W}); live registers hold K e_l: "
      f"{all(regs[l][l][coset_of[Kkey]] != 0 for l in range(W))}")
