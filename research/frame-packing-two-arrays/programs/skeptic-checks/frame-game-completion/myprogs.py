"""myprogs.py (v-pk2-complete): MY OWN frame-game programs, built so that two walkers reaching the SAME symplectic
element by DIFFERENT paths (different moves, different free letters) are added in one register BEFORE further moves.
In the engine, a naive lift gives (1 + n) dir.. with n a Pauli, which no free left factor repairs; Lambda must be
path independent for the lift to work.  Readouts are known by construction (no solver).
  conv1 (W=1): a = T_e1 y; b = p T_z (p0 y) with p T_z p0 = T_e1; d = a + b; d = T_e2..T_em d; read d/2.   M = m + 1
  conv2 (W=2): d = a + b + T_e1 y2 (same a, b); chain d; e = plain chain of y2; read y1: (d - e)/2, y2: e.  M = 2m + 1
  conv3 (W=1): as conv1 but b goes through TWO extra moves (b = p T_z2 p1 T_z1 p0 y, product = T_e1), and the
               register d = 3a - b is used, read d/2.                                                     M = m + 2
usage: myprogs.py m seed  -> progs_mine_m{m}.json"""
import sys, json, random
import numpy as np
from veng import *
m = int(sys.argv[1]); rng = random.Random(int(sys.argv[2]))
E = VEng(m); E.build_Ptable(); n = E.n
keys = list(E.Ptab.keys())
def Sp_of(k): return E.gbar(E.Ptab[k][0])
def rows(M):   # their row-bitmask format: row i, bit j = M[i, j]
    return [int(sum(int(M[i, j]) << j for j in range(2 * m))) for i in range(2 * m)]
def inv2(M):  # symplectic inverse over F2 via J M^T J
    J = np.block([[np.zeros((m, m), int), np.eye(m, dtype=int)], [np.eye(m, dtype=int), np.zeros((m, m), int)]])
    return J @ M.T @ J % 2
I2 = np.eye(2 * m, dtype=np.int64)
T = {z: E.gbar(E.dirm(z)) for z in range(1, n)}
for z in T: assert ((inv2(T[z]) @ T[z]) % 2 == I2).all()
def find1():
    while True:
        p0 = Sp_of(rng.choice(keys)); z = rng.randrange(1, n)
        if z == 1: continue
        p = T[1] @ inv2(T[z] @ p0 % 2) % 2
        if E.in_P(p): return p0, z, p
def find2():
    while True:
        p0 = Sp_of(rng.choice(keys)); p1 = Sp_of(rng.choice(keys)); z1 = rng.randrange(1, n); z2 = rng.randrange(1, n)
        F = T[z2] @ p1 @ T[z1] @ p0 % 2
        p = T[1] @ inv2(F) % 2
        if E.in_P(p) and not E.in_P(p1 @ T[z1] @ p0 % 2): return p0, z1, p1, z2, p
progs = []
def chain_rest(ops, r):
    for j in range(1, m): ops.append(['move', 1 << j, r]); r = len(ops) - 1
    return r
# conv1
p0, z, p = find1()
ops = [['in', 0]]
ops.append(['move', 1, 0]); a = 1
ops.append(['mix', [['1', rows(p0), 0]]]); ops.append(['move', z, 2]); ops.append(['mix', [['1', rows(p), 3]]]); b = 4
ops.append(['mix', [['1', rows(I2), a], ['1', rows(I2), b]]]); d = chain_rest(ops, 5)
progs.append(dict(name='conv1', W=1, ops=ops, read=[[[d, rows(I2), '1/2']]]))
# conv2
ops = [['in', 0], ['in', 1]]
ops.append(['move', 1, 0]); a = 2
ops.append(['mix', [['1', rows(p0), 0]]]); ops.append(['move', z, 3]); ops.append(['mix', [['1', rows(p), 4]]]); b = 5
ops.append(['move', 1, 1]); c = 6
ops.append(['mix', [['1', rows(I2), a], ['1', rows(I2), b], ['1', rows(I2), c]]]); d = chain_rest(ops, 7)
ops.append(['move', 1, 1]); e = chain_rest(ops, len(ops) - 1)
progs.append(dict(name='conv2', W=2, ops=ops, read=[[[d, rows(I2), '1/2'], [e, rows(I2), '-1/2']], [[e, rows(I2), '1']]]))
# conv3
p0, z1, p1, z2, p = find2()
ops = [['in', 0]]
ops.append(['move', 1, 0]); a = 1
ops.append(['mix', [['1', rows(p0), 0]]]); ops.append(['move', z1, 2]); ops.append(['mix', [['1', rows(p1), 3]]])
ops.append(['move', z2, 4]); ops.append(['mix', [['1', rows(p), 5]]]); b = 6
ops.append(['mix', [['3', rows(I2), a], ['-1', rows(I2), b]]]); d = chain_rest(ops, 7)
progs.append(dict(name='conv3', W=1, ops=ops, read=[[[d, rows(I2), '1/2']]]))
json.dump(progs, open(f'progs_mine_m{m}.json', 'w'))
print('ok', m, [(q['name'], sum(1 for o in q['ops'] if o[0] == 'move')) for q in progs])
