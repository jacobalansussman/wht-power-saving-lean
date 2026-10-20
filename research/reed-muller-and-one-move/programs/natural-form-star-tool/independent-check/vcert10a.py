import numpy as np, vf
ks = [0,3,6,8,9,18,20,22,28,30]; r = vf.PAR[:, ks].sum(axis=1).astype(int); s = vf.WT - r
F = 17170; A = int(vf.MASK[F]); print("F", F, "|D_F|", vf.WT[F], "r_F", r[F], "s_F", s[F])
cells = np.zeros((32, 32), int); tot = 0
for l in open('../cert_10a_17170.txt'):
    t = l.split()
    if not t or t[0] != "P": continue
    R, C, b, _, G = map(int, t[1:])
    assert (R | C) & ~A == 0 and (R | C) & ~int(vf.MASK[G]) == 0
    bound = min(bin(R).count('1'), bin(C).count('1'), s[G])
    assert bound == b, (bound, b)
    tot += bound
    for x in range(1, 32):
        for y in range(1, 32):
            if R >> (x-1) & 1 and C >> (y-1) & 1: cells[x, y] += 1
D = [z for z in range(1, 32) if A >> (z-1) & 1]
part = all(cells[x, y] == 1 for x in D for y in D) and cells.sum() == len(D) ** 2
print("pieces partition D_F x D_F:", part, "sum of bounds", tot, "< s_F", s[F], ":", tot < s[F])
