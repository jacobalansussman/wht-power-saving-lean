"""F_2 ranks of invariant candidates for the BW16 orthogonality pattern: N = I + (some orthogonal pairs)."""
import numpy as np, sys, time
L = np.load('bw16_L.npy'); v = len(L)
G = (L @ L.T) // 2
def f2rank(M):
    M = (M & 1).astype(np.uint8)
    P = np.packbits(M, axis=1)
    rows = [int.from_bytes(r.tobytes(), 'big') for r in P]
    basis = {}
    for x in rows:
        while x:
            h = x.bit_length() - 1
            if h in basis: x ^= basis[h]
            else: basis[h] = x; break
    return len(basis)
I = np.eye(v, dtype=np.int64)
O = (G == 0).astype(np.int64)
A1 = (np.abs(G) == 1).astype(np.int64); A2 = (np.abs(G) == 2).astype(np.int64)
t = time.time()
print('rank F2 of I+A_orth:', f2rank(I + O), '%.1fs' % (time.time()-t))
print('rank F2 of G mod 2 (=A1):', f2rank(A1), ' rank of (G/2 on even) A2:', f2rank(A2))
print('rank F2 of J + G mod 2:', f2rank(1 + G))
# refine orthogonal pairs by counting common |ip|=2 neighbours and common |ip|=1 neighbours
C2 = A2 @ A2; C1 = A1 @ A1; C12 = A1 @ A2
keys = {}
iu = np.where(O == 1)
inv = np.stack([C2[iu], C1[iu], C12[iu], C12.T[iu]], axis=1)
u, cnt = np.unique(inv, axis=0, return_counts=True)
print('orthogonal-pair classes by (common ip2, common ip1, ip1-ip2, ip2-ip1):')
for a, c in zip(u.tolist(), cnt.tolist()): print('  ', a, c)
np.save('inv_orth.npy', inv)
# same for |ip|=1 and |ip|=2 pairs
for name, A in (('ip1', A1), ('ip2', A2)):
    iu2 = np.where(A == 1)
    inv2 = np.stack([C2[iu2], C1[iu2]], axis=1)
    u2, c2 = np.unique(inv2, axis=0, return_counts=True)
    print(name, 'classes:', list(zip(u2.tolist(), c2.tolist())))
