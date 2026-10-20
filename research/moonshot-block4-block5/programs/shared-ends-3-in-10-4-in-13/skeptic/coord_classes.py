"""v-o3-shared: in the coordinate world, how many end-element classes E ~ r E b (r in R, b any Z4-monomial) are there
among the one-move coordinate end elements D_e0..D_e3, for R = signs (16 translations), o3 = <s012> x translations,
o4 = <c0123> x translations, addr?  (o3-shared claims 'every end element is one class' in the coordinate world.)"""
import sys, itertools
import numpy as np, warnings
warnings.filterwarnings('ignore')
n = 16; m = 4; A = (1 + 1j) / 2; B = (1 - 1j) / 2; I = np.eye(n, dtype=complex)
def X(z):
    M = np.zeros((n, n), complex)
    for x in range(n): M[x ^ z, x] = 1
    return M
def D(z): return A * I + B * X(z)
def perm_mat(perm):
    sig = [sum(((x >> j) & 1) << perm[j] for j in range(m)) for x in range(n)]
    M = np.zeros((n, n), complex)
    for y in range(n): M[sig[y], y] = 1
    return M
def is_mono_z4(M):
    nz = np.abs(M) > 1e-9
    if not ((nz.sum(0) == 1).all() and (nz.sum(1) == 1).all()): return False
    return all(min(abs(v - z) for z in (1, 1j, -1, -1j)) < 1e-9 for v in M[nz])
perms = {'sign': [(0, 1, 2, 3)], 'o3': [(0, 1, 2, 3), (1, 2, 0, 3), (2, 0, 1, 3)],
         'o4': [(0, 1, 2, 3), (1, 2, 3, 0), (2, 3, 0, 1), (3, 0, 1, 2)], 'addr': list(itertools.permutations(range(4)))}
Es = [D(1 << k) for k in range(4)]
for name, ps in perms.items():
    R = [perm_mat(p) @ X(v) for p in ps for v in range(n)]
    # E_i ~ E_j iff exists r: (r E_i)^-1 E_j monomial
    cls = []
    for j, Ej in enumerate(Es):
        for c in cls:
            Ei = Es[c[0]]
            if any(is_mono_z4(np.linalg.inv(r @ Ei) @ Ej) for r in R):
                c.append(j); break
        else:
            cls.append([j])
    print(name, '|R| =', len(R), 'classes of D_e0..D_e3:', cls)
# the selective-advance example: g1 = D_e1 D_e2 D_e3, g2 = D_e0 D_e2 D_e3, E = D_e0
K = D(1) @ D(2) @ D(4) @ D(8)
g1 = D(2) @ D(4) @ D(8); g2 = D(1) @ D(4) @ D(8)
E = K @ np.linalg.inv(g1)
dep = lambda M: int(round(np.log2((np.abs(M) > 1e-9).sum(1).max())))
print('E = K g1^-1 equals D_e0:', np.allclose(E, D(1)), '| depth(E g1) =', dep(E @ g1), 'E g1 == K:', np.allclose(E @ g1, K),
      '| depth(E g2) =', dep(E @ g2), 'E g2 == X_e0 D_e2 D_e3:', np.allclose(E @ g2, X(1) @ D(4) @ D(8)))
