# c15_cal105.py -- calibration of the per-irreducible filters (hook and 105a) on REAL symmetric programs:
# close the plain chain (h1 on e, h2 on h1 e, h3 on h2 h1 e; one fixed array) under H (exact symmetry
# tau_h: letter m -> c m c^-1, content g -> c g c', c' = K^-1 c^-1 K). The closure is an H-symmetric program computing K,
# so Theorem R holds in every irreducible psi of H:  sum_orbits <psi, Ind_T^H cap> >= <psi, demand>  for both certificates.
import numpy as np, time
from lib514 import *
import c7_symscreen2 as S7
import c14_105a as C14
t0 = time.time()
Kinv = inv(K)
def levi_of_perm(p):
    A = np.zeros((3, 3), dtype=int)
    for i in range(3): A[p[i], i] = 1
    return levi(A)
l2 = levi_of_perm([1, 0, 2]); l3 = levi_of_perm([2, 1, 0])
assert np.array_equal(comp(comp(l2, h1), inv(l2)), h2) and np.array_equal(comp(comp(l3, h1), inv(l3)), h3)
steps = [(h1, ID, ID), (h2, h1, l2), (h3, comp(h2, h1), l3)]     # (letter, content element, b with letter = b d b^-1)
def check(H, name):
    n = len(H); cls, rep, sizes, chars, clsof, idx = S7.char_table(H)
    ip = lambda a, b: np.real(np.sum(sizes * a * np.conj(b)) / n)
    tot_h = np.zeros(len(chars)); tot_5 = np.zeros(len(chars)); nmoves = 0
    for (m, g, b) in steps:
        T = []
        for c in H:
            cp = comp(comp(Kinv, inv(c)), K)
            if np.array_equal(comp(comp(c, m), inv(c)), m) and np.array_equal(comp(comp(c, g), cp), g): T.append(c)
        bi = inv(b); xs = [comp(comp(bi, t), b) for t in T]
        assert all(key(x) in S7.Ekeys for x in xs)
        ch_h = [S7.cap_trace(x) for x in xs]; ch_5 = [C14.capN_trace(x) for x in xs]
        tot_h += np.array([ip(S7.induced(H, rep, T, ch_h), c) for c in chars])
        tot_5 += np.array([ip(S7.induced(H, rep, T, ch_5), c) for c in chars])
        nmoves += n // len(T)
    dem_h = np.array([ip(S7.induced(H, rep, H, [S7.yperp_trace(x) for x in H]), c) for c in chars])
    dem_5 = np.array([ip(S7.induced(H, rep, H, [C14.demN_trace(x) for x in H]), c) for c in chars])
    ok_h = np.all(tot_h - dem_h > -1e-6); ok_5 = np.all(tot_5 - dem_5 > -1e-6)
    print(f'{name}: symmetrised plain program has {nmoves} moves for 1 fixed array; hook per-irrep caps {np.round(tot_h,3)} vs demand {np.round(dem_h,3)} -> {ok_h};'
          f' 105a caps {np.round(tot_5,3)} vs demand {np.round(dem_5,3)} -> {ok_5}', flush=True)
    return ok_h and ok_5
Lsubs = S7.subgroups_upto_conj(LEVI)
S3c = gen_closure([levi_of_perm([1, 0, 2]), levi_of_perm([1, 2, 0])])
allok = check(S3c, 'S3 (coordinate permutations)')
# the address-111 S4 and others: use the conjugate that contains the coordinate S3 when possible
S4_111 = gen_closure([g for g in LEVI if int(g[7]) == 7])           # address maps fixing x = 111 (x bits 0..2 -> v=7)
print('|S4_111| =', len(S4_111))
allok &= check(S4_111, 'S4 fixing address 111')
for H in Lsubs:
    if len(H) in (1, 168): continue
    allok &= check(H, f'class rep of order {len(H)}')
print('ALL CALIBRATIONS PASS:', allok, '[%.0fs]' % (time.time() - t0))
