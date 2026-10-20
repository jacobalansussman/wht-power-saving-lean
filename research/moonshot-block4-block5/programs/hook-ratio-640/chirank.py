"""o3-ratio (A'): TRUE hook values f_k(A) = sum over all 16 translation characters chi of rank(M_chi), for
translation-equivariant A (products of dirs).  M_chi[S, T] = sum_c chi(c) A[c,0] lam_T[S + c] eps(c, S) (rows: out reps
S not containing 0; cols: N reps T containing 0); chi = 0 is blkk.block0.  For even k all blocks have equal rank (Z_c
maps chi to chi + c); for odd k they need not.  Characters are grouped into orbits of the address maps in GL(4,2) that
fix A's set of directions (pi D(z) pi^-1 = D(pi z)), so only one block per orbit is ranked.
usage: HOOKK=k PRIME=p python3 -I chirank.py name1,name2   (d1, d2, d3, d3b, d4)  [ALLCHI=1: rank every chi, no orbits]"""
import sys, os, time, subprocess, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import blkk as B
here = os.path.dirname(os.path.abspath(__file__))
p, k, n = B.p, B.k, 16
ELTS = {'d1': (1,), 'd2': (1, 2), 'd3': (1, 2, 4), 'd3b': (3, 5, 9), 'd4': (1, 2, 4, 8)}
def apply(Mcols, v):  # Mcols = images of e0..e3 (ints), v int
    r = 0
    for j in range(4):
        if v >> j & 1: r ^= Mcols[j]
    return r
GL = []
for cols in itertools.product(range(1, 16), repeat=4):
    if len({apply(cols, v) for v in range(16)}) == 16: GL.append(cols)
assert len(GL) == 20160
def chi_orbits(zs):
    st = [g for g in GL if sorted(apply(g, z) for z in zs) == sorted(zs)]
    # dual action on characters: chi -> chi o g^-1 ; chi(x) = popcount(chi & x) mod 2. Orbit of chi = {chi o g : g in st}
    seen, orbs = set(), []
    for chi in range(16):
        if chi in seen: continue
        orb = set()
        for g in st:
            # chi o g : x -> chi . (g x); as a vector: bit j = chi . g(e_j)
            orb.add(sum(((bin(chi & g[j]).count('1') & 1) << j) for j in range(4)))
        seen |= orb; orbs.append((chi, len(orb)))
    return orbs, len(st)
def block_chi(A, chi):
    rows, cols = B.Outrows, B.Ncols
    M = np.zeros((len(rows), len(cols)), dtype=np.int64)
    a0 = A[:, 0]; sgn = np.array([1 if bin(chi & c).count('1') % 2 == 0 else p - 1 for c in range(n)], dtype=np.int64)
    w = (a0 * sgn) % p
    for j, ci in enumerate(cols):
        lam = B.wedge_cols(A, B.Sk[ci])
        G = (w[:, None] * lam[B.tr] % p * B.eps) % p
        M[:, j] = G.sum(axis=0)[rows] % p
    return M
def rank(M, tag):
    fn = os.path.join(here, f"C_{tag}.bin")
    with open(fn, 'wb') as fh:
        np.array([M.shape[0], M.shape[1], p], dtype=np.uint32).tofile(fh); M.astype(np.uint32).tofile(fh)
    out = subprocess.run([os.path.join(here, 'rankp'), fn], capture_output=True, text=True).stdout
    os.remove(fn); return int(out.split()[1])
for name in sys.argv[1].split(','):
    zs = ELTS[name]; A = B.prod(zs); t0 = time.time()
    orbs, nst = chi_orbits(zs) if not os.environ.get('ALLCHI') else ([(c, 1) for c in range(16)], 0)
    if os.environ.get('CHIS'): orbs = [(c, s) for c, s in orbs if c in [int(x) for x in os.environ['CHIS'].split(',')]]
    tot = 0; parts = []
    for chi, sz in orbs:
        r = rank(block_chi(A, chi), f"k{k}_{name}_{chi}"); tot += sz * r; parts.append(f"chi={chi}:{r}x{sz}")
    print(f"k={k} p={p} {name} dirs {zs}: f = {tot}  [{' '.join(parts)}] (stabiliser {nst}; {time.time()-t0:.0f}s)", flush=True)
