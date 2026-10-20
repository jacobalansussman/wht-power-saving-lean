"""quad.py -- w5-twoshared (resumed): is the "quad lemma" geometric?  For private cores, count over ALL sigma in
Stab(L0) (20,643,840 symplectic engine letters e') the solutions of WEAKER conditions than Lemma D, and how many of
them survive round 2 (both tops' readings of w_J at layer <= 2).  Per top e (j = 1 - e), with I(h) = S_h(L0):
  K   : lambda(w_k) <= 2 for k = 0, 1         <=> dim(L_D cap sigma I(q_k)) >= 2           (B = sigma^-1 L_D)
  SJ  : w_J of injection shape (I(w_J) ~ L_D)  <=> dim(sigma L_D cap I(J)) >= 3
  SK  : w_k of injection shape (I(w_k) ~ L_D)  <=> dim(sigma I(q_k) cap S_(y_e) L_D) >= 3
  R2  : lambda(N1) <= 2 and lambda(y_j w_J) <= 2 <=> dim(sigma A_e cap I(J)) >= 2 and dim(sigma A_j cap I(J)) >= 2
        (A_x = Lag(y_x)).
Prints #K, #K&R2, #K&SJ&SK, #K&SJ&SK&R2 per top (and the exact Lemma-D count for reference is in rd2batch).
Usage: quad.py ncores seed [seed ...]"""
import sys, os, random, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from esearch import basis, img_L0, GL, SYM, PX, PTZ, BT
import sp
from sweep import core
D = sp.D; comp = sp.comp; inv = sp.inv; lag = sp.lag
LD = img_L0(D)
# rank of four vectors of F2^4, table over 16^4
R4 = np.zeros(1 << 16, dtype=np.int8)
for idx in range(1 << 16):
    vs = [(idx >> (4 * t)) & 15 for t in range(4)]
    sset = {0}
    for v in vs: sset |= {w ^ v for w in sset}
    R4[idx] = len(sset).bit_length() - 1

def quot_table(T):
    """4-bit coordinates of F2^8 / T (T a 4-dim subspace)"""
    Tb = basis(T); comp_ = []
    span_ = set(T)
    for v in range(1, 256):
        if v not in span_:
            comp_.append(v); span_ = span_ | {w ^ v for w in span_}
        if len(comp_) == 4: break
    tab = np.zeros(256, dtype=np.int64)
    for c in range(16):
        base = 0
        for t in range(4):
            if (c >> t) & 1: base ^= comp_[t]
        for w in T: tab[base ^ w] = c
    return tab

def images(A):
    """sigma(a) for the basis vectors a of A, over all of Stab(L0): list of (nP, 1024) int arrays"""
    out = []
    for a in basis(A):
        ax, az = a & 15, a >> 4
        px = PX[:, ax]; pz = PTZ[:, az]
        z2 = pz[:, None] ^ BT[:, px].T
        out.append(px[:, None] | (z2 << 4))
    return out

def dim_cap(A, T):
    """dim(sigma A cap T) over all sigma"""
    tab = quot_table(T)
    ims = images(A)
    idx = tab[ims[0]] | (tab[ims[1]] << 4) | (tab[ims[2]] << 8) | (tab[ims[3]] << 12)
    return 4 - R4[idx]

ncores = int(sys.argv[1]); tot = {}
for seed in [int(a) for a in sys.argv[2:]]:
    rng = random.Random(seed)
    for c in range(ncores):
        t0 = time.time()
        cr = core(rng, 'priv')
        if cr is None: continue
        R0, R1, _, _, tg = cr; R = [R0, R1]
        Iq = [img_L0(R[k]['q']) for k in (0, 1)]
        line = []
        for e in (0, 1):
            j = 1 - e
            y = comp(D, R[e]['mu'][0]); yj = comp(D, R[j]['mu'][0])
            J = comp(y, R[j]['q']); IJ = img_L0(J)
            Ae, Aj = lag(y), lag(yj)
            K = (dim_cap(Iq[0], LD) >= 2) & (dim_cap(Iq[1], LD) >= 2)
            R2 = (dim_cap(Ae, IJ) >= 2) & (dim_cap(Aj, IJ) >= 2)
            SJ = dim_cap(LD, IJ) >= 3
            YL = img_L0(comp(y, D))           # S_(y_e) L_D = S_(y_e) S_D L0
            SK = (dim_cap(Iq[0], YL) >= 3) & (dim_cap(Iq[1], YL) >= 3)
            r = (int(K.sum()), int((K & R2).sum()), int((K & SJ & SK).sum()), int((K & SJ & SK & R2).sum()))
            line.append(f"top {e}: K {r[0]}, K&R2 {r[1]}, K&shapes {r[2]}, K&shapes&R2 {r[3]}")
            for nm, val in zip(('K', 'K&R2', 'K&shapes', 'K&shapes&R2'), r):
                tot[nm] = tot.get(nm, 0) + (val > 0)
        print(f"seed {seed} core {c}: " + "; ".join(line) + f" ({time.time()-t0:.0f}s)", flush=True)
print("tops with a nonzero count:", tot, flush=True)
