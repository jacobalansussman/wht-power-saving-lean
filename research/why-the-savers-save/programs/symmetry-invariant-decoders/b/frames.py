"""frames.py (d-R7): pools of subspaces of F_2^9 (dims 6, 7, 8) as port masks on P112, octad Hall test for a
family of totals, and a local search for families (multiset of dims) that pass Hall on all 1260 octads."""
import sys, itertools, time; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
import warnings; warnings.filterwarnings('ignore')
from model import *
import numpy as np
LABA = np.array(LAB, dtype=np.int64)
PT = np.array([[par(h & int(u)) for u in LABA] for h in range(512)], dtype=np.uint8)   # 512 x 112
OCTA = np.array(OCTP, dtype=np.int64)          # 1260 x 8

def duals(k):
    """reduced echelon bases (k vectors) of all k-dim subspaces of F_2^9"""
    out = []
    for piv in itertools.combinations(range(8, -1, -1), k):
        frees = [[j for j in range(p) if j not in piv] for p in piv]
        for bits in itertools.product(*[range(1 << len(f)) for f in frees]):
            out.append(tuple((1 << p) | sum(1 << j for i, j in enumerate(f) if b >> i & 1) for p, f, b in zip(piv, frees, bits)))
    return out

def pool(d, minports=0):
    """all d-dim subspaces (as duals) with their P-port masks; returns (duals list, bool array n x 112)"""
    ds = duals(9 - d)
    D = np.array(ds, dtype=np.int64)
    M = np.ones((len(ds), NP), dtype=bool)
    for i in range(9 - d): M &= (PT[D[:, i]] == 0)
    keep = M.sum(1) >= minports
    return [ds[i] for i in np.nonzero(keep)[0]], M[keep]

def hall_fail(S):
    """S: K x 112 bool support.  Returns bool array (1260,) True where the octad's K x 8 pattern has no matching of
    size 8 (columns = ports)."""
    K = S.shape[0]
    A = S[:, OCTA]                      # K x 1260 x 8
    st = np.zeros((OCTA.shape[0], 1 << K), dtype=bool); st[:, 0] = True
    bits = 1 << np.arange(K)
    idx = np.arange(1 << K)
    for p in range(8):
        new = np.zeros_like(st)
        for k in range(K):
            src = idx[(idx & bits[k]) == 0]
            can = st[:, src] & A[k, :, p][:, None]
            new[:, src | bits[k]] |= can
        st = new
    return ~st.any(1)

POPC = np.array([bin(i).count('1') for i in range(1 << 12)], dtype=np.int64)
SUBSZ = np.array([bin(q).count('1') for q in range(256)], dtype=np.int64)
def hall_deficit(S):
    """fast Hall: S K x 112 bool (K <= 12).  Returns array (1260,) of max_Q (|Q| - |N(Q)|) (> 0 means Hall fails)."""
    K = S.shape[0]
    w = (S.astype(np.int64) << np.arange(K)[:, None]).sum(0)     # 112 port masks
    pm = w[OCTA]                                                 # 1260 x 8
    N = np.zeros((256, OCTA.shape[0]), dtype=np.int64)
    for q in range(1, 256):
        low = (q & -q).bit_length() - 1
        N[q] = N[q & (q - 1)] | pm[:, low]
    deficit = SUBSZ[:, None] - POPC[N]
    return deficit.max(0)
def hall_fail(S):
    return hall_deficit(S) > 0
