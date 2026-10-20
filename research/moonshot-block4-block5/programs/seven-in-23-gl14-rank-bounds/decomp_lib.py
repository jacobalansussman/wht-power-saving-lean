"""decomp_lib.py -- w8-bpgl14: the Lift, the slot module D4, the joint sigma and the two twists p (F1's, from q) and p2
(O O F F O's second twist, ooffo.py's q2 rule), exactly as in decomp.py, without the commutant."""
import sys, os, random
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d4gen import *
_v16 = np.ones(16) / 4; _U16, _, _ = np.linalg.svd(np.eye(16) - np.outer(_v16, _v16)); _Q15 = _U16[:, :15].T
def sig_15p(b): return _Q15 @ sig_aff16(b) @ _Q15.T
def get_sig(name):
    if name == 'trivial': return sig_triv, 1
    if name == '15p': return sig_15p, 15
    if name == 'V': return sig_nat, 16
    return make_sig(name)
def load_twists(seed, psi, names):
    L = Lift(seed=seed, psi=psi)
    sigl = [get_sig(nm) for nm in names]
    dss = [ds for (_, ds) in sigl]; offs = np.cumsum([0] + dss); DS = int(offs[-1])
    def sigJ(b):
        M = np.zeros((DS, DS), dtype=complex)
        for (s, ds), o in zip(sigl, offs[:-1]): M[o:o + ds, o:o + ds] = s(b)
        return M
    mi = minv(L.m); D = D4(L)
    rng = random.Random(seed + 1)
    base = [b for (b, s) in L.Bm]; Bm = list(base); seen = set(key(b) for b in Bm)
    while len(Bm) < 400:
        w = mm(rng.choice(base), rng.choice(Bm))
        if key(w) not in seen: seen.add(key(w)); Bm.append(w)
    Pi = lambda k: D.R(k, sig_triv, 1).real
    def exact7(b):
        s = sigma(b); o = L.sp_order(s)
        if o <= 0 or o % 7: return None
        qq = snap(np.linalg.matrix_power(b, o // 7)); M7 = np.linalg.matrix_power(qq, 7); a2 = 0
        while not np.allclose(M7, I16):
            M7 = snap(M7 @ M7); a2 += 1
            if a2 > 10: return None
        e = (2 ** a2) * pow(2 ** a2, -1, 7)
        qf = snap(np.linalg.matrix_power(qq, e))
        return qf if np.allclose(np.linalg.matrix_power(qf, 7), I16) and not np.allclose(qf, I16) else None
    q1 = L.q; pow1 = set(key(x) for x in L.qp); q2 = None
    for b in Bm:
        qf = exact7(b)
        if qf is not None and key(qf) not in pow1 and not any(np.allclose(sigma(qf), sigma(x)) for x in L.qp):
            P1 = Pi(mm(L.m, q1, mi)); P2 = Pi(mm(L.m, qf, mi))
            if not any(np.allclose(np.linalg.matrix_power(P1, j), P2) for j in range(7)):
                q2 = qf; break
    q2p = [I16]
    for j in range(1, 7): q2p.append(mm(q2p[-1], q2))
    p2 = [mm(L.m, x, mi) for x in q2p]
    return L, D, sigJ, DS, L.p, p2
