# vlift.py -- v-w7-mixed3: is w7-mixed3's "two letters can first-touch a two-row array at exactly 1.0 in the 105
# (Lemma Q3 attained)" real?  vn105.py reports min over j in {nbar, 2 nbar, 4 nbar} (j = least power of the relating
# letter g that is a scalar ON V); the value 1.0 needs j = 2.  Here: find every pair (b, g), g in P0 of order 2, with
# A_{bg} proportional to A_b (rank 2; then delta_2 = 1), and test whether g has a Clifford lift C (any Pauli multiple of
# one lift) with C^2 scalar on V = C^8.  Lifts: A-generators t_(0|w) -> diag(i^(x.w mod 2)); Levi words -> CNOT-type
# permutations, each lift checked against the symplectic action on Pauli labels.
import itertools, sys, collections, numpy as np
src = open(__file__.rsplit('/', 1)[0] + '/v105.py').read()
exec(src.split('# ---------------- exact part, mod p ----------------')[0])
# ---- Pauli labels and Clifford lifts ----
def pauli(v):
    x, z = v & 7, v >> 3; M = np.zeros((8, 8), complex)
    for b in range(8): M[b ^ x, b] = (-1) ** bin(b & z).count('1')     # X^x Z^z : first Z^z then X^x
    return M
PA = [pauli(v) for v in range(64)]
def sympl_of(C):
    Ci = np.linalg.inv(C); tab = []
    for v in range(64):
        M = C @ PA[v] @ Ci; hit = None
        for w in range(64):
            ip = np.vdot(PA[w], M)
            if abs(abs(ip) - 8) < 1e-6: hit = w; break
        assert hit is not None; tab.append(hit)
    return tuple(tab)
def table(f): return tuple(f(v) for v in range(64))
def diag_lift(w):     # t_(0|w): phase i^(x.w mod 2)
    return np.diag([1j ** (bin(b & (w >> 3)).count('1') & 1) for b in range(8)])
A_lift = [diag_lift(w[0]) for w in A_ws]
for Cl, w in zip(A_lift, A_ws):
    assert sympl_of(Cl) == table(word_fn(w)), "A lift mismatch"
def perm_lift_for(tab):   # find a permutation-type lift x -> M x (x in F2^3) matching the symplectic table on x-part
    for cols in itertools.product(range(8), repeat=3):
        P = np.zeros((8, 8), complex)
        ok = True
        imgs = []
        for b in range(8):
            y = 0
            for j in range(3):
                if b >> j & 1: y ^= cols[j]
            imgs.append(y)
        if len(set(imgs)) != 8: continue
        for b in range(8): P[imgs[b], b] = 1
        for sgn in itertools.product(range(4), repeat=0):
            if sympl_of(P) == tab: return P
    return None
L_lift = []
for w in L_ws:
    tab = table(word_fn(w)); P = perm_lift_for(tab); assert P is not None, w; L_lift.append(P)
print('lifts of the 6 A-generators and 6 Levi generators verified against the symplectic action', flush=True)
# ---- P0 with lifts: element = a(bits) . l(BFS word) ----
Aall_l = []
for bits in range(64):
    Cm = np.eye(8, dtype=complex)
    for j in range(6):
        if bits >> j & 1: Cm = Cm @ A_lift[j]
    Aall_l.append(Cm)
Lev_l = [np.eye(8, dtype=complex)]; Lev7b = [np.eye(7, dtype=np.int64)]; kk = {Lev7b[0].tobytes()}; fr = [0]
while fr:
    nx = []
    for i in fr:
        for l7, ll in zip(L7, L_lift):
            g7 = l7 @ Lev7b[i]
            if g7.tobytes() not in kk: kk.add(g7.tobytes()); Lev7b.append(g7); Lev_l.append(ll @ Lev_l[i]); nx.append(len(Lev7b) - 1)
    fr = nx
assert len(Lev7b) == 168
P7 = [a @ l for a in Aall7 for l in Lev7b]; PL = [a @ l for a in Aall_l for l in Lev_l]
idx = {g.tobytes(): i for i, g in enumerate(P7)}; assert len(idx) == 10752
# spot-check the lift homomorphism on random products
rng = np.random.default_rng(5)
for _ in range(20):
    i, j = rng.integers(10752, size=2); k = idx[(P7[i] @ P7[j]).tobytes()]
    M = PL[i] @ PL[j] @ np.linalg.inv(PL[k])     # must be a Pauli times a scalar
    assert any(abs(abs(np.vdot(PA[w], M)) - 8) < 1e-6 for w in range(64))
print('lift map spot-checked (products agree up to Paulis and scalars)', flush=True)
# ---- float defect maps, two-row slot, orthonormal model (from v105's float part) ----
fl = src.split('# ---------------- float part: stable ranks, orthonormalised model ----------------')[1].split('for nm, E_ in Ef_.items():')[0]
fl = fl.replace('assert np.abs(Zf - Zf.T).max() < 1e-8', '')
# O28 was computed in the mod-p part; recompute here from characters (same orbit rule)
Aidx = {a.tobytes(): b for b, a in enumerate(Aall7)}
def act(l, s):
    li = np.round(np.linalg.inv(l)).astype(np.int64); out = 0
    for j in range(6):
        c = Aidx[(li @ A7[j] @ l).tobytes()]; out |= (bin(s & c).count('1') & 1) << j
    return out
orbs, seen = [], set()
for s in range(64):
    if s in seen: continue
    orb = {s}; frn = [s]
    while frn:
        nx = []
        for x in frn:
            for l in L7:
                y = act(l, x)
                if y not in orb: orb.add(y); nx.append(y)
        frn = nx
    seen |= orb; orbs.append(sorted(orb))
O28 = [o for o in orbs if len(o) == 28][0]
exec(fl)
ER_ = Ef_['two-row']
Aflat = np.zeros((10752, D * ER_.shape[1])); rks = np.zeros(10752, int)
W1P = [W1f(g) for g in P7]
for i in range(10752):
    W = hf @ W1P[i] @ ER_; W = W - PN @ W
    s = np.linalg.svd(W, compute_uv=False); rks[i] = int((s > 1e-7).sum()); Aflat[i] = W.ravel()
print('two-row defect counts', dict(collections.Counter(rks.tolist())), flush=True)
nrm = np.linalg.norm(Aflat, axis=1); An = Aflat / np.where(nrm > 0, nrm, 1)[:, None]
r2 = np.nonzero(rks == 2)[0]
def order7(g):
    P = np.eye(7, dtype=np.int64)
    for e in range(1, 30):
        P = P @ g
        if (P == np.eye(7, dtype=np.int64)).all(): return e
inv = [i for i in range(10752) if order7(P7[i]) == 2]
print('involutions in P0:', len(inv), flush=True)
found = collections.Counter(); ex = {}
for gi in inv:
    tgt = np.array([idx[(P7[b] @ P7[gi]).tobytes()] for b in r2])
    ok = (rks[tgt] == 2) & (np.abs(np.einsum('ij,ij->i', An[r2], An[tgt])) > 1 - 1e-9)
    if ok.any(): found[gi] = int(ok.sum()); ex[gi] = int(r2[np.nonzero(ok)[0][0]])
print('order-2 relating letters g with some b: A_{bg} proportional to A_b (delta_2 = 1):', len(found), 'g, pairs', sum(found.values()), flush=True)
# lift test: is there a Pauli P with (C_g P)^2 scalar?
res = collections.Counter(); good = []
for gi in found:
    Cg = PL[gi]; okl = False
    for w in range(64):
        M = Cg @ PA[w]; M2 = M @ M
        if np.allclose(M2, M2[0, 0] * np.eye(8), atol=1e-9): okl = True; break
    res[okl] += 1
    if okl: good.append(gi)
print('of these g, lifts with C^2 scalar exist for:', res[True], '; none for:', res[False], flush=True)
if good:
    gi = good[0]; b = ex[gi]
    print('example: b index', b, 'g index', gi, '-> Lemma Q3 bound 2/2 = 1.0 attained by D(b + c b g) with a j = 2 lift')
