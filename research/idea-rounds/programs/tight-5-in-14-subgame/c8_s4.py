# c8_s4.py -- the order-24 survivors: which S4 of L, its 9 shapes in detail; then the same screen inside the singular
# subgame (letters b d b' with b, b' in B' = AGL(3,2)) for the hook AND for other members Y(T) of the tight family.
import numpy as np, itertools, time
from lib514 import *
import c7_symscreen2 as S7
from libs8 import NS
t0 = time.time()
Lsubs = S7.subgroups_upto_conj(LEVI)
S4s = [H for H in Lsubs if len(H) == 24]
def fixed_subspaces(H):
    """which nonzero vectors / which planes of L_Z = z-space does H fix (z-part of the Levi action)"""
    pts = [z for z in range(1, 8) if all(int(h[z << 3]) == (z << 3) for h in H)]
    planes = []
    for a, b in itertools.combinations(range(1, 8), 2):
        P = {0, a, b, a ^ b}
        if all({int(h[z << 3]) >> 3 for z in P} == P for h in H): planes.append(tuple(sorted(P)))
    return pts, sorted(set(planes))
if __name__ == '__main__':
    for H in S4s:
        info = S7.analyse(H, 'S4')
        pts, planes = fixed_subspaces(H)
        res, lines = S7.describe(info, 5, 14, tight=True, gauge=False)
        print(f'S4 fixing z-vectors {pts} and z-planes {planes}: {sum(len(s) for _, _, s in res)} exact shapes')
        if res:
            print('   irreps', info['dims'])
            for a in info['arr']: print('   array orbit |S|=%d size %d demand %s' % (a['S'], a['size'], a['mult']))
            for m in info['mov']: print('   move orbit |T|=%d size %d cap %s exact=%s' % (m['T'], m['size'], m['mult'], m['exact']))
            for l in lines: print('     ', l)
print('[%.0fs]' % (time.time() - t0), flush=True)

# ---- subgame version: B' instead of P0, and the whole tight family Y(t7, t21)
R = S7.R
alt = [g for g in gen_closure(A_GENS) if all(q0(int(g[v])) == q0(v) for v in range(NV))]
Bp_gens = [levi([[1, 1, 0], [0, 1, 0], [0, 0, 1]]), levi([[0, 0, 1], [1, 0, 0], [0, 1, 0]])] + [a for a in alt if not np.array_equal(a, ID)][:3]
Bp = gen_closure(Bp_gens); assert len(Bp) == 1344
d = h1; dinv = inv(d)
dLZ = frozenset(int(d[v]) for v in LZ)
Ep = [b for b in Bp if frozenset(int(b[v]) for v in dLZ) == dLZ]; Epk = {key(b) for b in Ep}
print("|B'| =", len(Bp), " |E'| =", len(Ep))
# S8-adapted pieces of V56: A=(2,1^6), B=(3,1^5), C=(2,2,2,1,1) via the transposition class sum (eigenvalues -20,-12,-10)
Zt = sum(R.r56(TRANS[v]) for v in NS); Zt = (Zt + Zt.T) / 2
w, V = np.linalg.eigh(Zt)
Ab = V[:, np.abs(w + 20) < 1e-6]; Bb = V[:, np.abs(w + 12) < 1e-6]; Cb = V[:, np.abs(w + 10) < 1e-6]
print('dims A, B, C:', Ab.shape[1], Bb.shape[1], Cb.shape[1])
# B'-isotypic split of C into C7, C21, and B'-isomorphisms Phi7: A -> C7, Phi21: B -> C21
mats = [R.r56(b) for b in Bp]
def proj_type(Q, ref):
    """component of span(Q) isomorphic to the B'-module span(ref): use the B'-averaged intertwiners"""
    rng = np.random.default_rng(3)
    X = rng.standard_normal((Q.shape[1], ref.shape[1]))
    T = sum((Q.T @ M @ Q) @ X @ (ref.T @ M @ ref).T for M in mats) / len(mats)
    return T
T7 = proj_type(Cb, Ab); T21 = proj_type(Cb, Bb)
u7, s7, _ = np.linalg.svd(T7); u21, s21, _ = np.linalg.svd(T21)
C7 = Cb @ u7[:, :7]; C21 = Cb @ u21[:, :21]
print('C7 ⟂ C21:', np.abs(C7.T @ C21).max() < 1e-8)
def intertwiner(Q1, Q2):
    eqs = [np.kron(np.eye(Q1.shape[1]), (Q1.T @ M @ Q1).T) - np.kron(Q2.T @ M @ Q2, np.eye(Q1.shape[1])) for M in [R.r56(g) for g in Bp_gens]]
    _, sv, vt = np.linalg.svd(np.vstack(eqs)); X = vt[-1].reshape(Q1.shape[1], Q1.shape[1])
    X = X / np.sqrt(np.trace(X @ X.T) / Q1.shape[1]); assert sv[-1] < 1e-8 and np.allclose(X @ X.T, np.eye(Q1.shape[1]), atol=1e-8)
    return Q2 @ X          # columns: image of Q1's basis in span(Q2), B'-equivariantly
C7a = intertwiner(Ab, C7); C21b = intertwiner(Bb, C21)
def Yfam(t7, t21):
    return np.hstack([Ab + t7 * C7a, Bb + t21 * C21b]) / np.sqrt(1 + 0.0)
# locate the hook in the family
Yh = S7.Y
for nm, (P, Q) in {'7': (Ab, C7a), '21': (Bb, C21b)}.items():
    a = np.linalg.norm(P.T @ Yh); c = np.linalg.norm(Q.T @ Yh)
    M1 = P.T @ Yh; M2 = Q.T @ Yh
    t = np.sum(M1 * M2) / np.sum(M1 * M1)
    print(f'hook type {nm}: t = {t:.6f}, residual {np.linalg.norm(M2 - t * M1):.2e}')
def check_fam(t7, t21):
    Y = Yfam(t7, t21); Q, _ = np.linalg.qr(Y)
    f = lambda g: S7.rank(np.hstack([Q, R.r56(g) @ Q])) - 28
    return f(h1), f(K)
print('family members (cap, demand):', [check_fam(a, b) for a, b in [(0.3, -1.7), (2.0, 0.5), (-0.4, -0.4), (1.0, 1.0)]])

def cap_char_family(Q, xs):
    """trace of r(x) on (Q + dQ)/Q for x in E' (Q orthonormal 28-dim B'-invariant)"""
    PQ = Q @ Q.T
    u, s, _ = np.linalg.svd((np.eye(56) - PQ) @ R.r56(d) @ Q); Wq = u[:, :int((s > 1e-8).sum())]
    QW = np.hstack([Q, Wq]); out = []
    for x in xs:
        M = R.r56(x); C = np.linalg.lstsq(QW, M @ QW, rcond=None)[0]
        assert np.abs(QW @ C - M @ QW).max() < 1e-8
        out.append(np.trace(C[28:, 28:]))
    return out, Wq.shape[1]

def analyse_subgame(H, Q):
    """like S7.analyse but letters from B' and certificate Q (B'-invariant 28-dim); demand = V56/Q ~ complement"""
    n = len(H); cls, rep, sizes, chars, clsof, idx = S7.char_table(H)
    ip = lambda a, b: np.real(np.sum(sizes * a * np.conj(b)) / n)
    Tcls = S7.subgroups_upto_conj(H)
    PQ = Q @ Q.T
    def comp_trace(x): M = R.r56(x); return np.trace(M) - np.trace(Q.T @ M @ Q)
    arr = []
    for S in Tcls:
        ind = S7.induced(H, rep, S, [comp_trace(x) for x in S])
        arr.append(dict(S=len(S), size=n // len(S), mult=np.round([ip(ind, c) for c in chars], 6)))
    mov = []
    for T in Tcls:
        caps = {}; seenconj = set()
        for b in Bp:
            bi = inv(b); xs = [comp(comp(bi, t), b) for t in T]
            kx = tuple(key(x) for x in xs)
            if kx in seenconj: continue
            seenconj.add(kx)
            if not all(key(x) in Epk for x in xs): continue
            chv, wdim = cap_char_family(Q, xs)
            ch = tuple(np.round(chv, 6))
            tg = [comp(comp(dinv, x), d) for x in xs]
            ex = caps.get(ch, False)
            if not ex:
                ex = any(all(np.array_equal(comp(comp(bp, x), inv(bp)), y) for x, y in zip(xs, tg)) for bp in Bp)
            caps[ch] = ex
        for ch, ex in caps.items():
            ind = S7.induced(H, rep, T, list(ch))
            mov.append(dict(T=len(T), size=n // len(T), mult=np.round([ip(ind, c) for c in chars], 6), exact=ex))
    dims = [int(round(np.real(c[[i for i, cc in enumerate(cls) if np.array_equal(H[cc[0]], ID)][0]]))) for c in chars]
    return dict(n=n, dims=dims, arr=arr, mov=mov)

if __name__ == '__main__':
    fam = {"hook": S7.Y}
    rng = np.random.default_rng(5)
    for i in range(3):
        t7, t21 = rng.normal(size=2) * 1.3
        Qf, _ = np.linalg.qr(Yfam(t7, t21)); fam[f'Y({t7:.2f},{t21:.2f})'] = Qf
    for H in [h for h in Lsubs if len(h) in (24, 12, 8, 6)]:
        shapes = {}
        for nm, Q in fam.items():
            info = analyse_subgame(H, Q)
            res, lines = S7.describe(info, 5, 14, tight=True, gauge=False)
            shapes[nm] = set(lines)
        common = set.intersection(*shapes.values())
        pts, planes = fixed_subspaces(H)
        print(f'SUBGAME H order {len(H)} (fixes z-vectors {pts}, z-planes {planes}): exact shapes per certificate',
              {k: len(v) for k, v in shapes.items()}, ' common to all:', len(common), '[%.0fs]' % (time.time() - t0), flush=True)
        for l in sorted(common)[:10]: print('      ', l)
