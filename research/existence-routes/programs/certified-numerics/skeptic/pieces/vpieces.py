# v-ex2-certnum (skeptic): INDEPENDENT recomputation of ex2-certnum Computations 12-14 in a different model of Sp(6,2):
# W(E7)^+ via e7model.py (written by skeptic v-meta2-finite-group-compute; imports nothing from c2-finite-group or
# ex2-certnum).  Float linear algebra with reported singular-value gaps (evidence, not proof); I check: the irreducible
# chi as the Z2a-eigenspace of a tensor module (multiplicity one), its P_0-decomposition (A-weight orbits, commutant
# splitting), irreducibility of each piece (commutant dim 1), and the block ranks of h (distance 1), h2 (2), k (3).
import sys, itertools, numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from e7model import *
np.set_printoptions(linewidth=180)
A_ws = [[z_(1)], [z_(2)], [z_(3)], [z_(1) | z_(2)], [z_(1) | z_(3)], [z_(2) | z_(3)]]
def table(f): return tuple(f(v) for v in range(64))
def levi_fn(i, j):
    def f(v):
        x, z = v & 7, v >> 3
        if (x >> j) & 1: x ^= 1 << i
        if (z >> i) & 1: z ^= 1 << j
        return x | (z << 3)
    return f
L_ws = []
for i, j in itertools.permutations(range(3), 2):
    w = [1 << i, 1 << (3 + j), (1 << i) | (1 << (3 + j))]; assert table(word_fn(w)) == table(levi_fn(i, j)); L_ws.append(w)
L0 = frozenset(v for v in range(64) if v & 7 == 0)
for w in A_ws + L_ws: assert frozenset(word_fn(w)(v) for v in L0) == L0
words = {'h': [ee(1)], 'h2': [ee(1), ee(2)], 'k': [ee(1), ee(2), ee(3)]}
for nm, w in words.items():
    print(nm, 'dim(L0 cap g L0) =', len(frozenset(word_fn(w)(v) for v in L0) & L0).bit_length() - 1)
pairs = list(itertools.combinations(range(7), 2))
def wedge2(R):
    out = np.zeros((21, 21), dtype=np.int64)
    for a, (i, j) in enumerate(pairs):
        for b, (k, l) in enumerate(pairs): out[a, b] = R[i, k] * R[j, l] - R[i, l] * R[j, k]
    return out
def sym2(R):
    pr = [(i, j) for i in range(7) for j in range(i, 7)]; out = np.zeros((28, 28), dtype=np.int64)
    for a, (i, j) in enumerate(pr):
        for b, (k, l) in enumerate(pr):
            out[a, b] = R[i, k] * R[j, l] + (R[j, k] * R[i, l] if i != j else 0)
    return out
MODS = {'L2S2': lambda R: np.kron(wedge2(R), sym2(R)), '77': lambda R: np.kron(R, R), '777': lambda R: np.kron(np.kron(R, R), R),
        'L2L2': lambda R: np.kron(wedge2(R), wedge2(R)), '7777': lambda R: np.kron(np.kron(R, R), np.kron(R, R))}
mod, dimchi, chi2a = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
rho = MODS[mod]
z = 63 * chi2a / dimchi; assert abs(z - round(z)) < 1e-12; z = round(z)
Dm = rho(np.eye(7, dtype=np.int64)).shape[0]
Z = np.zeros((Dm, Dm))
for w in range(1, 64): Z += rho(R7_t(w))
ev = np.round(np.linalg.eigvals(Z).real).astype(int); vals, cnt = np.unique(ev, return_counts=True)
print('module', mod, 'dim', Dm, 'Z2a spectrum', dict(zip(vals.tolist(), cnt.tolist())))
if dict(zip(vals.tolist(), cnt.tolist())).get(z) == dimchi:
    U, s, Vh = np.linalg.svd(Z - z * np.eye(Dm)); X = Vh[-dimchi:].conj().T   # null space basis (orthonormal, Dm x d)
    print('null-space gap of Z - zI: s[-d-1] = %.3e, s[-d] = %.3e' % (s[-dimchi - 1], s[-dimchi]))
else:   # refine by the class sum of 2b (t_u t_v t_{u+v}, {0,u,v,u+v} totally isotropic): pick the joint eigenspace of dim d
    mz = dict(zip(vals.tolist(), cnt.tolist())).get(z, 0)
    U, s, Vh = np.linalg.svd(Z - z * np.eye(Dm)); Y = Vh[-mz:].conj().T
    print('Z2a eigenspace for %d has dim %d; refining by Z2b' % (z, mz))
    lines = set()
    for u in range(1, 64):
        for v in range(1, 64):
            if u < v and Bl(u, v) == 0 and (u ^ v) > v: lines.add((u, v))
    print('isotropic lines:', len(lines))
    Zb = np.zeros((Dm, Dm))
    for (u, v) in lines: Zb += rho(word_R7([u, v, u ^ v]))
    Zb_Y = np.linalg.pinv(Y) @ Zb @ Y
    evb = np.round(np.linalg.eigvals(Zb_Y).real, 6); vb, cb = np.unique(np.round(evb).astype(int), return_counts=True)
    print('Z2b spectrum on it:', dict(zip(vb.tolist(), cb.tolist())))
    zb = [v for v, c in zip(vb.tolist(), cb.tolist()) if c == dimchi]; assert len(zb) == 1, 'no unique d-dim joint eigenspace'
    U2, s2, Vh2 = np.linalg.svd(Zb_Y - zb[0] * np.eye(mz)); X = Y @ Vh2[-dimchi:].conj().T
    X, _ = np.linalg.qr(X)
    print('joint eigenspace (z2a, z2b) = (%d, %d), dim %d; gap %.3e vs %.3e' % (z, zb[0], dimchi, s2[-dimchi - 1], s2[-dimchi]))
Xp = np.linalg.pinv(X)
def on_X(Mbig):
    Y = Mbig @ X; Mx = Xp @ Y; res = np.abs(X @ Mx - Y).max(); assert res < 1e-8, res; return Mx
gen = {nm: on_X(rho(word_R7(w)).astype(float)) for nm, w in words.items()}
Ar = [on_X(rho(word_R7(w)).astype(float)) for w in A_ws]; Lr = [on_X(rho(word_R7(w)).astype(float)) for w in L_ws]
tr_t = np.trace(on_X(rho(R7_t(ee(1))).astype(float))); print('character of a transvection on X: %.6f (want %d)' % (tr_t, chi2a))
d = dimchi; I = np.eye(d)
def Pi(sv):
    P = I.copy()
    for i in range(6): P = P @ (I + (-1 if sv >> i & 1 else 1) * Ar[i]) / 2
    return P
# A-character orbits under the Levi (via conjugation action on A)
def compose_tab(a, b): return tuple(a[b[x]] for x in range(64))
Atab = {}
for c in range(64):
    g = tuple(range(64))
    for i in range(6):
        if c >> i & 1: g = compose_tab(g, table(word_fn(A_ws[i])))
    Atab[g] = c
def inv_tab(a):
    r = [0] * 64
    for x in range(64): r[a[x]] = x
    return tuple(r)
Ltabs = [table(word_fn(w)) for w in L_ws]
def act(l, sv):
    li = inv_tab(l); out = 0
    for j in range(6):
        aj = table(word_fn(A_ws[j])); c = Atab[compose_tab(li, compose_tab(aj, l))]
        out |= (bin(sv & c).count('1') & 1) << j
    return out
orbits, seen = [], set()
for sv in range(64):
    if sv in seen: continue
    orb = {sv}; fr = [sv]
    while fr:
        nx = []
        for x in fr:
            for l in Ltabs:
                y = act(l, x)
                if y not in orb: orb.add(y); nx.append(y)
        fr = nx
    seen |= orb; orbits.append(sorted(orb))
def colbasis(P, tol=1e-8):
    U, s, _ = np.linalg.svd(P); r = int((s > tol * max(1, s[0])).sum()); return U[:, :r], r, s
P0 = Ar + Lr
def commutant(B):      # commutant of P_0 restricted to the invariant subspace spanned by columns of B (d x m)
    Bp = np.linalg.pinv(B); m = B.shape[1]
    gs = []
    for g in P0:
        G = Bp @ g @ B; assert np.abs(B @ G - g @ B).max() < 1e-8, 'not invariant'; gs.append(G)
    Eq = np.concatenate([np.kron(np.eye(m), G.T) - np.kron(G, np.eye(m)) for G in gs], 0)
    U, s, Vh = np.linalg.svd(Eq); nul = int((s < 1e-8 * s[0]).sum()) + max(0, m * m - len(s))
    gap = (s[len(s) - nul - 1] if nul < len(s) else None, s[len(s) - nul] if nul else None)
    basis = [Vh[-i - 1].conj().reshape(m, m) for i in range(nul)]
    return nul, gap, basis, gs
pieces = []
for orb in orbits:
    PO = sum(Pi(sv) for sv in orb)
    B, r, _ = colbasis(PO)
    if r == 0: continue
    nul, gap, basis, gs = commutant(B)
    print('orbit size %d: part dim %d, commutant dim %d (singular values around the gap: %s)' % (len(orb), r, nul, gap))
    if nul == 1: pieces.append(B); continue
    rng = np.random.default_rng(1); C = sum(rng.normal() * b for b in basis)
    w, V = np.linalg.eig(C); order = np.argsort(w.real); w = w[order]; V = V[:, order]
    groups = []; st = 0
    for i in range(1, r + 1):
        if i == r or abs(w[i] - w[i - 1]) > 1e-6: groups.append((st, i)); st = i
    print('   splits by a generic commutant element into dims', [b - a for a, b in groups])
    for a, b in groups:
        Bi = B @ V[:, a:b]; Q, _ = np.linalg.qr(Bi)
        n2, g2, _, _ = commutant(Q); print('   piece dim %d: commutant dim %d (gap %s)' % (b - a, n2, g2)); pieces.append(Q)
for Q in pieces:
    if Q.shape[1] and commutant(Q)[0] != 1: print('   !! piece of dim', Q.shape[1], 'NOT irreducible')
pieces.sort(key=lambda Q: Q.shape[1])
dims = [Q.shape[1] for Q in pieces]; print('pieces (dims):', dims, 'sum', sum(dims))
T = np.concatenate(pieces, 1); Ti = np.linalg.inv(T)
off = np.cumsum([0] + dims); idx = [list(range(off[i], off[i + 1])) for i in range(len(dims))]
def rk(M):
    if M.size == 0: return 0, None
    s = np.linalg.svd(M, compute_uv=False); r = int((s > 1e-8 * max(1, s[0])).sum())
    return r, (s[r - 1] if r else None, s[r] if r < len(s) else 0.0)
for nm in ('h', 'h2', 'k'):
    Mt = Ti @ gen[nm] @ T
    br = [[rk(Mt[np.ix_(idx[a], idx[b])])[0] for b in range(len(dims))] for a in range(len(dims))]
    rows = []; gaps = []
    for a in range(len(dims)):
        others = [c for b in range(len(dims)) if b != a for c in idx[b]]
        r_, g_ = rk(Mt[np.ix_(idx[a], others)]); rows.append(r_); gaps.append(g_)
    print('%s: block ranks %s; rank of [tau, others] per piece %s; gaps %s' % (nm, br, rows,
          ['(%.2e, %.2e)' % g for g in gaps if g and g[0] is not None]))
