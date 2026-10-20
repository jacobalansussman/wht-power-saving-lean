"""t9_galg.py (session 2) -- O O F O F (Z7, 7 in 23) with a middle content that is a GROUP-ALGEBRA element:
F1 content y = sum_i a_i l_i t (letters l_i in Stab(L_Z) cap Stab(Lam_m) on the one O2 term t), so F1's top part in character
chi is V_chi = Y W_chi with Y = sum_i a_i X_i in A = span(perm matrices of G on nbrs(Lam_F)) and W_chi free in S_chi
(span of the 4 regular <pt>-orbit vectors, weights w^(a h), a = chi - psi1).  The relaxation needs, for every effective
character a in Xi:  e_K = Y w_a + f_a,  w_a in S_a,  f_a in Fr (free span on the 30 slots).
Test: random w_a, then a linear solve for (Y, f_a).  Free span: one O1 cell (main) -> nbrs(Lam_F) cap nbrs(Lam_1) plus
the fixed slot Y_line (reached by O2's items)."""
import sys, os, random, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from sp8 import *
from collections import Counter
t0 = time.time()
P = Polar(); D = P.D; iZ = P.idx[LZ]
iK = P.idx[img_mask(table(inv(sigma(G.K))), LZ)]
rng = random.Random(3); nrng = np.random.default_rng(5)
gensM = [sigma(G.D[z]) for z in range(1, 16)] + [sigma(G.rand_engine(rng)) for _ in range(30)]
trans = {iZ: IDM}; fr = [iZ]
while fr:
    nf = []
    for i in fr:
        for g in gensM:
            s2 = mul(g, trans[i]); j = P.idx[span_mask([app(s2, 1 << (4 + k)) for k in range(4)])]
            if j not in trans: trans[j] = s2; nf.append(j)
    fr = nf
letters = [sigma(G.rand_engine(rng)) for _ in range(4000)]
F = [j for j in P.nb[iK] if P.dZ[j] == 3][0]
nbF = P.nb[F]; ni = {Y: k for k, Y in enumerate(nbF)}
T = [Y for Y in nbF if P.dZ[Y] == 2 and D[Y, iK] == 2]
t = T[int(sys.argv[1]) if len(sys.argv) > 1 else 0]
def pm(a, b): return tuple(a[b[i]] for i in range(len(b)))
def pinv(a):
    r = [0] * len(a)
    for x, y in enumerate(a): r[y] = x
    return tuple(r)
s0 = trans[t]; s0i = inv(s0)
Fp = P.idx[img_mask(table(s0i), P.L[F])]
gens = []
for b in letters:
    if P.idx[img_mask(table(b), P.L[Fp])] == Fp:
        g = mul(s0, mul(b, s0i)); Tg = table(g)
        gens.append(tuple(ni[P.idx[img_mask(Tg, P.L[Y])]] for Y in nbF))
    if len(gens) >= 60: break
idp = tuple(range(30)); grp = {idp}; fr = [idp]
while fr:
    nf = []
    for g in fr:
        for s in gens:
            h = pm(s, g)
            if h not in grp: grp.add(h); nf.append(h)
    fr = nf
grpl = sorted(grp)
sev = []
for g in grpl:
    h = g; k = 1
    while h != idp and k < 8: h = pm(g, h); k += 1
    if k == 7: sev.append(g)
cls = {}; reps = []
for g in sev:
    if g in cls: continue
    reps.append(g)
    for x in grpl: cls.setdefault(pm(x, pm(g, pinv(x))), len(reps) - 1)
# basis of A: greedy group elements with rank increase
def pmat(g):
    M = np.zeros((30, 30)); M[list(g), list(range(30))] = 1; return M
order = list(grpl); random.Random(7).shuffle(order)
basis = []; Bflat = np.zeros((0, 900))
for g in order:
    v = pmat(g).reshape(-1)
    cand = np.vstack([Bflat, v])
    if np.linalg.matrix_rank(cand, tol=1e-9) > Bflat.shape[0]:
        Bflat = cand; basis.append(g)
    if len(basis) >= 900: break
dA = len(basis)
print('content cell', t, ': |G| =', len(grp), ' dim A (span of the 30x30 permutation matrices of G) =', dA, ' %.0fs' % (time.time() - t0), flush=True)
fixed = [s for s in range(30) if all(g[s] == s for g in gens)]
print('slots fixed by G:', [nbF[s] for s in fixed], ' input distances', [int(P.dZ[nbF[s]]) for s in fixed], flush=True)
k0 = ni[iK]
w = np.exp(2j * np.pi / 7)
def ovecs(g, a):
    gi = pinv(g); seen = set(); vs = []
    for N in range(30):
        if N in seen: continue
        orb = [N]; cur = N
        for _ in range(6): cur = gi[cur]; orb.append(cur)
        seen |= set(orb)
        if len(set(orb)) == 1: continue
        v = np.zeros(30, dtype=complex)
        for h, c in enumerate(orb): v[c] += w ** (a * h)
        vs.append(v)
    return np.array(vs).T
L1c = [Y for Y in P.nb[t] if P.dZ[Y] == 1]
Bm = [pmat(g) for g in basis]
e = np.zeros(30, dtype=complex); e[k0] = 1
for lam1 in L1c:
    Fr_slots = sorted(set([ni[c] for c in nbF if D[c, lam1] == 1]) | set(fixed))
    for c_ in range(len(reps)):
        pt = reps[c_]
        for excl in [None] + list(range(1, 7)):
            Xi = [a for a in range(1, 7) if a != excl]
            best = 9e9
            for trial in range(3):
                rows = []; rhs = []
                nF = len(Fr_slots)
                ncol = dA + nF * len(Xi)
                Mbig = np.zeros((30 * len(Xi), ncol), dtype=complex); bb = np.zeros(30 * len(Xi), dtype=complex)
                for ia, a in enumerate(Xi):
                    S = ovecs(pt, a)
                    wa = S @ (nrng.standard_normal(S.shape[1]) + 1j * nrng.standard_normal(S.shape[1]))
                    for k, B in enumerate(Bm):
                        Mbig[30 * ia:30 * ia + 30, k] = B @ wa
                    for j, s in enumerate(Fr_slots):
                        Mbig[30 * ia + s, dA + nF * ia + j] = 1
                    bb[30 * ia:30 * ia + 30] = e
                x, *_ = np.linalg.lstsq(Mbig, bb, rcond=None)
                res = float(np.linalg.norm(Mbig @ x - bb))
                best = min(best, res)
            print(f'  O1 cell {lam1} class {c_} characters {Xi}: residual {best:.2e}', flush=True)
print('done %.0fs' % (time.time() - t0))
