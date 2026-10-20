"""t7c_pairchars.py (session 2) -- t7 per character, for EVERY pair of eigencharacters (psi1, psi2), with a positive control.
Setting as t7_pairorbits.py: two single-term shared moves F1, F2 on the same O2 content cell Lam_t and the same F cell
Lam_F (Z7, 7 in 23), twists pt1, pt2 = all labelled order-7 elements of Stab(Lam_t) cap Stab(Lam_F) (image on the 30
neighbours of Lam_F), O1 cell Lam_1 = each of the 3 neighbours of Lam_t at input distance 1.
In character chi the top's items on F_i give, per <pt_i>-orbit, the vector sum_h c_i(h) e_{pt_i^-h N} with a free
per-chi coefficient, where c_i = chi - psi_i is F_i's effective character (nonzero, since chi is outside E).  So chi is
served iff e_{L_1} is in span(orbit vectors of pt1 at index a, of pt2 at index b) + free span, (a, b) = (chi - psi1,
chi - psi2).  Galois relabelling h -> kh maps (pt1, pt2, a, b) to (pt1^k, pt2^k, a/k, b/k) and permutes the labelled
twists, so a = 1 and b = 1..6 cover every (psi1, psi2, chi) (b = 1 is t7's case psi1 = psi2).
mode 'main': free span = e_c, c in nbrs(Lam_F) cap nbrs(Lam_1) (O1's and O2's items; t7's model).
mode 'ctl' : free span = every slot of nbrs(Lam_F) at input distance <= 3 (positive control: then only the distance-4
             junk must cancel, which the shift-matched pairs of t5 do)."""
import sys, os, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sp8 import *
from collections import Counter
t0 = time.time()
part = int(sys.argv[1]); nparts = int(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else 'main'
P = Polar(); D = P.D; iZ = P.idx[LZ]
iK = P.idx[img_mask(table(inv(sigma(G.K))), LZ)]
rng = random.Random(3)
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
t = T[0]
def pm(a, b): return tuple(a[b[i]] for i in range(len(b)))
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
sev = []
for g in sorted(grp):
    h = g; k = 1
    while h != idp and k < 8: h = pm(g, h); k += 1
    if k == 7: sev.append(g)
print('mode', mode, '; content cell', t, 'F cell', F, 'image order', len(grp), 'labelled twists', len(sev), ' %.0fs' % (time.time() - t0), flush=True)
k0 = ni[iK]
w = np.exp(2j * np.pi / 7)
def orbit_vecs(g, chi):
    gi = [0] * 30
    for x, y in enumerate(g): gi[y] = x
    seen = set(); vecs = []
    for N in range(30):
        if N in seen: continue
        orb = [N]; cur = N
        for _ in range(6): cur = gi[cur]; orb.append(cur)
        seen |= set(orb)
        if len(set(orb)) == 1: continue
        v = np.zeros(30, dtype=complex)
        for h, c in enumerate(orb): v[c] += w ** (chi * h)
        vecs.append(v)
    return np.array(vecs).T
L1c = [Y for Y in P.nb[t] if P.dZ[Y] == 1]
V = {chi: np.stack([orbit_vecs(g, chi) for g in sev]) for chi in range(1, 7)}
U = {chi: V[chi] / np.sqrt(7) for chi in V}
e = np.zeros(30, dtype=complex); e[k0] = 1
dprof = Counter(int(P.dZ[c]) for c in nbF)
print('input distances of the 30 slots:', dict(sorted(dprof.items())), '; O1 cells', L1c, flush=True)
res = {b: Counter() for b in range(1, 7)}; wit = {}
lams = L1c if mode == 'main' else [None]
for lam1 in lams:
    if mode == 'main':
        Rcols = [ni[c] for c in nbF if D[c, lam1] == 1]
    else:
        Rcols = [ni[c] for c in nbF if P.dZ[c] <= 3]
    Rm = np.zeros((30, len(Rcols)), dtype=complex)
    for j, c in enumerate(Rcols): Rm[c, j] = 1
    for ai in range(part, len(sev), nparts):
        A = np.concatenate([V[1][ai], Rm], axis=1)
        Qa, Ra = np.linalg.qr(A)
        keep = np.abs(np.diag(Ra)) > 1e-9
        Qa = Qa[:, keep] if keep.any() else Qa[:, :0]
        r = e - Qa @ (Qa.conj().T @ e)
        for b in range(1, 7):
            U2 = U[b]
            Wb = U2 - np.einsum('ij,njk->nik', Qa, np.einsum('ji,njk->nik', Qa.conj(), U2))
            Gb = np.conj(np.transpose(Wb, (0, 2, 1))) @ Wb
            cb = np.einsum('nji,j->ni', np.conj(Wb), r)
            Gp = np.linalg.pinv(Gb, rcond=1e-9, hermitian=True)
            proj2 = np.real(np.einsum('ni,nij,nj->n', np.conj(cb), Gp, cb))
            res2 = np.real(np.vdot(r, r)) - proj2
            ok = res2 < 1e-9
            res[b][int(ok.sum())] += 0
            res[b]['served'] += int(ok.sum()); res[b]['tested'] += len(ok)
            if ok.any() and b not in wit: wit[b] = (lam1, ai, int(np.argmax(ok)))
for b in range(1, 7):
    print(f'a = 1, b = {b}: (twist1, twist2, O1 cell) tested {res[b]["tested"]}, chi served in {res[b]["served"]}; witness {wit.get(b)}', flush=True)
print('part', part, 'of', nparts, 'done  %.0fs' % (time.time() - t0), flush=True)
