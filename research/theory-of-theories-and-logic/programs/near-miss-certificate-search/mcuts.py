# mcuts.py -- cert-07 (a'): the hook's own module M = Hom_Pauli(V, V (x) Lam^8 V), V = C^16 (block 4), basis f_T
# (T an 8-subset of F_2^4, f_T(e_0) = e_sigma(T) (x) e_T, sigma(T) = sum of T); X = V (x) M, and every P_0-invariant
# subspace of X is V (x) N_M.  P_0 acts monomially on M (A diagonally, GL(4,2) by T -> AT), so unions of GL(4,2)-orbits
# of 8-subsets are P_0-invariant ("coordinate cuts" of M); the hook is N_M = {f_T : sigma(T) in T}.
# Depth-j element C = h^(x)j (x) I (qubits 1..j = address bits 0..j-1):
#   (C.f_T)(e_0) = C f_T(C^-1 e_0) = 2^(-j/2) sum_{a in W_j} (C X^a (x) Lam^8(C X^a)) (e_sigma(T) (x) e_T), so
#   G_j[T', T] = 2^(-j/2) sum_a (C X^a)[sigma(T'), sigma(T)] * Lam^8(C X^a)[T', T].
# C X^a is block diagonal over the 2^(4-j) contiguous address blocks of size 2^j, so Lam^8(C X^a)[T',T] =
# prod_z det(block[T'_z, T_z]) and G_j is block diagonal by the composition (|T_z|)_z.  psi_j(cut S) = sum over
# composition classes of rank G_j[class - S, class cap S].   usage: mcuts.py validate | mcuts.py all OUT.json [part k of K]
import itertools, sys, json, time
import numpy as np
h = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
SUBS = list(itertools.combinations(range(16), 8)); IDX = {s: i for i, s in enumerate(SUBS)}
NT = len(SUBS)
SIG = np.array([np.bitwise_xor.reduce(np.array(s)) for s in SUBS])
def orbits():
    def matact(M, x): return sum((((sum(M[i][j] * ((x >> j) & 1) for j in range(4))) & 1) << i) for i in range(4))
    E12 = [[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]; CYC = np.roll(np.eye(4, dtype=int), 1, axis=0).tolist()
    P = [np.array([matact(M, x) for x in range(16)]) for M in (E12, CYC)]
    S = np.array(SUBS); parent = list(range(NT))
    def find(a):
        while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for p in P:
        im = np.sort(p[S], axis=1)
        for i, r in enumerate(im):
            a, b = find(i), find(IDX[tuple(r)])
            if a != b: parent[a] = b
    roots = {}
    for i in range(NT): roots.setdefault(find(i), []).append(i)
    orb = sorted(roots.values(), key=lambda v: (len(v), min(v)))
    lab = np.empty(NT, dtype=int)
    for k, v in enumerate(orb): lab[v] = k
    return orb, lab
def blocks_for(j):
    n = 2 ** j; Z = 2 ** (4 - j)
    Hj = np.array([[1.0]])
    for _ in range(j): Hj = np.kron(Hj, h)
    Ba = []
    for a in range(n):
        Xa = np.zeros((n, n)); Xa[np.arange(n) ^ a, np.arange(n)] = 1
        Ba.append(Hj @ Xa)
    # minor tables Mn[a][c', c] for subsets (bitmasks) of an n-set
    masks = list(range(1 << n)); pc = [bin(c).count('1') for c in masks]
    Mn = np.zeros((n, 1 << n, 1 << n))
    for a in range(n):
        for c1 in masks:
            r = [i for i in range(n) if c1 >> i & 1]
            for c2 in masks:
                if pc[c2] != pc[c1]: continue
                cc = [i for i in range(n) if c2 >> i & 1]
                Mn[a, c1, c2] = np.linalg.det(Ba[a][np.ix_(r, cc)]) if r else 1.0
    S = np.array(SUBS)
    codes = np.zeros((NT, Z), dtype=np.int64)
    for z in range(Z):
        inb = (S >> j) == z
        codes[:, z] = ((inb * (1 << (S & (n - 1)))).sum(axis=1))
    comp = [tuple(bin(int(c)).count('1') for c in row) for row in codes]
    classes = {}
    for i, k in enumerate(comp): classes.setdefault(k, []).append(i)
    out = {}
    for k, rows in classes.items():
        rows = np.array(rows); c = codes[rows]
        sl = SIG[rows] & (n - 1)
        G = np.zeros((len(rows), len(rows)))
        for a in range(n):
            T = Ba[a][sl[:, None], sl[None, :]]
            for z in range(Z): T = T * Mn[a][c[:, z][:, None], c[:, z][None, :]]
            G += T
        G /= np.sqrt(n)
        out[k] = (rows, G)
    return out
def psi_of(blocks, inS, tol=1e-8, want_phi=True):
    psi = 0; phi = 0.0
    for k, (rows, G) in blocks.items():
        m = inS[rows]
        if m.all() or (~m).all(): continue
        M = G[np.ix_(~m, m)]
        s = np.linalg.svd(M, compute_uv=False)
        psi += int((s > tol).sum()); phi += float((s ** 2).sum())
    return psi, phi
if __name__ == '__main__':
    t0 = time.time()
    orb, lab = orbits()
    print("orbits", [len(o) for o in orb], flush=True)
    BL = {j: blocks_for(j) for j in (1, 2, 3)}
    for j in (1, 2, 3):
        err = max(np.abs(G @ G.T - np.eye(len(r))).max() for r, G in BL[j].values())
        print(f"j={j}: {len(BL[j])} composition classes, max size {max(len(r) for r, G in BL[j].values())}, orth err {err:.1e} ({time.time()-t0:.0f}s)", flush=True)
    hook = np.array([SIG[i] in SUBS[i] for i in range(NT)])
    hk = sorted(set(lab[hook])); print("hook orbits", hk, "dim", hook.sum(), flush=True)
    if sys.argv[1] == 'validate':
        for j in (1, 2, 3):
            print("hook psi_%d, phi_%d =" % (j, j), psi_of(BL[j], hook), f"({time.time()-t0:.0f}s)", flush=True)
        sys.exit()
    # all cuts: unions of orbits not containing the last orbit (complement pairs give equal psi)
    K = len(orb); part, nparts = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (0, 1)
    res = []
    cuts = [b for b in range(1, 1 << (K - 1))]
    cuts = cuts[part::nparts]
    for b in cuts:
        inS = np.array([(b >> int(lab[i])) & 1 == 1 for i in range(NT)])
        p1, f1 = psi_of(BL[1], inS)
        # cheap bound: psi3 <= sum over j=3 classes of min(|S cap c|, |c - S|)
        ub = sum(min(int(inS[r].sum()), len(r) - int(inS[r].sum())) for r, G in BL[3].values())
        rec = dict(bits=b, orbits=[k for k in range(K) if b >> k & 1], dimN=int(inS.sum()), psi1=p1, phi1=f1, ub3=ub)
        if p1 > 0 and ub / p1 > 2.3:
            p3, f3 = psi_of(BL[3], inS); rec.update(psi3=p3, phi3=f3, ratio3=p3 / p1, spec3=f3 / f1 if f1 else None)
        res.append(rec)
        if len(res) % 20 == 0:
            json.dump(res, open(sys.argv[2], "w")); print(len(res), "cuts", f"({time.time()-t0:.0f}s)", flush=True)
    json.dump(res, open(sys.argv[2], "w"))
    best = sorted([r for r in res if 'ratio3' in r], key=lambda r: -r['ratio3'])[:10]
    for r in best: print(r, flush=True)
    print(f"done {len(res)} cuts ({time.time()-t0:.0f}s)")
