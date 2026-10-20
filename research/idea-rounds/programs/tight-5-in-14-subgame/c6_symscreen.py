# c6_symscreen.py -- which symmetry groups can a TIGHT W-in-M block-3 unit carry? (exact symmetry, t2-block3's
# definition: tau_h(g e_l) = c_h g c'_h e_{h.l}, one global c: H -> L for programs with a fixed array.)
# Necessary conditions used (all from V_56, own construction):
#  (1) per-irreducible hook balance (Theorem R restricted to rho-isotypic inputs, tight => equality in EVERY rho):
#      sum over move orbits O of <rho, Ind_{T_O}^H W_O> = sum over array orbits of <rho, Ind_S^H Yperp|_S>
#      where W = (Y + dY)/Y (10-dim cap module) with T acting through x_t = b^-1 c_t b in E = B cap dBd^-1, and the
#      demand of an array with stabiliser S is Yperp = V56/Y ~ KY with S acting through c (derivation in NOTE.md);
#  (2) t2-block3 S5' (upheld): number of move orbits >= number of array orbits + 2;
#  (3) sum of orbit sizes = M.
# Calibration: Z7 (Singer) all arrays fixed must fail 5-in-14 by (2) only; S3 must pass plain 6 in 18 and 5 in 15
# and fail 2 in 5 (and 5 in 13).
import numpy as np, itertools, time, sys
from fractions import Fraction
from collections import defaultdict
from lib514 import *
t0 = time.time()
R = Irreps()
dat = np.load('v56_basis.npz'); Vw = dat['V']; Y = Vw[:, dat['o28']]; Yp = Vw[:, dat['o7']]
d = h1
P0 = gen_closure(A_GENS + [levi([[1, 1, 0], [0, 1, 0], [0, 0, 1]]), levi([[0, 0, 1], [1, 0, 0], [0, 1, 0]])])
dLZ = frozenset(int(d[v]) for v in LZ)
E = [b for b in P0 if frozenset(int(b[v]) for v in dLZ) == dLZ]
print('|P0| =', len(P0), ' |E| =', len(E))
Ekeys = {key(b) for b in E}
# cap module basis: columns of r(d) Y orthogonal to Y
PY = Y @ Y.T
Wraw = (np.eye(56) - PY) @ R.r56(d) @ Y
u, s, vt = np.linalg.svd(Wraw); Wb = u[:, :int((s > 1e-8).sum())]
assert Wb.shape[1] == 10
QYW = np.hstack([Y, Wb])                       # basis of Y + dY (orthonormal: Y and Wb orthogonal)
def cap_trace(x):
    """trace of r(x) on (Y + dY)/Y for x in E (x keeps Y and Y + dY)"""
    M = R.r56(x)
    # matrix of M on Y+dY in basis QYW, then the quotient block
    C = np.linalg.lstsq(QYW, M @ QYW, rcond=None)[0]
    assert np.abs(QYW @ C - M @ QYW).max() < 1e-8
    return np.trace(C[28:, 28:])
def yperp_trace(x): return np.trace(Yp.T @ R.r56(x) @ Yp)
# ---------------- subgroups of L up to conjugacy
Lkey = {key(l): l for l in LEVI}
def subgroup(gens): return gen_closure(gens)
subs = {}
for a, b in itertools.combinations_with_replacement(range(168), 2):
    S = subgroup([LEVI[a], LEVI[b]]); k = frozenset(key(x) for x in S)
    subs[k] = S
for a in range(168):
    S = subgroup([LEVI[a]]); subs[frozenset(key(x) for x in S)] = S
classes = []
seen = set()
for k, S in sorted(subs.items(), key=lambda kv: len(kv[1])):
    if k in seen: continue
    conj = set()
    for l in LEVI:
        li = inv(l); conj.add(frozenset(key(comp(comp(l, x), li)) for x in S))
    seen |= conj
    classes.append(S)
print('subgroup classes of L:', [len(S) for S in classes], '[%.0fs]' % (time.time() - t0))
def name(S):
    n = len(S); orders = sorted(order(x) for x in S)
    from collections import Counter
    oc = Counter(orders)
    fixed_hyp = sum(1 for U in range(1, 8) if all(True for _ in [0]))  # placeholder
    return f"order {n} elt-orders {dict(oc)}"
# ---------------- character table of a small permutation-represented group (by class functions)
def conj_classes(S):
    ks = [key(x) for x in S]; idx = {k: i for i, k in enumerate(ks)}
    cls = []; done = set()
    for x in S:
        if key(x) in done: continue
        c = {key(comp(comp(g, x), inv(g))) for g in S}; done |= c; cls.append(sorted(idx[k] for k in c))
    return cls
def char_table(S):
    """irreducible characters via Burnside/Dixon-free method: eigenvectors of class-multiplication matrices (numerical)"""
    cls = conj_classes(S); n = len(S); k = len(cls)
    ks = [key(x) for x in S]; idx = {kk: i for i, kk in enumerate(ks)}
    clsof = np.zeros(n, dtype=int)
    for ci, c in enumerate(cls):
        for i in c: clsof[i] = ci
    # class multiplication coefficients a_{ij}^l: (C_i C_j) = sum a_ijl C_l
    Mats = []
    rep = [S[c[0]] for c in cls]
    for i in range(k):
        A = np.zeros((k, k))
        for j in range(k):
            cnt = np.zeros(k)
            for xi in cls[i]:
                for xj in cls[j]:
                    cnt[clsof[idx[key(comp(S[xi], S[xj]))]]] += 1
            # number of times each element of C_l appears
            for l in range(k): A[j, l] = cnt[l] / len(cls[l])
        Mats.append(A)
    rng = np.random.default_rng(0)
    Mix = sum(rng.normal() * A for A in Mats)
    w, V = np.linalg.eig(Mix)
    chars = []
    sizes = np.array([len(c) for c in cls])
    for j in range(k):
        v = V[:, j]                     # v_l proportional to |C_l| chi(g_l)/chi(1)
        v = v / v[[ci for ci, c in enumerate(cls) if any(np.array_equal(S[t], ID) for t in c)][0]]
        omega = v                        # central character omega_l = |C_l| chi_l / chi(1)
        chi_over = omega / sizes        # chi_l/chi(1)
        # chi(1)^2 = n / sum_l |C_l| |chi_l/chi1|^2
        d2 = n / np.sum(sizes * np.abs(chi_over) ** 2)
        d1 = np.sqrt(d2.real); chars.append(np.round(chi_over * d1, 8))
    return cls, rep, sizes, chars
def inner(sizes, a, b, n): return np.sum(sizes * a * np.conj(b)) / n
# ---------------- the screen
def screen(H, W, M, label, max_orbits=None, verbose=True, exact_only=False):
    n = len(H); cls, rep, sizes, chars = char_table(H)
    Hk = {key(x): i for i, x in enumerate(H)}
    clsof = {}
    for ci, c in enumerate(cls):
        for i in c: clsof[key(H[i])] = ci
    # demand: all arrays fixed (orbit structure 1^W) -> W * Yperp|_H
    dem = np.array([yperp_trace(r) for r in rep])
    dmult = [W * inner(sizes, dem, ch, n) for ch in chars]
    # subgroups T of H up to H-conjugacy
    Hsubs = {}
    for a, b in itertools.combinations_with_replacement(range(n), 2):
        S = gen_closure([H[a], H[b]]); Hsubs[frozenset(key(x) for x in S)] = S
    Tcls = []; seenT = set()
    for kk, S in sorted(Hsubs.items(), key=lambda kv: len(kv[1])):
        if kk in seenT: continue
        conj = set()
        for g in H:
            gi = inv(g); conj.add(frozenset(key(comp(comp(g, x), gi)) for x in S))
        seenT |= conj; Tcls.append(S)
    # orbit types: (T, cap character class) for every b in P0 with b^-1 T b in E
    types = []
    for T in Tcls:
        caps = {}; seenconj = set()
        for b in P0:
            bi = inv(b)
            xs = [comp(comp(bi, t), b) for t in T]
            kx = tuple(key(x) for x in xs)
            if kx in seenconj: continue
            seenconj.add(kx)
            if all(key(x) in Ekeys for x in xs):
                ch = tuple(np.round([cap_trace(x) for x in xs], 6))
                if ch not in caps:
                    # exactness: exists b' with b' x b'^-1 = d^-1 x d for all x (move letter d b' commutes with x)
                    ex = None
                    if exact_only or verbose:
                        di = inv(d)
                        tg = [comp(comp(di, x), d) for x in xs]
                        ex = any(all(np.array_equal(comp(comp(bp, x), inv(bp)), y) for x, y in zip(xs, tg)) for bp in P0)
                    caps[ch] = (xs, ex)
        for ch, (xs, ex) in caps.items():
            if exact_only and not ex: continue
            # induced character Ind_T^H of the class function t -> ch on T: value at class rep g:
            # (1/|T|) sum_{y in H} [y g y^-1 in T] ch(y g y^-1)
            Tk = {key(t): i for i, t in enumerate(T)}
            indv = []
            for g in rep:
                s = 0.0
                for y in H:
                    z = comp(comp(y, g), inv(y)); kz = key(z)
                    if kz in Tk: s += ch[Tk[kz]]
                indv.append(s / len(T))
            indv = np.array(indv)
            mult = [inner(sizes, indv, chi, n) for chi in chars]
            types.append(dict(T=len(T), size=n // len(T), mult=np.round(np.real(mult), 6), exact=ex))
    if verbose:
        print(f'  [{label}] |H|={n}, irreps dims {[int(round(c[[i for i,c0 in enumerate(cls) if any(np.array_equal(H[t], ID) for t in c0)][0]].real)) for c in chars]}')
        print(f'    demand multiplicities (W={W}):', np.round(np.real(dmult), 4))
        for t in types: print(f'    move-orbit type: |T|={t["T"]} orbit size {t["size"]} cap mults {t["mult"]} exact={t["exact"]}')
    # integer search: choose counts of orbit types
    target = np.round(np.real(dmult), 6)
    sols = []
    nt = len(types)
    def rec(i, remaining_moves, acc, counts, norb):
        if remaining_moves < 0: return
        if np.any(acc - target > 1e-6): return
        if i == nt:
            if remaining_moves == 0 and np.all(np.abs(acc - target) < 1e-6): sols.append((tuple(counts), norb))
            return
        t = types[i]; c = 0
        while c * t['size'] <= remaining_moves:
            rec(i + 1, remaining_moves - c * t['size'], acc + c * t['mult'], counts + [c], norb + c)
            c += 1
    rec(0, M, np.zeros(len(target)), [], 0)
    good = [s for s in sols if s[1] >= W + 2]
    if verbose:
        print(f'    shapes balancing every irreducible with {M} moves: {len(sols)}; with >= {W+2} move orbits (S5\'): {len(good)}')
        for s in good[:6]: print('      counts', s[0], 'orbits', s[1])
    return sols, good, types

if __name__ == '__main__':
    # calibration: Singer Z7 (all arrays fixed) and S3 (coordinate permutations)
    Z7 = [S for S in classes if len(S) == 7][0]
    print('--- calibration Z7, all arrays fixed, 5 in 14 (t2-block3: excluded by the orbit count only)')
    screen(Z7, 5, 14, 'Z7 5in14')
    S3 = [S for S in classes if len(S) == 6]
    S3c = gen_closure([levi([[0, 1, 0], [1, 0, 0], [0, 0, 1]]), levi([[0, 0, 1], [1, 0, 0], [0, 1, 0]])])
    print('--- calibration S3 (coordinate permutations): plain 6 in 18 must pass, 2 in 5 must fail')
    s, g, _ = screen(S3c, 6, 18, 'S3 6in18', verbose=False); print('   S3 6 in 18: shapes', len(s), 'with orbit count', len(g))
    s, g, _ = screen(S3c, 2, 5, 'S3 2in5', verbose=False); print('   S3 2 in 5: shapes', len(s), 'with orbit count', len(g))
    s, g, _ = screen(S3c, 5, 15, 'S3 5in15', verbose=False); print('   S3 5 in 15: shapes', len(s), 'with orbit count', len(g))
    s, g, _ = screen(S3c, 5, 13, 'S3 5in13', verbose=False); print('   S3 5 in 13: shapes', len(s), 'with orbit count', len(g))
    print('[%.0fs]' % (time.time() - t0))
    print('=== 5 in 14, all arrays fixed, every subgroup class of L')
    for S in classes:
        if len(S) == 1: continue
        s, g, ty = screen(S, 5, 14, f'|H|={len(S)}', verbose=(len(S) >= 6))
        print(f'  H of order {len(S)} (element orders {sorted(set(order(x) for x in S))}): balanced shapes {len(s)}, with orbit count {len(g)}', flush=True)
    print('[%.0fs]' % (time.time() - t0))
