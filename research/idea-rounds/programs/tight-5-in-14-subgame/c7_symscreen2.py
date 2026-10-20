# c7_symscreen2.py -- symmetry screen, general array orbits (H <= L, exact symmetry in t2-block3's sense).
# For a W-in-M program symmetric under H (global c: H -> L, arrays in orbits H/S, moves in orbits H/T):
#   (1) per-irreducible Theorem R on V_56: for every irreducible psi of H
#         sum_{move orbits} <psi, Ind_T^H W_(T,b)>  >=  sum_{array orbits} <psi, Ind_S^H Yperp|_S>
#       with equality in every psi when the program is tight (10 M = 28 W);
#       W_(T,b) = cap module (Y + dY)/Y with T acting through b^-1 c(t) b in E (b in P0), Yperp = V56/Y (via c);
#   (2) S5' (t2-block3, upheld): #move orbits >= #array orbits + 2;   (3) sum of move-orbit sizes = M.
# A move-orbit type is EXACT if the letter can be chosen to commute with c(T) (the t2-block3 definition);
# non-exact types are admitted only in the 'gauge' column (cap module still T-stable; see NOTE.md).
import numpy as np, itertools, time, sys
from collections import Counter
from lib514 import *
t0 = time.time()
R = Irreps()
dat = np.load('v56_basis.npz'); Vw = dat['V']; Y = Vw[:, dat['o28']]; Yp = Vw[:, dat['o7']]
d = h1
P0 = gen_closure(A_GENS + [levi([[1, 1, 0], [0, 1, 0], [0, 0, 1]]), levi([[0, 0, 1], [1, 0, 0], [0, 1, 0]])])
dLZ = frozenset(int(d[v]) for v in LZ)
E = [b for b in P0 if frozenset(int(b[v]) for v in dLZ) == dLZ]
Ekeys = {key(b) for b in E}
PY = Y @ Y.T
u, s, vt = np.linalg.svd((np.eye(56) - PY) @ R.r56(d) @ Y); Wb = u[:, :10]
QYW = np.hstack([Y, Wb])
def cap_trace(x):
    M = R.r56(x); C = np.linalg.lstsq(QYW, M @ QYW, rcond=None)[0]
    assert np.abs(QYW @ C - M @ QYW).max() < 1e-8
    return np.trace(C[28:, 28:])
def yperp_trace(x): return np.trace(Yp.T @ R.r56(x) @ Yp)
dinv = inv(d)

def subgroups_upto_conj(H):
    n = len(H); subs = {}
    for a, b in itertools.combinations_with_replacement(range(n), 2):
        S = gen_closure([H[a], H[b]]); subs[frozenset(key(x) for x in S)] = S
    out = []; seen = set()
    for k, S in sorted(subs.items(), key=lambda kv: len(kv[1])):
        if k in seen: continue
        conj = set()
        for g in H:
            gi = inv(g); conj.add(frozenset(key(comp(comp(g, x), gi)) for x in S))
        seen |= conj; out.append(S)
    return out

def char_table(S):
    n = len(S); ks = [key(x) for x in S]; idx = {kk: i for i, kk in enumerate(ks)}
    cls = []; done = set()
    for x in S:
        if key(x) in done: continue
        c = {key(comp(comp(g, x), inv(g))) for g in S}; done |= c; cls.append(sorted(idx[k] for k in c))
    k = len(cls); clsof = np.zeros(n, dtype=int)
    for ci, c in enumerate(cls):
        for i in c: clsof[i] = ci
    Mats = []
    for i in range(k):
        A = np.zeros((k, k))
        for j in range(k):
            cnt = np.zeros(k)
            for xi in cls[i]:
                for xj in cls[j]: cnt[clsof[idx[key(comp(S[xi], S[xj]))]]] += 1
            for l in range(k): A[j, l] = cnt[l] / len(cls[l])
        Mats.append(A)
    rng = np.random.default_rng(0); Mix = sum(rng.normal() * A for A in Mats)
    w, V = np.linalg.eig(Mix)
    sizes = np.array([len(c) for c in cls]); idc = [ci for ci, c in enumerate(cls) if np.array_equal(S[c[0]], ID)][0]
    chars = []
    for j in range(k):
        om = V[:, j] / V[idc, j]; co = om / sizes
        d1 = np.sqrt((n / np.sum(sizes * np.abs(co) ** 2)).real); chars.append(co * d1)
    rep = [S[c[0]] for c in cls]
    # sanity: orthonormality
    G = np.array([[np.sum(sizes * a * np.conj(b)) / n for b in chars] for a in chars])
    assert np.allclose(G, np.eye(k), atol=1e-6), G
    return cls, rep, sizes, chars, clsof, idx

def induced(H, rep, T, chvals):
    """Ind_T^H of the class function chvals (list aligned with T) evaluated at class reps"""
    Tk = {key(t): i for i, t in enumerate(T)}; out = []
    for g in rep:
        s = 0.0
        for y in H:
            kz = key(comp(comp(y, g), inv(y)))
            if kz in Tk: s += chvals[Tk[kz]]
        out.append(s / len(T))
    return np.array(out)

def analyse(H, label):
    n = len(H); cls, rep, sizes, chars, clsof, idx = char_table(H)
    ip = lambda a, b: np.real(np.sum(sizes * a * np.conj(b)) / n)
    Tcls = subgroups_upto_conj(H)
    # array-orbit types: stabiliser S (any subgroup class) -> demand multiplicities
    arr = []
    for S in Tcls:
        ch = [yperp_trace(x) for x in S]
        ind = induced(H, rep, S, ch)
        arr.append(dict(S=len(S), size=n // len(S), mult=np.round([ip(ind, c) for c in chars], 6)))
    # move-orbit types
    mov = []
    for T in Tcls:
        caps = {}; seenconj = set()
        for b in P0:
            bi = inv(b); xs = [comp(comp(bi, t), b) for t in T]
            kx = tuple(key(x) for x in xs)
            if kx in seenconj: continue
            seenconj.add(kx)
            if not all(key(x) in Ekeys for x in xs): continue
            ch = tuple(np.round([cap_trace(x) for x in xs], 6))
            tg = [comp(comp(dinv, x), d) for x in xs]
            ex = caps.get(ch, False)
            if not ex:
                ex = any(all(np.array_equal(comp(comp(bp, x), inv(bp)), y) for x, y in zip(xs, tg)) for bp in P0)
            caps[ch] = ex
        for ch, ex in caps.items():
            ind = induced(H, rep, T, list(ch))
            mov.append(dict(T=len(T), size=n // len(T), mult=np.round([ip(ind, c) for c in chars], 6), exact=ex))
    dims = [int(round(np.real(c[[i for i, cc in enumerate(cls) if np.array_equal(H[cc[0]], ID)][0]]))) for c in chars]
    return dict(n=n, dims=dims, arr=arr, mov=mov)

def solve(info, W, M, tight, gauge, max_sols=50):
    """enumerate array-orbit multisets (sizes summing to W) and move-orbit multisets (sizes summing to M)"""
    arr, mov = info['arr'], [m for m in info['mov'] if gauge or m['exact']]
    out = []
    def arr_rec(i, rem, acc, cnt, norb):
        if rem < 0: return
        if i == len(arr):
            if rem == 0: yield acc, cnt, norb
            return
        c = 0
        while c * arr[i]['size'] <= rem:
            yield from arr_rec(i + 1, rem - c * arr[i]['size'], acc + c * arr[i]['mult'], cnt + [c], norb + c)
            c += 1
    for dem, acnt, aorb in arr_rec(0, W, np.zeros(len(info['dims'])), [], 0):
        sols = []
        def mov_rec(i, rem, acc, cnt, norb):
            if len(sols) >= max_sols: return
            if rem < 0: return
            if tight and np.any(acc - dem > 1e-6): return
            if i == len(mov):
                if rem == 0 and norb >= aorb + 2:
                    ok = np.all(np.abs(acc - dem) < 1e-6) if tight else np.all(acc - dem > -1e-6)
                    if ok: sols.append((tuple(cnt), norb))
                return
            c = 0
            while c * mov[i]['size'] <= rem:
                mov_rec(i + 1, rem - c * mov[i]['size'], acc + c * mov[i]['mult'], cnt + [c], norb + c)
                c += 1
        mov_rec(0, M, np.zeros(len(info['dims'])), [], 0)
        if sols: out.append((tuple(acnt), aorb, sols))
    return out

def describe(info, W, M, tight, gauge):
    res = solve(info, W, M, tight, gauge)
    arr = info['arr']; mov = [m for m in info['mov'] if gauge or m['exact']]
    lines = []
    for acnt, aorb, sols in res:
        a = '+'.join(f"{c}x{arr[i]['size']}" for i, c in enumerate(acnt) if c)
        for cnt, norb in sols[:3]:
            m = '+'.join(f"{c}x{mov[i]['size']}(T{mov[i]['T']}{'' if mov[i]['exact'] else '*'})" for i, c in enumerate(cnt) if c)
            lines.append(f"arrays {a} | moves {m}")
    return res, lines

if __name__ == '__main__':
    Lsubs = subgroups_upto_conj(LEVI)
    print('subgroup classes of L (orders):', [len(S) for S in Lsubs], '[%.0fs]' % (time.time() - t0))
    S3c = gen_closure([levi([[0, 1, 0], [1, 0, 0], [0, 0, 1]]), levi([[0, 0, 1], [1, 0, 0], [0, 1, 0]])])
    info = analyse(S3c, 'S3')
    print('S3 (coordinate permutations): irreps', info['dims'])
    for a in info['arr']: print('   array orbit |S|=%d size %d demand %s' % (a['S'], a['size'], a['mult']))
    for m in info['mov']: print('   move orbit |T|=%d size %d cap %s exact=%s' % (m['T'], m['size'], m['mult'], m['exact']))
    print('CALIBRATION (slack mode = inequality):')
    for (W, M) in [(6, 18), (5, 15), (2, 5), (5, 13), (3, 9)]:
        res, lines = describe(info, W, M, tight=False, gauge=False)
        print(f'   S3 exact, {W} in {M}: {len(res)} array structures pass; e.g.', lines[:2])
    print('S3, 5 in 14 tight, EXACT types: ')
    res, lines = describe(info, 5, 14, tight=True, gauge=False)
    for l in lines: print('     ', l)
    print('S3, 5 in 14 tight, GAUGE types allowed: ')
    res, lines = describe(info, 5, 14, tight=True, gauge=True)
    for l in lines[:12]: print('     ', l)
    print('[%.0fs]' % (time.time() - t0), flush=True)
    summary = []
    for H in Lsubs:
        if len(H) == 1: continue
        info = analyse(H, str(len(H)))
        rex, lex = describe(info, 5, 14, tight=True, gauge=False)
        rga, lga = describe(info, 5, 14, tight=True, gauge=True)
        eo = sorted(set(order(x) for x in H))
        nshape_ex = sum(len(s) for _, _, s in rex); nshape_ga = sum(len(s) for _, _, s in rga)
        print(f'H order {len(H)} elt orders {eo}: exact shapes {nshape_ex} (array structures {len(rex)}), gauge shapes {nshape_ga} (array structures {len(rga)})', flush=True)
        for l in lex[:4]: print('      exact:', l)
        for l in lga[:3]: print('      gauge:', l)
        summary.append((len(H), eo, nshape_ex, nshape_ga))
    print('SUMMARY (order, element orders, #exact shapes, #gauge shapes):', summary)
    print('[%.0fs]' % (time.time() - t0))
