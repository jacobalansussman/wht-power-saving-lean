# c14_105a.py -- second filter on the symmetric shapes: Theorem M's certificate N in the irreducible 105a
# (N = sum of the radical weight spaces in the size-28 orbit; profile f(h)=20, f(h1h2)=38, f(K)=49; E = N cap K^-1 N).
# Per-irreducible Theorem R inequality for H-symmetric programs (slack mode):
#   sum_{move orbits} <psi, Ind_T^H Wn_(T,b)> >= sum_{array orbits} <psi, Ind_S^H (N/E)|_S>
# with the same move-orbit types (T, b) as the hook screen. A shape survives if the SAME move-orbit multiset that balances
# the hook also satisfies these inequalities.
import numpy as np, itertools, time
from collections import Counter
from lib514 import *
import c7_symscreen2 as S7
t0 = time.time()
R = S7.R
B7, B27 = R.B7, R.B27
Z189 = sum(np.kron(R.r7(TRANS[w]), R.r27(TRANS[w])) for w in range(1, NV)); Z189 = (Z189 + Z189.T) / 2
w, V = np.linalg.eigh(Z189); B105 = V[:, np.abs(w + 21) < 1e-6]; assert B105.shape[1] == 105
cache = {}
def r105(g):
    k = key(g)
    if k not in cache: cache[k] = B105.T @ np.kron(R.r7(g), R.r27(g)) @ B105
    return cache[k]
print('trace check 105a: id', round(np.trace(r105(ID))), ' transvection', round(np.trace(r105(TRANS[9]))), '(table 105, -35)')
A = [r105(a) for a in A_GENS]
rng = np.random.default_rng(4); Mix = sum(rng.normal() * a for a in A)
ww, VV = np.linalg.eigh((Mix + Mix.T) / 2)
wts = [tuple(int(np.round(VV[:, i] @ a @ VV[:, i])) for a in A) for i in range(105)]
# Levi orbits on weights (by transporting weight vectors)
wset = sorted(set(wts)); widx = {x: i for i, x in enumerate(wset)}
par = list(range(len(wset)))
def find(i):
    while par[i] != i: par[i] = par[par[i]]; i = par[i]
    return i
Lg = [levi([[1, 1, 0], [0, 1, 0], [0, 0, 1]]), levi([[0, 0, 1], [1, 0, 0], [0, 1, 0]])]
for l in Lg:
    rl = r105(l)
    for i in range(105):
        v = rl @ VV[:, i]; wv = tuple(int(np.round((v @ (a @ v)) / (v @ v))) for a in A)
        par[find(widx[wts[i]])] = find(widx[wv])
orb = Counter(find(widx[x]) for x in wset)
print('Levi orbits on the distinct radical weights of 105a:', sorted(orb.values()))
o28 = [r for r, c in orb.items() if c == 28][0]
N = VV[:, [i for i in range(105) if find(widx[wts[i]]) == o28]]
fN = lambda g: rank(np.hstack([N, r105(g) @ N])) - N.shape[1]
print('dim N', N.shape[1], ' profile', fN(h1), fN(comp(h1, h2)), fN(K), '(Theorem M: 20, 38, 49)')
# E = N cap K^-1 N
M = np.hstack([N, -np.linalg.inv(r105(K)) @ N]); u, s, vt = np.linalg.svd(M)
nul = vt[int((s > 1e-8).sum()):].T; E = N @ nul[:N.shape[1]]
Eb = np.linalg.qr(E)[0] if E.shape[1] else E
print('dim E =', Eb.shape[1])
d = h1; dinv = inv(d)
PN = N @ N.T
uu, ss, _ = np.linalg.svd((np.eye(105) - PN) @ r105(d) @ N); Wn = uu[:, :int((ss > 1e-8).sum())]
QNW = np.hstack([N, Wn]); print('cap module dim', Wn.shape[1])
def capN_trace(x):
    Mx = r105(x); C = np.linalg.lstsq(QNW, Mx @ QNW, rcond=None)[0]
    assert np.abs(QNW @ C - Mx @ QNW).max() < 1e-7
    return np.trace(C[N.shape[1]:, N.shape[1]:])
NK = r105(K) @ N; EK = r105(K) @ Eb if Eb.shape[1] else Eb
def demN_trace(x):
    """demand character on the c-side: trace on KN/KE for x in L (x keeps KN and KE)"""
    Mx = r105(x); C = np.linalg.lstsq(NK, Mx @ NK, rcond=None)[0]
    assert np.abs(NK @ C - Mx @ NK).max() < 1e-7
    t = np.trace(C)
    if EK.shape[1]:
        CE = np.linalg.lstsq(EK, Mx @ EK, rcond=None)[0]; t -= np.trace(CE)
    return t
P0 = S7.P0; Ekeys = S7.Ekeys
def analyse105(H):
    n = len(H); cls, rep, sizes, chars, clsof, idx = S7.char_table(H)
    ip = lambda a, b: np.real(np.sum(sizes * a * np.conj(b)) / n)
    Tcls = S7.subgroups_upto_conj(H)
    arr = []
    for S in Tcls:
        ind = S7.induced(H, rep, S, [demN_trace(x) for x in S])
        arr.append(np.round([ip(ind, c) for c in chars], 6))
    mov = []
    for T in Tcls:
        caps = {}; seen = set()
        for b in P0:
            bi = inv(b); xs = [comp(comp(bi, t), b) for t in T]
            kx = tuple(key(x) for x in xs)
            if kx in seen: continue
            seen.add(kx)
            if not all(key(x) in Ekeys for x in xs): continue
            ch_hook = tuple(np.round([S7.cap_trace(x) for x in xs], 6))
            ch_105 = tuple(np.round([capN_trace(x) for x in xs], 6))
            tg = [comp(comp(dinv, x), d) for x in xs]
            key2 = (ch_hook, ch_105)
            if not caps.get(key2, False):
                caps[key2] = any(all(np.array_equal(comp(comp(bp, x), inv(bp)), y) for x, y in zip(xs, tg)) for bp in P0)
        for (chh, ch5), ex in caps.items():
            indh = S7.induced(H, rep, T, list(chh)); ind5 = S7.induced(H, rep, T, list(ch5))
            mov.append(dict(T=len(T), size=n // len(T), hook=np.round([ip(indh, c) for c in chars], 6),
                            m105=np.round([ip(ind5, c) for c in chars], 6), exact=ex))
    arr_h = []
    for S in Tcls:
        ind = S7.induced(H, rep, S, [S7.yperp_trace(x) for x in S])
        arr_h.append(np.round([ip(ind, c) for c in chars], 6))
    return dict(n=n, Tsizes=[len(S) for S in Tcls], arr_h=arr_h, arr5=arr, mov=mov, k=len(chars))
def screen_both(info, W=5, Mv=14):
    n = info['n']; mov = [m for m in info['mov'] if m['exact']]; k = info['k']
    arrsz = [n // s for s in info['Tsizes']]
    results = []
    def arr_rec(i, rem, acch, acc5, cnt, norb):
        if rem < 0: return
        if i == len(arrsz):
            if rem == 0: yield acch, acc5, cnt, norb
            return
        c = 0
        while c * arrsz[i] <= rem:
            yield from arr_rec(i + 1, rem - c * arrsz[i], acch + c * info['arr_h'][i], acc5 + c * info['arr5'][i], cnt + [c], norb + c)
            c += 1
    for dh, d5, acnt, aorb in arr_rec(0, W, np.zeros(k), np.zeros(k), [], 0):
        sols = []
        def rec(i, rem, ah, a5, cnt, norb):
            if rem < 0 or np.any(ah - dh > 1e-6): return
            if i == len(mov):
                if rem == 0 and norb >= aorb + 2 and np.all(np.abs(ah - dh) < 1e-6):
                    sols.append((tuple(cnt), bool(np.all(a5 - d5 > -1e-6))))
                return
            c = 0
            while c * mov[i]['size'] <= rem:
                rec(i + 1, rem - c * mov[i]['size'], ah + c * mov[i]['hook'], a5 + c * mov[i]['m105'], cnt + [c], norb + c)
                c += 1
        rec(0, Mv, np.zeros(k), np.zeros(k), [], 0)
        if sols: results.append((tuple(acnt), sols))
    return results
if __name__ == '__main__':
    Lsubs = S7.subgroups_upto_conj(LEVI)
    for H in Lsubs:
        if len(H) not in (24, 12, 8, 6): continue
        info = analyse105(H)
        res = screen_both(info)
        nh = sum(len(s) for _, s in res); n5 = sum(sum(1 for x in s if x[1]) for _, s in res)
        print(f'H order {len(H)}: hook-balanced exact shapes {nh}; of these also passing the 105a inequalities: {n5}  [{time.time()-t0:.0f}s]', flush=True)
