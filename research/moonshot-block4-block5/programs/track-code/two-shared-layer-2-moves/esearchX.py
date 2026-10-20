"""esearchX.py -- w5-twoshared: Lemma D' (the extension of Lemma D to SHARED layer 2) tested exhaustively over all
20,643,840 symplectic engine letters, on cross cores (both merges carry both arrays, F diagonal; sweep.core 'cross').
Lemma D' [proved in NOTE s4]: a layer-2 junk J = y_e q_kj of top e that lies in NO frame of a u main term (B'p_k'j)
can only be cancelled through eps_e v; then ONE letter e' in supp(eps_e) must satisfy the J-condition (s_j -> v, or
x_j -> u_l) AND, for EVERY layer-3 term q_k'j' of v, a production condition for e' q_k'j' inside t_e:
   s_j' -> v type      S(S_(q_k'j' (D b_j')^-1) L0) = S_(D mu_e D) L0
   x_j' -> u_l type    S(S_(q_k'j') L0)             = S_(D mu_e D g_l D) L0
   u_k'' -> t_e type   S(S_(D g_k' p_k'j' p_k''j'^-1) L0) = S_D L0     (k'' != k': a second merge's term on array j')
A junk that lies in a u frame is reported separately (not covered).  Usage: esearchX.py ncores seed"""
import sys, os, random, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from esearch import maps, basis, img_L0, GL, SYM
GLidx = {c: i for i, c in enumerate(GL)}
def symidx(Bm):
    t = 0; bits = 0
    for i in range(4):
        for j in range(i, 4):
            bits |= Bm[i][j] << t; t += 1
    return bits
def single(S):
    """boolean (nP, 1024) marking the one engine letter with symplectic S (if S stabilises L0)"""
    out = np.zeros((len(GL), 1024), dtype=bool)
    if sp.lag(S) != sp.L0: return out
    P = tuple(S[i] & 15 for i in range(4))
    Pm = {}
    for v in range(16):
        w = 0
        for j in range(4):
            if (v >> j) & 1: w ^= P[j]
        Pm[v] = w
    Pinv = {w: v for v, w in Pm.items()}
    # B P e_i = S[i] >> 4  ->  B e_j = (S[Pinv-col] >> 4)
    BPe = [S[i] >> 4 for i in range(4)]
    Bm = [[0] * 4 for _ in range(4)]
    for j in range(4):
        x = Pinv[1 << j]          # P x = e_j
        val = 0
        for i in range(4):
            if (x >> i) & 1: val ^= BPe[i]
        for r in range(4): Bm[r][j] = (val >> r) & 1
    assert sp.symp_from(list(P), Bm) == tuple(S), "single(): reconstruction failed"
    out[GLidx[P], symidx(Bm)] = True
    return out
import sp
from sweep import core
ncores, seed = int(sys.argv[1]), int(sys.argv[2])
rng = random.Random(seed)
D = sp.D
tot = {'excluded': 0, 'open': 0, 'ufr': 0}
for c in range(ncores):
    t0 = time.time()
    R0, R1, a11, a20, tg = core(rng, 'cross')
    b = {0: R0['b'][0], 1: R1['b'][0]}; s = {0: R0['s'], 1: R1['s']}
    g = {1: R0['g'][0], 2: R1['g'][0]}
    a = {(1, 0): R0['a'][0], (2, 1): R1['a'][0], (1, 1): a11[0], (2, 0): a20[0]}
    p = {kj: sp.mv(a[kj], s[kj[1]]) for kj in a}
    q = {kj: sp.mv(g[kj[0]], p[kj]) for kj in a}
    mu = {0: R0['mu'][0], 1: R1['mu'][0]}; route = {0: (1, 0), 1: (2, 1)}
    lay3 = [kj for kj in q if sp.lam(q[kj]) == 3]
    line = [f"core {c}: tags {tg}, layer-3 terms of v {sorted(lay3)}"]
    tops_ok = 0; ufr = False
    for e in (0, 1):
        y = sp.comp(D, mu[e]); X = sp.comp(D, sp.comp(mu[e], D))
        # production conditions for every layer-3 term of v
        allk = None
        for (k1, j1) in lay3:
            cnd = maps(basis(img_L0(sp.comp(q[(k1, j1)], sp.inv(s[j1])))), img_L0(X))
            for l in (1, 2):
                W_ = sp.comp(D, sp.comp(mu[e], sp.comp(D, sp.comp(g[l], D))))
                cnd |= maps(basis(img_L0(q[(k1, j1)])), img_L0(W_))
            for k2 in (1, 2):
                if k2 == k1 or (k2, j1) not in p: continue
                Y = sp.comp(D, sp.comp(g[k1], sp.comp(p[(k1, j1)], sp.inv(p[(k2, j1)]))))
                cnd |= maps(basis(img_L0(Y)), img_L0(D))
            for kj2 in q:
                if kj2[1] == j1 and sp.lam(q[kj2]) == 2:
                    cnd |= single(sp.comp(y, sp.comp(q[kj2], sp.inv(q[(k1, j1)]))))
            allk = cnd if allk is None else (allk & cnd)
        top_ok = True; info = []
        for kj in q:
            if kj == route[e]: continue
            J = sp.comp(y, q[kj]); L = sp.lam(J)
            if L != 2: info.append(f"{kj}:lam{L}"); continue
            if any(sp.lag(J) == sp.lag(p[(k2, kj[1])]) for k2 in (1, 2) if (k2, kj[1]) in p):
                info.append(f"{kj}:u-frame"); ufr = True; continue
            j = kj[1]
            cj = maps(basis(img_L0(D)), img_L0(sp.comp(J, sp.inv(s[j]))))
            for l in (1, 2):
                cj |= maps(basis(img_L0(sp.comp(D, sp.comp(g[l], D)))), img_L0(J))
            for kj2 in q:
                if kj2[1] == j and sp.lam(q[kj2]) == 2:
                    cj |= single(sp.comp(J, sp.inv(q[kj2])))
            n = int((cj & allk).sum())
            info.append(f"{kj}:{n}")
            if n == 0: top_ok = False
        tops_ok += top_ok
        line.append(f"top {e}: {'passes' if top_ok else 'EXCLUDED'} [{' '.join(info)}]")
    v = 'open' if tops_ok == 2 else 'excluded'
    tot[v] += 1; tot['ufr'] += ufr
    print("\n   ".join(line), f"-> {v} ({time.time()-t0:.0f}s)", flush=True)
print("TOTAL", tot, flush=True)
