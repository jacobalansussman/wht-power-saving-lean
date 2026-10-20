"""rd2S.py -- (resumed) Lemma D' + D2 / D2' on SIGMA-PAIR cores (sweep_sig.core_sig, variants sig1 / sig0): the shape
BRIEF5 names.  Usage: rd2S.py variant ncores seed ...  (replays sweep_sig's cores: same seeds, same rng sequence).
Difference from cross cores: v's array-0 layer-3 part is (c1 + c2 sigma) q_10 (q_20 = sigma q_10, one frame), so letters
of eps_e can cancel each other's images there.  Since c1 + c2 sigma is invertible (F_00 is), every coset e' <sigma>
that meets supp(eps_e) must contain a letter E with E q_10 produced in t_e: the production condition for that frame
is relaxed to "e' sigma^k meets the production set for some k" (OR over the sigma-orbit, computed exactly by right
multiplication on all of Stab(L0)).  Any other frame holding two or more layer-3 terms of v is dropped (no condition),
which keeps the test a necessary one.
rd2X.py -- (resumed) Lemma D2 / D2' on the esearchX-open CROSS cores.  Built on esearchX.py: w5-twoshared: Lemma D' (the extension of Lemma D to SHARED layer 2) tested exhaustively over all
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

# Lemma D2 / D2' for cross cores [NOTE s3-s4]: for every Lemma-D' solution e' of a junk J = y_e q_kj (array a), put
# w_J = e'^-1 J (a term of v), N1 = y_e w_J.  Dead if lambda(N1) = 3 and N1 has no partner: Lag(w_J) is no u main frame
# Lag(p_k'a) (content cancellation / u -> t_e) and Lag(N1) is no frame Lag(q_k'a) of a layer-3 term of v (eps_e).
# D2': the other top reads w_J: G = y_e' w_J; dead if lambda(G) = 3, Lag(w_J) no u frame, Lag(G) no layer-3 V frame.

from sweep_sig import core_sig

def letter_PB(S):
    """(GL index, Bm) of a symplectic engine letter S (in Stab(L0))"""
    P = tuple(S[i] & 15 for i in range(4))
    Pm = {}
    for v in range(16):
        w = 0
        for j in range(4):
            if (v >> j) & 1: w ^= P[j]
        Pm[v] = w
    Pinv = {w: v for v, w in Pm.items()}
    BPe = [S[i] >> 4 for i in range(4)]
    Bm = [[0] * 4 for _ in range(4)]
    for j in range(4):
        x = Pinv[1 << j]; val = 0
        for i in range(4):
            if (x >> i) & 1: val ^= BPe[i]
        for r in range(4): Bm[r][j] = (val >> r) & 1
    assert sp.symp_from(list(P), Bm) == tuple(S)
    return GLidx[P], Bm

def right_mult_perm(S):
    """for the letter tau = S: index arrays (pi, c) with e'_(i,b) tau = e'_(pi[i], b ^ c[i])"""
    ip2, B2 = letter_PB(S); P2 = GL[ip2]
    pi = np.zeros(len(GL), dtype=np.int64); cc = np.zeros(len(GL), dtype=np.int64)
    for i, P1 in enumerate(GL):
        def ap(cols, x):
            r = 0
            for j in range(4):
                if (x >> j) & 1: r ^= cols[j]
            return r
        P = tuple(ap(P1, P2[j]) for j in range(4))
        pi[i] = GLidx[P]
        inv1 = {ap(P1, v): v for v in range(16)}
        Q = [[(inv1[1 << c] >> r) & 1 for c in range(4)] for r in range(4)]      # Q = P1^-1, Q[r][c]
        C = [[sum(Q[a][r] * B2[a][b] * Q[b][c] for a in range(4) for b in range(4)) & 1 for c in range(4)] for r in range(4)]
        cc[i] = symidx(C)
    return pi, cc

def shift(mask, pi, cc):
    """out[e'] = mask[e' tau]"""
    b = np.arange(1024)[None, :]
    return mask[pi[:, None], b ^ cc[:, None]]

variant = sys.argv[1]; ncores = int(sys.argv[2]); seeds = [int(x) for x in sys.argv[3:]]
D = sp.D
TALLY = {}
for seed in seeds:
  rng = random.Random(seed)
  for c in range(ncores):
    t0 = time.time()
    cr = core_sig(rng, variant)
    if cr is None: print("core failed", flush=True); continue
    R0, R1, a11, a20, tg, o = cr
    b = {0: R0['b'][0], 1: R1['b'][0]}; s = {0: R0['s'], 1: R1['s']}
    g = {1: R0['g'][0], 2: R1['g'][0]}
    a = {(1, 0): R0['a'][0], (2, 1): R1['a'][0], (2, 0): a20[0]}
    if a11 is not None: a[(1, 1)] = a11[0]
    p = {kj: sp.mv(a[kj], s[kj[1]]) for kj in a}
    q = {kj: sp.mv(g[kj[0]], p[kj]) for kj in a}
    mu = {0: R0['mu'][0], 1: R1['mu'][0]}; route = {0: (1, 0), 1: (2, 1)}
    lay3 = [kj for kj in q if sp.lam(q[kj]) == 3]
    sig = sp.comp(q[(2, 0)], sp.inv(q[(1, 0)]))
    assert sp.lag(sig) == sp.L0, "sigma is not a free letter"
    ns = 1; T = sig
    while T != sp.I8: T = sp.comp(T, sig); ns += 1
    pi, cc = right_mult_perm(sig)
    _r = random.Random(7)
    for _ in range(30):
        i_, b_ = _r.randrange(len(GL)), _r.randrange(1024)
        assert sp.symp_from(list(GL[pi[i_]]), SYM[b_ ^ cc[i_]]) == sp.comp(sp.symp_from(list(GL[i_]), SYM[b_]), sig), "right_mult_perm check failed"
    groups = {}
    for kj in lay3: groups.setdefault((kj[1], sp.lag(q[kj])), []).append(kj)
    print(f"seed {seed} core {c} ({variant}): sigma exact order {o}, symplectic order {ns}, tags {tg}, layer-3 terms of v {sorted(lay3)}, frame groups {sorted(len(v) for v in groups.values())}", flush=True)
    core_dead = False
    for e in (0, 1):
        y = sp.comp(D, mu[e]); X = sp.comp(D, sp.comp(mu[e], D)); yo = sp.comp(D, mu[1 - e])
        def prod_cond(k1, j1):
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
            return cnd
        allk = np.ones((len(GL), 1024), dtype=bool); notes = []
        for key_, members in groups.items():
            if len(members) == 1:
                allk &= prod_cond(*members[0]); continue
            if sorted(members) == [(1, 0), (2, 0)]:
                base = prod_cond(1, 0); acc = base.copy(); cur = base
                for _k in range(1, ns):
                    cur = shift(cur, pi, cc); acc |= cur
                allk &= acc; notes.append(f"sigma-orbit relaxation on {members} (orbit {ns})")
            else:
                notes.append(f"group {members} dropped")
        top_dead = False
        for kj in q:
            if kj == route[e]: continue
            J = sp.comp(y, q[kj]); L = sp.lam(J)
            if L != 2: print(f"   top {e} junk {kj}: layer {L}{' (kernel frame)' if sp.lag(J) == sp.LK else ''} (not tested)", flush=True); continue
            arr = kj[1]
            pfr = {sp.lag(p[(k2, arr)]) for k2 in (1, 2) if (k2, arr) in p}
            if sp.lag(J) in pfr: print(f"   top {e} junk {kj}: in a u frame (not tested)", flush=True); continue
            qfr = {sp.lag(q[kj2]) for kj2 in lay3 if kj2[1] == arr}
            cj = maps(basis(img_L0(D)), img_L0(sp.comp(J, sp.inv(s[arr]))))
            for l in (1, 2):
                cj |= maps(basis(img_L0(sp.comp(D, sp.comp(g[l], D)))), img_L0(J))
            for kj2 in q:
                if kj2[1] == arr and sp.lam(q[kj2]) == 2:
                    cj |= single(sp.comp(J, sp.inv(q[kj2])))
            sols = np.argwhere(cj & allk)
            st = {}
            for (ip, ib) in sols:
                Se = sp.symp_from(list(GL[ip]), SYM[ib])
                wJ = sp.comp(sp.inv(Se), J)
                N1 = sp.comp(y, wJ); G2 = sp.comp(yo, wJ)
                uf = sp.lag(wJ) in pfr
                d2 = sp.lam(N1) == 3 and not uf and sp.lag(N1) not in qfr
                d2p = sp.lam(G2) == 3 and not uf and sp.lag(G2) not in qfr
                k_ = (sp.lam(N1), sp.lam(G2), 'u' if uf else '-', 'dead' if (d2 or d2p) else 'ALIVE')
                st[k_] = st.get(k_, 0) + 1
            alive = sum(v for k_, v in st.items() if k_[3] == 'ALIVE')
            print(f"   top {e} junk {kj}: {len(sols)} solutions; (lam N1, lam y_other w_J, u-frame, verdict): {st}; alive {alive}", flush=True)
            if alive == 0: top_dead = True
        print(f"   top {e}: {'EXCLUDED' if top_dead else 'passes'} {notes}", flush=True)
        core_dead |= top_dead
    v = 'EXCLUDED' if core_dead else 'OPEN after round 2'
    TALLY[v] = TALLY.get(v, 0) + 1
    print(f"seed {seed} core {c}: {v} ({time.time()-t0:.0f}s)", flush=True)
print("TALLY", TALLY, flush=True)
