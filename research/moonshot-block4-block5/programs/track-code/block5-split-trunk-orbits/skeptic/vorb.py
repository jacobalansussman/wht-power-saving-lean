"""vorb.py (v-w7-b5, own code): decoupled 2-in-8, profile (4,2)/(2,4), lam = (3,3).  For sigma in Stab(M0) of three kinds
(fix: sigma N_a = N_a, sigma N_b = N_b; swap: N_a <-> N_b; cycle: N_a -> N_b but N_b -> X != N_a, so N_a, N_b lie in one
sigma-cycle of length l >= 3), decide whether BOTH outputs' own-branch Lemma-Q rows have an alive non-class producer among the
orbit's streams (Lemma E1-sigma: stream i+m produces stream i's row of output o through branch k' iff sigma^m N_c(o) = N_k';
in the 2-in-8 the producer's only main rows sit at Y = sigma^m A_o, predecessor N_k', for both arrays).  Status rules of
w7-b5 NOTE sec. 5 (my implementation).  usage: vorb.py seed nconf nlist"""
import sys, random
from vgeo import *

def build_config(rng):
    M0 = frozenset(range(32))
    for _ in range(200):
        L0 = rng.choice(sorted(neighbours(M0), key=lambda X: min(X ^ M0)))
        N1 = rng.choice([X for X in sorted(neighbours(L0), key=lambda X: min(X ^ L0)) if dist(M0, X) == 2])
        nm = sorted(neighbours(M0), key=lambda X: min(X ^ M0))
        anch = [X for X in nm if dist(X, N1) == 3]
        br = [X for X in sorted(neighbours(N1), key=lambda X: min(X ^ N1)) if dist(M0, X) == 3]
        rng.shuffle(anch)
        for P1 in anch:
            for Q in anch:
                if Q == P1:
                    continue
                Na = [X for X in br if dist(P1, X) == 4 and dist(Q, X) == 2]
                Nb = [X for X in br if dist(P1, X) == 2 and dist(Q, X) == 4]
                if not Na or not Nb:
                    continue
                Na = rng.choice(Na); Nb = rng.choice(Nb)
                A1s = [X for X in neighbours(Na) if dist(P1, X) == 5 and dist(Q, X) == 3]
                A2s = [X for X in neighbours(Nb) if dist(Q, X) == 5 and dist(P1, X) == 3]
                if A1s and A2s:
                    return dict(M0=M0, L0=L0, N1=N1, P1=P1, Q=Q, Na=Na, Nb=Nb, A1=rng.choice(A1s), A2=rng.choice(A2s))
    return None

def row_status(C, l, Y, Npred):
    anc = C['P1'] if l == 1 else C['Q']
    A = C['A1'] if l == 1 else C['A2']
    lam = dist(anc, Y); nu = dist(anc, Npred)
    br = [C['Na'], C['Nb']]
    if lam == 5:
        return 'FREE' if Y == A else 'DEAD'
    if lam == 4:
        if Y in br and dist(anc, Y) == 4:
            return 'AUTO'
        return 'NEED' if nu == 3 else 'DEAD'
    if lam == 3:
        if Y == C['N1'] or (Y in br and dist(anc, Y) == 3):
            return 'AUTO'
        return 'NEED'
    return 'LOW'

def orbit_alive(C, sig):
    """returns (out1_alive_nonclass, out2_alive_nonclass, info)"""
    Na, Nb, A1, A2 = C['Na'], C['Nb'], C['A1'], C['A2']
    seq = [(Na, Nb, A1, A2)]
    cur = seq[0]
    while True:
        cur = tuple(gframe(sig, X) for X in cur)
        if cur == seq[0]:
            break
        seq.append(cur)
        if len(seq) > 3000:
            return None
    per = len(seq)
    res = {}
    for o, (c, A) in ((1, ('a', A1)), (2, ('b', A2))):
        Nc = Na if c == 'a' else Nb
        alive_nc = []; twins = 0; dead = 0
        for m in range(per):
            na, nb, a1, a2 = seq[m]
            img = na if c == 'a' else nb          # sigma^m N_c
            Y = a1 if o == 1 else a2              # sigma^m A_o
            for kp, Nk in (('a', Na), ('b', Nb)):
                if img != Nk:
                    continue
                if m == 0 and kp == c:
                    continue                      # the stream's own class column
                if kp == c and Y == A:
                    twins += 1; continue          # twin (class family member)
                st = (row_status(C, 1, Y, Nk), row_status(C, 2, Y, Nk))
                if 'DEAD' in st:
                    dead += 1
                else:
                    alive_nc.append((m, kp, st))
        res[o] = (alive_nc, twins, dead)
    return res, per

def main():
    seed, nconf, nlist = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rng = random.Random(seed)
    tot = {}
    for ci in range(nconf):
        C = build_config(rng)
        assert C is not None
        for k in ('L0', 'N1', 'P1', 'Q', 'Na', 'Nb', 'A1', 'A2'):
            assert isotropic(C[k])
        G = []
        for _ in range(nlist):
            g = rand_stab_M0(rng)
            G.append((g, gframe(g, C['Na']), gframe(g, C['Nb'])))
        assert all(is_symplectic(G[i][0]) for i in range(5)) and all(gframe(G[i][0], C['M0']) == C['M0'] for i in range(5))
        byA = {}; byB = {}; byAB = {}
        for idx, (g, ga, gb) in enumerate(G):
            byA.setdefault(ga, []).append(idx); byB.setdefault(gb, []).append(idx); byAB.setdefault((ga, gb), []).append(idx)
        sigs = {'fix': [], 'swap': [], 'cycle': []}
        for idx, (g, ga, gb) in enumerate(G):
            for j in byB.get(ga, ()):        # h N_b = g N_a  -> sigma = h^-1 g maps N_a -> N_b
                h = G[j][0]
                if G[j][1] == gb:
                    sigs['swap'].append((idx, j))
                elif len(sigs['cycle']) < 400:
                    sigs['cycle'].append((idx, j))
            for j in byAB.get((ga, gb), ()):
                if j != idx and len(sigs['fix']) < 400:
                    sigs['fix'].append((idx, j))
        for kind in ('fix', 'swap', 'cycle'):
            for (i, j) in sigs[kind][:400]:
                sig = gmul(ginv(G[j][0]), G[i][0])
                assert gframe(sig, C['M0']) == C['M0']
                if kind == 'fix':
                    assert gframe(sig, C['Na']) == C['Na'] and gframe(sig, C['Nb']) == C['Nb']
                else:
                    assert gframe(sig, C['Na']) == C['Nb']
                r = orbit_alive(C, sig)
                if r is None:
                    continue
                res, per = r
                a1 = len(res[1][0]) > 0; a2 = len(res[2][0]) > 0
                key = (kind, 'out1 ' + ('alive' if a1 else 'none'), 'out2 ' + ('alive' if a2 else 'none'))
                tot[key] = tot.get(key, 0) + 1
                if kind == 'cycle' and a1 and a2 and tot[key] <= 3:
                    print('  example cycle sigma conf %d: period %d, out1 producers %s, out2 producers %s' %
                          (ci, per, res[1][0][:3], res[2][0][:3]), flush=True)
        print('conf %d: sigma counts fix %d swap %d cycle %d' % (ci, len(sigs['fix']), len(sigs['swap']), len(sigs['cycle'])), flush=True)
    for k in sorted(tot):
        print('   ', k, tot[k])

main()
