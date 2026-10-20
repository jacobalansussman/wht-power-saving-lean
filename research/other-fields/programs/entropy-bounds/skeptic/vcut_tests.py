# vcut_tests.py -- v-nc-gain: recompute nc-gain's quoted cut values with my own graph (vgamma) and LPs (vcut),
# in nc-gain's frame indexing (mapped by subspace equality onto dpg.py's order).
import sys, json, random, time, itertools, glob
import numpy as np
sys.path.insert(0, '.')
import vgamma, vcut, vlemmas

import os
HERE = os.path.dirname(os.path.abspath(__file__))
NC_CALC = os.path.join(HERE, '..', 'lp-and-layouts')
NC_RUNS = os.path.join(HERE, '..', '..', '..', 'data', 'entropy-bounds', 'colab-runs')

class Mapped:
    """my Gamma_m relabelled into dpg's indices."""
    def __init__(self, m):
        g = vgamma.VG(m)
        sys.path.insert(0, NC_CALC)
        import dpg
        verts, idx, adj, _ = dpg.build(m) if m < 4 else vlemmas.dpg_build_fast(m)
        perm = np.array([idx[S] for S in g.L])
        inv = np.argsort(perm)
        self.m = m; self.n = g.n
        self.D = g.D[np.ix_(inv, inv)]
        self.edges = sorted((min(int(perm[u]), int(perm[v])), max(int(perm[u]), int(perm[v]))) for u, v in g.edges)
        self.adj = [np.nonzero(self.D[i] == 1)[0] for i in range(self.n)]

def image(g, us):
    prof = np.stack([g.D[u] for u in us], axis=1)
    keys = {}; lab = np.zeros(g.n, int)
    for v in range(g.n): lab[v] = keys.setdefault(tuple(prof[v].tolist()), len(keys))
    E = sorted({(min(lab[a], lab[b]), max(lab[a], lab[b])) for a, b in g.edges if lab[a] != lab[b]})
    return len(keys), E, lab

def run(which):
    t0 = time.time()
    if '2' in which:
        g = Mapped(2)
        # (a) all-cuts exact vs flow form, random layouts W = 3..6 (shared frames allowed), Lemma W prediction
        rng = random.Random(7); rows = None
        stats = {}
        for W in (3, 4, 5, 6):
            k = 0; mism = 0; lw_bad = 0; gaps = 0; nst = 0
            for tr in range(40):
                P = []
                for _ in range(W):
                    s = rng.randrange(15); t = int(rng.choice(list(np.nonzero(g.D[s] == 2)[0]))); P.append((s, t))
                fl, lb, _ = vcut.cut_flow(g.n, g.edges, P)
                R, rhs = vcut.all_cuts_rows(g.n, g.edges, P)
                ex_f, ex = vcut.cut_all(g.n, g.edges, P, (R, rhs))
                if abs(fl - ex_f) > 1e-6: mism += 1
                st = vlemmas.two_class(g, P) is not None; nst += st
                if st and ex < 2 * W: lw_bad += 1
                if ex_f < 2 * W - 1e-6: gaps += 1
                k += 1
            stats[W] = (k, mism, nst, lw_bad, gaps)
            print(f"m=2 W={W}: {k} layouts; flow-form vs all-cuts mismatches {mism}; two-class stretchable {nst}, Lemma W violated {lw_bad}; cut gaps {gaps}", flush=True)
        P = [(0, 7), (3, 1), (4, 10), (0, 10), (10, 7)]
        print(f"m=2 session-1 gap layout {P}: directed {vcut.cut_flow(g.n, g.edges, P, True)[:2]}, undirected {vcut.cut_flow(g.n, g.edges, P)[:2]}, "
              f"all-cuts exact {vcut.cut_all(g.n, g.edges, P)}; distances {[int(g.D[s,t]) for s,t in P]}", flush=True)
    if '3' in which:
        g = Mapped(3)
        for P, claim in [([(4, 123), (25, 33), (78, 48)], 9), ([(34, 15), (65, 28), (126, 113)], 9), ([(114, 97), (17, 41), (16, 29)], 9),
                         ([(0, 95), (0, 80), (80, 95), (84, 60)], 11)]:
            ok = all(g.D[s, t] == 3 for s, t in P)
            u = vcut.cut_flow(g.n, g.edges, P); d = vcut.cut_flow(g.n, g.edges, P, True)
            print(f"m=3 {P}: opposite {ok}; e(D) {vlemmas.edef(g, P)[0] if len(P)==3 else '-'}; two-class {vlemmas.two_class(g, P) is not None}; "
                  f"undirected primal {u[0]:.6f} LB {u[1]:.6f}; directed primal {d[0]:.6f} LB {d[1]:.6f}; claimed (undirected) {claim}  ({time.time()-t0:.0f}s)", flush=True)
        for P, claim in [([(0, 95), (0, 80), (80, 95), (12, 80), (84, 60)], 13), ([(0, 128), (0, 124), (124, 128), (8, 124), (76, 44)], 14)]:
            d = vcut.cut_flow(g.n, g.edges, P, True)
            print(f"m=3 W=5 {P}: directed primal {d[0]:.6f} LB {d[1]:.6f}; claimed {claim}  ({time.time()-t0:.0f}s)", flush=True)
        P = [(4, 123), (25, 33), (78, 48)]
        best = 0
        for o in itertools.product((0, 1), repeat=3):
            us = [P[i][o[i]] for i in range(3)]
            nH, EH, lab = image(g, us)
            v = vcut.cut_flow(nH, EH, [(lab[s], lab[t]) for s, t in P])
            best = max(best, v[1])
        print(f"m=3 {P}: best three-end image LB {best:.6f} (claimed 8.5)", flush=True)
    if '4' in which:
        g = Mapped(4)
        files = sorted(glob.glob(NC_RUNS + '/job_w3s4*__w3screen_m4_s*.json'))
        recs = [r for f in files for r in json.load(open(f))['image']]
        rng = random.Random(3); sample = rng.sample(recs, 12)
        lows = sorted(recs, key=lambda r: r['value_float'])[:6]
        for r in lows + sample:
            P = [tuple(p) for p in r['P']]; o = r['orient']
            us = [P[i][o[i]] for i in range(3)]
            nH, EH, lab = image(g, us)
            v = vcut.cut_flow(nH, EH, [(lab[s], lab[t]) for s, t in P])
            comp = vlemmas.edef(g, P)[0]
            print(f"m=4 {P} orient {o}: e(D) {comp}, opposite {all(g.D[s,t]==4 for s,t in P)}; image {nH} pts {len(EH)} edges; "
                  f"my primal {v[0]:.6f} LB {v[1]:.6f}; nc-gain exact {r['value_float']:.6f}; agree {abs(v[0]-r['value_float'])<1e-5}  ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    run(sys.argv[1])
