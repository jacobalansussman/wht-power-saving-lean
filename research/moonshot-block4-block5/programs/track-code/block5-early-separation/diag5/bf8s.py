"""bf8s.py (mw2-b5beyond, resumed 14 Oct ~15:10): 2 in 8 at block 5 = the special-position butterfly of bf7b PLUS a second
merge register carrying v-ms-fano's separation, so that the separating (defect) stream is processed to the SAME depth as
the main stream (raw-input injections are one or two direction letters short and never collide: bf7b, bf7r).
  s_j = D b_j x_j;  t1 = D(a1 s1 + a2 s2);  t1' = D(a1 s1 + al * tau a2 s2)  (tau from a menu: odd-order Levi elements and
  powers, the order-3 coordinate cycle, Paulis, Clifford phases);  t2 = D(c2 t1 + al' c2' t1');  t3 = D(c3 t2);
  e_k = D(beta_k t3 + n_side side terms).
(t1, c2, c3, beta_k) is a bf7b special-position skeleton; c2' is a JOINT progress letter for t1''s leads (p, q').
Lemma Delta with r = q^-1 p: Delta(t1') = D (1 - tau-twisted factor) ... is a full-depth defect (PROOF 3d), so the lemma no
longer forbids layer 3.  Layered LM d = 5..1.
usage: bf8s.py ntrials seed n_side nq ntau"""
import sys, random, time
from b5eng import *
from lsearch2 import rand_free, dstep
from census import geo
from bf7 import sigma3, mul, solve
from bf7b import find


def menu(E):
    Lt = Letters(E); m = E.m; out = []
    lv = levi_menu(m)
    for name, Pm in sorted(lv.items()):
        o = gl_order(Pm, m)
        for k in (1, 2):
            if k < o:
                out.append(('%s^%d' % (name, k), Lt.addr(gl_pow(Pm, k, m))))
    out.append(('sigma3', sigma3(E)))
    out.append(('sigma3^2', sigma3(E) @ sigma3(E)))
    for z in (1, 2, 4, 3, 7):
        out.append(('X%d' % z, E.Xs[z])); out.append(('Z%d' % z, Lt.zp(z))); out.append(('S%d' % z, Lt.sphase(z)))
    return out


def main():
    ntr, seed, n_side, nq, ntau = [int(a) for a in sys.argv[1:6]]
    rng = random.Random(seed); t0 = time.time(); hist = {}
    for tr in range(ntr):
        E = Eng(5, seed=seed * 1000 + tr)
        sk, stat = find(E, rng, nq)
        if sk is None:
            print("trial %d: no special-position skeleton %s" % (tr, stat), flush=True); continue
        g1 = E.intern(E.snap(E.D @ sk['b1'])); g2 = E.intern(E.snap(E.D @ sk['b2']))
        p = mul(E, sk['a1'], g1)
        M = menu(E); rng.shuffle(M); tried = 0
        print("trial %d: skeleton %s" % (tr, stat), flush=True)
        for (tname, tau) in M:
            if tried >= ntau:
                break
            a2t = tau @ sk['a2']
            qt = mul(E, a2t, g2)
            # joint progress letter for t1''s leads (p, q')
            c2b = None
            for _ in range(3000):
                b = rand_free(E, rng)
                if dstep(E, b, p) == 1 and dstep(E, b, qt) == 1:
                    c2b = b; break
            if c2b is None:
                continue
            tried += 1
            P = Prog(E, 2)
            s1 = P.add_move([(0, sk['b1'])]); s2 = P.add_move([(1, sk['b2'])])
            t1 = P.add_move([(s1, sk['a1']), (s2, sk['a2'])])
            t1b = P.add_move([(s1, sk['a1']), (s2, a2t)])
            t2 = P.add_move([(t1, sk['c2']), (t1b, c2b)])
            t3 = P.add_move([(t2, sk['c3'])])
            for k in (0, 1):
                summ = [(t3, sk['e%d' % (k + 1)])]
                cand = [0, 1, s1, s2, t1, t1b, t2]
                for _ in range(n_side):
                    summ.append((cand[rng.randrange(len(cand))], rand_free(E, rng)))
                P.add_move(summ)
            if not P.setup_final(cap=3000, maxcells=4e7, maxterms=6000, dmin=1):
                print("  tau %s: refused %s" % (tname, P.why), flush=True); continue
            nf, res = solve(P, rng)
            deep = min([d for d in res if res[d] < 1e-10], default=6)
            hist[deep] = hist.get(deep, 0) + 1
            print("  tau %s: layer(q'^-1 p) %d unknowns %d terms %s closure %s | %s (%.0fs)" % (
                tname, E.layer(E.intern(E.snap(np.conj(E.MATS[qt]).T @ E.MATS[p]))), nf, [len(g) for g in P.gens][2:],
                [(f[0], f[1]) for f in P.fin], " ".join("d%d=%.1e" % (d, res[d]) for d in sorted(res, reverse=True)),
                time.time() - t0), flush=True)
    print("SUMMARY bf8s n_side %d: deepest layer passed (6 = none): %s" % (n_side, sorted(hist.items())), flush=True)


if __name__ == "__main__":
    main()
