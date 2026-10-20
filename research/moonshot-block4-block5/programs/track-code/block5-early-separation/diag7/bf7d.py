"""bf7d.py (mw2-b5beyond, resumed 14 Oct ~15:30): COLLISION CLOSURE on special-position 2-in-7 butterflies.
Round 0: bf7b skeleton + bf7c's designed injections.  Each round: layered solve; at the failing layer d take the rows that keep
a residual (coordinate l, group element g, in the readout frame) and try to give each one a colliding partner:
  (i)  an injection of x_l at t1 or t2 whose stream reaches an end in the coset B g  (W D beta in B g),
  (ii) a side term into an end from an earlier generator term h with D beta h in B g.
beta is found by frame matching against a fixed pool of random free letters (batched numpy).  Up to `rounds` rounds.
usage: bf7d.py ntrials seed nq npool rounds"""
import sys, random, time
from b5eng import *
from lsearch2 import rand_free, dstep
from bf7 import mul, solve
from bf7b import find
from bf7c import word


class Prog2(Prog):
    def setup_final(self, cap=3000, maxcells=4e7, maxterms=6000, dmin=1):
        E = self.E
        if sum(len(T) for T in self.gens) > maxterms:
            self.why = 'terms'; return False
        tok = []; gd = []
        for i, T in enumerate(self.gens):
            d = {}
            for j, t in enumerate(T):
                d[t] = len(tok) + 1; tok.append((i, j))
            gd.append(d)
        sizes = [len(T) for T in self.gens]
        self.off = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)
        self.fin = []; self.rowkeys = []
        for l in range(self.W):
            cl = E.closure_layered(gd, (l, E.KID), cap, dmin)
            if cl is None:
                self.why = 'cap'; return False
            rows, cols, ent = cl
            if len(rows) * max(1, len(cols)) > maxcells:
                self.why = 'cells'; return False
            lay = np.zeros(len(rows), dtype=np.int64); keys = [None] * len(rows)
            for t, r in rows.items():
                lay[r] = E.layer(t[1]); keys[r] = t
            Ea = np.array([(r, c, tok[co - 1][0], tok[co - 1][1]) for (r, c, co) in ent], dtype=np.int64).reshape(-1, 4)
            self.fin.append((len(rows), len(cols), Ea[:, 0], Ea[:, 1], self.off[Ea[:, 2]] + Ea[:, 3], lay))
            self.rowkeys.append(keys)
        return True


def residual_rows(P, alpha, d, tol=1e-6):
    fl = P.flat(alpha); out = []
    for a, (nr, nc, R, C, F, lay) in enumerate(P.fin):
        keep = lay >= d
        idx = -np.ones(nr, dtype=np.int64); idx[keep] = np.arange(keep.sum()); back = np.nonzero(keep)[0]
        sel = keep[R]
        R2, C2, F2 = idx[R[sel]], C[sel], F[sel]
        nr2 = int(keep.sum())
        A = np.zeros((nr2, nc), dtype=complex)
        np.add.at(A, (R2, C2), fl[F2])
        t = np.zeros(nr2, dtype=complex); t[idx[0]] = 1
        r = (A @ np.linalg.lstsq(A, t, rcond=1e-10)[0] - t) if nc else -t
        for j in np.nonzero(np.abs(r) > tol)[0]:
            out.append(P.rowkeys[a][back[j]])
    return out


def solve_keep(P, rng):
    """layered LM like bf7.solve, also returning the best alpha at the failing layer"""
    fi = P.free_idx(); nf = len(fi); res = {}; xbest = None
    for d in range(5, 0, -1):
        best = (np.inf, None)
        for st in range(3):
            x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
            x, cst = P.lm(d, x0, fi, iters=80) if nf else (x0, P.vp(np.ones(P.nal, dtype=complex), d, [], want_J=False)[0])
            if not np.isscalar(cst):
                cst = float(np.sum(np.abs(cst) ** 2))
            if cst < best[0]:
                best = (cst, x)
            if best[0] < 1e-20:
                break
        res[d] = best[0]; xbest = best[1]
        if best[0] > 1e-10:
            return res, d, P.full_alpha(xbest)
    return res, None, None


def main():
    ntr, seed, nq, npool, rounds = [int(a) for a in sys.argv[1:6]]
    rng = random.Random(seed); t0 = time.time(); hist = {}
    for tr in range(ntr):
        E = Eng(5, seed=seed * 1000 + tr)
        sk, stat = find(E, rng, nq)
        if sk is None:
            print("trial %d: no skeleton %s" % (tr, stat), flush=True); continue
        pool = np.stack([rand_free(E, rng) for _ in range(npool)])          # (npool, n, n)
        Ws = {(1, 0): word(E, [sk['e1'], sk['c3'], sk['c2']]), (1, 1): word(E, [sk['e2'], sk['c3'], sk['c2']]),
              (2, 0): word(E, [sk['e1'], sk['c3']]), (2, 1): word(E, [sk['e2'], sk['c3']])}

        def match(Wmat, target_mat):
            """pool letters b with Wmat D b target^-1 monomial"""
            G = np.einsum('ij,njk->nik', Wmat @ E.D, pool) @ np.conj(target_mat).T
            nzr = (np.abs(G) > 1e-9).sum(axis=2)
            ok = np.nonzero((nzr == 1).all(axis=1))[0]
            for i in ok:
                if E.is_mono(E.snap(G[i])):
                    return pool[i]
            return None
        inj = {1: [], 2: []}; extra = [[], []]; tried = set(); last = None
        for rd in range(rounds + 1):
            P = Prog2(E, 2)
            s1 = P.add_move([(0, sk['b1'])]); s2 = P.add_move([(1, sk['b2'])])
            t1 = P.add_move([(s1, sk['a1']), (s2, sk['a2'])] + inj[1])
            t2 = P.add_move([(t1, sk['c2'])] + inj[2])
            t3 = P.add_move([(t2, sk['c3'])])
            ends = [P.add_move([(t3, sk['e1'])] + extra[0]), P.add_move([(t3, sk['e2'])] + extra[1])]
            if not P.setup_final():
                print("  round %d: refused %s" % (rd, P.why), flush=True); break
            res, dfail, alpha = solve_keep(P, rng)
            print("trial %d round %d: injections t1 %d t2 %d, end side terms %s, closure %s | %s (%.0fs)" % (
                tr, rd, len(inj[1]), len(inj[2]), [len(x) for x in extra], [(f[0], f[1]) for f in P.fin],
                " ".join("d%d=%.1e" % (d, res[d]) for d in sorted(res, reverse=True)), time.time() - t0), flush=True)
            last = 1 if dfail is None else dfail + 1
            if dfail is None or rd == rounds:
                break
            rows = residual_rows(P, alpha, dfail)
            added = 0
            for (l, g) in rows[:12]:
                if (l, g) in tried:
                    continue
                tried.add((l, g)); gm = E.MATS[g]
                done = False
                for (site, k) in ((2, 0), (2, 1), (1, 0), (1, 1)):          # injections of x_l, stream through end k
                    b = match(Ws[(site, k)], gm)
                    if b is not None:
                        inj[site].append((l, b)); added += 1; done = True; break
                if done:
                    continue
                for k in (0, 1):                                            # side term into end k from a term h
                    for i in (s1, s2, t1, t2):
                        for (l2, h) in P.gens[i]:
                            if l2 != l:
                                continue
                            b = match(E.I, gm @ np.conj(E.MATS[h]).T)      # D b h in B g
                            if b is not None:
                                extra[k].append((i, b)); added += 1; done = True; break
                        if done:
                            break
                    if done:
                        break
            print("  round %d: %d leftover rows at layer %d, %d partners added" % (rd, len(rows), dfail, added), flush=True)
            if added == 0:
                break
        hist[last] = hist.get(last, 0) + 1
    print("SUMMARY bf7d: deepest layer passed after collision closure (6 = none): %s" % sorted(hist.items()), flush=True)


if __name__ == "__main__":
    main()
