#!/usr/bin/env python3
"""t_e10.py: Theorem E10 (bordered rank identity) on genuine engine units and random engine words, in rho_f = g^(x)f;
the zero-move fake (T_0 = P_K I_W, M = 0) against every rank identity; K conjugate to the free letter P_K."""
import random, sys, time
import numpy as np
from engE import *

def lam_list(rng): return [1, 1j, -1, -1j, 0.37 + 0.81j, 0]

def check_word(name, E, word, R, W, outs, f, rng, unit):
    rep = lambda g: kpow(g, f)
    d = E.N ** f
    tr = transfer(word, R, W, d, rep, outs)
    tgt = np.kron(np.eye(W), rep(E.K))
    is_unit = np.abs(tr["total"] - tgt).max() < 1e-9
    worst = 0
    for lam in lam_list(rng):
        lhs = rank(bordered(tr, lam))
        rhs = sum(tr["ks"]) + rank(tr["total"] - lam * np.eye(W * d))
        worst = max(worst, abs(lhs - rhs))
    T0dev = float(np.abs(tr["T0"] - np.eye(W * d)).max())
    print("%-34s m=%d f=%d d=%4d M=%2d sum k=%5d  unit=%s (expected %s)  transfer-form err=%.1e  "
          "max|rank N - (sum k + rank(Z-lam))|=%d  |T0 - I|max=%.3f" % (name, E.m, f, d, len(tr["ks"]), sum(tr["ks"]),
          is_unit, unit, tr["tf_err"], worst, T0dev))
    return worst == 0 and tr["tf_err"] < 1e-8 and is_unit == unit

def words(E, rng):
    m, N = E.m, E.N
    out = []
    # plain, W = 2
    w = [("move", l, E.dir[1 << j]) for l in range(2) for j in range(m)]
    out.append(("plain W=2", w, 2, 2, [0, 1], True))
    # pair trick (theory-guide t_hb.py), W = 2, R = 4
    w = [("mix", 2, [(0, 1), (1, 1)]), ("mix", 3, [(0, 1), (1, -1)])]
    for r in (2, 3): w += [("move", r, E.dir[1 << j]) for j in range(m)]
    w += [("mix", 0, [(2, .5), (3, .5)]), ("mix", 1, [(2, .5), (3, -.5)])]
    out.append(("pair trick W=2", w, 4, 2, [0, 1], True))
    # wasteful unit with T_0 = X^{e1} != I: dir e1 . dir e1 . X . dir e1 = dir e1^5 = dir e1, then the other directions
    w = [("move", 0, E.dir[1]), ("free", 0, E.free["X1"]), ("move", 0, E.dir[1]), ("move", 0, E.dir[1])]
    w += [("move", 0, E.dir[1 << j]) for j in range(1, m)]
    out.append(("wasteful unit W=1 (T0=X)", w, 1, 1, [0], True))
    # address-map conjugated plain: g^-1 dir(g z) g = dir(z)
    Gm = E.rand_gl(rng); Gi = np.linalg.inv(Gm)
    # find images of basis vectors under the permutation x -> Gx
    img = {x: int(np.argmax(np.abs(Gm[:, x]))) for x in range(N)}
    w = [("free", 0, Gm)] + [("move", 0, E.dir[img[1 << j]]) for j in range(m)] + [("free", 0, Gi)]
    out.append(("address-map conjugated plain W=1", w, 1, 1, [0], True))
    # random words with phases, shifts, address maps, mixes (not units)
    names = list(E.free.keys())
    for t in range(3):
        w = []
        for _ in range(14):
            u = rng.random(); r = rng.randrange(4)
            if u < 0.3: w.append(("move", r, E.dir[rng.randrange(1, N)]))
            elif u < 0.45: w.append(("free", r, E.free[rng.choice(names)]))
            elif u < 0.55: w.append(("free", r, E.rand_gl(rng)))
            else:
                src = rng.sample(range(4), rng.randint(1, 3))
                w.append(("mix", r, [(s_, complex(rng.gauss(0, 1), rng.gauss(0, 1))) for s_ in src]))
        out.append(("random word %d" % t, w, 4, 2, [0, 1], False))
    return out

if __name__ == "__main__":
    rng = random.Random(7)
    allok = True
    for m, fs in ((2, (1, 2, 3)), (3, (1, 2))):
        E = Eng(m)
        W_ = words(E, rng)
        for f in fs:
            for (name, w, R, W, outs, unit) in W_:
                allok &= check_word(name, E, w, R, W, outs, f, rng, unit)
    print("ALL E10 CHECKS PASS" if allok else "SOME CHECK FAILED")
    # zero-move fake
    print("\nzero-move fake: T_0 = P_K (x) I_W with P_K = prod_j P_{e_j,0} = diag(i^{|x|}) (a free letter), M = 0")
    for m, fs in ((2, (1, 2, 3, 4)), (3, (1, 2, 3))):
        E = Eng(m)
        PK = np.eye(E.N, dtype=complex)
        for j in range(m): PK = PK @ E.free["P%d_0" % (1 << j)]
        Hm = kpow(np.array([[1, 1], [1, -1]]) / np.sqrt(2), m)
        conj = np.abs(Hm @ PK @ Hm - E.K).max()
        # H in G up to a scalar: S . dir(e1) . S at m = 1 equals omega * H
        E1 = Eng(1); S = E1.free["P1_0"]; Hs = S @ E1.dir[1] @ S
        ratio = Hs / (np.array([[1, 1], [1, -1]]) / np.sqrt(2))
        for f in fs:
            rp = lambda g: kpow(g, f)
            d = E.N ** f
            ok = all(rank(rp(PK) - lam * np.eye(d)) == rank(rp(E.K) - lam * np.eye(d)) for lam in (1, 1j, -1, -1j))
            print("  m=%d f=%d: rank(rho(P_K) - lam) = rank(rho(K) - lam) for lam in {1,i,-1,-i}: %s" % (m, f, ok))
        print("  m=%d: |H^(x)m P_K H^(x)m - K|max = %.1e;  S.dir(e1).S / H entries = %s (constant => H in G up to a scalar)"
              % (m, conj, np.round(ratio.flatten(), 6)))
