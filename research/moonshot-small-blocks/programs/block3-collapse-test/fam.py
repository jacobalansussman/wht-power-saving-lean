"""fam.py -- ms-xw3.  Shape families for 2 in 5 / 3 in 7 / 3 in 8 at block 3, solved by als3.bcd.

usage: python3 -I fam.py <shape> <menu> <pattern> <starts> <sweeps> <seed> <target K|plant> [FS e.g. 3,4]
shapes:
  SHT  (2 in 5, shared top):  y0 = M0(v0), y1 = M1(G1 y0), y2 = M2(v1), y3 = M3(G3 y2),
                              y4 = D_{z4}( GEN over menu on v0, v1, y0..y3 )
                              out_l = GEN over menu on v0, v1, y0..y4            (blocks: c4, R)
       pattern = "z0,C0;G1;z1,C1;z2,C2;G3;z3,C3;z4"
  SHB  (2 in 5, shared bottom): y0 = D_{z0}( GEN on v0, v1 ), y1 = M1(G1 y0 + v0-word), y2 = M2(G2 y0 + v1-word),
                              y3 = D_{z3}( GEN on y0, y1, y2, inputs ), y4 = D_{z4}( GEN on y0..y3, inputs )
                              out_l = GEN on everything     (blocks: c0, c3, c4, R)
       pattern = "z0;G1;z1,C1;G2;z2,C2;z3;z4"
  CH5  (2 in 5, chain): y_i = D_{z_i}^{C_i}( GEN over menu on inputs and all earlier y ) for i = 0..4;
                              out_l = GEN on everything (blocks c0..c4, R).   pattern = "z0,C0;z1,C1;...;z4,C4"
"""
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from xw3lib import Game, menu as MENU
from als3 import Prog, bcd, targets_K, residual, per_output, fl


def build(shape, pat, Mw, W=2):
    P = [q for q in pat.split(";")]
    def zc(q):
        a = q.split(","); return int(a[0]), (a[1] if len(a) > 1 else "1")
    moves = []; outs = []
    if shape == "SHT":
        (z0, C0), G1, (z1, C1), (z2, C2), G3, (z3, C3), z4 = zc(P[0]), P[1], zc(P[2]), zc(P[3]), P[4], zc(P[5]), int(P[6])
        moves = [dict(z=z0, C=C0, content=[("in", "1", 0, 1.0)]),
                 dict(z=z1, C=C1, content=[("y", G1, 0, 1.0)]),
                 dict(z=z2, C=C2, content=[("in", "1", 1, 1.0)]),
                 dict(z=z3, C=C3, content=[("y", G3, 2, 1.0)]),
                 dict(z=z4, C="1", content=[("Gin", ("c4", "in"), Mw)] + [("Gy", ("c4", j), Mw, j) for j in range(4)])]
        nmov = 5
    elif shape == "SHB":
        z0, G1, (z1, C1), G2, (z2, C2), z3, z4 = int(P[0]), P[1], zc(P[2]), P[3], zc(P[4]), int(P[5]), int(P[6])
        moves = [dict(z=z0, C="1", content=[("Gin", ("c0", "in"), Mw)]),
                 dict(z=z1, C=C1, content=[("y", G1, 0, 1.0), ("in", "1", 0, 1.0)]),
                 dict(z=z2, C=C2, content=[("y", G2, 0, 1.0), ("in", "1", 1, 1.0)]),
                 dict(z=z3, C="1", content=[("Gin", ("c3", "in"), Mw)] + [("Gy", ("c3", j), Mw, j) for j in range(3)]),
                 dict(z=z4, C="1", content=[("Gin", ("c4", "in"), Mw)] + [("Gy", ("c4", j), Mw, j) for j in range(4)])]
        nmov = 5
    elif shape == "ST":
        # generalized shared top: "<chain_0>|...|<chain_{W-1}>|<top>", chain = "z,C;G;z,C;..." (private moves on
        # array l, each fed by the previous one through the fixed word G), top = "z;z;..." shared GENERAL moves
        parts = pat.split("|"); chains, top = parts[:-1], parts[-1]
        assert len(chains) == W
        for l, ch in enumerate(chains):
            q = ch.split(";")
            prev = None
            for k in range(0, len(q), 2):
                z, C = zc(q[k])
                if prev is None: cont = [("in", "1", l, 1.0)]
                else: cont = [("y", q[k - 1], prev, 1.0)]
                moves.append(dict(z=z, C=C, content=cont)); prev = len(moves) - 1
        for k, zt in enumerate(top.split(";")):
            i = len(moves); blk = "t%d" % k
            moves.append(dict(z=int(zt.split(",")[0]), C=(zt.split(",")[1] if "," in zt else "1"),
                              content=[("Gin", (blk, "in"), Mw)] + [("Gy", (blk, j), Mw, j) for j in range(i)]))
        nmov = len(moves)
    elif shape.startswith("CH"):
        nmov = len(P)
        for i in range(nmov):
            z, C = zc(P[i])
            moves.append(dict(z=z, C=C, content=[("Gin", ("c%d" % i, "in"), Mw)] +
                              [("Gy", ("c%d" % i, j), Mw, j) for j in range(i)]))
    else:
        raise ValueError(shape)
    for l in range(W):
        outs.append([("Gin", ("R", l, "in"), Mw)] + [("Gy", ("R", l, j), Mw, j) for j in range(nmov)])
    return moves, outs


def main():
    shape, mk, pat = sys.argv[1], sys.argv[2], sys.argv[3]
    starts, sweeps, seed = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    target = sys.argv[7] if len(sys.argv) > 7 else "K"
    FS = [int(a) for a in (sys.argv[8] if len(sys.argv) > 8 else "3,4").split(",")]
    W = int(os.environ.get("MSW", "2"))
    Mw = MENU(mk)
    G = Game(FS, W=W, S=1, seed=1000 + seed)
    moves, outs = build(shape, pat, Mw, W)
    P = Prog(G, moves, outs, W)
    rng = np.random.default_rng(seed)
    t0 = time.time()
    if target == "plant":
        pv = P.init_vals(rng, scale=0.5, density=0.15)
        o, _, _ = P.evaluate(pv, None)
        T = [o[l].c for l in range(W)]
    else:
        T = targets_K(G, W)
    print(f"# {shape} menu={mk} ({len(Mw)} words) pat={pat} FS={FS} target={target} seed={seed} "
          f"blocks={P.blocks} unknowns={sum(s[2] for s in P.slots.values())}", flush=True)
    best = 9.0
    for s in range(starts):
        vals = P.init_vals(rng)
        res, hist = bcd(P, T, vals, sweeps=sweeps)
        per = per_output(P, vals, T)
        best = min(best, res)
        print(f"  start {s}: residual {res:.4e} per-output {['%.3f' % q for q in per]} sweeps {len(hist)} "
              f"hist {['%.2e' % h for h in hist[::max(1, len(hist) // 6)]]} ({time.time() - t0:.0f}s)", flush=True)
        if res < 1e-7:
            print("CANDIDATE " + json.dumps(dict(shape=shape, menu=mk, pat=pat, seed=seed, start=s, res=res)), flush=True)
    print(f"# DONE best {best:.4e} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
