#!/usr/bin/env python3
"""sweep.py: for each toy (G, P), every paid set made of one double coset P h P (with its inverse), every target
double coset at the largest distance m >= 2: look for two-array units with M = 2m-1 and 2m-2 moves (and calibrate on
M = 2m) by Levenberg-Marquardt from random starts.  Prints one line per case."""
import sys, time, warnings, itertools
warnings.filterwarnings("ignore")
import numpy as np
from toylab import *

def hyperoct(m):
    n = 2 * m
    gens = [cyc(n, [0, 1])] + [cyc(n, [2 * i, 2 * i + 2], [2 * i + 1, 2 * i + 3]) for i in range(m - 1)]
    G = perm_group(gens, f"B{m}")
    P = G.subgroup([G.idx[g] for g in gens[1:]])
    return G, P, None

def gl3():
    E = lambda i, j: tuple((1 << r) | ((1 << j) if r == i else 0) for r in range(3))
    gens = [E(0, 1), E(1, 2), E(1, 0), E(2, 1)]
    G = f2_group(gens, "GL(3,2)")
    P = G.subgroup([G.idx[E(0, 1)], G.idx[E(1, 2)]])    # upper unitriangular = Borel (order 8)
    return G, P, None

def gl3_pt():
    E = lambda i, j: tuple((1 << r) | ((1 << j) if r == i else 0) for r in range(3))
    gens = [E(0, 1), E(1, 2), E(1, 0), E(2, 1)]
    G = f2_group(gens, "GL(3,2)")
    # stabiliser of the line spanned by e0 under v -> A v: columns... use the action on vectors
    def apply(A, v): return sum((bin(A[i] & v).count("1") & 1) << i for i in range(3))
    P = [i for i, A in enumerate(G.elems) if apply(A, 1) == 1]
    return G, P, None

def a5_sub(k):
    G = perm_group([cyc(5, [0, 1, 2]), cyc(5, [0, 1, 2, 3, 4])], "A5")
    if k == "Z2": P = G.subgroup([G.idx[cyc(5, [0, 1], [2, 3])]])
    elif k == "Z3": P = G.subgroup([G.idx[cyc(5, [0, 1, 2])]])
    elif k == "V4": P = G.subgroup([G.idx[cyc(5, [0, 1], [2, 3])], G.idx[cyc(5, [0, 2], [1, 3])]])
    elif k == "S3": P = G.subgroup([G.idx[cyc(5, [0, 1, 2])], G.idx[cyc(5, [0, 1], [3, 4])]])
    elif k == "Z5": P = G.subgroup([G.idx[cyc(5, [0, 1, 2, 3, 4])]])
    elif k == "D5": P = G.subgroup([G.idx[cyc(5, [0, 1, 2, 3, 4])], G.idx[cyc(5, [1, 4], [2, 3])]])
    return G, P, None

def s4_sub(k):
    G = perm_group([cyc(4, [0, 1]), cyc(4, [0, 1, 2, 3])], "S4")
    if k == "Z2": P = G.subgroup([G.idx[cyc(4, [0, 1])]])
    elif k == "Z2b": P = G.subgroup([G.idx[cyc(4, [0, 1], [2, 3])]])
    elif k == "Z3": P = G.subgroup([G.idx[cyc(4, [0, 1, 2])]])
    elif k == "Z4": P = G.subgroup([G.idx[cyc(4, [0, 1, 2, 3])]])
    elif k == "V4n": P = G.subgroup([G.idx[cyc(4, [0, 1])], G.idx[cyc(4, [2, 3])]])
    return G, P, None

def s5_sub(k):
    G = perm_group([cyc(5, [0, 1]), cyc(5, [0, 1, 2, 3, 4])], "S5")
    if k == "D5": P = G.subgroup([G.idx[cyc(5, [0, 1, 2, 3, 4])], G.idx[cyc(5, [1, 4], [2, 3])]])
    elif k == "S3": P = G.subgroup([G.idx[cyc(5, [0, 1, 2])], G.idx[cyc(5, [0, 1])]])
    elif k == "S4": P = G.subgroup([G.idx[cyc(5, [0, 1, 2, 3])], G.idx[cyc(5, [0, 1])]])
    elif k == "F20": P = G.subgroup([G.idx[cyc(5, [0, 1, 2, 3, 4])], G.idx[cyc(5, [1, 2, 4, 3])]])
    elif k == "S2S3": P = G.subgroup([G.idx[cyc(5, [0, 1])], G.idx[cyc(5, [2, 3, 4])], G.idx[cyc(5, [2, 3])]])
    return G, P, None

GAMES = {
    "D5": lambda: dihedral(5)[:2] + (None,), "D7": lambda: dihedral(7)[:2] + (None,), "D9": lambda: dihedral(9)[:2] + (None,),
    "J42": lambda: johnson(4, 2)[:2] + (None,), "J52": lambda: johnson(5, 2)[:2] + (None,),
    "J62": lambda: johnson(6, 2)[:2] + (None,), "J63": lambda: johnson(6, 3)[:2] + (None,),
    "B3": lambda: hyperoct(3), "B4": lambda: hyperoct(4), "GL3B": gl3, "GL3pt": gl3_pt,
    "A5Z2": lambda: a5_sub("Z2"), "A5Z3": lambda: a5_sub("Z3"), "A5V4": lambda: a5_sub("V4"), "A5S3": lambda: a5_sub("S3"),
    "A5Z5": lambda: a5_sub("Z5"), "A5D5": lambda: a5_sub("D5"),
    "S4Z2": lambda: s4_sub("Z2"), "S4Z2b": lambda: s4_sub("Z2b"), "S4Z3": lambda: s4_sub("Z3"), "S4Z4": lambda: s4_sub("Z4"),
    "S5D5": lambda: s5_sub("D5"), "S5S3": lambda: s5_sub("S3"), "S5F20": lambda: s5_sub("F20"), "S5S2S3": lambda: s5_sub("S2S3"),
    "Sp4": lambda: sp4()[:2] + (None,),
}

def run(name, restarts=4, Mextra=(1, 2), max_cases=12, calib=True):
    t0 = time.time()
    G, P, _ = GAMES[name]()
    reps = irreps(G)
    fg0 = FrameGame(G, P, [])  if False else None
    # double cosets of P
    Pset = set(P)
    seen = -np.ones(G.n, dtype=np.int32); dcs = []
    for g in range(G.n):
        if seen[g] < 0:
            k = len(dcs); dcs.append(g)
            for p in P:
                for q in P: seen[G.m(G.m(p, g), q)] = k
    cases = 0
    done_paid = set()
    for k, h in enumerate(dcs):
        if h in Pset: continue
        kinv = int(seen[G.inv[h]])
        key = tuple(sorted({k, kinv}))
        if key in done_paid: continue
        done_paid.add(key)
        hreps = [dcs[x] for x in key]
        fg = FrameGame(G, P, hreps)
        if fg.ncos and min(fg.dist) < 0:  # not connected
            print(f"{name}: paid {key} not connected; skip", flush=True); continue
        m = fg.diam
        if m < 2:
            print(f"{name}: |G|={G.n} |P|={len(P)} frames={fg.ncos} paid {key} diam {m}: skip", flush=True); continue
        tops = sorted({int(seen[g]) for g in range(G.n) if fg.ell[g] == m})
        S = Solver(fg, reps, W=2)
        for tk in tops:
            T = dcs[tk]
            for Mm in ([0] if calib else []) + list(Mextra):
                M = 2 * m - Mm
                if M < m: continue
                best = (np.inf, None, None)
                for seq in itertools.product(range(len(hreps)), repeat=M):
                    if len(hreps) > 1 and seq[0] != 0: continue
                    hs = [hreps[i] for i in seq]
                    res, Fs, allr = search_lm(S, hs, T, restarts=(2 if Mm == 0 else restarts), seed=M * 7 + tk)
                    if res < best[0]: best = (res, seq, allr)
                    if res < 1e-10: break
                tag = "UNIT" if best[0] < 1e-10 else "none"
                print(f"{name}: |G|={G.n} |P|={len(P)} frames={fg.ncos} paid {key} m={m} target dc{tk} M={M} (2m-{Mm}): "
                      f"best rel.res {best[0]:.2e} seq {best[1]} -> {tag}   [{time.time()-t0:.0f}s]", flush=True)
            cases += 1
            if cases >= max_cases: return


# ---- more games (added 12:45)
def s6_sub(k):
    G = perm_group([cyc(6, [0, 1]), cyc(6, list(range(6)))], "S6")
    if k == "S3wrS2":
        P = G.subgroup([G.idx[cyc(6, [0, 1])], G.idx[cyc(6, [0, 1, 2])], G.idx[cyc(6, [0, 3], [1, 4], [2, 5])]])
    elif k == "S2wrS3":
        P = G.subgroup([G.idx[cyc(6, [0, 1])], G.idx[cyc(6, [0, 2], [1, 3])], G.idx[cyc(6, [0, 2, 4], [1, 3, 5])]])
    elif k == "S5":
        P = G.subgroup([G.idx[cyc(6, [0, 1])], G.idx[cyc(6, [0, 1, 2, 3, 4])]])
    elif k == "S3xS3":
        P = G.subgroup([G.idx[cyc(6, [0, 1])], G.idx[cyc(6, [0, 1, 2])], G.idx[cyc(6, [3, 4])], G.idx[cyc(6, [3, 4, 5])]])
    return G, P, None

def a6_sub(k):
    G = perm_group([cyc(6, [0, 1, 2]), cyc(6, [1, 2, 3, 4, 5])], "A6")
    if k == "A5": P = G.subgroup([G.idx[cyc(6, [0, 1, 2])], G.idx[cyc(6, [0, 1, 2, 3, 4])]])
    elif k == "S4": P = G.subgroup([G.idx[cyc(6, [0, 1, 2, 3], [4, 5])], G.idx[cyc(6, [0, 1], [4, 5])]])
    elif k == "33": P = G.subgroup([G.idx[cyc(6, [0, 1, 2])], G.idx[cyc(6, [3, 4, 5])], G.idx[cyc(6, [0, 1], [3, 4])]])
    return G, P, None

def prod_game(a, b):
    """direct product of two dihedral games D_a x D_b (frames: product of cycles), as permutations on a + b points"""
    n = a + b
    ra = cyc(n, list(range(a))); sa = tuple([(-i) % a for i in range(a)] + list(range(a, n)))
    rb = cyc(n, list(range(a, n))); sb = tuple(list(range(a)) + [a + ((-(i - a)) % b) for i in range(a, n)])
    G = perm_group([ra, sa, rb, sb], f"D{a}xD{b}")
    P = G.subgroup([G.idx[sa], G.idx[sb]])
    return G, P, None

GAMES.update({
    "S6S3wr2": lambda: s6_sub("S3wrS2"), "S6S5": lambda: s6_sub("S5"), "S6S3xS3": lambda: s6_sub("S3xS3"),
    "Match3all": lambda: s6_sub("S2wrS3"),
    "A6A5": lambda: a6_sub("A5"), "A6S4": lambda: a6_sub("S4"), "A633": lambda: a6_sub("33"),
    "D5xD5": lambda: prod_game(5, 5), "D3xD5": lambda: prod_game(3, 5), "D5xD7": lambda: prod_game(5, 7),
})

if __name__ == "__main__":
    for name in sys.argv[1:]:
        run(name)
