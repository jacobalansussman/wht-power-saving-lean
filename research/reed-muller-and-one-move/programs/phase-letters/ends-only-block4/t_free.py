"""t_free.py (rm-theory, after the resume): checks for the free-direction lemma.

Lemma (free directions). In the multi-target slice, if the DISTINCT move directions y have F_2-independent rank-ones
y y^T, then deletion alone forces M >= sum_l rank(T_l) (= mW for targets congruent to I).

Part 1 (dependency structure): the F_2-dependencies among the rank-ones z z^T (z != 0) are the words of the dual of the
shortened RM(2,m), i.e. punctured RM(m-3,m). Computed here from scratch (kernel of the (2^m - 1) x m(m+1)/2 matrix):
its minimum weight and the number and shape of the minimum-weight words (expected: 7, the 3-dim subspaces).

Part 2 (independent check of the lemma's consequence, exhaustive): for W = 2 arrays at block m and M moves, every
M-multiset of directions up to GL(m,2), every pair of targets congruent to I (invertible, non-alternating symmetric),
the multi-target exact tests: (i) T_l in the F_2-span of the moves' rank-ones; (ii) deletion |D_Q| >= r_Q for all Q.
Reports the surviving (multiset class, number of target pairs) and whether the survivor's distinct directions are
dependent. Calibration: block 3 2 in 5 (impossible: wider-game Theorem W2), 2 in 6 (plain), 2 in 7 (U_3 exists in the
slice); block 4 2 in 8 (plain) must survive.

Usage: python3 t_free.py deps m | pairs m M
"""
import sys, itertools, time
import numpy as np


def forms(m):
    pairs = [(i, j) for i in range(m) for j in range(i, m)]
    nn = 1 << m
    zb = np.array([[(z >> i) & 1 for i in range(m)] for z in range(nn)], dtype=np.int64)
    mon = np.stack([zb[:, i] * zb[:, j] for (i, j) in pairs], axis=1)          # nn x P : (z z^T)_{ij}, i <= j
    P = len(pairs)
    coeffs = np.array(list(itertools.product([0, 1], repeat=P)), dtype=np.int64)  # all forms
    vals = (coeffs @ mon.T) % 2                                                    # forms x nn : Q(z)
    return pairs, mon, coeffs, vals


def rank2(rows):
    rows = [r for r in rows]; r = 0; basis = []
    for v in rows:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v); r += 1
    return r


def deps(m):
    pairs, mon, coeffs, vals = forms(m)
    nn = 1 << m; n = nn - 1; P = len(pairs)
    # rank-one of z as a P-bit integer
    ro = [int("".join(str(int(b)) for b in mon[z][::-1]), 2) for z in range(1, nn)]
    print(f"m={m}: {n} rank-ones in a space of dim {P}; rank of all = {rank2(ro)}")
    # kernel over F_2 of the n x P matrix (row z = rank-one): Gaussian elimination tracking combinations
    piv = {}            # pivot bit -> (vector, combination mask over z)
    kern = []
    for idx, v in enumerate(ro):
        c = 1 << idx
        while v:
            hb = v.bit_length() - 1
            if hb in piv:
                pv, pc = piv[hb]; v ^= pv; c ^= pc
            else:
                piv[hb] = (v, c); break
        if v == 0:
            kern.append(c)
    k = len(kern)
    print(f"  dependency space dimension {k} (expected 2^m - 1 - m(m+1)/2 = {n - P})")
    if k > 20:
        print("  too many words to list; skipping weight enumeration"); return
    # enumerate all words
    words = [0]
    for b in kern:
        words += [w ^ b for w in words]
    wts = {}
    for w in words[1:]:
        wts[bin(w).count("1")] = wts.get(bin(w).count("1"), 0) + 1
    print("  weight distribution of the dependencies:", dict(sorted(wts.items())))
    dmin = min(wts)
    mins = [w for w in words[1:] if bin(w).count("1") == dmin]
    shapes = {"3-subspace": 0, "other": 0}
    for w in mins:
        pts = [z + 1 for z in range(n) if (w >> z) & 1]
        s = set(pts)
        closed = all((a ^ b) in s for a in pts for b in pts if a != b)
        span_dim = rank2(pts)
        shapes["3-subspace" if (closed and span_dim == 3 and len(pts) == 7) else "other"] += 1
    nsub3 = {3: 1, 4: 15, 5: 155}.get(m)
    print(f"  minimum dependency size {dmin}; minimum words: {shapes} (number of 3-dim subspaces of F_2^{m}: {nsub3})")


def gl_elements(m):
    nn = 1 << m
    out = []
    def rec(cols):
        if len(cols) == m:
            out.append(tuple(cols)); return
        for v in range(1, nn):
            if rank2(cols + [v]) == len(cols) + 1:
                rec(cols + [v])
    rec([])
    return out


def apply(h, z):
    r = 0; i = 0
    while z:
        if z & 1: r ^= h[i]
        z >>= 1; i += 1
    return r


def pairs_check(m, M):
    t0 = time.time()
    pairs, mon, coeffs, vals = forms(m)
    nn = 1 << m; P = len(pairs)
    # targets congruent to I: invertible symmetric with a nonzero diagonal (non-alternating)
    targets = []
    for bits in range(1 << P):
        S = [[0] * m for _ in range(m)]
        for idx, (i, j) in enumerate(pairs):
            if (bits >> idx) & 1: S[i][j] = S[j][i] = 1
        if not any(S[i][i] for i in range(m)): continue
        rows = [sum(S[i][j] << j for j in range(m)) for i in range(m)]
        if rank2(rows) != m: continue
        targets.append(np.array([S[i][j] for (i, j) in pairs], dtype=np.int64))
    Tm = np.stack(targets)                            # t x P
    par = (coeffs @ Tm.T % 2).T.astype(np.int8)       # t x forms : <Q, T>
    print(f"m={m} M={M}: {len(targets)} targets congruent to I; {len(coeffs)} forms")
    G = gl_elements(m)
    print(f"  |GL({m},2)| = {len(G)}")
    seen = set(); reps = []
    for ms in itertools.combinations_with_replacement(range(1, nn), M):
        if ms in seen: continue
        orb = set()
        for h in G:
            orb.add(tuple(sorted(apply(h, z) for z in ms)))
        seen |= orb; reps.append((ms, len(orb)))
    print(f"  {len(reps)} classes of {M}-multisets of directions ({time.time() - t0:.1f} s)")
    surv = []
    for ms, osz in reps:
        cnt = np.zeros(vals.shape[0], dtype=np.int64)
        for z in ms: cnt += vals[:, z]                # |D_Q| with multiplicity
        Z = cnt == 0                                  # forms vanishing on every move
        # T in span  <=>  <Q,T> = 0 for every Q vanishing on all moves (Z also catches even multiplicities: cnt counts)
        # careful: span is about parities; Q with all Q(z_j) = 0 is exactly cnt == 0
        inV = np.where(~(par[:, Z].astype(bool)).any(axis=1))[0]
        if len(inV) == 0: continue
        sub = par[inV].astype(np.int64)               # k x forms
        good = 0; ex = None
        for a in range(len(inV)):
            r = sub[a][None, :] + sub[a:]             # pairs (a, b >= a)
            ok = (cnt[None, :] >= r).all(axis=1)
            nb = int(ok.sum()); good += nb
            if nb and ex is None:
                b = a + int(np.where(ok)[0][0]); ex = (int(inV[a]), int(inV[b]))
        if good:
            dist = sorted(set(ms))
            ro = [int("".join(str(int(x)) for x in mon[z][::-1]), 2) for z in dist]
            indep = rank2(ro) == len(dist)
            same = ex[0] == ex[1]
            surv.append((ms, osz, good, indep, ex))
    print(f"  surviving classes (span + deletion): {len(surv)}")
    for ms, osz, good, indep, ex in surv[:40]:
        print(f"    dirs {ms} (orbit {osz}): {good} unordered target pairs pass; distinct dirs independent: {indep}; "
              f"example pair {ex}{' (equal targets)' if ex[0] == ex[1] else ''}")
    nind = sum(1 for s in surv if s[3])
    print(f"  survivors with independent distinct directions: {nind} (the lemma says 0 when M < 2m)")
    print(f"  done in {time.time() - t0:.1f} s")


if __name__ == "__main__":
    if sys.argv[1] == "deps":
        deps(int(sys.argv[2]))
    else:
        pairs_check(int(sys.argv[2]), int(sys.argv[3]))
