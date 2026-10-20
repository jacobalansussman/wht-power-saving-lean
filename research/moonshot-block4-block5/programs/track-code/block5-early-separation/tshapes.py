#!/usr/bin/env python3
"""tshapes.py (mw2-b5beyond): Theorem T applied to dependency SHAPES (posets of moves).

Enumerates naturally labelled posets of T moves (as theory-guide/shapes.py), keeps those passing D2 both ways (>= W moves
of depth >= m, >= W of co-depth >= m), then applies the two-layer killer KT (PROOF.md section 1) forwards and on the
mirrored (transposed) poset:

  Top set  Tp = {k : depth(k) >= m}, input set S = {j : codepth(j) >= m}.  If |Tp| = |S| = W, the top system
  A_T Z Gamma_S = K I_W is square, so A_T, Z, Gamma_S are invertible (direct finiteness).  KT kills the shape if some
  j* in S satisfies
   (a) every move k with a chain of >= m-1 moves starting at j* (j* and k included) has depth >= m;
   (b) for every 'extra' input j (codepth(j) = m-1), every k with a chain of >= m-1 moves from j has depth >= m;
   (c) no top move p that is the second-to-last move of a >= m chain from j* to a top move is reached from an extra
       input j by a chain of >= m-1 moves.
  Then the next-layer equation forces rho_{p j*} = 0 for those p, hence the whole column j* of Z is 0: contradiction.

usage: tshapes.py W T m [maxprint]"""
import sys, itertools
W, T, m = map(int, sys.argv[1:4])
lim = int(sys.argv[4]) if len(sys.argv) > 4 else 50


def downsets(preds, n):
    out = [0]
    for i in range(n):
        new = []
        for D in out:
            if preds[i] & ~D == 0:
                new.append(D | (1 << i))
        out += new
    return out


def longest(preds, n):
    """L[i][j] = number of moves on the longest chain from i to j (i <= j), 0 if not comparable; L[i][i] = 1"""
    L = [[0] * n for _ in range(n)]
    for j in range(n):
        L[j][j] = 1
        for i in range(j):
            best = 0
            # chain i -> ... -> p -> j with p an immediate-or-not predecessor of j
            for p in range(i, j):
                if preds[j] >> p & 1 and L[i][p]:
                    best = max(best, L[i][p] + 1)
            L[i][j] = best
    return L


def depths(L, n):
    d = [max(L[i][j] for i in range(n)) for j in range(n)]
    c = [max(L[i][j] for j in range(n)) for i in range(n)]
    return d, c


def kt_kills(preds, n):
    L = longest(preds, n)
    d, c = depths(L, n)
    Tp = [k for k in range(n) if d[k] >= m]
    S = [j for j in range(n) if c[j] >= m]
    if len(Tp) != W or len(S) != W:
        return False, (len(Tp), len(S))
    extra = [j for j in range(n) if c[j] == m - 1]
    # (b)
    for j in extra:
        for k in range(n):
            if L[j][k] >= m - 1 and d[k] < m:
                return False, "b"
    for js in S:
        # (a)
        if any(L[js][k] >= m - 1 and d[k] < m for k in range(n)):
            continue
        # (c): second-to-last moves p of >= m chains js -> top k
        bad = False
        for k in Tp:
            for p in range(n):
                if preds[k] >> p & 1 and L[js][p] >= m - 1:
                    if any(L[j][p] >= m - 1 for j in extra):
                        bad = True
        if not bad:
            return True, js
    return False, "a/c"


def mirror(preds, n):
    """transposed poset, relabelled naturally (reverse order)"""
    succ = [0] * n
    for j in range(n):
        for i in range(n):
            if preds[j] >> i & 1:
                succ[i] |= 1 << j
    # new label of old i is n-1-i
    out = [0] * n
    for i in range(n):
        s = 0
        for j in range(n):
            if succ[i] >> j & 1:
                s |= 1 << (n - 1 - j)
        out[n - 1 - i] = s
    return out


def canon(preds, n):
    L = longest(preds, n)
    d, c = depths(L, n)
    up = [sum(1 for j in range(n) if preds[j] >> i & 1) for i in range(n)]
    dn = [bin(preds[i]).count("1") for i in range(n)]
    key = [(d[i], c[i], up[i], dn[i]) for i in range(n)]
    classes = {}
    for i in range(n):
        classes.setdefault(key[i], []).append(i)
    ks = sorted(classes)
    best = None
    for perms in itertools.product(*[itertools.permutations(classes[k]) for k in ks]):
        order = [x for p in perms for x in p]
        pos = {v: i for i, v in enumerate(order)}
        rel = tuple(sorted((pos[j], pos[i]) for i in range(n) for j in range(n) if preds[i] >> j & 1))
        if best is None or rel < best:
            best = rel
    return (tuple(key[i] for i in sorted(range(n), key=lambda i: key[i])), best)


passD2 = {}
total = 0


def rec(preds):
    global total
    n = len(preds)
    if n == T:
        total += 1
        L = longest(preds, n)
        d, c = depths(L, n)
        if sum(1 for x in d if x >= m) >= W and sum(1 for x in c if x >= m) >= W:
            cf = canon(preds, n)
            if cf not in passD2:
                passD2[cf] = list(preds)
        return
    for D in downsets(preds, n):
        rec(preds + [D])


rec([])
surv = []
killed = 0
for cf, preds in sorted(passD2.items()):
    k1, w1 = kt_kills(preds, T)
    k2, w2 = kt_kills(mirror(preds, T), T)
    if k1 or k2:
        killed += 1
    else:
        L = longest(preds, T); d, c = depths(L, T)
        Tp = sum(1 for x in d if x >= m); S = sum(1 for x in c if x >= m)
        surv.append((cf, Tp, S))


def hasse(rel, n):
    R = set(rel)
    return sorted((a, b) for (a, b) in R if not any((a, c) in R and (c, b) in R for c in range(n)))


print("W=%d T=%d m=%d: naturally labelled posets %d; classes passing D2 both ways %d; killed by KT (fwd or mirror) %d;"
      " SURVIVE %d" % (W, T, m, total, len(passD2), killed, len(surv)))
for (k, rel), Tp, S in surv[:lim]:
    print("  |Top|=%d |S|=%d (depth,codepth):" % (Tp, S), [x[:2] for x in k], " covers:", hasse(rel, T))
