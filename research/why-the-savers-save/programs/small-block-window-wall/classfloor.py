#!/usr/bin/env python3
"""an2-smallnet / classfloor.py (own code, stdlib; python3 -I).  A floor on the centre charge of ONE call that holds
for every decomposition of the index code, not only the hub's (proof in NOTE.md sec. 3):
  if a set Z of ports splits into groups Q_1..Q_p such that labels inside a group are pairwise non-orthogonal
  (B = 1) and labels in different groups are orthogonal (B = 0), then the centres of the call drop at least |Z|
  levels in all.  (Each group needs |Q| centres touching it on both sides, because the code must be the identity on
  a group; a centre touching c groups has a label pairing of F_2 rank >= c, and its drop is at least that rank.)
So a design whose ports admit such a Z with |Z| = n can never charge less than one level per port per call: it saves
nothing in any chain (gain 2 per pair, at least two calls per pair).  This program finds large Z (exact for n <= 20
by branch and bound, greedy with restarts above) for the designs that fit a block <= 16 (sigma <= 5)."""
import json, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import designs as D

def label_gram(L):
    n = len(L); norm = D.ip(L[0], L[0]); unit = norm // 2
    G = [[D.ip(a, b) for b in L] for a in L]
    return [[(1 + G[i][j] // unit) % 2 for j in range(n)] for i in range(n)]

def valid_add(B, groups, v):
    """can port v join the cluster structure? returns the group index it must join (or -1 for a new group), or None."""
    hit = [gi for gi, g in enumerate(groups) if any(B[v][u] == 1 for u in g)]
    if not hit:
        return -1
    if len(hit) > 1:
        return None
    g = groups[hit[0]]
    return hit[0] if all(B[v][u] == 1 for u in g) else None

def greedy(B, n, order):
    groups = []
    for v in order:
        r = valid_add(B, groups, v)
        if r is None: continue
        if r == -1: groups.append([v])
        else: groups[r].append(v)
    return groups

def exact(B, n):
    best = [0, None]
    def rec(i, groups, size):
        if size + (n - i) <= best[0]: return
        if i == n:
            best[0] = size; best[1] = [g[:] for g in groups]; return
        r = valid_add(B, groups, i)
        if r is not None:
            if r == -1:
                groups.append([i]); rec(i + 1, groups, size + 1); groups.pop()
            else:
                groups[r].append(i); rec(i + 1, groups, size + 1); groups[r].pop()
        rec(i + 1, groups, size)
    rec(0, [], 0)
    return best

def main():
    random.seed(15)
    small = [('D4', D.lines(D.roots_D(4))), ('D5', D.lines(D.roots_D(5))), ('D6', D.lines(D.roots_D(6))),
             ('A3', D.lines(D.roots_A(3))), ('A4', D.lines(D.roots_A(4))), ('A5', D.lines(D.roots_A(5)))]
    out = []
    for name, L in small:
        n = len(L); B = label_gram(L)
        if n <= 20:
            size, groups = exact(B, n); how = 'exact'
        else:
            size, groups = 0, None
            for _ in range(4000):
                order = list(range(n)); random.shuffle(order)
                g = greedy(B, n, order); s = sum(len(x) for x in g)
                if s > size: size, groups = s, g
            how = 'greedy, 4000 restarts (a lower bound on the best Z)'
        print(f"{name}: n = {n}; largest split set found |Z| = {size} ({how}); groups {sorted(len(g) for g in groups)};"
              f" centre charge per call >= {size}, per port >= {size / n:.3f}" + ("  -> saves nothing in any chain" if size == n else ""),
              flush=True)
        out.append(dict(design=name, n=n, Z=size, groups=sorted(len(g) for g in groups), method=how,
                        floor_per_port=size / n))
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'classfloor.json'), 'w'), indent=1)

if __name__ == '__main__':
    main()
