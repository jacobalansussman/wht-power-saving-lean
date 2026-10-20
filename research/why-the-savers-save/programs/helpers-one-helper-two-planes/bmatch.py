#!/usr/bin/env python3
"""Exact max b-matching in a general graph via the edge-gadget reduction to maximum matching and Edmonds' blossom
algorithm (classic O(V^3) BFS form). Node P with capacity b(P) -> b(P) copies; edge {a,b} -> two gadget nodes ea, eb
with ea-eb and a_i-ea, eb-b_j. Max matching = |E| + max b-matching.
Used by direct_vs_exact.py (imported) on the relaxed-H2 compatibility graph."""
from collections import deque
def max_matching(n, adj):
    match = [-1] * n; p = [-1] * n; base = list(range(n))
    def lca(a, b):
        used = [False] * n
        while True:
            a = base[a]; used[a] = True
            if match[a] == -1: break
            a = p[match[a]]
        while True:
            b = base[b]
            if used[b]: return b
            b = p[match[b]]
    def mark_path(v, b, ch, blossom):
        while base[v] != b:
            blossom[base[v]] = blossom[base[match[v]]] = True
            p[v] = ch; ch = match[v]; v = p[match[v]]
    def find_path(root):
        nonlocal p
        used = [False] * n; p = [-1] * n
        for i in range(n): base[i] = i
        used[root] = True; q = deque([root])
        while q:
            v = q.popleft()
            for to in adj[v]:
                if base[v] == base[to] or match[v] == to: continue
                if to == root or (match[to] != -1 and p[match[to]] != -1):
                    cur = lca(v, to); blossom = [False] * n
                    mark_path(v, cur, to, blossom); mark_path(to, cur, v, blossom)
                    for i in range(n):
                        if blossom[base[i]]:
                            base[i] = cur
                            if not used[i]: used[i] = True; q.append(i)
                elif p[to] == -1:
                    p[to] = v
                    if match[to] == -1: return to
                    used[match[to]] = True; q.append(match[to])
        return -1
    # greedy init
    for v in range(n):
        if match[v] == -1:
            for to in adj[v]:
                if match[to] == -1: match[to] = v; match[v] = to; break
    for v in range(n):
        if match[v] == -1:
            u = find_path(v)
            while u != -1:
                pv = p[u]; ppv = match[pv]; match[u] = pv; match[pv] = u; u = ppv
    return sum(1 for x in match if x != -1) // 2
def max_bmatching(cap, edges):
    nodes = sorted(cap); copies = {}; k = 0
    for P in nodes:
        copies[P] = list(range(k, k + cap[P])); k += cap[P]
    adj = []
    def node():
        nonlocal k
        adj.append([]); return None
    N = k + 2 * len(edges); adj = [[] for _ in range(N)]
    def link(a, b): adj[a].append(b); adj[b].append(a)
    for i, (a, b) in enumerate(edges):
        ea, eb = k + 2 * i, k + 2 * i + 1; link(ea, eb)
        for c in copies[a]: link(c, ea)
        for c in copies[b]: link(eb, c)
    return max_matching(N, adj) - len(edges)
if __name__ == '__main__':
    # tiny self-test: triangle with capacities 1 -> 1; capacities 2 -> 3
    print(max_bmatching({1: 1, 2: 1, 3: 1}, [(1, 2), (2, 3), (1, 3)]), max_bmatching({1: 2, 2: 2, 3: 2}, [(1, 2), (2, 3), (1, 3)]))
