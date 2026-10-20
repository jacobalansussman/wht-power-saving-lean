"""v-w3-blueprint: exhaustive mod-2 check of the Theorem R' condition in the ENGINE game.
Cell B g <-> Lagrangian sigma(g)^-1 L_Z; input distance = 4 - dim(Lambda cap L_Z).  For a middle single-term content
t = l D(z) a: M = sigma(t) L_Z (distance 1 from L_Z), N' = sigma(t) L_1 (L_1 = kernel Lagrangian), any complement of M.
g in S iff c_g N' = N';  junk g (not in S) at full distance iff c_g M cap N' = 0.
R' condition fails iff S != H and every junk g has c_g M cap N' != 0."""
import sys, random, itertools, numpy as np
np.seterr(all='ignore')
sys.path.insert(0, __import__('os').path.dirname(__file__))
from groups import G, z7, z2cube
n, m = 16, 4
def pauli(a, b):
    M = np.zeros((n, n), dtype=complex)
    for y in range(n): M[y ^ a, y] = (-1) ** bin(b & y).count('1')
    return M
def ident_pauli(P):
    a = int(np.argmax(np.abs(P[:, 0])))
    b = 0
    base = P[a, 0]
    for j in range(m):
        y = 1 << j
        v = P[y ^ a, y] / base
        if abs(v + 1) < 1e-6: b |= 1 << j
        elif abs(v - 1) > 1e-6: return None
    Q = pauli(a, b) * base
    return (a, b) if np.allclose(P, Q, atol=1e-6) else None
def sigma(g):
    gi = np.conj(g).T
    cols = []
    for j in range(8):
        a, b = (1 << j, 0) if j < 4 else (0, 1 << (j - 4))
        r = ident_pauli(g @ pauli(a, b) @ gi)
        assert r is not None
        cols.append(r[0] | (r[1] << 4))
    return cols  # image of basis vector j as 8-bit int (x bits 0-3, z bits 4-7)
def app(S, v):
    out = 0
    for j in range(8):
        if (v >> j) & 1: out ^= S[j]
    return out
def inv_map(S):
    img = {app(S, v): v for v in range(256)}
    return [img[1 << j] for j in range(8)]
def span(vs):
    s = {0}
    for v in vs: s |= {x ^ v for x in s}
    return frozenset(s)
def img(S, L): return frozenset(app(S, v) for v in L)
LZ = span([1 << (4 + j) for j in range(4)])
def cell(g): return img(inv_map(sigma(g)), LZ)
def dimcap(A, B): return int(round(np.log2(len(A & B))))
if __name__ == '__main__':
    rng = random.Random(5)
    # 1. validate the dictionary on random words: distance (row support) and frame equality vs Lagrangian
    ok = 0
    words = []
    for it in range(60):
        g = np.eye(n, dtype=complex)
        for k in range(rng.randrange(0, 6)):
            g = G.mul(G.rand_engine(rng), G.D[rng.randrange(1, 16)], g)
        g = G.mul(G.rand_engine(rng), g)
        words.append(g)
        assert G.dist_in(g) == 4 - dimcap(cell(g), LZ), 'distance mismatch'
        ok += 1
    same = 0; agree = 0
    for i in range(len(words)):
        for j in range(i + 1, len(words)):
            fr = G.is_mono_any(G.mul(words[j], np.conj(words[i]).T))
            lg = cell(words[i]) == cell(words[j])
            agree += (fr == lg); same += fr
    print('dictionary check: distances agree on', ok, 'words; B g = B g\' (mod scalars) iff same Lagrangian:', agree, 'of', len(words) * (len(words) - 1) // 2, 'pairs agree (', same, 'equal pairs)')
    L1 = cell(G.K)
    print('kernel Lagrangian distance', 4 - dimcap(L1, LZ))
    # all Lagrangians by orbit of L_Z under sigma of generators
    gens = [sigma(G.D[z]) for z in range(1, 16)] + [sigma(G.rand_engine(rng)) for _ in range(30)]
    seen = {LZ}; fr = [LZ]
    while fr:
        nf = []
        for L in fr:
            for S in gens:
                L2 = img(S, L)
                if L2 not in seen: seen.add(L2); nf.append(L2)
        fr = nf
    Lag = list(seen)
    print('Lagrangians found:', len(Lag), '(expected 2295)')
    Ms = [L for L in Lag if dimcap(L, LZ) == 3]
    print('distance-1 Lagrangians M:', len(Ms), '(expected 30)')
    for name, els in (('Z7', z7()[0]), ('Z2^3', z2cube()[0])):
        C = [sigma(c) for c in els]
        assert all(img(c, LZ) == LZ and img(c, L1) == L1 for c in C), 'c_g not in the Levi'
        from collections import Counter
        cnt = Counter(); cex = []
        for M in Ms:
            for Np in Lag:
                if len(M & Np) != 1: continue
                S = [g for g in range(len(C)) if img(C[g], Np) == Np]
                junk = [g for g in range(len(C)) if g not in S]
                full = [g for g in junk if len(img(C[g], M) & Np) == 1]
                cnt[(len(S), len(full), len(junk) - len(full))] += 1
                if junk and not full: cex.append((M, Np, S))
        print(name, ': pairs (M, N\') =', sum(cnt.values()), ' (|S|, #junk at full distance, #junk nearer):')
        for k, v in sorted(cnt.items()): print('    ', k, v)
        print('  R\' condition FAILS (S != H, no full-distance junk) for', len(cex), 'pairs')
