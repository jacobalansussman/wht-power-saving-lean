# Independent exact check of an interval lift at m = 3: own frame definition, own integer rank (Bareiss), mod-p rank as
# a second opinion. Input: JSON {name: 3 integer rows of length 6}; names "0", "F", "pXYZ" (point), "tABC..." (line).
import json, random
Y = json.load(open("e3_Y.json"))
def brank(rows):  # exact integer rank by fraction-free elimination
    M = [list(r) for r in rows]; r = 0; nc = len(M[0])
    for c in range(nc):
        p = next((i for i in range(r, len(M)) if M[i][c]), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        for i in range(r + 1, len(M)):
            if M[i][c]:
                M[i] = [M[r][c] * a - M[i][c] * b for a, b in zip(M[i], M[r])]
        r += 1
    return r
def prank(rows, P):
    M = [[v % P for v in r] for r in rows]; r = 0
    for c in range(len(M[0])):
        p = next((i for i in range(r, len(M)) if M[i][c]), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]; inv = pow(M[r][c], P - 2, P)
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c] * inv % P; M[i] = [(a - f * b) % P for a, b in zip(M[i], M[r])]
        r += 1
    return r
def bits(s): return int(s, 2)
pts = list(range(1, 8))
def subspace(name):
    if name == "0": return {0}
    if name == "F": return set(range(8))
    if name[0] == "p": return {0, bits(name[1:])}
    q = name[1:]; return {0} | {bits(q[i:i + 3]) for i in range(0, len(q), 3)}
def dot(a, b): return bin(a & b).count("1") & 1
def frame(s):
    ann = [z for z in range(8) if all(dot(z, x) == 0 for x in s)]
    return {(x << 3) | z for x in s for z in ann}
def red(v): return sum((v[i] & 1) << (5 - i) for i in range(6))
def rank2(vs):
    basis = []
    for v in vs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return len(basis), basis
for k, rows in Y.items():
    S = subspace(k); fr = frame(S)
    r2, _ = rank2([red(v) for v in rows])
    assert r2 == 3 and all(red(v) in fr for v in rows) and len(fr) == 8, k
names = list(Y)
P1, P2 = (1 << 61) - 1, 1000000007
E = 0; bad = []; nedge = 0; endE = 0
for a in names:
    for b in names:
        if a >= b: continue
        A, B = subspace(a), subspace(b)
        if not ((A < B and len(B) == 2 * len(A)) or (B < A and len(A) == 2 * len(B))): continue
        nedge += 1
        r = brank(Y[a] + Y[b]); assert r == prank(Y[a] + Y[b], P1) == prank(Y[a] + Y[b], P2)
        e = r - 4; assert e >= 0
        E += e
        if e:
            bad.append((a, b, e))
            if "0" in (a, b) or "F" in (a, b): endE += e
print("subspaces", len(Y), "covering edges", nedge, "E =", E, "excess on end edges", endE, "bad", bad)
