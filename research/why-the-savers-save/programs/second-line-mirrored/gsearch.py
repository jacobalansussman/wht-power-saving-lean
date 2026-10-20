#!/usr/bin/env python3
"""cx-mirror / gsearch.py (own code, stdlib).  Second lines beyond the twisted mirrors?  Write g_b = j + h_b: the
conditions (g_b even, u_b . g_b = 0, g_a . g_b = 1 on every junk pair) become h_b odd, u_b . h_b = 1, and
h_a . h_b = 0 whenever u_a . u_b = 0 (a != b).  The twisted mirrors are h_b = A u_b.  A min-conflicts local search
from random starts looks for solutions that are NOT of that form (h not the restriction of an orthogonal map), and
reports, for each solution found, whether it is linear, and the longest drop-free x chain inside each support
(randomised insertion; N[a][b] = g_a . u_b may be asymmetric here)."""
import os, json, random, time, importlib.util
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('tw', os.path.join(here, 'twist_lib.py')); tw = importlib.util.module_from_spec(spec); spec.loader.exec_module(tw)
U, Sk, v, dot = tw.U, tw.Sk, tw.v, tw.dot
J = (1 << 9) - 1
junk = {a: [b for b in range(v) if b != a and dot(U[a], U[b]) == 0] for a in range(v)}
cand = {b: [h for h in range(512) if bin(h).count('1') % 2 == 1 and dot(U[b], h) == 1] for b in range(v)}
rnd = random.Random(5)
def conflicts(H, a, h):
    return sum(1 for b in junk[a] if dot(h, H[b]) == 1)
def total(H): return sum(conflicts(H, a, H[a]) for a in range(v)) // 2
def linear(H):
    """is b -> h_b the restriction of a linear map? solve on a basis of ports, test the rest"""
    basis = []; imgs = []
    for b in range(v):
        x, y = U[b], H[b]
        for (bx, by) in zip(basis, imgs):
            if (x >> (bx.bit_length() - 1)) & 1: x ^= bx; y ^= by
        if x:
            basis.append(x); imgs.append(y)
            # keep reduced
        elif y: return False
    return True
sols = []; t0 = time.time(); starts = 0
while time.time() - t0 < 90 and len(sols) < 12:
    starts += 1
    H = [rnd.choice(cand[b]) for b in range(v)]
    cur = total(H)
    for it in range(4000):
        if cur == 0: break
        bad = [a for a in range(v) if conflicts(H, a, H[a])]
        a = rnd.choice(bad)
        best = min(cand[a], key=lambda h: (conflicts(H, a, h), rnd.random()))
        cur += conflicts(H, a, best) - conflicts(H, a, H[a]); H[a] = best
    if cur == 0:
        key = tuple(H)
        if key not in {tuple(s) for s in sols}: sols.append(H)
print('random starts', starts, '| solutions found', len(sols))
def chain_len(N, S, tries=150):
    best = 0
    for _ in range(tries):
        order = list(S); rnd.shuffle(order); ch = []
        for x in order:
            for pos in range(len(ch) + 1):
                if all(not N[x][a] for a in ch[:pos]) and all(not N[b][x] for b in ch[pos:]):
                    ch.insert(pos, x); break
        best = max(best, len(ch))
    return best
res = []
for H in sols:
    g = [J ^ h for h in H]
    N = [[0 if a == b else dot(g[a], U[b]) for b in range(v)] for a in range(v)]
    lin = linear(H); sym = all(N[a][b] == N[b][a] for a in range(v) for b in range(v))
    c = [chain_len(N, S) for S in Sk]
    res.append(dict(linear=lin, symmetric=sym, distinct_lines=len(set(H)), chains=c))
    print('solution: linear', lin, '| N symmetric', sym, '| distinct lines', len(set(H)), '| chains per support', c)
json.dump(dict(starts=starts, found=len(sols), results=res), open(os.path.join(here, 'gsearch.json'), 'w'), indent=1)
