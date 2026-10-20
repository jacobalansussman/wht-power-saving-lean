#!/usr/bin/env python3
"""cx-mirror / twist2.py (own code, stdlib).  Asymmetric twisted mirrors: A in O(9, F_2) built from random products
of Eichler maps and coordinate permutations, kept when u_b . A u_b = 1 for all 112 ports and A is NOT symmetric
(so N is asymmetric and a drop-free x chain is an acyclic set of the digraph a -> b iff N[a][b] = 1, ordered so that
every arc goes from an earlier to a later member).  For each kept A: the longest such chain found inside each
total's support by randomised insertion (a lower bound on the true maximum; the proved upper bound is 45)."""
import json, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
spec = importlib.util.spec_from_file_location('tw', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'twist_lib.py'))
tw = importlib.util.module_from_spec(spec); spec.loader.exec_module(tw)
U, Sk, v, dot, apply, compose, eichler, orth, ok, I9 = tw.U, tw.Sk, tw.v, tw.dot, tw.apply, tw.compose, tw.eichler, tw.orth, tw.ok, tw.I9
rnd = random.Random(11)
even = [w for w in range(1, 512) if bin(w).count('1') % 2 == 0]
def perm(p):
    return [1 << p[i] for i in range(9)]
def chain_len(N, S, tries=300):
    best = 0
    for _ in range(tries):
        order = list(S); rnd.shuffle(order); ch = []
        for x in order:
            # insert x at a position: no arc from x to an earlier member, no arc from a later member to x
            for pos in range(len(ch) + 1):
                if all(not N[x][a] for a in ch[:pos]) and all(not N[b][x] for b in ch[pos:]):
                    ch.insert(pos, x); break
        best = max(best, len(ch))
    return best
kept = []; tried = 0
while tried < 200000 and len(kept) < 25:
    tried += 1
    A = I9
    for _ in range(rnd.randint(2, 6)):
        if rnd.random() < 0.5:
            w = rnd.choice(even); z = rnd.choice([z for z in even if dot(w, z) == 0 and z != w])
            A = compose(eichler(w, z), A)
        else:
            p = list(range(9)); rnd.shuffle(p); A = compose(perm(p), A)
    sym = all(dot(A[i], 1 << j) == dot(A[j], 1 << i) for i in range(9) for j in range(9))
    if not sym and orth(A) and ok(A) and A not in kept:
        kept.append(A)
print('random products tried', tried, '| asymmetric A kept', len(kept))
res = []
for A in kept:
    N = [[0 if a == b else 1 ^ dot(U[b], apply(A, U[a])) for b in range(v)] for a in range(v)]
    c = [chain_len(N, S) for S in Sk]
    res.append(dict(A=A, chains=c))
    print('A', A, 'chains per support', c)
out = dict(tried=tried, kept=len(kept), results=res,
           best_worst_support=max((min(r['chains']) for r in res), default=None),
           best_any_support=max((max(r['chains']) for r in res), default=None))
print('best over A of the worst support:', out['best_worst_support'], '| best chain anywhere:', out['best_any_support'])
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'twist2.json'), 'w'), indent=1)
