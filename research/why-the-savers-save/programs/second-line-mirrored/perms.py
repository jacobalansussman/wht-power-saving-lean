#!/usr/bin/env python3
"""cx-mirror / perms.py (own code, stdlib).  Every coordinate permutation sigma of the 9-dim window as a twisted
mirror A = P_sigma (orthogonal): which pass u . A u = 1 on all 112 ports?  A non-identity passing sigma is asymmetric
unless sigma is an involution; for each passer, the longest drop-free x chain inside each total's support (exact for
symmetric N via maximum independent set; randomised insertion for asymmetric N)."""
import itertools, json, os, importlib.util, random
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('tw', os.path.join(here, 'twist_lib.py')); tw = importlib.util.module_from_spec(spec); spec.loader.exec_module(tw)
U, Sk, v, dot, apply, mis = tw.U, tw.Sk, tw.v, tw.dot, tw.apply, tw.mis
sup = [[i for i in range(9) if (u >> i) & 1] for u in U]
passers = []
for p in itertools.permutations(range(9)):
    good = True
    for s in sup:
        c = 0
        ss = set(s)
        for i in s:
            if p[i] in ss: c ^= 1
        if not c: good = False; break
    if good: passers.append(p)
print('permutations passing:', len(passers), '| involutions among them:', sum(1 for p in passers if all(p[p[i]] == i for i in range(9))))
rnd = random.Random(3)
def chain_len(N, S, tries=200):
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
for p in passers[:60]:
    A = [1 << p[i] for i in range(9)]
    N = [[0 if a == b else 1 ^ dot(U[b], apply(A, U[a])) for b in range(v)] for a in range(v)]
    sym = all(N[a][b] == N[b][a] for a in range(v) for b in range(v))
    if sym:
        adj = {a: {b for b in range(v) if N[a][b]} for a in range(v)}
        c = [len(mis(S, adj)) for S in Sk]
    else:
        c = [chain_len(N, S) for S in Sk]
    res.append(dict(perm=p, symmetric=sym, chains=c))
    print(p, 'symmetric N' if sym else 'asymmetric N', 'chains', c)
json.dump(dict(passers=len(passers), results=res), open(os.path.join(here, 'perms.json'), 'w'), indent=1)
