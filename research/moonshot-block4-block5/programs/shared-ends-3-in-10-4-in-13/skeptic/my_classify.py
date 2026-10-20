"""v-o3-shared: independent recount of rationed / minimal (W+1)-piece configurations over abelian R.
Each Galois orbit O gets S(O) = support of the unique relation (nonempty subset of rows 0..W).
minimal  : intersection of all S(O) empty
rationed : every pair {i,j} contains some S(O)
classes  : orbits of S_{W+1} (rows) x Aut(R) (acting on the orbit labels)."""
import itertools, sys
def subsets(n):
    return [frozenset(c) for k in range(1, n+1) for c in itertools.combinations(range(n), k)]
GROUPS = {  # name: (number of orbits, list of permutations of orbit labels induced by Aut(R))
 'Z2': (2, [(0,1)]),
 'Z3': (2, [(0,1)]),
 'Z4': (3, [(0,1,2)]),
 'Z2^2': (4, [p for p in itertools.permutations(range(4)) if p[0]==0]),  # S3 on nontrivial chars
 'Z6': (4, [(0,1,2,3)]),
 'Z8': (4, [(0,1,2,3)]),
}
def run(W, name, mode):
    L = W+1; no, auts = GROUPS[name]
    S = subsets(L)
    pairs = [frozenset(p) for p in itertools.combinations(range(L), 2)]
    good = []
    for a in itertools.product(S, repeat=no):
        if frozenset.intersection(*a):
            continue
        if mode == 'rat' and not all(any(s <= p for s in a) for p in pairs):
            continue
        good.append(a)
    seen = set(); ncls = 0
    rowperms = list(itertools.permutations(range(L)))
    for a in good:
        if a in seen: continue
        ncls += 1
        for rp in rowperms:
            b = tuple(frozenset(rp[i] for i in s) for s in a)
            for ap in auts:
                c = tuple(b[ap[k]] for k in range(no))
                seen.add(c)
    return len(good), ncls
for W in (3, 4):
    for name in GROUPS:
        if W == 4 and name in ('Z8',): pass
        out = []
        for mode in ('min', 'rat'):
            if W == 4 and mode == 'min' and GROUPS[name][0] >= 4: out.append(('min', 'skip')); continue
            out.append((mode,) + run(W, name, mode))
        print('W', W, name, out, flush=True)
