"""d-mcopy-code STEP 2e: value-free exact relay lower bound per group.
Pickup set at frame G = needed sources s whose content can reach (value-free, timing ignored, fixed chains,
no helper/source reads a target) some existing helper/source register r having G in its chain.
A new relay = nested chain of such frames inside one group target's U_t; min #relays >= min set cover of the 64
needed sources by such chains (content reaches group targets only through new relays: reach.py).
Usage: python3 -I cover_valuefree.py FULL-WORD.json.gz"""
import sys, json, functools
exec(open('chain_cover.py').read().split("gi = int(sys.argv[2])")[0])
Gs = json.load(open('mc_groups.json')); AG = json.load(open('all_groups.json'))
start, end, chain, paid, runs = frames(w)
prog = program(w)
full = basis([1 << i for i in range(9)])
out = []
for gi, A in enumerate(AG):
    T = A['targets']; mem = set(A['members']); allM = set(q for g in Gs for q in g['members'])
    rows = [{q: F(1)} for q in range(n)]; need = []; pk_of = {}
    for t, (ph, i, e) in enumerate(prog):
        if e[0] == 'add' and ph == 'B' and e[1] in mem and e[2] not in allM:
            xs = sorted(s for s in rows[e[2]] if s < v)
            if len(xs) == 8:
                j = len(set(pk_of.values()))
                for s in xs: pk_of[s] = j
        apply(e, rows)
    need = sorted(pk_of); idx = {s: k for k, s in enumerate(need)}
    ch = {}; at = defaultdict(list)
    for r in range(n):
        if r in mem: continue
        c = [tuple(U) for tag, U in chain[r] if U and tuple(U) != full]
        ch[r] = c
        for G in c: at[G].append(r)
    isT = lambda r: v <= r < 2*v
    # per-source reachability -> bitmask per node
    mask = defaultdict(int)
    for s in need:
        seen = set(); stack = [(s, k) for k in range(len(ch[s]))]; seen |= set(stack)
        while stack:
            r, k = stack.pop(); G = ch[r][k]
            for k2 in range(k + 1, len(ch[r])):
                if (r, k2) not in seen: seen.add((r, k2)); stack.append((r, k2))
            for r2 in at[G]:
                if r2 == r or (isT(r) and not isT(r2)): continue
                k2 = ch[r2].index(G)
                if (r2, k2) not in seen: seen.add((r2, k2)); stack.append((r2, k2))
        for (r, k) in seen:
            if not isT(r): mask[ch[r][k]] |= 1 << idx[s]
    Ut = {t: tuple(end[t]) for t in T}
    maxi = set()
    for t in T:
        FL = [G for G in mask if mask[G] and contained(G, Ut[t])]
        succ = {G: [H for H in FL if len(H) > len(G) and contained(G, H)] for G in FL}
        @functools.lru_cache(maxsize=None)
        def best(G):
            # all maximal masks of chains starting at G
            res = set()
            if not succ[G]: return frozenset([mask[G]])
            for H in succ[G]:
                for m in best(H): res.add(m | mask[G])
            return frozenset(res)
        for G in FL: maxi |= set(best(G))
    ml = sorted(maxi, key=lambda m: -bin(m).count('1')); mx = []
    for m in ml:
        if not any((m | x) == x for x in mx): mx.append(m)
    fullm = (1 << len(need)) - 1
    cov = defaultdict(list)
    for m in mx:
        for k in range(len(need)):
            if m >> k & 1: cov[k].append(m)
    def solve(limit):
        @functools.lru_cache(maxsize=None)
        def rec(c, left):
            if c == fullm: return True
            if left == 0: return False
            if max(bin(m & ~c).count('1') for m in mx) * left < bin(fullm & ~c).count('1'): return False
            k = ((~c & fullm) & -(~c & fullm)).bit_length() - 1
            return any(rec(c | m, left - 1) for m in cov[k])
        return rec(0, limit)
    kmin = next(k for k in range(1, 20) if solve(k))
    pks = [dict(Counter(pk_of[need[k]] for k in range(len(need)) if m >> k & 1)) for m in mx]
    print('group', gi, 'maximal value-free chain sets', len(mx), 'sizes', sorted(Counter(bin(m).count('1') for m in mx).items()), 'MIN RELAYS (value-free) =', kmin)
    out.append(dict(group=gi, maximal_sets=pks, min_relays=kmin))
json.dump(out, open('cover_valuefree.json', 'w'), indent=0)
