"""d-mcopy-code STEP 2d: value-free reachability (any coefficients, any re-valuation keeping every chain).
Nodes (r, G) for each register r (all but the group's own 8 M copies) and each frame G of its fixed chain.
Content moves r:G -> r:G' (later, larger frame of the same chain) and r:G <-> r':G when both chains contain G
(a gate at a common frame), except that a non-target may not read a target (audit rule); targets may read anything.
Timing ignored.  Question: can any needed source of the group reach any group target?  Frames at the full space
are dropped (no group target ever reaches them).
Usage: python3 -I reach.py FULL-WORD.json.gz"""
import sys, json
exec(open('chain_cover.py').read().split("gi = int(sys.argv[2])")[0])
Gs = json.load(open('mc_groups.json')); AG = json.load(open('all_groups.json'))
start, end, chain, paid, runs = frames(w)
prog = program(w)
full = basis([1 << i for i in range(9)])
res = []
for gi, A in enumerate(AG):
    T = set(A['targets']); mem = set(A['members']); allM = set(q for g in Gs for q in g['members'])
    rows = [{q: F(1)} for q in range(n)]; need = set()
    for t, (ph, i, e) in enumerate(prog):
        if e[0] == 'add' and ph == 'B' and e[1] in mem and e[2] not in allM:
            xs = [s for s in rows[e[2]] if s < v]
            if len(xs) == 8: need |= set(xs)
        apply(e, rows)
    at = defaultdict(list)   # frame -> registers having it in chain
    ch = {}
    for r in range(n):
        if r in mem: continue
        c = [tuple(U) for tag, U in chain[r] if U and tuple(U) != full]
        ch[r] = c
        for G in c: at[G].append(r)
    isT = lambda r: v <= r < 2*v
    seen = set(); stack = []
    for s in need:
        for k, G in enumerate(ch[s]):
            if (s, k) not in seen: seen.add((s, k)); stack.append((s, k))
    hitT = set()
    cutr = defaultdict(set)
    for e in w['cut']: cutr[e[2]].add(e[1])
    while stack:
        r, k = stack.pop(); G = ch[r][k]
        if r in T: hitT.add(r)
        for y in cutr.get(r, ()):   # CUT edge: target (at frame ()) reads centre r
            for k2 in range(len(ch[y])):
                if (y, k2) not in seen: seen.add((y, k2)); stack.append((y, k2))
        for k2 in range(k + 1, len(ch[r])):
            if (r, k2) not in seen: seen.add((r, k2)); stack.append((r, k2))
        for r2 in at[G]:
            if r2 == r: continue
            if isT(r) and not isT(r2): continue      # a helper/source may not read a target
            k2 = ch[r2].index(G)
            if (r2, k2) not in seen: seen.add((r2, k2)); stack.append((r2, k2))
    res.append(dict(group=gi, reached_group_targets=sorted(hitT), nodes_reached=len(seen)))
    print('group', gi, 'needed sources', len(need), 'group targets reachable without the group M copies (value-free):', sorted(hitT), 'nodes', len(seen))
json.dump(res, open('reach_withcut.json', 'w'), indent=0)
