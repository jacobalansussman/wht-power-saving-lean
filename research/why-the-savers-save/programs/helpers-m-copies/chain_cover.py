"""d-mcopy-code STEP 2b: exact lower bound on relays for one M group (support relaxation, timing ignored).
Givers = every source/helper register other than the M copies, at every frame G of its FIXED chain with G inside U_t
for some group target t; the giver's content there = union over its window of its needed-source support.
A relay is a nested chain of giver frames inside one U_t (its final reader's hyperplane); it can pick up only
content of givers at frames of its chain (a gate at G needs G in the giver's chain: no new giver run; the relay
climbs to G).  Every needed source (the 64 sources of the group's 8 packets) must enter some relay
(targets cannot read givers at a common chain frame carrying needed content: checked below).
min #chains covering the 64 sources = a lower bound on the relays of any delivery with all other chains fixed.
Usage: python3 -I chain_cover.py FULL-WORD.json.gz GROUP"""
import sys, json, itertools
from fractions import Fraction as F
from collections import Counter, defaultdict
sys.path.insert(0, '.')
from wordlib import *
w = load(sys.argv[1]); v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2*v + R
gi = int(sys.argv[2])
A = json.load(open('all_groups.json'))[gi]
Gs = json.load(open('mc_groups.json')); allM = set(q for g in Gs for q in g['members'])
start, end, chain, paid, runs = frames(w)
T = A['targets']; Ut = {t: tuple(end[t]) for t in T}
prog = program(w)
tagtime = {(ph, i): t for t, (ph, i, e) in enumerate(prog)}
def ttime(tag):
    if tag is None: return -1
    if tag[0] == 'final': return 10**9
    return tagtime[tuple(tag)]
mem = set(A['members'])
T_CLEAN = min(t for t, (ph, i, e) in enumerate(prog) if ph == 'cleanup')
rows = [{q: F(1)} for q in range(n)]
need = {}; pk_of = {}
for t, (ph, i, e) in enumerate(prog):
    if e[0] == 'add' and ph == 'B' and e[1] in mem and e[2] not in allM:
        xs = sorted(s for s in rows[e[2]] if s < v)
        if len(xs) == 8:
            j = len(set(pk_of.values())); 
            for s in xs: pk_of[s] = j
    apply(e, rows)
needed = sorted(pk_of); idx = {s: k for k, s in enumerate(needed)}
print('group', gi, 'needed sources', len(needed), 'packets', len(set(pk_of.values())))
# windows of every non-M register at frames inside some U_t
wins = defaultdict(list)
for g in range(n):
    if g in allM or v <= g < 2*v: continue   # helpers may not read targets (audit: non-target reads target)
    ch = chain[g]
    for k, (tag, U) in enumerate(ch):
        U = tuple(U)
        if not U: continue
        if not any(contained(U, Ut[t]) for t in T): continue
        a = ttime(tag); b = ttime(ch[k+1][0]) if k + 1 < len(ch) else 10**9
        b = min(b, T_CLEAN)   # cleanup gates act at the full frame: a register is never at a lower frame there
        wins[g].append((a, b, U))
rows = [{q: F(1)} for q in range(n)]
sup = defaultdict(int)      # (g, U) -> bitmask of needed sources held at some time in window
def upd(q, t):
    for (a, b, U) in wins.get(q, []):
        if a <= t < b:
            m = 0
            for s in rows[q]:
                if s in idx: m |= 1 << idx[s]
            sup[(q, U)] |= m
for q in wins: upd(q, -1)
for t, (ph, i, e) in enumerate(prog):
    apply(e, rows)
    for q in (e[1:3] if e[0] == 'add' else e[1:2]): upd(q, t)
fr = defaultdict(int); frg = defaultdict(set)
for (g, U), m in sup.items():
    if m: fr[U] |= m; frg[U].add(g)
frames_list = sorted(fr, key=len)
print('giver frames with needed content:', len(frames_list), Counter(len(U) for U in frames_list))
yreads = []
# targets (non-M) as direct readers at a common chain frame?  (helpers may not read targets; targets may read helpers)
for t in T:
    tf = set(tuple(U) for tag, U in chain[t])
    for U in frames_list:
        if U in tf: yreads.append((t, len(U), bin(fr[U]).count('1')))
print('group targets sharing a chain frame with a giver of needed content:', yreads)
# best chains per final reader: DP over containment DAG (all frames inside U_t), maximal coverage sets
best_cover = set()
for t in T:
    FL = [U for U in frames_list if contained(U, Ut[t])]
    # enumerate all nested chains (one frame per dim, increasing) -> coverage masks; keep Pareto-maximal masks
    succ = {U: [V for V in FL if len(V) > len(U) and contained(U, V)] for U in FL}
    masks = set()
    def dfs(U, m):
        m |= fr[U]
        nxt = succ[U]
        if not nxt: masks.add(m); return
        masks.add(m)
        for V in nxt: dfs(V, m)
    for U in FL: dfs(U, 0)
    best_cover |= masks
# keep maximal
ml = sorted(best_cover, key=lambda m: -bin(m).count('1'))
maxi = []
for m in ml:
    if not any((m | x) == x for x in maxi): maxi.append(m)
print('distinct maximal chain coverage sets:', len(maxi), ' sizes', Counter(bin(m).count('1') for m in maxi))
# packet coverage per maximal chain
def pk(m): return Counter(pk_of[needed[k]] for k in range(len(needed)) if m >> k & 1)
print(' largest:', [dict(pk(m)) for m in maxi[:6]])
full = (1 << len(needed)) - 1
# exact min set cover by iterative deepening with bitmask branch on lowest uncovered element
import functools
sys.setrecursionlimit(10000)
cover_of = defaultdict(list)
for m in maxi:
    for k in range(len(needed)):
        if m >> k & 1: cover_of[k].append(m)
def solve(limit):
    @functools.lru_cache(maxsize=None)
    def rec(covered, left):
        if covered == full: return True
        if left == 0: return False
        # lower bound: max cover size
        rem = bin(full & ~covered).count('1')
        mx = max(bin(m & ~covered).count('1') for m in maxi)
        if mx * left < rem: return False
        k = (~covered & full & -(~covered & full)).bit_length() - 1
        for m in cover_of[k]:
            if rec(covered | m, left - 1): return True
        return False
    return rec(0, limit)
kmin = None
for k in range(1, 12):
    if solve(k): kmin = k; break
print('EXACT minimum number of nested giver chains covering all needed sources (support relaxation):', kmin)
json.dump(dict(group=gi, needed=len(needed), giver_frames=len(frames_list), maximal_chain_sets=len(maxi),
               max_chain_cover=max(bin(m).count('1') for m in maxi), min_chains=kmin, target_direct=yreads),
          open(f'chain_cover_g{gi}.json', 'w'), indent=1)
