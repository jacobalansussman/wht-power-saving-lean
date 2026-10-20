"""openunion.py -- w5-twoshared (resumed): decide MANY Lemma-D injection-set pairs of an open private core at once.
Every coefficient is free, so a solution of a skeleton S extends by zeros to a solution of any skeleton S' that contains
S's summands (closure argument: S's closure is closed in S' once the extra coefficients are 0, and the extra rows then get
zero contributions).  Hence INCONSISTENT for the UNION of several injection sets decides every sub-collection, in
particular every pair (one set per top) inside it.
Unions tried, smallest first: for each pair of type combinations (c0 of top 0, c1 of top 1), the union of the first n
sets of c0 and the first n sets of c1, n = 1, 2, 4, 8, ..., up to all; then everything.  A union is decided only if both
closures fit in maxcol columns.
Usage: openunion.py seed core cap maxcol timeout"""
import sys, os, itertools, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
seed, cidx, cap, maxcol, tmo = [int(a) for a in sys.argv[1:6]]
sys.argv = ['opencase.py', str(seed), str(cidx), '1', str(cap)]
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'opencase.py')).read().split("# (i) each output alone")[0]
exec(compile(src, 'opencase_build', 'exec'))
from exact import Closure
os.makedirs('open_sing', exist_ok=True)
t0 = time.time()

def key(item):
    mi, s_, M = item
    return (mi, s_, G.intern(M))

by_combo = [{}, {}]
for e in (0, 1):
    for (tt, new) in inj_by_top[e]:
        by_combo[e].setdefault(tt, []).append(new)
for e in (0, 1):
    print(f"top {e}: combos " + ", ".join(f"{tt}: {len(v)} sets" for tt, v in by_combo[e].items()), flush=True)

def decide_union(sets, tag):
    seen = {}; extra = []
    for new in sets:
        for it in new:
            k = key(it)
            if k not in seen: seen[k] = 1; extra.append(it)
    mvs = moves_of(R0, R1, None, None, extra)
    P = build(2, mvs)
    Cs = [Closure(P, l, cap) for l in range(2)]
    if not all(C.ok for C in Cs):
        return 'CAP', len(extra), None
    sizes = [(C.nr, C.nc) for C in Cs]
    if max(C.nc for C in Cs) > maxcol:
        return 'BIG', len(extra), sizes
    scr, _ = singular_script(P, cap=cap)
    v, out = decide(scr, f'open_sing/U_s{seed}_c{cidx}.sing', timeout=tmo)
    if v == 'CONSISTENT':
        pickle.dump({'moves': mvs}, open(f'open_sing/U_CONS_s{seed}_c{cidx}_{tag}.pkl', 'wb')); print(out[-1500:], flush=True)
    return v, len(extra), sizes

tally = {}
covered_pairs = 0; total_pairs = len(inj_by_top[0]) * len(inj_by_top[1])
for c0, c1 in itertools.product(by_combo[0], by_combo[1]):
    A, B = by_combo[0][c0], by_combo[1][c1]
    n = 1; best = 0
    while True:
        nA, nB = min(n, len(A)), min(n, len(B))
        v, ne, sizes = decide_union(A[:nA] + B[:nB], f'{c0}_{c1}_{n}')
        print(f"  combos {c0} / {c1}: first {nA} x {nB} sets ({ne} distinct injections): closures {sizes} {v} ({time.time()-t0:.0f}s)", flush=True)
        tally[v] = tally.get(v, 0) + 1
        if v == 'INCONSISTENT': best = nA * nB
        if v != 'INCONSISTENT' or (nA == len(A) and nB == len(B)): break
        n *= 2
    covered_pairs += best
    print(f" combos {c0} / {c1}: {best} of {len(A) * len(B)} pairs decided INCONSISTENT by unions", flush=True)
print(f"COVERED {covered_pairs} of {total_pairs} pairs (one set per top); tally {tally} ({time.time()-t0:.0f}s)", flush=True)
