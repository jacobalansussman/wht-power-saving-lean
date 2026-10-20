import sys, itertools, time
import numpy as np
from vtests import Tester
m = int(sys.argv[1]); W = int(sys.argv[2])
t = Tester(m); N = t.N; Mo = t.Mo
DM = {z: Mo.dmask(z) for z in range(1, N + 1)}
bucket = {}
for S in itertools.combinations(range(1, N + 1), m):
    v = 0
    for z in S: v ^= DM[z]
    bucket.setdefault(v, []).append(tuple(sorted(S)))
fr = {k: sorted(bucket[Mo.smask(Mo.target(k))]) for k in range(N)}
def canon(ks):
    return min(tuple(sorted(((k * 2 ** s) + c) % N for k in ks)) for s in range(m) for c in range(N))
seen = set(); nplain = 0; npassSD = 0; out = []
t0 = time.time()
for rest in itertools.combinations_with_replacement(range(N), W - 1):
    ks = canon((0,) + rest)
    if ks in seen: continue
    seen.add(ks)
    done = set()
    for ch in itertools.product(*[fr[k] for k in ks]):
        ms = tuple(sorted(z for F in ch for z in F))
        if ms in done: continue
        done.add(ms); nplain += 1
        for z in sorted(set(ms)):
            dl = list(ms); dl.remove(z)
            res = t.run(dl, ks)
            if res["S"] and res["D"]:
                npassSD += 1
                out.append((ks, ms, z, res))
                print("PASS-SD", ks, "delete", z, "mult", ms.count(z), {k: v for k, v in res.items() if k not in ("S", "D")}, flush=True)
print(f"m={m} W={W}: classes {len(seen)}, distinct plain multisets {nplain}, S+D-passing single deletions {npassSD}; "
      f"also passing C: {sum(o[3]['C'] for o in out)}, C+Z: {sum(o[3]['C'] and o[3]['Z'] for o in out)}, "
      f"C+Z+TR: {sum(o[3]['C'] and o[3]['Z'] and o[3]['TR'] for o in out)}  ({time.time()-t0:.0f}s)", flush=True)
