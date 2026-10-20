"""block-m W-array 'plain minus one move' scan; S+D screened by the fast form, then every S+D pass re-checked directly
and run through C, Z, TR (vtests.py)."""
import sys, itertools, time
import numpy as np
from vtests import Tester
m = int(sys.argv[1]); W = int(sys.argv[2]); lim = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9
t = Tester(m); N = t.N; Mo = t.Mo
DM = {z: Mo.dmask(z) for z in range(1, N + 1)}
bucket = {}
for S in itertools.combinations(range(1, N + 1), m):
    v = 0
    for z in S: v ^= DM[z]
    bucket.setdefault(v, []).append(tuple(sorted(S)))
fr = {k: sorted(bucket[Mo.smask(Mo.target(k))]) for k in range(N)}
E = t.E.astype(np.int16)
def canon(ks):
    return min(tuple(sorted(((k * 2 ** s) + c) % N for k in ks)) for s in range(m) for c in range(N))
seen = set(); nplain = 0; out = []; t0 = time.time()
for rest in itertools.combinations_with_replacement(range(N), W - 1):
    ks = canon((0,) + rest)
    if ks in seen: continue
    seen.add(ks)
    n = sum(t.parity(k).astype(np.int16) for k in ks)
    done = set()
    for ch in itertools.product(*[fr[k] for k in ks]):
        c = np.zeros(N, dtype=np.int16)
        for F in ch:
            for z in F: c[z - 1] += 1
        key = c.tobytes()
        if key in done: continue
        done.add(key); nplain += 1
        D = E @ c
        if (D < n).any(): continue
        tight = (D == n) & (n > 0)
        bad = E[tight].any(axis=0)
        for z in range(1, N + 1):
            if c[z - 1] > 0 and not bad[z - 1]:
                ms = [y for y in range(1, N + 1) for _ in range(int(c[y - 1]))]
                ms.remove(z)
                res = t.run(ms, ks)
                assert res["S"] and res["D"]
                out.append((ks, z, int(c[z - 1]), res))
                print("PASS-SD", ks, "delete", z, "mult", int(c[z - 1]),
                      {k: v for k, v in res.items() if k not in ("S", "D")}, flush=True)
    if len(seen) >= lim: break
print(f"m={m} W={W}: classes {len(seen)}, distinct plain multisets {nplain}, S+D-passing single deletions {len(out)}; "
      f"passing C: {sum(o[3]['C'] for o in out)}, C+Z: {sum(o[3]['C'] and o[3]['Z'] for o in out)}, "
      f"C+Z+TR: {sum(o[3]['C'] and o[3]['Z'] and o[3]['TR'] for o in out)}  ({time.time()-t0:.0f}s)", flush=True)
