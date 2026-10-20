"""control: every plain multiset (W frames, no deletion) must pass S, D, C, Z, TR (they are realised programs)."""
import sys, itertools, random
from vtests import Tester
m = int(sys.argv[1]); W = int(sys.argv[2]); nsamp = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9
t = Tester(m); N = t.N; Mo = t.Mo
DM = {z: Mo.dmask(z) for z in range(1, N + 1)}
bucket = {}
for S in itertools.combinations(range(1, N + 1), m):
    v = 0
    for z in S: v ^= DM[z]
    bucket.setdefault(v, []).append(tuple(sorted(S)))
fr = {k: sorted(bucket[Mo.smask(Mo.target(k))]) for k in range(N)}
rng = random.Random(7); bad = 0; n = 0
allks = list(itertools.combinations_with_replacement(range(N), W))
rng.shuffle(allks)
for ks in allks[:nsamp]:
    ch = [rng.choice(fr[k]) for k in ks]
    ms = [z for F in ch for z in F]
    r = t.run(ms, ks); n += 1
    ok = r["S"] and r["D"] and r.get("C", True) and r.get("Z", True) and r.get("TR", True)
    if not ok: bad += 1; print("PLAIN FAILS", ks, r)
print(f"m={m} W={W}: {n} plain programs tested, {bad} fail any test")
