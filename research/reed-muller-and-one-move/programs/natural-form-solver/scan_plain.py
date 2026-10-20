"""scan_plain.py W: canonical target sets {0 = k_1 <= k_2 <= ...} of size W at block 5, up to the affine maps
k -> 2^s k + c; for each, the plain frames' exact tests and how many single deletions pass span + deletion."""
import sys, itertools, warnings; warnings.filterwarnings("ignore")
from inst import canonical, Instance
W = int(sys.argv[1]); lim = int(sys.argv[2]) if len(sys.argv) > 2 else 10**9
seen = set(); rows = []
def canon(ks):
    best = None
    for s in range(5):
        for c in range(31):
            t = tuple(sorted(((k * 2 ** s) + c) % 31 for k in ks))
            if best is None or t < best: best = t
    return best
for rest in itertools.combinations_with_replacement(range(31), W - 1):
    ks = canon((0,) + rest)
    if ks in seen: continue
    seen.add(ks)
    inst, fd = canonical(5, list(ks), dirs="frames", order="lex")
    et = inst.exact_tests(); d = inst.dirs
    npass = 0; ntight0 = 0
    for i in range(len(d)):
        e2 = Instance(5, inst.T, d[:i] + d[i + 1:]).exact_tests()
        if e2["span_ok"] and e2["violations"] == 0:
            npass += 1; ntight0 += e2["tight"] == 0
    rows.append((ks, et["tight"], len(set(d)), npass, ntight0))
    print(ks, "tight", et["tight"], "distinct", len(set(d)), "deletions passing", npass, "of which tight-free", ntight0, flush=True)
    if len(seen) >= lim: break
print("classes", len(seen), "with a passing deletion", sum(r[3] > 0 for r in rows))
