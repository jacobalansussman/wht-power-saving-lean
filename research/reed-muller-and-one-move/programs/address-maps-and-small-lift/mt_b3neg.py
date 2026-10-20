"""find block-3 multi-target instances that are impossible (below the engine hook 2.8 moves per array: 3 in 8,
4 in 11) yet pass the exact tests: these are the NEGATIVE controls for the LM probe."""
import itertools, sys
from mt_lm import Setup
from mt_cover import cover_test
W, M = int(sys.argv[1]), int(sys.argv[2])
found = []
targets = [",".join(map(str, t)) for t in itertools.combinations_with_replacement(range(7), W) if t[0] == 0]
for comp in itertools.combinations_with_replacement(range(1, 8), M):
    if len(set(comp)) < 7: continue          # use every direction at least once (else projection-type failures)
    for ts in targets:
        st = Setup(3, "k:" + ts, "list:" + ",".join(map(str, comp)), "lex")
        et = st.exact_tests()
        if et["excluded"]: continue
        ct = cover_test(st)
        if ct["excluded"]: continue
        found.append((ts, comp, et["tight"]))
        if len(found) >= 6: break
    if len(found) >= 6: break
for f in found: print("passes exact tests:", f)
print("found", len(found))
