"""vshape.py: check the shape the track states for the 5-move certificates, on MY local certificates: the tight
supports inside U are exactly the three pairs of Y, each with 2 moves and 2 odd arrays; F acts on exactly U with
r_F = 0; the two moves of U outside Y are two copies of one direction (35 sets) or two different directions (9)."""
import json
from collections import Counter
C = json.load(open("vcerts_local.json")); stats = Counter()
for c in C:
    Y = c["Y"]; U = c["U"]; sups = {tuple(D): p for D, p in c["sups"]}
    tight = [D for D, p in sups.items() if len(D) == sum(p) and sum(p) >= 1]
    pairs = sorted(tuple(sorted((a, b))) for i, a in enumerate(Y) for b in Y[i + 1:])
    stats["tight inside U == the 3 pairs of Y"] += sorted(tuple(sorted(t)) for t in tight) == pairs
    stats["each tight pair has 2 odd arrays"] += all(sum(sups[t]) == 2 for t in tight)
    stats["full U is a support with r = 0"] += (tuple(U) in sups and sum(sups[tuple(U)]) == 0)
    rest = [c["moves"][i] for i in U if i not in Y]
    stats["outside-Y pair is one repeated direction"] += (rest[0] == rest[1])
    stats["Y has 3 different directions"] += len({c["moves"][i] for i in Y}) == 3
print(len(C), dict(stats))
