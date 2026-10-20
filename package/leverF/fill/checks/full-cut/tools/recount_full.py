# fc-gate: independent recount on the raw design-B FILL (no chronology, no widths): per bin (replica o, child c,
# port p) of the four rank-48 children (8, 10 = hub pair 3, kind A4; 12, 13 = hub pair 0, kind F1) count the hosted
# roles (rows naming (o, c, p)). A bin with no row is fully idle (remainder 48).  Prints counts per child and the
# full-cut profile derived from them (heads = idle 8/10 bins, split = idle 12/13 bins).
# usage: python3 -I -B recount_full.py FILL.json.gz PROFILE.json T
import json, gzip, sys
from collections import Counter
rows = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))['rows']
prof = {int(k): v for k, v in json.load(open(sys.argv[2])).items()} if True else None
T = int(sys.argv[3])
hosted = Counter((o, c, p) for t, h, o, c, p, z in rows)
idle = Counter(); part = Counter(); hostedby = Counter()
for o in range(T):
    for c in (8, 10, 12, 13):
        for p in range(112):
            n = hosted.get((o, c, p), 0)
            if n == 0: idle[c] += 1
            else: part[c] += 1; hostedby[(c, n)] += 1
heads = idle[8] + idle[10]; non = idle[12] + idle[13]
base = {r: n for r, n in prof.items() if r != 48}
full = dict(base); full[24] = full.get(24, 0) + 2 * non
out = dict(rows=len(rows), idle_by_child=dict(sorted(idle.items())), hosted_bins_by_child=dict(sorted(part.items())),
           roles_per_hosted_bin={'%d:%d' % k: v for k, v in sorted(hostedby.items())},
           heads=heads, other_idle48=non, profile=prof, profile48_eq_idle=prof.get(48) == heads + non,
           full_cut_children=dict(sorted(full.items())), W=896 * T, H=heads, kc=25)
print(json.dumps(out, indent=1))
