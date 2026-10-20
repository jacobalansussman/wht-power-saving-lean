# lever-reprice: independent recount of fully unhosted rank-48 bins per child on a staged544 FILL (no chronology,
# no widths: a bin is fully idle iff no role at all names (owner, child, port)). Run: python3 -I -B tools/recount.py FILL T
import json, gzip, sys
from collections import Counter
rows = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))['rows']; T = int(sys.argv[2])
hosted = {(o, c, p) for t, h, o, c, p, z in rows}
idle = Counter(c for o in range(T) for c in (8, 10, 12, 13) for p in range(112) if (o, c, p) not in hosted)
print('fully idle rank-48 bins by child:', dict(sorted(idle.items())), 'pair3 (8,10):', idle[8] + idle[10], 'pair0 (12,13):', idle[12] + idle[13])
