# o-model: per child-pair group statistics of a literal owner FILL (our code; reads data only).
# For each of the nine child-pair classes (window (lo, hi] from the template, as in the hub engines), count bins by
# remainder, and list literal heads (empty bins of children 8/10). Args: FILL ORDER DISPATCH PROVIDERS|- TEMPLATE [WIDTHS-JSON]
import json, gzip, sys
from collections import Counter, defaultdict
rd = lambda p: json.loads(gzip.decompress(open(p, 'rb').read())) if p.endswith('.gz') else json.load(open(p))
fill, order, disp, prov, tmpl = (rd(sys.argv[i]) if sys.argv[i] != '-' else None for i in range(1, 6))
widths = json.loads(open(sys.argv[6]).read()) if len(sys.argv) > 6 else {pid: v['widths'] for pid, v in prov.items()}
disp = {(s, r): p for s, r, p in disp}
T = 1 + max(r for s, r in order)
rank = {c['child']: c['rank'] for c in tmpl['idle_children']}
used = Counter()
for t, h, o, c, p, z in fill['rows']:
    used[o, c, p] += widths[disp[tuple(order[t])]][h]
pairs = [(0, 2), (1, 3), (4, 6), (5, 7), (8, 10), (9, 11), (12, 13), (14, 15), (16, 17)]
out = {}
for cs in pairs:
    rem = Counter()
    for o in range(T):
        for c in cs:
            for p in range(112):
                rem[rank[c] - used.get((o, c, p), 0)] += 1
    out['%d/%d (rank %d)' % (cs[0], cs[1], rank[cs[0]])] = dict(sorted(rem.items()))
heads = Counter(c for o in range(T) for c in (8, 10) for p in range(112) if not used.get((o, c, p)))
print(json.dumps(dict(T=T, remainders_by_pair_class=out, literal_heads=dict(heads), total_heads=sum(heads.values())), indent=1))
