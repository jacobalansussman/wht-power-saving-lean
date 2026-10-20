import json, sys, glob, collections
for fn in sys.argv[1:]:
    rs = [json.loads(l) for l in open(fn) if l.strip()]
    c = sorted(r['cost'] for r in rs if r.get('cost') is not None)
    nd = sum(1 for r in rs if r.get('cost') is None)
    hist = collections.Counter(round(x, 1) for x in c)
    print('%s: %d trials (%d without a deletion-passing sequence); lowest %s; most common %s' % (fn, len(rs), nd, [round(x, 4) for x in c[:6]], hist.most_common(4)))
