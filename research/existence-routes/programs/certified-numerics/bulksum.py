# bulksum.py -- summary of out/bulk_*.jsonl, de-duplicated by equation content
import json, glob, hashlib, collections, sys
recs = {}
for f in sorted(glob.glob('out/bulk_*.jsonl')) + sorted(glob.glob('out/bulk2_*.jsonl')):   # second pass overrides
    for l in open(f):
        r = json.loads(l)
        if r['file'] in recs and not r['nulla'] and recs[r['file']]['nulla']: continue
        recs[r['file']] = r
seen = {}; c = collections.Counter(); bad = []; undec = []
for fn, r in sorted(recs.items()):
    h = hashlib.md5(json.dumps(json.load(open(fn))['eqs']).encode()).hexdigest()
    if h in seen: continue
    seen[h] = fn
    key = (r['singular'], r['nulla'] or 'none found')
    c[key] += 1
    if r['nulla'] and r['singular'] != 'infeasible': bad.append(fn)
    if not r['nulla'] and r['singular'] == 'infeasible': undec.append((fn, r['n'], r['N'], [t[3] for t in r.get('tries', [])][-1:]))
print('records', len(recs), 'distinct systems', len(seen))
for k, v in sorted(c.items(), key=lambda t: -t[1]): print('%4d  Singular: %-16s  certnum: %s' % (v, k[0], k[1]))
print('refused by certnum but feasible for Singular:', bad)
reasons = collections.Counter()
for fn, n, N, last in undec:
    r = recs[fn]; msgs = [t[3] for t in r.get('tries', [])]
    reasons['complex coefficients (first pass only)' if any('complex' in m for m in msgs) else
            'size cap' if any('too large' in m for m in msgs) else 'needs higher degree'] += 1
print('infeasible for Singular, no certificate at degree <= 2:', dict(reasons))
