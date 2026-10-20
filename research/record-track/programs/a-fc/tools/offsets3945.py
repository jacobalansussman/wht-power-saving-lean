# copied unchanged from wht45/reopt/o-model/probe/offsets3945.py (o-model)
# o-model probe: #3945's owner/build_fill.py step (keep every (t, h, o, c, p); recompute compact chronological offsets with the
# final provider widths), re-implemented on any staged544 T360 fill. Our code. Args: PKG3945_OWNER_DIR IN_FILL OUT_DIR
import json, gzip, sys, hashlib
from collections import Counter
from pathlib import Path
R, fin, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
rd = lambda p: json.loads(gzip.decompress(p.read_bytes()) if p.suffix == '.gz' else p.read_bytes())
def rank(rows):
    piv = {}
    for x in rows:
        while x:
            b = x.bit_length()
            if b in piv: x ^= piv[b]
            else: piv[b] = x; break
    return len(piv)
I = R / 'inputs'; order = rd(I / 'ORDER.json'); dispatch = {(s, r): pid for s, r, pid in rd(I / 'DISPATCH.json')}
fill = rd(fin); widths = {}; profiles = {}; full = [1 << i for i in range(9)]
for tag, pid in [('base', 'gauge44'), ('extra', 'extra')]:
    word = rd(I / f'WORD-{tag}.json.gz'); origin = rd(I / f'old-NATIVE-{tag}.json.gz')['helper_origin']
    widths[pid] = [rank(word['helper_final'].get(str(r), full)) - rank(word['start_extra'].get(str(r), [])) for r in origin]
    profiles[pid] = {int(r): n for r, n in word['report']['local_profile'].items()}
caps = {c['child']: c['rank'] for c in rd(I / 'TEMPLATE.json.gz')['idle_children']}; used = {}
bad = 0
for j, row in enumerate(fill['rows']):
    t, h, o, c, p, z = row; assert (t, h) == divmod(j, 408)
    key = (o, c, p); prefix = used.get(key, 0); w = widths[dispatch[tuple(order[t])]][h]
    if prefix + w > caps[c]: bad += 1
    row[5] = prefix; used[key] = prefix + w
assert not bad, 'owner capacity exceeded in %d roles' % bad
profile = Counter()
for o in range(360):
    for c, cap in caps.items():
        for p in range(112):
            idle = cap - used.get((o, c, p), 0)
            if idle: profile[idle] += 1
for call in order: profile.update(profiles[dispatch[tuple(call)]])
out.mkdir(parents=True, exist_ok=False)
payload = json.dumps(fill, separators=(',', ':')).encode(); (out / 'FILL.json.gz').write_bytes(gzip.compress(payload, compresslevel=6, mtime=0))
(out / 'PROFILE.json').write_text(json.dumps(dict(sorted(profile.items())), indent=2))
print(json.dumps({'status': 'OFFSETS_REGENERATED', 'roles': len(fill['rows']), 'payload_sha256': hashlib.sha256(payload).hexdigest(), 'profile': dict(sorted(profile.items()))}))
