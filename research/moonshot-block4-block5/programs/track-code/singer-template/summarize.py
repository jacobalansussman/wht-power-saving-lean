"""summarize.py (mw2-singer4): totals over fano_gate outputs (final line, or the last progress line of an unfinished run).
usage: python3 summarize.py file...   -> per (kind, W, k): prefixes built, top-augmentation rank histogram, Singer passes,
whether r <= nK1 held in every prefix that recorded nK1."""
import sys, ast, re, os
from collections import defaultdict
tot = defaultdict(lambda: {'built': 0, 'aug': defaultdict(int), 'singer': 0, 'cap': 0, 'files': 0, 'viol': 0, 'nk': 0, 'partial': 0})
for f in sys.argv[1:]:
    m = re.search(r'out_fg_(\w+?)_(\d+)_(\d+)_s(\d+)\.txt', os.path.basename(f))
    if not m:
        continue
    kind, W, k = m.group(1), int(m.group(2)), int(m.group(3))
    lines = [l for l in open(f) if l.startswith('fano_gate ') or l.startswith('progress ')]
    if not lines:
        continue
    last = lines[-1]
    d = ast.literal_eval(last[last.index('{'):last.rindex('}') + 1])
    t = tot[(kind, W, k)]
    t['files'] += 1; t['built'] += d['built']; t['singer'] += d['singer_pass']; t['cap'] += d.get('full_cap', 0)
    t['partial'] += last.startswith('progress')
    for r, c in d['augrank'].items():
        t['aug'][int(r)] += c
    for key, c in d.get('nK1_rank', {}).items():
        a, b = map(int, key.split(':')); t['nk'] += c
        if b > a:
            t['viol'] += c
for (kind, W, k), t in sorted(tot.items(), key=lambda x: (x[0][1], x[0][2], x[0][0])):
    M = k + W; plain = 4 * W
    print(f'{kind} {W} in {M} (k={k}{", plain" if M == plain else ""}): built {t["built"]} in {t["files"]} runs'
          f'{" (" + str(t["partial"]) + " unfinished)" if t["partial"] else ""}; top-aug rank {dict(sorted(t["aug"].items()))}; '
          f'max {max(t["aug"]) if t["aug"] else None} (need {W}); Singer passes {t["singer"]}; '
          f'r <= nK1 checked on {t["nk"]}, violations {t["viol"]}')
