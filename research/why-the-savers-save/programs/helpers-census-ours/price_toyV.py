"""Price toy V's counts on #942 (what-ifs): remove 112 factory/face helpers (two copies: 392 -> 280; M copies kept).
Usage: python3 -I price_toyV.py census942.json"""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import Counter
from census942 import price_delta, pct
C = json.load(open(sys.argv[1])); H = {int(q): v for q, v in C['registers'].items() if v['cls'] in ('factory', 'face')}
per = {}
for q in H:
    d = Counter()
    for r in H[q]['runs']: d[r] -= 1
    per[q] = pct(price_delta(dict(d)))
for name, qs in (('least valuable 112', sorted(H, key=lambda q: per[q])[:112]), ('most valuable 112', sorted(H, key=lambda q: -per[q])[:112]),
                 ('least valuable 104 (hub-compatible outputs: 392 -> 288)', sorted(H, key=lambda q: per[q])[:104])):
    d = Counter()
    for q in qs:
        for r in H[q]['runs']: d[r] -= 1
    b = price_delta(dict(d), certify=True); print(name, f'{float(b):.7e}', f'{pct(b):+.3f}%')
