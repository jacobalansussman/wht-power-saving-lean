"""Independent audit of a hosting assignment written by host_price.py (no code shared except the word file):
recomputes windows from the stage table, checks one owner per helper, <= one residual per role per call, every hosted
residual inside its child's window and outside the invocation's own pair/replica, sum of widths <= rank, then rebuilds the
global histogram from the remainders + the recorded bank pads and checks deficit 2560 and the recorded histogram.
usage: python3 -I audit_slots.py WORDS.json NAME RESULT.json"""
import sys, json, gzip
from collections import Counter, defaultdict
words, name, resf = sys.argv[1], sys.argv[2], sys.argv[3]
w = json.load(open(words))[name]; res = json.load(open(resf))
copies, defer = res['copies'], res['defer']
slots = json.load(gzip.open(resf + '.slots.gz', 'rt'))
PAIR = [1, 1, 1, 2, 2, 3, 3, 4, 4]
CH = {'A2': (16, 2, 'A', 1), 'B2': (8, 2, 'B', 2), 'A3': (32, 3, 'A', 3), 'B3': (8, 3, 'B', 4), 'A4': (48, 4, 'A', 5),
      'B4': (8, 4, 'B', 6), 'F1': (48, 1, 'F', 2), 'F2': (32, 2, 'F', 4), 'F3': (16, 3, 'F', 6)}
def inwin(rp, ch, k):
    r, pair, kind, pos = CH[ch]; run, st = divmod(k, 9)
    if kind == 'A':
        return run < rp or (run == rp and st <= pos)
    if kind == 'B':
        return run == rp and st == pos
    return (run == rp and st > pos) or (defer and run > rp)
wid = Counter({int(a): b for a, b in w['widths'].items()})
prof = Counter({int(a): b for a, b in w['profile'].items()})
if copies:
    wid[9] += 112; prof[1] += 112; prof[8] += 112
owner = set(); rolecall = set(); load = Counter()
for k, wd, h, (rp, ch), inst, off in slots:
    assert 0 <= h < wid[wd], 'helper index'
    assert (k, wd, h) not in owner, 'two owners'; owner.add((k, wd, h))
    pair = CH[ch][1]
    assert (rp, pair, inst, k) not in rolecall, 'two residuals of one call on one role'; rolecall.add((rp, pair, inst, k))
    assert inwin(rp, ch, k), 'outside window'
    assert not (k // 9 == rp and PAIR[k % 9] == pair), 'own pair'
    load[(rp, ch, inst)] += wd
hist = Counter({r: 72 * n for r, n in prof.items() if n})
for rp in range(8):
    for ch, (r, *_ ) in CH.items():
        for inst in range(224):
            rem = r - load[(rp, ch, inst)]
            assert rem >= 0, 'over rank'
            if rem:
                hist[rem] += 1
pads = Counter({int(a): b for a, b in res['pads'].items()})
for p, n in pads.items():
    if p:
        hist[p] += n
nb = sum(pads.values()); W = 7168 + nb
assert len(owner) + res['unhosted'] == 72 * sum(wid.values())
hosted_w = sum(s_[1] for s_ in slots); total_w = 72 * sum(a * b for a, b in wid.items())
assert sum((73 - p) * n for p, n in pads.items()) == total_w - hosted_w, 'bank loads do not hold the unhosted residuals'
D = 73 * W - sum(r * n for r, n in hist.items())
assert D == 2560, D
assert {str(r): n for r, n in sorted(hist.items())} == res['hist'], 'histogram differs'
print(f'AUDIT PASS {name} copies={copies} defer={defer}: {len(owner)} hosted residuals, W={W}, D={D}, histogram matches')
