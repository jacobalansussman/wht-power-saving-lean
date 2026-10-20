#!/usr/bin/env python3
"""d-bit-exits / windows9.py: the nine-stage idle-hosting windows of PR #510's chain (PR #463 template, PR #442 hosting),
re-derived from global_k.py's recipe and idle schedule (K = 4, m = 64) and hosting_core.windows (both re-implemented
here verbatim in logic; read from the hub's vendored files). Prints each idle child's rank and the stages whose helper
residuals it may host, then checks the frozen #589 selection's remainder against fastc.py's REM.
usage: python3 -I windows9.py SELECTION9.json"""
import json, sys
from collections import Counter
K = 4
def recipe(K):
    out = [dict(kind='local', stage=0, pair=0, reverse=False)]; stage = 1
    for j in range(K):
        out.append(dict(kind='local', stage=stage, pair=j, reverse=True)); stage += 1
        if j < K - 1: out.append(dict(kind='bridge', dst=j, src=j + 1, coefficient=1))
        out.append(dict(kind='local', stage=stage, pair=j, reverse=False)); stage += 1
        if j < K - 1: out.append(dict(kind='bridge', dst=j, src=j + 1, coefficient=-1))
    for j in reversed(range(K - 1)): out.append(dict(kind='bridge', dst=j + 1, src=j, coefficient=1))
    out.append(dict(kind='exchange'))
    return out
levels = [-1] * K; sched = []
for e in recipe(K):
    if e['kind'] == 'local':
        s, p = e['stage'], e['pair']; levels[p] = s
        sched.append(dict(kind='local', stage=s, pair=p, reverse=e['reverse']))
        if s % 2 == 1 and s <= 2 * K - 3:
            j = (s + 1) // 2; levels[j] = s; sched.append(dict(kind='idle', pair=j, rank=7 * (s + 1)))
        elif s % 2 == 0 and 2 <= s <= 2 * K - 2:
            j = s // 2; levels[j] = s; sched.append(dict(kind='idle', pair=j, rank=7))
        elif s == 2 * K:
            for j in range(K):
                old = levels[j]; sched.append(dict(kind='idle', pair=j, rank=7 * (2 * K - old))); levels[j] = 2 * K
    elif e['kind'] == 'bridge':
        sched.append(dict(kind='bridge', dst=e['dst'], src=e['src']))
idx = 0
for e in sched:
    if e['kind'] == 'idle': e['idx'] = idx; idx += 1
def windows(schedule):
    last = {}; out = {}; locals_ = []
    for pos, e in enumerate(schedule):
        k = e['kind']
        if k == 'local':
            locals_.append((pos, e['stage'], e['pair'])); last[e['pair']] = pos
        elif k == 'bridge':
            last[e['dst']] = pos; last[e['src']] = pos
        elif k == 'idle':
            p = e['pair']; lo = last.get(p, -1)
            st = frozenset(s for q, s, pp in locals_ if lo < q < pos and pp != p)
            out[e['idx']] = (p, e['rank'], st)
    return out
win = {k: x for k, x in windows(sched).items() if x[1] > 0 and x[2]}
allidle = [(e['idx'], e['rank']) for e in sched if e['kind'] == 'idle']
print('idle children (idx, rank):', allidle)
for i, (p, r, st) in sorted(win.items()): print('window', i, 'pair', p, 'rank', r, 'stages', sorted(st))
if len(sys.argv) > 1:
    sel = json.load(open(sys.argv[1]))['bit']
    rem = Counter()
    for i, pat, n in sel['patterns']:
        r = win[i][1] - sum(pat)
        if r: rem[r] += 32 * n
    for i, r in allidle:
        if r and i not in win: rem[r] += 32 * 222
    print('remainder from the frozen selection:', dict(sorted(rem.items())))
    assert dict(rem) == {2: 13920, 4: 14208, 6: 14208, 7: 21312, 42: 288}, rem
    print('equals fastc.py REM {2:13920, 4:14208, 6:14208, 7:21312, 42:288}')
json.dump({str(i): [p, r, sorted(st)] for i, (p, r, st) in win.items()}, open('windows9.json', 'w'), indent=1)
