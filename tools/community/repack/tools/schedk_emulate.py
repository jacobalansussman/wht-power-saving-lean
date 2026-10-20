"""Python emulation of the Lean checker GX.schedK (Work/Bridge401/SchedK.lean, verbatim semantics: 72 rows of
length R, entries < nb and distinct per row via a Nat bitset, nb pads < 73, R residuals with 72*sum < 2^20, and the
packed identity  sum_rows sum_q res_q * 2^(20*bank) = sum_b (73 - pad_b) * 2^(20*b)), plus the rank counts
(wsumK wGt top = 0, wEq r counts) that C588Sched.lean proves by `decide +kernel`.
Usage: schedk_emulate.py <schedule.json.gz> <nb expected>"""
import sys, json, gzip
from collections import Counter
s = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
rows, pad, res = s['rows'], s['pads'], s['widths']; nb = int(sys.argv[2]); R = 616
def skInj(row):
    seen = 0
    for x in row:
        if not (x < nb and not (seen >> x) & 1): return False
        seen |= 1 << x
    return True
def skRow(row): return sum(r << (20 * x) for x, r in zip(row, res))
def skTgt(): return sum((73 - p) << (20 * i) for i, p in enumerate(pad))
ok = (len(rows) == 72 and all(len(r) == R and skInj(r) for r in rows) and len(pad) == nb and all(p < 73 for p in pad)
      and len(res) == R and 72 * sum(res) < 2 ** 20 and sum(skRow(r) for r in rows) == skTgt())
bad = [list(r) for r in rows]; bad[0][1] = bad[0][0]
okbad = (all(skInj(r) for r in bad) and sum(skRow(r) for r in bad) == skTgt())
top = max(pad)
print('schedK', ok, '| duplicate-owner control accepted?', okbad, '| top pad', top, '| rank counts', sorted(Counter(pad).items()))
assert ok and not okbad
