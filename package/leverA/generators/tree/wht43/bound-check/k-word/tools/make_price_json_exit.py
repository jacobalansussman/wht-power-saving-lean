"""Write a price JSON (+ .schedule.json.gz beside it) in the format generate_global_lean.py (Shea's, via
d-lean-route) reads, from a pack_ceil.py table.  Usage: make_price_json_exit.py TABLE.schedule.json.gz PROFILE_JSON OUT.json [--copies N]  (k-word: N source copies, default 112)"""
import sys, os, json, gzip, shutil
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_moment import certify
from route_price import IDLE, price
from pack_exact import validate
tab, prof, out = sys.argv[1], {int(k): n for k, n in json.loads(sys.argv[2]).items()}, sys.argv[3]
s = json.loads(gzip.decompress(open(tab, 'rb').read()))
validate(s['rows'], s['pads'], s['widths'])
local = Counter(prof)
ncop = int(sys.argv[sys.argv.index('--copies') + 1]) if '--copies' in sys.argv else 112   # k-word: 0 / 5 for the stacked variants
if '--histogram-includes-copies' not in sys.argv: local[1] += ncop; local[8] += ncop
hist = Counter({r: 72 * n for r, n in local.items()})
for r, n in IDLE.items(): hist[r] += n * 112
hist.update(p for p in s['pads'] if p)
W = 64 * 112 + len(s['pads']); assert 73 * W - sum(r * n for r, n in hist.items()) == 2560
rate = certify(dict(hist), 73, W)
pr = price(dict(hist=dict(hist), W=W, D=2560))
json.dump(dict(status='EXACT_ARITHMETIC_CANDIDATE_TIER2_EXIT_FRAMES_PAD_MINIMAL_SCHEDULE', h=9, m=73, v=112,
               local_profile={str(k): v for k, v in sorted(local.items())}, local_residual_count=len(s['widths']),
               residual_rule=('r_q = dim(NATURAL exit frame: frame after the last gate, tail and cleanup dropped) - dim(entrance frame)' if '--histogram-includes-copies' in sys.argv else 'r_q = dim(exit frame U_q) - dim(entrance frame E_q); copies r = 9; source copies kept: %d' % ncop),
               banks=len(s['pads']), W=W, deficit=2560, profile={str(k): v for k, v in sorted(hist.items())},
               pad_hist={str(k): v for k, v in sorted(Counter(s['pads']).items())}, exact_rate=rate,
               kernel_A_s40=pr['kernel_A'], dft_delta_1e10=pr['dft']), open(out, 'w'), indent=2)
shutil.copy(tab, out[:-5] + '.schedule.json.gz')
print('W', W, 'banks', len(s['pads']), 'root', rate['root_lower_decimal'][:16], 'kernel', pr['kernel_A'], 'dft', pr['dft'])
