#!/usr/bin/env python3
"""g-heads: price from OUR lowering's G0 output (tools/g0_staged.py JSON): (1) b of OUR global profile by the hub's
moment engine (m-anatomy's exact re-implementation, 1e-24 grid) and by the mixed-target Decimal root, compared with the
hub's stated root digit for digit; (2) the lever (half cut, target 49) with OUR literal head count, through
lever-reprice's variants() (same formula as rl-plan's rescue.py; its T74 controls reproduce rl-plan).
usage: price_ours.py G0.json STATED_ROOT_LOWER [STATED_ROOT_UPPER]"""
import sys, json
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pricer import price_b, dec
from mixed import root, D
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from lever_price import variants          # (its module body prints the controls; silenced here, run separately)
g = json.load(open(sys.argv[1]))
lo = sys.argv[2]; hi = sys.argv[3] if len(sys.argv) > 3 else None
prof = {int(k): v for k, v in g['call_blocks']['global_profile'].items()}
T = g['fill']['replicas']; W = 896 * T
b_eng = dec(price_b(prof, 73, W), 25)
b_mix = root(prof, {73: W})
heads = g['G0']['heads_rank48_unhosted']
non = prof[48] - heads
out = dict(label=g['label'], W=W, profile=prof, profile_equals_hub=g['call_blocks']['global_profile_equals_hub'],
           b_hub_engine=str(b_eng), b_mixed_root='%.24e' % b_mix, stated_root=[lo, hi],
           b_reproduces_stated_lower_digit_for_digit=D(str(b_eng)) == D(lo),
           heads=heads, other_rank48=non)
with contextlib.redirect_stdout(io.StringIO()):
    v = variants(g['label'], prof, W, heads, non, quiet=True)
hub = v['uniform']['b']
for k in ('uniform', 'half_cut', 'full_cut'):
    out[k] = dict(b='%.10e' % v[k]['b'], over_uniform_pct='%+.3f' % (float((v[k]['b'] - hub) / hub) * 100), separated=v[k]['separated'])
out['breakeven_full_cut_heads'] = v['breakeven_full']
out['half_cut_breakeven_heads'] = 1 if heads >= 1 else None
print(json.dumps(out, indent=1))
