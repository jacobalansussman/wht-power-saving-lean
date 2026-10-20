"""Write the price JSON that Shea's generate_global_lean.py reads (and copy the schedule beside it as
<stem>.schedule.json.gz).  Usage: mk_price_json.py <schedule.json.gz> <out/stem.json>"""
import sys, os, json, gzip, shutil
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from price_sched import price, LOCAL
sched, out = sys.argv[1], sys.argv[2]
r = price(sched)
dst = out[:-len('.json')] + '.schedule.json.gz'
if os.path.abspath(dst) != os.path.abspath(sched):
    shutil.copyfile(sched, dst)
json.dump(dict(status='EXACT_ARITHMETIC_CANDIDATE_EXPLICIT_SCHEDULE', h=9, m=73, v=112,
               local_profile={str(k): v for k, v in LOCAL.items()}, local_residual_count=616, banks=r['banks'], W=r['W'],
               deficit=2560, profile=r['profile'], pad_hist=r['pad_hist'], exact_rate=r['rate'],
               kernel_numerator=r['kernel_A'], dft_numerator=r['dft']), open(out, 'w'), indent=2)
print(out, 'banks', r['banks'], 'W', r['W'], 'kernel', r['kernel_A'], 'dft', r['dft'])
