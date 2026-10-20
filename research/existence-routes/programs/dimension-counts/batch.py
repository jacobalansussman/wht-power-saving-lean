# batch driver: python3 batch.py OUT.jsonl SPECFILE.json  (a list of run() specs); appends one JSON line per spec
import sys, json, traceback, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dimcount import run
out = sys.argv[1]; specs = json.load(open(sys.argv[2]))
for sp in specs:
    try: o = run(**sp)
    except Exception as e: o = {'label': sp.get('label'), 'error': repr(e), 'tb': traceback.format_exc()[-400:]}
    with open(out, 'a') as f: f.write(json.dumps(o) + '\n')
    print(json.dumps(o), flush=True)
