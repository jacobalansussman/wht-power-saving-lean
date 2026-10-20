# markdown table of the representation runs (NOTE.md section 3)
import json, sys
rows = []
for fn in sys.argv[1:]:
    for l in open(fn):
        o = json.loads(l)
        if 'error' in o or 'r' not in o: continue
        rows.append(o)
print('| Run | P | Eeff | naive P − Eeff | r | codim Eeff − r | dominant | fibre P − r | orbit o | o ≤ r | dim R / max |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for o in rows:
    print(f"| {o['label']} | {o['P']} | {o['Eeff']} | {o['naive_eff']} | {o['r']} | {o['Eeff'] - o['r']} | {'yes' if o['dominant'] else 'no'} | {o['fibre']} | {o['orbit']} | {'yes' if o['orbit_ok'] else 'NO'} | {o['dimR']} / {o['Eeff'] // o['W']} |")
