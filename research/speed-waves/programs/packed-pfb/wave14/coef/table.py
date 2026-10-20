# Collect every quantiser result line (qA_*, qB_*, q2*_*.out) into one markdown table.
import json, glob
rows = []
for fn in sorted(glob.glob('q*.out')):
    proto = 'B' if '_B' in fn or fn.startswith('qB') or fn.startswith('q2B') else 'A'
    grp = '7' if 'g7' in fn else ('4' if 'g4' in fn else '14')
    for l in open(fn):
        try: r = json.loads(l)
        except Exception: continue
        if r['case'] == 'float' and grp != '14': continue
        rows.append((proto, r['case'], grp, r.get('quantiser', 'float (no quantisation)'), r.get('peak', '-'), r['ripple'],
                     r['stop1'], r['stop2'], r['med2'], r.get('frac80', float('nan')), r['errfloor']))
seen = set()
print("| proto | case | taps per packed group | quantiser | peak coef | ripple +-dB | worst dB [0.685,1.5) ch | worst dB >= 1.5 ch | median dB >= 1.5 ch | share of stop-band <= -80 dB | coef error floor dB |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    k = (r[0], r[1], r[2], r[3])
    if k in seen: continue
    seen.add(k)
    print(f"| {r[0]} | {r[1]} | {r[2] if r[1].startswith('packed') else '-'} | {r[3]} | {r[4]} | {r[5]:.3f} | {r[6]:.2f} | {r[7]:.2f} | {r[8]:.1f} | {100*r[9]:.1f} % | {r[10]:.1f} |")
