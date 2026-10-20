import re, json, sys
s = open('out_ctbl.txt').read().replace('\\\n', '')
s = re.sub(r'\n\s+', ' ', s)          # GAP wraps lists over lines with leading spaces
L = s.split('\n')
def grab(prefix):
    for l in L:
        if l.startswith(prefix): return l[len(prefix):].strip()
def lst(x): return json.loads(x.replace('"', '"'))
D = {}
D['classes'] = json.loads(grab('class names'))
D['degrees'] = json.loads(grab('degrees'))
p0 = grab('P0 order'); D['P0_degrees'] = json.loads(p0.split('irr degrees')[1])
D['P0_classes'] = json.loads(grab('P0 class names'))
D['P0_orders'] = json.loads(grab('P0 class orders'))
D['fusion'] = json.loads(grab('fusion P0->G'))
rows = []
i = [k for k, l in enumerate(L) if l.startswith('RESTRICTION')][0]
for l in L[i+1:i+31]:
    a, b = l.split(':', 1); v = json.loads(b.split('sumsq')[0]); rows.append(v)
D['restriction'] = rows
i = [k for k, l in enumerate(L) if l.startswith('Irr values table')][0]
D['irr'] = [json.loads(l) for l in L[i+1:i+31]]
D['centralizers'] = json.loads(grab('centralizer orders'))
D['orders'] = json.loads(grab('element orders'))
i = [k for k, l in enumerate(L) if l.startswith('P0 Irr values')][0]
D['P0_irr'] = [[x.strip() for x in l.strip()[1:-1].split(',')] for l in L[i+1:i+25]]
D['P0_centralizers'] = json.loads(grab('P0 centralizers'))
dd = grab('2.S6(2) degrees'); D['cover_degrees'] = json.loads(dd) if dd else None
json.dump(D, open('ctbl.json', 'w'))
print(len(D['irr']), len(D['P0_irr']), D['cover_degrees'])
print('classes', D['classes'])
print('2a..2d values per irrep:')
ci = [D['classes'].index(c) for c in ['1a','2a','2b','2c','2d']]
for chi, r in zip(D['irr'], D['restriction']):
    mf = max(r) <= 1
    print(chi[0], [chi[c] for c in ci], 'P0-mult-free' if mf else 'REPEATED', 'pieces', sum(r), 'dims', [D['P0_degrees'][j] for j, m in enumerate(r) for _ in range(m)])
