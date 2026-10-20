"""v-mcopy-code-rules probe 2: the cut channel. Which registers do cut gates touch; do the group's needed (packet)
sources reach the group targets through the cut; what frame data do cut gates carry; centre chains."""
import sys, json
sys.path.insert(0, '../helpers-m-copies')
from wordlib import *
w = load(sys.argv[1]); v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2*v + R
cls = lambda q: 'x' if q < v else 'y' if q < 2*v else 'h'
print('cut gates', len(w['cut']), 'example', w['cut'][:3])
print('cut (dst,src) classes', Counter((cls(e[1]), cls(e[2])) for e in w['cut'] if e[0]=='add'))
print('centers', w['centers'], [cls(q) for q in w['centers']], 'center frame dims', [len(U) for U in w['center_frames']])
print('target_levels count', len(w.get('target_levels', {})), 'dims', Counter(len(U) for U in w.get('target_levels', {}).values()))
G = json.load(open('../helpers-m-copies/mc_groups.json'))
start, end, chain, paid, runs = frames(w)
prog = program(w)
for gi in range(len(G)):
    g = G[gi]; mem = set(g['members']); T = set(g['targets'])
    allM = set(q for x in G for q in x['members'])
    rows = [{q: F(1)} for q in range(n)]; need = set()
    cutd = defaultdict(lambda: defaultdict(F))
    for ph, i, e in prog:
        if e[0] == 'add' and ph == 'B' and e[1] in mem and e[2] not in allM:
            xs = [s for s in rows[e[2]] if s < v]
            if len(xs) == 8: need |= set(xs)
        if e[0] == 'add' and ph == 'cut' and e[1] in T:
            for s, c in rows[e[2]].items():
                if s < v: cutd[e[1]][s] += F(e[3], e[4]) * c
        apply(e, rows)
    hit = {t: sorted(s for s, c in cutd[t].items() if c and s in need) for t in sorted(T)}
    cent = {t: sorted({e[2] for e in w['cut'] if e[1] == t}) for t in sorted(T)}
    if gi < 3 or True:
        print('group', gi, 'needed', len(need), 'needed sources delivered to each group target through the CUT:', {t: len(x) for t, x in hit.items()}, 'centres read', cent)
# are centre chain frames / cut frames present in d-mcopy-code reach chains?
c0 = w['centers'][0]
print('centre', c0, 'chain', [(t, len(U)) for t, U in chain[c0]])
