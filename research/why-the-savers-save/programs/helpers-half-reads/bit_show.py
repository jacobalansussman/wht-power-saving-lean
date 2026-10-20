import sys, json, gzip
from collections import defaultdict
w = json.load(gzip.open(sys.argv[1], 'rt')); V = w['v']; R = w['R']; N = 2*V+R
fr = w['frames']
T = 0; rows = []
for sec in ('D', 'A', 'B', 'C'):
    for g in w[sec]:
        rows.append((T, sec, g)); T += 1
    if sec == 'A': rows.append((T, 'scatter', None)); T += 1
for q in map(int, sys.argv[2:]):
    print('REG', q, 'start frame', w['start'][q], fr[w['start'][q]], 'final', w['final'][q], len(fr[w['final'][q]]))
    for T, sec, g in rows:
        if g and (g[1] == q or q in g[2]):
            print('   ', T, sec, 'frame', g[0], fr[g[0]][:12], 'target', g[1], 'sources', g[2])
