"""k-lower: the y-free twin of a gcert/1 certificate (mirror of Lean's GXD.nyRaw, Work/KLow/NoY.lean): every add whose
TARGET is a y register is dropped from phase B (an `out` gate keeps its other targets, an `in` gate into a y register
keeps its frame and loses its adds) and the scatter table has v empty rows.  The non-y rows of the total matrix are
unchanged (Lean: GXD.nyRaw_Mt); the zero checks run on the twin (y registers stay empty, no tag overflow).
usage: nyraw.py <cert.gcert.json.gz> <out.gcert.json.gz>"""
import sys, json, gzip
c = json.loads(gzip.open(sys.argv[1], 'rt').read())
v = c['v']
isy = lambda r: v <= r < 2 * v
B = []
for g in c['B']:
    if g[0] == 'out':
        B.append(['out', g[1], g[2], [x for x in g[3] if not isy(x[0])], g[4]])
    else:
        B.append(['in', g[1], g[2], [] if isy(g[2]) else g[3]])
c['B'] = B
c['scat'] = {'table': [[] for _ in range(v)]}
open(sys.argv[2], 'wb').write(gzip.compress(json.dumps(c, sort_keys=True).encode(), mtime=0))
