# v-ex2-certnum: re-derive every degree-0 (and, with arg 1, degree-1 all-variable) certificate claimed in ex2-certnum's
# bulk run, exactly over Q(i), own code.  A record passes if sum_k c_k f_k = 1 is found and re-verified by expansion.
import sys, json, glob, itertools, time
from fractions import Fraction as Fr
from collections import Counter
base = __import__('os').path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), '..', 'out') + '/'
nag = __import__('os').environ.get('A2NAG_DIR', 'a2-nag') + '/'
Z = (Fr(0), Fr(0))
def add(a, b): return (a[0] + b[0], a[1] + b[1])
def sub(a, b): return (a[0] - b[0], a[1] - b[1])
def mul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def inv(a):
    d = a[0] * a[0] + a[1] * a[1]; return (a[0] / d, -a[1] / d)
def load(f):
    d = json.load(open(f)); eqs = []
    for e in d['eqs']:
        p = {}
        for mon, c in e:
            k = tuple(sorted(mon)); v = (Fr(c[0]), Fr(c[1]))
            p[k] = add(p.get(k, Z), v)
        eqs.append({k: v for k, v in p.items() if v != Z})
    return d, eqs
def mons_upto(n, d):
    out = [()]
    for k in range(1, d + 1): out += list(itertools.combinations_with_replacement(range(n), k))
    return out
def cert(eqs, n, d):
    cols = []
    for k, e in enumerate(eqs):
        for mu in mons_upto(n, d):
            col = {}
            for m, v in e.items():
                key = tuple(sorted(m + mu)); col[key] = add(col.get(key, Z), v)
            cols.append({m: c for m, c in col.items() if c != Z})
    rows = {}
    for j, col in enumerate(cols):
        for m, c in col.items(): rows.setdefault(m, {})[j] = c
    if () not in rows: return False
    pivots = []
    for m, r0 in rows.items():
        r = dict(r0); b = (Fr(1 if m == () else 0), Fr(0))
        for pc, pr, pb in pivots:
            if pc in r:
                f = r[pc]
                for jj, v in pr.items():
                    nv = sub(r.get(jj, Z), mul(f, v))
                    if nv != Z: r[jj] = nv
                    else: r.pop(jj, None)
                b = sub(b, mul(f, pb))
        if not r:
            if b != Z: return False
            continue
        pc = min(r, key=lambda jj: (len(cols[jj]), jj)); iv = inv(r[pc])
        r = {jj: mul(v, iv) for jj, v in r.items()}; b = mul(b, iv)
        newp = []
        for qc, qr, qb in pivots:
            if pc in qr:
                f = qr[pc]
                for jj, v in r.items():
                    nv = sub(qr.get(jj, Z), mul(f, v))
                    if nv != Z: qr[jj] = nv
                    else: qr.pop(jj, None)
                qb = sub(qb, mul(f, b))
            newp.append((qc, qr, qb))
        pivots = newp + [(pc, r, b)]
    tot = {}
    for pc, pr, pb in pivots:
        if pb != Z:
            for m, c in cols[pc].items(): tot[m] = add(tot.get(m, Z), mul(pb, c))
    tot = {m: c for m, c in tot.items() if c != Z}
    return tot == {(): (Fr(1), Fr(0))}
recs = {}
for f in sorted(glob.glob(base + 'bulk_*.jsonl')) + sorted(glob.glob(base + 'bulk2_*.jsonl')):
    for l in open(f):
        r = json.loads(l)
        if r['file'] in recs and not r['nulla'] and recs[r['file']]['nulla']: continue
        recs[r['file']] = r
want = int(sys.argv[1])
cnt = Counter(); t0 = time.time(); fails = []
for fn, r in sorted(recs.items()):
    if not r['nulla'] or not r['nulla'].startswith('degree %d' % want): continue
    D, eqs = load(nag + fn[3:])
    ok = cert(eqs, D['nvar'], want); cnt[ok] += 1
    if not ok: fails.append(fn)
print('degree', want, 'claimed records:', sum(cnt.values()), 'certificate re-derived:', cnt[True], 'NOT re-derived:', cnt[False], fails[:10], '(%.0fs)' % (time.time() - t0))
# soundness: the Singular-feasible ones must have no certificate at this degree
fe = Counter()
for fn, r in sorted(recs.items()):
    if r['singular'].startswith('feasible'):
        D, eqs = load(nag + fn[3:]); fe[cert(eqs, D['nvar'], want)] += 1
print('Singular-feasible records: certificate at degree %d found for' % want, fe[True], 'of', sum(fe.values()))
