"""vball.py (skeptic v-w8-b5sparse, own code): from the agent's pickled layered systems, write the one-shot exact
ball test myself (no leaf peeling, no code of the agent), for Singular and msolve.

System: for each output l in OUTS: rows of the TRUNCATED d=3 pickle within BFS radius K of the seed column(s)
(rows containing a seed column = radius 1; rows containing any column of a radius-1 row = radius 2; ...),
all their entries; optional (D4=1) every row of the d=4 pickle with separate unknowns.
Normalisation: sum over columns of their K-term coefficient on the seed's K-row(s), weighted 1 (assert exactly one
seed K-row), = 1.  Nonzero alphas: a_k*b_k - 1 for every alpha that occurs.  Entries must be exact small integers.
usage: vball.py d4.pkl d3.pkl K char out_prefix [D4=0|1] [OUTS=12]
Writes out_prefix.sing (Singular std) and out_prefix.ms (msolve; char must be a prime, products of distinct vars only)."""
import sys, pickle, os


def intval(v):
    iv = int(round(float(v.real if isinstance(v, complex) else v)))
    assert abs(complex(v) - iv) < 1e-12, ('non-integer entry', v)
    return iv


def ballrows(ent, seeds, K):
    rows_of_col, cols_of_row = {}, {}
    for (r, c, cf) in ent:
        rows_of_col.setdefault(c, set()).add(r)
        cols_of_row.setdefault(r, set()).add(c)
    cols = set(seeds); rows = set()
    for _ in range(K):
        newr = set()
        for c in cols:
            newr |= rows_of_col.get(c, set())
        newr -= rows
        rows |= newr
        for r in newr:
            cols |= cols_of_row[r]
    return rows


def system(cl, seeds, rowsel, pre):
    rows, cols, ent, (klass, kterms) = cl
    polys = {}
    used_a = set(); used_x = set()
    for (r, c, cf) in ent:
        if rowsel is not None and r not in rowsel:
            continue
        for k, v in cf.items():
            iv = intval(v)
            if iv == 0:
                continue
            x = '%sc%d' % (pre, c); used_x.add(x)
            mon = (('a%d' % k,) if k else ()) + (x,)
            if k:
                used_a.add(k)
            d = polys.setdefault(r, {})
            d[mon] = d.get(mon, 0) + iv
    srows = set()
    for c in seeds:
        for (rk, cf) in kterms.get(c, []):
            srows.add(rk)
    assert len(srows) == 1, ('seed K-rows', srows)
    norm = {}
    for c, lst in kterms.items():
        for (rk, cf) in lst:
            if rk in srows:
                for k, v in cf.items():
                    iv = intval(v)
                    if iv:
                        x = '%sc%d' % (pre, c); used_x.add(x)
                        mon = (('a%d' % k,) if k else ()) + (x,)
                        if k:
                            used_a.add(k)
                        norm[mon] = norm.get(mon, 0) + iv
    norm[()] = -1
    out = [d for d in polys.values() if any(v for v in d.values())] + [norm]
    return out, used_a, used_x


def fmt(d, p):
    terms = []
    for mon, v in sorted(d.items()):
        if v == 0:
            continue
        vv = v % p if p else v
        if vv == 0:
            continue
        terms.append(('%d' % vv) + ''.join('*' + m for m in mon))
    return ('+'.join(terms)).replace('+-', '-') if terms else '0'


if __name__ == '__main__':
    d4f, d3f, K, ch, outp = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    D4 = int(os.environ.get('D4', '0')); OUTS = os.environ.get('OUTS', '12')
    P4 = pickle.load(open(d4f, 'rb')); P3 = pickle.load(open(d3f, 'rb'))
    allp = []; A = set(); X = []
    info = []
    for l in (0, 1):
        if str(l + 1) not in OUTS:
            continue
        if D4:
            cl, sc, cap = P4['sys'][l]
            assert not cap
            ps, a, x = system(cl, sc, None, 'w%d' % (l + 1)); allp += ps; A |= a; X += sorted(x)
            info.append('d4 out%d: %d polys, %d cols' % (l + 1, len(ps), len(x)))
        cl, sc, cap = P3['sys'][l]
        rs = ballrows(cl[2], sc, K)
        ps, a, x = system(cl, sc, rs, 'v%d' % (l + 1)); allp += ps; A |= a; X += sorted(x)
        info.append('d3 out%d ball K=%d: %d rows (%d polys incl. norm), %d cols' % (l + 1, K, len(rs), len(ps), len(x)))
    for k in sorted(A):
        allp.append({('a%d' % k, 'b%d' % k): 1, (): -1})
    V = X + ['b%d' % k for k in sorted(A)] + ['a%d' % k for k in sorted(A)]
    p = ch
    # Singular
    s = '// vball.py ' + ' ; '.join(info) + '\nring R = %d, (%s), dp;\nideal I =\n' % (ch, ','.join(V))
    s += ',\n'.join(fmt(d, 0 if ch == 0 else p) for d in allp) + ';\n'
    s += 'int t0 = timer; ideal G = std(I);\n'
    s += 'if (size(G) == 1 && deg(G[1]) == 0) { "VERDICT: 1 in ideal (inconsistent)"; } else { "VERDICT: 1 NOT in ideal; dim", dim(G); }\n"time", timer - t0;\nquit;\n'
    open(outp + '.sing', 'w').write(s)
    if ch:
        m = ','.join(V) + '\n%d\n' % ch + ',\n'.join(fmt(d, p) for d in allp) + '\n'
        assert '^' not in m
        open(outp + '.ms', 'w').write(m)
    print(outp, ' ; '.join(info), '| vars', len(V), 'polys', len(allp), 'alphas', len(A))
