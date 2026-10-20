"""v4gb.py (w8-b5sparse, resumed 15 Oct): EXACT V4 for layer-4 systems with three- and four-term rows (Singular), and
the d = 3 ball test on it, with no free-parameter assumption.

Layered system (endsolve format, from dumpsys.py pickles): rows r, columns c, entries (r, c, cf), value of an entry
sum_k v_k a_k (a_0 = 1, a_k = alpha_k).  Unknowns: alpha (shared by both outputs), one w per column and output.
  d = 4 equations: every row of both outputs' complete d = 4 closures; seed normalisation (the seed's K-value = 1,
  exactly as modsolve.build(norm='seed'), K-row weights fixed small integers).
  all alpha_k != 0: one extra variable t with t * prod(alpha_k) - 1 (Rabinowitsch).
  V4 = projection to alpha of that variety (elimination); printed with its minimal associated primes.
  d = 3 ball (mode 'ball'): rows within structural BFS distance k of the seed in the TRUNCATED d = 3 closure of each
  output (complete rows only, every column touching them: a valid necessary test at any cap), seed normalisation,
  fresh column unknowns, the SAME alpha.  1 in the ideal  <=>  no alpha with all alpha_k != 0 satisfies the d = 4 system
  and both balls: the design is excluded at layer 3 (over the algebraic closure of the coefficient field).
usage: v4gb.py v4   sys_d4.pkl out.sing [char]
       v4gb.py ball sys_d4.pkl sys_d3.pkl k out.sing [char] [both|1|2]
       v4gb.py ball4 sys_d4.pkl k out.sing [char]   (positive control: the d = 4 system's own ball, must be consistent)
char: 0 (rationals, default) or a prime < 2^31."""
import sys, pickle
from fractions import Fraction


def rat(v):
    return Fraction(v).limit_denominator(10 ** 6)


def pstr(cf):
    out = []
    for k, v in sorted(cf.items()):
        f = rat(v)
        if f == 0:
            continue
        c = ('%d/%d' % (f.numerator, f.denominator)) if f.denominator != 1 else str(f.numerator)
        out.append('(' + c + ')' + ('' if k == 0 else '*a(%d)' % k))
    return '+'.join(out) if out else '0'


def equations(cl, seedcols, var, rowsel=None):
    """polys (as strings) of one output's system; var(c) -> variable name; rowsel: set of row indices (None = all)"""
    rows, cols, ent, kl = cl
    klass, kterms = kl
    byrow = {}
    for (r, c, cf) in ent:
        if rowsel is not None and r not in rowsel:
            continue
        byrow.setdefault(r, []).append('(' + pstr(cf) + ')*' + var(c))
    eqs = ['+'.join(v) for r, v in sorted(byrow.items())]
    seedrows = sorted(set(rk for c in seedcols for (rk, cf) in kterms.get(c, [])))
    lam = {rk: i + 1 for i, rk in enumerate(seedrows)}
    nt = []
    used = set()
    for c in klass:
        for (rk, cf) in kterms.get(c, []):
            if rk in lam:
                nt.append('(%d)*(%s)*%s' % (lam[rk], pstr(cf), var(c))); used.add(c)
    eqs.append('+'.join(nt) + '-1')
    cset = set(c for (r, c, cf) in ent if rowsel is None or r in rowsel) | used
    return eqs, cset, len(seedrows)


def alphas_in(cl, rowsel=None):
    s = set()
    for (r, c, cf) in cl[2]:
        if rowsel is None or r in rowsel:
            s |= set(k for k in cf if k)
    for c, lst in cl[3][1].items():
        for (rk, cf) in lst:
            s |= set(k for k in cf if k)
    return s


def ball(cl, seedcols, k):
    rows, cols, ent, kl = cl
    rowsof, colsof = {}, {}
    for (r, c, cf) in ent:
        rowsof.setdefault(c, set()).add(r); colsof.setdefault(r, set()).add(c)
    R = set(); C = set(seedcols); fr = set(seedcols)
    for _ in range(k):
        nr = set()
        for c in fr:
            nr |= rowsof.get(c, set())
        nr -= R; R |= nr
        nc = set()
        for r in nr:
            nc |= colsof.get(r, set())
        nc -= C; C |= nc; fr = nc
    return R


def leafpeel(cl, seedcols, rowsel):
    """EXACT reduction: a column that occurs in exactly one remaining row, not in the normalisation, with an entry that
    is a single nonzero term (a constant or one alpha_k, nonzero since every alpha_k != 0), can always absorb that
    row: drop the row (and the column).  Repeat.  The reduced system is solvable iff the original is (for that alpha)."""
    rows, cols, ent, kl = cl
    klass, kterms = kl
    seedrows = set(rk for c in seedcols for (rk, cf) in kterms.get(c, []))
    prot = set(c for c in klass for (rk, cf) in kterms.get(c, []) if rk in seedrows)
    E = {}
    for (r, c, cf) in ent:
        if rowsel is not None and r not in rowsel:
            continue
        d = E.setdefault((r, c), {})
        for k, v in cf.items():
            d[k] = d.get(k, 0) + v
    rowsof, colsof = {}, {}
    for (r, c), d in E.items():
        if any(abs(v) > 1e-12 for v in d.values()):
            rowsof.setdefault(c, set()).add(r); colsof.setdefault(r, set()).add(c)
    alive = set(colsof)
    queue = [c for c in rowsof if len(rowsof[c]) == 1 and c not in prot]
    dropped = 0
    while queue:
        c = queue.pop()
        if c in prot or len(rowsof.get(c, ())) != 1:
            continue
        r = next(iter(rowsof[c]))
        d = {k: v for k, v in E[(r, c)].items() if abs(v) > 1e-12}
        if len(d) != 1:
            continue
        alive.discard(r); dropped += 1
        for c2 in colsof[r]:
            rowsof[c2].discard(r)
            if len(rowsof[c2]) == 1 and c2 not in prot:
                queue.append(c2)
    return alive, dropped


def header(char, nvars):
    return 'option(redSB);\nLIB "elim.lib";\nLIB "primdec.lib";\n'


if __name__ == "__main__":
    mode = sys.argv[1]
    D4 = pickle.load(open(sys.argv[2], 'rb')); nal = D4['nal']
    sysd4 = D4['sys']
    blocks = []                                   # (prefix, cl, seedcols, rowsel)
    import os as _o0
    for l in (0, 1):
        cl, sc, cap = sysd4[l]
        if cap:
            print('NOTE: d=4 closure of output %d truncated at the cap: its complete rows only (still a valid test)' % (l + 1))
        if not int(_o0.environ.get('NO4', '0')):      # NO4=1: leave the d=4 systems out (ball alone)
            blocks.append(('x%d' % (l + 1), cl, sc, None))
    if mode == 'v4':
        out, char = sys.argv[3], (sys.argv[4] if len(sys.argv) > 4 else '0')
    elif mode == 'ball':
        D3 = pickle.load(open(sys.argv[3], 'rb')); k = int(sys.argv[4]); out = sys.argv[5]
        char = sys.argv[6] if len(sys.argv) > 6 else '0'
        which = sys.argv[7] if len(sys.argv) > 7 else 'both'
        for l in (0, 1):
            if which != 'both' and int(which) != l + 1:
                continue
            cl, sc, cap = D3['sys'][l]
            blocks.append(('z%d' % (l + 1), cl, sc, ball(cl, sc, k)))
    elif mode == 'ball4':
        k = int(sys.argv[3]); out = sys.argv[4]; char = sys.argv[5] if len(sys.argv) > 5 else '0'
        for l in (0, 1):
            cl, sc, cap = sysd4[l]
            blocks.append(('z%d' % (l + 1), cl, sc, ball(cl, sc, k)))
    eqs = []; varnames = []; info = []
    usedA = set()
    import os as _o
    LEAF = int(_o.environ.get('LEAF', '1'))
    for (pre, cl, sc, rs) in blocks:
        cols = cl[1]
        if LEAF:
            nb = len(rs) if rs is not None else len(set(r for (r, c, cf) in cl[2]))
            rs, ndr = leafpeel(cl, sc, rs)
            info.append('%s leafpeel: %d -> %d rows' % (pre, nb, len(rs)))
        cmap = {}                                  # compact numbering of the columns that occur (ring size)
        e, cset, nsr = equations(cl, sc, lambda c, pre=pre, cmap=cmap: '%s(%d)' % (pre, cmap.setdefault(c, len(cmap) + 1)), rs)
        eqs += e
        varnames.append('%s(1..%d)' % (pre, len(cmap)))
        usedA |= alphas_in(cl, rs)
        info.append('%s: %d rows, %d columns used, %d seed K-rows' % (pre, len(e) - 1, len(cset), nsr))
    A = sorted(usedA)
    wvars = ','.join(varnames)
    s = header(char, 0)
    s += '// %s\n' % ' ; '.join(info)
    s += '// alphas used: %s\n' % A
    import os
    INV = int(os.environ.get('INV', '0'))           # 1: alpha_k != 0 via a(k)*b(k) = 1 (one inverse per alpha)
    if INV:
        s += 'ring R = %s, (t,%s,b(1..%d),a(1..%d)), dp;\n' % (char, wvars, nal, nal)
        s += 'ideal I =\n' + ',\n'.join(eqs) + ',\n' + ',\n'.join('a(%d)*b(%d)-1' % (k, k) for k in A) + ',\nt;\n'
    else:
        s += 'ring R = %s, (t,%s,a(1..%d)), dp;\n' % (char, wvars, nal)
        s += 'ideal I =\n' + ',\n'.join(eqs) + ',\nt*' + '*'.join('a(%d)' % k for k in A) + '-1;\n'
    s += 'int t0 = timer;\n'
    if mode == 'v4':
        wprod = 't*' + ('*'.join('b(%d)' % k for k in A) + '*' if INV else '') + '*'.join('%s' % v.split('(')[0] + '(%d)' % (i + 1) for v in varnames
                                for i in range(int(v.split('..')[1][:-1])))
        s += 'ideal J = eliminate(I, %s);\n' % wprod
        s += 'J = simplify(J, 2+8);\n'
        s += '"V4 ideal (in alpha): generators", size(J), "dim", dim(std(J));\n'
        s += 'print(J);\n'
        s += 'list L = minAssGTZ(J);\n"minimal associated primes:", size(L);\n'
        s += 'int i; for (i = 1; i <= size(L); i++) { "component", i, "dim", dim(std(L[i])); print(L[i]); }\n'
    else:
        s += 'ideal G = std(I);\n'
        s += 'if (size(G) == 1 && deg(G[1]) == 0) { "VERDICT: 1 in ideal -> EXCLUDED (no alpha, all alpha_k != 0)"; }\n'
        s += 'else { "VERDICT: consistent (1 not in ideal); dim", dim(G); }\n'
    s += '"time (s)", timer - t0;\nquit;\n'
    open(out, 'w').write(s)
    print('wrote %s: %d equations; %s; alphas used %d of %d' % (out, len(eqs) + 1, ' ; '.join(info), len(A), nal))
