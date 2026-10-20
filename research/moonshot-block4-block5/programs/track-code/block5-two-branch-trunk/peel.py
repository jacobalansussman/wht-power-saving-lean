"""peel.py (b5-twobranch): closure with EXACT forced-zero propagation ("peeling").

Every non-free row of layer >= d must have total coefficient 0.  A row whose COMPLETE producer set (all z/y columns on
t1, t2, t3 whose term lands on it; found by frame matching, as endsolve.closure) has exactly one live column forces that
column to 0; killing it removes it from its other rows, which may force more zeros (cascade).  Only rows of live columns
are expanded.  If every class (seed) column dies, the output has NO solution for this generator set: an exact certificate,
no numerics.  Otherwise the live core goes to the usual LM solve.
peel_closure(P, target_l, dmin, seeds, cap) -> dict(rows, cols, ent, kl, stats) in endsolve.closure's format (live part)."""
import numpy as np
from sp2 import M, Fstar


def peel_closure(P, target_l, dmin, seeds, cap=20000):
    E, Sy = P.E, P.S
    Dinv = np.conj(E.D).T
    idx = {}
    for u, T in P.gens.items():
        for (l, g, c) in T:
            idx.setdefault((l, Sy.frame(g)), []).append((u, g))
    Fs = [Sy.frame(P.s[0]), Sy.frame(P.s[1])]
    dcache = {}

    def dinv(gid):
        v = dcache.get(gid)
        if v is None:
            h = E.intern(E.snap(Dinv @ E.MATS[gid])); v = (h, Sy.frame(h)); dcache[gid] = v
        return v

    def free(l, g):
        if E.layer(g) == 0:
            return True
        F = Sy.frame(g)
        if F == Fs[l] or F == Fstar:
            return True
        h, Fh = dinv(g)
        return E.layer(h) == 0 or Fh == Fs[l]

    def prod_keys(row):
        (l, g) = row
        keys = set()
        for (u, gj) in idx.get((l, Sy.frame(g)), ()):
            bm = E.snap(E.MATS[g] @ np.conj(E.MATS[gj]).T)
            if E.is_mono(bm):
                keys.add(('z', u, E.intern(bm)))
        h, Fh = dinv(g)
        for (u, gj) in idx.get((l, Fh), ()):
            bm = E.snap(E.MATS[h] @ np.conj(E.MATS[gj]).T)
            if E.is_mono(bm):
                keys.add(('y', u, E.intern(bm)))
        return keys

    colrows = {}        # created column key -> list of (row, coefdict)
    kterm = {}          # class columns: key -> list of ((l, gid), coefdict)
    rowprod = {}        # row -> set of live producer keys (complete)
    dead = set()
    queue = []          # rows to expand
    stats = dict(created=0, killed=0, triggers=[])

    def kill(key, why=None):
        stack = [(key, why)]
        while stack:
            k, w = stack.pop()
            if k in dead:
                continue
            dead.add(k); stats['killed'] += 1
            if w is not None:
                stats['triggers'].append((w, k))
            for (r, cf) in colrows.get(k, ()):
                s = rowprod.get(r)
                if s is None:
                    continue
                s.discard(k)
                if len(s) == 1:
                    stack.append((next(iter(s)), r))

    def create(key):
        if key in colrows or key in dead:
            return
        typ, u, b = key
        bm = E.MATS[b]
        pre = E.D @ bm if typ == 'y' else bm
        rl = []
        for (l2, g2, cf) in P.gens[u]:
            g3 = E.intern(E.snap(pre @ E.MATS[g2]))
            if l2 == target_l and E.layer(g3) == M and Sy.frame(g3) == Fstar:
                kterm.setdefault(key, []).append(((l2, g3), cf))
            if E.layer(g3) < dmin or free(l2, g3):
                continue
            rl.append(((l2, g3), cf))
        colrows[key] = rl; stats['created'] += 1
        for (r, cf) in rl:
            if r not in rowprod:
                rowprod[r] = prod_keys(r) - dead
                queue.append(r)
        # forced zeros from this column's rows
        for (r, cf) in rl:
            s = rowprod[r]
            if len(s) == 1:
                kill(next(iter(s)), r)
                break

    for (u0, b0) in seeds:
        create(('y', u0, b0))
    while queue:
        r = queue.pop()
        s = rowprod.get(r)
        if not s or len(s) < 2:
            continue
        # expand only rows that still lie in some live created column
        if not any((k in colrows) for k in s):
            continue
        for k in list(s):
            if k not in colrows and k not in dead:
                if stats['created'] >= cap:
                    stats['cap'] = True
                    return finish(colrows, kterm, rowprod, dead, seeds, stats)
                create(k)
    return finish(colrows, kterm, rowprod, dead, seeds, stats)


def finish(colrows, kterm, rowprod, dead, seeds, stats):
    live = [k for k in colrows if k not in dead]
    cols = {k: i for i, k in enumerate(live)}
    rows = {}; ent = []
    for k in live:
        for (r, cf) in colrows[k]:
            if r not in rows:
                rows[r] = len(rows)
            ent.append((rows[r], cols[k], cf))
    klass = [cols[k] for k in live if k in kterm]
    kterms = {cols[k]: kterm[k] for k in live if k in kterm}
    stats['live'] = len(live); stats['rows'] = len(rows); stats['class_live'] = len(klass)
    stats['seeds_live'] = sum(1 for (u0, b0) in seeds if ('y', u0, b0) not in dead)
    stats['class_total'] = len(kterm)
    return rows, cols, ent, (klass, kterms), stats
