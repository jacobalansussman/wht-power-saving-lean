"""peel2.py (w8-b5sparse): peel.py's exact forced-zero closure, made COMPLETE, plus a valid truncated system at the cap.

peel.py (b5-twobranch, used by w7-b5 and v-w7-b5) has two gaps.  (1) A kill of a producer that was never CREATED is not
removed from the producer sets of the rows it would produce (peel.py discards via colrows only), so later single-producer
rows are missed: weaker, not invalid.  peel2 keeps a reverse index key -> rows.  (2) A completeness gap: a row r is popped from the queue; if at that
moment every created producer of r is dead, r is skipped and never re-queued.  If later another (live) producer of r is
created through a different row, r becomes a row of a live column but its remaining producers are never created, so the
assembled system has r with an incomplete producer set (too strong: a possible false 'fail').  peel2 re-queues such rows
(the `pending` set), and at the end CHECKS that every row of every live column has all its non-dead producers created.

At the cap, peel2 returns the live columns and the set of COMPLETE rows (rows whose non-dead producers are all created).
Restricting the identity to complete rows only is a VALID necessary test (every row in it has all its producers; killed
columns are exact forced zeros), whatever the cap; it is weaker than the full closure.

peel_closure2(P, target_l, dmin, seeds, cap, fifo) -> dict with rows, cols, ent, kl (endsolve format, rows = complete
rows only when capped, all rows of live columns otherwise), stats (incl. 'incomplete_rows' = rows of live columns with an
uncreated live producer, which must be 0 when the closure finished)."""
from collections import deque
import numpy as np
from sp2 import M, Fstar
DEBUG = {}
RECORD = []
EVENTS = []


def peel_closure2(P, target_l, dmin, seeds, cap=20000, fifo=True, verbose_every=0):
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

    colrows = {}; kterm = {}; rowprod = {}; dead = set(); expanded = set(); pending = set()
    del EVENTS[:]
    keyrows = {}        # producer key -> rows whose producer set contains it (so kills of UNCREATED keys propagate)
    queue = deque()
    stats = dict(created=0, killed=0, triggers=[], requeued=0)

    def kill(key, why=None):
        stack = [(key, why)]
        while stack:
            k, w = stack.pop()
            if k in dead:
                continue
            dead.add(k); stats['killed'] += 1
            if RECORD:
                EVENTS.append((k, w))
            if w is not None and len(stats['triggers']) < 5000:
                stats['triggers'].append((w, k))
            for r in keyrows.get(k, ()):
                s = rowprod.get(r)
                if s is None or k not in s:
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
                for kk in rowprod[r]:
                    keyrows.setdefault(kk, []).append(r)
                queue.append(r)
            elif r in pending:
                pending.discard(r); queue.append(r); stats['requeued'] += 1
        for (r, cf) in rl:
            s = rowprod[r]
            if len(s) == 1:
                kill(next(iter(s)), r)
                break

    for (u0, b0) in seeds:
        create(('y', u0, b0))
    capped = False
    while queue:
        r = queue.popleft() if fifo else queue.pop()
        if r in expanded:
            continue
        s = rowprod.get(r)
        if not s:
            continue
        if not any((k in colrows) for k in s):
            pending.add(r)
            continue
        for k in list(s):
            if k not in colrows and k not in dead:
                if stats['created'] >= cap:
                    capped = True; break
                create(k)
        if capped:
            break
        expanded.add(r)
        if verbose_every and stats['created'] >= stats.get('next_report', verbose_every):
            stats['next_report'] = stats['created'] + verbose_every
            ncl = sum(1 for k in kterm if k not in dead)
            print("       progress: created %d, killed %d, live %d, rows registered %d, expanded %d, pending %d, queue %d, class live %d"
                  % (stats['created'], stats['killed'], stats['created'] - stats['killed'], len(rowprod), len(expanded),
                     len(pending), len(queue), ncl), flush=True)
    stats['cap'] = capped
    DEBUG['state'] = (colrows, rowprod, dead, expanded, pending)
    return finish2(colrows, kterm, rowprod, dead, seeds, stats)


def finish2(colrows, kterm, rowprod, dead, seeds, stats):
    live = [k for k in colrows if k not in dead]
    liveset = set(live)
    cols = {k: i for i, k in enumerate(live)}
    rows = {}; ent = []
    incomplete = set(); allrows = set()
    for k in live:
        for (r, cf) in colrows[k]:
            allrows.add(r)
            s = rowprod.get(r, set())
            if any((kk not in colrows and kk not in dead) for kk in s):
                incomplete.add(r)
    for k in live:
        for (r, cf) in colrows[k]:
            if r in incomplete:
                continue
            if r not in rows:
                rows[r] = len(rows)
            ent.append((rows[r], cols[k], cf))
    klass = [cols[k] for k in live if k in kterm]
    kterms = {cols[k]: kterm[k] for k in live if k in kterm}
    # columns with NO complete row are unconstrained in the truncated system
    constrained = set(c for (r, c, cf) in ent)
    stats['live'] = len(live); stats['rows_all'] = len(allrows); stats['rows'] = len(rows)
    stats['incomplete_rows'] = len(incomplete); stats['class_live'] = len(klass)
    stats['class_unconstrained'] = sum(1 for c in klass if c not in constrained)
    stats['cols_unconstrained'] = len(live) - len(constrained)
    stats['seeds_live'] = sum(1 for (u0, b0) in seeds if ('y', u0, b0) not in dead)
    stats['class_total'] = len(kterm)
    return rows, cols, ent, (klass, kterms), stats


def audit_old(P, rows, cols, kl, dead_hint=None):
    """for a system produced by the OLD peel (or any closure): count rows whose full producer set (prod_keys, minus the
    columns of the system) contains keys that are neither in the system nor provably dead.  We cannot know the old
    peel's dead set from its output, so we report rows having producers outside the system (an upper bound on the
    incompleteness; dead producers are legitimately absent)."""
    raise NotImplementedError
