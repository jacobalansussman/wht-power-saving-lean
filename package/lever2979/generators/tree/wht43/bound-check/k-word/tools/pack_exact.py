"""x-pack: an exact (pad-minimal) bank table for Shea's C588 certificate (PR588 atm3 + 112 source copies).

Ceiling argument: the 72 rows carry 72 * 5068 = 364896 residual width; a bank holds at most 73, so at least
ceil(364896/73) = 4999 banks, total pad >= 4999*73 - 364896 = 31; the pad price p((73/p)^a - 1) is concave
with value 0 at p = 0, so ONE pad of rank 31 (all other banks full) is the cheapest possible pad multiset.

Construction (reaches the ceiling):
  * every row: the same 69 helper groups of total width exactly 73 (found by a randomised pattern split, as in
    tight_pack.split_exact), placed cyclically (group g, row t, member j -> bank 72 g + (t + j) mod 72), so every
    one of the 69*72 = 4968 banks gets load exactly 73;
  * the remaining helpers L (width 31 per row, 72 rows -> 2232 = 30*73 + 42) go to 31 more banks: 30 full banks
    and one of load 42 (pad 31).  The bank contents are integer patterns over L's widths; the row assignment is a
    proper 72-edge-colouring of the bipartite multigraph banks x helpers (every helper has degree 72, every bank
    degree <= 18), which exists by Koenig's theorem and is found by alternating-path recolouring.
Output: schedule JSON (rows, pads, widths) in bank_reprice.py's format, widths in Shea's helper order (schedRes).
"""
import sys, os, re, json, gzip, random, itertools
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(__file__))
from tight_pack import patterns

M, CALLS = 73, 72


def read_res(sched_lean):
    src = open(sched_lean).read()
    res = [int(x) for x in re.search(r'schedRes : List Nat := \[([0-9, ]*)\]', src).group(1).split(',')]
    assert len(res) == 616
    return res


def split_full(cnt, nbins, seed, trials=20000):
    """nbins groups of width exactly M from the multiset cnt (dict width->count) using all of it."""
    sizes = sorted(cnt, reverse=True)
    pats = patterns(sizes, M)
    rng = random.Random(seed)
    for _ in range(trials):
        rem = [cnt[s] for s in sizes]
        chosen = []
        ok = True
        for b in range(nbins):
            feas = [p for p in pats if all(p[i] <= rem[i] for i in range(len(sizes)))]
            if not feas:
                ok = False
                break
            def score(p):
                return sum(p[i] / (rem[i] + 1e-9) for i in range(len(sizes)) if p[i]) + rng.random() * 0.5
            p = max(feas, key=score) if rng.random() < 0.7 else rng.choice(feas)
            for i in range(len(sizes)):
                rem[i] -= p[i]
            chosen.append(p)
        if ok and not any(rem):
            return sizes, chosen
    return None


def leftover_patterns(Lc, nfull, last):
    """bank contents for the 72 copies of the leftover L (Counter width->multiplicity): nfull banks of load M and
    one of load `last`; returns list of Counters or None"""
    ws = sorted(Lc, reverse=True)
    need = [CALLS * Lc[w] for w in ws]
    pf = patterns(ws, M)
    pl = patterns(ws, last)
    for q in pl:
        rest = [need[i] - q[i] for i in range(len(ws))]
        if min(rest) < 0:
            continue
        # DFS: nfull patterns from pf (non-increasing index) summing to rest
        sol = []
        def dfs(k, start, rem):
            if k == 0:
                return not any(rem)
            for idx in range(start, len(pf)):
                p = pf[idx]
                if all(p[i] <= rem[i] for i in range(len(ws))):
                    # bound: remaining width must equal k*M
                    sol.append(p)
                    if dfs(k - 1, idx, [rem[i] - p[i] for i in range(len(ws))]):
                        return True
                    sol.pop()
            return False
        if dfs(nfull, 0, rest):
            return [Counter({ws[i]: p[i] for i in range(len(ws)) if p[i]}) for p in sol + [q]]
    return None


def edge_colour(banks_items, helpers, ncol=CALLS):
    """banks_items: list over banks of list of helper ids (multiset = multi-edges).  Returns colour per edge:
    dict (bank, k) -> colour, proper at banks and at helpers.  Alternating-path (Kempe) algorithm."""
    edges = [(b, h) for b, items in enumerate(banks_items) for h in items]
    atB = [dict() for _ in banks_items]          # colour -> edge index at bank
    atH = defaultdict(dict)                       # helper -> colour -> edge index
    col = [None] * len(edges)
    for e, (b, h) in enumerate(edges):
        fb = next(c for c in range(ncol) if c not in atB[b])
        fh = next(c for c in range(ncol) if c not in atH[h])
        if fb not in atH[h]:
            c = fb
        elif fh not in atB[b]:
            c = fh
        else:
            # flip the fb/fh alternating path starting at helper h with colour fb
            path = []
            side, node, cc = 'H', h, fb
            while True:
                tbl = atH[node] if side == 'H' else atB[node]
                if cc not in tbl:
                    break
                f = tbl[cc]
                path.append(f)
                fb2, fh2 = edges[f]
                if side == 'H':
                    side, node = 'B', fb2
                else:
                    side, node = 'H', fh2
                cc = fh if cc == fb else fb
            for f in path:
                bb, hh = edges[f]
                del atB[bb][col[f]]
                del atH[hh][col[f]]
            for f in path:
                bb, hh = edges[f]
                col[f] = fh if col[f] == fb else fb
                atB[bb][col[f]] = f
                atH[hh][col[f]] = f
            assert fb not in atH[h] and fb not in atB[b]
            c = fb
        col[e] = c
        atB[b][c] = e
        atH[h][c] = e
    return edges, col


def build(res, L_widths, seed=0):
    cnt = Counter(res)
    Lc = Counter(L_widths)
    rest = cnt - Lc
    assert sum(w * n for w, n in rest.items()) % M == 0
    nbins = sum(w * n for w, n in rest.items()) // M
    lw = sum(L_widths)
    tot = CALLS * lw
    nfull = tot // M
    last = tot - nfull * M
    pats = leftover_patterns(Lc, nfull, last)
    if pats is None:
        return None
    sp = split_full(dict(rest), nbins, seed)
    if sp is None:
        return None
    sizes, chosen = sp
    pools = defaultdict(list)
    for q, w in enumerate(res):
        pools[w].append(q)
    groups = []
    for p in chosen:
        g = []
        for i, s in enumerate(sizes):
            for _ in range(p[i]):
                g.append(pools[s].pop())
        groups.append(g)
    Lh = []
    for w in L_widths:
        Lh.append(pools[w].pop())
    assert all(not v for v in pools.values())
    rows = [[None] * len(res) for _ in range(CALLS)]
    pads = []
    for b, g in enumerate(groups):
        assert sum(res[q] for q in g) == M and len(g) <= CALLS
        pads.extend([0] * CALLS)
        for t in range(CALLS):
            for j, q in enumerate(g):
                rows[t][q] = b * CALLS + (t + j) % CALLS
    base = len(pads)
    # leftover banks: distribute the width-w items of each pattern over the helpers of width w (72 each)
    hw = defaultdict(list)
    for q in Lh:
        hw[res[q]].append(q)
    quota = {q: CALLS for q in Lh}
    banks_items = []
    for pat in pats:
        items = []
        for w, n in sorted(pat.items()):
            for _ in range(n):
                q = max(hw[w], key=lambda x: quota[x])   # spread evenly
                quota[q] -= 1
                items.append(q)
        banks_items.append(items)
    assert all(v == 0 for v in quota.values())
    edges, col = edge_colour(banks_items, Lh)
    for (bk, q), c in zip(edges, col):
        assert rows[c][q] is None
        rows[c][q] = base + bk
    for pat in pats:
        pads.append(M - sum(w * n for w, n in pat.items()))
    return rows, pads, dict(groups=groups, leftover_helpers=Lh, leftover_patterns=[dict(p) for p in pats])


def validate(rows, pads, widths):
    """Shea's bank_reprice.validate, plus his two mutation controls"""
    import shea_bank_reprice as sb
    sb.validate(rows, pads, widths)
    bad = json.loads(json.dumps(rows)); bad[0][1] = bad[0][0]
    try:
        sb.validate(bad, pads, widths)
    except AssertionError:
        pass
    else:
        raise AssertionError('duplicate-owner control accepted')
    bp = pads.copy(); bp[0] += 1
    try:
        sb.validate(rows, bp, widths)
    except AssertionError:
        pass
    else:
        raise AssertionError('padding control accepted')
    return True


if __name__ == '__main__':
    sched_lean, out = sys.argv[1], sys.argv[2]
    res = read_res(sched_lean)
    print('widths', sorted(Counter(res).items()), 'sum', sum(res), '72*sum', CALLS * sum(res),
          'ceil banks', -(-CALLS * sum(res) // M), 'min pad', -(-CALLS * sum(res) // M) * M - CALLS * sum(res))
    cands = set()
    ws = sorted(set(res), reverse=True)
    for k in range(2, 8):
        for comb in itertools.combinations_with_replacement(ws, k):
            if sum(comb) == 31 and all(Counter(comb)[w] <= Counter(res)[w] for w in comb):
                cands.add(tuple(sorted(comb, reverse=True)))
    for L in sorted(cands, key=lambda c: (len(c), c)):
        r = build(res, list(L))
        print('L', L, 'ok' if r else 'infeasible', flush=True)
        if r:
            rows, pads, info = r
            validate(rows, pads, res)
            print('banks', len(pads), 'pads', sorted(Counter(pads).items()), 'leftover patterns', info['leftover_patterns'])
            open(out, 'wb').write(gzip.compress(json.dumps(dict(rows=rows, pads=pads, widths=res), separators=(',', ':')).encode(), mtime=0))
            json.dump(dict(L=L, **info), open(out + '.info.json', 'w'))
            break
