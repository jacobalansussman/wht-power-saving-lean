"""x-exit-words: the pad-minimal ("ceiling") dedicated-bank table for a tier-2 (exit-frame) word.
Method = x-pack's pack_exact (cyclic exact-73 groups + a leftover multiset L whose 72 copies fill full banks and
ONE partial bank, rows by a Koenig edge-colouring), generalised to any residual list, with a memoised
pattern-count search for L.  Helper order = Shea's complex588 convention: the 112 source copies (width 9) first,
then the 504 word helpers in canonical order.
Usage: pack_ceil.py CANONICAL-WORD.json.gz OUT.schedule.json.gz [--profile JSON]  (profile override, e.g. U8)
k-word change (13 Oct): [--copies N] keeps only N width-9 source copies (0: mutable sources; 5: the centre copies)."""
import sys, os, json, gzip, itertools
from functools import lru_cache
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tight_pack import patterns
from pack_exact import split_full, edge_colour, validate, M, CALLS
from exit_price import read, word_data, IDLE
from route_price import price

def leftover_banks(Lc, nfull, last, budget=20000):
    ws = sorted(Lc, reverse=True)
    need = tuple(CALLS * Lc[w] for w in ws)
    pf = patterns(ws, M); pl = patterns(ws, last) if last else [tuple(0 for _ in ws)]
    for q in pl:
        rest = tuple(need[i] - q[i] for i in range(len(ws)))
        if min(rest) < 0: continue
        nodes = [0]
        @lru_cache(None)
        def dfs(k, start, rem):
            nodes[0] += 1
            if nodes[0] > budget: return None
            if k == 0: return () if not any(rem) else None
            if sum(r * w for r, w in zip(rem, ws)) != k * M: return None
            for idx in range(start, len(pf)):
                p = pf[idx]
                if all(p[i] <= rem[i] for i in range(len(ws))):
                    s = dfs(k - 1, idx, tuple(rem[i] - p[i] for i in range(len(ws))))
                    if s is not None: return (p,) + s
            return None
        sol = dfs(nfull, 0, rest)
        if sol is not None:
            return [Counter({ws[i]: p[i] for i in range(len(ws)) if p[i]}) for p in list(sol) + ([q] if last else [])]
    return None

def split_prop(cnt, nbins):
    """partition the multiset cnt (width->count, total width = 73 nbins) into nbins groups of width exactly 73:
    each step takes the exact-73 pattern closest to the remaining average composition (the last group is the
    remainder, automatically of width 73)."""
    sizes = sorted(cnt, reverse=True); pats = patterns(sizes, M)
    rem = [cnt[s] for s in sizes]; chosen = []
    for k in range(nbins, 1, -1):
        avg = [r / k for r in rem]
        feas = [q for q in pats if all(q[i] <= rem[i] for i in range(len(sizes)))]
        if not feas: return None
        q = min(feas, key=lambda q: sum((q[i] - avg[i]) ** 2 for i in range(len(sizes))))
        rem = [rem[i] - q[i] for i in range(len(sizes))]; chosen.append(q)
    assert sum(rem[i] * sizes[i] for i in range(len(sizes))) == M
    chosen.append(tuple(rem))
    assert all(sum(q) <= CALLS for q in chosen)
    return sizes, chosen

def build(res, L, seed=0):
    cnt = Counter(res); Lc = Counter(L); rest = cnt - Lc
    assert sum(w * n for w, n in rest.items()) % M == 0
    nbins = sum(w * n for w, n in rest.items()) // M
    tot = CALLS * sum(L); nfull, last = divmod(tot, M)
    pats = leftover_banks(Lc, nfull, last)
    if pats is None: return None
    sp = split_prop(dict(rest), nbins)
    if sp is None: return None
    sizes, chosen = sp
    pools = defaultdict(list)
    for q, w in enumerate(res): pools[w].append(q)
    groups = []
    for p in chosen:
        g = []
        for i, s in enumerate(sizes):
            for _ in range(p[i]): g.append(pools[s].pop())
        groups.append(g)
    Lh = [pools[w].pop() for w in L]
    assert all(not v for v in pools.values())
    rows = [[None] * len(res) for _ in range(CALLS)]; pads = []
    for b, g in enumerate(groups):
        assert sum(res[q] for q in g) == M and len(g) <= CALLS
        pads.extend([0] * CALLS)
        for t in range(CALLS):
            for j, q in enumerate(g): rows[t][q] = b * CALLS + (t + j) % CALLS
    base = len(pads)
    hw = defaultdict(list)
    for q in Lh: hw[res[q]].append(q)
    quota = {q: CALLS for q in Lh}; banks_items = []
    for pat in pats:
        items = []
        for w, n in sorted(pat.items()):
            for _ in range(n):
                q = max(hw[w], key=lambda x: quota[x]); quota[q] -= 1; items.append(q)
        banks_items.append(items)
    assert all(v == 0 for v in quota.values())
    edges, col = edge_colour(banks_items, Lh)
    for (bk, q), c in zip(edges, col):
        assert rows[c][q] is None; rows[c][q] = base + bk
    for pat in pats: pads.append(M - sum(w * n for w, n in pat.items()))
    return rows, pads, dict(L=list(L), leftover_patterns=[dict(p) for p in pats], n_groups=len(groups))

def main():
    wp, out = sys.argv[1], sys.argv[2]
    copies = 112
    if wp.endswith('.info.json'):      # natural-exit gcert info: residuals of ALL slots (copies included), full histogram
        d = json.load(open(wp)); res = list(d['residuals']); prof = {int(k): n for k, n in d['histogram'].items()}; copies = 0
        # Shea's order: the source copies are slots 504..615 of the lowering; put them first
        res = res[504:] + res[:504]
    else:
        c = read(wp); prof, es, us = word_data(c)
        if '--copies' in sys.argv: copies = int(sys.argv[sys.argv.index('--copies') + 1])   # k-word: 0 = mutable sources, 5 = centre copies only
        res = [9] * copies + [u - e for e, u in zip(es, us)]
    if '--profile' in sys.argv: prof = {int(k): n for k, n in json.loads(sys.argv[sys.argv.index('--profile') + 1]).items()}
    S = sum(res); B = -(-CALLS * S // M); P = B * M - CALLS * S; lw = S % M
    print('residuals', sorted(Counter(res).items()), 'S', S, 'banks', B, 'pad', P, 'leftover width', lw, flush=True)
    ws = sorted(set(res), reverse=True); cnt = Counter(res)
    cands = set()
    for k in range(1, 11):
        for comb in itertools.combinations_with_replacement(ws, k):
            if sum(comb) == lw and all(Counter(comb)[w] <= cnt[w] for w in comb):
                cands.add(tuple(sorted(comb, reverse=True)))
    # prefer L with many distinct widths (more bank patterns), then short
    order = sorted(cands, key=lambda c: (abs(len(set(c)) - 3), len(c), c))
    if '--L' in sys.argv: order = [tuple(json.loads(sys.argv[sys.argv.index('--L') + 1]))] + order
    for L in order:
        import time; t0 = time.time()
        r = build(res, list(L))
        print('L', L, 'ok' if r else 'no', '%.1fs' % (time.time() - t0), flush=True)
        if not r: continue
        rows, pads, info = r
        validate(rows, pads, res)
        assert len(pads) == B and sorted(p for p in pads if p) == ([P] if P else []), (len(pads), B)
        hist = Counter({rr: CALLS * n for rr, n in prof.items()})
        hist[1] += CALLS * copies; hist[8] += CALLS * copies
        for rr, n in IDLE.items(): hist[rr] += n * 112
        hist.update(p for p in pads if p)
        W = 64 * 112 + len(pads); D = M * W - sum(rr * n for rr, n in hist.items())
        assert D == 2560
        pr = price(dict(hist=dict(sorted(hist.items())), W=W, D=D))
        print('realised table: banks', len(pads), 'pads', sorted(Counter(pads).items()), 'W', W, 'root', pr['root_lower'], 'kernel', pr['kernel_A'], 'dft', pr['dft'], flush=True)
        open(out, 'wb').write(gzip.compress(json.dumps(dict(rows=rows, pads=pads, widths=res), separators=(',', ':')).encode(), mtime=0))
        json.dump(dict(word=wp, profile=prof, hist=dict(sorted(hist.items())), W=W, D=D, banks=len(pads), **pr, **info),
                  open(out + '.info.json', 'w'), indent=1)
        return
    raise SystemExit('no leftover multiset worked')

if __name__ == '__main__':
    main()
