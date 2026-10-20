"""relay.py (g-recombine) -- breaks H5/H17: targets that are already served act as side information for later
targets (successive decoding through target registers, 'relay').  Decoder level only (frames for the relay are
NOT modelled).  Pattern with an order: B[S,S] = 1, B[S,T] = 0 when S does not know T AND T is not served before S.
Exact minrank over F_2 by brute force on small patterns: the C5 side-information graph (scalar minrank 3, clique 2),
and its best order.  usage: python3 -I relay.py"""
import itertools
def rank2(rows):
    rows = list(rows); r = 0
    for bit in reversed(range(16)):
        piv = next((i for i in range(r, len(rows)) if rows[i] >> bit & 1), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> bit & 1: rows[i] ^= rows[r]
        r += 1
    return r
def minrank(v, free):
    cells = [(s, t) for s in range(v) for t in range(v) if s != t and free[s][t]]
    best = v
    for bits in range(1 << len(cells)):
        rows = [1 << s for s in range(v)]
        for k, (s, t) in enumerate(cells):
            if bits >> k & 1: rows[s] |= 1 << t
        best = min(best, rank2(rows))
    return best
def knows_cycle(v):
    return [[(t - s) % v in (1, v - 1) for t in range(v)] for s in range(v)]
for v in (5,):
    K = knows_cycle(v)
    print('C%d: scalar minrank (static side information):' % v, minrank(v, K))
    res = {}
    for order in ([0] + list(p) for p in itertools.permutations(range(1, v))):
        pos = {x: i for i, x in enumerate(order)}
        free = [[K[s][t] or pos[t] < pos[s] for t in range(v)] for s in range(v)]
        res.setdefault(minrank(v, free), order)
        if v == 7: break
    print('   with relays along an order: minrank values reached', sorted(res), ' e.g. order', res[min(res)])
