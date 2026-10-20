"""s-stack: tier 1 on top of tier 3b.  x-hosting's host_price.banks() puts the unhosted helper residuals into banks by
first-fit decreasing (<= one residual per call per bank, load <= 73, one pad block per bank).  Here the same items
are re-packed: banks are filled one at a time with an EXACT load-73 set of items from distinct calls when one exists
(DP over calls, preferring the calls with the most items left, since #banks >= max items of one call), the rest by
best fit.  Same rules as host_price.banks / audit_slots (one residual per call per bank, load <= 73).  Then the
invoice is re-priced exactly (x-pack/d-lean-route tools: rational moment root, kernel grid point n=3 and n=4).
Usage: python3 -I host_repack.py words.json result.json [seeds]"""
import sys, os, json, gzip, random
from collections import Counter, defaultdict
from fractions import Fraction as F
XP = os.path.join(os.environ.get('BC', '.'), 'x-pack/tools')  # as run: the agents' working folder; see generators/README.txt
sys.path.insert(0, XP)
from exact_moment import certify
import foldlib as fl
from rate_gap import best_A
M = 73


def leftovers(words, r):
    w = words[r['name']] if r['name'] in words else None
    return w


def items_of(word, r):
    wid = Counter({int(k): n for k, n in word['widths'].items()})
    if r['copies']: wid[9] += 112
    slots = json.load(gzip.open(r['_path'] + '.slots.gz', 'rt'))
    used = Counter((s[0], s[1]) for s in slots)
    items = []
    for k in range(72):
        for ww, n in wid.items():
            items += [(ww, k)] * (n - used[(k, ww)])
    return items


def pack(items, rng, noise):
    left = defaultdict(Counter)              # call -> Counter(width)
    for w, k in items: left[k][w] += 1
    banks = []
    while any(sum(c.values()) for c in left.values()):
        calls = sorted((k for k in left if sum(left[k].values())), key=lambda k: (-sum(left[k].values()) - noise * rng.random(), k))
        # DP over calls in priority order: reachable loads -> chosen list (first found = higher-priority calls used)
        dp = {0: []}
        for k in calls:
            nd = dict(dp)
            for load, ch in dp.items():
                for w in left[k]:
                    if left[k][w] and load + w <= M and (load + w) not in nd:
                        nd[load + w] = ch + [(w, k)]
            dp = nd
        best = dp[max(dp)]
        for w, k in best: left[k][w] -= 1
        banks.append(best)
    return banks


def price(hist, W):
    D = M * W - sum(k * v for k, v in hist.items()); assert D == 2560, D
    rate = certify(dict(hist), M, W)
    A0 = int(F(rate['root_lower']) * 10 ** 10)
    u = dict(m=M, W=W, D=D, H=dict(hist))
    A3 = A0
    while fl.plan_block(u, A3, 10, 40) is None: A3 -= 1
    return A0, A3, best_A(u, A0, 4)


if __name__ == '__main__':
    words = json.load(open(sys.argv[1])); path = sys.argv[2]; seeds = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    r = json.load(open(path)); r['_path'] = path
    word = words[r['name']]
    items = items_of(word, r)
    per = Counter(k for _, k in items); tot = sum(w for w, _ in items)
    print(path, 'items', len(items), 'width', tot, 'max per call', max(per.values()), 'ceil', -(-tot // M), 'old banks', r['banks'], flush=True)
    old = Counter({int(k): v for k, v in r['pads'].items()})
    best = None
    rng = random.Random(1)
    for s in range(seeds):
        banks = pack(items, rng, 0 if s == 0 else 3.0)
        # checks: each bank load <= 73, distinct calls, all items used once
        assert all(sum(w for w, _ in b) <= M and len({k for _, k in b}) == len(b) for b in banks)
        assert Counter(x for b in banks for x in b) == Counter(items)
        pads = [M - sum(w for w, _ in b) for b in banks]
        key = (len(banks), sum(p * 0 for p in pads))
        import math
        cost = sum(p * math.log(M / p) for p in pads if p)
        if best is None or cost < best[0]:
            best = (cost, banks, pads)
    cost, banks, pads = best
    hist = Counter({int(k): v for k, v in r['hist'].items()})
    for p, n in old.items():
        if p: hist[p] -= n
    hist.update(p for p in pads if p)
    hist = +hist
    W = r['W'] - r['banks'] + len(banks)
    A0, A3, A4 = price(hist, W)
    import math
    oldcost = sum(n * p * math.log(M / p) for p, n in old.items() if p)
    print(f"  old: banks {r['banks']} pad L {oldcost:.0f} dft3 {r['dft']}   new: banks {len(banks)} pad L {cost:.0f} pads {sorted(Counter(pads).items())}", flush=True)
    print(f"  new price: W={W} exact_grid={A0} dft3={A3-1} dft4={A4-1}", flush=True)
    out = path.replace('.json', '') if len(sys.argv) <= 4 else sys.argv[4]
    json.dump(dict(source=path, banks=[[list(x) for x in b] for b in banks], pads=pads, W=W, hist=dict(sorted(hist.items())),
                   exact_grid=A0, dft3=A3 - 1, dft4=A4 - 1), open(sys.argv[4] if len(sys.argv) > 4 else '/dev/null', 'w'))
