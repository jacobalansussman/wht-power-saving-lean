"""fs-price: consolidate the bank pads of a re-packed hosting result (host_repack2 output).  Same bank rules as
host_price.banks / audit_slots: load <= 73, at most one residual of a call in a bank, one pad per bank.  The pad price
g(p) = p ln(73/p) is concave, so for a fixed number of banks the best invoice puts the whole pad total into ONE bank.
Random item moves and swaps between banks, accepting every move that does not raise sum g(pad) (neutral moves let a pad
of 1 wander until it meets another pad and merges).  Re-prices exactly (host_repack.price) and writes the result.
With T0 > 0 (5th argument) it anneals: a move that raises the cost by d is accepted with probability exp(-d/T), T falls
linearly from T0 to 0; the best state seen is kept.
usage: python3 -I pad_sa.py REPACKED.json OUT.json [iterations] [seed] [T0]"""
import sys, os, json, math, random, shutil
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from host_repack import price
M = 73
g = lambda p: p * math.log(M / p) if p > 0 else 0.0
r = json.load(open(sys.argv[1])); out = sys.argv[2]
iters = int(sys.argv[3]) if len(sys.argv) > 3 else 3000000
rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 11)
T0 = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
import copy
banks = [[tuple(x) for x in b] for b in r['bank_items']]
items0 = Counter(x for b in banks for x in b)
nb = len(banks)
load = [sum(w for w, _ in b) for b in banks]
calls = [Counter(k for _, k in b) for b in banks]
cost = sum(g(M - l) for l in load)
tot = sum(M - l for l in load)
best = g(tot)
bcost, bstate = cost, copy.deepcopy(banks)
ok = lambda d, it: d <= 1e-9 or (T0 > 0 and rng.random() < math.exp(-d / max(1e-9, T0 * (1 - it / iters))))
print('start: banks %d, pad total %d, pad cost %.1f (one-pad optimum %.1f)' % (nb, tot, cost, best), flush=True)
padded = lambda: [i for i in range(nb) if load[i] < M]
for it in range(iters):
    if it % 20000 == 0:
        P = padded()
        if len(P) <= 1 and T0 == 0: break
    i = rng.choice(P) if rng.random() < 0.7 else rng.randrange(nb)
    j = rng.randrange(nb)
    if i == j or not banks[i] or not banks[j]: continue
    ai = rng.randrange(len(banks[i])); a = banks[i][ai]
    if rng.random() < 0.3:                       # move a: i -> j
        if load[j] + a[0] > M or calls[j][a[1]]: continue
        d = g(M - load[i] + a[0]) + g(M - load[j] - a[0]) - g(M - load[i]) - g(M - load[j])
        if not ok(d, it) or len(banks[i]) == 1: continue
        banks[i].pop(ai); banks[j].append(a); load[i] -= a[0]; load[j] += a[0]
        calls[i][a[1]] -= 1; calls[j][a[1]] += 1; cost += d
    else:                                        # swap a (in i) with b (in j)
        bj = rng.randrange(len(banks[j])); b = banks[j][bj]
        li, lj = load[i] - a[0] + b[0], load[j] - b[0] + a[0]
        if li > M or lj > M: continue
        if a[1] != b[1] and (calls[i][b[1]] or calls[j][a[1]]): continue
        d = g(M - li) + g(M - lj) - g(M - load[i]) - g(M - load[j])
        if not ok(d, it): continue
        banks[i][ai] = b; banks[j][bj] = a; load[i], load[j] = li, lj
        calls[i][a[1]] -= 1; calls[i][b[1]] += 1; calls[j][b[1]] -= 1; calls[j][a[1]] += 1; cost += d
    if cost < bcost - 1e-9:
        bcost, bstate = cost, copy.deepcopy(banks)
        if bcost < best + 1e-6: break
    if it % 500000 == 0:
        print('  it %d: pad cost %.1f, padded banks %d' % (it, cost, len(padded())), flush=True)
banks = bstate
assert all(sum(w for w, _ in b) <= M and len({k for _, k in b}) == len(b) for b in banks)
assert Counter(x for b in banks for x in b) == items0 and len(banks) == nb
pads = [M - sum(w for w, _ in b) for b in banks]
hist = Counter({int(k): v for k, v in r['hist'].items()})
for p, n in Counter({int(k): v for k, v in r['pads'].items()}).items():
    if p: hist[p] -= n
hist.update(p for p in pads if p); hist = +hist
A0, A3, A4 = price(hist, r['W'])
print('pad_sa: %d banks, pad cost %.1f, pads %s; W=%d dft3=%d dft4=%d' % (nb, sum(g(p) for p in pads), sorted(Counter(pads).items()), r['W'], A3 - 1, A4 - 1))
r.update(pads={str(k): v for k, v in sorted(Counter(pads).items())}, hist={str(k): v for k, v in sorted(hist.items())},
         dft=A3 - 1, dft_quartic=A4 - 1, kernel_A=A3, bank_items=banks, pad_sa_from=sys.argv[1])
json.dump(r, open(out, 'w'), indent=1); shutil.copyfile(sys.argv[1] + '.slots.gz', out + '.slots.gz')
