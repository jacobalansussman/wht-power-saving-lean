"""s-stack: second re-packer for the unhosted residuals (same rules as host_price.banks: <= one residual per call per
bank, load <= 73, one pad per bank).  For a target bank count B (from max-items-per-call upwards) banks are filled one
at a time; a call whose remaining items equal the remaining banks MUST place one item in the current bank (forced),
then the bank is filled to the largest load <= 73 (exact DP over calls, forced calls mandatory, others optional,
preferring calls with many items left).  Best (lowest first-order pad cost) feasible packing over B and seeds wins.
Usage: python3 -I host_repack2.py words.json result.json out.json"""
import sys, os, json, gzip, random, math
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from host_repack import items_of, price
M = 73


def attempt(items, B, rng, noise):
    left = defaultdict(Counter)
    for w, k in items: left[k][w] += 1
    banks = []
    for i in range(B):
        R = B - i
        cnt = {k: sum(c.values()) for k, c in left.items() if sum(c.values())}
        if not cnt: break
        if max(cnt.values()) > R: return None
        forced = [k for k in cnt if cnt[k] == R]
        opt = sorted((k for k in cnt if cnt[k] < R), key=lambda k: (-cnt[k] - noise * rng.random(), k))
        # DP over forced calls (mandatory), then optional calls; state load -> choice
        dp = {0: []}
        for k in forced:
            nd = {}
            for load, ch in dp.items():
                for w in sorted(left[k], key=lambda w: -w - noise * rng.random()):
                    if left[k][w] and load + w <= M and (load + w) not in nd:
                        nd[load + w] = ch + [(w, k)]
            dp = nd
            if not dp: return None
        for k in opt:
            nd = dict(dp)
            for load, ch in dp.items():
                for w in left[k]:
                    if left[k][w] and load + w <= M and (load + w) not in nd:
                        nd[load + w] = ch + [(w, k)]
            dp = nd
        best = dp[max(dp)]
        for w, k in best: left[k][w] -= 1
        banks.append(best)
    if any(sum(c.values()) for c in left.values()): return None
    return banks


if __name__ == '__main__':
    words = json.load(open(sys.argv[1])); path = sys.argv[2]; out = sys.argv[3]
    r = json.load(open(path)); r['_path'] = path
    items = items_of(words[r['name']], r)
    per = Counter(k for _, k in items); LB = max(max(per.values()), -(-sum(w for w, _ in items) // M))
    rng = random.Random(7); best = None
    for B in range(LB, r['banks'] + 1):
        for s in range(12):
            banks = attempt(items, B, rng, 0 if s == 0 else 2.0)
            if banks is None: continue
            assert all(sum(w for w, _ in b) <= M and len({k for _, k in b}) == len(b) for b in banks)
            assert Counter(x for b in banks for x in b) == Counter(items)
            pads = [M - sum(w for w, _ in b) for b in banks]
            cost = sum(p * math.log(M / p) for p in pads if p)
            if best is None or cost < best[0]: best = (cost, banks, pads, B)
        if best and best[3] == B and B > LB + 6: break
    cost, banks, pads, B = best
    old = Counter({int(k): v for k, v in r['pads'].items()})
    hist = Counter({int(k): v for k, v in r['hist'].items()})
    for p, n in old.items():
        if p: hist[p] -= n
    hist.update(p for p in pads if p); hist = +hist
    W = r['W'] - r['banks'] + len(banks)
    A0, A3, A4 = price(hist, W)
    print(f"{path}: LB {LB} old banks {r['banks']} -> {len(banks)} banks, pad L {cost:.0f}, pads {sorted(Counter(pads).items())}; W={W} dft3={A3-1} dft4={A4-1}", flush=True)
    r2 = {k: v for k, v in r.items() if k != '_path'}
    r2.update(pads={str(k): v for k, v in sorted(Counter(pads).items())}, hist={str(k): v for k, v in sorted(hist.items())}, W=W,
              banks=len(banks), dft=A3 - 1, dft_quartic=A4 - 1, kernel_A=A3, repacked_from=path, bank_items=banks)
    json.dump(r2, open(out, 'w'), indent=1)
    import shutil; shutil.copyfile(path + '.slots.gz', out + '.slots.gz')
