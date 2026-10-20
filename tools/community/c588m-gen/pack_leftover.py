#!/usr/bin/env python3
"""m-gen (13 Oct 2026): an OPTIONAL upgrade of the annealed bank table: the leftover bin, packed across calls.

The annealed table (x-mutable / d-lean-route sa_pack, K = 56) has 55 bins of load exactly 73 and one leftover bin of
load 45 (helpers of widths 4, 8, 8, 8, 8, 9), replicated over the 72 calls into 72 banks with pad 28 each.
Lean's bank checker GX.schedK (Work/Bridge401/SchedK.lean) asks only: in every call (row) each helper sits in its
own bank (injective rows), and every bank's total load over all calls is 73 - pad.  So the 72 x 6 leftover
(helper, call) items may share banks ACROSS calls (d-lean-route's tight_pack note).  Here they are packed into
44 banks of load 73 and one of load 28 (pad 45): 32 banks {9, 8 x 8}, 7 banks {9 x 5, 4 x 7}, 5 banks
{9, 8 x 6, 4 x 4}, 1 bank {8 x 2, 4 x 3}; within a bank all items come from different calls.
Banks: 55 x 72 + 45 = 4,005; W = 11,173.  The 55 full bins keep their cyclic rows (bank b*72 + (t+j) mod 72).
usage: pack_leftover.py <price-annealed.json> <out dir> [seed]"""
import sys, os, json, gzip, random, hashlib
from collections import Counter
from fractions import Fraction as Q

HERE = os.path.dirname(os.path.abspath(__file__))


def main(price, out, seed=0):
    sys.path.insert(0, HERE)
    from make_sched import schedK
    tk = os.path.join(HERE, "..", "toolkit")
    sys.path.insert(0, os.path.join(tk, "verification")); sys.path.insert(0, os.path.join(tk, "tools", "gen"))
    from exact_moment import certify
    import foldlib as fl
    p = json.load(open(price))
    s = json.loads(gzip.decompress(open(price[:-5] + ".schedule.json.gz", "rb").read()))
    rows, pads, res = s["rows"], s["pads"], s["widths"]
    R = len(res)
    nbin = len(pads) // 72
    left = [b for b in range(nbin) if pads[b * 72] != 0]
    assert len(left) == 1 and all(pads[b * 72] == 0 for b in range(nbin) if b not in left)
    L = left[0]
    G = sorted(q for q in range(R) if rows[0][q] // 72 == L)
    assert sorted(res[q] for q in G) == [4, 8, 8, 8, 8, 9]
    # renumber the full bins 0..54
    full = [b for b in range(nbin) if b != L]
    newbin = {b: i for i, b in enumerate(full)}
    nrow = [[None] * R for _ in range(72)]
    for t in range(72):
        for q in range(R):
            b = rows[t][q] // 72
            if b != L:
                nrow[t][q] = newbin[b] * 72 + rows[t][q] % 72
    base = len(full) * 72
    types = [(1, 8, 0)] * 32 + [(5, 0, 7)] * 7 + [(1, 6, 4)] * 5 + [(0, 2, 3)]
    byw = {9: [q for q in G if res[q] == 9], 8: [q for q in G if res[q] == 8], 4: [q for q in G if res[q] == 4]}
    rng = random.Random(seed)
    for attempt in range(2000):
        left_items = {w: {t: list(byw[w]) for t in range(72)} for w in byw}     # remaining helpers per call
        order = list(range(len(types)))
        rng.shuffle(order)
        assign = {}
        ok = True
        for bi in order:
            n9, n8, n4 = types[bi]
            used = set()
            for w, n in ((9, n9), (8, n8), (4, n4)):
                for _ in range(n):
                    cands = [t for t in range(72) if t not in used and left_items[w][t]]
                    if not cands:
                        ok = False
                        break
                    mx = max(len(left_items[w][t]) for t in cands)
                    t = rng.choice([t for t in cands if len(left_items[w][t]) == mx])
                    q = left_items[w][t].pop()
                    used.add(t)
                    assign[(q, t)] = base + bi
                if not ok:
                    break
            if not ok:
                break
        if ok and len(assign) == 72 * len(G):
            break
    else:
        raise SystemExit("no packing found")
    for (q, t), b in assign.items():
        nrow[t][q] = b
    npads = [0] * base + [73 - sum(w for w in [0]) for _ in types]
    load = Counter()
    for t in range(72):
        for q in range(R):
            load[nrow[t][q]] += res[q]
    npads = [73 - load[b] for b in range(base + len(types))]
    assert all(0 <= x < 73 for x in npads) and sorted(Counter(npads).items()) == [(0, base + 44), (45, 1)], Counter(npads)
    nb = len(npads)
    ok = schedK(nb, R, nrow, npads, res)
    bad = json.loads(json.dumps(nrow)); bad[0][1] = bad[0][0]
    badp = list(npads); badp[0] += 1
    local = {int(r): n for r, n in p["local_profile"].items()}
    hist = Counter({r: 72 * n for r, n in local.items()})
    for r, n in {48: 32, 32: 32, 16: 32, 8: 48}.items():
        hist[r] += n * 112
    hist.update(x for x in npads if x)
    W = 64 * 112 + nb
    assert 73 * W - sum(r * n for r, n in hist.items()) == 2560
    rate = certify(dict(hist), 73, W)
    u = dict(m=73, W=W, D=2560, H=dict(hist))
    A = int(Q(rate["root_lower"]) * 10 ** 10)
    while fl.plan_block(u, A, 10, 40) is None:
        A -= 1
    sp = os.path.join(out, "price-packed.schedule.json.gz")
    open(sp, "wb").write(gzip.compress(json.dumps(dict(rows=nrow, pads=npads, widths=res), separators=(",", ":")).encode(), mtime=0))
    rec = dict(p, schedule="packed", schedule_sha256=hashlib.sha256(open(sp, "rb").read()).hexdigest(), banks=nb, W=W,
               profile={str(r): n for r, n in sorted(hist.items())},
               pad_hist={str(x): n for x, n in sorted(Counter(npads).items())}, exact_rate=rate, kernel_numerator=A,
               dft_numerator=A - 1, schedK_mirror=ok, schedK_duplicate_owner_control=schedK(nb, R, bad, npads, res),
               schedK_wrong_pad_control=schedK(nb, R, nrow, badp, res), packing_attempt=attempt, packing_seed=seed)
    assert ok and not rec["schedK_duplicate_owner_control"] and not rec["schedK_wrong_pad_control"]
    json.dump(rec, open(os.path.join(out, "price-packed.json"), "w"), indent=2)
    print(json.dumps(dict(banks=nb, W=W, pads=rec["pad_hist"], root=rate["root_lower_decimal"][:17], kernel=A, dft=A - 1,
                          schedK=ok, attempt=attempt)))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 0)
