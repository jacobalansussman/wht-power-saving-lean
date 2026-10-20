#!/usr/bin/env python3
"""m-gen (13 Oct 2026): the bank schedules and price files of the copy-free #588 certificate `C588M`.

Two schedules over the 504 helpers (residual widths res[q] = 9 - dim(entrance frame of slot q), in HELPER ORDER,
exactly the list `schedRes` that Lean's `res_ok` compares with the start state of the label check):

* annealed (the build target): x-mutable's annealed packing (d-lean-route's sa_pack, K = 56 bins, seed 0:
  x-mutable/out/sa/588atm3_mutable_sources_no_copies.schedule.json.gz) whose bins are given in sorted-width
  positions; the bins are mapped onto the actual helper ids width class by width class (as d-lean-route's
  make_price_json.py does for Shea's control), then replicated cyclically (bin b, call t, member j -> bank
  b*72 + (t+j) mod 72), validated;
* cyclic (fallback): Ryan Shea's `cyclic_schedule` (bank_reprice.py) on the same widths.

For each: the global histogram (72 invocations of the kernel-counted local histogram, Bridge377 idle banks
{48:32, 32:32, 16:32, 8:48} per source, the pads), W = 64 v + banks, the deficit 73 W - mass = 2560, the exact
rational moment-root bracket (Shea's exact_moment.certify), the largest kernel numerator A <= 10^10 root with a
finite-fill rate certificate (Sussman's foldlib.plan_block, B = 10, s = 40), DFT numerator A - 1; an exact mirror
of Lean's GX.schedK (Work/Bridge401/SchedK.lean) on the literal tables, and the two negative controls (duplicate
owner in one call; a wrong pad) evaluated by that mirror.
Output: <out>/price-<variant>.json (bank_reprice format, read by generate_global_lean_z.py) and
<out>/price-<variant>.schedule.json.gz ({rows, pads, widths}).
usage: make_sched.py <cert.gcert.json.gz> <annealed.schedule.json.gz> <toolkit dir> <out dir>"""
import sys, os, json, gzip, hashlib
from collections import Counter, defaultdict
from fractions import Fraction as Q

LOCAL_IDLE = {48: 32, 32: 32, 16: 32, 8: 48}


def schedK(nb, R, rows, pad, res):
    """Lean's GX.schedK, literally"""
    def skInj(row):
        seen = 0
        for x in row:
            if not (x < nb) or (seen >> x) & 1:
                return False
            seen |= 1 << x
        return True
    if len(rows) != 72 or not all(len(row) == R and skInj(row) for row in rows):
        return False
    if len(pad) != nb or not all(p < 73 for p in pad) or len(res) != R:
        return False
    if not (72 * sum(res) < 2 ** 20):
        return False
    lhs = sum(sum(r << (20 * x) for x, r in zip(row, res)) for row in rows)
    rhs = sum((73 - p) << (20 * i) for i, p in enumerate(pad))
    return lhs == rhs


def cyclic_rows(widths, groups, count=72, m=73):
    rows = [[None] * len(widths) for _ in range(count)]
    pads = []
    for b, group in enumerate(groups):
        pads.extend([m - sum(widths[q] for q in group)] * count)
        for t in range(count):
            for j, q in enumerate(group):
                rows[t][q] = b * count + (t + j) % count
    return rows, pads


def main(cert, sa_sched, toolkit, out):
    sys.path.insert(0, os.path.join(toolkit, "experiments"))
    sys.path.insert(0, os.path.join(toolkit, "verification"))
    sys.path.insert(0, os.path.join(toolkit, "tools", "gen"))
    from bank_reprice import cyclic_schedule, validate
    from exact_moment import certify
    import foldlib as fl
    c = json.loads(gzip.open(cert, "rt").read())
    v, R, h = c["v"], c["R"], c["h"]
    dim = [len(f) for f in c["frames"]]
    res = [h - dim[c["start"][2 * v + q]] for q in range(R)]
    assert all(dim[c["final"][2 * v + q]] == h for q in range(R))
    local = Counter()
    for w in c["blocks"].values():
        for r, n in w.items():
            local[int(r)] += n
    local = dict(sorted(local.items()))
    s = json.loads(gzip.decompress(open(sa_sched, "rb").read()))
    w = s["widths"]
    assert Counter(w) == Counter(res), (Counter(w), Counter(res))
    validate(s["rows"], s["pads"], w)
    nb0 = len(s["pads"]) // 72
    groups = defaultdict(list)
    for q, b in enumerate(s["rows"][0]):
        groups[b // 72].append(q)
    pools = defaultdict(list)
    for q, x in enumerate(res):
        pools[x].append(q)
    G = [[pools[w[q]].pop() for q in sorted(groups[b])] for b in range(nb0)]
    assert all(not p for p in pools.values())
    variants = {"annealed": cyclic_rows(res, G), "cyclic": cyclic_schedule(res)}
    summary = {}
    for name, (rows, pads) in variants.items():
        validate(rows, pads, res)
        nb = len(pads)
        ok = schedK(nb, R, rows, pads, res)
        bad = json.loads(json.dumps(rows)); bad[0][1] = bad[0][0]          # Shea's duplicate-owner control
        badp = list(pads); badp[0] += 1                                      # wrong padding
        hist = Counter({r: 72 * n for r, n in local.items()})
        for r, n in LOCAL_IDLE.items():
            hist[r] += n * v
        hist.update(p for p in pads if p)
        W = 64 * v + nb
        mass = sum(r * n for r, n in hist.items())
        assert 73 * W - mass == 2560 == 64 * v - 72 * c["cst"], (73 * W - mass)
        rate = certify(dict(hist), 73, W)
        u = dict(m=73, W=W, D=2560, H=dict(hist))
        A = int(Q(rate["root_lower"]) * 10 ** 10)
        while fl.plan_block(u, A, 10, 40) is None:
            A -= 1
        sp = os.path.join(out, "price-%s.schedule.json.gz" % name)
        open(sp, "wb").write(gzip.compress(json.dumps(dict(rows=rows, pads=pads, widths=res), separators=(",", ":")).encode(),
                                           mtime=0))
        rec = dict(status="EXACT_ARITHMETIC_CANDIDATE_MUTABLE_SOURCES_NO_COPIES", source_copies=False, schedule=name,
                   source_sha256=hashlib.sha256(open(cert, "rb").read()).hexdigest(),
                   schedule_sha256=hashlib.sha256(open(sp, "rb").read()).hexdigest(),
                   h=h, m=73, v=v, local_profile={str(r): n for r, n in local.items()}, local_residual_count=R,
                   banks=nb, W=W, deficit=2560, profile={str(r): n for r, n in sorted(hist.items())},
                   pad_hist={str(p): n for p, n in sorted(Counter(pads).items())}, exact_rate=rate,
                   kernel_numerator=A, dft_numerator=A - 1,
                   schedK_mirror=ok, schedK_duplicate_owner_control=schedK(nb, R, bad, pads, res),
                   schedK_wrong_pad_control=schedK(nb, R, rows, badp, res))
        assert ok and not rec["schedK_duplicate_owner_control"] and not rec["schedK_wrong_pad_control"]
        json.dump(rec, open(os.path.join(out, "price-%s.json" % name), "w"), indent=2)
        summary[name] = dict(banks=nb, W=W, pads=rec["pad_hist"], root=rate["root_lower_decimal"][:17],
                             kernel=A, dft=A - 1, schedK=ok)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:5])
