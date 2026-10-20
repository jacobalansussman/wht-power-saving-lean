#!/usr/bin/env python3
"""k-lower (13 Oct 2026): the bank schedule and price file of an EXIT-FRAME certificate.

Residual of slot q (in slot order, the list `schedRes` that Lean's `res_ok` compares with the label states):
  dim(final frame) - dim(start frame)   for an exit slot (x-exit-lean's `fitsX`, `GLabX.res`),
  h - dim(start frame)                  otherwise.
The bank table is x-exit-words' pad-minimal table (pack_ceil.py: 4721 banks, one pad of 41) given in that tool's slot
order; it is mapped onto this certificate's slots width class by width class (loads unchanged), then checked with
the exact mirror of Lean's GX.schedK (m-gen's make_sched.schedK) and Shea's validate, with the two negative controls.
Then (as m-gen's make_sched.py / Shea's bank_reprice.py): global histogram (72 calls, Bridge377 idle banks, pads),
W = 64 v + banks, deficit 2560, exact moment root (Shea's exact_moment.certify), kernel numerator by Sussman's
foldlib.plan_block (B = 10, s = 40), DFT numerator = kernel - 1.
usage: make_sched_x.py <cert.gcert.json.gz> <xw.schedule.json.gz> <toolkit dir> <out prefix>"""
import sys, os, json, gzip, hashlib
from collections import Counter, defaultdict
from fractions import Fraction as Q
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_sched import schedK, LOCAL_IDLE


def main(cert, xsched, toolkit, outp):
    sys.path.insert(0, os.path.join(toolkit, "experiments"))
    sys.path.insert(0, os.path.join(toolkit, "verification"))
    sys.path.insert(0, os.path.join(toolkit, "tools", "gen"))
    from bank_reprice import validate
    from exact_moment import certify
    import foldlib as fl
    c = json.loads(gzip.open(cert, "rt").read())
    v, R, h = c["v"], c["R"], c["h"]
    dim = [len(f) for f in c["frames"]]
    E = set(c["exits"])
    res = [(dim[c["final"][2 * v + q]] if 2 * v + q in E else h) - dim[c["start"][2 * v + q]] for q in range(R)]
    local = Counter()
    for w in c["blocks"].values():
        for r, n in w.items():
            local[int(r)] += n
    local = dict(sorted(local.items()))
    s = json.loads(gzip.decompress(open(xsched, "rb").read()))
    w = s["widths"]
    assert Counter(w) == Counter(res), (Counter(w), Counter(res))
    validate(s["rows"], s["pads"], w)
    mine, theirs = defaultdict(list), defaultdict(list)
    for q, x in enumerate(res): mine[x].append(q)
    for q, x in enumerate(w): theirs[x].append(q)
    perm = {}
    for x in mine:
        for a, b in zip(mine[x], theirs[x]): perm[a] = b
    rows = [[row[perm[q]] for q in range(R)] for row in s["rows"]]
    pads = s["pads"]
    validate(rows, pads, res)
    nb = len(pads)
    ok = schedK(nb, R, rows, pads, res)
    bad = json.loads(json.dumps(rows)); bad[0][1] = bad[0][0]
    badp = list(pads); badp[0] += 1
    hist = Counter({r: 72 * n for r, n in local.items()})
    for r, n in LOCAL_IDLE.items():
        hist[r] += n * v
    hist.update(p for p in pads if p)
    W = 64 * v + nb
    mass = sum(r * n for r, n in hist.items())
    assert 73 * W - mass == 2560, (73 * W - mass)
    rate = certify(dict(hist), 73, W)
    u = dict(m=73, W=W, D=2560, H=dict(hist))
    A = int(Q(rate["root_lower"]) * 10 ** 10)
    while fl.plan_block(u, A, 10, 40) is None:
        A -= 1
    sp = outp + ".schedule.json.gz"
    open(sp, "wb").write(gzip.compress(json.dumps(dict(rows=rows, pads=pads, widths=res), separators=(",", ":")).encode(),
                                       mtime=0))
    rec = dict(status="EXACT_ARITHMETIC_CANDIDATE_EXIT_FRAMES_TAIL_IN_B", schedule="x-exit-words pad-minimal (remapped)",
               source_sha256=hashlib.sha256(open(cert, "rb").read()).hexdigest(),
               schedule_sha256=hashlib.sha256(open(sp, "rb").read()).hexdigest(),
               h=h, m=73, v=v, local_profile={str(r): n for r, n in local.items()}, local_residual_count=R,
               banks=nb, W=W, deficit=2560, profile={str(r): n for r, n in sorted(hist.items())},
               pad_hist={str(p): n for p, n in sorted(Counter(pads).items())}, exact_rate=rate,
               kernel_numerator=A, dft_numerator=A - 1,
               schedK_mirror=ok, schedK_duplicate_owner_control=schedK(nb, R, bad, pads, res),
               schedK_wrong_pad_control=schedK(nb, R, rows, badp, res))
    assert ok and not rec["schedK_duplicate_owner_control"] and not rec["schedK_wrong_pad_control"]
    json.dump(rec, open(outp + ".json", "w"), indent=2)
    print(json.dumps(dict(banks=nb, W=W, pads=rec["pad_hist"], root=rate["root_lower_decimal"][:17], kernel=A, dft=A - 1,
                          schedK=ok), indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:5])
