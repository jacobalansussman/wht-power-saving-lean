#!/usr/bin/env python3
"""an-build calibration (own code, stdlib). A replica-level frame model of the #401 nine-stage line as our Lean layout
states it (HDef.lean, u8xmhh overlay: gluing x_{k,k+1}: Phi_k(1) = Phi_{k+1}(t), y_{k,k+1}: Phi_k(t^perp) = Phi_{k+1}(0);
host chains hH0_*: side X starts at Phi_0(t) / Phi_k(t), side Y at Phi_0(0) / Phi_k(0); idle children table CH of
x-hosting/tools/host_price.py).  Frames are tracked as DIMENSIONS of the global flag of one family (a port path); the
local word enters only through its per-call profile.  Checks: every data array pays 72 of 73, D = 320 per replica, and
the exact-fill price of the U8 line, 1.6647327e-3 (an-ours knobs.out), with our fold (an-ours/anatomy_ours.price, the
Lean fold hypothesis, s = 40)."""
import sys, os, json
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../anatomy-our-lean-words')))
from anatomy_ours import price, rank, narrow, add, LOCAL

M, STAGES, V, T = 73, 9, 112, 8
PAIR = [1, 1, 1, 2, 2, 3, 3, 4, 4]
CH = [('A2', 16, 2, 1, 'A'), ('B2', 8, 2, 2, 'B'), ('A3', 32, 3, 3, 'A'), ('B3', 8, 3, 4, 'B'),
      ('A4', 48, 4, 5, 'A'), ('B4', 8, 4, 6, 'B'), ('F1', 48, 1, 2, 'F'), ('F2', 32, 2, 4, 'F'), ('F3', 16, 3, 6, 'F')]

def chain(h=9, fresh=8, shared_y=0):
    """global dims of the family flag: base(k) = dim Phi_k(0); x-chain Phi_k(t) = base+1 .. Phi_k(1) = base+h;
    y-chain Phi_k(0)+shared_y .. Phi_k(t^perp) = base+h-1.  Gluing: base(k+1) = dim Phi_k(t^perp) - shared_y."""
    base = [0]
    for k in range(STAGES - 1):
        base.append(base[-1] + h - 1 - shared_y)
    return base

def replica(h=9, shared_y=0, cst=64):
    """per data array: start dim, in-call levels, idle levels, end dim (X, Y) for each pair; the run's total paid."""
    base = chain(h, shared_y=shared_y)
    m = base[-1] + h
    rows, paid = [], 0
    for p in (1, 2, 3, 4):
        act = [k for k in range(STAGES) if PAIR[k] == p]
        for side in ('X', 'Y'):
            start = 1 if side == 'X' else shared_y            # X on its line; Y at 0 (one-sided) or on a line
            end = m - shared_y if side == 'X' else m - 1        # X full (or s^perp); Y t^perp
            incall = sum((h - 1 - shared_y) for _ in act)      # x: t -> full (h-1) ; y: 0 -> t^perp (h-1); two-sided: h-2
            idle = end - start - incall
            # the idle children of the CH table must cover exactly the non-active windows
            ch = sum(r for _, r, q, _, _ in CH if q == p) if (h, shared_y) == (9, 0) else None
            rows.append(dict(pair=p, side=side, start=start, end=end, incall=incall, idle=idle, CH=ch))
            paid += 2 * V // 2 * 0 + V * (end - start)          # V families per side and pair
    paid += STAGES * cst
    W = 2 * 4 * V
    return dict(m=m, base=base, rows=rows, paid=paid, W=W, D=m * W - paid)

if __name__ == '__main__':
    out = {}
    r = replica()
    print('one-sided #401 chain: m =', r['m'], ' window bases', r['base'])
    for row in r['rows']:
        print('  pair %(pair)d %(side)s: start dim %(start)d -> end dim %(end)d : in-call %(incall)d + idle %(idle)d (CH table %(CH)s)' % row)
        assert row['incall'] + row['idle'] == 72 and (row['CH'] is None or row['CH'] == row['idle'])
    print('  paid per replica', r['paid'], ' = 73*896 - D, D =', r['D'])
    assert r['m'] == 73 and r['D'] == 320
    out['one_sided'] = dict(m=r['m'], D_per_replica=r['D'], paid=r['paid'])
    # exact-fill price of the U8 line: 72 x local profile + idle remainder in rank-48 pieces (an-ours knobs.fill)
    U = {k: 72 * n for k, n in LOCAL['U8E1'].items()}
    W = T * r['W']; D = T * r['D']
    rem = M * W - D - rank(U); q, rr = divmod(rem, 48)
    H = add(U, {48: q, rr: 1} if rr else {48: q})
    a = price(H, W)
    print('U8 exact fill: W', W, 'D', M * W - rank(H), 'price %.7e' % a)
    assert M * W - rank(H) == 2560 and abs(float(a) - 1.6647327e-3) < 1e-10
    out['U8_exact_fill'] = dict(W=W, D=M * W - rank(H), a=float(a))
    # the Lean-checked word as a cross-check of the fold (an-ours knobs.py profile of U8X0HQ)
    HQ = {1: 185776, 2: 25200, 3: 32368, 4: 12560, 5: 5184, 8: 576, 12: 96, 48: 2688, 59: 1}
    aq = price(HQ, 7491)
    print('U8X0HQ profile, W 7491: D', M * 7491 - rank(HQ), 'price %.7e (Lean: 1.6256629e-3)' % aq)
    out['U8X0HQ'] = dict(D=M * 7491 - rank(HQ), a=float(aq))
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'calib.json'), 'w'), indent=1)
