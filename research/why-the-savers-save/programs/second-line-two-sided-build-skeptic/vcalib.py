#!/usr/bin/env python3
"""v-an-build calibration, own code (stdlib).  Rebuilds:
 (1) per-call frame bookkeeping of the U8X word from its gcert/1 file: dims of start/final frames of x, y; centre frames;
 (2) the replica deficit D = W_rep * m - paid with my own count (W_rep = 4 pairs x 2 sides x 112 = 896, each data array
     m - 1 = 72 paid, 9 calls x 64 centre charge, helpers neutral) and the per-word D = 8 * D_rep;
 (3) the Lean fold root, my own Decimal bisection of
        (2^40 - 1) * sum_r n_r (r/m)^(1-a) + W ((m-1)/m)^(1-a) + W (1/m)^(1-a) < 2^40 W
     (Work/BlockAccounting/Moment.lean T, moment; the hypothesis printed in RateU8X0HQx.lean), on
       - the U8X0HQ block list copied from fsp-build/hq/lean/Work/XHost/Data/RateU8X0HQx.lean (Lean K = 16256630),
       - the U8 exact fill of an-ours (72 x local U8E1 profile from fsp-price/data/words.json + rank-48 idle pieces),
     and the Lean grid check: the inequality holds at a = K/1e10 and fails at (K+1)/1e10.
usage: python3 -I vcalib.py"""
import gzip, json, os
from decimal import Decimal as D, getcontext
getcontext().prec = 70
HOME = os.path.expanduser('~')
BC = __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs') + '/'


def lhs_minus_rhs(H, W, m, a, s=40):
    e = 1 - a
    lm = D(m)
    tot = sum(D(n) * ((D(r) / lm).ln() * e).exp() for r, n in H.items())
    fill = ((D(m - 1) / lm).ln() * e).exp() + ((D(1) / lm).ln() * e).exp()
    two = D(2) ** s
    return (two - 1) * tot + D(W) * fill - two * D(W)


def root(H, W, m):
    lo, hi = D(0), D('0.05')
    assert lhs_minus_rhs(H, W, m, lo) < 0 and lhs_minus_rhs(H, W, m, hi) > 0
    for _ in range(120):
        mid = (lo + hi) / 2
        if lhs_minus_rhs(H, W, m, mid) < 0: lo = mid
        else: hi = mid
    return lo


def dim(frames, i):
    return len(frames[i])


def main():
    out = {}
    d = json.load(gzip.open(BC + 'U8X-parsed.json.gz', 'rt'))
    v, fr = d['v'], d['frames']
    xs = sorted(set((dim(fr, d['start'][r]), dim(fr, d['final'][r])) for r in range(v)))
    ys = sorted(set((dim(fr, d['start'][r]), dim(fr, d['final'][r])) for r in range(v, 2 * v)))
    cen = [dim(fr, f) for _, _, f in d['ret']]
    print('U8X word: x (start dim, final dim):', xs, '| y:', ys, '| centre frames dims:', cen, 'sum', sum(cen))
    out['x_start_final_dims'] = xs; out['y_start_final_dims'] = ys; out['centre_dims'] = cen
    # local profile of the call from words.json
    wj = json.load(open(BC + 'words.json'))
    U8E1 = {int(k): n for k, n in wj['869/870/890_U8_exit_frames(E1)']['profile'].items()}
    wid = {int(k): n for k, n in wj['869/870/890_U8_exit_frames(E1)']['widths'].items()}
    lr = sum(r * n for r, n in U8E1.items()); lw = sum(w * n for w, n in wid.items())
    print('local rank per call', lr, '; leftover widths', lw, '; difference', lr - lw, '(= 224*8 + 64 =', 224 * 8 + 64, ')')
    out['local_rank_minus_widths'] = lr - lw
    # replica bookkeeping (my count)
    m, Wrep, calls, cst = 73, 4 * 2 * 112, 9, 64
    paid = Wrep * (m - 1) + calls * cst
    Drep = Wrep * m - paid
    print('replica: W', Wrep, 'paid', paid, 'D', Drep, '; per word (8 replicas): W', 8 * Wrep, 'D', 8 * Drep)
    out['D_replica'] = Drep
    # two-sided bookkeeping at m = 74: each data array pays m - 2 = 72, plus crossings
    for xgb in (0, 112, 224):
        D2 = Wrep * 74 - (Wrep * 72 + calls * cst + 4 * xgb)
        print('  two-sided m=74, %3d crossings per GB call (4 GB calls): D per replica %d' % (xgb, D2))
    # fold on U8X0HQ (Lean block list)
    HQ = {1: 185776, 2: 25200, 3: 32368, 4: 12560, 5: 5184, 8: 576, 12: 96, 48: 2688, 59: 1}
    W = 7491
    assert sum(r * n for r, n in HQ.items()) == W * 73 - 2560
    a = root(HQ, W, 73)
    K = 16256630
    ok_K = lhs_minus_rhs(HQ, W, 73, D(K) / D(10) ** 10) < 0
    ok_K1 = lhs_minus_rhs(HQ, W, 73, D(K + 1) / D(10) ** 10) < 0
    print('U8X0HQ: root a* = %s ; floor(a* 1e10) = %d ; inequality at K=%d: %s, at K+1: %s ; DFT figure (K-1)/1e10 = %.7e'
          % (str(a)[:14], int(a * D(10) ** 10), K, ok_K, ok_K1, (K - 1) / 1e10))
    out['U8X0HQ_root'] = float(a); out['U8X0HQ_K_holds'] = ok_K; out['U8X0HQ_K1_holds'] = ok_K1
    # U8 exact fill
    Wf = 8 * Wrep
    U = {r: 72 * n for r, n in U8E1.items()}
    rem = 73 * Wf - 8 * Drep - sum(r * n for r, n in U.items())
    q, rr = divmod(rem, 48)
    H = dict(U); H[48] = H.get(48, 0) + q
    if rr: H[rr] = H.get(rr, 0) + 1
    af = root(H, Wf, 73)
    print('U8 exact fill: W', Wf, 'rank-48 pieces', q, 'remainder', rr, 'D', 73 * Wf - sum(r * n for r, n in H.items()),
          'root %.7e' % float(af))
    out['U8_exact_fill_root'] = float(af)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'vcalib.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
