#!/usr/bin/env python3
"""an-ours (15 Oct 2026): why our Lean-checked words save, and what caps them.

Own code; reads only block lists copied by hand from our own Rate*.lean files (listed in WORDS below, with their source
paths) and the per-call local profiles of fsp-price/data/words.json (our own data file).  No third-party code is run.

Price: the Lean fold hypothesis of RAM.engine_program_group_list (fs-build/tools/fsb_rate_check.py states it):
    (2^s - 1) * sum_r n_r (r/m)^(1-a) + W ((m-1)/m)^(1-a) + W (1/m)^(1-a) < 2^s W,   s = 40, m = 73
  The kernel numerator is the largest K with the inequality at a = K/10^10; the DFT figure is K - 1.
  Here: Decimal bisection (prec 60) of the same inequality -> the real root a*; K = floor(a* 10^10) (the certificates sit
  on that grid point when term_le4 is used; with term_le the rate lemma loses ~30 grid points).
Equivalent moment form: sum_r (r n_r / (m W)) (m/r)^a = 1 - O(2^-40).
First-order reading: with sum_r r n_r = m W - D,  a* ~= D / N,  N = sum_r r n_r ln(m/r)  ("narrow mass").
"""
import json, math, sys
from decimal import Decimal as Dm, getcontext
getcontext().prec = 60
M = 73
S = 40

def price(H, W, m=M, s=S):
    """real root a* of the fold inequality (Decimal bisection)."""
    H = {int(r): int(n) for r, n in H.items() if int(n)}
    lm = Dm(m)
    two = Dm(2) ** s
    terms = [(Dm(n), (Dm(r) / lm).ln()) for r, n in H.items()]
    fill = [((Dm(m - 1) / lm).ln()), ((Dm(1) / lm).ln())]
    Wd = Dm(W)
    def lhs(a):
        e = 1 - a
        tot = sum(n * (e * l).exp() for n, l in terms)
        return (two - 1) * tot + Wd * sum((e * l).exp() for l in fill)
    lo, hi = Dm(0), Dm('0.01')
    for _ in range(200):
        mid = (lo + hi) / 2
        if lhs(mid) < two * Wd: lo = mid
        else: hi = mid
    return lo

def rank(H): return sum(int(r) * int(n) for r, n in H.items())
def narrow(H, m=M): return sum(int(r) * int(n) * math.log(m / int(r)) for r, n in H.items())
def add(*Hs, k=None):
    out = {}
    for i, H in enumerate(Hs):
        c = 1 if k is None else k[i]
        for r, n in H.items(): out[int(r)] = out.get(int(r), 0) + c * int(n)
    return {r: n for r, n in sorted(out.items()) if n}
def sub(A, B):
    out = {int(r): int(n) for r, n in A.items()}
    for r, n in B.items(): out[int(r)] = out.get(int(r), 0) - int(n)
    return {r: n for r, n in sorted(out.items()) if n}

REL = 'community/'
BC = 'bound-check/'
WORDS = [  # name, W, Lean kernel numerator, block list (copied from the Rate*.lean file named), source
 ('Shea C588 (#588 atm3, 112 source copies, 5040 banks)', 12208, 12331686,
  {1:189360, 2:35424, 3:20160, 4:16704, 5:4032, 6:4032, 8:14088, 16:3584, 32:3584, 48:3584},
  'shea-src research/complex-transfer-588-596/ATTEMPT.md (H = 72 H_local + 112 H_idle + H_padding)'),
 ('C588 repack (5009 banks)', 12177, 12383292,
  {1:187636, 2:35424, 3:20160, 4:16704, 5:4032, 6:4032, 8:14016, 16:3584, 32:3584, 37:1, 48:3584}, REL+'c588-repack/.../RateC588x.lean'),
 ('C588M (mutable sources)', 11200, 13552607,
  {1:178848, 2:35424, 3:20160, 4:16704, 5:4032, 6:4032, 8:5952, 16:3584, 28:72, 32:3584, 48:3584}, REL+'c588m/.../RateC588Mx.lean'),
 ('C588MHH (mutable + same-run hosting)', 7279, 15699306,
  {1:180640, 2:35424, 3:20224, 4:17056, 5:5728, 6:4032, 7:320, 8:1584, 39:1, 48:1680}, REL+'c588mhh/.../RateC588MHHx.lean'),
 ('U8X (U8 word, exit frames, source copies, banks)', 11889, 12704024,
  {1:191520, 2:24336, 3:32112, 4:11664, 5:5184, 8:14016, 16:3584, 32:3584, 41:1, 48:3584}, REL+'u8x/.../RateU8Xx.lean'),
 ('U8XM (+ mutable sources, 5 centre copies)', 10939, 13868857,
  {1:183816, 2:24336, 3:32112, 4:11664, 5:5184, 8:6312, 16:3584, 27:1, 32:3584, 48:3584}, REL+'u8xm/.../RateU8XMx.lean'),
 ('U8XMHH (+ same-run hosting)', 7478, 16197706,
  {1:186137, 2:25216, 3:32368, 4:12560, 5:5184, 8:936, 12:64, 48:2608, 61:1}, REL+'u8xmhh/.../RateU8XMHHx.lean'),
 ('U8X0 (no centre copies, banks)', 10895, 13928619,
  {1:183456, 2:24336, 3:32112, 4:11664, 5:5184, 8:5952, 16:3584, 32:3584, 48:3584, 55:1}, BC+'fsp-cc/out/global-U8X0/.../RateU8X0x.lean'),
 ('U8X0HH (no copies + hosting)', 7491, 16256507,
  {1:185776, 2:25201, 3:32368, 4:12560, 5:5184, 8:576, 12:96, 48:2688, 57:1}, BC+'fsp-cc/stage-u8x0hh/.../RateU8X0HHx.lean'),
 ('U8X0HQ (one-pad table, term_le4) = our best', 7491, 16256630,
  {1:185776, 2:25200, 3:32368, 4:12560, 5:5184, 8:576, 12:96, 48:2688, 59:1}, BC+'fsp-build/hq/lean/.../RateU8X0HQx.lean'),
 ('T74E (hub T74 gauge44 word, our lowering, banks)', 10855, 14040061,
  {1:168480, 2:36576, 3:30240, 4:9072, 5:5904, 8:5952, 15:1, 16:3584, 32:3584, 48:3584}, BC+'lg-a1/out/global-T74E/.../RateT74Ex.lean'),
]

LOCAL = {  # per-call local profiles (fsp-price/data/words.json; Shea ATTEMPT.md for the copy version)
 '588atm3': {1:2484, 2:492, 3:280, 4:232, 5:56, 6:56, 8:8},
 '588atm3+copies': {1:2596, 2:492, 3:280, 4:232, 5:56, 6:56, 8:120},
 'U8E1': {1:2548, 2:338, 3:446, 4:162, 5:72, 8:8},
 'U8XM_5cc': {1:2553, 2:338, 3:446, 4:162, 5:72, 8:13},
}
WIDTHS = {'588atm3': {4:48, 5:8, 7:32, 8:140, 9:276}, 'U8E1': {4:80, 5:32, 7:54, 8:122, 9:216}}

def band(H):
    b = {'r=1': 0.0, 'r=2..8': 0.0, 'r=9..32': 0.0, 'r>=33': 0.0}
    for r, n in H.items():
        v = r * n * math.log(M / r)
        k = 'r=1' if r == 1 else 'r=2..8' if r <= 8 else 'r=9..32' if r <= 32 else 'r>=33'
        b[k] += v
    return b

def main():
    out = {'words': []}
    print(f"{'word':52s} {'W':>6s} {'banks':>5s} {'D':>5s} {'N':>10s} {'N/Wd':>7s} {'D/N':>10s} {'a* (ours)':>13s} {'Lean K':>9s}")
    for name, W, K, H, src in WORDS:
        D = M * W - rank(H)
        N = narrow(H)
        a = price(H, W)
        Kours = int(a * 10**10)
        Wd = 7168 if not name.startswith('T74') else 7168
        row = dict(name=name, W=W, banks=W - 7168, D=D, rank=rank(H), N=round(N, 1), N_per_data_array=round(N / 7168, 3),
                   first_order=D / N, a_star=str(a)[:16], K_ours=Kours, K_lean=K, bands={k: round(v, 1) for k, v in band(H).items()},
                   src=src)
        out['words'].append(row)
        print(f"{name:52s} {W:6d} {W-7168:5d} {D:5d} {N:10.1f} {N/7168:7.2f} {D/N:10.4e} {float(a):13.7e} {K:9d}  ours K={Kours}")
    # decomposition of each word: 72 x local + rest
    print('\nDecomposition: profile = 72 x local word + rest (idle remainders, bank/host pieces, pads, copies)')
    dec = {}
    for name, W, K, H, src in WORDS:
        if name.startswith('Shea') or name.startswith('C588 repack'): loc = LOCAL['588atm3+copies']
        elif name.startswith('C588'): loc = LOCAL['588atm3']
        elif name.startswith('U8XM') : loc = LOCAL['U8XM_5cc']
        elif name.startswith('U8X (') : loc = None
        elif name.startswith('U8'): loc = LOCAL['U8E1']
        else: loc = None
        if loc is None: continue
        L72 = add(loc, k=[72]) if False else {r: 72 * n for r, n in loc.items()}
        rest = sub(H, L72)
        neg = {r: n for r, n in rest.items() if n < 0}
        NL, NR = narrow(L72), narrow({r: n for r, n in rest.items() if n > 0})
        dec[name] = dict(local72_rank=rank(L72), local72_N=round(NL, 1), rest=rest, rest_rank=rank(rest), rest_N=round(NR, 1),
                         share_N_local=round(NL / (NL + NR), 4), negative=neg)
        print(f"{name:52s} local: rank {rank(L72):7d} N {NL:10.1f} | rest: rank {rank(rest):7d} N {NR:9.1f} ({100*NR/(NL+NR):5.2f} % of N)  rest={rest}")
    out['decomposition'] = dec
    # ceilings for the U8 line
    print('\nCeilings and what-ifs (price a* of the same fold; our own arithmetic)')
    U = {r: 72 * n for r, n in LOCAL['U8E1'].items()}
    cases = {}
    # (i) exact fill: W = 7168 data arrays, no banks, every residual hosted, idle remainder all rank 48 (+ one leftover piece)
    rem = 73 * 7168 - 2560 - rank(U)
    q, rr = divmod(rem, 48)
    H1 = add(U, {48: q} if not rr else {48: q, rr: 1})
    cases['U8 at an exact fill (W 7168, no banks, idle remainder {48: %d%s})' % (q, ', %d: 1' % rr if rr else '')] = (H1, 7168)
    # (ii) same, idle remainder as rank 72 pieces (the largest legal move), to see how much the idle split matters
    q2, r2 = divmod(rem, 72)
    H2 = add(U, {72: q2} if not r2 else {72: q2, r2: 1})
    cases['U8, idle remainder in rank-72 pieces (not legal in this network; scale)'] = (H2, 7168)
    # (iii) U8X0HQ with its 323 banks deleted but every bank piece kept (W 7168): is it the banks' W or their pieces?
    HQ = WORDS[9][3]
    cases['U8X0HQ profile, W forced to 7168 (bank pieces kept): infeasible rank, shown for scale'] = (HQ, 7168)
    # (iv) the local word alone: if the 72 x local climbs were the only narrow mass (rest at rank 72)
    for k, (H, W) in cases.items():
        D = 73 * W - rank(H)
        a = price(H, W)
        cases[k] = dict(W=W, D=D, N=round(narrow(H), 1), D_over_N=D / narrow(H), a_star=float(a))
        print(f"  {k}: W {W} D {D} N {narrow(H):.1f} D/N {D/narrow(H):.4e} a* {float(a):.7e}")
    out['ceilings'] = cases
    # knobs at fixed profile shape: deficit per data array is 5/14 on every word of the line (and on T74: 23680/66304)
    print('\nDeficit per data array: U8/C588 2560/7168 =', 2560 / 7168, '; T74 (hub, T = 37) 23680/66304 =', 23680 / 66304)
    print('Network geometry (pr717 P4): m = 1 + 8(2P+1), W = 2PTv, D = T[2Pv - (2P+1)cst]; P=4, T=8, v=112, cst=64 ->',
          1 + 8 * 9, 2 * 4 * 8 * 112, 8 * (8 * 112 - 9 * 64))
    # sensitivities on U8X0HQ: one more deficit unit; one fewer bank (its pieces hosted); narrow mass -1 %
    base = price(HQ, 7491)
    sens = {}
    sens['D + 72 (one unit per call; e.g. a centre climb 8 -> 7)'] = price(sub(HQ, {1: 72}), 7491)  # 72 rank-1 moves removed
    # one bank fewer: remove its 73 rank units of pad/pieces: model as 73 rank-1 -> impossible; use the bank's real mix
    sens['one port fewer (v 111): W -64, D -64, local scaled 111/112'] = None
    for k, v in sens.items():
        if v is not None:
            print(f"  sensitivity {k}: a* {float(v):.7e} ({100*(float(v)/float(base)-1):+.3f} %)")
    out['sensitivity'] = {k: (float(v) if v is not None else None) for k, v in sens.items()}
    json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else 'anatomy_ours.json', 'w'), indent=1, default=str)

if __name__ == '__main__':
    main()
