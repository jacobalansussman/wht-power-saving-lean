#!/usr/bin/env python3
"""an-ours knobs (own code): what-if prices on the U8 line with the same fold as anatomy_ours.price."""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from anatomy_ours import price, rank, narrow, add, sub, LOCAL
M = 73
U = {r: 72 * n for r, n in LOCAL['U8E1'].items()}
def fill(local, W, D=2560, piece=48):
    rem = M * W - D - rank(local)
    q, rr = divmod(rem, piece)
    H = add(local, {piece: q}) if not rr else add(local, {piece: q, rr: 1})
    return H
res = {}
def show(k, H, W):
    D = M * W - rank(H); a = float(price(H, W)); N = narrow(H)
    res[k] = dict(W=W, D=D, N=round(N, 1), D_over_N=D / N, a=a, mean_log=N / rank(H))
    print(f"{k:78s} W {W:5d} D {D:5d} N {N:10.1f} lbar {N/rank(H):.3f} a* {a:.7e}")
HQ = {1:185776, 2:25200, 3:32368, 4:12560, 5:5184, 8:576, 12:96, 48:2688, 59:1}
show('U8X0HQ (our best, Lean-checked)', HQ, 7491)
show('U8 at exact fill (ceiling of the U8 line in this network)', fill(U, 7168), 7168)
# centres free: remove the 8 rank-8 centre climbs per call (cst 64 -> 0): D = W
Uc = sub(U, {8: 576})
show('exact fill, central totals free (cst 0, D = W: one saved direction per array)', fill(Uc, 7168, D=7168), 7168)
# U8X0HQ with centres free (same rest)
show('U8X0HQ rest, central totals free', sub(HQ, {8: 576}), 7491)
# local word with every rank-1 climb merged into rank-2 climbs (hypothetical: halves the count, same rank)
U2 = dict(U); n1 = U2.pop(1); U2[2] = U2.get(2, 0) + n1 // 2
if n1 % 2: U2[1] = 1
show('exact fill, every local rank-1 climb paired into rank-2 (hypothetical)', fill(U2, 7168), 7168)
# small-unit special case: all moves rank 1, one saved move per array at block m
for m in (3, 4, 5, 6, 7, 9, 16, 73):
    print(f"  one saved move per array, all moves rank 1, block m={m:3d}: s = ln(m/(m-1))/ln m = {math.log(m/(m-1))/math.log(m):.4f}")
    res[f'one_saved_move_rank1_m{m}'] = math.log(m/(m-1))/math.log(m)
# per-data-array deficit across results (README section 6 table and the Rate files)
deltas = {'U8/C588 line (hub P112 network, P=4, T=8)': (2560, 7168), 'hub T74 (T=37)': (23680, 66304),
          'our 4th result (E8 unit, bridged 5-stage, 4v data)': (120, 480), 'our 2nd/3rd (paired cube p=11)': (3080, 4*1320),
          'our 1st (NStar3)': (1035, 4*560)}
for k, (D, Wd) in deltas.items():
    print(f"  deficit per data array, {k}: {D}/{Wd} = {D/Wd:.4f}")
    res['delta ' + k] = D / Wd
json.dump(res, open('knobs.json', 'w'), indent=1)
