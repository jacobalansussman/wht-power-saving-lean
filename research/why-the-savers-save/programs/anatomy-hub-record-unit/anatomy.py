"""an-hub anatomy: price hub complex profiles with our own few lines, decompose the saving, try levers.

Model (hub's own, read in research/complex-local-word-ceiling/ceiling.py and the bottleneck brief):
  bill(a) = sum_r n_r (r/m)^(1-a) / W ;  saving a = largest a with bill(a) <= 1.
  W = live data families (helpers hosted), m = 73.
Our language: "W in M at block m" with rank-1 moves gives s = ln(mW/M)/ln m; with children of rank r,
  a ~= -ln(bill(0)) / ln(m / r_eff), r_eff = exp(mass-weighted mean of ln r).
Standard library only; floats + bisection (relative accuracy ~1e-12 is enough for anatomy).
"""
import json, math, sys

def bill(prof, m, W, a):
    return sum(n * (r / m) ** (1 - a) for r, n in prof.items()) / W

def saving(prof, m, W):
    lo, hi = 0.0, 0.2
    if bill(prof, m, W, 0.0) >= 1:
        return float('nan')
    for _ in range(200):
        mid = (lo + hi) / 2
        if bill(prof, m, W, mid) <= 1: lo = mid
        else: hi = mid
    return lo

def decompose(prof, m, W):
    mass = sum(r * n for r, n in prof.items())
    b0 = mass / (m * W)
    lnr = sum(r * n * math.log(r) for r, n in prof.items()) / mass
    return dict(mass=mass, bill0=b0, deficit_frac=1 - b0, reff=math.exp(lnr),
                a_rank1=-math.log(b0) / math.log(m), a_first_order=-math.log(b0) / (math.log(m) - lnr))

def add(p, q, s=1):
    out = dict(p)
    for r, n in q.items():
        out[r] = out.get(r, 0) + s * n
        if out[r] == 0: del out[r]
    return out

def scale(p, k):
    return {r: n * k for r, n in p.items()}

def report(name, prof, m, W):
    d = decompose(prof, m, W)
    a = saving(prof, m, W)
    d.update(a=a, name=name)
    print(f"{name:58s} a={a:.10e}  mass={d['mass']:.6g}  1-bill0={d['deficit_frac']:.6e}  r_eff={d['reff']:.4f}"
          f"  a(all rank1)={d['a_rank1']:.6e}  a(1st order)={d['a_first_order']:.6e}")
    return d

if __name__ == '__main__':
    m = 73
    res = {}
    # T74 (#1839 line): global profile as shipped in the #1424-era PROFILE.json (rl-price/hub1839)
    t74 = {1: 1559100, 2: 339000, 3: 277392, 4: 84912, 5: 54612, 8: 5328, 48: 22748}
    Wt74, calls74 = 74 * 896, 74 * 9
    res['T74'] = report('T74 profile (PROFILE.json, PR1424 era)', t74, m, Wt74)
    mass = sum(r * n for r, n in t74.items())
    print('  coverage check: mass', mass, '== 72*W + 64*calls =', 72 * Wt74 + 64 * calls74)
    w = Wt74 / calls74
    print(f'  per call: w={w:.4f} families, paid mass/call={mass/calls74:.2f}, plain m*w={m*w:.2f}, '
          f'deficit/call={m*w-mass/calls74:.4f} = w - 64 = {w-64:.4f}')
    # anatomy in our "W in M" words, per call
    print(f'  our words: block {m}, W={w:.3f} arrays/call, M={mass/calls74:.1f} paid mass/call, '
          f'rank-1 saving ln(mW/M)/ln m = {math.log(m*w/(mass/calls74))/math.log(m):.6e}')
    idle = {48: 22748}
    local = add(t74, idle, -1)
    print('  T74 local mass/call', sum(r*n for r,n in local.items())/calls74, ' idle mass/call', 48*22748/calls74)
    # Levers on T74 (one premise at a time; favourable relaxations, realisation cost omitted)
    lv = {}
    # (i) idle at rank 51 (halving-guard cap) instead of 48
    im = 48 * 22748
    lv['idle at 51 (fractional)'] = saving(add(local, {51: im / 51}), m, Wt74)
    # (ii) all local mass of rank<=8 merged into rank-8 / rank-16 children (mass kept)
    lm = sum(r * n for r, n in local.items())
    lv['local mass as rank-8'] = saving({8: lm / 8, 48: 22748}, m, Wt74)
    lv['local mass as rank-16'] = saving({16: lm / 16, 48: 22748}, m, Wt74)
    lv['everything at rank 51'] = saving({51: mass / 51}, m, Wt74)
    # (iii) centre path free: remove 64 per call of rank-8 mass (8 children of rank 8 per call)
    lv['centre free (drop 8x rank-8 per call)'] = saving(add(t74, {8: 8 * calls74}, -1), m, Wt74)
    lv['centre halved (drop 4x rank-8 per call)'] = saving(add(t74, {8: 4 * calls74}, -1), m, Wt74)
    # (iv) two unpaid coordinates per family: remove W units of rank-1 mass
    lv['2 unpaid coords/family (drop W rank-1)'] = saving(add(t74, {1: Wt74}, -1), m, Wt74)
    # (v) more families per call at the same centre: scale everything except centre by k
    for k in (2, 4):
        p = add(scale(add(t74, {8: 8 * calls74}, -1), k), {8: 8 * calls74})
        lv[f'{k}x families per call, same centre'] = saving(p, m, Wt74 * k)
    for k, v in lv.items():
        print(f'  lever {k:45s} a={v:.6e}  x{v/res["T74"]["a"]:.3f}')
    res['T74_levers'] = lv
    # #1042 (T16, 144 calls): local profile (same as #1094's) + derived idle (hub ceiling README)
    loc1094 = {1: 2412, 2: 418, 3: 432, 4: 142, 5: 82, 8: 8}
    idle1042 = {1: 2249, 2: 1, 3: 400, 7: 3, 16: 47, 48: 4850}
    p1042 = add(scale(loc1094, 144), idle1042)
    res['1042'] = report('#1042 (T16): 144 x local + derived idle', p1042, m, 16 * 896)
    print('  coverage check #1042:', sum(r*n for r,n in p1042.items()), '==', 72*16*896 + 64*144)
    # #1094 packing ceiling (all idle at 51), the hub's figure 0.001690054
    idle_mass = 1646 * 144
    res['1094ceil'] = report('#1094 local + all idle mass at rank 51 (hub ceiling)', add(scale(loc1094, 144), {51: idle_mass / 51}), m, 16 * 896)
    # what the saving would be if the SAME deficit were spent with rank-1 children only (our plain-move game)
    print('\nShare of the T74 saving due to child grouping (r_eff>1):',
          1 - res['T74']['a_rank1'] / res['T74']['a'])
    json.dump(res, open(sys.argv[1], 'w'), indent=1, default=float) if len(sys.argv) > 1 else None
