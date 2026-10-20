#!/usr/bin/env python3
"""an2-smallnet / ceiling.py (own code, stdlib; python3 -I).  Ceilings for a small-block crossing network, solved
exactly from the moment equation  sum_children (r/m)^(1-a) = W  (the hub's per-child pricing), under deliberately
generous assumptions: helpers and junk cleaning free (no arrays, no levels), the centre charge at the hub's own
accounting t(k-1) per call (an assumption for D6; the true floor is unknown), a long chain (every pair in 2 calls),
every data array's levels in blocks of sigma (the window it takes up per call), centre levels in blocks of sigma.
u = free directions per data array: 1 (one-sided), 1.5 (the build's ceiling for a second line), 2 (two-sided with
the second line free, which the build refuted at our word).  Calibration: proposal 3's own assumptions (u = 2,
free bridges, r_eff = 4, m = 16) must give an-synth's 0.096."""
import json, math, os

def solve(children, W):
    """children: list of (rank, count); m implied by caller via ratio list (r/m, count)."""
    lo, hi = -5.0, 0.999999
    f = lambda a: sum(c * x ** (1 - a) for x, c in children) - W
    # f decreasing in a? x<1: x^(1-a) increases with a. so f increasing in a.
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > 0: hi = mid
        else: lo = mid
    return (lo + hi) / 2

def ceiling(m, sigma, cst_per_port, u, calls_per_pair=2):
    """per pair: 2 data arrays paying m - u levels in blocks of sigma; centre: calls_per_pair * cst_per_port levels in
    blocks of sigma; W = 2."""
    data_levels = 2 * (m - u)
    centre_levels = calls_per_pair * cst_per_port
    ch = [(sigma / m, data_levels / sigma), (sigma / m, centre_levels / sigma)]
    return solve(ch, 2.0), 2 * u - centre_levels

out = {}
# calibration: proposal 3 (an-synth 4.3): u = 2 per array, bridges free, r_eff = 4, m = 16
a_cal = solve([(4 / 16, 2 * (16 - 2) / 4)], 2.0)
print(f"calibration, proposal 3's own assumptions (m 16, u 2, free bridges, r_eff 4): a = {a_cal:.4f} (an-synth: 0.096)")
out['calibration_prop3'] = a_cal
rows = [  # design, sigma, hub-type centre levels per port per call, smallest block 1 + 3 sigma
    ('D6 (hub-type charge 0.8 assumed, not shown)', 4, 0.8, 13),
    ('E7', 6, 0.6667, 19),
    ('P112', 8, 0.5714, 25),
    ('E8 (120 ports)', 8, 0.5333, 25),
]
for name, sg, c, m in rows:
    for u in (1, 1.5, 2):
        a, Dpair = ceiling(m, sg, c, u)
        print(f"{name:45s} m={m:2d} u={u:3.1f}: D per pair {Dpair:5.2f} (D/W {Dpair/2:5.3f}), ceiling a = {a:.4f}")
        out[f"{name}|u={u}"] = dict(m=m, sigma=sg, cst_per_port=c, u=u, D_per_pair=Dpair, a=a)
# the same designs at the block proposal 3 wanted (m = 16) if the wall did not exist (windows shared, sigma kept)
print('-- if stages could share windows (the premise of proposal 3; excluded by the wall):')
for name, sg, c, m in rows:
    for u in (1, 1.5):
        mm = max(sg + 2, 10)
        a, Dpair = ceiling(mm, sg, c, u)
        print(f"{name:45s} m={mm:2d} u={u:3.1f}: ceiling a = {a:.4f}")
        out[f"shared|{name}|u={u}"] = dict(m=mm, a=a)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ceiling.json'), 'w'), indent=1)

# --- the designs that fit a block <= 16, at their PROVED centre floors (classfloor.py), all three u, block 1 + 3 sigma
#     (+1 for a second line's coordinate when u > 1); helpers free, ideal grouping: the most generous reading.
print('-- blocks <= 16, proved centre floors (classfloor.json), generous everywhere else:')
fl = {d['design']: d['floor_per_port'] for d in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'classfloor.json')))}
sig = {'D4': 2, 'A3': 2, 'D5': 4, 'D6': 4, 'A4': 4, 'A5': 4}
best = 0
for name in ('D4', 'A3', 'D5', 'D6', 'A4', 'A5'):
    for u in (1, 1.5, 2):
        m = 1 + 3 * sig[name] + (1 if u > 1 else 0)
        a, Dpair = ceiling(m, sig[name], fl[name], u)
        best = max(best, a)
        print(f"{name}: floor {fl[name]:.2f}/port/call, m={m:2d}, u={u:3.1f}: D per pair <= {Dpair:5.2f}, ceiling a <= {a:.4f}")
        out[f"small|{name}|u={u}"] = dict(m=m, floor=fl[name], u=u, D_per_pair=Dpair, a=a)
print(f"largest generous ceiling among the designs that fit m <= 16: {best:.4f}")
out['small_best'] = best
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ceiling.json'), 'w'), indent=1)

# --- the block floors used in NOTE.md: proved 1 + 2 sigma (an array that is a source in two amortised calls),
#     argued 1 + 3 sigma (its third call too); two-sided lines with a mirrored second window (argued): 2 + 3 sigma.
print('-- at the proved floor 1 + 2 sigma, the argued floor 1 + 3 sigma, and two-sided with mirrored windows 2 + 3 sigma:')
for name, sg, c in (('E7', 6, 0.6667), ('P112', 8, 0.5714), ('E8 (120 ports)', 8, 0.5333)):
    for label, m, us in (('1+2sigma (proved)', 1 + 2 * sg, (1, 1.5, 2)), ('1+3sigma (argued)', 1 + 3 * sg, (1, 1.5, 2)),
                         ('2+3sigma two-sided', 2 + 3 * sg, (2,))):
        for u in us:
            a, Dpair = ceiling(m, sg, c, u)
            print(f"{name:16s} {label:20s} m={m:2d} u={u:3.1f}: ceiling a = {a:.4f}")
            out[f"floors|{name}|{label}|u={u}"] = dict(m=m, u=u, a=a)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ceiling.json'), 'w'), indent=1)

# --- the small designs at the PROVED block floor 1 + 2 sigma (+1 for a second line), with the proved centre floors
#     and with the hub-type charges t(k-1)/n: this is where the rigorous reading leaves a gap (NOTE sec. 6).
print('-- small designs at the proved floor 1 + 2 sigma:')
hub = {'D4': 0.6667, 'A3': 1.0, 'D5': 1.0, 'D6': 0.8, 'A4': 1.6, 'A5': 1.3333}
for name in ('D4', 'A3', 'D5', 'D6', 'A4', 'A5'):
    for lab, c in (('proved centre floor', fl[name]), ('hub-type charge', hub[name])):
        for u in (1, 1.5, 2):
            m = 1 + 2 * sig[name] + (1 if u > 1 else 0)
            a, Dpair = ceiling(m, sig[name], c, u)
            print(f"{name}: {lab:20s} {c:.3f}/port/call, m={m:2d}, u={u:3.1f}: D per pair {Dpair:5.2f}, ceiling a = {a:.4f}")
            out[f"proved_floor|{name}|{lab}|u={u}"] = dict(m=m, c=c, u=u, D_per_pair=Dpair, a=a)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ceiling.json'), 'w'), indent=1)
