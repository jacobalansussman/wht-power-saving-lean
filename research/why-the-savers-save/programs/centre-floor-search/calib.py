#!/usr/bin/env python3
"""ctr-search / calib.py: calibration decoders.  (1) the hub's 64: the frame-root coordinate decoder in Hamming
coordinates, f = e_i (the removed frame lines);  (2) the D8-type decoder 72: orthogonal basis f = e_{2a} +- e_{2a+1}
(the 'butterfly pairs' of pr717 D4; cx-medium's e_i in D8 coordinates);  (3) a few other orthogonal bases.
Every decoder: Phi_f(b) = f.r_b, Psi_f(a) = f.r_a / (4 |f|^2), exact check in rationals, both costs."""
import sys, os, json, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p112 import *

def frame_decoder(F):
    Phi, Psi = [], []
    for f in F:
        nf = ip(f, f)
        Phi.append([ip(f, r) for r in LINES])
        Psi.append([Fraction(ip(f, r), 4 * nf) for r in LINES])
    return Phi, Psi

def report(name, Phi, Psi):
    ok, bad = check_exact(Phi, Psi)
    c = costs(Phi, Psi)
    j = junk_map(Phi, Psi)
    vals = sorted(set(str(v) for v in j.values()))
    print(f"{name}: exact {ok} (bad {bad}); n {c['n']}; pair {c['pair']} {c['pair_terms']}; star_gather {c['star_gather']} "
          f"{c['star_gather_terms']}; star_min {c['star_min']}; junk entries {len(j)} values {vals[:6]}")
    return dict(name=name, exact=ok, **c, junk_entries=len(j))

if __name__ == '__main__':
    out = []
    F1 = [tuple(1 if j == i else 0 for j in range(8)) for i in range(8)]
    out.append(report('hub64: f = e_i (frame roots)', *frame_decoder(F1)))
    F2 = []
    for a in range(4):
        for s in (1, -1):
            v = [0] * 8; v[2 * a] = 1; v[2 * a + 1] = s; F2.append(tuple(v))
    out.append(report('butterfly (e_2a +- e_2a+1)', *frame_decoder(F2)))
    F3 = []
    for a in range(4):
        for s in (1, -1):
            v = [0] * 8; v[a] = 1; v[a + 4] = s; F3.append(tuple(v))
    out.append(report('halves (e_j +- e_j+4)', *frame_decoder(F3)))
    # Hadamard 8x8 rows
    H = [tuple((-1) ** bin(i & j).count('1') for j in range(8)) for i in range(8)]
    out.append(report('WHT rows', *frame_decoder(H)))
    json.dump(out, open(os.path.join(HERE, 'calib.json'), 'w'), indent=1)
