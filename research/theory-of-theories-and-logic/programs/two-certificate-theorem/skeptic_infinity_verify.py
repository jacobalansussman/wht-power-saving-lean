"""infinity_verify.py -- v-ml-fooling: top-degree forms of the K3 and W3 coordinate slices (from the certs), Bezout
number, and an exact check of the hand decomposition of the locus at infinity (see VERDICT.md)."""
import json, sys, random
from fractions import Fraction
for fn in ['K3_coord.json', 'W3_coord.json']:
    C = json.load(open('../../data/two-certificate-theorem/certs/' + fn)); names = C['unknowns']
    tops = []; bez = 1
    for f in C['equations']:
        d = max(sum(e) for e, c in f); bez *= d
        tops.append([(e, Fraction(c)) for e, c in f if sum(e) == d])
    print(fn, 'Bezout', bez)
    for t in tops:
        print('   top:', ' + '.join(f"{c}*" + '*'.join(f"{names[i]}" + (f"^{k}" if k > 1 else '') for i, k in enumerate(e) if k) for e, c in t))
    # family B: u0=u1=a, v0=v1=b, w1 = -w0, others free; family A: w0 = w1 = 0
    def ev(t, x): 
        s = 0
        for e, c in t:
            p = c
            for i, k in enumerate(e): p *= x[i] ** k
            s += p
        return s
    rng = random.Random(1)
    for fam in 'AB':
        ok = True
        for _ in range(50):
            x = {v: Fraction(rng.randint(-99, 99)) for v in names}
            if fam == 'B':
                x['u1_1'] = x['u0_1']; x['v1_1'] = x['v0_1']
                for v in names:
                    if v.startswith('w1_'): x[v] = -x['w0_' + v[3:]]
            else:
                for v in names:
                    if v.startswith('w0_') or v.startswith('w1_'): x[v] = Fraction(0)
            xx = [x[v] for v in names]
            ok &= all(ev(t, xx) == 0 for t in tops)
        free = len(names) - (2 + sum(v.startswith('w1_') for v in names)) if fam == 'B' else len(names) - sum(v.startswith('w0_') or v.startswith('w1_') for v in names)
        print(f"   family {fam}: all top forms vanish on 50 random points: {ok}; affine cone dim {free} -> projective dim {free-1}")
