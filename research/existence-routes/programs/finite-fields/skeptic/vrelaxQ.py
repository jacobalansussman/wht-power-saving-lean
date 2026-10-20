"""vrelaxQ.py -- v-ex2-modp: build each relaxation EXACTLY (nag-rerun's myrelax.Skeleton, integer coefficients, no prime)
and compare it with the system ex2-modp used (raw mod-p export reconstructed over Q).  Run with numpy (uv)."""
import sys, os, json
NR = os.environ.get('NAG_RERUN_DIR', 'nag-rerun')
sys.path.insert(0, NR); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from myrelax import Skeleton
from vlib import *
from fractions import Fraction
M = os.environ.get('MODP_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
def canonQ(names, P):
    return tuple(sorted((tuple(sorted((names[i], e) for i, e in enumerate(k) if e)), Fraction(c)) for k, c in P.items()))
names_list = sorted({f.split('.u.p1.ms')[0] for f in os.listdir(f'{M}/data/ms') if f.endswith('.u.p1.ms')})
tot = ok = 0
for name in names_list:
    rp = f'{NR}/recipes/{name}.json'
    if not os.path.exists(rp):
        print(name, 'no recipe'); continue
    rec = json.load(open(rp)); sk = Skeleton(rec, normalise=False)
    rows, Kkey = sk.relaxation(); polys, nf = sk.polys(rows, Kkey)
    nb = len(sk.gens); vn = [f"a{k+1}" for k in range(nf)] + [f"b{k+1}" for k in range(nb)]
    ex = set(); maxc = 0
    for key, P in polys.items():
        if not P: continue
        Q = {}
        for mono, c in P.items():
            assert c.imag == 0 and c.real == int(c.real)
            e = [0] * len(vn)
            for v in mono: e[v] += 1
            Q[tuple(e)] = Q.get(tuple(e), 0) + int(c.real); maxc = max(maxc, abs(int(c.real)))
        Q = {k: v for k, v in Q.items() if v}
        if Q: ex.add(canonQ(vn, Q))
    rn, F, _ = to_Q(read_ms(f'{M}/data/ms/{name}.u.p1.ms'), read_ms(f'{M}/data/ms/{name}.u.p2.ms'))
    raw = set(canonQ(rn, P) for P in F)
    tot += 1; same = (raw == ex) and rn == vn; ok += same
    print(f"{name}.u: exact rows {len(ex)} raw rows {len(F)} (distinct {len(raw)}) max|coef| {maxc} -> {'IDENTICAL' if same else 'DIFFERENT'}", flush=True)
print(f'TOTAL {ok}/{tot} relaxation exports identical to the exact build')
