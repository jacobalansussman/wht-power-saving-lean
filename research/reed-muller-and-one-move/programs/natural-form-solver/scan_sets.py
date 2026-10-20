"""scan_sets.py: exact tests (span, deletion, tight forms) for structured target sets at block 5, one move per
direction (U_5 with per-array canonical Singer targets T_k)."""
import itertools, warnings
import numpy as np
warnings.filterwarnings("ignore")
from inst import canonical, Field
fd = Field(5)
F = fd.F
def test(ks):
    inst, _ = canonical(5, list(ks), dirs="all", order="lex")
    return inst.exact_tests()
cos = []
seen = set()
for a in range(31):
    if a in seen: continue
    c = sorted({(a * 2 ** s) % 31 for s in range(5)}); cos.append(c); seen |= set(c)
print("cyclotomic cosets:", cos)
# Frobenius-invariant sets
for r in range(1, 3):
    for cs in itertools.combinations(cos[1:], r):
        for z0 in (0, 1, 2, 3):
            ks = [0] * z0 + [k for c in cs for k in c]
            if 8 <= len(ks) <= 10:
                print("frob", ks, test(ks), flush=True)
# additive: 3-dim F2-subspaces V of GF(32) inside {alpha^(2k)} u {0}: targets T_k with alpha^(2k) in V
elts = list(range(1, 32))
subs = set()
for b1, b2, b3 in itertools.combinations(elts, 3):
    V = {0}
    for c in itertools.product([0, 1], repeat=3):
        v = (b1 if c[0] else 0) ^ (b2 if c[1] else 0) ^ (b3 if c[2] else 0); V.add(v)
    if len(V) == 8: subs.add(frozenset(V))
print("3-dim subspaces:", len(subs))
def k_of(beta):   # alpha^(2k) = beta  -> k = log(beta)/2 mod 31
    return (F.log[beta] * 16) % 31
done = 0
for V in sorted(subs, key=sorted)[:6]:
    ks7 = sorted(k_of(b) for b in V if b)
    for extra in [ks7[0], (ks7[0] + 1) % 31]:
        ks = ks7 + [extra]
        print("addsub", ks, test(ks), flush=True)
