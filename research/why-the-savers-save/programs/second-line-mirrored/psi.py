#!/usr/bin/env python3
"""cx-mirror / psi.py (own code, stdlib).  The mirror map psi on F_2^10 (window bits 0..8, s = bit 9):
psi(w + e s) = w + wt(w) (s + j) + e j.  Checks: linear, an involution, an isometry of the dot product, and
psi(u_b) = v_b = s + j + u_b, psi(v_b) = u_b for all 112 ports."""
import gzip, json, os
d = json.load(gzip.open(os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs/U8X-parsed.json.gz'))))
U = d['ports']; S = 1 << 9; J = (1 << 9) - 1
dot = lambda a, b: bin(a & b).count('1') & 1
def psi(x):
    w, e = x & J, (x >> 9) & 1
    y = w ^ ((S | J) if bin(w).count('1') % 2 else 0)
    return y ^ (J if e else 0)
X = range(1 << 10)
lin = all(psi(a ^ b) == psi(a) ^ psi(b) for a in X for b in range(0, 1 << 10, 37))
inv = all(psi(psi(x)) == x for x in X)
iso = all(dot(psi(a), psi(b)) == dot(a, b) for a in X for b in range(0, 1 << 10, 29))
sw = all(psi(u) == S | J ^ u and psi(S | J ^ u) == u for u in U)
print('psi linear', lin, '| involution', inv, '| isometry', iso, '| swaps u_b and v_b for all ports', sw)
