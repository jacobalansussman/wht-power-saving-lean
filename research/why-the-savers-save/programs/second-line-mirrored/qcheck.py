#!/usr/bin/env python3
"""cx-mirror / qcheck.py (own code, stdlib).  Is the 9-dim space Q of quadratic forms that vanish on the 112 ports the
same as the space of forms that vanish on ALL 256 odd vectors of F_2^9?  If so, u . A u = 1 on the ports forces
u . A u = 1 on every odd u, i.e. <e, A e> = 0 on the whole even space E, so A + A^-1 = 0: every twisted mirror is an
involution, A = A^T, and N is symmetric (a drop-free x chain is an independent set of {N = 1})."""
import os, importlib.util, json
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('tw', os.path.join(here, 'twist_lib.py')); tw = importlib.util.module_from_spec(spec); spec.loader.exec_module(tw)
pairs, nullspace, bit = tw.pairs, tw.nullspace, tw.bit
odd = [u for u in range(512) if bin(u).count('1') % 2 == 1]
rows = [[bit(u, i) for i in range(9)] + [bit(u, i) & bit(u, j) for i, j in pairs] for u in odd]
Qall = nullspace(rows, 45)
print('dim of forms vanishing on all 256 odd vectors:', len(Qall), '| on the 112 ports:', len(tw.Q))
def vanish_all(q):
    return all(sum(c * r for c, r in zip(q, row)) % 2 == 0 for row in rows)
same = len(Qall) == len(tw.Q) and all(vanish_all(q) for q in tw.Q)
print('the two spaces coincide:', same)
json.dump(dict(dim_Q_ports=len(tw.Q), dim_Q_all_odd=len(Qall), coincide=same), open(os.path.join(here, 'qcheck.json'), 'w'))
# the extra form Q*: an element of Q that does not vanish on every odd vector
U = tw.U
ext = next(q for q in tw.Q if not vanish_all(q))
val = lambda q, u: (sum(q[i] * bit(u, i) for i in range(9)) + sum(c * bit(u, i) * bit(u, j) for c, (i, j) in zip(q[9:], pairs))) % 2
ones = [u for u in odd if val(ext, u)]
print('Q*: nonzero on', len(ones), 'odd vectors; on ports:', sum(1 for u in U if val(ext, u)),
      '| weights of those odd vectors:', sorted({bin(u).count("1") for u in ones}))
# its polar (cross part) restricted to E = even vectors: B(e, f) = Q*(j+e+f) + Q*(j+e) + Q*(j+f) + Q*(j)
E = [e for e in range(512) if bin(e).count('1') % 2 == 0]
Eb = [1 << i | 1 << 8 for i in range(8)]   # a basis of E: e_i + e_8
B = [[(val(ext, ((1 << 9) - 1) ^ a ^ b) + val(ext, ((1 << 9) - 1) ^ a) + val(ext, ((1 << 9) - 1) ^ b) + val(ext, (1 << 9) - 1)) % 2 for b in Eb] for a in Eb]
rk = []; 
for r in [int(''.join(map(str, row)), 2) for row in B]:
    for x in rk: r = min(r, r ^ x)
    if r: rk.append(r)
print('rank of the polar form of Q* on E:', len(rk))
