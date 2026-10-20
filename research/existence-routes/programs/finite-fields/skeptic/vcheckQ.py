"""vcheckQ.py CORE.Q.txt CERT.txt : exact check of 1 = sum g_i f_i on an exact core file (rational cofactors, char 0)."""
import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vlib import *
from fractions import Fraction
L = open(sys.argv[1]).read().split('\n'); names = [v for v in L[0].split(',') if v]
F = [parse_poly_str(c, names) for c in '\n'.join(L[2:]).split(',\n') if c.strip()]
blocks = re.findall(r'COF (\d+)\n(.*?)(?=\nCOF |\Z)', open(sys.argv[2]).read(), re.S)
T = {}
for i, s in blocks:
    add_into(T, mul(parse_poly_str(s.strip(), names), F[int(i) - 1]))
print('OK  exact identity 1 = sum g_i f_i' if T == {tuple([0] * len(names)): Fraction(1)} else 'FAIL', len(blocks), 'cofactors', sys.argv[1])
