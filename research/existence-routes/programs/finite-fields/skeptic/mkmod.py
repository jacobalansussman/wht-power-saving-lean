"""mkmod.py -- v-ex2-modp: write a Singular script for a beta-free (module) certificate of an EXACT core file (fq/*.Q.txt)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vlib import parse_poly_str
from fractions import Fraction
src, out, ch = sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else '0'
L = open(src).read().split('\n'); names = [v for v in L[0].split(',') if v]
polys = [parse_poly_str(c, names) for c in '\n'.join(L[2:]).split(',\n') if c.strip()]
ai = [i for i, v in enumerate(names) if v[0] == 'a']; bi = [i for i, v in enumerate(names) if v[0] != 'a']
def sp(P):
    t = []
    for k, c in P.items():
        mon = '*'.join(f"{names[i]}^{k[i]}" if k[i] > 1 else names[i] for i in ai if k[i])
        cs = f"({c.numerator}/{c.denominator})" if c.denominator != 1 else f"({c.numerator})"
        t.append(cs + ('*' + mon if mon else ''))
    return '+'.join(t) if t else '0'
gens = []
for P in polys:
    comps = []
    for b in bi:
        Q = {k: c for k, c in P.items() if k[b] == 1}; assert all(k[b] <= 1 for k in P)
        comps.append(sp(Q))
    comps.append(sp({k: c for k, c in P.items() if all(k[b] == 0 for b in bi)}))
    gens.append('[' + ','.join(comps) + ']')
an = [names[i] for i in ai]
s = ["short = 0;", f"ring r = {ch}, ({','.join(an)}), dp;", "option(redSB);", "module M = " + ",\n ".join(gens) + ";",
     "module E = [" + ','.join(['0'] * len(bi) + ['1']) + "];", "module G = std(M);",
     "if (size(reduce(E, G)) == 0) {", " matrix T = lift(M, E);",
     f" int i; for (i = 1; i <= {len(polys)}; i++) {{ print(\"COF \" + string(i)); print(string(T[i,1])); }}",
     "} else { print(\"NOT IN MODULE\"); }", "quit;"]
open(out, 'w').write('\n'.join(s) + '\n'); print('wrote', out, len(polys), 'rows', len(bi), 'betas', an)
