import re, ast, sys
from fractions import Fraction as Fr
from vnulla0 import load, evaluate
txt = open('../out/cal1.txt').read()
blocks = re.split(r'^=== ', txt, flags=re.M)[1:]
n = 0
for b in blocks:
    name = b.split()[0]
    m = re.search(r'exact_point: (\{.*\})', b)
    v = re.search(r'verdict: (.*)', b)
    if not m: print(name, 'no exact point;', v.group(1)[:60] if v else ''); continue
    pt = ast.literal_eval(m.group(1)); d, eqs = load('../sys/' + name)
    P = {int(k[1:]): (Fr(a), Fr(c)) for k, (a, c) in pt.items()}
    assert sorted(P) == list(range(d['nvar']))
    bad = evaluate(eqs, P); n += bad == 0
    print(name, 'exact point: %d of %d equations nonzero' % (bad, len(eqs)))
print('points solving every equation exactly:', n)
