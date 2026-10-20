"""combine_cert.py -- ex2-modp: a full-chart system that has no beta-free certificate on its own may have one once
the chart's RELAXATION rows are added (every unit satisfies both: a2-nag Lemma 1).  This replaces nag-rerun's
"extras" (alpha-only elements of the relaxation's Groebner basis mod p, which carry no proof) by the relaxation
itself (exact rows).  Writes PREFIX_OUT.p1.ms / .p2.ms (nag-rerun format) with
   relaxation rows of c_<chart> (alphas a*, betas b*)  +  full-core rows of full_c_<chart>.a<l> minus its extras
   (alphas a*, betas c*),
after substituting every alpha fixed by a relaxation row of the form  k * a_i + d  (an exact equivalence).
Both blocks are linear in their betas, so lift_cert.py makem applies with beta = (b*, c*).

  python3 combine_cert.py DECDIR CHART ARRAY OUTPREFIX
"""
import sys, os, re
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zsys import load_Q, subst

dec, chart, arr, outp = sys.argv[1:5]
rn, rel = load_Q(os.path.join(dec, f'c_{chart}'))
cn, core = load_Q(os.path.join(dec, f'full_c_{chart}.a{arr}'))
k = int(re.search(rf"p1 array {arr}:.*extras (\d+)", open(os.path.join(dec, f'full_c_{chart}.log')).read()).group(1))
core = core[:len(core) - k] if k else core
names = list(rn) + [v for v in cn if v not in rn]
def remap(P, src):
    out = {}
    for e, c in P.items():
        f = [0] * len(names)
        for i, x in enumerate(e):
            f[names.index(src[i])] += x
        out[tuple(f)] = c
    return out
polys = [remap(P, rn) for P in rel] + [remap(P, cn) for P in core]
# substitute alphas fixed by single-variable linear rows
changed = True; fixed = {}
while changed:
    changed = False
    for P in polys:
        vs = {i for e in P for i, x in enumerate(e) if x}
        if len(vs) == 1:
            i = vs.pop()
            if names[i][0] == 'a' and max(e[i] for e in P) == 1:
                a1 = sum(c for e, c in P.items() if e[i] == 1); a0 = sum(c for e, c in P.items() if e[i] == 0)
                val = Fraction(-a0) / a1; fixed[names[i]] = val
                polys = [subst(Q, i, val) for Q in polys]
                changed = True; break
polys = [P for P in polys if P]
if any(set(P) == {tuple([0] * len(names))} for P in polys):
    print('a row became a nonzero constant: empty by substitution alone')
used = sorted({i for P in polys for e in P for i, x in enumerate(e) if x})
nm = [names[i] for i in used]
polys = [{tuple(e[i] for i in used): c for e, c in P.items()} for P in polys]
for j, p in ((1, 1073741789), (2, 1073741741)):
    lines = []
    for P in polys:
        t = []
        for e, c in P.items():
            cm = c.numerator * pow(c.denominator, -1, p) % p
            mon = '*'.join(f"{nm[i]}^{x}" if x > 1 else nm[i] for i, x in enumerate(e) if x)
            t.append(f"{cm}" + ('*' + mon if mon else ''))
        lines.append('+'.join(t))
    open(f"{outp}.p{j}.ms", 'w').write(','.join(nm) + f"\n{p}\n" + ',\n'.join(lines) + '\n')
print(f"{chart}.a{arr}: {len(polys)} rows, {sum(1 for v in nm if v[0]=='a')} alphas, {sum(1 for v in nm if v[0]!='a')} betas; fixed {fixed}")
