"""subst_cert.py -- ex2-modp: certificates for nag-rerun's --linsub step.  fullsys.py --linsub substitutes, before
eliminating, the degree-1 alpha-only elements of the CHART RELAXATION's Groebner basis mod p ("a8 = a7 - 1", ...).
That is sound only if the relaxation really implies the relation over C.  For each such relation x = L(alpha) this
writes the system  (chart relaxation) + {t (x - L) - 1}  (mod p1 and p2, nag-rerun's format) into OUTDIR, so that
lift_cert.py can certify it empty (=> every point of the relaxation has x = L).

  python3 subst_cert.py DECDIR OUTDIR
"""
import sys, os, re, glob, ast
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zsys import parse_ms, ratrec

P1, P2 = 1073741789, 1073741741
dec, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
seen = set()
for log in sorted(glob.glob(os.path.join(dec, 'full_c_*.log'))):
    chart = os.path.basename(log)[len('full_c_'):-len('.log')]
    txt = open(log).read()
    m1 = re.search(r"p1: linear substitutions from the relaxation: (\{.*\})", txt)
    m2 = re.search(r"p2: linear substitutions from the relaxation: (\{.*\})", txt)
    if not m1: continue
    s1, s2 = ast.literal_eval(m1.group(1)), ast.literal_eval(m2.group(1))
    for x, L1 in s1.items():
        L2 = s2[x]
        # rebuild L over Q: {var or '1': coef}
        L = {}
        for v in set(L1) | set(L2):
            r = ratrec(L1.get(v, 0), P1)
            assert r is not None and (r.numerator - r.denominator * L2.get(v, 0)) % P2 == 0, (log, x, v)
            if r: L[v] = r
        key = (chart, x, tuple(sorted(L.items())))
        if key in seen: continue
        seen.add(key)
        for j, p in ((1, P1), (2, P2)):
            names, pp, polys = parse_ms(os.path.join(dec, f'c_{chart}.p{j}.ms'))
            body = open(os.path.join(dec, f'c_{chart}.p{j}.ms')).read().split('\n', 2)[2].strip().rstrip(',')
            # row t*(x - L) - 1
            terms = [f"1*t*{x}"]
            for v, c in L.items():
                cm = (-c.numerator * pow(c.denominator, -1, p)) % p
                terms.append(f"{cm}*t" + ("" if v == '1' else f"*{v}"))
            terms.append(f"{p - 1}")
            row = "+".join(terms)
            tag = f"sub_{chart}_{x}"
            open(os.path.join(out, f"{tag}.p{j}.ms"), 'w').write(",".join(names + ['t']) + f"\n{p}\n" + body + ",\n" + row + "\n")
        print(f"{chart}: {x} = " + " + ".join(f"{c}*{v}" for v, c in L.items()) if L else f"{chart}: {x} = 0")
