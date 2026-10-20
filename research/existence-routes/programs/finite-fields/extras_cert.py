"""extras_cert.py -- ex2-modp: certify nag-rerun's appended 'extras' (alpha-only elements of the chart relaxation's
Groebner basis mod p).  For each extra e (the last k rows of full_c_<chart>.a<l>), write the system
c_<chart> + {t e - 1} (nag-rerun format, both primes): if lift_cert certifies it empty, then e vanishes at every point
of the relaxation, hence at every unit (a2-nag Lemma 1), and the extras may be used as exact rows.

  python3 extras_cert.py DECDIR CHART ARRAY OUTDIR"""
import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zsys import parse_ms
dec, chart, arr, out = sys.argv[1:5]
os.makedirs(out, exist_ok=True)
k = int(re.search(rf"p1 array {arr}:.*extras (\d+)", open(os.path.join(dec, f'full_c_{chart}.log')).read()).group(1))
for j in (1, 2):
    txt = open(os.path.join(dec, f'full_c_{chart}.a{arr}.p{j}.ms')).read().strip().split('\n')
    names = txt[0].split(','); p = int(txt[1]); rows = '\n'.join(txt[2:]).split(',\n')
    extras = [r.strip().rstrip(',') for r in rows[len(rows) - k:]]
    rel = open(os.path.join(dec, f'c_{chart}.p{j}.ms')).read().strip().split('\n')
    rnames = rel[0].split(','); rbody = '\n'.join(rel[2:]).strip().rstrip(',')
    for q, e in enumerate(extras):
        # t * e - 1 : multiply every term of e by t
        terms = [t + '*t' for t in e.split('+') if t]
        row = '+'.join(terms) + f'+{p - 1}'
        open(os.path.join(out, f'ext_{chart}_a{arr}_{q}.p{j}.ms'), 'w').write(','.join(rnames + ['t']) + f'\n{p}\n' + rbody + ',\n' + row + '\n')
print(chart, 'extras:', k)
