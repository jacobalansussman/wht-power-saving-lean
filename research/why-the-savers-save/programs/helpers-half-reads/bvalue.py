"""d-halfread: after the merge a += k*b at P, how long does the partner b keep holding exactly beta (+ its dirt)?
For each of the 168 planes: time of b's first write after the merge vs. the times at which h3's readers read h3.
Usage: python3 -I bvalue.py WORD copies.json"""
import sys, json
from collections import Counter
sys.path.insert(0, __import__('os').path.dirname(__file__))
from wordlib import *
w = load(sys.argv[1]); v, R = w['report']['v'], w['report']['R']; n = 2*v+R
prog = program(w); C = json.load(open(sys.argv[2]))
first_write_after = {}
st = Counter(); st2 = Counter()
for r in C:
    Tm, phm, a, b, k = r['merge']
    fw = next((T for T, (ph, e) in enumerate(prog) if T > Tm and e[1] == b and ph != 'compensation'), None)
    rd = [T for T, ph, c, cl, qq, d in r['reads'] if T > r['copy'][0] and T < r['Tend']]
    ph_fw = prog[fw][0] if fw is not None else None
    cls_b = 'x' if b < v else 'h'
    ok = [T for T in rd if fw is None or T < fw]
    st[(cls_b, 'b rewritten in ' + str(ph_fw), 'reads of h3 before that: %d of %d' % (len(ok), len(rd)))] += 1
for s, x in sorted(st.items(), key=lambda z: -z[1]): print(x, s)
