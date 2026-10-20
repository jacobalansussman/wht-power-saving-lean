"""vcert.py -- v-ex2-modp: independent check of an ex2-modp Nullstellensatz certificate against the RAW exports.
  python3 vcert.py SYSPREFIX CERT.txt [KEEP_EXTRAS]
SYSPREFIX.p1.ms / .p2.ms are nag-rerun's (or ex2-modp's derived) mod-p exports.  If SYSPREFIX looks like
<dir>/full_c_X.a<l> and <dir>/full_c_X.log reports 'extras N' for array l at p1, the last N rows are dropped
(as ex2-modp claims to do) unless KEEP_EXTRAS is given.
Checks: (1) system rebuilt over Q by MY reconstruction (support equal at both primes, values agree);
(2) cofactors parsed by MY parser; if Singular ran in char p1, each coefficient reconstructed by MY ratrec;
(3) sum g_i f_i == 1 EXACTLY over Q; (4) independently, sum g_i f_i == 1 mod p2 on the RAW p2 export.
Prints one line: OK/FAIL with sizes, integrality and coefficient bounds."""
import sys, os, re
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vlib import *

pre, cert = sys.argv[1], sys.argv[2]
keep = len(sys.argv) > 3 and sys.argv[3] == 'KEEP_EXTRAS'
A = read_ms(pre + '.p1.ms'); B = read_ms(pre + '.p2.ms')
names, F, (mn, md) = to_Q(A, B)
PB = B[2]
drop = 0
m = re.match(r'(.*)\.a(\d+)$', pre)
if m and os.path.exists(m.group(1) + '.log') and not keep:
    mm = re.search(rf"p1 array {m.group(2)}:.*extras (\d+)", open(m.group(1) + '.log').read())
    if mm:
        drop = int(mm.group(1))
if drop:
    F = F[:-drop]; PB = PB[:-drop]
txt = open(cert).read()
blocks = re.findall(r'COF (\d+)\n(.*?)(?=\nCOF |\Z)', txt, re.S)
if not blocks:
    print(f'FAIL no cofactors in {cert}'); sys.exit(1)
sing = cert[:-4] + '.sing'
char0 = os.path.exists(sing) and re.search(r'ring r = 0,', open(sing).read()) is not None
G = {}
for i, s in blocks:
    g = parse_poly_str(s.strip(), names)
    if not char0:
        g2 = {}
        for k, v in g.items():
            assert v.denominator == 1
            r = ratrec(int(v), P1)
            if r is None:
                print('FAIL reconstruction of a cofactor coefficient'); sys.exit(1)
            g2[k] = r
        g = g2
    G[int(i)] = g
if max(G) > len(F):
    print(f'FAIL cofactor index {max(G)} > rows {len(F)} (drop={drop})'); sys.exit(1)
T = {}
for i, g in G.items():
    add_into(T, mul(g, F[i - 1]))
one = {tuple([0] * len(names)): Fraction(1)}
okQ = (T == one)
# mod p2 on the raw export
T2 = {}
for i, g in G.items():
    gp = {k: v.numerator * pow(v.denominator, -1, P2) % P2 for k, v in g.items()}
    for k, v in mul(gp, PB[i - 1]).items():
        T2[k] = (T2.get(k, 0) + v) % P2
T2 = {k: v for k, v in T2.items() if v}
okP2 = (T2 == {tuple([0] * len(names)): 1})
den = max((v.denominator for g in G.values() for v in g.values()), default=1)
used = sum(1 for g in G.values() if g)
print(f"{'OK  ' if okQ and okP2 else 'FAIL'} Q-identity={okQ} p2-raw-identity={okP2} rows={len(F)} drop={drop} "
      f"nonzero_cofactors={used} cof_maxden={den} sys_maxnum={mn} sys_maxden={md} char0={char0}  {pre}")
