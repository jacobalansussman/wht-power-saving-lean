"""How far the Lean-plannable rate (plan_block, fill s, Dd log decimals) sits below the exact moment root."""
import sys, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Q
from exact_moment import certify
import foldlib as fl
from exit_price import ceil_unit
prof = {1:2548,2:338,3:446,4:162,5:72,8:8}
w = [4]*80+[5]*32+[7]*54+[8]*122+[9]*216
u = ceil_unit(prof, w)
root = Q(certify(u['hist'], 73, u['W'])['root_lower'])
print('U8 ceil W', u['W'], 'D', u['D'], 'root', float(root))
uu = dict(m=73, W=u['W'], D=u['D'], H=u['hist'])
for s, Dd in [(40,12),(40,16),(60,12),(60,16),(80,20)]:
    A = int(root*10**10)
    while fl.plan_block(uu, A, 10, s, Dd) is None: A -= 1
    print('s', s, 'Dd', Dd, 'kernel A', A, 'gap to root grid', int(root*10**10) - A)
