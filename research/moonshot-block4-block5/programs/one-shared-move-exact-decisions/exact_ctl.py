"""exact_ctl.py -- w3-asym: positive controls for the Singular pipeline: the plain 2 x 4 program (must be CONSISTENT),
and the plain program with an extra summand in array 0's last move carrying a free coefficient (still CONSISTENT: a = 0
works), and the plain program minus array 0's last move (must be INCONSISTENT)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_decide import *
os.makedirs('sing_ctl', exist_ok=True)
P = plain(2)
s, sz = singular_script(P); open('sing_ctl/plain.sing', 'w').write(s); print("plain", sz)
rng = random.Random(3)
P2 = Prog(2)
for mv in P.moves[:3]:
    P2.add_move([(src, b) for (src, b, pidx) in mv])
mv = P.moves[3]
P2.add_move([(src, b) for (src, b, pidx) in mv] + [(1, rand_engine(rng))])
for mv in P.moves[4:]:
    P2.add_move([(src, b) for (src, b, pidx) in mv])
s, sz = singular_script(P2); open('sing_ctl/plain_extra.sing', 'w').write(s); print("plain+extra", sz, "params", P2.np)
P3 = Prog(2)
for k, mv in enumerate(P.moves):
    if k == 3: continue
    P3.add_move([(src if src < 2 + 3 else src - 1, b) for (src, b, pidx) in mv])
s, sz = singular_script(P3); open('sing_ctl/plain_minus.sing', 'w').write(s); print("plain-minus", sz)
