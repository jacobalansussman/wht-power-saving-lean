# Widened game at block 2: add the one-point phase CS = diag(1,1,1,i) to the free letters; recheck VV ranks mod p.
import sys, os, json, numpy as np; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common, jac
orig = common.letters
def letters_w(m, p):
    free, Wr = orig(m, p); I, h = common.setup(p)
    free = free + [np.diag([I if x == (1 << m) - 1 else 1 for x in range(1 << m)]).astype(np.int64)]
    return free, Wr
jac.letters = letters_w
for p in (2000029, 1500041):
    for W, n in ((1, 1), (1, 2), (2, 3), (2, 4)):
        Pr = jac.Prog(2, 'VV', W, n, p, 13)
        print(json.dumps({'game': 'widened', 'p': p, 'shape': f'{W}in{n}', '|B_w|': Pr.nB, 'dimA': Pr.s, 'gU': Pr.gU, 'Eeff': Pr.Eeff,
                          'rank_modp': Pr.jac_rank(Pr.random_point())}), flush=True)
