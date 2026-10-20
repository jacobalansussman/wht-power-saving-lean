"""levi.py -- test of the shared-move coset lemma: for the shared first move v = delta sum_l tau^l(b0 x0), the K-frame
cosets of different arrays coincide when b0^-1 h^l b0 lies in B^K = K B K^-1 (then a unit of shape one-shared-move +
private chains is impossible: res^2 >= W-1).  Count, over progress letters b0 for fixed directions delta, how many are
Levi-breaking (b0^-1 h b0 not in B^K for every power), per symmetry group."""
import sys
from symsolve import *
for (m, W, Hn) in [(4, 2, 'Z2swap'), (4, 2, 'Z2swapNL'), (4, 3, 'Z3rot'), (4, 3, 'Z3T_rot'), (4, 4, 'Z4rot'), (4, 7, 'Z7rot'), (4, 6, 'Z6rot'), (3, 2, 'Z2swap'), (3, 3, 'Z3rot')]:
    R = Run(m, W, Hn, 1); G = R.G; rng = random.Random(11)
    Ki = np.conj(G.K).T
    fd = fixed_dirs(R, rng, 4000)
    tot = brk = 0
    for d in fd:
        pl = G.progress_letters(G.ID, rng, want=24, Dmat=d, steps=1)
        for b in pl:
            tot += 1
            ok = True
            for (c, pi, cp) in R.Hc[1:]:
                x = G.mul(np.conj(b).T, c, b)
                if G.is_monomial(G.mul(Ki, x, G.K)):       # x in K^-1 B K = B^K
                    ok = False; break
            brk += ok
    print(f"m={m} W={W} H={Hn}: fixed dirs {len(fd)}, progress letters on inputs {tot}, Levi-breaking {brk}", flush=True)
