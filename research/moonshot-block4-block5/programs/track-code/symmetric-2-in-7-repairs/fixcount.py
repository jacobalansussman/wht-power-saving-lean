"""fixcount.py -- number of tau-fixed rank-1 move directions (frames of delta = b D1 b^-1, b free, with
c delta c^-1 delta^-1 free for all c in H), sampled; and how many of them make progress on an input term."""
import sys
from symsolve import *
for (m, W, Hn) in [(4, 2, 'Z7'), (4, 2, 'Z2swap'), (4, 2, 'Z2swapNL'), (4, 2, 'Z3fixT'), (4, 3, 'Z3rot'), (4, 3, 'Z3T_rot'), (4, 4, 'Z4rot'), (4, 4, 'V4'), (4, 2, 'GL32')]:
    try:
        R = Run(m, W, Hn, 1)
    except Exception as e:
        print(Hn, "skip", e); continue
    G = R.G; rng = random.Random(7)
    fd = fixed_dirs(R, rng, int(sys.argv[1]) if len(sys.argv) > 1 else 3000)
    prog = sum(1 for d in fd if G.progress_letters(G.ID, rng, want=2, Dmat=d, steps=1))
    print(f"m={m} H={Hn} |H|={len(R.Hc)}: fixed rank-1 directions (frames of delta^-1 L0) found {len(fd)}, with progress on an input: {prog}", flush=True)
