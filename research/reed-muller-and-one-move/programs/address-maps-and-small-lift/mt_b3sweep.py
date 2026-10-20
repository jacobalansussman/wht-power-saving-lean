"""mt_b3sweep.py (rm-theory): block 3, the whole no-phase game (Theorem B: targets = any of the 28 matrices congruent
to I). Random 5-in-14 instances (random targets, random direction multisets using every direction, random order)
that pass the exact tests, each with a few LM restarts; plus one positive control (plain 5 in 15 with Singer
targets) and one negative control (3 in 8 below the hook) per batch. Prints every instance's best residual."""
import sys, random
import numpy as np
from mt_lm import Setup, MTSystem
from mt_cover import cover_test
from mt_b3 import runinst
N = int(sys.argv[1]); restarts = int(sys.argv[2]); iters = int(sys.argv[3]); seed = int(sys.argv[4])
rng = random.Random(seed)
runinst("POSITIVE control", 3, "k:0,1,2,3,4", "frames", "rand%d" % seed, 3, iters, seed)
runinst("NEGATIVE control (hook)", 3, "k:0,4,5", "list:1,1,2,3,4,5,6,7", "rand%d" % seed, 3, iters, seed)
done = 0; tried = 0; floors = []
while done < N and tried < 50 * N:
    tried += 1
    ts = [rng.randrange(28) for _ in range(5)]
    comp = list(range(1, 8)) + [rng.randrange(1, 8) for _ in range(7)]
    rng.shuffle(comp)
    st = Setup(3, "gl:" + ",".join(map(str, ts)), "list:" + ",".join(map(str, comp)), "lex")
    if st.exact_tests()["excluded"] or cover_test(st)["excluded"]: continue
    b = runinst("5 in 14 GL", 3, "gl:" + ",".join(map(str, ts)), "list:" + ",".join(map(str, comp)), "lex", restarts, iters, seed + done)
    floors.append(b); done += 1
print(f"instances {done} (of {tried} drawn); best residuals: min {min(floors):.3e}, sorted head {sorted(floors)[:8]}")
