import sys
from symsolve import *
m, W, Hn = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
R = Run(m, W, Hn, 1); G = R.G; rng = random.Random(5)
fd = fixed_dirs(R, rng, 3000)
print(Hn, "|H|", len(R.Hc), "fixed dirs found:", len(fd))
for d in fd[:6]:
    pl = G.progress_letters(G.ID, rng, want=3, Dmat=d, steps=1)
    dist_d = G.dist(G.intern(d))
    print("  dist(d)=", dist_d, "progress letters on input:", len(pl))
