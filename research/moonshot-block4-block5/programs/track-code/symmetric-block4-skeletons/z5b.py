from symsolve import *
R = Run(4, 5, 'Z5rot', 1); G = R.G; rng = random.Random(3)
fd = fixed_dirs(R, rng, 20000)
print("|H| =", len(R.Hc), " Z5 fixed rank-1 directions in 20000 samples:", len(fd))
