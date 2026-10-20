import sys; sys.path.insert(0, "."); from indep import *
m = 3; Ls = lagrangians(m); tot = 0; it = 0
gauss = {0: 1, 1: 2, 2: 5, 3: 16}
want = 0
for X in Ls:
    for Y in Ls:
        d = dist(X, Y, m); want += gauss[d]
        for Z in Ls:
            it += 1
            if dist(X, Z, m) + dist(Z, Y, m) == d: tot += 1
print("iterations", it, "interval members", tot, "expected sum of #subspaces(F_2^d)", want)
