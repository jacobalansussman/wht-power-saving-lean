"""classes.py: orbit representatives of W-subsets of Z_n (canonical Singer copy indices k, target mult(alpha^2k))
under k -> 2^j k + c (Singer cycle and Frobenius; they act on the canonical targets and fix the direction set)."""
import itertools, sys
def classes(m, W):
    n = (1 << m) - 1
    maps = [(pow(2, j, n), c) for j in range(m) for c in range(n)]
    seen = set(); reps = []
    for S in itertools.combinations(range(n), W):
        if S in seen: continue
        orb = set()
        for (a, c) in maps:
            orb.add(tuple(sorted((a * k + c) % n for k in S)))
        seen |= orb; reps.append((min(orb), len(orb)))
    return reps
if __name__ == "__main__":
    m, W = int(sys.argv[1]), int(sys.argv[2])
    reps = classes(m, W)
    print(len(reps), "classes; orbit sizes sum", sum(o for _, o in reps))
    for r, o in reps: print(",".join(map(str, r)), o)
