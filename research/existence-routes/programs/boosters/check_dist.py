# e3-boost check_dist.py: brute-force check of the distance-set certificate profiles (dualpolar.py) on the actual
# frames at m = 2, 3, 4 (all Lagrangians, built by transvections in certB.py): |S(L0) \ S(L_d)| for every S.
import itertools
import certB, dualpolar
for m in (2, 3, 4):
    frames = certB.all_frames(m)
    def dist(U, V): return m - (len(U & V).bit_length() - 1)
    L = [certB.L_d(m, d) for d in range(m + 1)]
    dist0 = {U: dist(L[0], U) for U in frames}
    distd = [{U: dist(L[d], U) for U in frames} for d in range(m + 1)]
    P, ks = dualpolar.intersection_numbers(m)
    bad = 0
    for r in range(1, m + 2):
        for S in itertools.combinations(range(m + 1), r):
            brute = [sum(1 for U in frames if dist0[U] in S and distd[d][U] not in S) for d in range(1, m + 1)]
            if brute != dualpolar.dist_cert(m, S, P, ks): bad += 1; print("MISMATCH", m, S, brute)
    print(f"m={m}: {len(frames)} frames; all {2 ** (m + 1) - 1} distance sets checked; mismatches {bad}")
