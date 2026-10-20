# Own implementation of Corollary 10' (even test): positions (i, j) of moves lying in a common tight odd support are
# forced zeros of E; an even support whose bipartite graph D x D minus forced zeros has no perfect matching excludes.
import geo, sys, numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
def even_test(m, W, dirs):
    forms = geo.all_form_supports(m)
    M = len(dirs)
    Z = np.zeros((M, M), bool)
    for S, tr in forms:
        D = [i for i, z in enumerate(dirs) if z in S]
        if tr == 1 and len(D) < W:
            return "EXCLUDED (deletion)"
        if tr == 1 and len(D) == W:
            Z[np.ix_(D, D)] = True
    fails = 0; first = None
    for S, tr in forms:
        if tr: continue
        D = [i for i, z in enumerate(dirs) if z in S]
        if not D: continue
        allowed = ~Z[np.ix_(D, D)]
        mt = maximum_bipartite_matching(csr_matrix(allowed.astype(np.int8)), perm_type="column")
        if (mt >= 0).sum() < len(D):
            fails += 1
            if first is None: first = sorted(dirs[i] for i in D)
    return ("EXCLUDED: %d even supports without perfect matching, first %s" % (fails, first)) if fails else "passes"
U4 = list(range(1, 16))
cases = [("U4 + {7,8}  (4 in 17)", 4, 4, U4 + [7, 8]), ("U4 + {4,11} (4 in 17)", 4, 4, U4 + [4, 11]),
         ("U4 + {2,13} (4 in 17)", 4, 4, U4 + [2, 13]), ("U4 + {1,14} (4 in 17)", 4, 4, U4 + [1, 14]),
         ("U4 + {3,12} (4 in 17, not in their list)", 4, 4, U4 + [3, 12]),
         ("4 x (7,11,13,14) + 15 (4 in 17, exists)", 4, 4, [7, 11, 13, 14] * 4 + [15]),
         ("U3 lifted + e4 e4 (2 in 9, exists: own solver)", 4, 2, [1, 2, 3, 4, 5, 6, 7, 8, 8]),
         ("plain 3 in 12", 4, 3, [1, 2, 4, 8] * 3),
         ("survivor W=12 M=47: 3xU4 + {7,8}", 4, 12, U4 * 3 + [7, 8])]
for name, m, W, dirs in cases:
    print("%-48s %s" % (name, even_test(m, W, dirs)))
