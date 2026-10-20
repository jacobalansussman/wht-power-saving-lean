import numpy as np, itertools
from vb import *
for m in (2, 4):
    S = [t for t in sym_nonsing(m) if not t.diagonal().any()][0]
    pts = [vec(z, m) for z in range(1, 2 ** m)]
    forms = set()
    for a in range(1, 2 ** m):
        for b in range(1, 2 ** m):
            va, vb_ = vec(a, m), vec(b, m)
            if (va @ S @ vb_) % 2 == 1:
                forms.add(tuple(int((va @ p) % 2 * (vb_ @ p) % 2) for p in pts))
    forms = list(forms); F = np.array(forms)
    wts = set(F.sum(1).tolist()); cover = F.sum(0)
    # parity <Q,S> via a decomposition S = sum v v^T: find a subset of points whose rank-ones sum to S
    dec = None
    for k in range(1, m + 3):
        for sub in itertools.combinations(range(len(pts)), k):
            if (sum(np.outer(pts[i], pts[i]) for i in sub) % 2 == S).all(): dec = sub; break
        if dec: break
    par = set(int(F[q][list(dec)].sum() % 2) for q in range(len(F)))
    nQ = len(F); c = cover[0]
    print(f"m={m}: {nQ} distinct forms (a.z)(b.z) with a^T S b = 1; weights {wts}; parities on S {par}; "
          f"min decomposition {len(dec)} rank-ones; direction cover {set(cover.tolist())}; "
          f"deletion bound per array = nQ/cover = {nQ / c} vs (2^m-1)/2^(m-2) = {(2 ** m - 1) / 2 ** (m - 2)}")
