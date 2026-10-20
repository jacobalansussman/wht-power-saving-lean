# closed form of N_r(s), extracted verbatim from nrs.py (PROOF.md Lemma 2)
def N_closed(m, r):
    """number of nonzero p in F_2^{2m} with d(pi_p L, pi_p M) = s, for d(L, M) = r (1 <= r <= m)."""
    out = {}
    def add(s, c):
        if c: out[s] = out.get(s, 0) + c
    add(r, (2 ** (m - r) - 1) + (4 ** m - 2 ** (m + r)))            # p in L cap M, or p outside L + M
    add(r - 1, 2 * (2 ** m - 2 ** (m - r)) + (2 ** r - 1) * 2 ** (r - 1) * 2 ** (m - r))   # p in one of L, M; or p = l + m', <l,m'> = 1
    add(r - 2, (2 ** r - 1) * (2 ** (r - 1) - 1) * 2 ** (m - r))      # p = l + m', <l,m'> = 0
    assert sum(out.values()) == 4 ** m - 1
    return dict(sorted(out.items()))

