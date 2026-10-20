"""exact integer rank (fraction-free Bareiss), shared by rm2-phase scripts."""
def rank_int(rows):
    M = [list(r) for r in rows if any(r)]
    if not M: return 0
    n = len(M[0]); r = 0; prev = 1
    for c in range(n):
        p = next((k for k in range(r, len(M)) if M[k][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        for k in range(r + 1, len(M)):
            M[k] = [(M[r][c] * M[k][j] - M[k][c] * M[r][j]) // prev for j in range(n)]
        prev = M[r][c]; r += 1
        if r == len(M): break
    return r

