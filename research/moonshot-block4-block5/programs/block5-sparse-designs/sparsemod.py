"""sparsemod.py (w8-b5sparse): exact SPARSE consistency test of A w = e over F_p (structured Gaussian elimination).

A is given as {(r, c): v} (v in F_p), e as {r: v}.  Reductions, all exact and consistency-preserving:
  * a column with a single nonzero (in row r): row r can always be met by that column -> delete row r and the column;
  * a row with a single nonzero (column c): w_c is determined -> substitute into the other rows of c, delete both;
  * an empty row: inconsistent iff its rhs is nonzero;
then general elimination with the Markowitz pivot (min (row count - 1) * (col count - 1)).
Returns (consistent, rank, stats)."""
import heapq


def consistent_sparse(A, e, p, max_fill=50_000_000, dense_limit=4e7):
    rows = {}; cols = {}
    for (r, c), v in A.items():
        v %= p
        if v:
            rows.setdefault(r, {})[c] = v; cols.setdefault(c, set()).add(r)
    rhs = {r: v % p for r, v in e.items() if v % p}
    for r in rhs:
        rows.setdefault(r, {})
    rank = 0; fill = 0
    inv = lambda x: pow(x, p - 2, p)

    def drop_row(r):
        for c in rows[r]:
            s = cols.get(c)
            if s is not None:
                s.discard(r)
        del rows[r]; rhs.pop(r, None)

    def drop_col(c):
        for r in cols.get(c, ()):
            rows[r].pop(c, None)
        cols.pop(c, None)

    def pivot(r, c):
        """use row r to eliminate column c from every other row"""
        nonlocal fill
        R = rows[r]; a = R[c]; ia = inv(a); br = rhs.get(r, 0)
        for r2 in list(cols[c]):
            if r2 == r:
                continue
            R2 = rows[r2]; f = R2[c] * ia % p
            for c2, v in R.items():
                if c2 == c:
                    continue
                nv = (R2.get(c2, 0) - f * v) % p
                if nv:
                    if c2 not in R2:
                        fill += 1; cols[c2].add(r2)
                    R2[c2] = nv
                elif c2 in R2:
                    del R2[c2]; cols[c2].discard(r2)
            del R2[c]
            if br:
                nb = (rhs.get(r2, 0) - f * br) % p
                if nb:
                    rhs[r2] = nb
                else:
                    rhs.pop(r2, None)
        cols[c] = {r}
        drop_row(r); cols.pop(c, None)

    changed = True
    while rows:
        if fill > max_fill:
            return None, rank, dict(fill=fill, left_rows=len(rows), left_cols=len(cols))
        # empty rows
        for r in [r for r, R in rows.items() if not R]:
            if rhs.get(r, 0):
                return False, rank, dict(fill=fill, reason='empty row with nonzero rhs', left_rows=len(rows))
            drop_row(r)
        if not rows:
            break
        # singleton columns: free -> drop row + column
        sc = [c for c, s in cols.items() if len(s) == 1]
        if sc:
            for c in sc:
                s = cols.get(c)
                if s is None or len(s) != 1:
                    continue
                r = next(iter(s))
                drop_row(r); cols.pop(c, None); rank += 1
            continue
        for c in [c for c, s in cols.items() if not s]:
            cols.pop(c, None)
        # singleton rows
        sr = [r for r, R in rows.items() if len(R) == 1]
        if sr:
            for r in sr:
                if r not in rows or len(rows[r]) != 1:
                    continue
                c = next(iter(rows[r])); pivot(r, c); rank += 1
            continue
        # small enough: finish with FLINT on the dense core
        if dense_limit and len(rows) * (len(cols) + 1) <= dense_limit:
            import flint
            rl = list(rows); cl_ = [c for c in cols if cols[c]]
            ri = {r: i for i, r in enumerate(rl)}; ci = {c: j for j, c in enumerate(cl_)}
            Ma = flint.nmod_mat(len(rl), max(1, len(cl_)), p); Mt = flint.nmod_mat(len(rl), len(cl_) + 1, p)
            for r, R in rows.items():
                for c, v in R.items():
                    Ma[ri[r], ci[c]] = v; Mt[ri[r], ci[c]] = v
                if rhs.get(r):
                    Mt[ri[r], len(cl_)] = rhs[r]
            ra = Ma.rank() if cl_ else 0; rt = Mt.rank()
            return ra == rt, rank + ra, dict(fill=fill, dense_core=(len(rl), len(cl_)))
        # pivot on a column of minimum count (lazy heap), row of minimum length in it; a batch of pivots per pass
        heap = [(len(s), c) for c, s in cols.items() if s]
        heapq.heapify(heap)
        npiv = 0
        while heap and npiv < 2000:
            cnt, c = heapq.heappop(heap)
            s = cols.get(c)
            if not s:
                continue
            if len(s) != cnt:
                heapq.heappush(heap, (len(s), c)); continue
            if len(s) == 1:
                break                       # back to the cheap reductions
            r = min(s, key=lambda rr: len(rows[rr]))
            pivot(r, c); rank += 1; npiv += 1
            if fill > max_fill:
                break
    return True, rank, dict(fill=fill)


if __name__ == "__main__":
    import random
    p = 1000003
    rng = random.Random(1)
    # self-test against dense elimination
    for t in range(200):
        nr, nc = rng.randrange(1, 9), rng.randrange(1, 9)
        A = {(r, c): rng.randrange(p) for r in range(nr) for c in range(nc) if rng.random() < 0.4}
        e = {r: rng.randrange(p) for r in range(nr) if rng.random() < 0.5}
        if rng.random() < 0.5:   # make it consistent
            w = [rng.randrange(p) for _ in range(nc)]
            e = {}
            for (r, c), v in A.items():
                e[r] = (e.get(r, 0) + v * w[c]) % p
        try:
            import flint  # noqa
            dl = 4e7
        except ImportError:
            dl = 0
        ok, rk, st = consistent_sparse(dict(A), dict(e), p, dense_limit=dl if t % 2 else 0)
        # dense check
        M = [[A.get((r, c), 0) for c in range(nc)] + [e.get(r, 0)] for r in range(nr)]
        def rank(Mx, ncol):
            Mx = [row[:ncol] for row in Mx]; rk_ = 0; R = len(Mx)
            for c in range(ncol):
                pr = next((i for i in range(rk_, R) if Mx[i][c] % p), None)
                if pr is None:
                    continue
                Mx[rk_], Mx[pr] = Mx[pr], Mx[rk_]; iv = pow(Mx[rk_][c], p - 2, p)
                Mx[rk_] = [x * iv % p for x in Mx[rk_]]
                for i in range(R):
                    if i != rk_ and Mx[i][c] % p:
                        f = Mx[i][c]; Mx[i] = [(x - f * y) % p for x, y in zip(Mx[i], Mx[rk_])]
                rk_ += 1
            return rk_
        okd = rank([r[:] for r in M], nc) == rank([r[:] for r in M], nc + 1)
        assert ok == okd, (t, ok, okd, A, e)
    print("sparsemod self-test: 200 random systems agree with dense elimination")
