"""Fixed-order frame MILP with a widened candidate pool: P/Q/current frames of both chains, all pairwise joins and
meets of those that fall in [P(k),Q(k)], and (for gaps <= 3) every subspace between P(k) and Q(k)."""
import sys, itertools
sys.path.insert(0,'.')
import fw, search, ilp
from search import le, join, meet
def all_between(Pk, Qk):
    # all subspaces U with Pk <= U <= Qk: enumerate by adding vectors of Qk outside Pk
    out = {Pk}; frontier = [Pk]
    els = [x for x in range(1, 512) if (search.mask(Qk) >> x) & 1]
    while frontier:
        nxt = []
        for U in frontier:
            mU = search.mask(U)
            for x in els:
                if not (mU >> x) & 1:
                    V = join(U, (x,))
                    if V not in out: out.add(V); nxt.append(V)
        frontier = nxt
    return out
if __name__ == "__main__":
    W = search.Word(fw.load('K8-FULL-WORD.json.gz'))
    P, Q = ilp.bounds_PQ(W)
    N = len(W.seq)
    extra = {}
    big = 0
    for k in range(N):
        if P[k] == Q[k]: continue
        gap = len(Q[k]) - len(P[k])
        s = set()
        if gap <= 3:
            s = all_between(P[k], Q[k])
        else:
            big += 1
            pool = set()
            for r in W.regs[k]:
                for k2 in W.pos[r]: pool.update((P[k2], Q[k2], W.f[k2]))
            pool = [U for U in pool if le(P[k], join(P[k], U))]
            for U in pool:
                for V in pool:
                    for X in (join(P[k], join(U, V)), meet(Q[k], join(P[k], meet(U, V)))):
                        if le(P[k], X) and le(X, Q[k]): s.add(X)
        extra[k] = s
    print('positions with widened pools', len(extra), 'gap>3', big, 'max pool', max(len(s) for s in extra.values()), flush=True)
    fr, info = ilp.solve(W, extra_pool=extra, time_limit=900)
    print(info)
    if fr:
        W.f = fr
        print('cost', W.total_cost(), dict(sorted(W.profile().items())))
