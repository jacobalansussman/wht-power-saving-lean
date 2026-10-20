"""t_units_cases.py: the frame-game test programs of t_units.py (same construction, same random seeds), as a module.
cases(E, rng) -> list of (name, FrameProg)."""
from lift import *

def cases(E, rng):
    S = E.S; m = E.m
    def kchain(P, r, order=None):
        for j in (order or range(m)): r = P.move(1 << j, r)
        return r
    def rand_path(P, r):
        """random free letters between m moves, steered: accept only if it ends at Kbar up to a free letter"""
        for _ in range(20000):
            ps = [S.rand_P(rng) for _ in range(m + 1)]; zs = [rng.randrange(1, E.n) for _ in range(m)]
            F = ps[0]
            for t in range(m): F = S.mul(ps[t + 1], S.mul(S.T(zs[t]), F))
            if S.in_P(S.mul(S.K, S.inv(F))): break
        else: raise RuntimeError("no path")
        c = P.mix([(1, ps[0], r)])
        for t in range(m):
            c = P.move(zs[t], c); c = P.mix([(1, ps[t + 1], c)])
        return c

    def cases():
        out = []
        P = FrameProg(S, 1); y = P.inp(0); kchain(P, y); out.append(("plain1", P))
        P = FrameProg(S, 2); y1, y2 = P.inp(0), P.inp(1); kchain(P, y1); kchain(P, y2); out.append(("plain2", P))
        sw = S.addr([2, 1] + [1 << k for k in range(2, m)])          # swap e0 <-> e1
        assert S.mul(sw, S.K) == S.mul(S.K, sw)
        P = FrameProg(S, 2); y1, y2 = P.inp(0), P.inp(1)
        r1 = P.mix([(1, S.I, y1), (1, sw, y2)]); r2 = P.mix([(1, S.I, y2), (-1, sw, y1)])
        kchain(P, r1); kchain(P, r2, order=list(range(m))[::-1]); out.append(("xchg", P))
        if m >= 3:
            s3 = S.addr([1 << ((k + 1) % 3) if k < 3 else 1 << k for k in range(m)])
            P = FrameProg(S, 2); y1, y2 = P.inp(0), P.inp(1)
            r1 = P.mix([(1, S.I, y1), (1, S.I, y2)]); r2 = P.mix([(1, S.I, y1), (-1, s3, y2)])
            kchain(P, r1); kchain(P, r2); out.append(("o3", P))
        P = FrameProg(S, 2); y1, y2 = P.inp(0), P.inp(1)
        r1 = P.mix([(1, S.I, y1), (1, S.I, y2)]); r2 = P.mix([(1, S.I, y1), (-2, sw, y2)])
        rand_path(P, r1); rand_path(P, r2); out.append(("mix13", P))
        for k in range(2):
            P = FrameProg(S, 1); y = P.inp(0); rand_path(P, y); out.append((f"randW1_{k}", P))
        return out

    return cases()
