"""rc_run.py (rm2-exact): build a multi-target system, run rcdp (guillotine rank-cover bound) on chosen forms,
and verify every reported violation INDEPENDENTLY (partition of D_F x D_F, each piece's bound re-derived with a
second parity routine through field evaluation, sum < s_F).

usage: rc_run.py <case> [dmax] [maxforms] [order]
cases: see CASES below.  Output: one summary line + certificate counts.
"""
import sys, os, subprocess, struct, time, random, itertools
import numpy as np
from xbase import Field, Forms, System, canon_targets, identity

HERE = os.path.dirname(os.path.abspath(__file__))


def frame_dirs(Fd, beta):
    """coordinates of beta*f_j (the frame of target mult(beta^2))"""
    return [Fd.coord(Fd.mul(beta, f)) for f in Fd.basis]


def sqrt_el(Fd, g):
    return Fd.exp[(Fd.log[g] * ((Fd.n + 1) // 2)) % Fd.n]


def build(case):
    """returns (m, targets as matrices, counts dict, label, expected)"""
    kind = case.split(":")[0]
    if kind == "canon":            # canon:m:k1,k2,...  one move per direction, canonical targets
        _, m, ks = case.split(":"); m = int(m); ks = [int(x) for x in ks.split(",")]
        Fd = Field(m); return m, canon_targets(Fd, ks), None, "open"
    if kind == "slice":            # slice:m:W:r   r copies of every direction, all targets I
        _, m, W, r = case.split(":"); m, W, r = int(m), int(W), int(r)
        return m, [identity(m)] * W, {z: r for z in range(1, 1 << m)}, "known"
    if kind == "plain":            # plain:m:k1,...  plain schedule for canonical targets (exists)
        _, m, ks = case.split(":"); m = int(m); ks = [int(x) for x in ks.split(",")]
        Fd = Field(m); cnt = {}
        for k in ks:
            beta = Fd.exp[k % Fd.n]          # T_k = mult(alpha^(2k)) = mult(beta^2), beta = alpha^k
            for z in frame_dirs(Fd, beta): cnt[z] = cnt.get(z, 0) + 1
        return m, canon_targets(Fd, ks), cnt, "exists"
    if kind == "dirs":             # dirs:m:z1,z2,...:W  explicit direction multiset, targets I
        _, m, zs, W = case.split(":"); m, W = int(m), int(W); cnt = {}
        for z in zs.split(","): cnt[int(z)] = cnt.get(int(z), 0) + 1
        return m, [identity(m)] * W, cnt, "exists"
    if kind == "gen":              # gen:m:z1,z2,...:t1,t2,...  rm-theory t_small target indices (any targets ~ I)
        _, m, zs, ts = case.split(":"); m = int(m); cnt = {}
        for z in zs.split(","): cnt[int(z)] = cnt.get(int(z), 0) + 1
        allT = targets_list(m)
        return m, [allT[int(t)] for t in ts.split(",")], cnt, "open"
    raise ValueError(case)


def targets_list(m):
    """rm-theory t_small.targets_of order: bits over the pairs (i<=j), nonzero diagonal, full rank over F_2"""
    pairs = [(i, j) for i in range(m) for j in range(i, m)]; out = []
    for bits in range(1 << len(pairs)):
        S = [[0] * m for _ in range(m)]
        for idx, (i, j) in enumerate(pairs):
            if (bits >> idx) & 1: S[i][j] = S[j][i] = 1
        if not any(S[i][i] for i in range(m)): continue
        rows = [sum(S[i][j] << j for j in range(m)) for i in range(m)]; basis = []
        for v in rows:
            for b in basis: v = min(v, v ^ b)
            if v: basis.append(v)
        if len(basis) != m: continue
        out.append(S)
    return out


_FRAME = {}
def frame_of(m, T):
    """columns v_1..v_m with sum v_j v_j^T = T (brute force, cached)"""
    key = (m, tuple(map(tuple, T)))
    if key in _FRAME: return _FRAME[key]
    for vs in itertools.combinations_with_replacement(range(1, 1 << m), m):
        S = [[0] * m for _ in range(m)]
        for v in vs:
            for j in range(m):
                for k in range(m): S[j][k] ^= ((v >> j) & 1) & ((v >> k) & 1)
        if S == [list(r) for r in T]:
            _FRAME[key] = list(vs); return list(vs)
    raise ValueError("no frame")


def second_parity(Fd, fo, q, T):
    """<Q,T> recomputed as sum_j Q(g e_j) with T = g g^T, g = mult(sqrt) for canonical / identity targets."""
    m = fo.m
    cols = frame_of(m, T)
    def Qv(z):
        v = 0
        for p, (j, k) in enumerate(fo.PAIRS):
            if (q >> p) & 1: v ^= ((z >> j) & 1) & ((z >> k) & 1)
        return v
    return sum(Qv(c) for c in cols) % 2


def verify(Fd, fo, sy, targets, F, pieces):
    m = fo.m; nd = len(fo.dirs); cnt = sy.counts
    w = lambda msk: sum(cnt.get(z, 0) for z in range(1, nd + 1) if msk >> (z - 1) & 1)
    def supp(q):
        msk = 0
        for z in range(1, nd + 1):
            v = 0
            for p, (j, k) in enumerate(fo.PAIRS):
                if (q >> p) & 1: v ^= ((z >> j) & 1) & ((z >> k) & 1)
            if v and cnt.get(z, 0): msk |= 1 << (z - 1)
        return msk
    def slack(q):
        return w(supp(q)) - sum(second_parity(Fd, fo, q, T) for T in targets)
    DF = supp(F); sF = slack(F)
    cells = set(); total = 0
    for (R, C, b, how, G) in pieces:
        assert R & ~DF == 0 and C & ~DF == 0 and R and C
        for x in range(nd):
            if R >> x & 1:
                for y in range(nd):
                    if C >> y & 1:
                        assert (x, y) not in cells, "overlap"; cells.add((x, y))
        if how == 0: assert b == min(w(R), w(C)), "dim bound"
        else:
            assert (R | C) & ~supp(G) == 0, "piece not inside D_G"
            assert b == slack(G), "slack mismatch"
        total += b
    dl = [x for x in range(nd) if DF >> x & 1]
    assert len(cells) == len(dl) ** 2, "not a partition"
    assert total < sF, "no violation"
    return total, sF


def main():
    case = sys.argv[1]; dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    maxforms = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 9
    order = sys.argv[4] if len(sys.argv) > 4 else "srank"
    m, targets, counts, expected = build(case)
    Fd = Field(m); fo = Forms(m); sy = System(fo, targets, counts)
    nd = len(fo.dirs)
    viol = sy.violations()
    if len(viol):
        print("CASE %s expected=%s : DELETION VIOLATION at %d forms (excluded by deletion)" % (case, expected, len(viol)))
        return
    # masks restricted to moves present (directions with count > 0)
    present = sum(1 << i for i, z in enumerate(fo.dirs) if sy.counts.get(z, 0) > 0)
    masks = (fo.mask & present).astype(np.int64)
    ndirs_in = np.array([bin(int(x)).count("1") for x in masks])
    cand = [F for F in range(1, fo.nf) if 0 < ndirs_in[F] <= dmax and sy.s[F] > 0]
    if order == "srank":   # most demanding first: large s_F relative to support
        cand.sort(key=lambda F: (-(sy.s[F] / max(1, sy.d[F])), ndirs_in[F]))
    elif order == "random":
        random.Random(1).shuffle(cand)
    cand = cand[:maxforms]
    inp = os.path.join(HERE, "tmp_rc_%d.bin" % os.getpid())
    with open(inp, "wb") as f:
        f.write(struct.pack("<5i", m, nd, fo.nf, len(cand), dmax))
        f.write(np.array([sy.counts.get(z, 0) for z in fo.dirs], dtype="<i4").tobytes())
        f.write(masks.astype("<u4").tobytes())
        f.write(sy.s.astype("<i4").tobytes())
        f.write(np.array(cand, dtype="<i4").tobytes())
    t0 = time.time()
    out = subprocess.run([os.path.join(HERE, "rcdp"), inp], capture_output=True, text=True).stdout
    os.remove(inp)
    lines = out.splitlines(); nv = 0; nok = 0; tight_margin = {}
    i = 0; first = None
    while i < len(lines):
        L = lines[i]
        if L.startswith("F ") and "bound" in L:
            p = L.split(); s = int(p[5]); bnd = int(p[7])
            tight_margin[bnd - s] = tight_margin.get(bnd - s, 0) + 1
        if L.startswith("CERT"):
            F = int(L.split()[1]); pieces = []; i += 1
            while lines[i] != "END":
                p = lines[i].split(); pieces.append((int(p[1]), int(p[2]), int(p[3]), int(p[4]), int(p[5]))); i += 1
            tot, sF = verify(Fd, fo, sy, targets, F, pieces); nv += 1
            if first is None: first = (F, tot, sF, len(pieces))
        i += 1
    print("CASE %s expected=%s W=%d M=%d : tested %d forms (dmax %d) in %.1fs; VERIFIED violations %d%s; margin hist %s"
          % (case, expected, sy.W, sum(sy.counts.values()), len(cand), dmax, time.time() - t0, nv,
             (" first F=%d bound %d < s_F %d (%d pieces)" % first) if first else "",
             sorted(tight_margin.items())[:6]))


if __name__ == "__main__":
    main()
