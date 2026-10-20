"""tile_alt.py (rm2-phase, restart session 15 Oct 2026): rm2-exact's rank-cover tiling test (rcdp.c, Theorem 10
with a target per array, rm-theory 8.5) run on per-array targets that may be ALTERNATING (Lemma A of ../PROOF.md:
phase letters at the two ends allow every nonsingular symmetric target).

Why the test is valid for alternating targets: rm-theory 8.5 uses the targets only through p_l(Q) = <Q, T_l>, the
sign the program owes array l at the form Q; nothing there needs T_l = h h^T.

Input: a survivor file with lines  "z1,...,zM|t1,...,tW"  (directions as ints, bit j = coordinate j; target indices in
t_small_alt.py's order with ALT=1: the non-alternating nonsingular symmetric matrices first, then the alternating ones).
Steps:
  1. canonical form of each (direction multiset, target multiset) under GL(m,2) (z -> A z, T -> A T A^T), dedupe;
  2. per class: exact deletion check, then rcdp on every form F with s_F > 0 (|D_F| <= dmax), and every reported
     violation re-verified independently (the partition of D_F x D_F is checked cell by cell and every piece's bound
     re-derived with a second parity routine: <Q, T> = sum_j Q(v_j) for a decomposition T = sum_j v_j v_j^T found by
     brute force; alternating T need m + 1 vectors).
Usage: python3 tile_alt.py m survivors.txt [dmax] [maxclasses] [selftest]
"""
import sys, os, itertools, struct, subprocess, time
import numpy as np
os.environ.setdefault("ALT", "1")
from xbase import Forms, System
from t_free import forms as tf_forms, gl_elements, apply as gl_apply
import t_small_alt

HERE = os.path.dirname(os.path.abspath(__file__))


def target_matrices(m):
    t_small_alt.ALT = True
    pairs, mon, coeffs, vals = tf_forms(m)
    Tm = t_small_alt.targets_of(m, pairs)
    mats = []
    for row in Tm:
        S = [[0] * m for _ in range(m)]
        for idx, (i, j) in enumerate(pairs):
            S[i][j] = S[j][i] = int(row[idx])
        mats.append(S)
    return mats, t_small_alt.NONALT


def mat_key(S):
    m = len(S); return sum(S[i][j] << (i * m + j) for i in range(m) for j in range(m))


def congr(A_cols, S):
    """A S A^T with A given by its columns (ints)"""
    m = len(S)
    A = [[(A_cols[j] >> i) & 1 for j in range(m)] for i in range(m)]
    AS = [[sum(A[i][k] & S[k][j] for k in range(m)) & 1 for j in range(m)] for i in range(m)]
    return [[sum(AS[i][k] & A[j][k] for k in range(m)) & 1 for j in range(m)] for i in range(m)]


_DEC = {}
def decomposition(S):
    """vectors v_1..v_k (k minimal, k <= m + 2) with sum v v^T = S"""
    m = len(S); key = mat_key(S)
    if key in _DEC: return _DEC[key]
    for k in range(1, m + 3):
        for vs in itertools.combinations(range(1, 1 << m), k):
            T = [[0] * m for _ in range(m)]
            for v in vs:
                for i in range(m):
                    if (v >> i) & 1:
                        for j in range(m):
                            T[i][j] ^= (v >> j) & 1
            if T == S:
                _DEC[key] = list(vs); return list(vs)
    raise ValueError("no decomposition")


def qval(fo, q, z):
    v = 0
    for p, (j, k) in enumerate(fo.PAIRS):
        if (q >> p) & 1: v ^= ((z >> j) & 1) & ((z >> k) & 1)
    return v


def second_parity(fo, q, S):
    return sum(qval(fo, q, v) for v in decomposition(S)) % 2


def verify(fo, counts, targets, F, pieces):
    nd = len(fo.dirs)
    w = lambda msk: sum(counts.get(z, 0) for z in range(1, nd + 1) if msk >> (z - 1) & 1)
    def supp(q):
        msk = 0
        for z in range(1, nd + 1):
            if counts.get(z, 0) and qval(fo, q, z): msk |= 1 << (z - 1)
        return msk
    def slack(q):
        return w(supp(q)) - sum(second_parity(fo, q, T) for T in targets)
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


def run_case(fo, counts, targets, dmax=16, tag="", maxforms=None):
    """returns (verdict, info): 'deletion', 'excluded', or 'passes'"""
    sy = System(fo, targets, counts)
    nd = len(fo.dirs)
    viol = sy.violations()
    if len(viol):
        return "deletion", int(viol[0])
    present = sum(1 << i for i, z in enumerate(fo.dirs) if counts.get(z, 0) > 0)
    masks = (fo.mask & present).astype(np.int64)
    nin = np.array([bin(int(x)).count("1") for x in masks])
    cand = [F for F in range(1, fo.nf) if 0 < nin[F] <= dmax and sy.s[F] > 0]
    if maxforms is not None:      # quick pass: most demanding forms first (rm2-exact's 'srank' order)
        cand.sort(key=lambda F: (-(sy.s[F] / max(1, sy.d[F])), nin[F])); cand = cand[:maxforms]
    inp = os.path.join(HERE, "tmp_ta_%d.bin" % os.getpid())
    with open(inp, "wb") as f:
        f.write(struct.pack("<5i", fo.m, nd, fo.nf, len(cand), dmax))
        f.write(np.array([counts.get(z, 0) for z in fo.dirs], dtype="<i4").tobytes())
        f.write(masks.astype("<u4").tobytes())
        f.write(sy.s.astype("<i4").tobytes())
        f.write(np.array(cand, dtype="<i4").tobytes())
    out = subprocess.run([os.path.join(HERE, "rcdp"), inp], capture_output=True, text=True).stdout
    os.remove(inp)
    lines = out.splitlines(); nv = 0; first = None; margins = {}
    i = 0
    while i < len(lines):
        L = lines[i]
        if L.startswith("F ") and "bound" in L:
            p = L.split(); s = int(p[5]); bnd = int(p[7]); margins[bnd - s] = margins.get(bnd - s, 0) + 1
        if L.startswith("CERT"):
            F = int(L.split()[1]); pieces = []; i += 1
            while lines[i] != "END":
                p = lines[i].split(); pieces.append(tuple(int(x) for x in p[1:6])); i += 1
            tot, sF = verify(fo, counts, targets, F, pieces); nv += 1
            if first is None: first = (F, tot, sF, len(pieces))
        i += 1
    if nv: return "excluded", (nv, first, len(cand))
    return "passes", (len(cand), sorted(margins.items())[:4])


def perm_table(Gs, mats, keyidx, m):
    """perm[g, i] = index of A_g S_i A_g^T (vectorised in chunks; spot-checked against congr)"""
    Mt = np.array(mats, dtype=np.int64)                                     # T x m x m
    wts = (1 << (np.arange(m)[:, None] * m + np.arange(m)[None, :])).astype(np.int64)
    lut = np.full(1 << (m * m), -1, dtype=np.int64)
    for k, i in keyidx.items(): lut[k] = i
    perm = np.zeros((len(Gs), len(mats)), dtype=np.int32)
    for c0 in range(0, len(Gs), 512):
        As = np.array([[[(h[j] >> i) & 1 for j in range(m)] for i in range(m)] for h in Gs[c0:c0 + 512]], dtype=np.int64)
        AS = np.einsum('gik,tkj->gtij', As, Mt) & 1
        ASA = np.einsum('gtik,gjk->gtij', AS, As) & 1
        keys = (ASA * wts[None, None]).sum(axis=(2, 3))
        perm[c0:c0 + 512] = lut[keys]
    assert (perm >= 0).all()
    for gi in (0, len(Gs) // 2, len(Gs) - 1):
        for i in (0, len(mats) // 3, len(mats) - 1):
            assert perm[gi, i] == keyidx[mat_key(congr(list(Gs[gi]), mats[i]))]
    return perm


def canon(G, dirs, tset, tkeys, img_cache):
    """canonical key of (direction multiset, target multiset) under the group G (list of column tuples)"""
    best = None
    for gi, h in enumerate(G):
        dd = tuple(sorted(gl_apply(h, z) for z in dirs))
        tt = tuple(sorted(img_cache[(gi, t)] for t in tset))
        key = (dd, tt)
        if best is None or key < best: best = key
    return best


def main():
    m = int(sys.argv[1]); path = sys.argv[2]
    dmax = int(sys.argv[3]) if len(sys.argv) > 3 else 16
    maxcl = int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 9
    t0 = time.time()
    mats, nonalt = target_matrices(m)
    keyidx = {mat_key(S): i for i, S in enumerate(mats)}
    fo = Forms(m)
    rows = []
    for L in open(path):
        L = L.strip()
        if not L or "|" not in L: continue
        a, b = L.split("|")
        rows.append((tuple(int(x) for x in a.split(",")), tuple(int(x) for x in b.split(","))))
    print(f"[tile_alt m={m} file={os.path.basename(path)}: {len(rows)} survivor lines; {len(mats)} targets "
          f"({nonalt} non-alternating)]", flush=True)
    # group: stabiliser of each direction multiset (computed per multiset)
    G_all = gl_elements(m)
    classes = {}
    by_dirs = {}
    for d, t in rows: by_dirs.setdefault(d, []).append(t)
    for d, tl in by_dirs.items():
        dset = tuple(sorted(d))
        Gs = [h for h in G_all if tuple(sorted(gl_apply(h, z) for z in dset)) == dset]
        # target index permutation under each stabiliser element
        perm = perm_table(Gs, mats, keyidx, m)
        T = np.array(tl, dtype=np.int64)                  # n x W
        W = T.shape[1]; base = len(mats)
        best = None
        for gi in range(len(Gs)):
            im = np.sort(perm[gi][T], axis=1)
            code = np.zeros(len(T), dtype=np.int64)
            for c in range(W): code = code * base + im[:, c]
            best = code if best is None else np.minimum(best, code)
        for n_, code in enumerate(best.tolist()):
            tt = []
            c = code
            for _ in range(W): tt.append(c % base); c //= base
            key = (dset, tuple(sorted(tt)))
            classes[key] = classes.get(key, 0) + 1
        print(f"  dirs {dset}: stabiliser order {len(Gs)}; {len(tl)} tuples -> {len(set(best.tolist()))} classes "
              f"({time.time()-t0:.1f} s)", flush=True)
    keys = sorted(classes)
    print(f"  total classes {len(keys)}", flush=True)
    tally = {}
    nalt_hist = {}
    SH = os.environ.get("SHARD", "")          # "i/n": only classes with index = i mod n
    si, sn = (int(SH.split("/")[0]), int(SH.split("/")[1])) if SH else (0, 1)
    for ci, (dset, tt) in enumerate(keys[:maxcl]):
        if ci % sn != si: continue
        counts = {}
        for z in dset: counts[z] = counts.get(z, 0) + 1
        targets = [mats[t] for t in tt]
        QUICK = int(os.environ.get("QUICK", "0"))
        v, info = run_case(fo, counts, targets, dmax, maxforms=QUICK) if QUICK else ("passes", None)
        if v == "passes": v, info = run_case(fo, counts, targets, dmax)
        na = sum(1 for t in tt if t >= nonalt)
        ntight = len(System(fo, targets, counts).tight())
        v = v + ("" if ntight else "-notight")
        tally[v] = tally.get(v, 0) + 1
        nalt_hist.setdefault(na, {}); nalt_hist[na][v] = nalt_hist[na].get(v, 0) + 1
        if not v.startswith("excluded") or ci < 5:
            print(f"  class {ci}: dirs {dset} targets {tt} ({na} alternating; {classes[(dset, tt)]} tuples): {v} {info}",
                  flush=True)
    print(f"  verdicts: {tally}; by number of alternating targets: {nalt_hist}")
    print(f"  done in {time.time()-t0:.1f} s")


def selftest():
    """calibration: plain programs (must pass), lone alternating target with m moves (must fail by deletion)."""
    m = 4; fo = Forms(m)
    mats, nonalt = target_matrices(m)
    alts = [S for S in mats[nonalt:]]
    I = [[int(i == j) for j in range(m)] for i in range(m)]
    out = []
    # (a) plain 4 arrays: two alternating targets (5 moves each) + two targets I (4 moves each): 18 moves; must pass
    for trial in range(3):
        A1, A2 = alts[trial], alts[(7 * trial + 5) % len(alts)]
        counts = {}
        for S in (A1, A2, I, I):
            for v in decomposition(S): counts[v] = counts.get(v, 0) + 1
        out.append(("plain alt,alt,I,I", run_case(fo, counts, [A1, A2, I, I])))
    # (b) plain 2 alternating targets sharing nothing: 10 moves; must pass
    A1, A2 = alts[0], alts[1]; counts = {}
    for S in (A1, A2):
        for v in decomposition(S): counts[v] = counts.get(v, 0) + 1
    out.append(("plain alt,alt", run_case(fo, counts, [A1, A2])))
    # (c) negative: one alternating target on its own 4 rank-ones of a frame (no decomposition with 4): must be refused
    counts = {1: 1, 2: 1, 4: 1, 8: 1}
    out.append(("alt on 4 moves (impossible)", run_case(fo, counts, [alts[0]])))
    # (d) negative: alternating target with m + 1 moves NOT summing to it: refused
    dec = decomposition(alts[0]); bad = dec[:-1] + [15 if dec[-1] != 15 else 14]
    counts = {}
    for v in bad: counts[v] = counts.get(v, 0) + 1
    out.append(("alt on a wrong 5-set", run_case(fo, counts, [alts[0]])))
    for name, res in out: print("  selftest", name, "->", res)


def one_case():
    """tile_alt.py case m z1,z2,...  t1,t2,...  [dmax]"""
    m = int(sys.argv[2]); fo = Forms(m); mats, nonalt = target_matrices(m)
    counts = {}
    for z in sys.argv[3].split(","): counts[int(z)] = counts.get(int(z), 0) + 1
    tt = [int(t) for t in sys.argv[4].split(",")]
    dmax = int(sys.argv[5]) if len(sys.argv) > 5 else 16
    print("  case dirs", sys.argv[3], "targets", tt, "(alternating:", sum(t >= nonalt for t in tt), ") ->",
          run_case(fo, counts, [mats[t] for t in tt], dmax))


if __name__ == "__main__":
    if sys.argv[1] == "selftest": selftest()
    elif sys.argv[1] == "case": one_case()
    else: main()
