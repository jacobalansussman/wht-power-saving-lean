#!/usr/bin/env python3
"""twist.py (t-sandwich) -- interval lemma and the twisted-offset (interference-alignment) pair search on G_lag(m).

Frames: Lagrangians of F_2^(2m), built by ../lib/frames.py (v = x | z<<m, L_0 = {(0|z)}, L_F = F L_0,
subspace frame L_U = {(u | u+v): u in U, v in U^perp}).  d(L,L') = m - dim(L & L').
I(X,Y) = {Z : d(X,Z)+d(Z,Y) = d(X,Y)} (geodesic interval).
Port u (odd u, the hub's complex labels):  X_u = I(L_<u>, L_F) (source register X_u: <u> -> F, drop-free),
                                           Y_u = I(L_0, L_{u^perp}) (target register Y_u: 0 -> u^perp, drop-free).
Old knows:  S knows T  iff  X_T & Y_S nonempty  (= u_S.u_T = 0 on subspace frames).
Twisted pair (wider game, see NOTE): one total R carries x_S and x_T with relative symplectic twist phi.
  a0 in X_S, phi a0 in X_T            (R reads both sources at their own frames at time t1)
  a* in Y_S, phi a* in Y_T            (R delivers x_S to Y_S and x_T to Y_T with exact offset)
  mu in Y_S cap I(a*, eps_S), phi mu in X_T cap I(phi a0, L_F)   (Y_S cancels its x_T term with X_T itself)
  nu in X_S cap I(a0, L_F), phi nu in Y_T cap I(phi a*, eps_T)   (Y_T cancels its x_S term with X_S itself)
  copy cost of R: l = d(a0, a*).
usage: python3 twist.py M   (M = 2 or 3; prints and writes out/twist_mM.txt)
"""
import os, sys, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__))
TG = os.path.normpath(os.path.join(HERE, "../lib"))
sys.path.insert(0, TG)
import numpy as np
from frames import lag_graph, L_of, span, perp, popc, omega

m = int(sys.argv[1]) if len(sys.argv) > 1 else 2
out_lines = []
def out(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out_lines.append(s)

G = lag_graph(m); n = G["n"]; frames = G["frames"]; idx = G["idx"]; i0, iF = G["i0"], G["iF"]
masks = [sum(1 << v for v in L) for L in frames]
D = np.array([[m - (popc(masks[i] & masks[j]).bit_length() - 1) for j in range(n)] for i in range(n)], dtype=np.int8)
def interval(a, b):
    return np.nonzero(D[a, :] + D[:, b] == D[a, b])[0]
def nsub(d):
    return sum(1 for k in range(d + 1) for _ in [0] for __ in range(1)) and \
        sum(int(np.prod([(2 ** (d - i) - 1) for i in range(k)]) // max(1, int(np.prod([(2 ** (k - i) - 1) for i in range(k)])))) for k in range(d + 1))

# ---------------- 1. interval lemma ----------------
subs = set()
for k in range(m + 1):
    for vs in itertools.combinations(range(1, 1 << m), k):
        subs.add(span(vs))
subframes = {idx[L_of(U, m)]: U for U in subs}
Sub = set(interval(i0, iF).tolist())
assert Sub == set(subframes), "I(L_0,L_F) != subspace frames"
srcs = range(n) if m <= 3 else [i0, iF] + list(range(0, n, 97))
bad = 0; cnt = 0
for a in srcs:
    for b in range(n):
        I = interval(a, b); cnt += 1
        d = int(D[a, b])
        if len(I) != nsub(d):
            bad += 1
        # structure: Z in I  iff  Z = (Z&a) + (Z&b) (as subspaces)
        La, Lb = frames[a], frames[b]
        for z in range(n):
            Lz = frames[z]
            gen = span(list((Lz & La) | (Lz & Lb)))
            if (z in set(I.tolist())) != (gen == Lz):
                bad += 1
out("m=%d frames %d, subspace frames %d = I(L_0,L_F): yes; interval lemma |I(X,Y)| = #subspaces of F_2^d(X,Y) and"
    " I(X,Y) = {Z : Z = (Z&X)+(Z&Y)} on %d ordered pairs: %d failures" % (m, n, len(Sub), cnt, bad))
exc = {z: int(D[i0, z] + D[z, iF] - m) for z in range(n)}
hist = {}
for z in range(n): hist[exc[z]] = hist.get(exc[z], 0) + 1
out("excess d(L_0,Z)+d(Z,L_F)-m over all frames:", dict(sorted(hist.items())), "(0 exactly on the subspace frames)")

# ---------------- 2. Sp(2m,2) as permutations of the frames ----------------
N2 = 2 * m
def apply_lin(cols, v):
    r = 0; b = 0
    while v:
        if v & 1: r ^= cols[b]
        v >>= 1; b += 1
    return r
def sp_elements():
    """all symplectic maps as column tuples, via a symplectic basis e_i=(e_i|0), f_i=(0|e_i)."""
    basis = [1 << i for i in range(m)], [1 << (m + i) for i in range(m)]
    res = []
    def rec(k, used_span_vecs, imgs):
        if k == m:
            res.append(tuple(imgs)); return
        Wp = [v for v in range(1, 1 << N2) if all(omega(v, w, m) == 0 for w in used_span_vecs)]
        for e in Wp:
            if e in span(used_span_vecs): continue
            for f in Wp:
                if omega(e, f, m) == 1:
                    rec(k + 1, used_span_vecs + [e, f], imgs + [(e, f)])
    rec(0, [], [])
    return res
t = time.time()
els = sp_elements()
out("|Sp(%d,2)| enumerated: %d (%.1fs)" % (N2, len(els), time.time() - t))
# column images: col for (e_i|0) is imgs[i][0], for (0|e_i) is imgs[i][1]; vector bit b<m is x_b, bit m+b is z_b
NE = len(els)
tab = np.zeros((NE, 1 << N2), dtype=np.uint8)
cols = np.zeros((NE, N2), dtype=np.uint8)
for k, im in enumerate(els):
    for i in range(m):
        cols[k, i] = im[i][0]; cols[k, m + i] = im[i][1]
for v in range(1, 1 << N2):
    acc = np.zeros(NE, dtype=np.uint8)
    for b in range(N2):
        if (v >> b) & 1: acc ^= cols[:, b]
    tab[:, v] = acc
mask_index = {mk: i for i, mk in enumerate(masks)}
sorted_masks = np.array(sorted(masks), dtype=object)
P = np.zeros((NE, n), dtype=np.uint8 if n < 256 else np.int16)
for fi, L in enumerate(frames):
    mk = np.zeros(NE, dtype=object)
    imgs = [tab[:, v] for v in L]
    # frame image as a frozenset -> index; use sum of 1<<img as python ints (2m <= 6 -> 64 bits fits int64 except bit 63)
    acc = np.zeros(NE, dtype=np.uint64)
    for a in imgs: acc |= (np.uint64(1) << a.astype(np.uint64))
    lut = {np.uint64(mk_): i for i, mk_ in enumerate(masks)}
    P[:, fi] = np.array([lut[x] for x in acc.tolist()], dtype=np.int16) if NE < 5000 else \
        np.searchsorted(np.array(sorted(np.uint64(x) for x in masks), dtype=np.uint64), acc)
if NE >= 5000:   # map sorted positions back to frame indices
    order = np.argsort(np.array(masks, dtype=np.uint64))
    P = order[P].astype(P.dtype)
assert all(sorted(P[k].tolist()) == list(range(n)) for k in range(0, NE, max(1, NE // 50)))
# identity present, F present
ident = [k for k in range(NE) if (P[k] == np.arange(n)).all()]
assert len(ident) == 1
out("frame permutations built (%.1fs)" % (time.time() - t))

# ---------------- 3. ports and the twisted pair search ----------------
odd = [u for u in range(1, 1 << m) if popc(u) & 1]
def port(u):
    a = idx[L_of(span([u]), m)]; e = idx[L_of(perp(span([u]), m), m)]
    return a, e
inI = {}
def I_mask(a, b):
    key = (a, b)
    if key not in inI:
        z = np.zeros(n, dtype=bool); z[interval(a, b)] = True; inI[key] = z
    return inI[key]
leqF = np.array([I_mask(a, iF) for a in range(n)])            # leqF[a, z]: z in I(a, L_F)
def run_pair(uS, uT, phis=None, timed=False):
    """timed=False: the hub's generous convention (side information free of timing, u-floor hypothesis 5):
         mu in Y_S with phi mu in X_T;  nu in X_S with phi nu in Y_T.
       timed=True: chains are monotone in time; R reads at t1 < delivers at t0; each cancellation time t may be
         before t1, between, or after t0 (pre-cancellation is allowed: offsets are invariant)."""
    aS, eS = port(uS); aT, eT = port(uT)
    XS = np.nonzero(I_mask(aS, iF))[0]; YS = np.nonzero(I_mask(i0, eS))[0]
    XT = I_mask(aT, iF); YT = I_mask(i0, eT)
    Pp = P if phis is None else P[phis]
    best = None; found = []
    for a0 in XS:
        for ast in YS:
            l = int(D[a0, ast])
            base = XT[Pp[:, a0]] & YT[Pp[:, ast]]
            if not base.any(): continue
            okmu = np.zeros(len(Pp), dtype=bool)
            for mu in YS:
                c = XT[Pp[:, mu]]
                if timed:
                    pa0, pmu = Pp[:, a0], Pp[:, mu]
                    xt_up = leqF[pa0, pmu]          # X_T reaches phi mu after phi a0  (t > t1)
                    xt_dn = leqF[pmu, pa0]          # X_T at phi mu before phi a0     (t < t1)
                    ys_up = bool(I_mask(ast, eS)[mu])   # Y_S at mu after a*            (t > t0)
                    ys_dn = bool(I_mask(mu, eS)[ast])   # Y_S at mu before a*           (t < t0)
                    c = c & ((xt_dn & ys_dn) | (xt_up & ys_dn) | (xt_up & ys_up))
                okmu |= c
            oknu = np.zeros(len(Pp), dtype=bool)
            for nu in XS:
                c = YT[Pp[:, nu]]
                if timed:
                    pas, pnu = Pp[:, ast], Pp[:, nu]
                    xs_up = bool(leqF[a0, nu]); xs_dn = bool(leqF[nu, a0])
                    yt_up = I_mask_rows_eT(eT)[pas, pnu]; yt_dn = I_mask_rows_eT(eT)[pnu, pas]
                    c = c & ((xs_dn & yt_dn) | (xs_up & yt_dn) | (xs_up & yt_up))
                oknu |= c
            ok = base & okmu & oknu
            if ok.any():
                k = int(np.nonzero(ok)[0][0])
                found.append((l, int(a0), int(ast), k))
                if best is None or l < best: best = l
    return best, found, (aS, eS, aT, eT)
_rows = {}
def I_mask_rows_eT(eT):
    if eT not in _rows:
        _rows[eT] = np.array([I_mask(b, eT) for b in range(n)])
    return _rows[eT]

out("\nTwisted pair search, odd labels, m=%d (all phi in Sp(%d,2)):" % (m, N2))
idn = ident[0]
reps = [(uS, uT) for uS in odd for uT in odd if uS != uT]
if m == 3:   # coordinate permutations act transitively on each of these three classes
    reps = [(0b001, 0b010), (0b001, 0b111), (0b111, 0b001)]
for uS, uT in reps:
    if True:
        dot = popc(uS & uT) & 1
        old = "knows (u_S.u_T=0)" if dot == 0 else "does NOT know (u_S.u_T=1)"
        for timed in (False, True):
            best, found, pr = run_pair(uS, uT, timed=timed)
            bid, fid, _ = run_pair(uS, uT, phis=[idn], timed=timed)
            out("  u_S=%s u_T=%s  %s, %s: twisted min copy cost l = %s (a0,a* configs %d); phi = identity: %s"
                % (format(uS, "0%db" % m), format(uT, "0%db" % m), old, "timed" if timed else "generous",
                   best, len(found), bid))
os.makedirs(os.path.join(HERE, "out"), exist_ok=True)
open(os.path.join(HERE, "out", "twist_m%d.txt" % m), "w").write("\n".join(out_lines) + "\n")
