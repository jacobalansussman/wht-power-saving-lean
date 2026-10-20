"""fam_enum.py (rm2-phase, restart session): exhaustive S+D+C+Z pruning for ONE direction multiset with a large
stabiliser, using the stabiliser to fix one target (t_small_alt.py enumerates every target multiset; this does the
same job for one class, faster, and dedupes the survivors up to the stabiliser).

Completeness: let Gs be the stabiliser of the direction multiset in GL(m,2) (z -> A z, T -> A T A^T). Every target
multiset containing a target t is Gs-equivalent to one containing the representative r of t's Gs-orbit. So
enumerating, for each orbit representative r (only alternating orbits with ALTONLY=1), the multisets {r} + (W-1 others
from all targets in span) covers every class; the survivors are deduplicated by their canonical form under Gs.
The tests (S), (D), (C), (Z) are those of t_small_alt.py (leaf_tests imported unchanged); they are Gs-invariant.

usage: python3 fam_enum.py m W z1,z2,...,zM  [dumpfile]      env ALT=1 (alternating targets allowed), ALTONLY=1
"""
import sys, os, time
import numpy as np
from t_free import forms as tf_forms, gl_elements, apply as gl_apply
import t_small_alt
from t_small_alt import leaf_tests
from tile_alt import target_matrices, mat_key, perm_table

ALT = os.environ.get("ALT", "1") == "1"; ALTONLY = os.environ.get("ALTONLY", "0") == "1"
_FS = os.environ.get("FE_SHARD", "0/1"); FSI, FSN = int(_FS.split("/")[0]), int(_FS.split("/")[1])   # split the first free target


def main():
    m, W = int(sys.argv[1]), int(sys.argv[2]); ms = tuple(sorted(int(z) for z in sys.argv[3].split(",")))
    dump = sys.argv[4] if len(sys.argv) > 4 else None
    t0 = time.time()
    t_small_alt.ALT = ALT
    pairs, mon, coeffs, vals = tf_forms(m)
    Tm = t_small_alt.targets_of(m, pairs); nonalt = t_small_alt.NONALT
    mats, nonalt2 = target_matrices(m)
    if not ALT: mats = mats[:nonalt]
    assert nonalt2 == nonalt and len(mats) == len(Tm)
    for i in (0, nonalt - 1, len(Tm) - 1):          # same indexing as t_small_alt
        assert [mats[i][a][b] for (a, b) in pairs] == list(Tm[i])
    par = (coeffs @ Tm.T % 2).T.astype(np.int16)      # targets x forms
    cnt = np.zeros(vals.shape[0], dtype=np.int16)
    for z in ms: cnt += vals[:, z].astype(np.int16)
    Z = cnt == 0
    inV = np.where(~(par[:, Z].astype(bool)).any(axis=1))[0]
    Zv = vals[:, list(ms)].astype(bool)
    G = gl_elements(m)
    Gs = [h for h in G if tuple(sorted(gl_apply(h, z) for z in ms)) == ms]
    keyidx = {mat_key(S): i for i, S in enumerate(mats)}
    perm = perm_table(Gs, mats, keyidx, m)
    # orbits of Gs on targets in span
    inset = set(int(x) for x in inV); orbit_of = {}; reps = []
    for t in sorted(inset):
        if t in orbit_of: continue
        orb = set(int(x) for x in perm[:, t]); reps.append((t, len(orb)))
        for u in orb: orbit_of[u] = t
    print(f"[fam_enum m={m} W={W} dirs {ms} ALT={int(ALT)} ALTONLY={int(ALTONLY)}] stabiliser {len(Gs)}; targets in span "
          f"{len(inV)} of {len(Tm)}; orbits {[(r, s, 'alt' if r >= nonalt else 'na') for r, s in reps]} "
          f"({time.time()-t0:.1f} s)", flush=True)
    sub = par[inV].astype(np.int16); k = len(inV)
    base = len(mats)
    surv = {}; stats = {"leaves": 0, "excl_C": 0, "excl_Z": 0, "pass": 0}
    for r, osz in reps:
        if ALTONLY and r < nonalt: continue
        ri = int(np.where(inV == r)[0][0])
        s0 = sub[ri].copy()
        if not (cnt >= s0).all(): continue
        def dfs(tup, s, depth, a0):
            cand = s[None, :] + sub[a0:]
            ok = np.where((cnt[None, :] >= cand).all(axis=1))[0]
            for o in ok:
                j = a0 + int(o); t2 = tup + (j,)
                if depth == 0 and FSN > 1 and j % FSN != FSI: continue
                if depth + 1 == W - 1:
                    stats["leaves"] += 1
                    res = leaf_tests(Zv, cnt, cand[o], W)
                    if res[0]:
                        stats["pass"] += 1
                        tt = np.array([r] + [int(inV[x]) for x in t2], dtype=np.int64)
                        im = np.sort(perm[:, tt], axis=1)
                        code = np.zeros(len(Gs), dtype=np.int64)
                        for c in range(W): code = code * base + im[:, c]
                        cmin = int(code.min())
                        if cmin not in surv: surv[cmin] = (res[1], res[2])
                    elif res[3] is None: stats["excl_C"] += 1
                    else: stats["excl_Z"] += 1
                else:
                    dfs(t2, cand[o], depth + 1, j)
        if W == 1:
            stats["leaves"] += 1
            res = leaf_tests(Zv, cnt, s0, W)
            if res[0]: surv[r] = (res[1], res[2])
        else:
            dfs((), s0, 0, 0)
        print(f"  rep {r} ({'alt' if r >= nonalt else 'non-alt'}, orbit {osz}) done: {stats}; classes so far "
              f"{len(surv)} ({time.time()-t0:.1f} s)", flush=True)
    print(f"  SURVIVOR CLASSES (up to the stabiliser): {len(surv)}; leaves {stats['leaves']}, passing {stats['pass']}, "
          f"excluded by cliques {stats['excl_C']}, by the zero-block test {stats['excl_Z']}")
    if dump:
        with open(dump, "w") as fh:
            for code in sorted(surv):
                tt = []; c = code
                for _ in range(W): tt.append(c % base); c //= base
                fh.write(",".join(map(str, ms)) + "|" + ",".join(map(str, sorted(tt))) + "\n")
    nalt = {}
    for code in surv:
        c = code; na = 0
        for _ in range(W):
            na += (c % base) >= nonalt; c //= base
        nalt[na] = nalt.get(na, 0) + 1
    print(f"  classes by number of alternating targets: {dict(sorted(nalt.items()))}; with no tight form: "
          f"{sum(1 for v in surv.values() if v[0] == 0)}")
    print(f"  done in {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()
