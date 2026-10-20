"""vgen.py (v-r3-3in14): generic-m version of my tests (deletion, zero-block with local certificate, supports inside
U), run on the track's must-fail controls (block-4 3 in 11, block-3 3 in 8) and their plain parents, using the
track's own matrices and move coordinates (my code is coordinate-free).  Output: vcerts_ctl.json, out_vgen.txt."""
import json, itertools
import numpy as np
from vcore import max_clique
TRK = "../"
D = json.load(open(TRK + "certs_min.json"))
out = []; res = []
for inst in D:
    if not inst["tag"].startswith("ctl-b") and inst["tag"] not in ("ctl-U3", "ctl-plain1in5"): continue
    m = inst["m"]; mons = [(i, j) for i in range(m) for j in range(i, m)]; P = len(mons); NF = 1 << P
    lam = np.arange(NF, dtype=np.int64)
    par = np.zeros(NF, dtype=np.uint8)
    for b in range(P): par ^= ((lam >> b) & 1).astype(np.uint8)
    def vec(T): return sum(1 << t for t, (i, j) in enumerate(mons) if T[i][j] & 1)
    def zz(z): return vec([[((z >> i) & 1) & ((z >> j) & 1) for j in range(m)] for i in range(m)])
    Qm = np.stack([par[lam & zz(z)] for z in inst["moves"]]); Pt = np.stack([par[lam & vec(T)] for T in inst["Ts"]])
    M, W = len(inst["moves"]), len(inst["Ts"])
    cnt = Qm.sum(0).astype(int); r = Pt.sum(0).astype(int)
    Dm = (Qm.astype(np.int64) * (np.int64(1) << np.arange(M, dtype=np.int64))[:, None]).sum(0)
    Pm = (Pt.astype(np.int64) * (np.int64(1) << np.arange(W, dtype=np.int64))[:, None]).sum(0)
    ok = bool((cnt >= r).all())
    ts = sorted(set(int(Dm[q]) for q in np.nonzero((cnt == r) & (r >= 1))[0]))
    J = [[any(t >> i & 1 and t >> j & 1 for t in ts) for j in range(M)] for i in range(M)]
    live = [i for i in range(M) if J[i][i]]
    cl = max_clique(live, J); best = None
    for q in np.nonzero(cnt + r < 2 * len(cl))[0]:
        d = int(Dm[q]); Y = max_clique([i for i in live if d >> i & 1], J)
        inside = all(any(t >> i & 1 and t >> j & 1 and t & ~d == 0 for t in ts) for i in Y for j in Y)
        if cnt[q] + r[q] < 2 * len(Y) and inside and (best is None or cnt[q] < best[1]):
            best = (int(q), int(cnt[q]), int(r[q]), Y, d)
    verdict = "EXCLUDED" if (not ok or best or len(cl) > W) else "passes"
    out.append(f"{inst['tag']} idx {inst['idx']}: m={m} W={W} M={M} deletion {'ok' if ok else 'FAILS'} clique {len(cl)} "
               f"my verdict {verdict} (track: {inst['verdict']}) local cert {best and best[1:4]}")
    if (best and inst["expect"] == "fail") or inst["tag"] in ("ctl-U3", "ctl-plain1in5"):
        d = best[4] if inst["expect"] == "fail" else (1 << M) - 1; U = [i for i in range(M) if d >> i & 1]; sups = {}
        for q in range(1, NF):
            dm = int(Dm[q])
            if dm and dm & ~d == 0: sups[dm] = [int(Pm[q]) >> l & 1 for l in range(W)]
        res.append({"n": f"{inst['tag']}-{inst['idx']}", "U": U, "sups": [[[i for i in U if dm >> i & 1], p] for dm, p in sorted(sups.items())]})
json.dump([x for x in res if not x["n"].startswith(("ctl-U3", "ctl-plain"))], open("vcerts_ctl.json", "w")); json.dump([x for x in res if x["n"].startswith(("ctl-U3", "ctl-plain"))], open("vcerts_pass.json", "w")); open("out_vgen.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
