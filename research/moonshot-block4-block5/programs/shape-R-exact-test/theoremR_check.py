"""theoremR_check.py -- o3-shapeR: computational check of THEOREM R (PROOF.md 6d).  Each input of a shape-R prefix is
typed by u3's frames other than O and N_l (one term each, else 'multi'): type 0 = no depth-1 frame; type R = some
depth-1 frame next to at most one depth-2 frame of u3; type M = every depth-1 frame next to both of its possible
depth-2 frames, both in u3.  Theorem R: (0, 0) -> dim eps(M4), dim eps(M5') <= 1; every other single-term pair ->
eps(M4) = 0.  Builds prefixes with run_fam.py's builder and checks every one.
usage: python3 -I theoremR_check.py FAMILY SEED NPREFIX"""
import sys, os, collections
sys.argv = [sys.argv[0], sys.argv[1], "0", sys.argv[2], "/dev/null"] + sys.argv[3:]
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_fam.py")).read().split("stats = collections.Counter()")[0])
NP = int(sys.argv[5]) if len(sys.argv) > 5 else 100
st = collections.Counter()
def itype(u3, l, Nl):
    frames = collections.Counter(); rep = {}
    for (ll, M, q) in u3:
        if ll != l: continue
        f = fk(M)
        if f in (W.L0, Nl): continue
        frames[f] += 1; rep[f] = M
    if any(v > 1 for v in frames.values()):
        return "multi"
    d2 = [f for f in frames if abs(W.dep.get(f, W.st(rep[f])[0]) - 2) < 1e-9 or W.st(rep[f])[0] == 2]
    d1 = [f for f in frames if W.st(rep[f])[0] == 1]
    if not d1:
        return "0"
    for f in d1:
        nb = set(W.nbrs(rep[f])[1])
        if sum(1 for h in d2 if h in nb) <= 1:
            return "R"
    return "M"
for k in range(NP):
    b = build()
    if b is None or isinstance(b, str):
        continue
    gens, info = b
    u3 = gens[4]
    ty = tuple(itype(u3, l, fk(gens[2 + l][0][1])) for l in (0, 1))
    res, n1, n2 = test_prefix(W, gens)
    d4 = max((x["d4"][0] for x in res), default=0); d5 = max((x["d5"][0] for x in res), default=0)
    key = "multi" if "multi" in ty else ("00" if ty == ("0", "0") else "other")
    ok = (d4 <= 1 and d5 <= 1) if key == "00" else (d4 == 0 if key == "other" else None)
    st[(fam, "types", tuple(sorted(ty)), "dim eps M4 max", d4, "M5' max", d5, {True: "AS THEOREM R", False: "VIOLATES", None: "outside"}[ok])] += 1
for kk, v in sorted(st.items(), key=lambda kv: str(kv[0])):
    print(kk, v)
