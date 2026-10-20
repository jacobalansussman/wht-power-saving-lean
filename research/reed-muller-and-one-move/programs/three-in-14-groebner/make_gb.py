"""make_gb.py: write the Groebner inputs (gb.py) for every instance in certs.json into gbin/.
- Excluded instances: the subsystem of the forms in its zero-block certificate (tight witnesses + F).
- Passing controls: the same kind of subsystem for the form closest to failing (smallest |D_F| + r_F - 2|Y| over
  cliques Y inside D_F), so the must-pass controls go through the same pipeline.
- Plain controls: also substitute the explicit plain program (XOR chains, time order = move order) into every
  polynomial of the subsystem, exactly (Fractions): every one must vanish.  This checks the equation generator."""
import json, os, itertools
from fractions import Fraction
from f2q import Field, Config, popc, all_frames
from ztest import links, cliques
from gb import build, write

FIELDS = {m: Field(m) for m in (3, 4, 5)}
os.makedirs("gbin", exist_ok=True)
certs = json.load(open("certs.json"))


def closest(cfg):
    tight, wit = links(cfg)
    cl = [y for y in cliques(cfg, wit) if y]
    best = None; cmax = max(len(y) for y in cl) if cl else 0
    cl = sorted(cl, key=len, reverse=True)
    for Q in range(1, len(cfg.Dmask)):
        d = cfg.Dmask[Q]; r = popc(cfg.pmask[Q]); sd = popc(d)
        if best is not None and sd + r - 2 * cmax > best[0][0]: continue
        for y in cl:
            ym = sum(1 << i for i in y)
            if d & ym == ym:
                key = (sd + r - 2 * len(y), sd, -len(y))
                if best is None or key < best[0]: best = (key, Q, y)
    if best is None: return None, None
    _, Q, y = best
    forms = sorted({wit[(min(a, b), max(a, b))][1] for a in y for b in y}) + [Q]
    return forms, best[0]


def plain_solution(cfg, ks):
    """assign moves to arrays so each array's moves are a frame of T_k; return (a, b, H) exact."""
    F = cfg.F; W, M = cfg.W, cfg.M
    frs = [all_frames(F, k) for k in ks]
    target = sorted(cfg.moves)
    for ch in itertools.product(*frs):
        if sorted(z for fr in ch for z in fr) == target:
            owner = [None] * M; used = [False] * M
            for l, fr in enumerate(ch):
                for z in fr:
                    i = next(i for i in range(M) if cfg.moves[i] == z and not used[i])
                    used[i] = True; owner[i] = l
            a = [[Fraction(int(owner[i] == l)) for i in range(M)] for l in range(W)]
            H = [[Fraction(2) if (i < j and owner[i] == owner[j]) else Fraction(0) for j in range(M)] for i in range(M)]
            return a, a, H
    return None


def inv_exact(Mx):
    n = len(Mx); A = [row[:] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(Mx)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[p] = A[p], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


def substitute(S, cfg, sol):
    a, b, H = sol
    val = {}
    for name in S.vars:
        if name[0] == "a": l, i = map(int, name[1:].split("_")); val[name] = a[l][i]
        elif name[0] == "b": l, i = map(int, name[1:].split("_")); val[name] = b[l][i]
        elif name[0] == "g": i, j = map(int, name[1:].split("_")); val[name] = H[i][j]
    # n variables: N_D = (I + H_D)^-1
    tags = {}
    for name in S.vars:
        if name[0] == "n":
            tag, s, t = name[1:].split("_"); tags.setdefault(tag, set()).update([int(s), int(t)])
    for tag, mv in tags.items():
        D = sorted(mv); Mx = [[Fraction(int(i == j)) + (H[i][j] if i != j else 0) for j in D] for i in D]
        N = inv_exact(Mx)
        for x, s in enumerate(D):
            for y, t in enumerate(D): val[f"n{tag}_{s}_{t}"] = N[x][y]
    bad = 0
    for p in S.polys:
        tot = Fraction(0)
        for mono, c in p.items():
            term = Fraction(c)
            for vi in mono: term *= val[S.vars[vi]]
            tot += term
        bad += (tot != 0)
    return bad


manifest = []
for c in certs:
    F = FIELDS[c["m"]]; cfg = Config(F, c["Ts"], c["moves"])
    if c["cert"] and c["cert"]["kind"] == "Z":
        forms = sorted({p[2] for p in c["cert"]["pairs"]}) + [c["cert"]["F"]]; note = "zero-block certificate"
    else:
        forms, key = closest(cfg); note = f"closest-to-failing form, slack {key}"
    S = build(cfg, forms)
    base = f"gbin/i{c['idx']:03d}"
    write(S, base)
    sub = None
    if c["tag"].endswith("-plain") or c["tag"] == "ctl-plain1in5":
        sol = plain_solution(cfg, c["ks"])
        sub = "no frame partition found" if sol is None else substitute(S, cfg, sol)
    manifest.append(dict(idx=c["idx"], tag=c["tag"], expect=c["expect"], z_verdict=c["verdict"], base=base,
                         nvars=len(S.vars), npolys=len(S.polys), nforms=len(set(forms)), note=note,
                         plain_substitution_nonzero=sub))
    print(manifest[-1], flush=True)
json.dump(manifest, open("gbin/manifest.json", "w"), indent=0)
