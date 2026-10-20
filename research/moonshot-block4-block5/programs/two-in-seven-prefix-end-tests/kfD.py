"""kfD.py -- o3-finish2 (14 Oct 2026): Criterion D of the augmentation lemma (widest end test, every end move).

eps(M_D1) (o3-finish; seeds = every row of K's frame, treated as sinks here) lies in the span of the outputs o(c) of
the columns c that touch a seed and are NOT dead (dead = coefficient forced to 0 in every solution of the non-seed
equations).  Seed-touching columns, up to the symmetries (B' for span-only components, Lambda = B' cap D1 B' D1^-1 for
the rest; both permute solutions and keep outputs):
  span class t   : b u_t, u_t with a term in K's frame; output sigma(u_t) for every b.  All b u_t are dead iff
                   (pure) u_t is dead in the span-only system AND every special translate c1^-1 u_t (one per
                   Lambda-class, from r g'^-1 = c1 D1 c2, r a row of the span-only component N_t, g' a generator term)
                   is dead in the D1 system.  Needs N_t complete.
  direct class   : D1 b4 u_s with K g^-1 = b3 D1 b4 for a term g of u_s; output sigma(D1 b4 u_s); one Lambda-orbit.
Criterion D: if for some direction v every class whose output is not parallel to v is dead, dim eps(M_D1) <= 1 and the
prefix is excluded for every end move E6, E7 and every Q in GL_2(Q[B']).
Deadness = the start column's own coordinate is forced to 0 by the equations of the processed rows of a BFS from it
(sound relaxation: fewer equations can only free it), checked at row checkpoints; linear algebra = fulltest.solve_dimM
(exact singleton peeling + dense rank mod BOTH primes 1048573, 2097133; 'dead' needs dead mod both).
usage: python3 kfD.py PREFIXES.jsonl OUT.jsonl [first] [step] [caprows] [only RUN:IT,...]"""
import sys, os, json, time, collections
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault('O3F_PRIMES', '1048573,2097133')
FT = os.path.join(HERE, 'fulltest_dir')
if not os.path.isdir(FT):
    FT = HERE   # fulltest.py and fulltestG.py sit in this folder
sys.path.insert(1, FT)
from kfG import *
import fulltest as FTM

CHECKS = (100, 400, 1500, 5000, 15000, 40000)

EXACT = os.environ.get('KFD_EXACT', '0') == '1'


def qrank(rowsd):
    """exact rank over Q of a sparse matrix given as list of dict(col -> Fraction)"""
    piv = {}          # pivot col -> row dict (normalised, pivot = leading entry in column order)
    r = 0
    for row in rowsd:
        row = {c: v for c, v in row.items() if v != 0}
        while row:
            c = min(row)
            if c in piv:
                f = row[c]
                for c2, v2 in piv[c].items():
                    nv = row.get(c2, 0) - f * v2
                    if nv == 0:
                        row.pop(c2, None)
                    else:
                        row[c2] = nv
            else:
                inv = 1 / row[c]
                piv[c] = {c2: v2 * inv for c2, v2 in row.items()}
                r += 1
                break
    return r


def exact_dimT(nrows, ent, tset):
    """as fulltest.solve_dimM, exactly over Q (single T column): 1 if the T coordinate is free, 0 if forced"""
    rows = collections.defaultdict(dict)
    for (r, c, co) in ent:
        nv = rows[r].get(c, Fraction(0)) + Fraction(co)
        if nv == 0:
            rows[r].pop(c, None)
        else:
            rows[r][c] = nv
    R = [d for d in rows.values() if d]
    (t,) = tuple(tset)
    rA = qrank([dict(d) for d in R])
    rX = qrank([{c: v for c, v in d.items() if c != t} for d in R])
    return 1 - rA + rX


def exact_parallel(comp, v):
    """image of ker(A) under O lies in span(v)?  <=>  rank([A; v2 o1 - v1 o2]) = rank(A), exactly over Q"""
    keep = set(comp["processed"])
    rows = collections.defaultdict(dict); w = {}
    for (rk, c, co) in comp["ent"]:
        co = Fraction(co)
        if rk[0] == 'seed':
            f = v[1] if rk[1] == 0 else -v[0]
            w[c] = w.get(c, Fraction(0)) + f * co
        elif rk in keep:
            rows[rk][c] = rows[rk].get(c, Fraction(0)) + co
    R = [dict(d) for d in rows.values()]
    rA = qrank([dict(d) for d in R])
    return qrank(R + [w]) == rA



def sigma(S, gen):
    s = [Fraction(0), Fraction(0)]
    for (l, gid), co in gen.items():
        if S.G.frm(gid) == S.Kfrm:
            s[l] += Fraction(co)
    return s


def dead_test(S, start, withE, caprows):
    """BFS from the start column; 'dead' if its coordinate is forced to 0 by the processed rows' equations.
    Returns (status, nproc, ncols, processed_rows, comp) with comp = kfG-style component when the BFS closed."""
    G = S.G
    cols = {}; ent = []; rows = {}; went = []
    q = collections.deque(); processed = []
    def addcol(kind, i, bmat):
        ck = (kind, i, G.intern(bmat))
        if ck in cols:
            return cols[ck]
        c = len(cols); cols[ck] = c
        for (rk, co) in S.col_entries(kind, i, bmat, withE):
            went.append((rk, c, co))
            if rk[0] == 'seed':
                continue
            if rk not in rows:
                rows[rk] = len(rows); q.append(rk)
            ent.append((rows[rk], c, co))
        return c
    c0 = addcol(*start)
    checks = [x for x in CHECKS if x < caprows] + [caprows]
    ci = 0
    while True:
        if not q or (ci < len(checks) and len(processed) >= checks[ci]):
            keep = set(rows[rk] for rk in processed)
            use = [(r, c, co) for (r, c, co) in ent if r in keep]
            d = exact_dimT(len(rows), use, {c0}) if EXACT else FTM.solve_dimM(len(rows), use, {c0})
            if d == 0:
                return "dead", len(processed), len(cols), list(processed), None
            if not q:
                comp = {"cols": cols, "rows": rows, "ent": went, "processed": processed, "capped": False}
                return ("alive" if d is not None else "undecided"), len(processed), len(cols), list(processed), comp
            ci += 1
            if ci >= len(checks):
                return "undecided", len(processed), len(cols), None, None
        if len(G.MATS) > MATCAP:
            return "undecided-mem", len(processed), len(cols), None, None
        rk = q.popleft(); processed.append(rk)
        (l, gid) = rk
        Mg = G.MATS[gid]
        for (i, gterm) in S.index.get((G.frm(gid), l), ()):
            if i >= 2:
                addcol('S', i, G.mul(Mg, np.conj(G.MATS[gterm]).T))
        if withE:
            hm = G.mul(S.D1i, Mg)
            for (i, gterm) in S.index.get((G.frame_key(hm), l), ()):
                if i >= 2:
                    addcol('E', i, G.mul(hm, np.conj(G.MATS[gterm]).T))


def parallel(a, b):
    return a[0] * b[1] - a[1] * b[0] == 0


def test(rec, caprows=40000):
    t0 = time.time()
    S = Sys(rec); G = S.G
    gens = S.gens
    out = {"run": rec.get("run", "").replace(".jsonl", ""), "it": rec.get("it")}
    classes = []     # dict(kind, label, out, start(s))
    for t in range(2, len(gens)):
        v = sigma(S, gens[t])
        if any(G.frm(gid) == S.Kfrm for (l, gid) in gens[t]):
            classes.append({"kind": "span", "t": t, "out": v})
    seen = set()
    for s in range(2, len(gens)):
        for (l, gid) in gens[s]:
            d = decomp1(G.mul(G.K, np.conj(G.MATS[gid]).T))
            if d is None:
                continue
            key = (s, G.frame_key(G.mul(G.D1, d[1])))
            if key in seen:
                continue
            seen.add(key)
            v = sigma(S, G.lmul(G.mul(G.D1, d[1]), gens[s]))
            classes.append({"kind": "dir", "s": s, "out": v, "b4": d[1]})
    outs = [c["out"] for c in classes if c["out"] != [0, 0]]
    out["outputs"] = [[c["kind"], c.get("t", c.get("s")), [str(x) for x in c["out"]]] for c in classes]
    # candidate directions: none (all dead) + each distinct output direction
    dirs = [None]
    for v in outs:
        if not any(w is not None and parallel(v, w) for w in dirs):
            dirs.append(v)
    def bad(v):
        return [k for k, c in enumerate(classes) if c["out"] != [0, 0] and (v is None or not parallel(c["out"], v))]
    # order: fewest bad classes, then fewest span classes among them
    dirs.sort(key=lambda v: (len(bad(v)), sum(classes[k]["kind"] == "span" for k in bad(v))))
    if len(dirs) == 1 or len(bad(dirs[0])) == 0:
        out.update(verdict="excluded", how="crit0", sec=round(time.time() - t0, 1))
        return out
    status = {}
    def comp_images(comp):
        if EXACT:
            return ("EXACT", comp)
        ims = []
        for pr in PRIMES:
            im = S.image(comp, pr)
            if im is None:
                return None
            ims.append([v for v in im if tuple(v) != (0, 0)])
        return ims
    def run_one(start, withE):
        nonlocal S
        if len(S.G.MATS) > MATCAP // 2:
            S = Sys(rec)      # memory: rebuild the Game between tests (starts are plain matrices)
        r = dead_test(S, start, withE, caprows)
        if os.environ.get('KFD_VERBOSE') == '1':
            print("   test", start[0], start[1], "withE", withE, "->", r[0], r[1], r[2], "mats", len(S.G.MATS),
                  round(time.time() - t0, 1), file=sys.stderr, flush=True)
        if r[0] == "dead":
            return ("dead", r[1:3], r[3])
        if r[0] == "alive":
            ims = comp_images(r[4])
            if ims is not None:
                return ("alive", r[1:3], ims)
        return (r[0], r[1:3], None)
    def class_status(k):
        """list of (status, size, data) over the class's symmetry classes of columns; stops at the first failure"""
        if k in status:
            return status[k]
        c = classes[k]
        if c["kind"] == "dir":
            st = run_one(('E', c["s"], c["b4"]), True)
            status[k] = [st if st[0] != "dead" else (st[0], st[1], None)]
            return status[k]
        t = c["t"]
        pure = run_one(('S', t, np.eye(16, dtype=complex)), False)
        res = [("pure-" + pure[0], pure[1], pure[2] if pure[0] == "alive" else None)]
        if pure[0] not in ("dead", "alive"):
            status[k] = res; return res
        R0 = pure[2] if pure[0] == "dead" else None
        if R0 is None:
            # alive and complete: its rows are the whole component
            R0 = None
        rowsR0 = R0 if R0 is not None else []
        if pure[0] == "alive":
            # pure component complete: translate starts from all its rows (the BFS closed)
            rr = dead_test(S, ('S', t, np.eye(16, dtype=complex)), False, caprows)
            rowsR0 = rr[3] or []
        starts = {}
        rowmats = [(l, S.G.MATS[gid]) for (l, gid) in rowsR0]
        for (l, Mr) in rowmats:
            for s2 in range(len(S.gens)):
                for (l2, g2) in S.gens[s2]:
                    if l2 != l:
                        continue
                    dd = decomp1(S.G.mul(Mr, np.conj(S.G.MATS[g2]).T))
                    if dd is None:
                        continue
                    bt = np.conj(dd[0]).T
                    starts.setdefault(S.G.frame_key(S.G.mul(S.D1i, bt)), bt)
        res.append(("translates", len(rowsR0), len(starts)))
        for key, bt in starts.items():
            st = run_one(('S', t, bt), True)
            res.append(st if st[0] != "dead" else ("dead", st[1], None))
            if st[0] not in ("dead", "alive"):
                break
        status[k] = res
        return res
    def ok_for(k, v):
        for st in class_status(k):
            s0 = st[0].replace("pure-", "")
            if s0 in ("dead", "translates"):
                continue
            if s0 != "alive" or st[2] is None:
                return False
            if v is None:
                return False
            if EXACT:
                if not exact_parallel(st[2][1], v):
                    return False
                continue
            for pi, ims in enumerate(st[2]):
                pr = PRIMES[pi]
                vp = (int(v[0].numerator) * pow(int(v[0].denominator), pr - 2, pr) % pr,
                      int(v[1].numerator) * pow(int(v[1].denominator), pr - 2, pr) % pr)
                for w in ims:
                    if (vp[0] * w[1] - vp[1] * w[0]) % pr:
                        return False
        return True
    tried = []
    for v in dirs:
        ok = True
        for k in bad(v):
            if not ok_for(k, v):
                ok = False; break
        tried.append([None if v is None else [str(x) for x in v], ok])
        if ok:
            out.update(verdict="excluded", how="dead", direction=tried[-1][0])
            break
    else:
        out["verdict"] = "open"
    out["tried"] = tried
    out["class_status"] = {str(k): [[x[0], list(x[1]) if isinstance(x[1], (list, tuple)) else x[1]] for x in v]
                           for k, v in status.items()}
    out["sec"] = round(time.time() - t0, 1)
    return out


if __name__ == "__main__":
    inp, outp = sys.argv[1], sys.argv[2]
    first = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    step = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    caprows = int(sys.argv[5]) if len(sys.argv) > 5 else 40000
    only = None
    if len(sys.argv) > 6 and sys.argv[6]:
        only = {(x.split(':')[0], int(x.split(':')[1])) for x in sys.argv[6].split(',')}
    recs = [json.loads(l) for l in open(inp)]
    fo = open(outp, "a")
    for kk in range(first, len(recs), step):
        rec = recs[kk]
        key = (rec.get("run", "").replace(".jsonl", ""), rec.get("it"))
        if only is not None and key not in only:
            continue
        r = test(rec, caprows)
        r["k"] = kk; r["exact"] = EXACT
        fo.write(json.dumps(r, default=str) + "\n"); fo.flush()
        print(kk, r["run"], r["it"], r["verdict"], r.get("how"), r.get("direction"), r.get("class_status"),
              r["sec"], "s", flush=True)
