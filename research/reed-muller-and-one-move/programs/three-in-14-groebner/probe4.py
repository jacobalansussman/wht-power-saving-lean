"""probe4.py WHICH [maxU]: BY-PRODUCT, outside the brief. Does the local exact algebra of gb2.py see past the
zero-block test on the smallest test-blind sharing shapes next to plain (four arrays)?
WHICH = b4 (my block-4 4-in-15 survivor, out_z4b4.txt) or b5:<line> (line of v-rm2-solve's survivors_b5w4.txt).
For every move set U = D_Q of a form with 4 <= |D_Q| <= maxU (and the unions of two such supports up to maxU):
the time-ordered system with all forms supported in U (gb2.system), Singular std over F_32003 with a degree bound
of 12 for the first two orders (a 1 found below the bound is a proof of 1 in the ideal mod p; a miss proves
nothing), then full std (degBound = 0) for every order.  A U is reported when
EVERY relative order (up to copies and reversal) is UNIT mod p.  Mod p only: anything reported needs Singular over Q
and an exact lift check before it counts.
Env: TARGETED=1 (targeted candidates), MAXC=n (16 Oct, after launch of the uncapped runs: try only the first n)."""
import sys, re, json, os, subprocess, random, tempfile, itertools
from f2q import Field, Config, popc
from gb2 import supports_in, orders, Ring, system, msolve_text, pstr

which = sys.argv[1]; maxU = int(sys.argv[2]) if len(sys.argv) > 2 else 6
if which in ("b4", "b4plain"):
    F = Field(4)
    line = [l for l in open("out_z4b4.txt") if l.startswith("SURVIVES")][0]
    ks = [int(v) for v in re.search(r"ks=\(([0-9, ]+)\)", line).group(1).split(",")]
    moves = [int(v) for v in re.search(r"moves=\[([0-9, ]+)\]", line).group(1).split(",")]
    if which == "b4plain":  # calibration: restore the deleted copy of direction 8 (a genuine plain 4 in 16)
        moves = sorted(moves + [8])
elif which.startswith("w3:"):  # calibration: an excluded 3-in-14 instance from certs_min.json
    c = json.load(open("certs_min.json"))[int(which.split(":")[1])]
    F = Field(c["m"]); ks = c["ks"]; moves = c["moves"]
else:
    F = Field(5); k = int(which.split(":")[1])
    line = [l for l in open("survivors_b5w4.txt") if l.startswith("SURVIVES")][k]
    ks = [int(v) for v in re.search(r"\(([0-9, ]+)\)", line).group(1).split(",")]
    moves = [int(v) for v in re.search(r"dirs \[([0-9, ]+)\]", line).group(1).split(",")]
cfg = Config(F, [F.target(k) for k in ks], moves)
print("instance", which, ks, moves, "span/deletion", cfg.span_deletion(), flush=True)
tmp = tempfile.mkdtemp(prefix="p4_", dir=".")


def sing_trunc(polysets, deg=12, char=32003, tl=60):
    """Singular std with the degree bound degBound over F_p: every element it returns lies in the ideal, so a 1 proves
    UNIT mod p; anything else means only 'not found below the degree bound' (reported as NOT-FOUND)."""
    NV = max(max(v for p in ps for m in p for v in m) for ps in polysets) + 1
    txt = [f"ring r = {char}, (x(1..{NV})), dp;", f"degBound = {deg};", "ideal I; ideal G;"]
    for k, ps in enumerate(polysets):
        txt.append("I = " + ",\n".join(pstr(p, "sing") for p in ps) + ";")
        txt.append(f'G = std(I); if (G[1] == 1) {{ "S {k} UNIT"; }} else {{ "S {k} NOT-FOUND"; }}')
    txt.append("quit;")
    fi = os.path.join(tmp, "t.sing"); open(fi, "w").write("\n".join(txt) + "\n")
    try:
        out = subprocess.run(["Singular", "-q", fi], capture_output=True, text=True, timeout=tl).stdout
    except subprocess.TimeoutExpired:
        out = ""
    return ["UNIT" if f"S {k} UNIT" in out else "NOT-FOUND" for k in range(len(polysets))]


def ms(polys):
    fi = os.path.join(tmp, "m.ms"); fo = os.path.join(tmp, "m.out")
    open(fi, "w").write(msolve_text(polys, 65521))
    if os.path.exists(fo): os.remove(fo)
    subprocess.run(["msolve", "-f", fi, "-o", fo, "-g", "2"], capture_output=True, timeout=600)
    s = open(fo).read() if os.path.exists(fo) else ""
    return "UNIT" if re.search(r"length of basis:\s+1 element", s) and re.search(r"^\[1\]:", s, re.M) else \
           ("NOT-UNIT" if "#length of basis" in s else "ERROR")


cands = set()
if os.environ.get("TARGETED"):
    # targeted: U = D_F for the forms F closest to failing the zero-block test (smallest |D_F| + r_F - 2|Y| over
    # cliques Y inside D_F), plus the tight supports linking Y, |U| <= maxU; ordered by that slack, then size
    from ztest import links, cliques
    tight, wit = links(cfg); cl = [y for y in cliques(cfg, wit) if len(y) >= 2]
    sc = {}
    for Q in range(1, len(cfg.Dmask)):
        d = cfg.Dmask[Q]; r = popc(cfg.pmask[Q]); sd = popc(d)
        if sd > maxU: continue
        for y in cl:
            ym = sum(1 << i for i in y)
            if d & ym != ym: continue
            u = d
            for a in y:
                for b in y: u |= wit[(min(a, b), max(a, b))][2]
            if popc(u) > maxU: continue
            key = (sd + r - 2 * len(y), popc(u))
            if u not in sc or key < sc[u]: sc[u] = key
    cands = sorted(sc, key=lambda u: (sc[u], u))
    print("targeted: slack of the first candidates", [sc[u] for u in cands[:10]], flush=True)
else:
    small = sorted({cfg.Dmask[Q] for Q in range(1, len(cfg.Dmask)) if 4 <= popc(cfg.Dmask[Q]) <= maxU})
    cands |= set(small)
    for a, b in itertools.combinations(small, 2):
        u = a | b
        if popc(u) <= maxU: cands.add(u)
    cands = sorted(cands, key=lambda u: (popc(u), u))
maxc = int(os.environ.get("MAXC", "0"))  # 16 Oct: optional cap, the first MAXC candidates only (closest first)
if maxc: print("capped to the first", maxc, "of", len(cands), "candidates", flush=True); cands = cands[:maxc]
print("candidate move sets", len(cands), flush=True)
found = []
for ci, u in enumerate(cands):
    if ci % 10 == 0: print(f"progress: tried {ci} of {len(cands)} candidates", flush=True)
    U = [i for i in range(cfg.M) if u >> i & 1]
    sup = supports_in(cfg, U); R = Ring(cfg.W, U)
    ords = list(orders(cfg, U))
    first = [ords[0], random.Random(u).choice(ords)]
    if any(v != "UNIT" for v in sing_trunc([system(cfg, R, o, sup) for o in first], tl=20)): continue
    allunit = True
    for b in range(0, len(ords), 20):
        if any(v != "UNIT" for v in sing_trunc([system(cfg, R, o, sup) for o in ords[b:b + 20]], deg=0, tl=600)):
            allunit = False; break
    print(f"U={U} dirs={[moves[i] for i in U]} supports={len(sup)} orders={len(ords)} first-two UNIT, all orders UNIT mod p: {allunit}", flush=True)
    if allunit:
        found.append(U); break
print("RESULT", which, "excluding move set mod p:", found[0] if found else None, flush=True)
