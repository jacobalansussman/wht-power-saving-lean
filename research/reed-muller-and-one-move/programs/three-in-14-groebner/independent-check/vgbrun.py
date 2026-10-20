"""vgbrun.py (v-r3-3in14) -- runs on a Colab.  For MY 44 sets (vcerts_local.json: U = a 5-move zero-block
certificate support whose tight pairs lie inside U), EVERY one of the 120 permutations of U (no symmetry
reduction), my own equations (vpoly.equations: X_D = A (I+G)^-1 B^T = diag p for every form with support in U):
  Singular std over Q -> UNIT / NOT-UNIT for every order;
  for LIFTN sampled orders per set: lift(I, 1) and MY exact check sum c_i f_i = 1 (vpoly.parse_sing);
  for MSN sampled orders per set: msolve -g 2 over 65521 and 1073741827 (must print [1]).
Must-pass control (mode plain): the plain 3 in 15 (deleted move restored), U+ = U plus every copy of the deleted
direction, 2 random orders per set: must NOT give 1.
Usage: python3 vgbrun.py <mode main|plain> <shard> <nshards>  -> res_<mode>_<shard>.jsonl"""
import json, sys, random, itertools, subprocess, os, time
from vpoly import VRing, equations, to_sing, to_msolve, parse_sing, pmul, padd

mode, shard, nsh = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
LIFTN, MSN = int(os.environ.get("LIFTN", "20")), 5
certs = json.load(open(os.environ.get("CERTS", "vcerts_local.json")))
plain = json.load(open("vplain_sups.json")) if mode == "plain" else None
random.seed(1000 + shard)
TAG = os.environ.get("TAG", "")
outf = open(f"res_{mode}{TAG}_{shard}.jsonl", "a")


def singular(script, tmo):
    fn = f"tmp_{mode}{TAG}_{shard}.sing"; open(fn, "w").write(script)
    try:
        r = subprocess.run(["Singular", "-q", fn], capture_output=True, text=True, timeout=tmo)
        if os.environ.get("RAW"): open(f"raw_{mode}{TAG}_{shard}.txt", "a").write(f"=== rc {r.returncode}\n" + r.stdout[-3000:] + "\n--- stderr\n" + r.stderr[-3000:] + "\n")
        return r.stdout
    except subprocess.TimeoutExpired:
        return "TIMEOUT"


def msolve(polys, prime):
    uv = sorted({v for p in polys for m in p for v in m}); ren = {v: k for k, v in enumerate(uv)}
    txt = ",".join(f"x{k + 1}" for k in range(len(uv))) + f"\n{prime}\n" + ",\n".join(to_msolve(p, ren) for p in polys) + "\n"
    fi, fo = f"tmp_{mode}{TAG}_{shard}.ms", f"tmp_{mode}{TAG}_{shard}.msout"
    open(fi, "w").write(txt)
    if os.path.exists(fo): os.remove(fo)
    r = subprocess.run(["msolve", "-g", "2", "-t", "1", "-f", fi, "-o", fo], capture_output=True, timeout=600)
    txt = open(fo).read() if os.path.exists(fo) else ""
    body = "".join(l for l in txt.splitlines() if not l.startswith("#")).strip()
    if body == "[1]:": return True
    hdr = [l for l in txt.splitlines() if l.startswith("#length of basis")]
    # anything but a complete basis header + body is reported as a string, never as a plain False
    return f"BASIS {hdr[0].split(':')[1].strip()}" if (hdr and body.endswith("]:")) else f"BAD rc={r.returncode} len={len(txt)}"


for c in certs:
    n = c["n"]
    if (n if isinstance(n, int) else 0) % nsh != shard: continue
    if isinstance(n, int) and n >= int(os.environ.get("NMAX", "999")): continue
    t0 = time.time()
    if mode == "main":
        U = c["U"]; sups = [(D, tuple(p)) for D, p in c["sups"]]
        orders = [list(U)] if os.environ.get("FIXORDER") else [list(o) for o in itertools.permutations(U)]
    else:
        P = plain[str(n)]; U = P["U"]; sups = [(D, tuple(p)) for D, p in P["sups"]]
        orders = [random.sample(U, len(U)) for _ in range(2)]
    W = len(sups[0][1]); R = VRing(W, U)
    systems = [equations(R, o, sups) for o in orders]
    for f in systems:
        for p in f:
            for m in p: assert len(set(m)) == len(m)      # multilinear
    lift_idx = set(random.sample(range(len(orders)), min(LIFTN, len(orders)))) if mode == "main" and not os.environ.get("NOLIFT") else set()
    lines = [f"ring r = 0, (x(1..{len(R.names)})), dp;", "option(redSB);", "short = 0;", "ideal I; ideal G; matrix T; int i;"]
    for k, polys in enumerate(systems):
        lines.append("I = " + ",\n".join(to_sing(p) for p in polys) + ";")
        lt = ""
        if k in lift_idx:
            lt = (f' T = lift(I, ideal(1)); "LIFT {k} BEGIN"; for (i = 1; i <= ncols(I); i++) '
                  f'{{ "C " + string(i) + " " + string(T[i, 1]); }} "LIFT {k} END";')
        lines.append(f'G = std(I); if (G[1] == 1) {{ "ORDER {k}: UNIT";{lt} }} else {{ "ORDER {k}: NOT-UNIT dim " + string(dim(G)); }}')
    lines.append("quit;")
    so = "SKIPPED" if os.environ.get("SKIPSING") else singular("\n".join(lines) + "\n", 1500 if mode == "main" else 300)
    verdict = {}; lifts = {}; cur = None
    for line in so.splitlines():
        if line.startswith("ORDER "):
            k = int(line.split()[1].rstrip(":")); verdict[k] = "UNIT" if "NOT-UNIT" not in line else line.split(": ")[1]
        elif line.startswith("LIFT ") and line.endswith("BEGIN"):
            cur = int(line.split()[1]); lifts[cur] = []
        elif line.startswith("LIFT ") and line.endswith("END"):
            cur = None
        elif cur is not None and line.startswith("C "):
            _, i, poly = line.split(" ", 2); lifts[cur].append((int(i), poly))
    # exact check of the cofactor identities with MY polynomials
    lift_ok = {}
    for k, cof in lifts.items():
        polys = systems[k]; assert len(cof) == len(polys), (len(cof), len(polys))
        s = {}
        for (i, ps), f in zip(cof, polys): s = padd(s, pmul(parse_sing(ps), f))
        lift_ok[k] = (s == {(): 1})
    ms = {}
    if mode == "main" and not os.environ.get("NOMS"):
        for k in random.sample(range(len(orders)), min(MSN, len(orders))):
            ms[k] = [msolve(systems[k], 65521), msolve(systems[k], 1073741827)]
    else:
        for k in range(len(orders)):
            try:
                ms[k] = [msolve(systems[k], 65521)]
            except Exception as e:
                ms[k] = [f"ERR {type(e).__name__}"]
    rec = {"n": n, "mode": mode, "orders": len(orders), "unit": sum(v == "UNIT" for v in verdict.values()),
           "decided": len(verdict), "timeout": so == "TIMEOUT", "lifts_checked": len(lift_ok),
           "lifts_ok": sum(lift_ok.values()), "msolve": {str(k): v for k, v in ms.items()},
           "nonunit": {str(k): v for k, v in verdict.items() if v != "UNIT"}, "secs": round(time.time() - t0, 1),
           "neq": [len(f) for f in systems[:1]], "nvars_used": len({v for p in systems[0] for m in p for v in m})}
    outf.write(json.dumps(rec) + "\n"); outf.flush()
    print(json.dumps(rec), flush=True)
print("DONE", flush=True)
