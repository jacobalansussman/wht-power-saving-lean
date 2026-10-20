"""run_gb2.py SHARD NSHARDS [--count]: the exact-algebra decision (gb2.py) for every instance in certs.json.

Per instance, the move set U and the orders:
- excluded instances (zero-block certificate): U = the certificate's moves; every relative order of U;
- plain controls (W in 5W, a genuine program in every time order): U' = the same directions as the matching
  excluded instance's U plus the restored copy of the deleted direction; 5 sampled relative orders (must all be
  NOT-UNIT: this shows that the algebra sees the deletion and does not refuse real programs);
- plain 1 in 5: all 5 moves, every order (must be NOT-UNIT);
- U_3 = 2 in 7: all 7 moves, lex order (known to exist there; must be NOT-UNIT).
Each (instance, order): Singular std over Q (all orders of an instance in one Singular process) and msolve with -g 2
over 65521 and 1073741827.  For every order that Singular calls UNIT it also writes lift(I, 1) over Q, and the
cofactor identity sum c_i f_i = 1 is checked exactly in Python (verify_lift).
Output: res_gb2_SHARD.jsonl (one line per instance)."""
import json, sys, os, subprocess, random, time, re, tempfile
from fractions import Fraction
from f2q import Field, Config, popc
from gb2 import supports_in, orders, Ring, system, pstr, msolve_text, singular_text

shard, nsh = int(sys.argv[1]), int(sys.argv[2]); COUNT = "--count" in sys.argv
FIELDS = {m: Field(m) for m in (3, 4, 5)}
certs = json.load(open("certs_min.json"))
PR = (65521, 1073741827)


def pick(moves, dirs):
    """indices in moves realising the multiset of directions dirs (first copies)."""
    used = set(); out = []
    for z in dirs:
        i = next(i for i in range(len(moves)) if moves[i] == z and i not in used); used.add(i); out.append(i)
    return sorted(out)


jobs = []
for c in certs:
    if c["cert"] and c["cert"]["kind"] == "Z":
        U = [i for i in range(len(c["moves"])) if c["cert"]["U"] >> i & 1]
        jobs.append(dict(c=c, U=U, mode="all", expect="UNIT"))
    elif c["tag"].endswith("-plain"):
        ex = certs[c["idx"] - 1]; assert ex["tag"] + "-plain" == c["tag"] and ex["cert"]["kind"] == "Z"
        z = int(re.search(r"delete dir (\d+)", ex["src"]).group(1))
        dU = [ex["moves"][i] for i in range(len(ex["moves"])) if ex["cert"]["U"] >> i & 1]
        jobs.append(dict(c=c, U=pick(c["moves"], dU + [z]), mode="sample5", expect="NOT-UNIT"))
    elif c["tag"] == "ctl-plain1in5":
        jobs.append(dict(c=c, U=list(range(5)), mode="all", expect="NOT-UNIT"))
    elif c["tag"] == "ctl-U3":
        jobs.append(dict(c=c, U=list(range(7)), mode="lex", expect="NOT-UNIT"))
import math
def cost(j):
    from math import factorial
    return 5 if j["mode"] == "sample5" else (1 if j["mode"] == "lex" else factorial(len(j["U"])))
load = [0] * nsh; assign = []
for j in sorted(jobs, key=cost, reverse=True):          # greedy balance by (an upper bound on) the order count
    s_ = min(range(nsh), key=lambda t: load[t]); load[s_] += cost(j); assign.append((s_, j))
# excluded instances first (the result), then the controls; U_3 (the slowest) last
mine = sorted([j for s_, j in assign if s_ == shard],
              key=lambda j: (j["expect"] != "UNIT", j["mode"] == "lex", j["mode"] == "sample5", j["c"]["idx"]))
TQ = lambda j: 2400 if j["expect"] == "UNIT" else 600        # Singular time limit per instance (s)
TM = lambda j: 1200 if j["expect"] == "UNIT" else 300        # msolve time limit per order (s)


def parse_ms(path):
    s = open(path).read()
    if re.search(r"length of basis:\s+1 element", s) and re.search(r"^\[1\]:", s, re.M): return "UNIT"
    if "[-1]" in s: return "UNIT"
    if "#length of basis" in s: return "NOT-UNIT"
    return "ERROR"


def parse_poly(s, nv):
    """Singular output polynomial (short=0) -> dict mono(tuple of exps) -> Fraction."""
    s = s.replace(" ", "")
    if s == "0": return {}
    terms = re.findall(r"[+-]?[^+-]+", s); out = {}
    for t in terms:
        sign = -1 if t[0] == "-" else 1; t = t.lstrip("+-")
        coef = Fraction(1); exps = [0] * nv
        for f in t.split("*"):
            mm = re.fullmatch(r"x\((\d+)\)(?:\^(\d+))?", f)
            if mm: exps[int(mm.group(1)) - 1] += int(mm.group(2) or 1)
            else: coef *= Fraction(f)
        key = tuple(exps); out[key] = out.get(key, 0) + sign * coef
    return {k: v for k, v in out.items() if v}


def verify_lift(cof_strs, polys, nv):
    """exact check of sum cof_i * f_i == 1."""
    tot = {}
    for cs, f in zip(cof_strs, polys):
        cp = parse_poly(cs, nv)
        if not cp: continue
        fp = {}
        for mono, c in f.items():
            e = [0] * nv
            for v in mono: e[v] += 1
            fp[tuple(e)] = Fraction(c)
        for m1, c1 in cp.items():
            for m2, c2 in fp.items():
                k = tuple(a + b for a, b in zip(m1, m2)); tot[k] = tot.get(k, 0) + c1 * c2
    tot = {k: v for k, v in tot.items() if v}
    return tot == {tuple([0] * nv): Fraction(1)}


out = open(f"res_gb2_{shard}.jsonl", "a")
tmp = tempfile.mkdtemp(prefix=f"gb2s{shard}_", dir=".")
for j in mine:
    c = j["c"]; t0 = time.time()
    cfg = Config(FIELDS[c["m"]], c["Ts"], c["moves"]); U = j["U"]
    sup = supports_in(cfg, U); R = Ring(cfg.W, U)
    if j["mode"] == "lex": ords = [sorted(U)]
    else:
        ords = list(orders(cfg, U))
        if j["mode"] == "sample5" and len(ords) > 5:
            ords = random.Random(77 + c["idx"]).sample(ords, 5)
    if COUNT:
        print(c["idx"], c["tag"], "U", len(U), "orders", len(ords), "supports", len(sup), flush=True); continue
    systems = [system(cfg, R, o, sup) for o in ords]
    # Singular over Q, all orders
    sp = os.path.join(tmp, "q.sing")
    open(sp, "w").write(singular_text(R, systems, 0, f"i{c['idx']}", lift=(j["expect"] == "UNIT")))
    try:
        so = subprocess.run(["Singular", "-q", sp], capture_output=True, text=True, timeout=TQ(j)).stdout
    except subprocess.TimeoutExpired:
        so = "SINGULAR-TIMEOUT"
    qv = {}
    for mm in re.finditer(r"ORDER (\d+): (UNIT|NOT-UNIT)", so): qv[int(mm.group(1))] = mm.group(2)
    # msolve over two primes, every order
    mv = {p: {} for p in PR}
    for k, sy in enumerate(systems):
        for p in PR:
            fi = os.path.join(tmp, "m.ms"); fo = os.path.join(tmp, "m.out")
            open(fi, "w").write(msolve_text(sy, p))
            if os.path.exists(fo): os.remove(fo)
            try:
                subprocess.run(["msolve", "-f", fi, "-o", fo, "-g", "2", "-t", "1"], capture_output=True, timeout=TM(j))
                mv[p][k] = parse_ms(fo) if os.path.exists(fo) else "ERROR"
            except subprocess.TimeoutExpired:
                mv[p][k] = "TIMEOUT"
    # lift certificates over Q for EVERY order that Singular calls UNIT; each identity checked exactly here
    lifts = dict(checked=0, ok=0, failed=[], denominators=set())
    NV = len(R.names)
    for k in range(len(ords)):
        if qv.get(k) != "UNIT" or j["expect"] != "UNIT": continue
        mm = re.search(rf"LIFT {k} BEGIN\n(.*?)\nLIFT {k} END", so, re.S)
        lifts["checked"] += 1
        body = mm.group(1).strip().split("\n") if mm else []
        good = bool(mm) and len(body) == len(systems[k]) and verify_lift(body, systems[k], NV)
        if good: lifts["ok"] += 1
        else: lifts["failed"].append(k)
        if mm: lifts["denominators"] |= {Fraction(x).denominator for x in re.findall(r"(?<![\w(])\d+/\d+", mm.group(1))}
    lifts["denominators"] = sorted(lifts["denominators"])
    res = dict(idx=c["idx"], tag=c["tag"], src=c["src"], expect=j["expect"], U=U, dirsU=[cfg.moves[i] for i in U],
               nsupports=len(sup), norders=len(ords),
               Q=dict(unit=sum(v == "UNIT" for v in qv.values()), not_unit=sum(v == "NOT-UNIT" for v in qv.values()),
                      missing=len(ords) - len(qv)),
               msolve={str(p): dict(unit=sum(v == "UNIT" for v in mv[p].values()),
                                    not_unit=sum(v == "NOT-UNIT" for v in mv[p].values()),
                                    error=sum(v == "ERROR" for v in mv[p].values())) for p in PR},
               agree=all(qv.get(k) == mv[p][k] for k in range(len(ords)) for p in PR),
               singular_timeout=(so == "SINGULAR-TIMEOUT"),
               msolve_timeouts=sum(v == "TIMEOUT" for p in PR for v in mv[p].values()),
               lifts=lifts, secs=round(time.time() - t0, 1))
    allv = [qv.get(k) for k in range(len(ords))] + [mv[p][k] for p in PR for k in range(len(ords))]
    res["verdict"] = ("UNIT" if all(v == "UNIT" for v in allv) else
                      "NOT-UNIT" if all(v == "NOT-UNIT" for v in allv) else "MIXED")
    res["as_expected"] = (res["verdict"] == j["expect"]) and (j["expect"] != "UNIT" or lifts["ok"] == len(ords))
    out.write(json.dumps(res) + "\n"); out.flush()
    print(json.dumps(res), flush=True)
print("SHARD DONE", shard, flush=True)
