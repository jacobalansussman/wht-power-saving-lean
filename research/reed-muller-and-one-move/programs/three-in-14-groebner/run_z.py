"""run_z.py: span, deletion, clique and zero-block tests (my code) on
  (1) my own 44 block-5 3-in-14 sets (out_scan3_b5.txt, my coordinates),
  (2) rm2-solve's 44 (data_rm2_44.json, its coordinates),
  (3) must-fail controls through the same code: block-4 3 in 11 (both sets) and block-3 3 in 8 (my scans),
  (4) must-pass controls: every plain multiset behind (1) and (3) (W in 5W, a genuine program),
      U_3 = 2 in 7 (block 3, targets I, I, all seven directions; real sharing), plain 1 in 5 with T_0.
Writes certs.json (one certificate per instance) and prints a table."""
import json, re, sys, time
from f2q import Field, Config, popc
from ztest import zero_block, check_cert

FIELDS = {m: Field(m) for m in (3, 4, 5)}
inst = []
def from_scan(fn, m, tag):
    F = FIELDS[m]
    for line in open(fn):
        if not line.startswith("PASS"): continue
        ks = [int(v) for v in re.search(r"ks=\(([0-9, ]+)\)", line).group(1).split(",")]
        moves = [int(v) for v in re.search(r"moves=\[([0-9, ]+)\]", line).group(1).split(",")]
        z = int(re.search(r"delete dir (\d+)", line).group(1))
        inst.append(dict(tag=tag, m=m, ks=ks, Ts=[F.target(k) for k in ks], moves=moves, src=line.split(" moves=")[0],
                         expect="fail"))
        inst.append(dict(tag=tag + "-plain", m=m, ks=ks, Ts=[F.target(k) for k in ks], moves=sorted(moves + [z]),
                         src="plain multiset of: " + line.split(" moves=")[0], expect="pass"))
from_scan("out_scan3_b5.txt", 5, "mine-b5")
for d in json.load(open("data_rm2_44.json")):
    inst.append(dict(tag="rm2-b5", m=5, ks=d["ks"], Ts=d["Ts"], moves=d["moves"], src=d["line"], expect="fail"))
from_scan("out_scan3_b4.txt", 4, "ctl-b4-3in11")
from_scan("out_scan3_b3.txt", 3, "ctl-b3-3in8")
I3 = [[int(i == j) for j in range(3)] for i in range(3)]
inst.append(dict(tag="ctl-U3", m=3, ks=[0, 0], Ts=[I3, I3], moves=list(range(1, 8)), src="U_3 = 2 in 7", expect="pass"))
inst.append(dict(tag="ctl-plain1in5", m=5, ks=[0], Ts=[FIELDS[5].target(0)], moves=FIELDS[5].frame_canonical(0),
                 src="plain 1 in 5", expect="pass"))

t0 = time.time(); certs = []; tally = {}
for k, d in enumerate(inst):
    F = FIELDS[d["m"]]
    cfg = Config(F, d["Ts"], d["moves"])
    span, dele = cfg.span_deletion()
    res = zero_block(cfg) if (span and dele) else dict(clique=None, C=None, Z=None)
    excluded = (not span) or (not dele) or (res.get("C") is False) or (res.get("Z") is False)
    verdict = "EXCLUDED" if excluded else "passes"
    ok = (verdict == "EXCLUDED") == (d["expect"] == "fail")
    if res.get("cert") and res["cert"]["kind"] == "Z":
        check_cert(cfg, res["cert"])
    key = (d["tag"], verdict); tally[key] = tally.get(key, 0) + 1
    c = res.get("cert")
    desc = ""
    if c and c["kind"] == "Z":
        desc = (f"Y={c['Y']} |D_F|={popc(c['DF'])} r_F={c['rF']} moves-in-cert={popc(c['U'])} "
                f"pairs-forms={len({p[2] for p in c['pairs']})}")
    print(f"[{k}] {d['tag']} {d['src']} | span={span} del={dele} clique={res.get('clique')} C={res.get('C')} "
          f"Z={res.get('Z')} nfail={res.get('nfail')} -> {verdict} {'(as expected)' if ok else '(UNEXPECTED)'} {desc}",
          flush=True)
    certs.append(dict(idx=k, tag=d["tag"], m=d["m"], ks=d["ks"], Ts=d["Ts"], moves=d["moves"], src=d["src"],
                      expect=d["expect"], verdict=verdict, span=span, deletion=dele, clique=res.get("clique"),
                      cert=c))
json.dump(certs, open("certs.json", "w"))
print("TALLY", sorted(tally.items()), f"{time.time() - t0:.0f}s")
