import json, sys, os
from f2q import Field, Config, popc
from gb2 import supports_in, orders, Ring, system, msolve_text, singular_text
FIELDS = {m: Field(m) for m in (3, 4, 5)}
certs = json.load(open("certs.json"))
os.makedirs("jobs/gbt2", exist_ok=True)
for idx in [int(v) for v in sys.argv[1:]]:
    c = certs[idx]; cfg = Config(FIELDS[c["m"]], c["Ts"], c["moves"])
    U = [i for i in range(cfg.M) if c["cert"]["U"] >> i & 1]
    sup = supports_in(cfg, U); R = Ring(cfg.W, U)
    ords = list(orders(cfg, U))
    print(idx, c["tag"], "U", U, "dirs", [cfg.moves[i] for i in U], "supports in U", len(sup),
          "sizes", sorted(popc(d) for d in sup), "orders", len(ords))
    sysl = [system(cfg, R, o, sup) for o in ords[:3]]
    print(" polys per order", [len(s) for s in sysl], "vars used", len({v for p in sysl[0] for m in p for v in m}))
    open(f"jobs/gbt2/c{idx}_o0_65521.ms", "w").write(msolve_text(sysl[0], 65521))
    open(f"jobs/gbt2/c{idx}_Q3.sing", "w").write(singular_text(R, sysl, 0, f"c{idx}"))
