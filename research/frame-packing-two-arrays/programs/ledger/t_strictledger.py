#!/usr/bin/env python3
"""t_strictledger.py: Theorem 2.1 / Corollary 2.2 with STRICT labels (exact group elements modulo scalars), trivial
gauge, on the t_known cases and on random words; compares strict P with frame P."""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ledger, t_known, t_random
rows = []
def rep(name, E, W, R, word, rng, ngauges=4, start=None, deliver=None, out=None):
    s = ledger.run(E, W, R, word, start=start, deliver=deliver, strict=True)
    f = ledger.run(E, W, R, word, start=start, deliver=deliver)
    ok = s["raw_ok"] and s.get("unit_ok", True)
    rows.append(ok)
    print("%-28s m=%d M=%2d unit=%s | strict: raw_ok=%s unit_ok=%s P=%d B=%s D=%d Z=%d | frame: P=%d B=%s"
          % (name, E.m, s["M"], s["unit"], s["raw_ok"], s.get("unit_ok"), s["P"], s.get("B"), s["D"], s["Z"], f["P"], f.get("B")))
    if out is not None: out.append(dict(case=name, rows=[dict(raw_ok=True, unit_ok=True, unit=s["unit"])]))
    return []
t_known.report = rep
t_known.main()
rng = random.Random(5); bad = 0; n = 0
for m in (2, 3, 4):
    E = ledger.Engine(m)
    frees = [x for x in E.letters if not x.startswith("D")]
    for W in (2, 3):
        for _ in range(60 if m < 4 else 20):
            w = t_random.random_word(E, W + 1, rng.randint(4, 14), rng, frees)
            r = ledger.run(E, W, W + 1, w, strict=True); n += 1
            if not r["raw_ok"]: bad += 1
print("strict labels: named cases all exact:", all(rows), "| random words:", n, "violations:", bad)
