"""mincert.py: for every excluded instance in certs.json, the zero-block certificate whose move set U (D_F plus the
tight supports used) needs the fewest time orders (multiset permutations of U's directions, up to reversal).
Each new certificate is re-checked by ztest.check_cert. Writes certs_min.json (same records, 'cert' replaced,
'cert_old_U' kept)."""
import json, itertools, math
from collections import Counter
from f2q import Field, Config, popc
from ztest import links, cliques, check_cert

FIELDS = {m: Field(m) for m in (3, 4, 5)}
certs = json.load(open("certs.json"))


def norders(cfg, U):
    cnt = Counter(cfg.moves[i] for i in U)
    n = math.factorial(len(U))
    for v in cnt.values(): n //= math.factorial(v)
    return n


for c in certs:
    if not (c["cert"] and c["cert"]["kind"] == "Z"): continue
    cfg = Config(FIELDS[c["m"]], c["Ts"], c["moves"])
    tight, wit = links(cfg)
    cand = {}                                   # pair -> list of (D, Q) tight supports containing it
    for D, (Q, p) in tight.items():
        mv = [i for i in range(cfg.M) if D >> i & 1]
        for a in range(len(mv)):
            for b in range(a, len(mv)):
                cand.setdefault((mv[a], mv[b]), []).append((D, Q))
    cl = [y for y in cliques(cfg, wit) if len(y) >= 1]
    best = None
    for Qf in range(1, len(cfg.Dmask)):
        d = cfg.Dmask[Qf]; r = popc(cfg.pmask[Qf]); sd = popc(d)
        for y in cl:
            ym = sum(1 << i for i in y)
            if d & ym != ym or sd + r >= 2 * len(y): continue
            pairs = [(y[a], y[b]) for a in range(len(y)) for b in range(a, len(y))]
            opts = []
            for pr in pairs:
                inside = [(D, Q) for (D, Q) in cand[pr] if D & ~d == 0]
                opts.append(inside[:1] if inside else sorted(cand[pr], key=lambda t: popc(t[0]))[:6])
            for choice in itertools.product(*opts):
                U = d
                for (D, Q) in choice: U |= D
                Ul = [i for i in range(cfg.M) if U >> i & 1]
                key = (norders(cfg, Ul), popc(U), sd)
                if best is None or key < best[0]:
                    best = (key, dict(kind="Z", Y=list(y), F=Qf, DF=d, rF=r,
                                      pairs=[(pr[0], pr[1], ch[1], ch[0]) for pr, ch in zip(pairs, choice)], U=U))
    check_cert(cfg, best[1])
    c["cert_old_U"] = c["cert"]["U"]; c["cert"] = best[1]
    print(c["idx"], c["tag"], "orders", best[0][0], "|U|", best[0][1], "|D_F|", best[0][2],
          "(old |U|", popc(c["cert_old_U"]), ")", flush=True)
json.dump(certs, open("certs_min.json", "w"))
print("done")
