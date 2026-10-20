"""v-an-synth: count direct x-y gates in our gcert/1 lowerings and test them under an-synth's own mirrored-v model
(v_S . v_T = u_S . u_T; a gate joining x_T and y_S is legal two-sided iff u_S . u_T = 0). Own code; reads our own files.
Usage: python3 -I calc/xyprobe.py FILE..."""
import gzip, json, sys
from collections import Counter
def par(x): return bin(x).count('1') & 1
res = {}
for p in sys.argv[1:]:
    d = json.load(gzip.open(p, 'rt'))
    v, ports = d['v'], d['ports']
    def regs(g): return [g[2]] + [t for t, _, _ in g[3]] + (list(g[4]) if g[0] == 'out' else [])
    out = {}
    for ph in ('A', 'B'):
        n_xy = 0; bad = 0; same = 0; kinds = Counter()
        for g in d[ph]:
            rs = regs(g)
            xs = [r for r in rs if r < v]; ys = [r - v for r in rs if v <= r < 2 * v]
            if xs and ys:
                n_xy += 1
                for x in xs:
                    for y in ys:
                        if x == y: same += 1
                        if par(ports[x] & ports[y]): bad += 1
                kinds[g[0]] += 1
        out[ph] = dict(gates=len(d[ph]), direct_xy_gates=n_xy, xy_pairs_nonorth=bad, same_family_pairs=same, kinds=dict(kinds))
    res[p.split('/')[-1]] = out
    print(p.split('/')[-1], json.dumps(out))
json.dump(res, open('calc/xyprobe.json', 'w'), indent=1)
