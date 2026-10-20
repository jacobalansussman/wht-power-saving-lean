"""shrink.py (d-relay) -- find a SMALL sub-pattern of bit P111 (hub numbering) on which rank-8 relay decoders with
L layers and NO closure are still UNSAT.  Patterns are hereditary (a principal submatrix of a fitting relay decoder
fits the induced layering), so UNSAT on a sub-pattern implies UNSAT on P111 and P120 for the same L.
Greedy: drop ports one at a time (keeping one 8-clique K for the WLOG gauge X_K = I) while rank 8 stays UNSAT at L.
usage: venv/bin/python -I shrink.py L budget_seconds [seed] [start.json]   -> writes out/shrink_L{L}.json"""
import sys, json, time, random
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from patterns import instance, max_clique_nonknow
from relaysat import build
from pysat.solvers import Cadical153
L = int(sys.argv[1]); budget = float(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
knows, G = instance('P111'); v = len(knows)
K = max_clique_nonknow(knows, size=8)
out = 'out/shrink_L%d.json' % L
cur = list(range(v))
if len(sys.argv) > 4:
    d = json.load(open(sys.argv[4])); cur = d['ports']; K = d['K']
tried = set(json.load(open(out)).get('tried', [])) if len(sys.argv) > 4 and sys.argv[4] == out else set()
def unsat(ports):
    sub = [[knows[a][b] for b in ports] for a in ports]
    Ks = [ports.index(k) for k in K]
    pool, cl, X, Y, ge = build(sub, 8, L, 'none', Ks, True)
    s = Cadical153(bootstrap_with=cl); r = s.solve(); s.delete(); return not r
t0 = time.time(); rng = random.Random(seed)
assert unsat(cur), 'start pattern not UNSAT'
cand = [p for p in cur if p not in K and p not in tried]; rng.shuffle(cand)
for p in cand:
    if time.time() - t0 > budget: break
    trial = [q for q in cur if q != p]
    if unsat(trial): cur = trial
    tried.add(p)
    json.dump({'L': L, 'ports': cur, 'K': K, 'tried': sorted(tried), 'note': 'P111 hub numbering indices into live list'}, open(out, 'w'))
print('L=%d: UNSAT sub-pattern of P111 with %d ports (tried %d removals, %.0fs)' % (L, len(cur), len(tried), time.time() - t0))
