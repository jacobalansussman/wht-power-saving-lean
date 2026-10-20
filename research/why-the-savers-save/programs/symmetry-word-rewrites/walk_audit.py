"""Take K random cost-neutral commuting reorders (frames re-optimised), write the word, and audit it with the hub's
audit_local (py39 port): confirms that the search's legality model matches the hub's checker."""
import sys, random, json, subprocess
sys.path.insert(0, '.')
import fw, search, reorder
random.seed(int(sys.argv[1])); K = int(sys.argv[2]); out = sys.argv[3]
W = search.Word(fw.load('K8-FULL-WORD.json.gz')); O = reorder.OW(W)
regs = [r for r in O.chain if len(O.chain[r]) > 1]
done = 0; moved = 0
while done < K:
    r = random.choice(regs); gs = O.chain[r]; i = random.randrange(len(gs) - 1); g1, g2 = gs[i], gs[i + 1]
    if O.sec[g1] != O.sec[g2] or not O.commute(g1, g2): continue
    a, b = O.posof[g1], O.posof[g2]; between = O.order[a + 1:b]
    if not all(O.commute(g2, h) for h in between): continue
    d = O.try_move(g2, a, accept=lambda d: abs(d) < 1e-9)
    if d is not None: done += 1
nd = sum(1 for p, g in enumerate(O.order) if g != p)
print('accepted moves', done, 'gates displaced', nd, 'cost', O.total_cost(), dict(sorted(O.profile().items())))
w2 = O.to_word(); fw.save(w2, out)
print('saved', out)
