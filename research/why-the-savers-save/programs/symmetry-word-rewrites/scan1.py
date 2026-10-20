import sys, time, random, json
sys.path.insert(0, '.')
import fw, search, reorder
from collections import Counter
import os; W = search.Word(fw.load(os.environ.get('WORD', 'K8-FULL-WORD.json.gz')))
O = reorder.OW(W)
c0 = O.total_cost(); print('cost0', c0, 'feasible', O.feasible, flush=True)
mv = O.candidate_moves()
print('candidate moves', len(mv), Counter(O.sec[g] for g,_ in mv), flush=True)
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)
random.shuffle(mv)
t = time.time(); stats = Counter(); deltas = Counter(); n=0
lim = int(sys.argv[2]) if len(sys.argv)>2 else 300
for g, p in mv[:lim]:
    # evaluate without accepting: accept never
    rec = {}
    def acc(d): rec['d'] = d; return False
    O.try_move(g, p, accept=acc)
    if 'd' in rec: stats['evaluated'] += 1; deltas[round(rec['d'],3)] += 1
    else: stats['infeasible'] += 1
    n += 1
print(stats, 'secs', round(time.time()-t,1))
print(sorted(deltas.items()))
