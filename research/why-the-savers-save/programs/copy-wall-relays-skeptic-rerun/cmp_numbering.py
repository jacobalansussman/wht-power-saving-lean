import sys
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'hubcopy'))
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-relays'))
from e8units import hamming_e8
from patterns import hamming_e8_hub, instance, max_clique_nonknow
Lh, Gh = hamming_e8()
Lm, Gm = hamming_e8_hub()
print('lines identical:', [tuple(int(a) for a in r) for r in Lh] == [tuple(r) for r in Lm])
print('gram identical:', all(int(Gh[i][j]) == Gm[i][j] for i in range(120) for j in range(120)))
for name in ('P120', 'P111', 'P111ours'):
    k, G = instance(name)
    v = len(k)
    K = max_clique_nonknow(k, size=8)
    ok = all(not k[a][b] for a in K for b in K if a != b)
    sym = all(k[a][b] == k[b][a] for a in range(v) for b in range(v))
    print(name, 'v', v, 'known per row min/max', min(sum(r) for r in k), max(sum(r) for r in k), 'clique', len(K), 'pairwise non-knowing', ok, 'symmetric', sym)
# max clique overall (size 9 would break clique bound 8)
k, G = instance('P120')
print('P120 clique of 9 exists?', len(max_clique_nonknow(k, size=9)) >= 9)
