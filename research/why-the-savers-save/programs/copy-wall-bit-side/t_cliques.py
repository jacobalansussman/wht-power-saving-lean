import sys, time, collections; sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from e8geo import *
for name, P in (('P120', []), ('P111-s2', PUNCT_S2), ('P111-word', PUNCT_WORD)):
    live = [i for i in range(NL) if i not in P]
    t = time.time()
    C8 = cliques(live, 8)
    C9 = cliques(live, 9) if name == 'P120' else []
    per = collections.Counter(T for c in C8 for T in c)
    print(name, 'v', len(live), '8-cliques', len(C8), '9-cliques', len(C9), 'per-line counts', sorted(collections.Counter(per.values()).items()), round(time.time() - t, 1), flush=True)
