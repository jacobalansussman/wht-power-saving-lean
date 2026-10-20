import sys
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-relays'))
from collections import Counter
from patterns import e8_ours, hamming_e8_hub
def prof(G, P):
    live = [i for i in range(120) if i not in P]
    c = Counter(sum(1 for p in P if G[i][p] != 0) for i in live)
    d = Counter(sum(1 for q in P if q != p and G[p][q] != 0) for p in P)
    return sorted(c.items()), sorted(d.items())
L, G = e8_ours()
print('P111ours  :', prof(G, [8, 9, 10, 12, 16, 17, 24, 32, 64]))
Lg, Gg = hamming_e8_hub()
print('#589 true :', prof(Gg, [9, 31, 37, 52, 66, 76, 78, 111, 119]))
print('builder P :', prof(Gg, [39, 52, 57, 71, 76, 78, 103, 109, 114]))
