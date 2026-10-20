"""common.py (g-small): E8 lines and the bit/complex patterns, reusing u-floor/copyfloor.py."""
import sys
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../copy-wall-floor-checker'))
import copyfloor as cf

L, G = cf.e8_lines()          # 120 lines (norm 4 vectors), G = <a,b>/2 in {0,+-1,+-2}
V = 120
NONORTH = [[S != T and G[S][T] != 0 for T in range(V)] for S in range(V)]
ORTH = [[S != T and G[S][T] == 0 for T in range(V)] for S in range(V)]

def max_cliques(adj, size, verts=None):
    verts = list(range(len(adj))) if verts is None else verts
    out = []
    def rec(cur, cand):
        if len(cur) == size:
            out.append(tuple(cur)); return
        for i, c in enumerate(cand):
            if len(cur) + len(cand) - i < size: break
            rec(cur + [c], [d for d in cand[i + 1:] if adj[c][d]])
    rec([], verts)
    return out
