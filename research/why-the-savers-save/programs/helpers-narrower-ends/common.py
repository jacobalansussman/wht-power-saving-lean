"""d-hpe-multigate: shared helpers. F2^9 subspaces as reduced echelon bases (sorted tuple of ints)."""
import json, gzip
from collections import Counter, defaultdict
WORD = __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/hub-inputs/U8-CANONICAL-WORD.json.gz')
PHASES = ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']
def load(path=WORD):
    return json.loads(gzip.decompress(open(path, 'rb').read()))
def rref(vecs):
    bb = []
    for y in vecs:
        for z in bb: y = min(y, y ^ z)
        if y:
            bb = [min(z, z ^ y) for z in bb]; bb.append(y)
    # full reduction
    bb.sort(reverse=True)
    out = []
    for y in bb:
        for z in out: y = min(y, y ^ z)
        if y: out.append(y)
    # back-reduce
    out.sort(reverse=True)
    res = []
    for i, y in enumerate(out):
        for z in out:
            if z != y and (y ^ z) < y and (z.bit_length() < y.bit_length()):
                pass
        res.append(y)
    # canonical: reduced echelon via pivot elimination
    piv = {}
    for y in sorted(res, reverse=True):
        for p, z in piv.items():
            if (y >> p) & 1: y ^= z
        if y:
            p = y.bit_length() - 1
            for q in list(piv):
                if (piv[q] >> p) & 1: piv[q] ^= y
            piv[p] = y
    return tuple(sorted(piv.values(), reverse=True))
def join(a, b): return rref(list(a) + list(b))
def inside(a, b): return join(a, b) == b
def dim(a): return len(a)
FULL = rref([1 << i for i in range(9)])
def cls(r, v): return 'x' if r < v else ('y' if r < 2*v else 'h')
