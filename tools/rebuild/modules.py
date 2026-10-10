"""modules.py: the three positive query modules of the paired-cube circuit, loaded from data files.

  triple module  : inputs = 3-subsets I of [p];   root J = sum of the inputs I disjoint from J
  pair module    : inputs = 2-subsets of [n];     root {i,j} = sum of the inputs disjoint from {i,j}
  all-but-one    : inputs 0..n-1;                 root i = sum of the inputs other than i
Every addition joins disjoint supports.  The modules are frozen data files of the outside repository (JSON
lists of additions).  They are loaded as data, and their contract (every root has exactly the support above,
every addition joins disjoint supports) is checked here.
"""
from itertools import combinations
import json


def _check(args, n_in, roots, want, one_based):
    off = 1 if one_based else 0
    supp = ([0] if one_based else []) + [1 << i for i in range(n_in)]
    for x in range(n_in + off, len(args)):
        a, b = args[x]
        assert off <= a < x and off <= b < x and not supp[a] & supp[b]
        supp.append(supp[a] | supp[b])
    assert len(roots) == len(want)
    for r, w in zip(roots, want):
        assert supp[r] == w, "module contract"
    return True


def triple_want(p):
    L = list(combinations(range(p), 3))
    return [sum(1 << i for i, I in enumerate(L) if not set(I) & set(J)) for J in L]


def pair_want(n):
    L = list(combinations(range(n), 2))
    return [sum(1 << k for k, q in enumerate(L) if not set(q) & set(ij)) for ij in L]


def load_triple(path, p):
    d = json.load(open(path))
    n = len(list(combinations(range(p), 3)))
    assert d["input_count"] == n and all(list(a) == [0, 0] for a in d["args"][:n + 1])
    _check(d["args"], n, d["roots"], triple_want(p), True)
    return dict(input_count=n, args=[list(a) for a in d["args"]], roots=list(d["roots"]))


def load_zero_based(path, n_in, want):
    d = json.load(open(path))
    assert d["input_count"] == n_in and all(a is None for a in d["args"][:n_in])
    _check(d["args"], n_in, d["roots"], want, False)
    return dict(input_count=n_in, args=[None if a is None else list(a) for a in d["args"]], roots=list(d["roots"]))
