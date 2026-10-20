# d-bit-exits shared frame helpers: exact Q-spans of E8 root lines and closure in the live labels, as xcheck_p.py does.
from fractions import Fraction as _Q


def _rank(vecs):
    B = {}
    for rr in vecs:
        row = list(map(_Q, rr))
        for c in sorted(B):
            q = row[c]
            if q: row = [a - q * b for a, b in zip(row, B[c])]
        c = next((i for i, x in enumerate(row) if x), None)
        if c is not None:
            q = row[c]; B[c] = [x / q for x in row]
    return B


def _inspan(B, vec):
    row = list(map(_Q, vec))
    for c in sorted(B):
        q = row[c]
        if q: row = [a - q * b for a, b in zip(row, B[c])]
    return not any(row)


class Frames:
    def __init__(self, w):
        self.w = w; self.V = w['v']; self.L = [tuple(r) for r in w['lines']]
        self.mask = [sum(1 << i for i in f) for f in w['frames']]
        self.dim = [len(_rank([self.L[i] for i in f])) for f in w['frames']]
        self.key = {m: k for k, m in enumerate(self.mask)}
        self.FULL = next(k for k, d in enumerate(self.dim) if d == 8)
        self.ZERO = next(k for k, d in enumerate(self.dim) if d == 0)
        self.new = []

    def le(self, a, b):
        return self.mask[a] & self.mask[b] == self.mask[a]

    def frame_of_labels(self, labels):
        B = _rank([self.L[i] for i in labels])
        m = sum(1 << i for i in range(self.V) if _inspan(B, self.L[i]))
        if m in self.key: return self.key[m]
        k = len(self.mask); self.mask.append(m); self.dim.append(len(B)); self.key[m] = k; self.new.append(k)
        return k

    def join(self, fs):
        m = 0
        for f in fs: m |= self.mask[f]
        if m in self.key: return self.key[m]
        return self.frame_of_labels([i for i in range(self.V) if m >> i & 1])

    def labels(self, f):
        return [i for i in range(self.V) if self.mask[f] >> i & 1]
