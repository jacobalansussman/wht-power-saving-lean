"""myxw.py -- v-ex2-limits independent re-implementation of the block-2 XW letters (written from the spec in the
xw2lib docstring, not from its code).  A monomial word is a pair (t, p): basis vector e_v -> i^p[v] e_t[v].
Address v = sum_j c_j 4^j, column c_j in 0..3, bit0 = row 0, bit1 = row 1.  Words act LEFT TO RIGHT."""
import numpy as np

def cols_of(f):
    v = np.arange(4 ** f, dtype=np.int64)
    return np.array([(v >> (2 * j)) & 3 for j in range(f)])          # shape (f, N)

def addr(C):
    out = np.zeros(C.shape[1], np.int64)
    for j in range(C.shape[0]):
        out += C[j].astype(np.int64) << (2 * j)
    return out

def bit(c, r): return (c >> r) & 1

# per-column maps as functions on column value arrays: (newcol(c), phase(c))
def col_letter(name):
    tab = {
        'X1': (lambda c: c ^ 1, None), 'X2': (lambda c: c ^ 2, None), 'X3': (lambda c: c ^ 3, None),
        'Z1': (None, lambda c: 2 * bit(c, 0)), 'Z2': (None, lambda c: 2 * bit(c, 1)),
        'Z3': (None, lambda c: 2 * (bit(c, 0) ^ bit(c, 1))),
        'S1': (None, lambda c: bit(c, 0)), 'S2': (None, lambda c: bit(c, 1)),
        'S3': (None, lambda c: bit(c, 0) ^ bit(c, 1)),
        'Om': (None, lambda c: np.ones_like(c)), 'Neg': (None, lambda c: 2 * np.ones_like(c)),
        'P': (lambda c: (bit(c, 0) << 1) | bit(c, 1), None),
        'C01': (lambda c: c ^ (bit(c, 0) << 1), None),        # row1 ^= row0
        'C10': (lambda c: c ^ bit(c, 1), None),               # row0 ^= row1
        'CZ': (None, lambda c: 2 * (bit(c, 0) & bit(c, 1))),
        'CS': (None, lambda c: bit(c, 0) & bit(c, 1)),
        'CSi': (None, lambda c: 3 * (bit(c, 0) & bit(c, 1))),
    }
    return tab[name]

def letter_map(name, f):
    """return function (C, ph) -> (C', ph') acting on all basis points at once (C = columns of current target)."""
    if name in ('1', ''):
        return lambda C, ph: (C, ph)
    def colwise(nm, js):
        perm, phase = col_letter(nm)
        def g(C, ph):
            C = C.copy(); ph = ph.copy()
            for j in js:
                if phase is not None: ph = ph + phase(C[j])
                if perm is not None: C[j] = perm(C[j])
            return C, ph % 4
        return g
    try:
        col_letter(name); return colwise(name, range(f))
    except KeyError:
        pass
    if name[0] in 'lf' and len(name) > 1:
        try:
            col_letter(name[1:]); return colwise(name[1:], [f - 1] if name[0] == 'l' else [0])
        except KeyError:
            pass
    if name in ('Pr0', 'Pr1'):
        r = int(name[2])
        return lambda C, ph: (C, (ph + (sum(bit(C[j], r) for j in range(f)) % 2)) % 4)
    def chain(r, s, mirror):           # row r of column j ^= row s of the PREVIOUS column (old values)
        order = list(range(f)) if not mirror else list(range(f - 1, -1, -1))
        def g(C, ph):
            old = C.copy(); C = C.copy()
            for t in range(1, f):
                C[order[t]] = C[order[t]] ^ (bit(old[order[t - 1]], s) << r)
            return C, ph
        return g
    def prefix(r, mirror):             # row r of column j := XOR of rows r of columns up to j (in order)
        order = list(range(f)) if not mirror else list(range(f - 1, -1, -1))
        def g(C, ph):
            C = C.copy(); acc = np.zeros(C.shape[1], np.int64)
            for t in range(f):
                j = order[t]; acc = acc ^ bit(C[j], r)
                C[j] = (C[j] & (3 ^ (1 << r))) | (acc << r)
            return C, ph
        return g
    def seq(*gs):
        def g(C, ph):
            for h in gs: C, ph = h(C, ph)
            return C, ph
        return g
    if name == 'A': return seq(chain(0, 0, False), chain(1, 1, False))
    if name == 'B': return seq(chain(0, 0, True), chain(1, 1, True))
    if name == 'Ai': return seq(prefix(0, False), prefix(1, False))
    if name == 'Bi': return seq(prefix(0, True), prefix(1, True))
    if name in ('a0', 'a1'): return chain(int(name[1]), int(name[1]), False)
    if name in ('b0', 'b1'): return chain(int(name[1]), int(name[1]), True)
    if name in ('a0i', 'a1i'): return prefix(int(name[1]), False)
    if name in ('b0i', 'b1i'): return prefix(int(name[1]), True)
    if name in ('E01', 'E10'): return chain(int(name[1]), int(name[2]), False)
    raise KeyError(name)

def word(w, f):
    """free word -> (t, p) with e_v -> i^p[v] e_t[v]"""
    C = cols_of(f); ph = np.zeros(C.shape[1], np.int64)
    for nm in [x for x in w.split('.') if x and x != '1']:
        C, ph = letter_map(nm, f)(C, ph)
    return addr(C), ph % 4
