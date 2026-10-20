"""t-egraph: frame-word library for the k8 word (hub audit_local semantics, re-implemented for search).

A word is the hub FULL-WORD dict.  The climbed sequence is A, then (centre paths), then B, then the tail
(early_cleanup gates with their replacement frames), then every register's endpoint.  Frames are canonical
reduced bases (tuples), as in audit_local.basis.  Nothing here is trusted for a claim: every candidate is
re-checked by the hub's audit_local (py39 port, only int.bit_count replaced).
"""
import gzip, json, math
from collections import Counter, defaultdict

def basis(vs):
    piv = {}
    for a in vs:
        while a:
            bit = a.bit_length() - 1
            if bit in piv: a ^= piv[bit]
            else:
                piv[bit] = a; break
    for k in sorted(piv):
        for j in piv:
            if j > k and piv[j] >> k & 1: piv[j] ^= piv[k]
    return tuple(piv[k] for k in sorted(piv, reverse=True))

def contained(a, b): return basis(a + b) == b
def join(a, b): return basis(a + b)

def popc(x): return bin(x).count('1')

def load(path):
    return json.loads(gzip.decompress(open(path, 'rb').read()))

def save(w, path):
    data = json.dumps(w, separators=(',', ':')).encode()
    with open(path, 'wb') as f: f.write(gzip.compress(data, mtime=0))

def registers(e): return e[1:3] if e[0] == 'add' else e[1:2]

def geometry(w):
    v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2 * v + R
    full = basis([1 << i for i in range(h)])
    end = [full] * v + [basis([x for x in range(1, 1 << h) if popc(x & p) % 2 == 0]) for p in w['ports']] + [full] * R
    for q, U in w.get('helper_final', {}).items(): end[int(q)] = basis(U)
    start = [(p,) for p in w['ports']] + [()] * (v + R)
    for q, U in w.get('start_extra', {}).items(): start[int(q)] = basis(U)
    return v, h, R, n, full, start, end

def climbed(w):
    """The climbed gate list: [(section, index, gate-with-frame)] in climb order (A, B, tail)."""
    out = []
    for i, e in enumerate(w['A']): out.append(('A', i, e))
    for i, e in enumerate(w['B']): out.append(('B', i, e))
    for j, (i, U) in enumerate(w.get('early_cleanup', [])):
        out.append(('T', j, w['cleanup'][i][:-1] + [U]))
    return out

def F(d, m=73): return d * math.log(m / d) if d else 0.0

def paid_profile(w, frames=None):
    """Exact re-run of audit_local's paid Counter (climb sizes). frames: optional override list aligned with climbed()."""
    v, h, R, n, full, start, end = geometry(w)
    seq = climbed(w)
    cur = list(start); paid = Counter()
    def climb(q, U):
        U = basis(U)
        assert contained(cur[q], U), ('not monotone', q, cur[q], U)
        d = len(U) - len(cur[q])
        if d: paid[d] += 1
        cur[q] = U
    nA = len(w['A'])
    for k, (sec, i, e) in enumerate(seq):
        if k == nA:
            lev = {int(t): basis(U) for t, U in w.get('target_levels', {}).items()}
            assert all(cur[v + t] == lev.get(t, ()) for t in range(v))
            for q, U in zip(w['centers'], w['center_frames']):
                U = basis(U); assert cur[q] == U, ('center', q)
                served = {lev.get(g[1] - v, ()) for g in w['cut'] if g[2] == q}
                for T in sorted(served, key=lambda x: -len(x)):
                    if T != U: paid[len(U) - len(T)] += 1
                    U = T
        U = frames[k] if frames is not None else e[-1]
        for q in registers(e): climb(q, U)
    for q, U in enumerate(end): climb(q, U)
    return paid

def chains(w, frames=None):
    """Per register: list of (k, frame) over climbed positions k touching it (A,B,tail), plus start/end."""
    v, h, R, n, full, start, end = geometry(w)
    seq = climbed(w)
    ch = defaultdict(list)
    for k, (sec, i, e) in enumerate(seq):
        U = basis(frames[k] if frames is not None else e[-1])
        for q in registers(e): ch[q].append((k, U))
    return ch, start, end

def cost_of_profile(p, m=73):
    return sum(n * F(d, m) for d, n in p.items())
