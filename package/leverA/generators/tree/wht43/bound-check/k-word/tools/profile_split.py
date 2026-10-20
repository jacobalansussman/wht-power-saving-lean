"""k-word: where the gcert lowering's histogram differs from the hub's paid local profile.
Recounts the hub's climb profile (hpe_audit_and_transport's rule: every gate climbs its registers to its frame, centres
pay their served descents at the cut, final climbs) on the canonical word (a) as is, (b) without the tail gates.
usage: profile_split.py CANONICAL-WORD.json.gz"""
import sys, json, gzip
from collections import Counter
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
v = c['v']; F = c['frames']
def basis(vs):
    piv = {}
    for x in vs:
        for k in sorted(piv, reverse=True):
            if x >> k & 1: x ^= piv[k]
        if x:
            k = x.bit_length() - 1
            for j in piv:
                if piv[j] >> k & 1: piv[j] ^= x
            piv[k] = x
    return tuple(piv[k] for k in sorted(piv, reverse=True))
def prof(use_tail, served=True):
    cur = list(c['start']); paid = Counter()
    lev = {int(t): f for t, f in c['target_levels'].items()}
    def climb(r, f):
        o = cur[r]
        if o != f:
            assert basis(list(F[o]) + list(F[f])) == tuple(F[f])
            paid[len(F[f]) - len(F[o])] += 1; cur[r] = f
    for part in ('A', 'B'):
        if part == 'B' and served:
            for r, fr in zip(c['centers'], c['center_frames']):
                sv = {lev.get(e[1] - v, 0) for e in c['cut'] if e[2] == r}
                old = fr
                for f in sorted(sv, key=lambda f: -len(F[f])):
                    if old != f: paid[len(F[old]) - len(F[f])] += 1; old = f
        for e in c[part]:
            for r in [e[1]] + ([e[2]] if e[0] == 'add' else []): climb(r, e[-1])
    if use_tail:
        for e in c['tail']:
            for r in [e[1]] + ([e[2]] if e[0] == 'add' else []): climb(r, e[-1])
    for r, f in enumerate(c['final']): climb(r, f)
    return dict(sorted(paid.items()))
print('hub local_profile        ', dict(sorted((int(k), n) for k, n in c['local_profile'].items())))
print('recount with tail        ', prof(True))
print('recount without tail     ', prof(False))
print('without tail, no served  ', prof(False, served=False))
