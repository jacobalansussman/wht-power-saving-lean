#!/usr/bin/env python3
"""cx-mirror / counts.py (own code, stdlib).  an2-word's per-line crossing formulas, applied call by call.

  GF (y += x):  paid climbs of v >= #{b active : v in Y_b(call start), v never in X_b during the call}
  GB (x -= y):  paid drops  of v >= #{b active : v in Y_b(call start), v not in X_b(call end)}

Lines v_b = s + g_b with the natural mirror g_b = j + u_b (lines.py: passes the junk orthogonality test).  The
local space of a call is F_2^10 (window bits 0..8, s = bit 9).  Every v in F_2^10 is tested, not only the 112 lines.

Two readings of the template's frame chain between calls:
  G (gauged): our Lean gluing is per family (Phi_k(1) = Phi_{k+1}(t), Phi_k(t^perp) = Phi_{k+1}(0) with t the
    family's own port; an-build calib.py).  Two-sided, family b ends a call with X_b = C_b + <u_b>,
    Y_b = C_b + <v_b>, C_b = u_b^perp & v_b^perp, and the next call re-identifies C_b with the background and u_b, v_b
    with its own local lines.  So EVERY call starts with X_b on <u_b>, Y_b on <v_b>, X_b inside v_b^perp throughout.
  F (fixed s): s is one global coordinate that is never re-identified (the shared-s reading of an-build / an2-word).
    Y_b ends its first call on u_b^perp, which contains s; afterwards the window part g_b of its line is background,
    so in each later call Y_b starts holding s locally and X_b never holds it.
Schedule: GF at stages 0,2,4,6,8, GB at 1,3,5,7 (template Gates.lean, as read by v-an2-word); PAIR = 1,1,1,2,2,3,3,4,4."""
import gzip, json, os
P = os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs/U8X-parsed.json.gz'))
d = json.load(gzip.open(P)); v = d['v']; U = d['ports']
S = 1 << 9; J9 = (1 << 9) - 1
dot = lambda a, b: bin(a & b).count('1') & 1
V = [S | (J9 ^ u) for u in U]
PAIR = [1, 1, 1, 2, 2, 3, 3, 4, 4]; KIND = ['GF' if k % 2 == 0 else 'GB' for k in range(9)]

def span(vs):
    out = {0}
    for x in vs:
        out |= {y ^ x for y in out}
    return out

def counts(ystart, xperp_line):
    """ystart[b]: set of vectors in Y_b's local start frame; X_b lies in xperp_line[b]^perp throughout the call
    (so a vector w is never in X_b iff w . xperp_line[b] = 1).  The same set answers both formulas here, because
    X_b's end frame is exactly xperp_line[b]^perp."""
    c = {}
    for b in range(v):
        for w in ystart[b]:
            if w and dot(w, xperp_line[b]) == 1:
                c[w] = c.get(w, 0) + 1
    return c

rows = []; first = {}
for k in range(9):
    p = PAIR[k]
    for model in ('G', 'F'):
        if model == 'G' or k == PAIR.index(p):
            ys = [span([V[b]]) for b in range(v)]; xl = V
        else:                                   # F, a later call of the pair: Y holds s, X avoids s
            ys = [span([S]) for b in range(v)]; xl = [S] * v
        c = counts(ys, xl)
        rows.append(dict(stage=k, call=KIND[k], pair=p, model=model, lines_with_count=len(c),
                         max_count=max(c.values()), total=sum(c.values())))
for r in rows:
    print('stage %(stage)d %(call)s pair %(pair)d  model %(model)s : %(lines_with_count)3d lines need a crossing, '
          'largest count per line %(max_count)3d, sum %(total)d' % r)
perrep = {m: sum(r['max_count'] for r in rows if r['model'] == m) for m in 'GF'}
print('per replica, the sum over calls of the LARGEST per-line count (a floor on paid register-crossings):', perrep)
print('gain of the second line per replica: 896.  Model F floor per replica:',
      sum(r['max_count'] for r in rows if r['model'] == 'F'), '(kills: the shared-s lemma returns in every later call)')
json.dump(rows, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'counts.json'), 'w'), indent=1)
