"""v-x-exit-words-reprice: does the Lean engine's derived undo gate (g8 Mti: "slots := rows of the inverse",
Work/Carrier/Live.lean, used in GLab.fwd_erasedE at OLD (full) frames) touch the exit slots of a gcert/1 file?
Replays the word symbolically (every register its own symbol), then the inverse word, to get Mt and the slot rows
of Mti.  For each slot q: dirty (Mt row q != e_q)?  partners j (Mti(q, j) != 0, j != q).  For an exit slot
(final frame != full), the frame the undo of q needs: span(end-of-B frame of q, end-of-B frames of its partners).
Reports how many exit slots are dirty and whether their exit frame contains what the undo needs.
usage: undo_rows.py GCERT.json.gz gx.py"""
import sys, json, gzip, importlib.util
from fractions import Fraction as Fr
from collections import Counter
cp, gxp = sys.argv[1], sys.argv[2]
spec = importlib.util.spec_from_file_location('gx', gxp); gx = importlib.util.module_from_spec(spec); spec.loader.exec_module(gx)
c = json.loads(gzip.decompress(open(cp, 'rb').read()))
v, R, h = c['v'], c['R'], c['h']; n = 2 * v + R
fr = [tuple(f) for f in c['frames']]; dim = [len(f) for f in fr]
cls = lambda r: 'x' if r < v else 'y' if r < 2 * v else 's'

def axpy(reg, t, src, co):
    d = reg[t]
    for k, u in src.items():
        w = d.get(k, 0) + co * u
        if w: d[k] = w
        else: d.pop(k, None)

def run(reg, gates, sign=1, rev=False):
    for g in (reversed(gates) if rev else gates):
        if g[0] == 'neg':
            reg[g[2]] = {k: -u for k, u in reg[g[2]].items()}
        elif g[0] == 'out':
            src = dict(reg[g[2]])
            assert g[2] not in [t for t, _, _ in g[3]]
            for t, a, b in g[3]: axpy(reg, t, src, sign * Fr(a, b))
        else:
            assert g[2] not in [s for s, _, _ in g[3]]
            for s, a, b in g[3]: axpy(reg, g[2], dict(reg[s]), sign * Fr(a, b))

def scatter(reg, sign=1):
    tot = {k: dict(reg[s]) for k, s, _ in c['ret']}
    for t in range(v):
        for k, co in gx.scat_row(c, t): axpy(reg, v + t, tot[k], sign * co)

# forward: Mt rows
reg = [{r: Fr(1)} for r in range(n)]
run(reg, c['A']); scatter(reg); run(reg, c['B'])
Mt = reg
# x rows unit? y rows?
assert all(Mt[t] == {t: 1} for t in range(v)), 'x not restored'
# inverse: Mti rows
reg = [{r: Fr(1)} for r in range(n)]
run(reg, c['B'], -1, True); scatter(reg, -1); run(reg, c['A'], -1, True)
Mti = reg
# end-of-B frame of every register (before the final climb): replay frames
cur = list(c['start'])
Bcore = c['B'][:-n] if all(g[0] == 'out' and not g[3] and not g[4] for g in c['B'][-n:]) and [g[2] for g in c['B'][-n:]] == list(range(n)) else c['B']
print('stripped final frame-only gates:', len(c['B']) - len(Bcore))
for g in c['A'] + Bcore:
    for r in gx.gate_regs(g)[0]: cur[r] = g[1]
endB = cur
def span_dim(vecs):
    basis = []
    for x in vecs:
        for b in basis: x = min(x, x ^ b)
        if x: basis.append(x)
    return len(basis)
out = Counter(); bad = []; needd = Counter(); extra = 0; extra_hist = Counter()
for q in range(2 * v, n):
    dirty = Mt[q] != {q: 1}
    part = [j for j, a in Mti[q].items() if j != q]
    assert all(cls(j) != 'y' for j in part), 'undo row reads y'
    ex = c['final'][q]
    full = (ex == 1)
    out[('full' if full else 'exit', 'dirty' if dirty else 'clean')] += 1
    if not full and dirty:
        need = list(fr[endB[q]]) + [u for j in part for u in fr[endB[j]]]
        dn = span_dim(need + list(fr[ex]))
        ok = dn == dim[ex]
        needd[(dim[ex], dn)] += 1
        if not ok:
            bad.append((q, dim[c['start'][q]], dim[ex], dn, len(part), Counter(cls(j) for j in part)))
            extra += dn - dim[ex]; extra_hist[dn - dim[ex]] += 1
print('file', cp)
print('slots by (exit kind, dirty after the word):', dict(out))
print('dirty exit slots: (dim natural/given exit, dim of span(exit, end-of-B frames of undo partners)):', dict(sorted(needd.items())))
print('dirty exit slots whose exit frame does NOT contain the undo partners\' frames:', len(bad), ' total missing rank', extra, dict(extra_hist))
for b in bad[:8]: print('  slot', b[0], 'dim ent', b[1], 'dim exit', b[2], 'needs', b[3], 'partners', b[4], dict(b[5]))

# ---- corrected exits under the Lean undo constraint (GLabX / GLab): an exit slot's exit frame must contain the
# end-of-B frames of its undo partners (J1, per row), and of every role in its tail component (J2: closure through
# exit-slot partners; here every slot is treated as an exit slot, as the natural-exit pricing does).
import os
def span_basis(vecs):
    basis = []
    for x in vecs:
        for b in basis: x = min(x, x ^ b)
        if x: basis.append(x)
    return basis
part_of = {q: [j for j in Mti[q] if j != q] for q in range(2 * v, n)}
J1 = {}; J2 = {}
for q in range(2 * v, n):
    J1[q] = len(span_basis(list(fr[c['final'][q]]) + [u for j in part_of[q] for u in fr[endB[j]]]))
# closure (undirected components of the partner graph restricted so that propagation goes through slots)
adj = {q: set() for q in range(n)}
for q in range(2 * v, n):
    for j in part_of[q]: adj[q].add(j); adj[j].add(q)
seen = {}
for q in range(2 * v, n):
    if q in seen: continue
    comp = []; st = [q]; seen[q] = q
    while st:
        a = st.pop(); comp.append(a)
        for b in adj[a]:
            if b not in seen and b >= 2 * v: seen[b] = q; st.append(b)
            elif b < 2 * v: comp.append(b)
    d = len(span_basis([u for a in comp for u in fr[endB[a] if a < 2 * v or True else 0]] + [u for a in comp if a >= 2*v for u in fr[c['final'][a]]]))
    for a in comp:
        if a >= 2 * v: J2[a] = d
ent = {q: dim[c['start'][q]] for q in range(2 * v, n)}
nat = {q: dim[c['final'][q]] for q in range(2 * v, n)}
blocks = c['blocks']; tot = Counter()
for p in blocks.values(): tot.update({int(k): v_ for k, v_ in p.items()})
res = {}
for lab, J in (('J1_row', J1), ('J2_component', J2)):
    hist = Counter(tot)
    for q in range(2 * v, n):
        if J[q] > nat[q]: hist[J[q] - nat[q]] += 1
    widths = [J[q] - ent[q] for q in range(2 * v, n)]
    res[lab] = dict(hist=dict(sorted(hist.items())), widths=widths, exit_dims=dict(sorted(Counter(J.values()).items())))
    print(lab, 'exit dims', res[lab]['exit_dims'], 'residual sum', sum(widths), '(natural', sum(nat[q] - ent[q] for q in nat), ')')
if len(sys.argv) > 3:
    json.dump(res, open(sys.argv[3], 'w'))
