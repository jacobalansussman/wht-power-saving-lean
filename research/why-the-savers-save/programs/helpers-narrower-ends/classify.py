#!/usr/bin/env python3
"""d-hpe-multigate STEP 1: exact classification of the U8 FULL-ended helpers for hyperplane ends.
For every helper q that ends at FULL in U8 (378, of which 272 make their last climb 8->9 in cleanup):
  - its cleanup gates (FULL cleanup after the 171 early/tail gates);
  - single-tail test of #591 (Kothe): exactly one cleanup gate, neither operand touched by an earlier cleanup gate;
  - multi-gate commute test: all its cleanup gates can move to right after B together with the least set of earlier
    cleanup gates they depend on (closure in the read/write conflict order; additions into one register commute);
  - U1 = join of the operands' frames at the first of its gates, U_any = join of q's and every partner's frame at the
    start of the cleanup (a lower bound on q's end for ANY order and frame choice of the SAME gate set);
  - Umin from cones.py (lower bound for ANY cleanup at all, given A + cut + B).
Usage: python3 -I classify.py FULL-WORD.json.gz cones.json out.json"""
import sys, gzip, json
from collections import Counter, defaultdict
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from cones_lib import basis, contained
w = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read())); cones = json.load(open(sys.argv[2]))['helpers']
v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2*v+R
full = basis([1 << i for i in range(h)])
def registers(e): return e[1:3] if e[0] == 'add' else e[1:2]
start = [(p,) for p in w['ports']]+[()]*(v+R)
for q, U in w['start_extra'].items(): start[int(q)] = basis(U)
cur = start.copy(); last = {}
def climb(q, U, ph):
    U = basis(U)
    if U != cur[q]: last[q] = (len(cur[q]), len(U), ph)
    cur[q] = U
for section in ('A', 'B'):
    for e in w[section]:
        for q in registers(e): climb(q, e[-1], section)
early = w['early_cleanup']; earlyids = {i for i, U in early}
for i, U in early:
    for q in registers(w['cleanup'][i]): climb(q, U, 'tail')
pre = cur.copy()   # frames at the start of the FULL cleanup
CL = [e for i, e in enumerate(w['cleanup']) if i not in earlyids]
cg = defaultdict(list)
for j, e in enumerate(CL):
    for q in registers(e): cg[q].append(j)
writes_at = defaultdict(list); reads_at = defaultdict(list); scales_at = defaultdict(list); pred = []
for j, e in enumerate(CL):
    if e[0] == 'add':
        a, b = e[1], e[2]
        pred.append(set(writes_at[b]) | set(scales_at[b]) | set(reads_at[a]) | set(scales_at[a]))
        writes_at[a].append(j); reads_at[b].append(j)
    else:
        a = e[1]; pred.append(set(writes_at[a]) | set(reads_at[a])); scales_at[a].append(j)
def closure(S):
    seen = set(S); st = list(S)
    while st:
        for i in pred[st.pop()]:
            if i not in seen: seen.add(i); st.append(i)
    return seen
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
lowend = {int(q) for q in w['helper_final']}
rows = {}; tab = Counter()
for q in range(2*v, n):
    if q in lowend: continue
    G = cg[q]; touched_before = set()
    for j in range(min(G)):
        touched_before |= set(registers(CL[j]))
    partners = [r for j in G for r in registers(CL[j]) if r != q]
    clo = closure(G)
    e0 = CL[G[0]]
    U1 = basis([x for r in registers(e0) for x in pre[r]])
    Uany = basis([x for r in [q]+partners for x in pre[r]])
    # same gate set, closure moved first, every gate at the join of its operands' current frames
    f = {}
    for j in sorted(clo):
        e = CL[j]
        if e[0] == 'add':
            jn = basis(list(f.get(e[1], pre[e[1]]))+list(f.get(e[2], pre[e[2]]))); f[e[1]] = jn; f[e[2]] = jn
    Uclo = f.get(q, pre[q])
    single = len(G) == 1 and not (set(registers(CL[G[0]])) & touched_before)
    rec = dict(start_dim=len(start[q]), cold=len(start[q]) == 0, dim_at_cleanup=len(pre[q]),
               last_climb=last.get(q), eight_to_nine_in_cleanup=len(pre[q]) == 8,
               cleanup_gates=len(G), partner_classes=dict(Counter(cls(r) for r in partners)),
               partners_touched_by_earlier_cleanup=sum(1 for j in G for r in registers(CL[j]) if r != q and any(r in registers(CL[i]) for i in range(j))),
               single_tail_591=single, commutes_alone=(clo == set(G)), closure_size=len(clo),
               U_first_gate_dim=len(U1), U_any_order_dim=len(Uany), U_closure_first_dim=len(Uclo),
               Umin_cone_dim=cones[str(q)]['Umin_dim'], Umin_row_dim=cones[str(q)]['Umin_row_dim'], Umin_col_dim=cones[str(q)]['Umin_col_dim'])
    rows[q] = rec
    tab[(rec['eight_to_nine_in_cleanup'], rec['cold'], rec['single_tail_591'], rec['commutes_alone'], rec['U_any_order_dim'], rec['Umin_cone_dim'])] += 1
S = [q for q in rows if rows[q]['eight_to_nine_in_cleanup']]
print('FULL-ended helpers', len(rows), ' of which 8->9 in cleanup', len(S))
print('single-tail (#591 rule) candidates:', sum(rows[q]['single_tail_591'] for q in rows))
print('cleanup gates per helper: min', min(rows[q]['cleanup_gates'] for q in rows), 'max', max(rows[q]['cleanup_gates'] for q in rows))
print('helpers whose gates commute to right after B with no other gate:', sum(rows[q]['commutes_alone'] for q in rows))
print('U of the first gate alone, dim:', dict(Counter(rows[q]['U_first_gate_dim'] for q in S)), '(the 272)')
print('U over all partners (any order, same gates), dim:', dict(Counter(rows[q]['U_any_order_dim'] for q in rows)))
print('U with closure moved first at joins, dim:', dict(Counter(rows[q]['U_closure_first_dim'] for q in rows)))
print('Umin from the cones (any cleanup), dim:', dict(Counter(rows[q]['Umin_cone_dim'] for q in rows)))
print('table (8->9, cold, single, commutes, Uany, Umin):')
for k, x in sorted(tab.items()): print('  ', k, x)
json.dump(dict(rows={str(q): r for q, r in rows.items()}, table={str(k): x for k, x in tab.items()}), open(sys.argv[3], 'w'), indent=0)
