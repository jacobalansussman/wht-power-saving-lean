"""d-halfread: exact obstruction census for retiring a U8 second copy by pre-merge half reads.
(1) targets: each target's single rank-2 stop is already a pre-merge half read of a plane outside its cap;
    for every (plane P, target t with P inside u_t^perp) pair, is t's rank-2 line a half line of P?
(2) per copy: is the partner b still holding beta when the copy's readers read it; reader classes.
Usage: python3 -I obstruction.py WORD copies.json screen.json out.json"""
import sys, json
from collections import Counter, defaultdict
sys.path.insert(0, __import__('os').path.dirname(__file__))
from wordlib import *
w = load(sys.argv[1]); v, R = w['report']['v'], w['report']['R']; n = 2*v+R
prog = program(w); chain, bad = frame_replay(w, prog); C = json.load(open(sys.argv[2])); S = json.load(open(sys.argv[3]))
start, end = ends(w)
out = {}
# (1) targets
half_of = {}  # target -> (plane, 'a'/'b', reg)
for r in C:
    Tm, _, a, b, k = r['merge']
    for T, (ph, e) in enumerate(prog[:Tm]):
        if e[0] == 'add' and e[2] in (a, b) and v <= e[1] < 2*v and ph == 'B' and len(basis(e[-1])) == 2:
            half_of.setdefault(e[1], []).append((tuple(r['P']), 'a' if e[2] == a else 'b', e[2], T))
rank2 = {t: [U for T, U in chain[t] if len(U) == 2] for t in range(v, 2*v)}
caps = {t: end[t] for t in range(v, 2*v)}
pairs = 0; line_is_half = 0; own_half_in_cap = 0
for r in C:
    P = tuple(r['P']); La = tuple(r['La']); Lb = tuple(r['Lb'])
    for t in range(v, 2*v):
        if contained(P, caps[t]):
            pairs += 1
            L = rank2[t][0]
            if L == La or L == Lb: line_is_half += 1
for t, hs in half_of.items():
    for P, ab, reg, T in hs:
        if contained(P, caps[t]): own_half_in_cap += 1
out['targets_with_one_rank2_stop'] = sum(1 for t in rank2 if len(rank2[t]) == 1)
out['targets_whose_rank2_stop_is_a_premerge_half_read'] = sum(1 for t in half_of if len(half_of[t]) == 1)
out['halves_read_by_targets'] = dict(Counter(ab for hs in half_of.values() for _, ab, _, _ in hs))
out['half_reads_of_a_plane_inside_the_readers_cap'] = own_half_in_cap
out['plane_target_pairs_with_P_in_cap'] = pairs
out['of_which_target_rank2_line_is_a_half_line_of_P'] = line_is_half
# (2) copies: partner b intact during the copy's reads, and reader classes in that case
gl = defaultdict(list)
for T, (ph, e) in enumerate(prog):
    if ph != 'compensation':
        for q in regs(e): gl[q].append(T)
cls = lambda q: 'target' if v <= q < 2*v else ('source' if q < v else ('entering helper' if start[q] else 'cold helper'))
bstat = Counter(); readers_when_b_intact = Counter(); plane_list = []
for r, s in zip(C, S):
    Tm, phm, a, b, k = r['merge']
    fw = next((T for T in gl[b] if T > Tm and prog[T][1][1] == b), None)
    rd = [(T, c) for T, ph, c, cl, qq, d in r['reads'] if r['copy'][0] < T < r['Tend']]
    intact = all(fw is None or T < fw for T, c in rd)
    bstat[('merge in ' + phm, 'b keeps beta through the copy reads' if intact else 'b rewritten before a copy read (in %s)' % prog[fw][0])] += 1
    if intact:
        readers_when_b_intact[tuple(sorted(cls(c) for c in {c for _, c in rd}))] += 1
        plane_list.append(dict(P=r['P'], h3=r['copy'][2], readers=sorted({c for _, c in rd}), verdicts=s['verdicts']))
out['partner_b'] = {' / '.join(k): x for k, x in bstat.items()}
out['reader_classes_where_b_keeps_beta'] = {', '.join(k): x for k, x in readers_when_b_intact.items()}
json.dump(dict(out, planes_b_intact=plane_list), open(sys.argv[4], 'w'), indent=1)
print(json.dumps(out, indent=1))
