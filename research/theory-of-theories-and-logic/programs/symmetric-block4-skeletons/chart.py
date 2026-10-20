"""chart.py NAME TMO : exact decision of the full system of array 0 by charts (my own enumeration).
Forced zeros (out/NAME.fz, exact over Q) dropped.  For every move: option k = 'first nonzero remaining summand is k'
(earlier remaining summands 0, summand k = 1: WLOG by rescaling the move and absorbing the factor in later alphas,
which stay free), or 'whole move zero'.  Every chart: relaxation std over Q (GB = 1 -> pruned, Lemma 1), else full
closure system, constant-pivot elimination, std over Q.  Prints a per-chart verdict and a final count."""
import sys, os, json, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import vnr
from vnr_run import sing
name, tmo = sys.argv[1], int(sys.argv[2])
rec = json.load(open(os.path.join(HERE, 'recipes', name + '.json')))
fz = set(json.load(open(os.path.join(HERE, 'out', name + '.fz'))))
OUT = os.path.join(HERE, 'out', 'ch_' + name); os.makedirs(OUT, exist_ok=True)
KJ = vnr.Skel(rec).keptj
moves = []; a0 = 0
for (kind, summ, Dm) in rec['moves']:
    idx = [a0 + k for k in range(len(summ)) if a0 + k not in fz]
    moves.append(idx); a0 += len(summ)
opts = [list(range(len(ix))) + ['z'] for ix in moves]
tot = {'pruned': 0, 'excluded': 0, 'open': 0, 'FEASIBLE': 0}
for ch in itertools.product(*opts):
    zero = set(fz); one = set()
    for ix, o in zip(moves, ch):
        if o == 'z': zero.update(ix)
        else: zero.update(ix[:o]); one.add(ix[o])
    tag = ''.join(str(o) for o in ch)
    sk = vnr.Skel(rec, zero=zero, one=one, keptj=KJ)
    txt, av, bv = vnr.relax_singular(sk, tag)
    r, dt = sing(txt, os.path.join(OUT, tag + '.rel.sing'), 300)
    if r == 'GB=1':
        tot['pruned'] += 1; continue
    cl = sk.closure(0, 60000)
    if cl is None:
        print('CHART', name, tag, 'relax', r, 'closure > cap -> OPEN', flush=True); tot['open'] += 1; continue
    A, nc = cl
    core = vnr.const_eliminate(A, nc)
    if any(set(R) == {'r'} and set(R['r']) == {()} for R in core):
        res = 'GB=1'; dt2 = 0
    else:
        t2, av2, bv2 = vnr.core_singular(core, tag)
        res, dt2 = sing(t2, os.path.join(OUT, tag + '.full.sing'), tmo)
    k = 'excluded' if res == 'GB=1' else ('FEASIBLE' if res == 'GB!=1' else 'open')
    tot[k] += 1
    print('CHART', name, tag, 'relax', r, 'closure', len(A), 'x', nc, 'core', len(core), 'full', res, round(dt2, 1), flush=True)
print('CHARTS', name, 'forced', sorted(fz), tot, flush=True)
