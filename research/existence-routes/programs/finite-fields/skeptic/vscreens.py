"""vscreens.py -- v-ex2-modp: is every screen / chart-relaxation system that ex2-modp certified exactly decide.py's
system: the relaxation (u export, checked = exact build) plus the rows decide.py adds for that step?"""
import os, sys, json, re, ast, glob, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vlib import *
M = os.environ.get('MODP_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
def canon(names, P):
    return tuple(sorted((tuple(sorted((names[i], e) for i, e in enumerate(k) if e)), c) for k, c in P.items()))
tot = good = 0
for D in sorted(glob.glob(f'{M}/data/dec/*')) + sorted(glob.glob(f'{M}/data/hard/dec/*')):
    sk = os.path.basename(D)
    if not os.path.exists(f'{D}/log.txt') and not os.path.exists(f'{D}/verdict.json'): continue
    rec = json.load(open(os.path.join(os.environ.get('NAG_RERUN_DIR', 'nag-rerun'), 'recipes', f'{sk}.json')))
    ns = [len(s) for (_, s, _) in rec['moves']]; idx = []; b = 0
    for n in ns: idx.append(list(range(b, b + n))); b += n
    if os.path.exists(f'{D}/verdict.json'):
        Z = json.load(open(f'{D}/verdict.json'))['forced_zero']
    else:
        Z = ast.literal_eval(re.search(r'forced zeros \(0-based original alpha indices\): (\[.*?\])', open(f'{D}/log.txt').read()).group(1))
    K = [[k for k in range(n) if idx[t][k] not in Z] for t, n in enumerate(ns)]
    def chart_eqs(t, opt):
        eqs = [('z', i) for i in Z]
        if opt == 'zero': eqs += [('z', idx[t][k]) for k in K[t]]
        else:
            for k in K[t]:
                if k == opt: break
                eqs.append(('z', idx[t][k]))
            eqs.append(('o', idx[t][opt]))
        return eqs
    for j, p in ((1, P1), (2, P2)):
        un, _, U = read_ms(f'{M}/data/ms/{sk}.u.p{j}.ms')
        for f in sorted(glob.glob(f'{D}/*.p{j}.ms')):
            b = os.path.basename(f)[:-6]
            if b.startswith('full_'): continue
            if b.startswith('fz'): want = [('t', int(b[2:]) - 1)]
            elif b.startswith('m'):
                t, o = b[1:].split('_', 1); want = chart_eqs(int(t), o if o == 'zero' else int(o))
            elif b.startswith('c_'):
                combo = b[2:].split('_'); want = sorted(set(sum((chart_eqs(t, o if o == 'zero' else int(o)) for t, o in enumerate(combo)), [])))
            elif b == 'rel_u': want = []
            else: continue
            names, pp, P = read_ms(f)
            exp_names = un + (['t'] if any(w[0] == 't' for w in want) else [])
            rows = []
            for kind, i in want:
                e = [0] * len(exp_names); e[exp_names.index(f'a{i+1}')] = 1
                if kind == 'z': rows.append({tuple(e): 1})
                elif kind == 'o': rows.append({tuple(e): 1, tuple([0] * len(exp_names)): p - 1})
                else:
                    e[exp_names.index('t')] = 1; rows.append({tuple(e): 1, tuple([0] * len(exp_names)): p - 1})
            Uc = [canon(un, Q) for Q in U]
            PC = [canon(names, Q) for Q in P]
            ok = names == exp_names and PC[:len(Uc)] == Uc and sorted(PC[len(Uc):]) == sorted(canon(exp_names, Q) for Q in rows)
            tot += 1; good += ok
            if not ok: print('MISMATCH', sk, b, j)
print(f'TOTAL {good}/{tot} screen/chart-relaxation files equal relaxation + decide.py rows')
